"""`database_identity` must be able to FAIL: every way two databases can differ, and every way it can look.

The instrument is the gate for "this change moves no number", so each test here builds two small
schemas in the probe database that differ in exactly one way and asserts the instrument says so —
and, for the three things it deliberately ignores, that it ignores exactly those.
"""

from __future__ import annotations

import pytest
from sqlalchemy import create_engine, text

from bacatlas_backend.instruments.database_identity import (
    BUILD_STAMPS,
    CoverageMismatch,
    compare_databases,
)
from tests.conftest import PROBE_URL

REFERENCE, CANDIDATE = "identity_reference", "identity_candidate"

#: Two tables: one ordinary, one carrying the build stamp and a timestamp.
LOCUS_ROWS = [(1, "a", b"\x00\x01"), (2, "b", b"\x02"), (3, None, None)]
PANGENOME_ROWS = [(1, "run-1", "abc123", "2026-01-01T00:00:00Z")]


@pytest.fixture(scope="module")
def probe():
    engine = create_engine(PROBE_URL, future=True)
    try:
        with engine.connect() as connection:
            connection.execute(text("SELECT 1"))
    except Exception as error:  # noqa: BLE001
        pytest.skip(f"probe database unavailable at {PROBE_URL}: {error}")
    yield engine
    with engine.begin() as connection:
        for schema in (REFERENCE, CANDIDATE):
            connection.execute(text(f"DROP SCHEMA IF EXISTS {schema} CASCADE"))


def build(engine, schema, *, locus=LOCUS_ROWS, pangenome=PANGENOME_ROWS, label_type="text", with_locus=True):
    with engine.begin() as connection:
        connection.execute(text(f"DROP SCHEMA IF EXISTS {schema} CASCADE"))
        connection.execute(text(f"CREATE SCHEMA {schema}"))
        connection.execute(text(
            f"CREATE TABLE {schema}.pangenome (pangenome_id int, run_id text, git_sha text, ingested_at timestamptz)"
        ))
        for row in pangenome:
            connection.execute(
                text(f"INSERT INTO {schema}.pangenome VALUES (:a, :b, :c, :d)"),
                dict(zip("abcd", row, strict=True)),
            )
        if with_locus:
            connection.execute(text(f"CREATE TABLE {schema}.locus (locus_id int, label {label_type}, block bytea)"))
            for row in locus:
                connection.execute(
                    text(f"INSERT INTO {schema}.locus VALUES (:a, :b, :c)"), dict(zip("abc", row, strict=True))
                )


def compare(engine, **keywords):
    return compare_databases(
        engine, engine, reference_schema=REFERENCE, candidate_schema=CANDIDATE,
        must_hold_rows=frozenset({"locus", "pangenome"}), **keywords,
    )  # fmt: skip


def test_two_identical_schemas_agree_and_the_coverage_is_counted(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE)
    report = compare(probe)
    assert report.identical
    assert [table.table for table in report.tables] == ["locus", "pangenome"]
    assert report.column_count == 7
    assert sum(table.reference_rows for table in report.tables) == 4
    assert "0 table(s) differ, 0 required table(s) empty" in report.render()


def test_one_changed_value_is_one_row_on_each_side(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, locus=[(1, "a", b"\x00\x01"), (2, "B", b"\x02"), (3, None, None)])
    report = compare(probe)
    assert [table.table for table in report.differing] == ["locus"]
    locus = report.differing[0]
    assert (locus.only_in_reference, locus.only_in_candidate) == (1, 1)
    assert "b" in locus.reference_examples[0] and "B" in locus.candidate_examples[0]
    assert not report.identical


def test_one_changed_byte_in_a_bytea_block_is_a_difference(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, locus=[(1, "a", b"\x00\x02"), (2, "b", b"\x02"), (3, None, None)])
    assert [table.table for table in compare(probe).differing] == ["locus"]


def test_a_null_and_an_empty_string_are_different_rows(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, locus=[(1, "a", b"\x00\x01"), (2, "b", b"\x02"), (3, "", None)])
    assert [table.table for table in compare(probe).differing] == ["locus"]


def test_a_duplicated_row_is_a_difference_and_is_counted_once(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, locus=[*LOCUS_ROWS, LOCUS_ROWS[0]])
    locus = compare(probe).differing[0]
    assert (locus.reference_rows, locus.candidate_rows) == (3, 4)
    assert (locus.only_in_reference, locus.only_in_candidate) == (0, 1)


def test_the_order_rows_were_inserted_in_is_not_a_difference(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, locus=LOCUS_ROWS[::-1])
    assert compare(probe).identical


def test_a_timestamp_and_a_build_stamp_are_excluded_BY_NAME_and_nothing_else_is(probe):
    assert ("pangenome", "git_sha") in BUILD_STAMPS
    build(probe, REFERENCE)
    build(probe, CANDIDATE, pangenome=[(1, "run-1", "fff999", "2027-06-06T00:00:00Z")])
    report = compare(probe)
    assert report.identical
    pangenome = next(table for table in report.tables if table.table == "pangenome")
    assert pangenome.excluded_columns == ("git_sha", "ingested_at")
    assert pangenome.compared_columns == ("pangenome_id", "run_id")
    # the same table, differing in a column that is NOT excluded
    build(probe, CANDIDATE, pangenome=[(1, "run-2", "abc123", "2026-01-01T00:00:00Z")])
    assert [table.table for table in compare(probe).differing] == ["pangenome"]


def test_a_missing_table_is_refused_before_any_row_is_compared(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, with_locus=False)
    with pytest.raises(CoverageMismatch, match="different tables"):
        compare(probe)


def test_a_column_of_another_type_is_refused(probe):
    build(probe, REFERENCE)
    build(probe, CANDIDATE, label_type="varchar(8)")
    with pytest.raises(CoverageMismatch, match="locus"):
        compare(probe)


def test_two_EMPTY_required_tables_agree_and_are_still_a_failure(probe):
    build(probe, REFERENCE, locus=[])
    build(probe, CANDIDATE, locus=[])
    report = compare(probe)
    assert not report.differing
    assert report.empty_required == ["locus"]
    assert not report.identical
    assert "EMPTY   locus" in report.render()


def test_a_schema_with_no_table_is_refused_rather_than_reported_identical(probe):
    with probe.begin() as connection:
        for schema in (REFERENCE, CANDIDATE):
            connection.execute(text(f"DROP SCHEMA IF EXISTS {schema} CASCADE"))
            connection.execute(text(f"CREATE SCHEMA {schema}"))
    with pytest.raises(CoverageMismatch, match="no table"):
        compare(probe)
