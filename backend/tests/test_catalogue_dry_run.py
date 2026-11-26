"""The dry run must be able to say NOT CLEAN — for a row the schema will not hold, and for a refusal.

An instrument that builds a catalogue, sends nothing and reports "clean" is worth exactly as much as the
proof that it would have said otherwise. So: each kind of value Postgres would reject is put through
`SchemaFit` and named; the published Kp probe runs clean with the very row counts its real load writes;
and one genome's `.loci` is corrupted by a single ordinal and the run ends 1, naming the genome.
"""

from __future__ import annotations

import enum
import os
import struct

import pytest
from sqlalchemy import (
    ARRAY,
    Boolean,
    Column,
    Enum,
    Float,
    Index,
    Integer,
    LargeBinary,
    MetaData,
    SmallInteger,
    String,
    Table,
    Text,
)

from bacatlas_backend.instruments import catalogue_dry_run
from bacatlas_backend.instruments.catalogue_dry_run import DrySession, NeedsDatabase, SchemaFit
from tests.conftest import (
    NUNA_DATA_ROOT,
    artifacts_for,
    published_catalogue,
    redesigned_audit_artifacts,
    requires_artifacts,
)


class Colour(enum.Enum):
    RED = "red"
    BLUE = "blue"


TABLE = Table(
    "probe_rows",
    MetaData(),
    Column("id", Integer, primary_key=True, autoincrement=True),
    Column("label", String(4), nullable=False),
    Column("rank", SmallInteger, nullable=True),
    Column("score", Float, nullable=True),
    Column("block", LargeBinary, nullable=True),
    Column("colour", Enum(Colour), nullable=True),
    Column("members", ARRAY(Integer), nullable=True),
)
COLUMNS = ("label", "rank", "score", "block", "colour", "members")
GOOD = ("abcd", 7, 0.5, b"\x00\x01", Colour.RED, [1, 2, 3])


def _kinds(fit: SchemaFit) -> dict:
    return {(column, kind): count for (_table, column, kind), count in fit.violations.items()}


def test_rows_that_fit_are_counted_and_nothing_is_reported():
    fit = SchemaFit()
    assert fit.copy_rows(None, TABLE, COLUMNS, [GOOD, ("ab", None, float("nan"), None, None, [])]) == 2
    assert fit.rows == {"probe_rows": 2}
    assert not fit.violations, "NaN is coerced to NULL exactly as the real loader coerces it"


@pytest.mark.parametrize(
    "row,column,kind",
    [
        (("abcde", 7, 0.5, b"x", Colour.RED, [1]), "label", "longer than String(4)"),
        ((None, 7, 0.5, b"x", Colour.RED, [1]), "label", "NULL in a NOT NULL column"),
        (("ab\x00", 7, 0.5, b"x", Colour.RED, [1]), "label", "a NUL byte, which text cannot hold"),
        ((5, 7, 0.5, b"x", Colour.RED, [1]), "label", "not text in a text column"),
        (("ab", 40_000, 0.5, b"x", Colour.RED, [1]), "rank", "outside ±32,768"),
        (("ab", 1.5, 0.5, b"x", Colour.RED, [1]), "rank", "not an integer in an integer column"),
        (("ab", 7, float("inf"), b"x", Colour.RED, [1]), "score", "NaN or infinity"),
        (("ab", 7, 0.5, b"x", "GREEN", [1]), "colour", "not one of the enum's 2 names"),
        (("ab", 7, 0.5, b"x", Colour.RED, [1, None]), "members", "NULL element in an array"),
        (("ab", 7, 0.5, b"x", Colour.RED, [2**31]), "members", "array element out of range"),
        (("ab", 7, 0.5, b"x", Colour.RED, [1.5]), "members", "array element that is not an integer"),
        (("ab", 7, "0.5", b"x", Colour.RED, [1]), "score", "not a number in a numeric column"),
        (("ab", 7, 0.5, b"x", Colour.RED, "1,2"), "members", "not a list in an array column"),
    ],
)
def test_each_value_postgres_would_reject_is_named_with_its_column(row, column, kind):
    fit = SchemaFit()
    fit.copy_rows(None, TABLE, COLUMNS, [GOOD, row, GOOD])
    assert _kinds(fit) == {(column, kind): 1}, fit.violations
    assert fit.rows["probe_rows"] == 3, "a row that does not fit is still counted: the run goes on to find the rest"


TEXTS = Table(
    "probe_texts",
    MetaData(),
    Column("id", Integer, primary_key=True, autoincrement=True),
    Column("term", Text, nullable=True),
    Column("note", Text, nullable=True),
    Column("blurb", Text, nullable=True),
    Column("flag", Boolean, nullable=True),
    Column("names", ARRAY(String), nullable=True),
    Index("ix_probe_texts__term", "term"),
    Index("ix_probe_texts__blurb_trigram", "blurb", postgresql_using="gin"),
)


def test_unbounded_text_is_limited_ONLY_where_a_btree_holds_it():
    """`locus_annotation_entry.term_value` is Text under a btree: 2,704 bytes an entry, and COPY refuses past it."""
    fit = SchemaFit()
    long = "x" * 3000
    fit.copy_rows(None, TEXTS, ("term", "note", "blurb"), [("short", long, long), (long, "short", "short")])
    assert _kinds(fit) == {("term", "may not fit a btree index entry (2,704 bytes)"): 1}, fit.violations
    # bytes, not characters: 1,000 three-byte characters are 3,000 bytes
    fit = SchemaFit()
    fit.copy_rows(None, TEXTS, ("term",), [("語" * 1000,), ("語" * 800,)])
    assert _kinds(fit) == {("term", "may not fit a btree index entry (2,704 bytes)"): 1}


@pytest.mark.parametrize(
    "columns,row,column,kind",
    [
        (("flag",), ("yes",), "flag", "not a boolean in a boolean column"),
        (("names",), (["a", 3],), "names", "array element that is not text"),
    ],
)
def test_booleans_and_text_arrays_are_typed(columns, row, column, kind):
    fit = SchemaFit()
    fit.copy_rows(None, TEXTS, columns, [row])
    assert _kinds(fit) == {(column, kind): 1}


def test_text_in_a_bytea_column_is_refused_as_the_real_loader_refuses_it():
    with pytest.raises(TypeError):
        SchemaFit().copy_rows(None, TABLE, COLUMNS, [("ab", 7, 0.5, "not bytes", Colour.RED, [1])])


def test_a_column_the_table_does_not_have_is_refused():
    with pytest.raises(KeyError, match="no column"):
        SchemaFit().copy_rows(None, TABLE, ("label", "nope"), [("ab", 1)])


def test_a_question_only_a_database_can_answer_is_named_not_answered():
    session = DrySession()
    session.add(object())
    session.flush()
    assert session.added == {"object": 1}
    with pytest.raises(NeedsDatabase):
        session.execute("select 1")


# ── the whole thing, on the published Kp probe ──────────────────────────────────────────────────
def _arguments(nseq_dir, scratch) -> list[str]:
    """⚠ The data root is `redesigned_audit_artifacts`' copy, not the mirror: the mirrored audit predates
    2026-11-26 and is refused as it stands, and these tests are about the dry run, not the audit facts."""
    entry = published_catalogue("kp")
    data_root = redesigned_audit_artifacts("kp", scratch / "root").data_root
    return [
        "--data-root", str(data_root), "--set-key", entry.set_key, "--species-key", entry.species_key,
        "--model-label", entry.model_label, "--run-id", entry.run_id, "--nseq-dir", str(nseq_dir),
    ]  # fmt: skip


SEQ = NUNA_DATA_ROOT.parent.parent / "bacatlas" / "data" / "seq"


@requires_artifacts
@pytest.mark.skipif(not SEQ.is_dir(), reason=f"the probe's .nseq/.loci are not at {SEQ}")
def test_the_published_kp_probe_runs_CLEAN_with_the_row_counts_its_real_load_writes(capsys, tmp_path):
    assert catalogue_dry_run.main(_arguments(SEQ, tmp_path)) == 0
    out = capsys.readouterr().out
    assert "DRY RUN CLEAN: 0 refusal(s), 0 column(s)" in out
    assert "the .loci oracle compared 100 genomes byte for byte; 0 have no .loci" in out
    assert "gene counts: 100 genomes from their .nseq gene table, 0 from their _meta row count" in out
    for table, rows in (
        ("locus", "15,670"), ("locus_arrangement", "57,831"), ("locus_offset_occupant", "228,281"),
        ("locus_annotation_entry", "56,831"), ("genome_locus_assignment", "100"), ("locus_similarity", "31,340"),
    ):  # fmt: skip
        assert any(line.split() == [table, rows] for line in out.splitlines()), (table, rows)


@requires_artifacts
@pytest.mark.skipif(not SEQ.is_dir(), reason=f"the probe's .nseq/.loci are not at {SEQ}")
def test_ONE_wrong_ordinal_in_one_genomes_loci_ends_the_run_1_and_names_the_genome(tmp_path, capsys):
    from bacatlas_backend.ingest.ingest_genome_collection_roster import read_genome_vocabulary

    sample = read_genome_vocabulary(artifacts_for("kp"))[0][0]  # a genome OF THIS catalogue, not any .loci
    os.symlink(SEQ / f"{sample}.nseq", tmp_path / f"{sample}.nseq")
    block = bytearray((SEQ / f"{sample}.loci").read_bytes())
    (first,) = struct.unpack_from("<i", block, 0)
    struct.pack_into("<i", block, 0, first + 1)
    (tmp_path / f"{sample}.loci").write_bytes(bytes(block))

    code = catalogue_dry_run.main(_arguments(tmp_path, tmp_path))
    out = capsys.readouterr().out
    assert code == 1
    assert "REFUSED genome locus assignments + the .loci oracle" in out and sample in out
    assert "DRY RUN NOT CLEAN: 1 refusal(s)" in out


@requires_artifacts
@pytest.mark.skipif(not SEQ.is_dir(), reason=f"the probe's .nseq/.loci are not at {SEQ}")
def test_a_loci_of_the_WRONG_LENGTH_is_refused_because_the_count_comes_from_the_nseq(tmp_path, capsys):
    """The gene count used to be the `.loci`'s own length, so its length check compared the file with itself."""
    from bacatlas_backend.ingest.ingest_genome_collection_roster import read_genome_vocabulary

    sample = read_genome_vocabulary(artifacts_for("kp"))[0][0]
    os.symlink(SEQ / f"{sample}.nseq", tmp_path / f"{sample}.nseq")
    (tmp_path / f"{sample}.loci").write_bytes((SEQ / f"{sample}.loci").read_bytes() + struct.pack("<i", -1))

    assert catalogue_dry_run.main(_arguments(tmp_path, tmp_path)) == 1
    out = capsys.readouterr().out
    assert "REFUSED genome locus assignments + the .loci oracle" in out and sample in out
