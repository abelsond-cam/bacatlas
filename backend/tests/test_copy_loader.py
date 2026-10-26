"""The COPY loader, tested against the REAL tables — including the four coercions that are traps.

⛔ Each of these is a value that would land in Postgres looking entirely plausible and meaning
something else: a NaN that compares equal to itself, an enum written by value instead of name, an
empty array collapsed to NULL, a numpy scalar the driver cannot adapt.
"""

from __future__ import annotations

import math

import pytest
from sqlalchemy import func, select

from bacatlas_backend.ingest.staging_table_loader import copy_rows, replace_rows_for
from bacatlas_backend.models.enumerations import PrevalenceBand, SampleIdentifierKind
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.locus import Locus
from tests.conftest import TEN_ZEROS

GENE_COLUMNS = (
    "genome_id", "flat_index", "pathogen_species_id", "contig_index",
    "start_position", "end_position", "length_nt", "strand", "strand_is_observed",
    "gff_phase", "is_five_prime_partial",
)


def _gene_row(genome_id, species_id, flat_index, **overrides):
    row = dict(
        genome_id=genome_id, flat_index=flat_index, pathogen_species_id=species_id,
        contig_index=0, start_position=1, end_position=300, length_nt=300,
        strand="+", strand_is_observed=True, gff_phase=0, is_five_prime_partial=False,
    )
    row.update(overrides)
    return tuple(row[name] for name in GENE_COLUMNS)


@pytest.fixture()
def clean_genes(session, _seed_ids):
    session.execute(Gene.__table__.delete())
    session.commit()
    yield
    session.execute(Gene.__table__.delete())
    session.commit()


def test_it_copies_rows_and_reports_what_it_wrote(session, _seed_ids, clean_genes):
    written = copy_rows(
        session, Gene.__table__, GENE_COLUMNS,
        (_gene_row(_seed_ids["genome"], _seed_ids["species"], i) for i in range(1_000)),
    )
    session.commit()

    assert written == 1_000
    assert session.execute(select(func.count()).select_from(Gene)).scalar_one() == 1_000


def test_a_nan_measurement_becomes_null_and_never_reaches_the_check(session, _seed_ids):
    """⚠ Postgres defines `NaN = NaN` as TRUE, so a NaN that got in is indistinguishable from data.

    ⚠ **This test moved from `gene.gc_percent` to `locus.interest_score` on 2026-10-26**, because
    dropping that column left `gene` with no `measurement()` column and therefore no `nan_guards`.
    The subject is the coercion, not the table — but it has to run against a column that really
    carries the CHECK, or it proves nothing about either.
    """
    session.execute(Locus.__table__.delete())
    copy_rows(
        session, Locus.__table__,
        ("pangenome_id", "pathogen_species_id", "node_label", "catalogue_ordinal",
         "member_gene_count", "member_genome_count", "prevalence_band", "display_name",
         "display_name_source", "named_member_count", "context_observed_member_counts",
         "total_arrangement_count", "arrangement_member_gene_count",
         "arrangement_member_genome_count", "interest_score"),
        [(_seed_ids["pangenome"], _seed_ids["species"], "0", 0, 99, 97,
          PrevalenceBand.CORE, "wzi", "bakta_symbol", 63, TEN_ZEROS, 4, 99, 97, math.nan)],
    )
    session.commit()

    stored = session.execute(select(Locus.interest_score)).scalar_one()
    assert stored is None
    session.execute(Locus.__table__.delete())
    session.commit()


def test_an_enum_is_written_by_name_not_by_value(session, _seed_ids):
    """SQLAlchemy persists `.name`; the CHECK constraints were written against that spelling."""
    species_id = _seed_ids["species"]
    copy_rows(
        session, Genome.__table__,
        ("pathogen_species_id", "sample_id", "sample_id_kind", "strand_is_observed"),
        [(species_id, "SAMEA_COPY_TEST", SampleIdentifierKind.ASSEMBLY_STEM, False)],
    )
    session.commit()

    kind = session.execute(
        select(Genome.sample_id_kind).where(Genome.sample_id == "SAMEA_COPY_TEST")
    ).scalar_one()
    assert kind is SampleIdentifierKind.ASSEMBLY_STEM
    session.execute(Genome.__table__.delete().where(Genome.__table__.c.sample_id == "SAMEA_COPY_TEST"))
    session.commit()


def test_an_empty_array_is_not_null(session, _seed_ids):
    """⛔ `{}` is *annotated with nothing*; NULL is *not annotated*. Different claims."""
    genome_id, _ = _seed_ids["genome"], None
    session.execute(GeneFunctionalAnnotation.__table__.delete())
    copy_rows(
        session, GeneFunctionalAnnotation.__table__,
        ("genome_id", "flat_index", "gene_ontology_terms"),
        [(genome_id, 0, []), (genome_id, 1, None), (genome_id, 2, ["GO:0003677"])],
    )
    session.commit()

    stored = dict(
        session.execute(
            select(GeneFunctionalAnnotation.flat_index, GeneFunctionalAnnotation.gene_ontology_terms)
        ).all()
    )
    assert stored[0] == []
    assert stored[1] is None
    assert stored[2] == ["GO:0003677"]
    session.execute(GeneFunctionalAnnotation.__table__.delete())
    session.commit()


def test_a_numpy_scalar_is_adapted(session, _seed_ids, clean_genes):
    """`numpy.int64` has no text dumper; the failure otherwise names the driver, not the column."""
    numpy = pytest.importorskip("numpy")
    copy_rows(
        session, Gene.__table__, GENE_COLUMNS,
        [_gene_row(_seed_ids["genome"], _seed_ids["species"], numpy.int64(7),
                   start_position=numpy.int32(11))],
    )
    session.commit()

    row = session.execute(select(Gene.flat_index, Gene.start_position)).one()
    assert row == (7, 11)


def test_a_numpy_float_is_adapted_on_a_measurement_column(session, _seed_ids):
    """⚠ The float half of the numpy coercion, which left `gene` with `gc_percent`."""
    numpy = pytest.importorskip("numpy")
    session.execute(Locus.__table__.delete())
    copy_rows(
        session, Locus.__table__,
        ("pangenome_id", "pathogen_species_id", "node_label", "catalogue_ordinal",
         "member_gene_count", "member_genome_count", "prevalence_band", "display_name",
         "display_name_source", "named_member_count", "context_observed_member_counts",
         "total_arrangement_count", "arrangement_member_gene_count",
         "arrangement_member_genome_count", "interest_score"),
        [(_seed_ids["pangenome"], _seed_ids["species"], "0", 0, 99, 97,
          PrevalenceBand.CORE, "wzi", "bakta_symbol", 63, TEN_ZEROS, 4, 99, 97,
          numpy.float32(48.25))],
    )
    session.commit()

    assert session.execute(select(Locus.interest_score)).scalar_one() == pytest.approx(48.25)
    session.execute(Locus.__table__.delete())
    session.commit()


def test_a_mis_named_column_is_refused_before_a_row_is_sent(session, _seed_ids):
    with pytest.raises(KeyError) as failure:
        copy_rows(session, Gene.__table__, ("genome_id", "start"), [(1, 1)])
    assert "start" in str(failure.value)


def test_replace_refuses_an_unscoped_predicate(session, _seed_ids):
    """An unbounded `where` here is a table wipe the caller would see as a successful load."""
    with pytest.raises(ValueError, match="scoped predicate"):
        replace_rows_for(session, Gene.__table__, GENE_COLUMNS, [], where=None)


def test_replace_is_idempotent_at_its_own_scope(session, _seed_ids, clean_genes):
    """⭐ Re-running an ingest after a fix must not double the rows it already wrote."""
    genome_id, species_id = _seed_ids["genome"], _seed_ids["species"]
    rows = [_gene_row(genome_id, species_id, i) for i in range(50)]
    scope = Gene.__table__.c.genome_id == genome_id

    replace_rows_for(session, Gene.__table__, GENE_COLUMNS, rows, where=scope)
    session.commit()
    deleted, written = replace_rows_for(session, Gene.__table__, GENE_COLUMNS, rows, where=scope)
    session.commit()

    assert (deleted, written) == (50, 50)
    assert session.execute(select(func.count()).select_from(Gene)).scalar_one() == 50


def test_the_loader_issues_one_copy_statement_per_call(session, engine, _seed_ids, clean_genes):
    """⭐ The cost assertion: 5,000 rows must not be 5,000 statements."""
    from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle

    rows = [_gene_row(_seed_ids["genome"], _seed_ids["species"], i) for i in range(5_000)]
    with SqlCostOracle(engine) as report:
        copy_rows(session, Gene.__table__, GENE_COLUMNS, rows, batch_rows=1_000)
        session.commit()

    copies = [record for record in report.statements if record.category == "copy"]
    assert len(copies) == 1, report.statements


def test_a_locus_array_column_round_trips(session, _seed_ids):
    """`context_observed_member_counts` is ten integers in OFFSETS order, and order is the meaning."""
    counts = [97, 96, 95, 94, 93, 92, 91, 90, 89, 88]
    session.execute(Locus.__table__.delete())
    copy_rows(
        session, Locus.__table__,
        ("pangenome_id", "pathogen_species_id", "node_label", "catalogue_ordinal",
         "member_gene_count", "member_genome_count", "prevalence_band", "display_name",
         "display_name_source", "named_member_count", "context_observed_member_counts",
         "total_arrangement_count", "arrangement_member_gene_count",
         "arrangement_member_genome_count"),
        [(_seed_ids["pangenome"], _seed_ids["species"], "0", 0, 99, 97,
          PrevalenceBand.CORE, "wzi", "bakta_symbol", 63, counts, 4, 99, 97)],
    )
    session.commit()

    assert session.execute(select(Locus.context_observed_member_counts)).scalar_one() == counts
    session.execute(Locus.__table__.delete())
    session.commit()


# ── bytea through text-format COPY ────────────────────────────────────────────────────────────
# ⛔ **No bytea column exists in this schema today, so this path has never been exercised.** The
# gene-layer redesign (`nuna/docs/design/storing_the_gene_layer.md`) stores every per-gene field in
# `bytea` blocks, and `staging_table_loader`'s own docstring is the reason this is a gate rather
# than an assumption: the binary COPY path it rejected "fails silently by writing a plausible wrong
# value", and text format carries exactly the same risk for a type whose text representation is
# itself backslash-escaped. `\x00 \x09 \x0a \x0d \x5c` are the five bytes COPY TEXT gives meaning
# to — NUL, TAB, LF, CR and the escape character itself.
#
# ⚠ The table is created and dropped here, not added to the schema: this stage changes no schema.
# The CHECK is the design's own text. ⭐ An earlier revision wrote `block IS NULL OR
# octet_length(…) = …`; that disjunct is a semantic NO-OP, because a CHECK evaluating to NULL is
# satisfied in SQL, so it was describing a difference that does not exist. The one real divergence
# from the design is that `block` here is NULLABLE where `genome_genes.gene_spans` is NOT NULL —
# deliberately, because the empty-versus-NULL distinction below is a COPY coercion worth proving and
# the real column cannot hold either case.
BYTEA_PROBE_TABLE = "copy_bytea_probe"
BYTEA_PROBE_DDL = f"""
CREATE TABLE {BYTEA_PROBE_TABLE} (
    probe_id   int PRIMARY KEY,
    gene_count int  NOT NULL,
    block      bytea,
    CONSTRAINT block_length_is_gene_count_times_stride
        CHECK (octet_length(block) = gene_count * 16)
)
"""
BYTEA_PROBE_COLUMNS = ("probe_id", "gene_count", "block")


@pytest.fixture()
def bytea_probe(engine, session):
    """A real `bytea` column with the design's CHECK, created and dropped around one test.

    ⛔ **`session` is a parameter for a reason, not for convenience.** Each test reads the table
    back inside the session's transaction, which then holds an `ACCESS SHARE` lock on it. Fixtures
    tear down in reverse order of setup, so without rolling that transaction back FIRST the `DROP
    TABLE` below blocks on the test's own still-open read and the suite hangs rather than fails —
    measured, 2026-10-26, as `idle in transaction` against a waiting `DROP`.
    """
    from sqlalchemy import Column, Integer, LargeBinary, MetaData, Table, text

    with engine.begin() as connection:
        connection.execute(text(f"DROP TABLE IF EXISTS {BYTEA_PROBE_TABLE}"))
        connection.execute(text(BYTEA_PROBE_DDL))
    table = Table(
        BYTEA_PROBE_TABLE,
        MetaData(),
        Column("probe_id", Integer, primary_key=True),
        Column("gene_count", Integer),
        Column("block", LargeBinary),
    )
    yield table
    session.rollback()
    with engine.begin() as connection:
        connection.execute(text(f"DROP TABLE IF EXISTS {BYTEA_PROBE_TABLE}"))


def _store(session, table, rows):
    copy_rows(session, table, BYTEA_PROBE_COLUMNS, rows)
    session.commit()
    return dict(session.execute(select(table.c.probe_id, table.c.block)).all())


def test_every_byte_copy_text_gives_meaning_to_survives_a_round_trip(session, bytea_probe):
    """⛔ NUL, TAB, LF, CR and backslash — the five characters COPY TEXT escapes.

    Each is placed alone, leading, trailing and doubled, because an escaping bug that drops a
    trailing byte and one that drops a leading byte are different bugs.
    """
    dangerous = {0x00: "NUL", 0x09: "TAB", 0x0A: "LF", 0x0D: "CR", 0x5C: "backslash"}
    blobs = {}
    probe_id = 0
    for value, name in dangerous.items():
        byte = bytes([value])
        # ⚠ Every shape is exactly one 16-byte row, so the CHECK stays in force while the byte moves.
        for position, shape in (
            ("leading", byte + b"A" * 15),
            ("trailing", b"A" * 15 + byte),
            ("every byte", byte * 16),
            ("doubled mid-row", b"A" * 7 + byte * 2 + b"A" * 7),
            ("single mid-row", b"A" * 8 + byte + b"A" * 7),
        ):
            assert len(shape) == 16
            blobs[probe_id] = (shape, f"{name} {position}")
            probe_id += 1

    stored = _store(session, bytea_probe, [(key, 1, shape) for key, (shape, _) in blobs.items()])

    assert len(stored) == len(blobs) == 25
    for key, (shape, name) in blobs.items():
        assert stored[key] == shape, f"{name} ({shape!r}) came back as {stored[key]!r}"


def test_all_256_byte_values_round_trip_in_one_blob(session, bytea_probe):
    """⭐ The whole alphabet at once — not a sample of it, and in a known order."""
    blob = bytes(range(256))
    stored = _store(session, bytea_probe, [(1, 16, blob)])
    assert stored[1] == blob
    assert len(stored[1]) == 256


def test_postgres_octet_length_is_pythons_len_so_the_check_constraint_means_what_it_says(
    session, bytea_probe
):
    """⛔ The design's CHECK is `octet_length(block) = gene_count * 16`.

    If escaping inflated or deflated the stored value, that constraint would be measuring the
    escaped form and every block would be rejected — or worse, accepted at the wrong length.
    """
    from sqlalchemy import func

    blobs = {1: bytes(16), 2: b"\x5c" * 32, 3: b"\x00\x09\x0a\x0d\x5c" * 16}
    copy_rows(
        session,
        bytea_probe,
        BYTEA_PROBE_COLUMNS,
        [(key, len(blob) // 16, blob) for key, blob in blobs.items()],
    )
    session.commit()

    measured = dict(
        session.execute(select(bytea_probe.c.probe_id, func.octet_length(bytea_probe.c.block))).all()
    )
    assert measured == {key: len(blob) for key, blob in blobs.items()}


def test_a_real_gene_spans_block_round_trips_and_still_decodes(session, bytea_probe):
    """⭐ End to end: 5,000 genes encoded, COPYed, read back, decoded — the same coordinates.

    A genome-sized block, because a truncation that a 16-byte blob cannot show is exactly what a
    80,000-byte one does.
    """
    from bacatlas_backend.gene_blocks import GeneSpanTable, encode_gene_spans, pack_gene_flags

    gene_count = 5_000
    starts = [1 + 300 * j for j in range(gene_count)]
    ends = [start + 299 for start in starts]
    contigs = [j % 607 for j in range(gene_count)]
    flags = [
        pack_gene_flags(is_minus_strand=j % 2 == 0, is_five_prime_partial=j % 97 == 0, gff_phase=j % 3)
        for j in range(gene_count)
    ]
    block = encode_gene_spans(starts, ends, contigs, flags)
    assert len(block) == gene_count * 16

    stored = _store(session, bytea_probe, [(1, gene_count, block)])

    assert stored[1] == block
    table = GeneSpanTable(stored[1], gene_count=gene_count)
    assert [table.start_position(j) for j in range(gene_count)] == starts
    assert [table.end_position(j) for j in range(gene_count)] == ends
    assert [table.contig_index(j) for j in range(gene_count)] == contigs
    assert [table.flags(j) for j in range(gene_count)] == flags


def test_an_empty_block_is_not_null(session, bytea_probe):
    """⛔ Same rule as the empty array: `''::bytea` is *a genome with no genes*; NULL is *not loaded*."""
    stored = _store(session, bytea_probe, [(1, 0, b""), (2, 0, None)])
    assert stored[1] == b""
    assert stored[2] is None


def test_the_check_constraint_refuses_a_block_whose_length_disagrees(session, bytea_probe):
    """⭐ The invariant has to be enforceable in Postgres, not only in Python.

    ⚠ The exception is **psycopg's, not SQLAlchemy's**. `copy_rows` goes straight to the DBAPI
    cursor because `COPY … FROM STDIN` has no SQLAlchemy expression, so nothing wraps the error into
    `sqlalchemy.exc.IntegrityError`. Any caller that means to catch a constraint violation from this
    loader has to catch the psycopg class.
    """
    import psycopg

    with pytest.raises(psycopg.errors.CheckViolation, match="block_length_is_gene_count_times_stride"):
        copy_rows(session, bytea_probe, BYTEA_PROBE_COLUMNS, [(1, 3, bytes(2 * 16))])
    session.rollback()


def test_a_str_where_bytes_were_meant_is_refused_before_it_reaches_the_column(session, bytea_probe):
    """⛔ The one silent-wrong-value this path has, and the gate that now stops it.

    Measured before the gate existed: `copy_rows` accepted `"\\x" + "41" * 16` — 34 characters —
    and stored **16 bytes** of `A`. A legal block length, so the CHECK was satisfied and nothing
    raised. `_refuse_text_in_a_binary_column` turns that into a `TypeError` naming the column, on
    the first row, before anything is written.
    """
    hex_text = "\\x" + "41" * 16  # 34 characters that bytea's input function reads as 16 bytes

    with pytest.raises(TypeError, match="bytea column and was given a str"):
        copy_rows(session, bytea_probe, BYTEA_PROBE_COLUMNS, [(1, 1, hex_text)])
    session.rollback()

    assert session.execute(select(func.count()).select_from(bytea_probe)).scalar_one() == 0


def test_the_binary_type_gate_leaves_a_bytes_column_alone(session, bytea_probe):
    """⚠ The gate must not fire on the value it exists to protect, nor on NULL."""
    stored = _store(session, bytea_probe, [(1, 1, b"A" * 16), (2, 0, None)])
    assert stored == {1: b"A" * 16, 2: None}
