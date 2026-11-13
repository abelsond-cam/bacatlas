"""The genome-layer ingest, against real artifacts and a real database.

⛔ Two properties matter more than the counts, and neither is visible in a single successful run:
**a re-run must not double a row**, and **the cost per genome must not scale with its gene count**.
A loader that violates the first looks correct until someone re-runs it after a fix; one that
violates the second looks correct at 100 genomes and is unusable at 80,000.
"""

from __future__ import annotations

import os
import random
import shutil
import struct
from pathlib import Path

import pytest
from sqlalchemy import delete, func, select

from bacatlas_backend.gene_blocks import GENE_SPAN_STRIDE, SPAN_ENCODING_NSEQ_GENE_TABLE, read_nseq_gene_table
from bacatlas_backend.gff.gff_cds_parser import parse_genome_annotation
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.genome_gene_alignment import flatten_to_extractor_order
from bacatlas_backend.ingest.ingest_command_line import load_genome_layer, reconcile
from bacatlas_backend.ingest.ingest_genome_gene_table import (
    GenomeIngestError,
    gene_spans_from_alignment,
    ingest_one_genome,
)
from bacatlas_backend.ingest.ingest_pathogen_species import (
    SPECIES_KEY_BY_PARQUET_VALUE,
    ingest_pathogen_species,
    species_key_for_parquet_value,
)
from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation, GenomeNoncodingFeature
from bacatlas_backend.models.genome import Genome, GenomeAssembly, GenomeContig
from bacatlas_backend.models.genome_genes import GenomeGenes
from tests.conftest import PUBLISHED_SITE_CATALOGUE_DIR

pandas = pytest.importorskip("pandas")

#: The 200 committed `.nseq`/`.loci` — nuna's own encoding of 200 of the 280 probe GFFs.
SEQUENCE_DIR = PUBLISHED_SITE_CATALOGUE_DIR / "seq"

DATA_ROOT = Path(os.environ.get("BACATLAS_NUNA_DATA_ROOT", "~/developer/nuna/data")).expanduser()
MODEL_LABEL = "ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL"
RUN_ID = "ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0"

pytestmark = pytest.mark.skipif(
    not (DATA_ROOT / "raw" / "gff").is_dir(), reason=f"artifacts not mirrored at {DATA_ROOT}"
)


@pytest.fixture(scope="module")
def artifacts():
    return CatalogueArtifacts(
        data_root=DATA_ROOT, set_key="ecoli", model_label=MODEL_LABEL, run_id=RUN_ID
    )


@pytest.fixture(scope="module")
def one_genome(artifacts):
    """A real genome that exists in the store: `(sample_id, bakrep_dataset_id)`."""
    path = sorted(artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz"))[0]
    return path.parent.name, path.parent.parent.name


@pytest.fixture(scope="module")
def one_sequenced_genome(artifacts):
    """A real genome with BOTH a GFF in the mirror and a committed `.nseq`: `(sample_id, bakrep_dataset_id)`."""
    for path in sorted(artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz")):
        if (SEQUENCE_DIR / f"{path.parent.name}.nseq").is_file():
            return path.parent.name, path.parent.parent.name
    pytest.skip(f"no probe genome has both a GFF under {artifacts.annotation_root} and a .nseq under {SEQUENCE_DIR}")


def _wipe_genome_layer(session):
    for table in (
        GenomeGenes.__table__, GeneFunctionalAnnotation.__table__, GenomeNoncodingFeature.__table__,
        Gene.__table__, GenomeContig.__table__, GenomeAssembly.__table__, Genome.__table__,
    ):
        session.execute(delete(table))
    session.commit()


@pytest.fixture()
def clean_slate(session):
    """Ingest tests write real rows, so each starts and ends with an empty genome layer."""
    _wipe_genome_layer(session)
    yield
    _wipe_genome_layer(session)


def _species_key_of(artifacts, sample_id):
    """The probe `_meta` files carry `species`; the `.nseq`-only route must be TOLD it, so read it here."""
    frame = pandas.read_parquet(artifacts.per_genome(sample_id, "meta"), columns=["species"])
    value = frame["species"].iloc[0]
    return species_key_for_parquet_value(str(value))


def _load(session, artifacts, one_genome, *, write_noncoding=False, nseq_path=None):
    """⚠ `write_noncoding` defaults to False because the LOADER does — `genome_noncoding_feature` has
    no reader and became opt-in on 2026-10-29. A test that wants the table has to say so, which is
    the same discipline the flag itself now imposes. `nseq_path` likewise: None is the GFF-only route."""
    species_ids = ingest_pathogen_species(session)
    session.commit()
    sample_id, dataset_id = one_genome
    report = ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
        pathogen_species_ids=species_ids, write_noncoding_features=write_noncoding, nseq_path=nseq_path,
    )
    session.commit()
    return report


def test_the_counts_it_writes_are_the_counts_the_parquets_hold(session, artifacts, one_genome, clean_slate):
    """⛔ Reconciled against the SOURCE, not against the loader's own report."""
    sample_id, _ = one_genome
    report = _load(session, artifacts, one_genome, write_noncoding=True)

    meta_rows = len(pandas.read_parquet(artifacts.per_genome(sample_id, "meta"), columns=["flat_index"]))
    # ⚠ `start`, not `nc_index`: nuna narrowed the `_noncoding` schema to the six columns that have a
    # reader on 2026-10-29, and `nc_index` was one of the four dropped. The sidecars in the mirror
    # still carry it, so a test naming it would keep passing here and fail the day they are rebuilt.
    noncoding_rows = len(pandas.read_parquet(artifacts.per_genome(sample_id, "noncoding"), columns=["start"]))

    assert report.genes_written == meta_rows
    assert report.functional_annotations_written == meta_rows
    assert report.noncoding_features_written == noncoding_rows
    assert session.execute(select(func.count()).select_from(Gene)).scalar_one() == meta_rows


def test_the_noncoding_table_is_NOT_written_unless_it_is_asked_for(session, artifacts, one_genome, clean_slate):
    """⭐ The same run, without the opt-in, writes the gene layer and none of the 422 non-coding rows.

    ⛔ Zero here must be a DECISION, never a quiet failure of the load — so it is asserted beside a
    non-zero gene count from the same report, which is what tells the two apart.
    """
    report = _load(session, artifacts, one_genome)
    assert report.genes_written > 3_000, "the rest of the gene layer still loaded"
    assert report.noncoding_features_written == 0


def test_a_re_run_replaces_rather_than_doubles(session, artifacts, one_genome, clean_slate):
    """⭐ The property that only fails the SECOND time, which is when a fix is being re-applied."""
    first = _load(session, artifacts, one_genome)
    after_first = session.execute(select(func.count()).select_from(Gene)).scalar_one()

    second = _load(session, artifacts, one_genome)
    after_second = session.execute(select(func.count()).select_from(Gene)).scalar_one()

    assert after_first == after_second == first.genes_written
    assert second.genes_deleted == first.genes_written, "the re-run did not delete what it replaced"
    assert session.execute(select(func.count()).select_from(Genome)).scalar_one() == 1


def test_the_cost_does_not_scale_with_the_gene_count(session, engine, artifacts, one_genome, clean_slate):
    """⛔ A ratio against a logical minimum: three tables, delete + COPY each, plus the upserts.

    ⚠ The bound is a CONSTANT, not a multiple of anything — that is the whole assertion. A loader
    that issued one statement per gene would pass every correctness test in this file.
    """
    species_ids = ingest_pathogen_species(session)
    session.commit()
    sample_id, dataset_id = one_genome

    with SqlCostOracle(engine) as report:
        result = ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
            pathogen_species_ids=species_ids,
        )
        session.commit()

    assert result.genes_written > 3_000
    report.assert_at_most(16, what=f"loading {sample_id} ({result.genes_written:,} genes)")
    # ⚠ FOUR by default — contigs, genes, functional annotations, and the one-row `genome_genes`
    # block — five with `--with-noncoding-features`. The count moved when that table became opt-in
    # and again when the block arrived; what the test asserts is unchanged — a CONSTANT, not a
    # multiple of the gene count — and both polarities are pinned so neither can start issuing a
    # statement per gene.
    assert len([r for r in report.statements if r.category == "copy"]) == 4, report.statements


def test_the_symbol_column_holds_a_symbol_and_the_tag_column_holds_the_tag(
    session, artifacts, one_genome, clean_slate
):
    """⛔ `meta.gene_name` falls back to the genome-private locus tag on 16.6 % of rows.

    Loading it into `bakta_gene_symbol` puts a per-genome tag in a cross-genome symbol column, and a
    family holding only those would report one distinct symbol and read as unanimous.
    """
    import re

    _load(session, artifacts, one_genome)
    locus_tag_shape = re.compile(r"^[A-Z0-9]{4,}_\d+$")

    symbols = session.execute(
        select(Gene.bakta_gene_symbol).where(Gene.bakta_gene_symbol.isnot(None))
    ).scalars().all()
    tags = session.execute(select(Gene.locus_tag).where(Gene.locus_tag.isnot(None))).scalars().all()

    assert symbols, "no gene carries a Bakta symbol"
    offenders = [s for s in symbols if locus_tag_shape.fullmatch(s)]
    assert offenders == [], f"{len(offenders)} locus tags landed in the symbol column: {offenders[:5]}"
    assert len(tags) == session.execute(select(func.count()).select_from(Gene)).scalar_one(), (
        "the GFF carries a locus_tag on every CDS, so every gene must have one"
    )
    assert all(locus_tag_shape.fullmatch(tag) for tag in tags[:200])


def test_a_multi_term_ec_set_survives_as_a_list(session, artifacts, one_genome, clean_slate):
    """The 12.3 % case that would have raised against the original `String(32)`."""
    _load(session, artifacts, one_genome)
    rows = session.execute(
        select(GeneFunctionalAnnotation.ec_numbers).where(GeneFunctionalAnnotation.ec_numbers.isnot(None))
    ).scalars().all()

    assert rows, "no gene carries an EC number"
    assert all(isinstance(value, list) for value in rows)
    multi = [value for value in rows if len(value) > 1]
    assert multi, "no multi-term EC set in this genome — the case the array exists for is untested"
    assert max(len(",".join(value)) for value in rows) > 0


def test_not_annotated_stays_null_and_never_becomes_an_empty_list(session, artifacts, one_genome, clean_slate):
    """⛔ *Not annotated* and *annotated with nothing* are different claims."""
    _load(session, artifacts, one_genome)
    empty = session.execute(
        select(func.count()).select_from(GeneFunctionalAnnotation).where(
            (GeneFunctionalAnnotation.ec_numbers == []) | (GeneFunctionalAnnotation.gene_ontology_terms == [])
        )
    ).scalar_one()
    nulls = session.execute(
        select(func.count()).select_from(GeneFunctionalAnnotation).where(
            GeneFunctionalAnnotation.ec_numbers.is_(None)
        )
    ).scalar_one()

    assert empty == 0, f"{empty} rows hold an empty array, which the parquets never contain"
    assert nulls > 0, "no gene lacks an EC number — the NULL path is untested here"


def test_the_genome_row_records_what_was_measured(session, artifacts, one_genome, clean_slate):
    """⚠ `total_base_count` is EVERY contig; `contig_count_with_coding_sequence` is the CDS-bearing ones."""
    sample_id, dataset_id = one_genome
    _load(session, artifacts, one_genome)

    genome = session.execute(select(Genome).where(Genome.sample_id == sample_id)).scalar_one()
    contigs = session.execute(select(func.count()).select_from(GenomeContig)).scalar_one()

    assert genome.sample_id_kind is SampleIdentifierKind.BIOSAMPLE
    assert genome.bakrep_dataset_id == dataset_id
    assert genome.strand_is_observed is True
    assert genome.contig_count_with_coding_sequence == contigs
    assert genome.total_base_count > sum(
        row for row in session.execute(select(GenomeContig.length_bases)).scalars()
    ), "total_base_count is not larger than the CDS-bearing contigs — it should cover ALL of them"

    assembly = session.execute(select(GenomeAssembly)).scalar_one()
    assert assembly.annotation_file_path.endswith(f"{sample_id}.bakta.gff3.gz")
    # ⚠ NULL because the assembly manifest is not mirrored here — a fact about this machine.
    assert assembly.assembly_file_path is None


def test_the_species_vocabulary_is_translated_and_never_matched_directly():
    """⛔ `kpneumoniae` is not `kp`. A direct join loses every Klebsiella genome and looks fine."""
    assert species_key_for_parquet_value("kpneumoniae") == "kp"
    assert species_key_for_parquet_value("ecoli") == "ecoli"
    assert "kp" not in SPECIES_KEY_BY_PARQUET_VALUE, "the parquet vocabulary must not contain the browser key"
    with pytest.raises(KeyError, match="has no browser key"):
        species_key_for_parquet_value("saureus")


def test_a_genome_with_no_meta_parquet_is_refused_by_name(session, artifacts, clean_slate):
    """A refusal names the genome and the reason; it never writes a partial row."""
    species_ids = ingest_pathogen_species(session)
    session.commit()
    with pytest.raises(GenomeIngestError, match="no meta parquet"):
        ingest_one_genome(
            session, artifacts, sample_id="SAMEA_NOT_A_GENOME", bakrep_dataset_id="EA10",
            pathogen_species_ids=species_ids,
        )
    assert session.execute(select(func.count()).select_from(Genome)).scalar_one() == 0


def test_the_symbol_fallback_recovers_what_the_product_parquet_dropped(session, artifacts, clean_slate):
    """⭐ `product.gene` alone silently discards real symbols, and the fallback is what recovers them.

    `join_products` drops every `(contig_idx, start, end)` key that is not unique in the genome, so
    an ambiguous coordinate yields NaN in the product parquet while `meta.gene_name` still carries
    the symbol. Measured over 60 probe genomes: 25 such symbols against 249,797 both carry.

    ⚠ The fallback admits a value ONLY where it is not locus-tag-shaped, so it can never trade one
    silent error for the other.
    """
    import re

    # SAMEA104305638 is one of the genomes measured to carry the case (perC, yfdN).
    target = "SAMEA104305638"
    path = next(artifacts.annotation_root.glob(f"*/{target}/*.bakta.gff3.gz"), None)
    if path is None:
        pytest.skip(f"{target} is not in this store")

    products = pandas.read_parquet(artifacts.per_genome(target, "product"), columns=["flat_index", "gene"])
    meta = pandas.read_parquet(artifacts.per_genome(target, "meta"), columns=["flat_index", "gene_name"])
    joined = meta.merge(products, on="flat_index", how="left")
    tag_shape = re.compile(r"^[A-Z0-9]{4,}_\d+$")
    recoverable = joined[
        joined["gene"].isna()
        & joined["gene_name"].notna()
        & ~joined["gene_name"].fillna("").map(lambda value: bool(tag_shape.fullmatch(value)))
    ]
    assert len(recoverable) > 0, f"{target} no longer exhibits the case this test is about"

    species_ids = ingest_pathogen_species(session)
    session.commit()
    ingest_one_genome(
        session, artifacts, sample_id=target, bakrep_dataset_id=path.parent.parent.name,
        pathogen_species_ids=species_ids,
    )
    session.commit()

    for row in recoverable.itertuples(index=False):
        stored = session.execute(
            select(Gene.bakta_gene_symbol).where(Gene.flat_index == int(row.flat_index))
        ).scalar_one()
        assert stored == row.gene_name, f"flat_index {row.flat_index}: {stored!r} != {row.gene_name!r}"

    # ...and no locus tag came in with them.
    symbols = session.execute(
        select(Gene.bakta_gene_symbol).where(Gene.bakta_gene_symbol.isnot(None))
    ).scalars().all()
    assert [s for s in symbols if tag_shape.fullmatch(s)] == []


def test_cog_categories_are_split_into_the_set_they_always_were(session, artifacts, one_genome, clean_slate):
    """⛔ `CP` is two categories. Stored as a scalar it FITS, which is why the error was silent."""
    _load(session, artifacts, one_genome)
    rows = session.execute(
        select(GeneFunctionalAnnotation.cog_categories).where(
            GeneFunctionalAnnotation.cog_categories.isnot(None)
        )
    ).scalars().all()

    assert rows, "no gene carries a COG category"
    assert all(isinstance(value, list) for value in rows)
    assert all(len(letter) == 1 for value in rows for letter in value), "a value is not one letter"
    multi = [value for value in rows if len(value) > 1]
    assert multi, "no multi-category gene in this genome — the case the array exists for is untested"
    # The query that was silently wrong before: a single-letter match must now find them.
    from sqlalchemy import func as sql_func

    found = session.execute(
        select(sql_func.count()).select_from(GeneFunctionalAnnotation).where(
            GeneFunctionalAnnotation.cog_categories.any("C")
        )
    ).scalar_one()
    scalar_style = sum(1 for value in rows if value == ["C"])
    assert found > scalar_style, "the set query finds no more than an exact-scalar match would"


# ── the gene_spans block, and the three routes ────────────────────────────────────────────────
def _stored_block(session, sample_id):
    genome_id = session.execute(select(Genome.genome_id).where(Genome.sample_id == sample_id)).scalar_one()
    return session.execute(select(GenomeGenes).where(GenomeGenes.genome_id == genome_id)).scalar_one()


def _count(session, model):
    return session.execute(select(func.count()).select_from(model)).scalar_one()


def test_the_dual_route_copies_nunas_block_and_says_it_checked_it(
    session, artifacts, one_sequenced_genome, clean_slate
):
    """GFF and `.nseq` both held: today's load, plus the block — nuna's bytes, checked against the GFF."""
    sample_id, _ = one_sequenced_genome
    nseq = SEQUENCE_DIR / f"{sample_id}.nseq"
    report = _load(session, artifacts, one_sequenced_genome, nseq_path=nseq)
    published = read_nseq_gene_table(nseq, expect_sample=sample_id)

    assert (report.gene_spans_written, report.gene_spans_source) == (1, "nseq+encoded_checked")
    assert report.genes_written == published.gene_count > 3_000, "the rest of the gene layer still loaded"
    row = _stored_block(session, sample_id)
    assert row.gene_count == published.gene_count
    assert row.span_encoding == SPAN_ENCODING_NSEQ_GENE_TABLE
    assert len(row.gene_spans) == row.gene_count * GENE_SPAN_STRIDE
    assert bytes(row.gene_spans) == published.gene_table


def test_the_gff_only_route_encodes_the_block_and_it_is_still_nunas_bytes(
    session, artifacts, one_sequenced_genome, clean_slate
):
    """No `.nseq` named: the block is encoded from the GFF — and is STILL byte-identical to nuna's file."""
    sample_id, _ = one_sequenced_genome
    report = _load(session, artifacts, one_sequenced_genome)

    assert (report.gene_spans_written, report.gene_spans_source) == (1, "encoded")
    published = read_nseq_gene_table(SEQUENCE_DIR / f"{sample_id}.nseq", expect_sample=sample_id)
    row = _stored_block(session, sample_id)
    assert (row.gene_count, bytes(row.gene_spans)) == (published.gene_count, published.gene_table)


def _corrupted_copy(tmp_path: Path, sample_id: str, *, flat_index: int = 7, word: int = 0) -> Path:
    """A byte-exact copy of the committed `.nseq` with ONE gene-table byte changed — gene 7's `start`."""
    raw = bytearray((SEQUENCE_DIR / f"{sample_id}.nseq").read_bytes())
    (header_length,) = struct.unpack_from("<I", raw, 8)
    offset = 12 + header_length + flat_index * GENE_SPAN_STRIDE + word * 4
    raw[offset] = (raw[offset] + 1) % 256
    target = tmp_path / f"{sample_id}.nseq"
    target.write_bytes(raw)
    return target


def test_a_corrupted_nseq_gene_table_is_REFUSED_by_name_on_the_dual_route(
    session, artifacts, one_sequenced_genome, clean_slate, tmp_path
):
    """⛔ One byte off in one gene: the genome is refused, the refusal names the gene, nothing is written."""
    sample_id, dataset_id = one_sequenced_genome
    corrupted = _corrupted_copy(tmp_path, sample_id)
    species_ids = ingest_pathogen_species(session)
    session.commit()

    with pytest.raises(GenomeIngestError, match=rf"{sample_id}: the \.nseq at .* does not match the GFF") as refusal:
        ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
            pathogen_species_ids=species_ids, nseq_path=corrupted,
        )
    session.rollback()
    assert "gene table: first difference at flat_index 7" in str(refusal.value), str(refusal.value)
    assert _count(session, Genome) == 0 and _count(session, GenomeGenes) == 0


def test_a_corrupted_nseq_gene_table_is_REFUSED_by_name_on_the_nseq_only_route(
    session, artifacts, one_sequenced_genome, clean_slate, tmp_path
):
    """⛔ With no GFF the parquet is the witness: the same byte is caught against `meta.start`."""
    sample_id, _ = one_sequenced_genome
    corrupted = _corrupted_copy(tmp_path, sample_id)
    species_ids = ingest_pathogen_species(session)
    session.commit()

    disagreement = rf"{sample_id}: .* disagrees with the meta parquet on `start`"
    with pytest.raises(GenomeIngestError, match=disagreement) as refusal:
        ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
            pathogen_species_ids=species_ids, species_key=_species_key_of(artifacts, sample_id),
            nseq_path=corrupted,
        )
    session.rollback()
    assert "first at flat_index 7" in str(refusal.value), str(refusal.value)
    assert _count(session, Genome) == 0 and _count(session, GenomeGenes) == 0


def test_an_nseq_published_under_another_genomes_name_is_refused(
    session, artifacts, one_sequenced_genome, clean_slate, tmp_path
):
    """⭐ The header's `sample` key does its job at ingest: a mis-named file decodes perfectly and is refused anyway."""
    sample_id, dataset_id = one_sequenced_genome
    other = next(path for path in sorted(SEQUENCE_DIR.glob("*.nseq")) if path.stem != sample_id)
    impostor = tmp_path / f"{sample_id}.nseq"
    shutil.copyfile(other, impostor)
    species_ids = ingest_pathogen_species(session)
    session.commit()

    with pytest.raises(GenomeIngestError, match=rf"{sample_id}: the \.nseq at .* could not be read .* was asked for"):
        ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
            pathogen_species_ids=species_ids, nseq_path=impostor,
        )
    session.rollback()
    assert _count(session, Genome) == 0


def test_the_nseq_only_route_writes_the_genome_layer_and_NO_per_gene_rows(
    session, artifacts, one_sequenced_genome, clean_slate
):
    """The 71k shape: no GFF, a `.nseq`, a named species. Genome, contigs and the block — nothing per gene."""
    sample_id, _ = one_sequenced_genome
    nseq = SEQUENCE_DIR / f"{sample_id}.nseq"
    species_ids = ingest_pathogen_species(session)
    session.commit()
    report = ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
        pathogen_species_ids=species_ids, species_key=_species_key_of(artifacts, sample_id), nseq_path=nseq,
    )
    session.commit()
    published = read_nseq_gene_table(nseq, expect_sample=sample_id)

    assert (report.gene_spans_written, report.gene_spans_source) == (1, "nseq")
    assert (report.genes_written, report.functional_annotations_written, report.noncoding_features_written) == (0, 0, 0)
    assert report.contigs_written == len(published.contig_names) > 100
    assert _count(session, Gene) == 0 and _count(session, GeneFunctionalAnnotation) == 0
    assert _count(session, GenomeNoncodingFeature) == 0

    genome = session.execute(select(Genome).where(Genome.sample_id == sample_id)).scalar_one()
    assert genome.bakrep_dataset_id is None
    assert genome.strand_is_observed is True
    assert genome.coding_gene_count == published.gene_count
    assert genome.contig_count_with_coding_sequence == len(published.contig_names)
    assert genome.total_base_count == published.base_count
    assert session.execute(select(GenomeAssembly)).scalar_one().annotation_file_path is None
    row = _stored_block(session, sample_id)
    assert (row.gene_count, bytes(row.gene_spans)) == (published.gene_count, published.gene_table)


def test_the_nseq_only_route_REQUIRES_a_species_key(session, artifacts, one_sequenced_genome, clean_slate):
    """The cohort this route exists for has no `species` column to read, so the key is a precondition."""
    sample_id, _ = one_sequenced_genome
    species_ids = ingest_pathogen_species(session)
    session.commit()
    with pytest.raises(GenomeIngestError, match=f"{sample_id}: .*needs `species_key`"):
        ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
            pathogen_species_ids=species_ids, nseq_path=SEQUENCE_DIR / f"{sample_id}.nseq",
        )
    assert _count(session, Genome) == 0


def test_a_genome_with_neither_a_gff_nor_an_nseq_is_refused_by_name(
    session, artifacts, one_sequenced_genome, clean_slate
):
    sample_id, _ = one_sequenced_genome
    species_ids = ingest_pathogen_species(session)
    session.commit()
    with pytest.raises(GenomeIngestError, match=f"{sample_id}: in the roster but under no dataset directory"):
        ingest_one_genome(
            session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
            pathogen_species_ids=species_ids, species_key=_species_key_of(artifacts, sample_id), nseq_path=None,
        )
    assert _count(session, Genome) == 0


def test_a_re_run_that_changes_route_replaces_the_per_gene_rows_it_no_longer_writes(
    session, artifacts, one_sequenced_genome, clean_slate
):
    """⭐ GFF route, then `.nseq`-only for the same genome: the first load's gene rows must be gone."""
    sample_id, _ = one_sequenced_genome
    nseq = SEQUENCE_DIR / f"{sample_id}.nseq"
    first = _load(session, artifacts, one_sequenced_genome, nseq_path=nseq)
    assert first.genes_written > 3_000

    species_ids = ingest_pathogen_species(session)
    second = ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
        pathogen_species_ids=species_ids, species_key=_species_key_of(artifacts, sample_id), nseq_path=nseq,
    )
    session.commit()

    assert second.genes_deleted == first.genes_written, "the route change did not delete what it replaced"
    assert _count(session, Gene) == 0 and _count(session, GeneFunctionalAnnotation) == 0
    assert _count(session, Genome) == 1 and _count(session, GenomeGenes) == 1


def test_the_nseq_only_route_writes_the_same_genome_and_contig_rows_as_the_gff_route(
    session, artifacts, one_sequenced_genome, clean_slate
):
    """⛔ Field by field, on one real genome loaded both ways. Three fields differ BY DESIGN and are named.

    `bakrep_dataset_id` and `annotation_file_path` are facts about the GFF, which the second route never
    opened. ⚠ `total_base_count` is the one that is NOT merely absent: the GFF route stores every contig
    in the `##FASTA` block, the `.nseq` route stores `n_bases`, which is the CDS-bearing contigs only —
    the only contigs a `.nseq` holds. The difference is measured here, not assumed away.
    """
    sample_id, dataset_id = one_sequenced_genome
    nseq = SEQUENCE_DIR / f"{sample_id}.nseq"
    species_key = _species_key_of(artifacts, sample_id)
    species_ids = ingest_pathogen_species(session)
    session.commit()

    def snapshot():
        genome = session.execute(select(Genome).where(Genome.sample_id == sample_id)).scalar_one()
        contigs = session.execute(
            select(GenomeContig.seqid, GenomeContig.contig_index, GenomeContig.contig_name, GenomeContig.length_bases)
            .where(GenomeContig.genome_id == genome.genome_id)
            .order_by(GenomeContig.contig_index)
        ).all()
        block = _stored_block(session, sample_id)
        assembly = session.execute(
            select(GenomeAssembly).where(GenomeAssembly.genome_id == genome.genome_id)
        ).scalar_one()
        fields = {
            column.name: getattr(genome, column.name)
            for column in Genome.__table__.columns
            if column.name != "genome_id"
        }
        return (
            fields,
            [tuple(row) for row in contigs],
            (block.gene_count, block.span_encoding, bytes(block.gene_spans)),
            assembly.annotation_file_path,
        )

    ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
        pathogen_species_ids=species_ids, species_key=species_key, nseq_path=nseq,
    )
    session.commit()
    via_gff = snapshot()
    _wipe_genome_layer(session)
    ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=None,
        pathogen_species_ids=species_ids, species_key=species_key, nseq_path=nseq,
    )
    session.commit()
    via_nseq = snapshot()

    # The contig rows and the block: identical, every field, every row.
    assert len(via_gff[1]) > 100
    assert via_nseq[1] == via_gff[1]
    assert via_nseq[2] == via_gff[2]

    gff_fields, nseq_fields = via_gff[0], via_nseq[0]
    assert gff_fields["bakrep_dataset_id"] == dataset_id and nseq_fields["bakrep_dataset_id"] is None
    assert via_gff[3].endswith(f"{sample_id}.bakta.gff3.gz") and via_nseq[3] is None

    alignment = flatten_to_extractor_order(parse_genome_annotation(artifacts.genome_annotation(sample_id, dataset_id)))
    cds_contig_bases = sum(alignment.contig_lengths_by_seqid[seqid] for seqid in alignment.contig_names_in_index_order)
    assert nseq_fields["total_base_count"] == cds_contig_bases == sum(length for *_, length in via_nseq[1])
    assert gff_fields["total_base_count"] == alignment.total_base_count >= cds_contig_bases
    cds_free_contigs = len(alignment.contig_lengths_by_seqid) - len(alignment.contig_names_in_index_order)
    assert (gff_fields["total_base_count"] > cds_contig_bases) == (cds_free_contigs > 0), (
        f"{cds_free_contigs} CDS-free contigs, {gff_fields['total_base_count'] - cds_contig_bases} bases apart"
    )

    differing = {name for name in gff_fields if gff_fields[name] != nseq_fields[name]}
    assert differing <= {"bakrep_dataset_id", "total_base_count"}, differing
    assert nseq_fields["coding_gene_count"] == gff_fields["coding_gene_count"] > 3_000
    assert nseq_fields["contig_count_with_coding_sequence"] == gff_fields["contig_count_with_coding_sequence"]


# ── the GFF → block differential, on a sample of the store ────────────────────────────────────
def test_the_block_encoded_from_the_gff_is_the_block_nuna_published(artifacts):
    """⭐ `gene_spans_from_alignment` over 20 seeded probe GFFs against the committed `.nseq`, byte for byte.

    ⛔ Coverage first: how many genomes and genes were compared is asserted before any difference is.
    """
    index = {path.parent.name: path for path in artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz")}
    both = sorted(sample for sample in index if (SEQUENCE_DIR / f"{sample}.nseq").is_file())
    if len(both) < 20:
        pytest.skip(f"only {len(both)} probe genomes have both a GFF in the mirror and a committed .nseq")
    chosen = random.Random(2026).sample(both, 20)

    genomes = genes = 0
    differences: list[str] = []
    for sample in chosen:
        alignment = flatten_to_extractor_order(parse_genome_annotation(index[sample]))
        encoded = gene_spans_from_alignment(alignment)
        published = read_nseq_gene_table(SEQUENCE_DIR / f"{sample}.nseq", expect_sample=sample)
        genomes += 1
        genes += published.gene_count
        if encoded != published.gene_table:
            differences.append(f"{sample}: gene table")
        if tuple(alignment.contig_names_in_index_order) != published.contig_names:
            differences.append(f"{sample}: contig names")
        lengths = tuple(alignment.contig_lengths_by_seqid[s] for s in alignment.contig_names_in_index_order)
        if lengths != published.contig_lengths:
            differences.append(f"{sample}: contig lengths")
    print(f"\nGFF → gene_spans vs the committed .nseq: genomes {genomes} · genes {genes:,} · differences {differences}")
    assert genomes == 20 and genes > 60_000, (genomes, genes)
    assert differences == []


# ── reconcile, through the loader ─────────────────────────────────────────────────────────────
def test_reconcile_counts_the_gene_spans_blocks_beside_the_rest(session, artifacts, one_sequenced_genome, clean_slate):
    """⭐ Through `load_genome_layer` with the tree named: the dual route, the frame line, the render, the count."""
    sample_id, _ = one_sequenced_genome
    located = CatalogueArtifacts(
        data_root=DATA_ROOT, set_key="ecoli", model_label=MODEL_LABEL, run_id=RUN_ID, nseq_root=SEQUENCE_DIR
    )
    report = load_genome_layer(session, located, only={sample_id}, progress_every=0)

    assert (report.genomes_loaded, report.genomes_refused, report.gene_spans) == (1, 0, 1), report.refusals
    assert report.gene_spans_routes == {"nseq+encoded_checked": 1, "encoded": 0, "nseq": 0}
    frame = "280 with a GFF, 200 with a .nseq, 200 with both, 0 with neither"
    assert any(frame in note for note in report.notes), report.notes
    assert "gene_spans blocks      1  (nseq+encoded_checked 1 · encoded 0 · nseq 0)" in report.render()
    assert reconcile(session, report) == []

    session.execute(delete(GenomeGenes.__table__))
    session.commit()
    assert reconcile(session, report) == [
        "gene_spans: wrote 1 but the table holds 0 for the 1 genome(s) this run loaded"
    ]
