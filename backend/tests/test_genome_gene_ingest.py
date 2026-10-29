"""The genome-layer ingest, against real artifacts and a real database.

⛔ Two properties matter more than the counts, and neither is visible in a single successful run:
**a re-run must not double a row**, and **the cost per genome must not scale with its gene count**.
A loader that violates the first looks correct until someone re-runs it after a fix; one that
violates the second looks correct at 100 genomes and is unusable at 80,000.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from sqlalchemy import delete, func, select

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_genome_gene_table import GenomeIngestError, ingest_one_genome
from bacatlas_backend.ingest.ingest_pathogen_species import (
    SPECIES_KEY_BY_PARQUET_VALUE,
    ingest_pathogen_species,
    species_key_for_parquet_value,
)
from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation, GenomeNoncodingFeature
from bacatlas_backend.models.genome import Genome, GenomeAssembly, GenomeContig

pandas = pytest.importorskip("pandas")

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


@pytest.fixture()
def clean_slate(session):
    """Ingest tests write real rows, so each starts and ends with an empty genome layer."""
    def _wipe():
        for table in (
            GeneFunctionalAnnotation.__table__, GenomeNoncodingFeature.__table__,
            Gene.__table__, GenomeContig.__table__, GenomeAssembly.__table__, Genome.__table__,
        ):
            session.execute(delete(table))
        session.commit()

    _wipe()
    yield
    _wipe()


def _load(session, artifacts, one_genome, *, write_noncoding=False):
    """⚠ `write_noncoding` defaults to False because the LOADER does — `genome_noncoding_feature` has
    no reader and became opt-in on 2026-10-29. A test that wants the table has to say so, which is
    the same discipline the flag itself now imposes."""
    species_ids = ingest_pathogen_species(session)
    session.commit()
    sample_id, dataset_id = one_genome
    report = ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
        pathogen_species_ids=species_ids, write_noncoding_features=write_noncoding,
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
    # ⚠ THREE by default, four with `--with-noncoding-features`. The count moved when that table
    # became opt-in; what the test asserts is unchanged — a CONSTANT, not a multiple of the gene
    # count — and both polarities are pinned so neither can start issuing a statement per gene.
    assert len([r for r in report.statements if r.category == "copy"]) == 3, report.statements


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
