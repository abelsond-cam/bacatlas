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
        GenomeGenes.__table__, GenomeContig.__table__, GenomeAssembly.__table__, Genome.__table__,
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


def _load(session, artifacts, one_genome, *, nseq_path=None):
    """`nseq_path` None is the GFF-only route."""
    species_ids = ingest_pathogen_species(session)
    session.commit()
    sample_id, dataset_id = one_genome
    report = ingest_one_genome(
        session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
        pathogen_species_ids=species_ids, nseq_path=nseq_path,
    )
    session.commit()
    return report


def test_the_counts_it_writes_are_the_counts_the_parquets_hold(session, artifacts, one_genome, clean_slate):
    """⛔ Reconciled against the SOURCE, not against the loader's own report."""
    sample_id, _ = one_genome
    report = _load(session, artifacts, one_genome)

    meta_rows = len(pandas.read_parquet(artifacts.per_genome(sample_id, "meta"), columns=["flat_index"]))

    assert report.gene_count == meta_rows > 3_000
    assert _stored_block(session, sample_id).gene_count == meta_rows


def test_a_re_run_replaces_rather_than_doubles(session, artifacts, one_genome, clean_slate):
    """⭐ The property that only fails the SECOND time, which is when a fix is being re-applied."""
    first = _load(session, artifacts, one_genome)
    after_first = (_count(session, GenomeContig), _count(session, GenomeGenes))

    _load(session, artifacts, one_genome)
    after_second = (_count(session, GenomeContig), _count(session, GenomeGenes))

    assert after_first == after_second == (first.contigs_written, 1)
    assert _count(session, Genome) == 1


def test_the_cost_does_not_scale_with_the_gene_count(session, engine, artifacts, one_genome, clean_slate):
    """⛔ A ratio against a logical minimum: two tables, delete + COPY each, plus the upserts.

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

    assert result.gene_count > 3_000
    report.assert_at_most(16, what=f"loading {sample_id} ({result.gene_count:,} genes)")
    # ⚠ TWO — contigs and the one-row `genome_genes` block — since the per-gene tables went (D8).
    # What the test asserts is unchanged: a CONSTANT, not a multiple of the gene count.
    assert len([r for r in report.statements if r.category == "copy"]) == 2, report.statements


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
    assert report.gene_count == published.gene_count > 3_000
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


def test_the_nseq_only_route_writes_the_genome_layer(session, artifacts, one_sequenced_genome, clean_slate):
    """The 71k shape: no GFF, a `.nseq`, a named species. Genome, contigs and the block."""
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
    assert report.gene_count == published.gene_count
    assert report.contigs_written == len(published.contig_names) > 100

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
