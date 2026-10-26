"""The Sequence tab's endpoint, against real GFFs and the real database.

⛔ **The defect this file exists for is a flank taken from the wrong end.** On a minus-strand gene —
about half of them — the upstream flank sits at HIGHER contig coordinates and is reverse-complemented
with the gene. Slicing `start - flank` unconditionally returns the *downstream* flank instead, and it
renders as a perfectly plausible 100 bases of DNA. Nothing on the page could contradict it, so the
check below goes back to the contig and reads the bases itself.

⚠ These run against `BACATLAS_DATABASE_URL` with the catalogues loaded AND `BACATLAS_ROOT_GFF`
pointing at the annotation store. Missing either, they skip with the reason.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import TEXT, func, select
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.gff.gene_sequence_reader import reverse_complement
from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle
from bacatlas_backend.models.gene import Gene, GeneLocusMembership
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.sequence_source.contig_sequence_provider import (
    ContigSequenceProvider,
    GenomeSequenceRequest,
)
from bacatlas_backend.services.gene_sequence_service import SequenceUnavailable, load_gene_sequences

PANGENOME = 1


@pytest.fixture(scope="module")
def application():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set")
    root = os.environ.get("BACATLAS_ROOT_GFF")
    if not root:
        pytest.skip("BACATLAS_ROOT_GFF is not set — the sequence endpoint reads the original GFFs")
    return create_application(Configuration.from_environment())


@pytest.fixture(scope="module")
def client(application):
    return application.test_client()


@pytest.fixture(scope="module")
def engine(application):
    return application.extensions["bacatlas_database"].engine


@pytest.fixture(scope="module")
def gff_root(application):
    return application.config["BACATLAS"].artifact_roots["gff"]


@pytest.fixture()
def provider(gff_root):
    """The HELD-only chain — no ATB store, so these tests cannot reach the network.

    ⛔ Function-scoped on purpose: a module-scoped provider would carry its parsed-genome cache
    between tests, and the "one statement and one file" assertion below measures a cold read.
    """
    return ContigSequenceProvider(gff_root=gff_root, atb_store=None)


def _a_gene(session, *, strand: str):
    """One (sample_id, node_label) whose gene sits on `strand`, away from either contig end."""
    row = session.execute(
        select(Genome.sample_id, Locus.node_label)
        .select_from(Gene)
        .join(
            GeneLocusMembership,
            (GeneLocusMembership.genome_id == Gene.genome_id)
            & (GeneLocusMembership.flat_index == Gene.flat_index),
        )
        .join(Genome, Genome.genome_id == Gene.genome_id)
        .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
        .where(
            GeneLocusMembership.pangenome_id == PANGENOME,
            Gene.strand == strand,
            Gene.start_position > 500,
            Gene.length_nt > 900,
        )
        .order_by(Gene.genome_id, Gene.flat_index)
        .limit(1)
    ).one()
    return row.sample_id, row.node_label


# ── the flank orientation, which is the whole point ────────────────────────────────────────────
def test_a_MINUS_strand_gene_takes_its_upstream_flank_from_the_HIGHER_coordinates(engine, provider):
    """⛔⛔ **The bug of record, checked against the contig itself.**

    `load_meta_flanks` orients on strand and `app.js:4376-4386` states the convention. A reader
    following a guide designed against it would otherwise be pointed at the opposite end of the gene.
    """
    with Session(engine) as session:
        sample_id, label = _a_gene(session, strand="-")
        rows = load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    assert rows, f"{sample_id} has no gene at {label}"
    gene = rows[0]
    assert gene.strand == "-"

    contig = provider.contig_sequences(
        GenomeSequenceRequest(sample_id=sample_id, bakrep_dataset_id=_dataset_id(engine, sample_id))
    )[gene.seqid]
    higher = contig[gene.end_position : gene.end_position + 100]
    lower = contig[max(0, gene.start_position - 101) : gene.start_position - 1]

    assert gene.sequence.upstream_flank_sequence == reverse_complement(higher)
    assert gene.sequence.downstream_flank_sequence == reverse_complement(lower)
    # ⚠ And they are NOT each other: on a palindromic region this test would pass either way.
    assert gene.sequence.upstream_flank_sequence != gene.sequence.downstream_flank_sequence


def test_a_PLUS_strand_gene_takes_its_upstream_flank_from_the_LOWER_coordinates(engine, provider):
    with Session(engine) as session:
        sample_id, label = _a_gene(session, strand="+")
        rows = load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    gene = rows[0]
    contig = provider.contig_sequences(
        GenomeSequenceRequest(sample_id=sample_id, bakrep_dataset_id=_dataset_id(engine, sample_id))
    )[gene.seqid]
    assert gene.sequence.upstream_flank_sequence == contig[gene.start_position - 101 : gene.start_position - 1]
    assert gene.sequence.downstream_flank_sequence == contig[gene.end_position : gene.end_position + 100]


def _dataset_id(engine, sample_id: str) -> str:
    with Session(engine) as session:
        return session.execute(
            select(Genome.bakrep_dataset_id).where(Genome.sample_id == sample_id)
        ).scalar_one()


# ── the sliced values against the stored ones ──────────────────────────────────────────────────
def test_the_sliced_sequence_REPRODUCES_the_columns_ingest_wrote(engine, provider):
    """⭐ Two independent computations of the same thing, from the same coordinates.

    `protein_length_aa` is a column written at ingest; the sequence here is sliced live. They are
    computed by different code at different times, so agreeing is evidence — and a disagreement
    would mean the coordinate convention moved under one of them.

    ⚠ **This test used to compare GC the same way and no longer can.** `stored_gc_percent` was
    dropped 2026-10-26 with its column, and that comparison is what confirmed AllTheBacteria
    streaming was byte-correct (44.968553 from streamed bases, 44.968553 stored). What survives is
    stronger and needs no column: T6 compares the page's own `gcPct` over the page's own decoded
    bases against this path's computed value, so neither side reads the database.
    """
    checked = 0
    with Session(engine) as session:
        pairs = session.execute(
            select(Genome.sample_id, Locus.node_label)
            .select_from(Gene)
            .join(
                GeneLocusMembership,
                (GeneLocusMembership.genome_id == Gene.genome_id)
                & (GeneLocusMembership.flat_index == Gene.flat_index),
            )
            .join(Genome, Genome.genome_id == Gene.genome_id)
            .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
            .where(GeneLocusMembership.pangenome_id == PANGENOME, Gene.protein_length_aa.is_not(None))
            .order_by(Gene.genome_id, Gene.flat_index)
            .limit(40)
        ).all()
        for sample_id, label in pairs:
            for gene in load_gene_sequences(
                session,
                pangenome_id=PANGENOME,
                node_label=label,
                sample_id=sample_id,
                sequence_provider=provider,
            ):
                assert len(gene.sequence.protein_sequence) == gene.stored_protein_length_aa
                # ⚠ The span INCLUDES the stop codon, so the protein is `length/3 - 1` residues.
                assert len(gene.sequence.coding_sequence) == gene.length_nt
                checked += 1
    # ⛔ Coverage before the verdict: a query that matched nothing would report 40 green assertions.
    assert checked >= 40, f"only {checked} genes were examined"


# ── the three different answers ────────────────────────────────────────────────────────────────
def test_a_genome_with_NO_gene_here_is_an_ANSWER_and_not_a_failure(client, engine):
    """⛔ `genes: []` with a 200. A 404 would read as "that locus does not exist"."""
    with Session(engine) as session:
        label = session.execute(
            select(Locus.node_label)
            .where(Locus.pangenome_id == PANGENOME, Locus.member_genome_count == 1)
            .limit(1)
        ).scalar_one()
        absent = session.execute(
            select(Genome.sample_id)
            .where(
                Genome.genome_id.not_in(
                    select(GeneLocusMembership.genome_id)
                    .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
                    .where(Locus.node_label == label, Locus.pangenome_id == PANGENOME)
                )
            )
            .limit(1)
        ).scalar_one()
    response = client.get(f"/api/v1/species/ecoli/genomes/{absent}/loci/{label}/sequence")
    assert response.status_code == 200
    assert response.get_json()["genes"] == []


def test_an_unknown_genome_and_an_unknown_locus_are_NAMED_404s(client):
    missing_genome = client.get("/api/v1/species/ecoli/genomes/SAMEA000000/loci/2811/sequence")
    assert missing_genome.status_code == 404
    assert "SAMEA000000" in missing_genome.get_json()["detail"]

    missing_locus = client.get(
        "/api/v1/species/ecoli/genomes/SAMEA103923484/loci/not-a-locus/sequence"
    )
    assert missing_locus.status_code == 404
    assert "not-a-locus" in missing_locus.get_json()["detail"]


def test_an_unreadable_annotation_file_RAISES_rather_than_returning_no_genes(engine, tmp_path):
    """⛔ The distinction the whole endpoint turns on: *we could not read it* is not *there is none*."""
    with Session(engine) as session:
        sample_id, label = _a_gene(session, strand="+")
        with pytest.raises(SequenceUnavailable, match="not on this server"):
            load_gene_sequences(
                session,
                pangenome_id=PANGENOME,
                node_label=label,
                sample_id=sample_id,
                sequence_provider=ContigSequenceProvider(gff_root=tmp_path, atb_store=None),
            )


# ── the things that are easy to get plausibly wrong ────────────────────────────────────────────
def test_the_contig_is_named_by_its_NAME_and_not_by_its_index(client, engine):
    """⛔ `contig_index` enumerates contigs that HAVE a CDS, so `index + 1` is a different contig.

    ⚠ **And it coincides most of the time, which is exactly why this is dangerous.** Measured on the
    loaded catalogues: the derived name is wrong on **7,696 of 26,878 contigs (28.6 %)** and right on
    the rest — so a page printing `contig{index + 1}` looks correct on seven contigs in ten and names
    a real contig the gene is not on for the other three. The test therefore selects from the
    DIVERGING set rather than taking whatever comes first, which is how this assertion passed
    vacuously the first time it was written.
    """
    derived = func.concat("contig", func.lpad(func.cast(Gene.contig_index + 1, TEXT), 5, "0"))
    with Session(engine) as session:
        row = session.execute(
            select(Genome.sample_id, Locus.node_label, Gene.contig_index, GenomeContig.contig_name)
            .select_from(Gene)
            .join(
                GeneLocusMembership,
                (GeneLocusMembership.genome_id == Gene.genome_id)
                & (GeneLocusMembership.flat_index == Gene.flat_index),
            )
            .join(Genome, Genome.genome_id == Gene.genome_id)
            .join(
                GenomeContig,
                (GenomeContig.genome_id == Gene.genome_id)
                & (GenomeContig.contig_index == Gene.contig_index),
            )
            .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
            .where(
                GeneLocusMembership.pangenome_id == PANGENOME,
                GenomeContig.contig_name != derived,
            )
            .limit(1)
        ).first()
    assert row is not None, "no gene sits on a contig whose name differs from its index"
    payload = client.get(
        f"/api/v1/species/ecoli/genomes/{row.sample_id}/loci/{row.node_label}/sequence"
    ).get_json()
    gene = payload["genes"][0]
    assert gene["contig_name"] == row.contig_name
    assert gene["contig_name"] != f"contig{row.contig_index + 1:05d}"
    assert gene["seqid"].endswith(gene["contig_name"])


def test_rho_above_one_returns_EVERY_copy_with_its_ordinal(client, engine):
    """⛔ `COUNT(DISTINCT)` territory: one genome really can hold two genes at one locus."""
    with Session(engine) as session:
        found = session.execute(
            select(Genome.sample_id, Locus.node_label, func.count())
            .select_from(GeneLocusMembership)
            .join(Genome, Genome.genome_id == GeneLocusMembership.genome_id)
            .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
            .where(GeneLocusMembership.pangenome_id == PANGENOME)
            .group_by(Genome.sample_id, Locus.node_label)
            .having(func.count() > 1)
            .limit(1)
        ).first()
    if found is None:
        pytest.skip("no genome holds two genes at one locus in this catalogue")
    sample_id, label, copies = found
    genes = client.get(
        f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence"
    ).get_json()["genes"]
    assert len(genes) == copies
    assert [gene["copy_ordinal"] for gene in genes] == list(range(1, copies + 1))
    assert {gene["copy_count"] for gene in genes} == {copies}
    # ⚠ And they are different genes, not one row repeated.
    assert len({gene["flat_index"] for gene in genes}) == copies


def test_a_flank_that_RUNS_OFF_the_contig_says_so(client, engine):
    """⚠ Draft assemblies: mean contig 14.4 kb, ~14 genes, so this is common rather than exotic.

    A short flank that does not say it is short reads as a complete one.
    """
    with Session(engine) as session:
        row = session.execute(
            select(Genome.sample_id, Locus.node_label)
            .select_from(Gene)
            .join(
                GeneLocusMembership,
                (GeneLocusMembership.genome_id == Gene.genome_id)
                & (GeneLocusMembership.flat_index == Gene.flat_index),
            )
            .join(Genome, Genome.genome_id == Gene.genome_id)
            .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
            .where(GeneLocusMembership.pangenome_id == PANGENOME, Gene.start_position < 50)
            .limit(1)
        ).one()
    gene = client.get(
        f"/api/v1/species/ecoli/genomes/{row.sample_id}/loci/{row.node_label}/sequence"
    ).get_json()["genes"][0]
    truncated = (
        gene["upstream_flank_is_truncated_by_contig_end"]
        if gene["strand"] == "+"
        else gene["downstream_flank_is_truncated_by_contig_end"]
    )
    assert truncated is True
    assert len(gene["upstream_flank_sequence"]) < 100 or len(gene["downstream_flank_sequence"]) < 100


def test_the_sequence_view_is_ONE_statement_and_ONE_file(application, engine, provider):
    """⚠ However many copies the genome holds. The file is opened once and cached per process."""
    with Session(engine) as session:
        sample_id, label = _a_gene(session, strand="+")
    with Session(engine) as session, SqlCostOracle(engine) as report:
        load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    report.assert_at_most(1, what="one gene sequence view")


def test_each_flank_NAMES_the_contig_coordinates_it_came_from(client, engine):
    """⭐ So the page can say *"contig 1,388–1,487, reverse-complemented"* without re-deriving it.

    ⛔ A client working out which end a flank came from would have to re-implement the strand rule,
    which is the exact thing the reader exists to hold in one place — and getting it wrong there
    mislabels a correct sequence, which is worse than showing the wrong one.
    """
    with Session(engine) as session:
        sample_id, label = _a_gene(session, strand="-")
    gene = client.get(
        f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence"
    ).get_json()["genes"][0]

    upstream = gene["upstream_flank_span"]
    downstream = gene["downstream_flank_span"]
    # ⛔ On a MINUS gene the upstream flank is at the HIGHER coordinates.
    assert upstream[0] > gene["end_position"]
    assert downstream[1] < gene["start_position"]
    # The spans describe exactly the strings that were returned.
    assert upstream[1] - upstream[0] + 1 == len(gene["upstream_flank_sequence"])
    assert downstream[1] - downstream[0] + 1 == len(gene["downstream_flank_sequence"])


def test_an_empty_flank_has_a_NULL_span_rather_than_a_zero_length_one(engine):
    """⚠ `[1, 0]` is not a span; it is arithmetic leaking. A flank that does not exist says so."""
    from bacatlas_backend.gff.gene_sequence_reader import read_gene_sequence

    view = read_gene_sequence("ATGAAATAG", start_position=1, end_position=9, strand="+", flank_length=100)
    assert view.upstream_flank_span is None
    assert view.upstream_flank_sequence == ""
    assert view.upstream_flank_is_truncated_by_contig_end is True
