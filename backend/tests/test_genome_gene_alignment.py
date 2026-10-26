"""The alignment gate, run against the real cohort — the one check a later test cannot make.

⛔ A misaligned `flat_index` puts one gene's sequence, product, UniRef50 family and locus membership
under another gene's name, and **nothing downstream can see it**: every row is well-formed, every
join succeeds, and the page renders. It is caught here or not at all.

⚠ **Coverage is reported, always.** By default this walks 30 genomes (~10 s). Set
`BACATLAS_FULL_COHORT=1` for all 280 (~85 s), which is what an ingest run must pass. Either way the
test asserts the number it actually examined, because a suite that examined 30 of 280 must say 30.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from bacatlas_backend.gff.gff_cds_parser import parse_genome_annotation
from bacatlas_backend.ingest.genome_gene_alignment import (
    GeneAlignmentError,
    check_against_meta,
    flatten_to_extractor_order,
)

DATA_ROOT = Path(os.environ.get("BACATLAS_NUNA_DATA_ROOT", "~/developer/nuna/data")).expanduser()
FULL_COHORT = os.environ.get("BACATLAS_FULL_COHORT", "").lower() in {"1", "true", "yes"}
SAMPLED_GENOMES = 30

pandas = pytest.importorskip("pandas")

pytestmark = pytest.mark.skipif(
    not (DATA_ROOT / "raw" / "gff").is_dir(), reason=f"the GFF store is not mirrored at {DATA_ROOT}"
)


def _genomes():
    paths = sorted((DATA_ROOT / "raw" / "gff").glob("*/*/*.bakta.gff3.gz"))
    return paths if FULL_COHORT else paths[:SAMPLED_GENOMES]


@pytest.fixture(scope="module")
def alignment_report():
    """Align every genome in scope and gather the counts, so each test asserts on one pass."""
    report = {"genomes": 0, "genes": 0, "checks": {"contig_index": 0, "coordinates": 0, "protein": 0}}
    failures = []
    for gff in _genomes():
        sample_id = gff.parent.name
        meta = DATA_ROOT / "proc" / "embeddings" / "meta" / f"{sample_id}_meta.parquet"
        if not meta.is_file():
            failures.append((sample_id, "no meta parquet"))
            continue
        frame = pandas.read_parquet(
            meta, columns=["flat_index", "contig_idx", "start", "end", "protein_sequence"]
        ).sort_values("flat_index")
        try:
            aligned = flatten_to_extractor_order(parse_genome_annotation(gff))
            counts = check_against_meta(
                aligned, sample_id, frame["flat_index"], frame["contig_idx"],
                frame["start"], frame["end"], frame["protein_sequence"],
            )
        except GeneAlignmentError as error:  # noqa: PERF203 — one genome must not stop the survey
            failures.append((sample_id, str(error)[:300]))
            continue
        report["genomes"] += 1
        report["genes"] += counts["genes"]
        for key in report["checks"]:
            report["checks"][key] += counts[key]
    report["failures"] = failures
    return report


def test_every_genome_reproduces_the_numbering_the_parquets_are_keyed_on(alignment_report):
    """⭐ The whole cohort, on four independent things at once.

    Measured over all 280 genomes: **1,436,421 genes**, every one matching on gene count, contig
    index, coordinates AND a byte-exact protein string. That last one is the only check that can
    catch a frame shift, and a frame shift is invisible on a page.
    """
    assert alignment_report["failures"] == [], alignment_report["failures"][:5]
    expected = 280 if FULL_COHORT else SAMPLED_GENOMES
    assert alignment_report["genomes"] == expected, "some genomes were skipped without failing"
    assert alignment_report["genes"] > 3_000 * expected, f"only {alignment_report['genes']:,} genes"


def test_all_four_checks_ran_on_every_gene(alignment_report):
    """⛔ Coverage before the claim: three of the four are per-gene, and all three must have fired."""
    genes = alignment_report["genes"]
    assert alignment_report["checks"] == {
        "contig_index": genes, "coordinates": genes, "protein": genes
    }, "a per-gene check did not run on every gene"


def test_the_gff_carries_more_cds_lines_than_the_parquet_has_genes(alignment_report):
    """⚠ The reason the two cannot be zipped, asserted rather than remembered."""
    surpluses = []
    for gff in _genomes()[:10]:
        annotation = parse_genome_annotation(gff)
        aligned = flatten_to_extractor_order(annotation)
        surpluses.append(len(annotation.coding_features) - len(aligned.genes))

    assert all(surplus >= 0 for surplus in surpluses)
    assert any(surplus > 0 for surplus in surpluses), (
        "no genome lost a CDS to the translation-dependent filters — the chain is not being "
        "exercised, and a zip would pass this cohort by accident"
    )


def test_contig_index_enumerates_only_contigs_that_carry_a_cds(alignment_report):
    """⛔ It is a position over CDS-bearing contigs, not a contig number, and the gap is large."""
    gff = _genomes()[0]
    aligned = flatten_to_extractor_order(parse_genome_annotation(gff))

    with_cds = len(aligned.contig_names_in_index_order)
    in_assembly = len(aligned.contig_lengths_by_seqid)
    assert with_cds < in_assembly, f"{with_cds} of {in_assembly} contigs carry a CDS"
    assert {gene.contig_index for gene in aligned.genes} == set(range(with_cds))
    # ⚠ total_base_count is ALL of them — a genome property, not an annotation one.
    assert aligned.total_base_count > sum(
        aligned.contig_lengths_by_seqid[name] for name in aligned.contig_names_in_index_order
    )


def test_a_misalignment_is_refused_rather_than_written(alignment_report):
    """The gate's own failure path — an off-by-one must raise, and say which check caught it."""
    gff = _genomes()[0]
    sample_id = gff.parent.name
    aligned = flatten_to_extractor_order(parse_genome_annotation(gff))
    frame = pandas.read_parquet(
        DATA_ROOT / "proc" / "embeddings" / "meta" / f"{sample_id}_meta.parquet",
        columns=["flat_index", "contig_idx", "start", "end", "protein_sequence"],
    ).sort_values("flat_index")

    shifted = frame.iloc[1:].reset_index(drop=True)
    shifted["flat_index"] = range(len(shifted))
    with pytest.raises(GeneAlignmentError, match="does not match the one that wrote flat_index"):
        check_against_meta(
            aligned, sample_id, shifted["flat_index"], shifted["contig_idx"],
            shifted["start"], shifted["end"], shifted["protein_sequence"],
        )

    # Same length, one protein corrupted: the count check passes and the protein check must not.
    tampered = frame.copy()
    tampered.loc[tampered.index[5], "protein_sequence"] = "MMMM"
    with pytest.raises(GeneAlignmentError, match="translates to"):
        check_against_meta(
            aligned, sample_id, tampered["flat_index"], tampered["contig_idx"],
            tampered["start"], tampered["end"], tampered["protein_sequence"],
        )
