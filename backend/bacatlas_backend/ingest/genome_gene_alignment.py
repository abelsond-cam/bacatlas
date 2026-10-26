"""Reproduce the extractor's gene order, and PROVE it against the parquet before writing a row.

⛔ **`flat_index` is a running counter over the CDS the extractor KEPT, and it cannot be inferred
from coordinates.** Every embedding matrix, every `{BS}_*.parquet` and every assignment file is
keyed on it. Admit or drop one CDS the extractor did not and every index after it names a different
gene — which puts one gene's sequence, product, UniRef50 family and locus membership under another
gene's name, with nothing on any page looking wrong.

Measured on 25 probe genomes: the GFFs carry **6–18 more CDS lines than their meta parquets**
(mean 11.2). So the two cannot be zipped, and the filter chain has to be reproduced:

1. drop a `pseudo` CDS  ── in the parser
2. drop a CDS whose contig is absent from the `##FASTA` block
3. drop one whose translation is empty
4. drop one carrying an internal stop
5. bucket what survives by seqid **in encounter order**, then concatenate the buckets

⛔ Step 5 is why the contig order is an enumeration and not a vote. The first version of nuna's own
build inferred it from `bakta_products.contig_map`'s modal coordinate vote and **failed the gate on
27 of 100 E. coli genomes**: that vote keeps only `(start, end)` pairs unique within a genome, so a
short contig whose only CDS shares coordinates with another gets no vote and cannot be named at all.
Enumeration has no such hole — every kept CDS names its own contig, so the map is total.

⭐ **And none of that is trusted.** `check_against_meta` re-derives the answer from the parquet on
four independent things at once — gene count, contig index, coordinates, and the protein string —
and refuses the genome on any disagreement. A misaligned flatten is exactly the failure a downstream
test cannot see, so it is caught here or not at all.
"""

from __future__ import annotations

from dataclasses import dataclass

from bacatlas_backend.gff.gene_sequence_reader import (
    reverse_complement,
    translate_coding_sequence,
)
from bacatlas_backend.gff.gff_cds_parser import CodingFeature, ParsedGenomeAnnotation


class GeneAlignmentError(ValueError):
    """The flatten does not reproduce the numbering the parquets are keyed on."""


@dataclass(frozen=True)
class AlignedGene:
    """One gene at its own `flat_index`, with everything only the GFF can supply."""

    flat_index: int
    contig_index: int
    seqid: str
    start_position: int
    end_position: int
    strand: str
    phase: int
    is_five_prime_partial: bool
    locus_tag: str | None
    protein_sequence: str

    @property
    def length_nt(self) -> int:
        """1-based inclusive, and `end` includes the stop codon."""
        return self.end_position - self.start_position + 1


@dataclass(frozen=True)
class GenomeAlignment:
    """One genome's genes in extractor order, plus the contig map that ordering defines."""

    genes: tuple[AlignedGene, ...]
    #: seqid in `contig_index` order — ⛔ **only contigs carrying a kept CDS**, which is what
    #: `contig_index` enumerates. Measured on SAMEA103923484: 270 of the assembly's 348.
    contig_names_in_index_order: tuple[str, ...]
    #: Every contig in the `##FASTA` block, including those with no CDS — the genome's real extent.
    contig_lengths_by_seqid: dict[str, int]

    @property
    def total_base_count(self) -> int:
        """The genome's real extent, over every contig in the `##FASTA` block.

        ⚠ **Not** the sum of the contigs `contig_index` enumerates: the extractor keeps only those
        carrying a CDS, and the two differ on every genome measured — 5,025,405 bases over 348
        contigs against 4,880,150 over 270 on SAMEA103923484.
        """
        return sum(self.contig_lengths_by_seqid.values())


def flatten_to_extractor_order(annotation: ParsedGenomeAnnotation) -> GenomeAlignment:
    """The CDS the extractor kept, numbered as it numbered them.

    ⚠ A `compute_gc` flag used to let a caller skip the GC pass. GC is no longer computed here at
    all — its column was dropped 2026-10-26 — so the flag went with it rather than surviving as a
    parameter that selects between two identical behaviours.
    """
    if not annotation.carries_sequence:
        raise GeneAlignmentError(
            "the GFF has no ##FASTA block, so the extractor's translation-dependent filters cannot "
            "be reproduced and flat_index cannot be established. Every probe GFF carries one; a "
            "cohort that does not needs the assembly FASTA wired in as a second source."
        )

    contigs = annotation.contig_sequences
    order: list[str] = []
    per_contig: dict[str, list[AlignedGene]] = {}

    for feature in annotation.coding_features:
        sequence = contigs.get(feature.seqid)
        if sequence is None:
            continue  # the extractor's n_skipped_missing_contig
        coding = sequence[feature.start_position - 1 : feature.end_position]
        if feature.strand == "-":
            coding = reverse_complement(coding)
        protein = translate_coding_sequence(coding, feature.phase, feature.is_five_prime_partial)
        # `not protein` covers both of the extractor's short-sequence exits; a residual '*' is its
        # internal-stop skip, the terminal one having already been dropped by the translation.
        if not protein or "*" in protein:
            continue
        if feature.seqid not in per_contig:
            order.append(feature.seqid)
            per_contig[feature.seqid] = []
        per_contig[feature.seqid].append(
            _aligned(feature, len(order) - 1, protein)
        )

    genes: list[AlignedGene] = []
    for seqid in order:
        for gene in per_contig[seqid]:
            genes.append(_renumber(gene, len(genes)))

    return GenomeAlignment(
        genes=tuple(genes),
        contig_names_in_index_order=tuple(order),
        contig_lengths_by_seqid={seqid: len(sequence) for seqid, sequence in contigs.items()},
    )


def _aligned(feature: CodingFeature, contig_index: int, protein: str) -> AlignedGene:
    # ⚠ GC is no longer computed here. The column it fed was dropped 2026-10-26, and the rule it
    # protected — NULL rather than 0.0 for a span with no unambiguous base, because under this
    # schema 0.0 is a MEASURED zero — now lives only where GC still is: `gene_sequence_reader`
    # returns 0.0 for such a span, and the serving path reports it as computed rather than stored.
    return AlignedGene(
        flat_index=-1,  # assigned in the concatenation pass, which is where the number is decided
        contig_index=contig_index,
        seqid=feature.seqid,
        start_position=feature.start_position,
        end_position=feature.end_position,
        strand=feature.strand,
        phase=feature.phase,
        is_five_prime_partial=feature.is_five_prime_partial,
        locus_tag=feature.locus_tag,
        protein_sequence=protein,
    )


def _renumber(gene: AlignedGene, flat_index: int) -> AlignedGene:
    return AlignedGene(**{**gene.__dict__, "flat_index": flat_index})


def check_against_meta(
    alignment: GenomeAlignment,
    sample_id: str,
    meta_flat_index,
    meta_contig_index,
    meta_start,
    meta_end,
    meta_protein_sequence,
) -> dict[str, int]:
    """Refuse the genome unless the flatten reproduces the parquet on all four things at once.

    ⛔ **Four independent checks, not one.** Each alone can pass on a misalignment: the counts match
    whenever the same number of genes shifted; coordinates can coincide between neighbouring genes
    on a short contig; the contig index is constant across a whole contig's worth of rows. Together
    they cannot. Returns the counts checked, so a caller reports coverage rather than a bare pass.

    ⚠ `meta_*` are the parquet columns, and `flat_index` must be `0..n-1` contiguous — the gene
    table is positional, and a gap in it means the file was written by something else.
    """
    counts = {
        "genes": len(alignment.genes),
        "contig_index": 0,
        "coordinates": 0,
        "protein": 0,
    }
    flat = [int(value) for value in meta_flat_index]
    if flat != list(range(len(flat))):
        raise GeneAlignmentError(
            f"{sample_id}: the meta parquet's flat_index is not 0..{len(flat) - 1} contiguous — the "
            "gene table is positional, so nothing here can be aligned to it"
        )
    if len(alignment.genes) != len(flat):
        raise GeneAlignmentError(
            f"{sample_id}: the GFF's filter chain yields {len(alignment.genes):,} genes but the "
            f"parquet has {len(flat):,}. The chain does not match the one that wrote flat_index, so "
            "every index after the first difference names a different gene."
        )

    contig_index = [int(v) for v in meta_contig_index]
    start = [int(v) for v in meta_start]
    end = [int(v) for v in meta_end]
    protein = [str(v) for v in meta_protein_sequence]

    for gene in alignment.genes:
        j = gene.flat_index
        if gene.contig_index != contig_index[j]:
            raise GeneAlignmentError(
                f"{sample_id}: gene {j} is on contig_index {gene.contig_index} by enumeration and "
                f"{contig_index[j]} in the parquet — the contig ORDER differs, which shifts a whole "
                "contig's worth of genes"
            )
        counts["contig_index"] += 1
        if (gene.start_position, gene.end_position) != (start[j], end[j]):
            raise GeneAlignmentError(
                f"{sample_id}: gene {j} spans {gene.start_position}-{gene.end_position} in the GFF "
                f"and {start[j]}-{end[j]} in the parquet"
            )
        counts["coordinates"] += 1
        if gene.protein_sequence != protein[j]:
            raise GeneAlignmentError(
                f"{sample_id}: gene {j} ({gene.locus_tag}) translates to {len(gene.protein_sequence)} "
                f"aa and the parquet holds {len(protein[j])} aa. This is the check that catches a "
                "frame shift, which no page would show."
            )
        counts["protein"] += 1
    return counts
