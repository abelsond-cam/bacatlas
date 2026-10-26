"""One genome's gene(s) at one locus — DNA, flanks and protein, read from the original GFF.

⭐ **The database holds no bases, and that is the design.** `.nseq` — a custom 8-byte-magic,
2-bit-packed binary with its own schema doc, its own build job and 290 MB of derived files — existed
for exactly one reason: GitHub Pages applies `Range` to the *compressed* stream, so byte offsets
return plausible wrong bytes with a 206 and no error, and the browser therefore had to fetch whole
files and decode DNA itself. **A server has no such constraint.** So the format is retired rather
than ported, and a gene's nucleotides are sliced from the file Bakta actually wrote.

⚠ **Parsed on demand, cached per process.** A gzipped bacterial GFF with its FASTA block is ~1.5–3
MB and parses in tens of milliseconds — against a tab the reader has to click. The cache holds the
contig sequences only, not the CDS features: the features are already in the database, and keeping
them would double the footprint of the one thing here that is measured in megabytes.

⭐ **Where the bases come from is no longer this module's business** — it asks
`sequence_source.ContigSequenceProvider`, which tries the held Bakta GFF, then a disk cache, then
AllTheBacteria. The held copy still wins, so a genome whose GFF is on this server behaves exactly as
before and enabling streaming can only change the answer for genomes that previously failed.

⚠ **A genome with no `bakrep_dataset_id` used to be unreadable by definition** — that value names the
directory the GFF sits in, no parquet carries it, and it is captured from the path at ingest. Under
streaming its absence means only "we hold no copy of this one", which is the ordinary case for a
streamed genome rather than a fault. The old hard failure is therefore gone from here and the
provider decides; what must NOT come back is a message that calls a streamed genome's missing
directory an error.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.gff.gene_sequence_reader import GeneSequenceView, read_gene_sequence
from bacatlas_backend.models.gene import Gene, GeneLocusMembership
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.sequence_source.contig_sequence_provider import (
    ContigSequenceProvider,
    GenomeSequenceRequest,
    SequenceSourceUnavailable,
)

#: Bases either side of the gene. The published page's own window.
FLANK_LENGTH = 100


class SequenceUnavailable(LookupError):
    """The gene exists but its bases cannot be read — and the reason is always named.

    ⛔ Never conflated with "this genome has no gene here", which is an ANSWER. `app.js:4604`: *"a
    sequence panel that fails silently is one a reader will read as 'this genome has nothing here',
    which is a different claim and a false one."*

    ⚠ **`is_retryable` is new with assembly streaming, and it is not cosmetic.** A genome ATB does
    not have will never arrive; a fetch that timed out probably will. Both were a bare 503 while
    every genome was held locally and the only cause was a missing file. The endpoint turns this
    into `Retry-After`, so the distinction has to survive the raise.
    """

    def __init__(self, message: str, *, is_retryable: bool = False) -> None:
        super().__init__(message)
        self.is_retryable = is_retryable


@dataclass(frozen=True)
class GeneSequenceRow:
    """One gene copy, with everything the tab prints about where it came from."""

    #: ⭐ ρ > 1 puts one genome in a locus twice, so a copy is `n of m` and the page says so.
    copy_ordinal: int
    copy_count: int
    flat_index: int
    #: ⛔ The contig's NAME, never `contig_index + 1` — the index enumerates contigs that HAVE a
    #: CDS, so index 269 is `contig00324` on SAMEA103923484.
    contig_name: str
    seqid: str
    start_position: int
    end_position: int
    strand: str
    strand_is_observed: bool
    gff_phase: int
    is_five_prime_partial: bool
    length_nt: int
    #: What the database stored at ingest, so the page can show it without opening a file — and so a
    #: divergence from the freshly-sliced sequence is visible rather than assumed away.
    #: ⚠ `stored_gc_percent` was here until 2026-10-26 and is gone with the column; GC is computed
    #: from the sliced bases in `sequence.gc_percent`, which is the only GC any client ever read.
    stored_protein_length_aa: int | None
    bakta_gene_symbol: str | None
    bakta_product: str | None
    locus_tag: str | None
    sequence: GeneSequenceView


def load_gene_sequences(
    session: Session,
    *,
    pangenome_id: int,
    node_label: str,
    sample_id: str,
    sequence_provider: ContigSequenceProvider,
    flank_length: int = FLANK_LENGTH,
) -> list[GeneSequenceRow]:
    """Every copy this genome has at this locus, in genome order.

    ⛔ Returns an EMPTY LIST where the genome simply has no gene here — an answer, and the caller
    renders it as one. It RAISES `SequenceUnavailable` only where a gene exists and its bases could
    not be read, which is a different thing and must read differently on the page.

    ⚠ **One statement**, then one assembly. The join reaches `genome_contig` for the seqid rather
    than deriving it, because two GFF seqids can map to one `contig_index` — a short-read assembler
    emits byte-identical duplicate contigs and the extractor keeps only one — so the map is stored as
    measured at ingest and never re-derived per request.

    ⭐ **Where the bases come from is `sequence_provider`'s business, not this function's**, and the
    seqid is why that works: the ATB assembly's defline carries the same first token our GFF writes
    (`SAMN04393397.contig00001` on both, verified on all 280 probe genomes), so the stored
    `genome_contig.seqid` joins against either source unchanged.
    """
    rows = session.execute(
        select(
            Gene.flat_index,
            Gene.start_position,
            Gene.end_position,
            Gene.strand,
            Gene.strand_is_observed,
            Gene.gff_phase,
            Gene.is_five_prime_partial,
            Gene.length_nt,
            Gene.protein_length_aa,
            Gene.bakta_gene_symbol,
            Gene.bakta_product,
            Gene.locus_tag,
            GenomeContig.seqid,
            GenomeContig.contig_name,
            Genome.bakrep_dataset_id,
            Genome.sample_id,
        )
        .join(GeneLocusMembership, GeneLocusMembership.genome_id == Gene.genome_id)
        .join(Genome, Genome.genome_id == Gene.genome_id)
        .join(
            GenomeContig,
            (GenomeContig.genome_id == Gene.genome_id)
            & (GenomeContig.contig_index == Gene.contig_index),
        )
        .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
        .where(
            GeneLocusMembership.pangenome_id == pangenome_id,
            GeneLocusMembership.flat_index == Gene.flat_index,
            Locus.node_label == node_label,
            Locus.pangenome_id == pangenome_id,
            Genome.sample_id == sample_id,
        )
        .order_by(Gene.flat_index)
    ).all()

    if not rows:
        return []

    # ⚠ `bakrep_dataset_id` being absent is no longer fatal — it means only that this genome has no
    # HELD copy, which a streamed one never does. The provider decides whether that is readable.
    try:
        sequences = sequence_provider.contig_sequences(
            GenomeSequenceRequest(sample_id=sample_id, bakrep_dataset_id=rows[0].bakrep_dataset_id)
        )
    except SequenceSourceUnavailable as error:
        raise SequenceUnavailable(str(error), is_retryable=error.is_retryable) from error

    out: list[GeneSequenceRow] = []
    for ordinal, row in enumerate(rows, start=1):
        contig = sequences.get(row.seqid)
        if contig is None:
            # ⚠ Named, not skipped. A silently dropped copy turns "2 of 2" into "1 of 1".
            # ⛔ "assembly", not "annotation file": under streaming the bases and the coordinates
            # come from different places, so a seqid the assembly lacks is a mismatch between the
            # two — which is a far more interesting failure than a missing file and must not be
            # described as one.
            raise SequenceUnavailable(
                f"contig {row.seqid} is not in {sample_id}'s assembly, so this gene's bases "
                "cannot be read"
            )
        out.append(
            GeneSequenceRow(
                copy_ordinal=ordinal,
                copy_count=len(rows),
                flat_index=row.flat_index,
                contig_name=row.contig_name,
                seqid=row.seqid,
                start_position=row.start_position,
                end_position=row.end_position,
                strand=row.strand,
                strand_is_observed=row.strand_is_observed,
                gff_phase=row.gff_phase,
                is_five_prime_partial=row.is_five_prime_partial,
                length_nt=row.length_nt,
                stored_protein_length_aa=row.protein_length_aa,
                bakta_gene_symbol=row.bakta_gene_symbol,
                bakta_product=row.bakta_product,
                locus_tag=row.locus_tag,
                sequence=read_gene_sequence(
                    contig,
                    start_position=row.start_position,
                    end_position=row.end_position,
                    strand=row.strand,
                    phase=row.gff_phase,
                    is_five_prime_partial=row.is_five_prime_partial,
                    flank_length=flank_length,
                ),
            )
        )
    return out
