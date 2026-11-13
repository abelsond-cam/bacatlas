"""One genome's gene(s) at one locus — DNA, flanks and protein, read from the original GFF.

⭐ **The database holds no bases, and that is the design.** `.nseq` — a custom 8-byte-magic,
2-bit-packed binary with its own schema doc, its own build job and 290 MB of derived files — existed
for exactly one reason: GitHub Pages applies `Range` to the *compressed* stream, so byte offsets
return plausible wrong bytes with a 206 and no error, and the browser therefore had to fetch whole
files and decode DNA itself. **A server has no such constraint.** So the format is retired rather
than ported, and a gene's nucleotides are sliced from the file Bakta actually wrote.

⭐ **The coordinates come from the genome-grain blocks, not from a row per gene** (design of record:
`nuna/docs/design/storing_the_gene_layer.md` §3, *"How each question is answered"*, the Sequence-tab
row). The algorithm is the one the published `app.js` already used over the committed `.nseq` +
`.loci`: take this genome's `gene_locus_ordinals`, scan it for the locus's `catalogue_ordinal`, and
read the hits' rows out of `gene_spans`. ~4,900 `int32` scanned in place rather than an inverse index
the database would have to keep consistent — and because `gene_spans` is byte-identical to the
`.nseq` gene table, the page's decoder and this one are differential oracles for each other.

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

⛔ **No numpy here, and no hand-rolled `struct`.** This is the serving path, where numpy is not
installed; every byte is read through `bacatlas_backend.gene_blocks`, which is stdlib-only and is the
one decoder the parity suite has proved against the page's.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from bacatlas_backend.gene_blocks import (
    SPAN_ENCODING_NSEQ_GENE_TABLE,
    GeneBlockError,
    GeneSpanTable,
    decode_int32_vector,
    positions_of,
)
from bacatlas_backend.gff.gene_sequence_reader import GeneSequenceView, read_gene_sequence
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
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
    """One gene copy, with everything the tab prints about where it came from.

    ⚠ Four fields left on 2026-11-13 with the gene-grain tables: `stored_protein_length_aa`,
    `bakta_gene_symbol`, `bakta_product` and `locus_tag`. Nothing rendered any of them (David,
    2026-10-26 / 2026-10-29); the symbol and product live at the LOCUS as rollups, and the protein
    length is read off the sliced sequence. `stored_gc_percent` went the same way on 2026-10-26.
    """

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
    sequence: GeneSequenceView


def _genome_at_locus(session: Session, *, pangenome_id: int, node_label: str, sample_id: str):
    """ONE statement: the genome, its two blocks, the locus's ordinal, and the genome's contig map.

    ⭐ Inner joins on `genome_locus_assignment` and `locus`, so a genome this catalogue never placed
    and a locus this catalogue does not have both come back as NO ROW — and neither moves a block
    across the wire to learn it. ⚠ The join onto `locus` is deliberately a filtered cross join: the
    locus is addressed by `(pangenome_id, node_label)`, which is unique, so it contributes at most one
    row and nothing about it depends on the genome.

    ⚠ `genome_genes` is an OUTER join, the one place a missing row is a fault rather than an answer:
    a genome placed in a catalogue whose coordinates were never loaded is a database defect, and the
    caller refuses it by name instead of reporting "no gene here" — UNLESS `genome.gene_table_omitted_reason`
    says why there is none (the genome was named at ingest as one nuna's `.nseq` gate refused). That
    reason rides on this same statement, so the named state costs no second one.

    ⭐ The contig map rides along as a correlated `json_agg` over `genome_contig` — O(this genome's
    contigs with a CDS, ~300 on a draft assembly) and a fraction of the blocks' own size. Reading it
    here rather than in a second statement keeps the tab at one statement, which is the budget
    `test_the_sequence_view_is_ONE_statement_and_ONE_file` holds it to.
    """
    contig_map = (
        select(
            func.json_agg(
                func.json_build_object(
                    "contig_index",
                    GenomeContig.contig_index,
                    "seqid",
                    GenomeContig.seqid,
                    "contig_name",
                    GenomeContig.contig_name,
                )
            )
        )
        .where(GenomeContig.genome_id == Genome.genome_id)
        .correlate(Genome)
        .scalar_subquery()
    )
    return session.execute(
        select(
            Genome.genome_id,
            Genome.bakrep_dataset_id,
            Genome.sample_id_kind,
            Genome.bare_accession,
            Genome.strand_is_observed,
            Genome.gene_table_omitted_reason,
            GenomeGenes.gene_count,
            GenomeGenes.span_encoding,
            GenomeGenes.gene_spans,
            GenomeLocusAssignment.gene_count.label("assigned_gene_count"),
            GenomeLocusAssignment.gene_locus_ordinals,
            Locus.catalogue_ordinal,
            contig_map.label("contigs"),
        )
        .select_from(Genome)
        .join(
            GenomeLocusAssignment,
            (GenomeLocusAssignment.genome_id == Genome.genome_id)
            & (GenomeLocusAssignment.pangenome_id == pangenome_id),
        )
        .join(Locus, (Locus.pangenome_id == pangenome_id) & (Locus.node_label == node_label))
        .outerjoin(GenomeGenes, GenomeGenes.genome_id == Genome.genome_id)
        .where(Genome.sample_id == sample_id)
    ).one_or_none()


def _contigs_by_index(sample_id: str, contigs: list[dict] | None) -> dict[int, tuple[str, str]]:
    """`contig_index → (seqid, contig_name)` from the `json_agg` payload; `{}` for a genome with none.

    ⚠ One seqid per index. `ingest_genome_gene_table._contig_rows` enumerates
    `contig_names_in_index_order`, so the map is a function; a second row for one index would mean
    the loader changed shape, and it is refused rather than letting one of the two win silently.
    """
    out: dict[int, tuple[str, str]] = {}
    for contig in contigs or ():
        index = int(contig["contig_index"])
        if index in out:
            raise GeneBlockError(
                f"{sample_id}: genome_contig holds two rows for contig_index {index} "
                f"({out[index][0]!r} and {contig['seqid']!r}); the gene table can name only one."
            )
        out[index] = (str(contig["seqid"]), str(contig["contig_name"]))
    return out


def _refuse_named_genome_without_gene_table(sample_id: str, found) -> list[GeneSequenceRow]:
    """A genome NAMED as holding no gene table: no gene here is `[]`, a gene here is a named 503.

    `SequenceUnavailable`, NOT retryable, carrying the stored reason.

    ⭐ The assignment block is still this genome's — written at `genome.coding_gene_count`, its `_meta`
    row count — so *"this genome has no gene at this locus"* stays an answer it can give truthfully.
    Where it has one, the bases cannot be read because the coordinates were never held, and that is
    said in those words: a 503 without `Retry-After`, because no retry will ever produce a table.
    """
    ordinals = decode_int32_vector(
        found.gene_locus_ordinals, name="gene_locus_ordinals", expected_count=found.assigned_gene_count
    )
    copies = len(positions_of(ordinals, found.catalogue_ordinal))
    if not copies:
        return []
    raise SequenceUnavailable(
        f"{sample_id} has {'a gene' if copies == 1 else f'{copies} genes'} at this locus, but BacAtlas "
        f"holds no gene table for this genome — {found.gene_table_omitted_reason} — so "
        f"{'its' if copies == 1 else 'their'} coordinates and bases cannot be shown",
        is_retryable=False,
    )


def load_gene_sequences(
    session: Session,
    *,
    pangenome_id: int,
    node_label: str,
    sample_id: str,
    sequence_provider: ContigSequenceProvider,
    flank_length: int = FLANK_LENGTH,
) -> list[GeneSequenceRow]:
    """Every copy this genome has at this locus, in `flat_index` order.

    ⛔ Returns an EMPTY LIST where the genome simply has no gene here — an answer, and the caller
    renders it as one. That includes a genome this catalogue never placed (no
    `genome_locus_assignment` row) and a locus the catalogue lacks: the endpoint has already turned
    an unknown genome or locus into a named 404 before calling this. It RAISES `SequenceUnavailable`
    only where a gene exists and its bases could not be read, which is a different thing and must
    read differently on the page — including a genome NAMED as holding no gene table, which raises
    it non-retryably with the stored reason.

    ⛔ **A malformed row is REFUSED BY NAME, never served as "no gene here".** A `genome_genes` row
    that is missing, in a layout this decoder does not know, or whose length disagrees with its own
    `gene_count`, raises `GeneBlockError` — which the API answers as a 500 with the message in the
    server log. The one thing it must never become is `genes: []`.

    ⚠ **One statement**, then one assembly. The statement carries the contig map rather than
    deriving it, because `contig_index` enumerates contigs that HAVE a CDS and `contig_index + 1`
    is a real contig the gene is not on for 28.6 % of contigs; the map is stored as measured at
    ingest and never re-derived per request.

    ⭐ **Where the bases come from is `sequence_provider`'s business, not this function's**, and the
    seqid is why that works: the ATB assembly's defline carries the same first token our GFF writes
    (`SAMN04393397.contig00001` on both, verified on all 280 probe genomes), so the stored
    `genome_contig.seqid` joins against either source unchanged.
    """
    found = _genome_at_locus(session, pangenome_id=pangenome_id, node_label=node_label, sample_id=sample_id)
    if found is None:
        return []
    if found.gene_count is None and found.gene_table_omitted_reason is not None:
        return _refuse_named_genome_without_gene_table(sample_id, found)
    if found.gene_count is None:
        raise GeneBlockError(
            f"{sample_id} is placed in pangenome {pangenome_id} but has no genome_genes row, so its "
            "coordinates cannot be read. The genome layer was not loaded for it."
        )
    if found.span_encoding != SPAN_ENCODING_NSEQ_GENE_TABLE:
        raise GeneBlockError(
            f"{sample_id}: genome_genes.span_encoding is {found.span_encoding}, and this decoder knows "
            f"only {SPAN_ENCODING_NSEQ_GENE_TABLE} (nuna's .nseq gene-table layout)."
        )
    if found.assigned_gene_count != found.gene_count:
        raise GeneBlockError(
            f"{sample_id}: genome_locus_assignment says {found.assigned_gene_count} genes and "
            f"genome_genes says {found.gene_count}. The two blocks are indexed by the same flat_index, "
            "so they cannot describe this genome together."
        )

    # ⭐ The page's own algorithm: scan the ordinals for the locus, then read those rows of the spans.
    ordinals = decode_int32_vector(
        found.gene_locus_ordinals, name="gene_locus_ordinals", expected_count=found.gene_count
    )
    flat_indexes = positions_of(ordinals, found.catalogue_ordinal)
    if not flat_indexes:
        return []
    spans = GeneSpanTable(found.gene_spans, gene_count=found.gene_count)
    contigs = _contigs_by_index(sample_id, found.contigs)

    # ⚠ `bakrep_dataset_id` being absent is no longer fatal — it means only that this genome has no
    # HELD copy, which a streamed one never does. The provider decides whether that is readable.
    # ⭐ The key KIND picks the remote station (ATB for a BioSample, NCBI for an assembly stem), and the
    # genome's contig names let the provider re-key an NCBI assembly (`NZ_CP…` ↔ `CP…`) into the
    # namespace the lookup below uses — so that lookup is the same for every source. ⚠ ALL of the
    # genome's contigs, not just this locus's: the `NZ_` rule is proved on the whole map, as nuna did.
    try:
        sequences = sequence_provider.contig_sequences(
            GenomeSequenceRequest(
                sample_id=sample_id,
                bakrep_dataset_id=found.bakrep_dataset_id,
                sample_id_kind=found.sample_id_kind,
                bare_accession=found.bare_accession,
                wanted_seqids=frozenset(seqid for seqid, _ in contigs.values()),
            )
        )
    except SequenceSourceUnavailable as error:
        raise SequenceUnavailable(str(error), is_retryable=error.is_retryable) from error

    out: list[GeneSequenceRow] = []
    for ordinal, flat_index in enumerate(flat_indexes, start=1):
        contig_index = spans.contig_index(flat_index)
        if contig_index not in contigs:
            raise GeneBlockError(
                f"{sample_id}: gene {flat_index} sits on contig_index {contig_index}, which has no "
                "genome_contig row. The contig map and the gene table were not loaded together."
            )
        seqid, contig_name = contigs[contig_index]
        contig = sequences.get(seqid)
        if contig is None:
            # ⚠ Named, not skipped. A silently dropped copy turns "2 of 2" into "1 of 1".
            # ⛔ "assembly", not "annotation file": under streaming the bases and the coordinates
            # come from different places, so a seqid the assembly lacks is a mismatch between the
            # two — which is a far more interesting failure than a missing file and must not be
            # described as one.
            raise SequenceUnavailable(
                f"contig {seqid} is not in {sample_id}'s assembly, so this gene's bases cannot be read"
            )
        start_position = spans.start_position(flat_index)
        end_position = spans.end_position(flat_index)
        strand = spans.strand(flat_index)
        gff_phase = spans.gff_phase(flat_index)
        is_five_prime_partial = spans.is_five_prime_partial(flat_index)
        out.append(
            GeneSequenceRow(
                copy_ordinal=ordinal,
                copy_count=len(flat_indexes),
                flat_index=flat_index,
                contig_name=contig_name,
                seqid=seqid,
                start_position=start_position,
                end_position=end_position,
                strand=strand,
                strand_is_observed=found.strand_is_observed,
                gff_phase=gff_phase,
                is_five_prime_partial=is_five_prime_partial,
                length_nt=spans.length_nt(flat_index),
                sequence=read_gene_sequence(
                    contig,
                    start_position=start_position,
                    end_position=end_position,
                    strand=strand,
                    phase=gff_phase,
                    is_five_prime_partial=is_five_prime_partial,
                    flank_length=flank_length,
                ),
            )
        )
    return out
