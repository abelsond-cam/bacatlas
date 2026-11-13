"""The gene layer at the genome grain — one row per genome, the coordinates as one packed block.

⛔ **This replaces `gene`, and the grain is the whole point.** Postgres charges ~28 B of tuple header
per row regardless of columns, so one row per gene put the serving schema at ~347 GB at 72,801
genomes (`gene` alone 131 GB at 80,000). One row per genome holding the same facts as a 16 B/gene
block is ~1.2 GB for the same cohort. Design of record: `nuna/docs/design/storing_the_gene_layer.md`
§3A; the codec is `bacatlas_backend.gene_blocks`.

⭐ **`gene_spans` is byte-identical to nuna's published `.nseq` gene table** — four little-endian
`uint32` per gene (`start | end | contig_index | flags`), in `flat_index` order. Three consequences:
nuna's encoder is the only encoder; `app.js` and the server read the same bytes; and a symmetric
encode/decode bug cannot cancel in the parity suite, because the page side decodes with code we did
not write. Where a `.nseq` exists for the genome the block is COPIED from it and the sample name in
its header is checked; where only a GFF exists the block is encoded from the parsed CDS lines with
`encode_gene_spans`, which the 200-genome differential oracle proved identical (0 differences).

⚠ What a reader of this row cannot get any more, by decision (David, 2026-10-26 / 2026-10-29): the
Bakta symbol, product and locus tag per gene, the stored protein length, and GC. None of them was
rendered by any page; the symbol/product live at the LOCUS as rollups, GC is computed from the
sliced bases, and protein length is `(end − start + 1) / 3 − 1` where phase is 0.
"""

from __future__ import annotations

from sqlalchemy import CheckConstraint, ForeignKey, Integer, LargeBinary, SmallInteger
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.gene_blocks import GENE_SPAN_STRIDE, SPAN_ENCODING_NSEQ_GENE_TABLE


class GenomeGenes(Base):
    """One genome's whole gene table as a block. Keyed by `genome_id`; no other index exists or is needed."""

    __tablename__ = "genome_genes"
    __table_args__ = (
        # ⛔ The length CHECK is the only thing the database can verify about an opaque block, and it is
        # what turns a truncated write into a refusal rather than a table of short rows.
        CheckConstraint(
            f"octet_length(gene_spans) = gene_count * {GENE_SPAN_STRIDE}", name="gene_spans_is_16_bytes_per_gene"
        ),
        CheckConstraint("gene_count >= 0", name="gene_count_is_not_negative"),
    )

    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)
    #: Rows in the block — and the one number every consumer reads before touching the bytes.
    gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Which layout the block is in. Only `SPAN_ENCODING_NSEQ_GENE_TABLE` exists; the column is here so a
    #: second layout can never be read as the first by a decoder that did not know it existed.
    span_encoding: Mapped[int] = mapped_column(SmallInteger, nullable=False, default=SPAN_ENCODING_NSEQ_GENE_TABLE)
    #: `gene_count * 16` bytes: per gene `start u32 | end u32 | contig_index u32 | flags u32`, little-endian.
    #: flags: bit 0 minus strand, bit 1 five-prime partial, bits 2–3 the GFF phase.
    gene_spans: Mapped[bytes] = mapped_column(LargeBinary, nullable=False)
