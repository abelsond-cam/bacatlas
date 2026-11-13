"""Which locus each of a genome's genes landed in, for one pangenome — one row per (catalogue, genome).

⛔ **This replaces `gene_locus_membership`, and it is keyed by GENOME, not by locus.** Per-locus
parallel arrays are the obvious shape and they fail at scale: a core locus at 430,000 genomes holds
430,000 `(genome, flat_index)` pairs ≈ 5.2 MB in one row, and the one surviving query needs the
genome's *position* in that array — which a GIN index cannot supply — so it would detoast and
`unnest … WITH ORDINALITY` 430,000 entries to return one row. Keyed by genome, one row is ~39 kB and
there is no cliff. Design of record: `nuna/docs/design/storing_the_gene_layer.md` §3B.

⭐ **`gene_locus_ordinals` is nuna's `.loci` sidecar, byte for byte**: `gene_count` little-endian
`int32`, the catalogue ordinal of `flat_index` *j* or `-1` where the gene is not in the catalogue.
Where a `.loci` exists for the genome, the block written at ingest is compared against it; where
none exists it is built from the same memberships frame that used to feed the per-gene table, so the
two routes cannot describe different catalogues.

⚠ `gene_arrangement_ranks` keeps the per-gene link to `locus_arrangement.rank_within_locus` with **no
FK and no index**, so the `ON DELETE SET NULL` stall recorded on the old table (*"the second sat in
`DELETE FROM locus_arrangement` for over 8 minutes and had to be killed"*) is structurally impossible
rather than index-mitigated. `-1` is *this genome has a gene here with no recorded neighbourhood* —
a different sentence from *no gene here*, and the page says it.

⭐ The question *"which genomes carry locus L"* is answered by `locus.present_genome_ids` (GIN), not
by this table — this row answers the inverse, *"which loci does genome G occupy, and where"*.
"""

from __future__ import annotations

from sqlalchemy import CheckConstraint, ForeignKey, Integer, LargeBinary
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.gene_blocks import INT32_VECTOR_STRIDE


class GenomeLocusAssignment(Base):
    """One genome's locus ordinals and arrangement ranks for one pangenome, as two int32 blocks."""

    __tablename__ = "genome_locus_assignment"
    __table_args__ = (
        CheckConstraint(
            f"octet_length(gene_locus_ordinals) = gene_count * {INT32_VECTOR_STRIDE}",
            name="gene_locus_ordinals_is_4_bytes_per_gene",
        ),
        CheckConstraint(
            f"octet_length(gene_arrangement_ranks) = gene_count * {INT32_VECTOR_STRIDE}",
            name="gene_arrangement_ranks_is_4_bytes_per_gene",
        ),
        CheckConstraint(
            "placed_gene_count >= 0 AND placed_gene_count <= gene_count", name="placed_genes_within_gene_count"
        ),
        CheckConstraint(
            "distinct_locus_count >= 0 AND distinct_locus_count <= placed_gene_count",
            name="distinct_loci_within_placed_genes",
        ),
    )

    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), primary_key=True
    )
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)
    #: ⛔ Must equal `genome_genes.gene_count` for the same genome — the two blocks are indexed by the same
    #: `flat_index`. Checked by the loader and by `verify_pangenome_is_servable`, not by a constraint,
    #: because a cross-table equality is not expressible as a CHECK.
    gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Entries of `gene_locus_ordinals` that are not `-1`.
    placed_gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Distinct loci among the placed genes — what `pangenome_genome_locus_count.locus_count` is fed from.
    distinct_locus_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: `gene_count` little-endian int32: `locus.catalogue_ordinal`, or `-1` where unassigned. ≡ nuna's `.loci`.
    gene_locus_ordinals: Mapped[bytes] = mapped_column(LargeBinary, nullable=False)
    #: `gene_count` little-endian int32: `locus_arrangement.rank_within_locus`, or `-1` where the gene has
    #: no recorded neighbourhood (alone on its contig). ⚠ A value, never a wildcard — see the module note.
    gene_arrangement_ranks: Mapped[bytes] = mapped_column(LargeBinary, nullable=False)
