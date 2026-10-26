"""The JOINT ±5 view — whole neighbourhoods that some set of genomes actually has.

⛔ **Not the same thing as `locus_offset_occupant`, and not per gene.** An arrangement is a complete
ten-slot vector carried together. The marginal view reduces each offset independently, and where a
locus occurs in more than one arrangement *"the marginal describes none of them, and read
left-to-right it is a gene order nothing has."* The page draws the TRACK from here and the BAR under
each block from there.

⚠ **There is no per-gene link any more.** `gene_locus_membership.locus_arrangement_id` carried one
and was dropped 2026-10-26, measured unread by every query, serialiser and test. Which genomes an
arrangement holds is `member_genome_ids` here; which genes had no recorded neighbourhood is counted
at ingest from the frame's own `arrangement_rank`. The per-gene rank returns as
`genome_locus_assignment.gene_arrangement_ranks` — see
`nuna/docs/design/storing_the_gene_layer.md`.
"""

from __future__ import annotations

from sqlalchemy import BigInteger, Boolean, ForeignKey, Index, Integer, UniqueConstraint
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base

#: `catalogue_ordinal * 2 + same_strand`, or **−1 where the contig ends**.
#: ⛔ −1 is a VALUE, not a wildcard: a member truncated at a contig edge has a genuinely different
#: *observed* neighbourhood and is not folded into the untruncated group, which would claim an
#: observation never made. The API serialises it as `null` + `absence_reason: "contig_end"`, never
#: as a bare −1, because −1 means five different things across the payload.
CONTIG_END_SLOT_CODE = -1


class LocusArrangement(Base):
    """One distinct ±5 neighbourhood at one locus."""

    __tablename__ = "locus_arrangement"
    __table_args__ = (
        UniqueConstraint("locus_id", "rank_within_locus"),
        Index("ix_locus_arrangement__locus_id", "locus_id"),
        # ⭐ GIN over the genome membership array: this is what answers "which arrangements does the
        # anchored genome occupy at this locus" with `member_genome_ids @> ARRAY[?]`, replacing a
        # binary search over a ~490k-entry Int32Array in the browser.
        Index(
            "ix_locus_arrangement__member_genome_ids",
            "member_genome_ids",
            postgresql_using="gin",
        ),
    )

    locus_arrangement_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)
    #: Denormalised for LIST partitioning at scale — ~49 M rows at 80,000 genomes.
    pangenome_id: Mapped[int] = mapped_column(BigInteger, nullable=False)

    #: Rank order is stable and is what makes an uncapped export provably additive against a capped
    #: one. ⚠ The API compares selection on THIS integer, not on object identity — `arrShown` uses
    #: `a === sel` today, which a framework's recomputed array silently breaks.
    rank_within_locus: Mapped[int] = mapped_column(Integer, nullable=False)

    member_gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: ⚠ **Not the same as the gene count when ρ > 1**, and the page says "genomes".
    member_genome_count: Mapped[int] = mapped_column(Integer, nullable=False)

    #: `arr.flip` — this arrangement matches rank 1 better reverse-complemented than as observed.
    #: ⚠ INTRINSIC to the arrangement, not display mirroring: the two are different slot spaces and
    #: conflating them puts a count from one position under a gene from another.
    is_recorded_reverse_complement: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    #: The ten slots, packed, in `OFFSETS` order. Never queried per slot — 490 M rows to hold 490 M
    #: small ints would be a 10× storage multiplier for no query value.
    neighbour_slot_codes: Mapped[list[int]] = mapped_column(ARRAY(Integer), nullable=False)

    #: ⭐ Which genomes carry it, as `genome_id`, ascending. This is `arr.gid`, the single biggest
    #: block in the payload (486,717 entries at 100 genomes) and O(genes) rather than O(loci).
    #: ⛔ **There is deliberately NO uniqueness constraint on (locus, genome).** A genome at ρ > 1
    #: occupies TWO arrangements at one locus, so `sum(member_genome_count)` over a locus may exceed
    #: `locus.member_genome_count` while the UNION equals it. A schema forbidding that would assert
    #: something false, and `COUNT(*)` where `COUNT(DISTINCT genome_id)` is meant returns a
    #: plausible, slightly-too-large number.
    member_genome_ids: Mapped[list[int]] = mapped_column(ARRAY(Integer), nullable=False)

    def __repr__(self) -> str:
        return f"<LocusArrangement locus={self.locus_id} rank={self.rank_within_locus}>"
