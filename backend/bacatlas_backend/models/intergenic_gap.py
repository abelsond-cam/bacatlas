"""Intergenic regions — a property of an ADJACENCY, not of a position.

⛔ Keyed on the **canonical pair** of flanking loci, which is what makes the whole block additive:
nothing is inserted into any sequence, so the offsets, the marginal and the joint views are all
untouched. Canonical, because the track is drawn in the focal gene's reading frame and a
neighbourhood may be mirrored — an orientation-dependent key would miss on half of them.

⛔⛔ **THE CANONICAL ORDER IS BY `node_label`, WHICH IS TEXT — NOT BY ANY INTEGER.**
`intergenic.gene_adjacencies` sorts the pair with ``left <= right`` on the **node ids**
(`intergenic.py:106-118`), and node ids are text. Everything downstream inherits that ordering.

This is not a stylistic choice, and getting it wrong has already cost the live site. The published
page keys its gap lookup by payload INDEX —
``GAPAT[(x < y ? x : y) + "|" + (x < y ? y : x)]`` (`app.js:75`) — against a map built from the
label-sorted pairs the exporter shipped. Measured on the two published catalogues:

    ecoli  22,838 gaps · 7,379 (32.3 %) have a > b by index · 100 % are label-sorted
    kp     20,544 gaps · 6,108 (29.7 %) have a > b by index · 100 % are label-sorted

Those 7,379 and 6,108 gaps **cannot be found by the page**: it asks for ``"890|3819"`` and the map
holds ``"3819|890"``. `gapBetween` returns null, `gapSlot` returns null, and the track simply draws
no intergenic block — which renders as *"these two genes are adjacent with nothing between them"*,
a claim the data does not make. It is the exact failure mode this codebase names elsewhere: an
absence that reads as a measurement.

So the columns below are named `a`/`b` and **not** `low`/`high`, because there is no numeric
ordering to appeal to, and a `CHECK (low <= high)` on surrogate ids would silently impose a second,
different canonicalisation on top of the real one. The ordering invariant is asserted at ingest
against the labels, and the API resolves a pair by sorting its two labels before it queries.

⚠ **Three things this is not.** A gap is not a locus — no cluster, no embedding, not in the locus
count, not walkable. Its length is a **median over genomes that disagree**, not a distance. And the
page draws gaps only where an arrangement is selected: without one the track shows each position's
marginal mode, and two marginal modes need not be neighbours in any genome, so the space between
them is not an observed gap.
"""

from __future__ import annotations

from sqlalchemy import BigInteger, ForeignKey, Index, Integer, SmallInteger, String, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.column_types import measurement, nan_guards


class IntergenicGap(Base):
    """The region between two adjacent loci, summarised over the genomes that carry the adjacency."""

    __tablename__ = "intergenic_gap"
    __table_args__ = (
        UniqueConstraint("pangenome_id", "flanking_locus_id_a", "flanking_locus_id_b"),
        # ⭐ Both directions indexed, because the canonical order is by LABEL and a caller holding
        # two locus ids cannot tell which is `a` without resolving them.
        Index("ix_intergenic_gap__flanking_locus_id_a", "flanking_locus_id_a"),
        Index("ix_intergenic_gap__flanking_locus_id_b", "flanking_locus_id_b"),
        # ⛔ No `CHECK (a <= b)`. See the module docstring: the pair is canonicalised by node_label,
        # and a numeric CHECK on surrogate ids would assert a DIFFERENT ordering that happens to be
        # satisfiable — the schema would then look correct while holding the pairs the wrong way
        # round for a third of the table.
        #
        # ⛔⛔ **AND NO `CHECK (a <> b)` EITHER — A LOCUS CAN FLANK ITSELF, AND IT IS NOT AN ERROR.**
        # This table briefly carried one, on the reasonable-sounding ground that a gap sits *between
        # two* loci. The real ecoli load refused 12 rows against it, and every one turned out to be
        # a **tandem repeat of a multi-copy locus**: `hokE`, `ldrB`, `ldrD`, `zorO`, `tpr`, `icaT`
        # and a transposase — type I toxin-antitoxin systems, which occur in E. coli as several
        # near-identical copies that sit next to each other. Every one of the twelve has ρ > 1
        # (`size/genomes` from 1.02 to 2.00), which is the signature. kp carries three (all `ldrB`
        # -shaped), and the PUBLISHED payloads ship all fifteen: 12 of 22,838 and 3 of 20,544.
        #
        # So the region between one copy and the next is a real intergenic region with a real
        # measured length, and a constraint forbidding it would have silently dropped the tandem
        # arrays — which are the loci where the spacing actually carries information. Uniqueness of
        # the pair is still enforced by the UniqueConstraint above; that is the invariant that holds.
        *nan_guards("length_variance_score"),
    )

    intergenic_gap_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )
    #: The locus whose `node_label` sorts FIRST as text. ⛔ Not the smaller id — see the docstring.
    flanking_locus_id_a: Mapped[int] = mapped_column(
        ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False
    )
    #: The locus whose `node_label` sorts SECOND as text.
    flanking_locus_id_b: Mapped[int] = mapped_column(
        ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False
    )

    observed_genome_count: Mapped[int] = mapped_column(Integer, nullable=False)

    #: ⛔ **SIGNED. Negative means the two genes OVERLAP** — 18.8 % of adjacent pairs do, at a median
    #: of −4 bases (a shared stop/start codon). Clamping this at zero made *abuts exactly* and
    #: *overlaps by 190 bases* the same number and threw the question away before the page saw it.
    median_signed_length_nt: Mapped[int] = mapped_column(Integer, nullable=False)

    #: ⛔ `q1 == q3` does NOT mean every genome agrees — it means the MIDDLE HALF does. `[1,1,1,1,100]`
    #: has an interquartile range of zero, and an earlier page reported that shape as "identical in
    #: every genome". Only `mn == mx` certifies the strong claim, which is why both pairs exist.
    #: ⚠ Absent means NOT MEASURED, never a spread of zero.
    quartile1_signed_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)
    quartile3_signed_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)
    minimum_signed_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)
    maximum_signed_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: ⛔ **NULL = not measured; 0.0 = every genome agrees.** The payload's sparse `vi`/`vd`/`vmd`
    #: triple existed only because a JSON array cannot be sparse, and there "absent" meant
    #: agreement. Stored densely the meaning INVERTS, so it is spelled out here: white on the track
    #: has to mean "identical in every genome", never "small".
    #: ⛔ **`0.0` is the MAJORITY case and it means *every genome agrees*.** The payload ships this
    #: as a sparse triple (`gaps.vi/vd/vmd`) over the gaps that vary at all, and 86.1 % of ecoli's
    #: gaps are absent from it — every one of them a measured zero. `app.js::gapVarRaw` reads
    #: absence-from-the-index as `0` and absence-of-the-whole-block as `null`, and the two are not
    #: interchangeable. Stored dense, so NULL here means only *this run did not measure variance*.
    #: ⚠ The payload's `vd` is the score × 1000; this column holds the score itself.
    length_variance_score: Mapped[float | None] = measurement()
    #: ⚠ NULL for a gap absent from the sparse index — the mode was simply not carried. So ONE
    #: absence has two different fates: the variance becomes `0.0`, the mode stays NULL.
    modal_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: Distinct named features seen in it — the coverage behind the capped list below.
    distinct_named_feature_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)


class IntergenicGapFeature(Base):
    """A named non-coding feature observed inside a gap, in some number of genomes.

    Editorial rules recorded as data, decided from the cohort census (job 34060534):
    `crispr-repeat`/`crispr-spacer` fold into the `CRISPR` array containing them at ingest;
    `regulatory_region` is kept off the track (50.7 % sit inside a CDS); `gap` rows never arrive,
    being bookkeeping dropped at the parser.
    """

    __tablename__ = "intergenic_gap_feature"
    __table_args__ = (UniqueConstraint("intergenic_gap_id", "rank_within_gap"),)

    intergenic_gap_feature_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    intergenic_gap_id: Mapped[int] = mapped_column(
        ForeignKey("intergenic_gap.intergenic_gap_id", ondelete="CASCADE"), nullable=False
    )
    rank_within_gap: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    feature_label: Mapped[str] = mapped_column(String(256), nullable=False)
    feature_type: Mapped[str] = mapped_column(String(64), nullable=False)
    observed_genome_count: Mapped[int] = mapped_column(Integer, nullable=False)
