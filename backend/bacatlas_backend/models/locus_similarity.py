"""Set-to-set embedding similarity: the per-locus numbers, the five nearest, and the baseline.

⛔ **What this replaced, and why it is a different measurement rather than a renaming.**
``locus.{esm,bacformer}_{within,nearest}_medoid_distance`` measured a locus by reducing it to ONE
member — the gene with the highest cosine to its centroid — and then measuring that single point:
``within`` was its members' distance to that gene, ``nearest`` that gene's distance to another
locus's medoid. ``bac80_marker_tokens.md`` §4.4 measured that construction inverting the
inter/intra ordering; on the probe catalogue it did not invert but is biased optimistic — the
all-pairs median sits 0.009–0.017 below the medoid one and separation < 0 was understated at 2.5 %
against the true 4.2 % (rank agreement ρ 0.936–0.988, so the *ordering* was largely right).

Everything here is measured over the whole **set**, or anchored on a **point**. The three views the
card offers are NOT three pictures of one thing: flagging by each — worst decile of separation,
negative weak margin, own fraction below 1 — only **34 %** (ecoli/bacformer) and **18 %** (esm) of
flagged loci are flagged by all three, and of the median's worst decile only 50 % / 22 % have a
negative margin. Rank ρ 0.79–0.93: they agree overall and part at the tail, which is where flagging
happens. That measurement is why all three are stored rather than one being chosen.

⛔ **Every name here is kept under the length that makes a constraint name TRUNCATE.** Postgres caps an
identifier at 63 characters; past that SQLAlchemy substitutes a 4-hex hash, so the database gets a name
nobody wrote — and Alembic's autogenerate compares check constraints against the name that was *asked* for,
so such a table reports as drifted from its own migration **forever**, with no fix short of renaming. The
binding constraint is ``ck_{table}__{column}_is_not_nan``: with this table's name that leaves 31 characters
for a column, and the longest here uses 22.

⚠ **``weakest → own`` implies the own fraction by ARITHMETIC**: if the weakest member's nearest gene
is a stranger then the fraction cannot be 1. Those two agreeing is never evidence of anything.
"""

from __future__ import annotations

from sqlalchemy import (
    BigInteger,
    ForeignKey,
    Index,
    Integer,
    SmallInteger,
    String,
    UniqueConstraint,
)
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.column_types import measurement, nan_guards
from bacatlas_backend.models.enumerations import EmbeddingRepresentation

#: How many other loci the card lists per representation. The artifact is RAGGED — a locus whose
#: shortlist held fewer keeps fewer rows — so this is a maximum, never a row count to assert.
NEAREST_LOCUS_COUNT = 5


class PangenomeSimilarityBaseline(Base):
    """One (pangenome, representation) baseline — what makes a cosine on this card mean anything.

    ⛔ **This is NOT the retired ``locus_map_projection`` null.** That sampled random pairs of
    **medoids**; every number on this card is a median over **gene** pairs. On the published *E.
    coli* catalogue the two sit at **0.0651** and **0.0587** — close enough to look interchangeable
    and not be. Sampled over ~1 M random gene pairs by ``build_cluster_similarity``.
    """

    __tablename__ = "pangenome_similarity_baseline"
    __table_args__ = (
        UniqueConstraint("pangenome_id", "representation"),
        *nan_guards(
            "floor_median",
            "floor_p25",
            "floor_p75",
            "floor_p99",
        ),
    )

    pangenome_similarity_baseline_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )
    representation: Mapped[EmbeddingRepresentation] = mapped_column(nullable=False)

    #: ``raw`` or ``centred``. ⛔ Raw is what is published (David, 2026-09-24); the centred path stays
    #: behind a flag so the 2026-09-24 numbers reproduce. Stored because a centred number captioned
    #: as a raw one is indistinguishable from the real thing — ESM's floor moves 0.7417 → ~0.005.
    similarity_form: Mapped[str] = mapped_column(String(16), nullable=False, default="raw")

    #: The median cosine over the random gene-pair sample, and its spread. A median alone cannot say
    #: whether a locus's 0.41 sits far outside random or inside its shoulder, so the card draws a
    #: p25–p75 box with a whisker to p99 rather than a bare tick.
    floor_median: Mapped[float | None] = measurement()
    floor_p25: Mapped[float | None] = measurement()
    floor_p75: Mapped[float | None] = measurement()
    floor_p99: Mapped[float | None] = measurement()

    #: ⭐ **The other half of "p12 of 12,104 loci".** The percentiles on ``locus_similarity``
    #: are MIDRANKS over the loci where the quantity is measurable — not over the catalogue — and the
    #: card prints both halves of that sentence. Both are assertions, so the denominator needs a home
    #: at exactly this grain.
    #: ⚠ The measurable set is defined by the MEASUREMENT, not by the prevalence band: on the
    #: published ecoli catalogue it is exactly the 5,427 loci of size 1 that are excluded, while
    #: ``prevalence_band = RARE`` covers 5,458 — the extra 31 are paralogues inside one genome, which
    #: DO have a within-locus pair. Gating on the band would blank 31 real measurements.
    measurable_locus_count: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: The kNN width the shortlist and the point measures were computed at (200). Above the largest
    #: locus in either species (143 ecoli / 159 kp), so no locus can fill its members' neighbour
    #: lists with itself and come back with an empty shortlist.
    neighbour_knn_k: Mapped[int | None] = mapped_column(Integer, nullable=True)
    nearest_locus_count: Mapped[int] = mapped_column(
        SmallInteger, nullable=False, default=NEAREST_LOCUS_COUNT
    )
    source_csv_path: Mapped[str | None] = mapped_column(String(512), nullable=True)


class LocusSimilarity(Base):
    """One locus's set-to-set similarity in one representation — the eight numbers the card shows.

    Three pairs and two shares, in the order the card's three views offer them. ⛔ The two
    DIFFERENCES the views print — ``separation`` and ``margin`` — are deliberately **not** stored.
    Each is the difference of two columns that are, the API subtracts, and a separately-rounded
    difference drifting from the two rounded numbers printed beside it is the exact class of quiet
    disagreement this card was rebuilt to end.
    """

    __tablename__ = "locus_similarity"
    __table_args__ = (
        UniqueConstraint("locus_id", "representation"),
        *nan_guards(
            "within_similarity",
            "nearest_similarity",
            "weak_own_similarity",
            "weak_other_similarity",
            "own_neighbour_fraction",
            "separation_percentile",
            "weak_margin_percentile",
            "own_fraction_percentile",
        ),
    )

    locus_similarity_id: Mapped[int] = mapped_column(
        BigInteger, primary_key=True, autoincrement=True
    )
    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)
    representation: Mapped[EmbeddingRepresentation] = mapped_column(nullable=False)

    #: ⚠ NULL on a SINGLETON, which has no pair inside its locus — 5,427 of *E. coli*'s 17,531 loci.
    #: A 0.0 here would claim its members are unrelated to each other.
    within_similarity: Mapped[float | None] = measurement()
    #: The highest such median against another locus. ⚠ Present for a singleton: that question is
    #: well posed for one gene, which is why it is not NULL alongside the others.
    nearest_similarity: Mapped[float | None] = measurement()

    #: The member least attached to its own locus: its nearest neighbour INSIDE the locus, and then
    #: that SAME gene's nearest gene anywhere else. ⭐ The point is invariant where a median moves
    #: with the set — and the clustering joins points, not medians (David, 2026-09-24).
    weak_own_similarity: Mapped[float | None] = measurement()
    weak_other_similarity: Mapped[float | None] = measurement()

    #: The share of members whose nearest gene in the whole species is another member of this locus.
    #: ⚠ A share, NOT a cosine — it must never be drawn against the gene-pair baseline.
    own_neighbour_fraction: Mapped[float | None] = measurement()

    #: ⭐ Precomputed MIDRANKS over **measurable loci only**, one per view. ⚠ Singletons are NULL and
    #: must read *not measurable*, never ``0.000``.
    #: ⛔ ``own_fraction_percentile`` is a midrank inside a tie block covering **88.5 %** of the
    #: catalogue, so at exactly 1.0 it reads "p53" — *better than half the catalogue* — when it means
    #: *tied with nearly all of it*. Stored because below 1.0 the whole block is above and it means
    #: what it looks like; the card is what must decline to print it at 1.0.
    separation_percentile: Mapped[float | None] = measurement()
    weak_margin_percentile: Mapped[float | None] = measurement()
    own_fraction_percentile: Mapped[float | None] = measurement()


class LocusNearestLocus(Base):
    """The five nearest OTHER loci, per representation — ranked, and a row per neighbour.

    ⭐ Rows rather than an ordinal array (which is what the retired ``nearest_locus_ordinals`` was),
    because each carries its own similarity and the list is **ragged**: a locus whose shortlist held
    fewer than five others simply has fewer rows. An array would have needed a sentinel, and a
    sentinel read as an index names a real locus that looks entirely plausible.

    ⚠ ``rank`` is 1-based, as the artifact writes it, and rank 1 is BY CONSTRUCTION the locus whose
    similarity is ``locus_similarity.nearest_similarity``. The ingest asserts that
    rather than trusting it: a disagreement can only mean two different runs.
    """

    __tablename__ = "locus_nearest_locus"
    __table_args__ = (
        UniqueConstraint("locus_id", "representation", "rank"),
        # ⚠ The unique constraint indexes `locus_id`; `neighbour_locus_id` is a SECOND FK to the same
        # table and has no index without this one. Deleting a pangenome would then scan every row
        # per locus — the omission that cost an 8-minute stall on `gene_locus_membership`.
        Index("ix_locus_nearest_locus__neighbour_locus_id", "neighbour_locus_id"),
        *nan_guards("cross_similarity"),
    )

    locus_nearest_locus_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)
    representation: Mapped[EmbeddingRepresentation] = mapped_column(nullable=False)
    rank: Mapped[int] = mapped_column(SmallInteger, nullable=False)

    #: ⚠ A real FK, not a catalogue ordinal. The retired array stored ordinals and every reader had
    #: to resolve them; a neighbour outside the catalogue simply has no row here.
    neighbour_locus_id: Mapped[int] = mapped_column(
        ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False
    )
    #: The cosine over every gene pair spanning the two loci: a median on a full-audit build, the exact
    #: mean on a lean one (`pair_estimator` says which).
    cross_similarity: Mapped[float | None] = measurement()


#: ⛔ The clause the ingest writes into `pangenome.omitted_sections["similarity_estimator[<rep>]"]` when a
#: similarity artifact holds MEANS — nuna's lean producer, `rank_by=mean_cross_pair_cosine`. It is the one
#: stored fact that says which estimator a catalogue's within- and nearest-locus numbers are, so the writer
#: (`catalogue_frames._similarity_estimator_note`) and the reader (`pair_estimator`) share it here: the api
#: package imports nothing from ingest.
LEAN_MEAN_ESTIMATOR_CLAUSE = "its within-locus and nearest-locus numbers are MEANS on this build"


def pair_estimator(omitted_sections: dict | None, representation: str) -> str:
    """`"mean"` or `"median"` — what one representation's within- and nearest-locus numbers ARE.

    A full-audit artifact (`rank_by=sim`) holds a median over every enumerated pair and writes no estimator
    note. nuna's lean producer holds the exact mean cosine over every pair of embedded members, computed from
    each locus's summed unit vectors, and the ingest says so in a note the footer renders — this reads it.
    ⚠ The note's KEY alone is not the test: it is also written for a median artifact whose genes are not all
    embedded. The clause is.
    """
    note = (omitted_sections or {}).get(f"similarity_estimator[{representation}]") or ""
    return "mean" if LEAN_MEAN_ESTIMATOR_CLAUSE in note else "median"
