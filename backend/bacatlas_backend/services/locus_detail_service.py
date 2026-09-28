"""One locus, whole — the hot path, and the fan-out it exists to fix.

⭐ **Rendering one locus touches a median 15–19 other loci and up to 303**, because a neighbour's
*name* is a transitive read through its Pfam and product lists. `display_name` and `best_product` are
materialised at ingest, so that read becomes `WHERE locus_id = ANY($1)` over ~20 ids: **one index
scan, in the same round trip**, carried as `neighbour_display_rows`.

⛔ **The statement count must not grow with the neighbour count.** That is the property worth
asserting, and it is what a naive implementation breaks: fetching each neighbour's name on demand is
15–19 round trips typically and 303 at worst, every one of them correct. A per-request budget of "4"
would be a number to satisfy rather than a property to hold — this module issues **one statement per
table the response draws from**, and the cost oracle asserts that against the tables, not a snapshot.

⛔ **`occ(code)` is deleted.** A slot is `null` plus an explicit `absence_reason`, never a bare `-1`:
that packed form is exactly where the "−1 means five different things" trap lives.

⚠ **The four remainders stay four.** `observed − listed` is a *display* cut; `size − observed` is
*missing data* (contig ends); `total − listed` is arrangements not drawn; `size − listed` is members
with no recorded neighbourhood. A rewrite that computes "the rest" once collapses them, and the page
starts making a claim it cannot support.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from sqlalchemy import or_, select, tuple_
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.models.intergenic_gap import IntergenicGap
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry, LocusUnirefFamilyCrosstab
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant
from bacatlas_backend.models.locus_similarity import LocusNearestLocus, LocusSimilarity
from bacatlas_backend.models.reference_vocabulary import PfamFamily
from bacatlas_backend.services.projected_genome_service import load_placements_at_locus

#: The signed offsets, in display order. ⛔ `0` is absent — it is the focal locus.
SIGNED_OFFSETS = (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)

#: How many arrangements a locus response carries before the reader must page.
#: ⚠ `total_arrangement_count` always ships separately, and the anchored genome's arrangement is
#: offered even past this cap — *"otherwise the reader is told in words that their genome sits in #37
#: and has no button to go back to it"*.
ARRANGEMENT_PAGE_SIZE = 8

#: `LocusAnnotationEntry` kinds that belong to the **function** tab and are fetched on tab open, not
#: on every walk. Splitting them is what keeps the hot path ~14 kB.
FUNCTION_ANNOTATION_KINDS = (
    AnnotationKind.COG_ORTHOGROUP,
    AnnotationKind.GENE_ONTOLOGY_SLIM,
    AnnotationKind.EC_NUMBER,
    AnnotationKind.KEGG_ORTHOLOGY,
)

#: The kinds the locus card itself shows.
CARD_ANNOTATION_KINDS = (
    AnnotationKind.GENE_SYMBOL,
    AnnotationKind.PROTEIN_PRODUCT,
    AnnotationKind.PFAM_ARCHITECTURE,
)


class LocusNotFound(LookupError):
    """No locus with that label in this species' published pangenome."""


@dataclass
class NeighbourDisplayRow:
    """Everything the track needs to DRAW a neighbour block, without a second request.

    ⚠ Each field is needed for a different part of the block: the name is its label, the genome
    count sets its alpha, and the median length sets its width (`max(6, len_nt / 10)` px).

    ⭐ `best_product` is for the **map** legend rather than the track: *"the locus NUMBER is not what
    tells you whether a neighbour belongs here — the product is"* (`app.js:2749`). The track has no
    room for it and does not ask.
    """

    locus_id: int
    #: ⛔ The row's OTHER address. Carried on the row so nothing has to look it up again in the
    #: wrong index — see `NeighbourDisplayIndex`.
    catalogue_ordinal: int
    node_label: str
    display_name: str
    display_name_source: str
    best_product: str | None
    member_genome_count: int
    median_gene_length_nt: int | None
    prevalence_band: str
    #: ⛔ The two `*_within_medoid_distance` columns and `map_position` went on 2026-09-24. The first
    #: pair fed the map's ring radius and the second the whole-catalogue sprite, both of which went
    #: with the medoid geometry they were drawn from. A neighbour row now carries only what NAMES a
    #: locus; the similarity that relates it to the focal one lives on `locus_nearest_locus`, where
    #: it belongs — it is a property of the PAIR, not of either locus.


@dataclass
class NeighbourDisplayIndex:
    """The neighbour rows, addressed BOTH ways and never in one dict.

    ⛔ `by_locus_id` and `by_catalogue_ordinal` are different key spaces over the same small
    integers. Merging them is the mistake this type exists to make unwriteable.
    """

    by_locus_id: dict = field(default_factory=dict)
    by_catalogue_ordinal: dict = field(default_factory=dict)

    def all_rows(self) -> list:
        """Every distinct row, for a serialiser that wants to emit the block once.

        ⛔ BOTH indexes. A row reached only by catalogue ordinal — every arrangement slot occupant
        that is not also a marginal mode — is a real neighbour the response refers to, and omitting
        it from `neighbour_display_rows` leaves the track with a block it can name from the slot but
        cannot label, size or colour.
        """
        return list(
            {
                row.locus_id: row
                for row in (*self.by_locus_id.values(), *self.by_catalogue_ordinal.values())
            }.values()
        )


@dataclass
class LocusDetail:
    """One locus and everything a walk step needs, in one response."""

    locus: Locus
    card_annotations: dict = field(default_factory=dict)
    uniref_families: list = field(default_factory=list)
    arrangements: list = field(default_factory=list)
    arrangements_listed: int = 0
    offset_occupants: dict = field(default_factory=dict)
    neighbour_display_rows: NeighbourDisplayIndex = field(default_factory=NeighbourDisplayIndex)
    intergenic_gaps: list = field(default_factory=list)
    #: `{representation: LocusSimilarity}` — the eight numbers the card's three views show.
    similarity: dict = field(default_factory=dict)
    #: `{representation: [LocusNearestLocus, …]}` in rank order, 1-based, ragged.
    nearest_loci: dict = field(default_factory=dict)
    #: ⭐ How many DISTINCT other loci this response resolved. The fan-out, measured per request.
    resolved_neighbour_count: int = 0
    #: ⭐ Every Pfam family this response MENTIONS, resolved once — `{accession: PfamFamily}`. The
    #: same move as `neighbour_display_rows`: a bounded block in one round trip instead of the page
    #: carrying an 833 kB vendored reference to render a chip. ⚠ Keyed VERSION-STRIPPED.
    pfam_families: dict = field(default_factory=dict)
    #: ⛔ Which arrangement RANKS the anchored genome carries here — a **list**, because a genome at
    #: rho > 1 occupies two arrangements at one locus and there is no uniqueness constraint on
    #: (locus, genome) anywhere. Empty when there is no anchor, or when the anchored genome has no
    #: gene at this locus; those two are distinguished by `anchor_genome_id` being set at all.
    anchor_arrangement_ranks: list[int] = field(default_factory=list)
    #: Whether a genome was anchored at all. ⚠ Distinguishes "no anchor set" from "anchored, and
    #: this genome has no gene at this locus" — an empty rank list means the second only when this
    #: is true, and the two are different sentences on the page.
    is_anchored: bool = False
    #: ⛔ A PROJECTED genome's genes at this locus — placed after the model was built, never members.
    #: ⚠ `anchor_kind` is what stops the page saying "this genome is in arrangement #3" about a
    #: genome that is in no arrangement at all: a projected genome MATCHES an arrangement, which is
    #: a different relation from being counted in one, and the two must not share a sentence.
    anchor_kind: str | None = None
    projected_placements: list = field(default_factory=list)


def pfam_accessions_in(architecture: str | None) -> list[str]:
    """`"PF00126.29,PF03466"` → `["PF00126", "PF03466"]`.

    ⛔ **The version suffix must be cut before lookup.** `pfam_reference` strips it when it builds
    the table, so an annotation carrying `PF00126.29` looked up as written silently misses and the
    chip falls back to a bare accession — which reads as "this family has no name" rather than as a
    failed join. The published page cuts at the dot for exactly this reason (`app.js:3285`).
    """
    if not architecture:
        return []
    out = []
    for raw in architecture.split(","):
        accession = raw.split(".")[0].strip()
        if accession:
            out.append(accession)
    return out


def _resolve_pfam_families(session: Session, detail: LocusDetail) -> dict:
    """Every Pfam accession this response mentions, in ONE statement.

    ⚠ Bounded by construction: the top-5 architectures plus one modal architecture per listed
    UniRef50 family, each a handful of domains. It does not grow with the catalogue, which is what
    makes a single `= ANY(...)` the right shape rather than a per-chip lookup.
    """
    wanted: set[str] = set()
    for entry in detail.card_annotations.get(AnnotationKind.PFAM_ARCHITECTURE.value, []):
        wanted.update(pfam_accessions_in(entry.term_value))
    for family in detail.uniref_families:
        wanted.update(pfam_accessions_in(family.modal_pfam_architecture))
    if not wanted:
        # ⛔ No statement at all rather than `IN ()`. An empty IN is valid SQL and returns nothing,
        # but it still costs a round trip on every locus with no Pfam coverage — which is 22 % of
        # them — and the cost oracle would then be measuring a query that can never return a row.
        return {}
    return {
        family.pfam_accession: family
        for family in session.execute(
            select(PfamFamily).where(PfamFamily.pfam_accession.in_(sorted(wanted)))
        ).scalars()
    }


def _locus_by_label(session: Session, pangenome_id: int, node_label: str) -> Locus:
    locus = session.execute(
        select(Locus).where(Locus.pangenome_id == pangenome_id, Locus.node_label == node_label)
    ).scalar_one_or_none()
    if locus is None:
        raise LocusNotFound(
            f"no locus {node_label!r} in pangenome {pangenome_id}. ⚠ Node labels are TEXT — "
            "`0123` and `123` are different loci, and a numeric round-trip loses the distinction."
        )
    return locus


def listed_arrangement_condition(
    arrangement_limit: int,
    anchor_genome_id: int | None,
    matched_arrangement_ids: Sequence[int] | None = None,
):
    """Which arrangements a response carries: the commonest few, OR any the anchored genome sits in.

    ⭐ `member_genome_ids @> ARRAY[?]` is the GIN index doing the work a binary search over a
    ~490k-entry Int32Array did in the browser. ⛔ **The OR is the rule, not an optimisation** — the
    anchored arrangement is offered however rare it is, which is `arrShown`'s rule: *"otherwise the
    reader is told in words that their genome sits in #37 and has no button to go back to it"*.
    Measured on the probe catalogues, 15,643 (ecoli) and 10,322 (kp) (locus, genome) pairs carry an
    arrangement past the cap of 8 — ranks as deep as #84 — and without this clause every one of them
    is a genome anchored to a neighbourhood the response does not contain. Parity suite T4 drives
    this statement for every one of them.
    """
    condition = LocusArrangement.rank_within_locus < arrangement_limit
    if anchor_genome_id is not None:
        condition = or_(condition, LocusArrangement.member_genome_ids.contains([anchor_genome_id]))
    # ⭐ The same rule for a PROJECTED genome, and it has to be a different clause: a projected
    # genome is in no `member_genome_ids`, so the array test above can never find it. Without this
    # a reader is told their genome's neighbourhood is arrangement #37 and the response does not
    # contain #37.
    if matched_arrangement_ids:
        condition = or_(condition, LocusArrangement.locus_arrangement_id.in_(list(matched_arrangement_ids)))
    return condition


def load_listed_arrangements(
    session: Session,
    *,
    locus_id: int,
    arrangement_limit: int,
    anchor_genome_id: int | None,
    matched_arrangement_ids: Sequence[int] | None = None,
) -> list[LocusArrangement]:
    """The arrangements one locus response carries, in rank order. One statement.

    Separate from `load_locus_detail` so parity suite T4 can drive THIS statement for every anchored
    pair past the cap without paying for the rest of the view each time.
    """
    return list(
        session.execute(
            select(LocusArrangement)
            .where(
                LocusArrangement.locus_id == locus_id,
                listed_arrangement_condition(
                    arrangement_limit, anchor_genome_id, matched_arrangement_ids
                ),
            )
            .order_by(LocusArrangement.rank_within_locus)
        ).scalars()
    )


def anchor_arrangement_ranks(arrangements, anchor_genome_id: int) -> list[int]:
    """The ranks, among `arrangements`, that the anchored genome sits in — ascending, and a LIST.

    ⛔ **A list, never the first match.** A genome at ρ > 1 has two genes at one locus and can sit in
    two of its arrangements (1,032 ecoli / 801 kp (locus, genome) pairs; up to 14 at one locus), and
    there is no uniqueness constraint on (locus, genome) anywhere. The published page's
    `anchorRanks` returns every one, marks each with ⚓, and opens on the lowest.

    ⛔ **Recomputed from the rows, never inferred from the query that added them.** The anchored
    arrangement is often ALSO within the cap, in which case the OR in `listed_arrangement_condition`
    added nothing, and a client assuming "the extra row is the anchored one" would mark the wrong one.

    ⚠ Ranks are `rank_within_locus`, 0-based — the same integer as the page's index into its
    arrangement list, which it prints as `#rank + 1`.
    """
    return [
        arrangement.rank_within_locus
        for arrangement in arrangements
        if anchor_genome_id in (arrangement.member_genome_ids or ())
    ]


def membership_is_complete(locus: Locus) -> bool:
    """Whether EVERY genome present at this locus reaches some arrangement — `app.js::membershipComplete`.

    ⛔ The only thing that settles which of the anchor line's two sentences is true when the anchored
    genome carries no arrangement here: *"has no gene at this locus"* (complete) or *"has no recorded
    neighbourhood at this locus"* (not). At 6.26 % of ecoli loci the first would be false.

    ⚠ A GENOME question, so it compares the genome union (`arrangement_member_genome_count`) and never
    the gene counts: a genome at ρ > 1 can lose one gene's window and keep its arrangement.
    """
    return locus.arrangement_member_genome_count >= locus.member_genome_count


def load_locus_detail(
    session: Session,
    *,
    pangenome_id: int,
    node_label: str,
    anchor_genome_id: int | None = None,
    projected_genome_id: int | None = None,
    arrangement_limit: int = ARRANGEMENT_PAGE_SIZE,
) -> LocusDetail:
    """The whole locus view. One statement per table, and none of them per neighbour.

    ⚠ `anchor_genome_id` and `projected_genome_id` are mutually exclusive — one genome is anchored
    at a time, and the two are different relations to this pangenome. The route refuses both.
    """
    locus = _locus_by_label(session, pangenome_id, node_label)
    detail = LocusDetail(locus=locus)

    # ⭐ Resolved BEFORE the arrangements, because the arrangement a projected gene matched has to be
    # listed even when it sits past the display cap — the same rule the anchor has, through a
    # different clause.
    matched_arrangement_ids: list[int] = []
    if projected_genome_id is not None:
        detail.projected_placements = load_placements_at_locus(
            session,
            pangenome_id=pangenome_id,
            genome_id=projected_genome_id,
            locus_id=locus.locus_id,
        )
        matched_arrangement_ids = [
            placement.matched_locus_arrangement_id
            for placement in detail.projected_placements
            if placement.matched_locus_arrangement_id is not None
        ]

    # ── the card's own lists ───────────────────────────────────────────────────────────────────
    entries = session.execute(
        select(LocusAnnotationEntry)
        .where(
            LocusAnnotationEntry.locus_id == locus.locus_id,
            LocusAnnotationEntry.annotation_kind.in_(CARD_ANNOTATION_KINDS),
        )
        .order_by(LocusAnnotationEntry.annotation_kind, LocusAnnotationEntry.rank_within_locus)
    ).scalars()
    for entry in entries:
        detail.card_annotations.setdefault(entry.annotation_kind.value, []).append(entry)

    detail.uniref_families = list(
        session.execute(
            select(LocusUnirefFamilyCrosstab)
            .where(LocusUnirefFamilyCrosstab.locus_id == locus.locus_id)
            .order_by(LocusUnirefFamilyCrosstab.rank_within_locus)
        ).scalars()
    )

    # ── the joint view, capped, plus whichever the anchored genome carries ─────────────────────
    detail.arrangements = load_listed_arrangements(
        session,
        locus_id=locus.locus_id,
        arrangement_limit=arrangement_limit,
        anchor_genome_id=anchor_genome_id,
        matched_arrangement_ids=matched_arrangement_ids,
    )
    detail.arrangements_listed = len(detail.arrangements)
    detail.is_anchored = anchor_genome_id is not None or projected_genome_id is not None
    detail.anchor_kind = (
        "catalogue" if anchor_genome_id is not None else ("projected" if projected_genome_id is not None else None)
    )
    if anchor_genome_id is not None:
        detail.anchor_arrangement_ranks = anchor_arrangement_ranks(
            detail.arrangements, anchor_genome_id
        )

    detail.pfam_families = _resolve_pfam_families(session, detail)

    # ── the marginal view ──────────────────────────────────────────────────────────────────────
    occupants = session.execute(
        select(LocusOffsetOccupant)
        .where(LocusOffsetOccupant.locus_id == locus.locus_id)
        .order_by(LocusOffsetOccupant.signed_offset, LocusOffsetOccupant.rank_within_offset)
    ).scalars()
    for occupant in occupants:
        detail.offset_occupants.setdefault(occupant.signed_offset, []).append(occupant)

    # ── the set-to-set similarity, both representations ────────────────────────────────────────
    # ⚠ Loaded BEFORE the fan-out, deliberately: the card lists the five nearest loci in each
    # representation, and those are a **different set** from the track's ±5 neighbours — the two
    # representations do not even agree with each other (their separations correlate at rho ~0.47).
    # Resolving them in the same statement costs nothing; a second one would be an N+1 the cost
    # oracle is there to refuse.
    for similarity in session.execute(
        select(LocusSimilarity).where(LocusSimilarity.locus_id == locus.locus_id)
    ).scalars():
        detail.similarity[similarity.representation.value] = similarity

    # ⭐ One statement for BOTH representations' lists, ordered so the card never sorts. The list is
    # RAGGED — a locus whose shortlist held fewer than five others simply has fewer rows — so the
    # card reads what is there rather than padding to five.
    for nearest in session.execute(
        select(LocusNearestLocus)
        .where(LocusNearestLocus.locus_id == locus.locus_id)
        .order_by(LocusNearestLocus.representation, LocusNearestLocus.rank)
    ).scalars():
        detail.nearest_loci.setdefault(nearest.representation.value, []).append(nearest)

    # ── the fan-out, resolved in ONE statement ─────────────────────────────────────────────────
    # ⛔⛔ **TWO DIFFERENT INTEGER SPACES, KEPT APART.** `locus_offset_occupant.neighbour_locus_id`
    # is a surrogate `locus_id`; an arrangement slot code carries a **catalogue ordinal**
    # (`code // 2`). Both are small integers over the same range, so locus_id 5 and ordinal 5 are
    # different loci and a single lookup keyed by "the number" silently draws one on the other —
    # with a page that still looks entirely correct. They are collected separately and returned
    # separately, so a caller has to say which address it is holding.
    neighbour_locus_ids = {
        occupant.neighbour_locus_id
        for rows in detail.offset_occupants.values()
        for occupant in rows
    }
    neighbour_ordinals = {
        code // 2
        for arrangement in detail.arrangements
        for code in arrangement.neighbour_slot_codes
        if code >= 0
    }
    # ⭐ The card's five nearest, per representation — real `locus_id`s now, so they join the
    # locus-id space above rather than the ordinal one.
    # ⚠ This is why the retired `-1` sentinel is gone: a neighbour outside the catalogue has no row
    # at all, instead of an ordinal that had to be filtered before it named locus 0.
    neighbour_locus_ids.update(
        nearest.neighbour_locus_id
        for rows in detail.nearest_loci.values()
        for nearest in rows
    )
    detail.neighbour_display_rows = _neighbour_display_rows(
        session,
        pangenome_id=pangenome_id,
        locus_ids=neighbour_locus_ids,
        catalogue_ordinals=neighbour_ordinals,
        focal_locus=locus,
    )
    # ⚠ Every DISTINCT locus the block resolved, across BOTH key spaces — not just the marginal
    # occupants. It counted only `by_locus_id` while the block already carried arrangement occupants
    # reached by ordinal, so the number under-reported the fan-out it exists to measure; adding the
    # map's nearest loci made that visible rather than causing it.
    detail.resolved_neighbour_count = len(detail.neighbour_display_rows.all_rows())

    # ── the gaps the drawn track needs — every adjacent pair, not just the focal's two ─────────
    # ⛔⛔ **Every adjacent pair in every LISTED arrangement.** This read only the gaps the focal
    # locus is an endpoint of, so the track drew the two regions beside the focal gene and nothing
    # else — nine genes packed edge to edge with nothing between them, which reads as *"these genes
    # are adjacent"*: a claim the data does not make, on every locus, found only by LOOKING at the
    # rebuilt page beside the published one. Still one statement: the pairs are known from the slot
    # codes and the ordinals the fan-out just resolved.
    # ⛔ Both orders of each pair, because the canonical order is by node_label and a caller holding
    # two locus_ids cannot tell which side each is on without resolving the labels.
    detail.intergenic_gaps = _load_track_gaps(session, pangenome_id=pangenome_id, detail=detail)

    return detail


def _adjacent_locus_id_pairs(detail: LocusDetail) -> set[tuple[int, int]]:
    """Every pair of loci drawn side by side in some listed arrangement — both orders.

    The ten slot codes are the ±5 window in recorded order, `-5 … -1, +1 … +5`, with the focal locus
    between them. ⚠ A code of `-1` is *the contig ends here* — a real observation, and it BREAKS the
    adjacency rather than being skipped over: the genes either side of a contig end are not
    neighbours, and pairing them would ask for a region that does not exist.
    """
    locus_id_by_ordinal = {
        ordinal: row.locus_id for ordinal, row in detail.neighbour_display_rows.by_catalogue_ordinal.items()
    }
    locus_id_by_ordinal[detail.locus.catalogue_ordinal] = detail.locus.locus_id
    focal_code = detail.locus.catalogue_ordinal * 2
    pairs: set[tuple[int, int]] = set()
    for arrangement in detail.arrangements:
        codes = list(arrangement.neighbour_slot_codes)
        window = codes[:5] + [focal_code] + codes[5:]
        for left, right in zip(window, window[1:], strict=False):
            if left < 0 or right < 0:
                continue
            a = locus_id_by_ordinal.get(left // 2)
            b = locus_id_by_ordinal.get(right // 2)
            if a is None or b is None:
                continue
            pairs.add((a, b))
            pairs.add((b, a))
    return pairs


def _load_track_gaps(session: Session, *, pangenome_id: int, detail: LocusDetail) -> list:
    """The regions between every adjacent DRAWN pair — one statement, and nothing else.

    ⚠ Only drawn pairs. The old query also returned every gap the focal locus is an endpoint of,
    including gaps to neighbours that sit only in arrangements past the display cap — loci this
    response never resolves, so those gaps arrived with a `null` label, which a client cannot key
    and silently drops. A gap is useful here only between two loci the track can put side by side.
    And with no arrangement drawn there are no gaps at all: two marginal modes are not neighbours
    in any genome, so the region between them was never observed (`app.js:987`).
    """
    pairs = _adjacent_locus_id_pairs(detail)
    if not pairs:
        return []
    return list(
        session.execute(
            select(IntergenicGap).where(
                IntergenicGap.pangenome_id == pangenome_id,
                tuple_(IntergenicGap.flanking_locus_id_a, IntergenicGap.flanking_locus_id_b).in_(sorted(pairs)),
            )
        ).scalars()
    )


def _neighbour_display_rows(
    session: Session,
    *,
    pangenome_id: int,
    locus_ids: set[int],
    catalogue_ordinals: set[int],
    focal_locus: Locus,
) -> NeighbourDisplayIndex:
    """Every locus this response refers to, by either address — **one statement, always**."""
    index = NeighbourDisplayIndex()
    if not locus_ids and not catalogue_ordinals:
        return index
    # ⭐ **One statement, and now a plain one.** It used to carry two OUTER joins onto
    # `locus_embedding_geometry` for each neighbour's position on the catalogue sprite; the sprite
    # went on 2026-09-24 and the joins with it. A neighbour row is what NAMES a locus.
    rows = session.execute(
        select(
            Locus.locus_id,
            Locus.catalogue_ordinal,
            Locus.node_label,
            Locus.display_name,
            Locus.display_name_source,
            Locus.best_product,
            Locus.member_genome_count,
            Locus.median_gene_length_nt,
            Locus.prevalence_band,
        )
        .where(
            Locus.pangenome_id == pangenome_id,
            # ⛔⛔ `catalogue_ordinals` on the second branch, NOT `locus_ids`. This line read
            # `catalogue_ordinal.in_(locus_ids)` — the precise mistake the comment at the call site
            # warns about, made one screen below it: two key spaces over the same small integers,
            # merged.
            #
            # ⚠ **It mostly worked, which is why nothing caught it.** An arrangement occupant is
            # usually also a marginal mode, so its row came back on the FIRST branch and was then
            # indexed by its ordinal anyway. Only the occupants that are not marginal modes fell
            # through: measured at **36 of 5,423 named slots (0.7 %) over a random 200 loci, but
            # touching 11 of those 200** — and each one is a blank, unwalkable block sitting where a
            # real neighbour is. On loci with unusually many arrangements it is far worse: on the
            # three the API-contract fixture pins, 37 of 190.
            or_(Locus.locus_id.in_(locus_ids), Locus.catalogue_ordinal.in_(catalogue_ordinals)),
        )
    ).all()
    for locus_id, ordinal, label, name, source, product, genomes, length, band in rows:
        row = NeighbourDisplayRow(
            locus_id=locus_id,
            catalogue_ordinal=ordinal,
            node_label=label,
            display_name=name,
            display_name_source=source,
            best_product=product,
            member_genome_count=genomes,
            median_gene_length_nt=length,
            prevalence_band=band.value,
        )
        if locus_id in locus_ids:
            index.by_locus_id[locus_id] = row
        if ordinal in catalogue_ordinals:
            index.by_catalogue_ordinal[ordinal] = row
    # ⚠ The focal locus is its own neighbour on a tandem repeat (hokE, ldrB, zorO …) — a real case
    # and not a degenerate one, so it is always resolvable from its own response.
    focal = NeighbourDisplayRow(
        locus_id=focal_locus.locus_id,
        catalogue_ordinal=focal_locus.catalogue_ordinal,
        node_label=focal_locus.node_label,
        display_name=focal_locus.display_name,
        display_name_source=focal_locus.display_name_source,
        best_product=focal_locus.best_product,
        member_genome_count=focal_locus.member_genome_count,
        median_gene_length_nt=focal_locus.median_gene_length_nt,
        prevalence_band=focal_locus.prevalence_band.value,
    )
    index.by_locus_id.setdefault(focal_locus.locus_id, focal)
    index.by_catalogue_ordinal.setdefault(focal_locus.catalogue_ordinal, focal)
    return index


def load_function_block(session: Session, *, locus_id: int) -> dict:
    """The EggNOG tab — fetched on tab open, not on every walk. One statement."""
    entries = session.execute(
        select(LocusAnnotationEntry)
        .where(
            LocusAnnotationEntry.locus_id == locus_id,
            LocusAnnotationEntry.annotation_kind.in_(FUNCTION_ANNOTATION_KINDS),
        )
        .order_by(
            LocusAnnotationEntry.annotation_kind,
            LocusAnnotationEntry.gene_ontology_namespace,
            LocusAnnotationEntry.rank_within_locus,
        )
    ).scalars()
    grouped: dict[str, list] = {}
    for entry in entries:
        grouped.setdefault(entry.annotation_kind.value, []).append(entry)
    return grouped


def load_arrangement_page(
    session: Session, *, locus_id: int, offset: int = 0, limit: int = 50
) -> list[LocusArrangement]:
    """Arrangements past the display cut — the full scroller, paged. One statement."""
    return list(
        session.execute(
            select(LocusArrangement)
            .where(LocusArrangement.locus_id == locus_id)
            .order_by(LocusArrangement.rank_within_locus)
            .offset(offset)
            .limit(limit)
        ).scalars()
    )
