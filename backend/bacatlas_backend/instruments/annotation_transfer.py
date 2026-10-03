"""Inferring a function for a node that has none — the ladder, and the rate the reader judges it by.

Two mechanisms, and the page must never blur them (`nuna/docs/model_evaluation/COG_function_inference.md`):

* **internal inference** — some of the node's own genes carry a call, so the modal call covers the
  rest. Measured at 99.5 % / 99.9 % within-node unanimity; this is where the lift is.
  ⛔ A node with exactly ONE annotated gene is unanimous **by construction** — the absence of
  evidence, not evidence — so `own_support.checkable` is False and the page says so.
* **neighbour transfer** — no member gene carries a call, so the nearest ESM neighbour that does
  supplies one, at a depth set by how similar it is.

⭐ **The agreement rate is COMPUTED from the loaded catalogue, never hard-coded.** It is the thing a
reader judges a suggestion by, so a literal lifted out of the write-up would keep reading plausibly
long after the catalogue moved underneath it — the exact failure recorded in nuna `PROJECT_STATE.md`
§6, 2026-10-03 ("a frozen oracle cannot detect a drift it shares"). This module is also what
`scripts/measure_cog_function_inference.py` reports from, so the page and the document cannot drift
apart: there is one implementation of the measurement, not two.

⚠ **`chance` is EXACT here, not sampled.** For a focal node's claim set it is the share of annotated
nodes whose own set intersects it — computed from an inverted index, memoised on the claim set
(there are ~118 distinct COG category values in a catalogue, so the memo collapses the work). An
earlier sampled version drew 200 random partners per pair, which was both slower and gave the page a
number that differed in the third digit from the document's.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import cache

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import AnnotationKind, EmbeddingRepresentation
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry
from bacatlas_backend.models.locus_similarity import LocusNearestLocus

#: ⛔ **The ladder is a DECISION (David, 2026-10-03), not a measurement** — recorded in nuna
#: `PROJECT_STATE.md` §6. Half-open `[low, high)`, ordered strongest first. Below 0.96 is ONE tier
#: because 0.95-0.96 and 0.90-0.95 were statistically indistinguishable, and three indistinguishable
#: claims read worse than one; below 0.90 nothing is called at all ("definitely too remote to call").
TIERS: tuple[tuple[str, float, float], ...] = (
    (">= 0.99", 0.99, 1.01),
    ("0.98-0.99", 0.98, 0.99),
    ("0.97-0.98", 0.97, 0.98),
    ("0.96-0.97", 0.96, 0.97),
    ("0.90-0.96", 0.90, 0.96),
)
#: Not a tier: a cosine this low is shown as a neighbour and carries NO suggested function.
NOT_CALLED = "< 0.90"

#: The deepest level each tier may quote. EC levels are field counts (4 = `2.7.10.1`, 1 = `2.-.-.-`);
#: COG's two rungs are 2 = the orthogroup accession, 1 = its functional-category letters.
QUOTED_LEVEL: dict[AnnotationKind, dict[str, int]] = {
    AnnotationKind.EC_NUMBER: {">= 0.99": 4, "0.98-0.99": 3, "0.97-0.98": 1,
                              "0.96-0.97": 1, "0.90-0.96": 1},
    AnnotationKind.COG_ORTHOGROUP: {">= 0.99": 2, "0.98-0.99": 2, "0.97-0.98": 1,
                                   "0.96-0.97": 1, "0.90-0.96": 1},
}
LEVEL_LABEL: dict[AnnotationKind, dict[int, str]] = {
    AnnotationKind.EC_NUMBER: {4: "full EC code", 3: "EC sub-subclass", 2: "EC subclass",
                              1: "EC class"},
    AnnotationKind.COG_ORTHOGROUP: {2: "COG orthogroup", 1: "COG category"},
}
#: Below this a cell reports its pair count and NO rate. A rate off 7 pairs is not a rate.
MIN_PAIRS = 30
SUPPORTED_KINDS = tuple(QUOTED_LEVEL)


def tier_for(cosine: float | None) -> str:
    """Which tier a neighbour's cosine falls in, or `NOT_CALLED`."""
    if cosine is None:
        return NOT_CALLED
    for name, low, high in TIERS:
        if low <= cosine < high:
            return name
    return NOT_CALLED


def ec_levels(term_value: str) -> dict[int, frozenset[str]]:
    """Level → the prefixes this EC value RESOLVES to; empty means it does not state that level.

    ⚠ Two shapes a `split_part` would get wrong, both measured in the published catalogues:
    1,034 of 9,398 values are comma-joined SETS of codes (`1.6.5.9,7.1.1.-`), so splitting the whole
    string on `.` yields `9,7` — a code neither side holds; and 2,257 carry `-` placeholders, so a
    pair can be **incomparable** at a level without disagreeing about anything.
    """
    found: dict[int, set[str]] = {level: set() for level in (1, 2, 3, 4)}
    for code in (part.strip() for part in term_value.split(",")):
        fields = code.split(".")
        for level in (1, 2, 3, 4):
            if len(fields) >= level and all(field.isdigit() for field in fields[:level]):
                found[level].add(".".join(fields[:level]))
    return {level: frozenset(values) for level, values in found.items()}


def cog_levels(term_value: str | None, categories: list[str] | None) -> dict[int, frozenset[str]]:
    """The two COG rungs.

    ⚠ Level 1 is `locus.modal_cog_categories` — the modal CONCATENATED category set over the node's
    annotated members, not the category of the orthogroup at level 2. Most nodes carry one letter
    (4,421 of 4,993 in *E. coli*) but 572 carry two to four, so agreement at this rung is overlap.
    """
    return {
        2: frozenset([term_value]) if term_value else frozenset(),
        1: frozenset(categories or ()),
    }


@dataclass(frozen=True)
class Cell:
    """One (tier, level) cell of the calibration: how often a transfer at this depth was right."""

    tier: str
    level: int
    level_label: str
    pairs: int
    agreeing: int
    chance: float

    @property
    def agreement(self) -> float | None:
        """`None` below `MIN_PAIRS` — a cell with too few pairs reports its n and no rate."""
        return self.agreeing / self.pairs if self.pairs >= MIN_PAIRS else None

    @property
    def interval(self) -> tuple[float, float] | None:
        """A 95 % Wilson interval, which stays inside [0, 1] at the 99 %+ rates the top tiers reach."""
        if self.pairs < MIN_PAIRS:
            return None
        proportion, divisor = self.agreeing / self.pairs, 1 + 1.96**2 / self.pairs
        centre = (proportion + 1.96**2 / (2 * self.pairs)) / divisor
        half = 1.96 * ((proportion * (1 - proportion) / self.pairs
                        + 1.96**2 / (4 * self.pairs**2)) ** 0.5) / divisor
        return (max(0.0, centre - half), min(1.0, centre + half))

    @property
    def lift(self) -> float | None:
        """Observed ÷ chance — **not** an odds ratio.

        ⛔ David, 2026-10-03. Where chance is 0.2 % the two read 394x and 3,432x; the odds ratio is
        unbounded as chance → 0 and cannot be said as *"x more likely than a random annotated node"*,
        which is the sentence the page actually prints.
        """
        rate = self.agreement
        return rate / self.chance if rate is not None and self.chance > 0 else None

    def as_json(self) -> dict:
        """The cell as the page reads it — every rate accompanied by the n it rests on."""
        interval = self.interval
        return {
            "tier": self.tier,
            "level": self.level,
            "level_label": self.level_label,
            "pairs": self.pairs,
            "agreement": self.agreement,
            "interval_low": interval[0] if interval else None,
            "interval_high": interval[1] if interval else None,
            "chance": self.chance,
            "lift": self.lift,
        }


@dataclass(frozen=True)
class Calibration:
    """Every (tier, level) cell for one catalogue, one representation and one vocabulary."""

    annotation_kind: AnnotationKind
    representation: EmbeddingRepresentation
    cells: dict[tuple[str, int], Cell]
    annotated_locus_count: int
    locus_count: int

    def cell(self, tier: str, level: int) -> Cell | None:
        """The cell a transfer at this tier and depth is judged by, or None if nothing measured it."""
        return self.cells.get((tier, level))

    def as_json(self) -> dict:
        """The whole ladder, strongest tier and deepest rung first — the page draws it as a table."""
        return {
            "annotation_kind": self.annotation_kind.value,
            "representation": self.representation.value,
            "annotated_locus_count": self.annotated_locus_count,
            "locus_count": self.locus_count,
            "min_pairs": MIN_PAIRS,
            "cells": [self.cells[key].as_json() for key in sorted(
                self.cells, key=lambda key: ([name for name, _, _ in TIERS].index(key[0]), -key[1]))],
        }


def _claims(session: Session, pangenome_id: int, kind: AnnotationKind):
    """Per locus, its top-rank call folded onto every level of the ladder.

    ⚠ The top-rank entry is already a **modal vote over the node's member genes**
    (`export_payload.top_counts`), with `member_gene_count` as its support — not one gene's call.
    """
    rows = session.execute(
        select(Locus.locus_id, Locus.member_gene_count, Locus.modal_cog_categories,
               LocusAnnotationEntry.term_value, LocusAnnotationEntry.term_name,
               LocusAnnotationEntry.member_gene_count)
        .join(
            LocusAnnotationEntry,
            (LocusAnnotationEntry.locus_id == Locus.locus_id)
            & (LocusAnnotationEntry.annotation_kind == kind)
            & (LocusAnnotationEntry.rank_within_locus == 0),
            isouter=True,
        )
        .where(Locus.pangenome_id == pangenome_id)
    ).all()
    levels: dict[int, dict[int, frozenset[str]]] = {}
    terms: dict[int, tuple[str, str | None, int]] = {}
    locus_count = 0
    for locus_id, _members, categories, term, term_name, support in rows:
        locus_count += 1
        folded = (ec_levels(term) if (kind is AnnotationKind.EC_NUMBER and term)
                  else cog_levels(term, categories) if kind is not AnnotationKind.EC_NUMBER
                  else None)
        if folded and any(folded.values()):
            levels[locus_id] = folded
            if term is not None:
                terms[locus_id] = (term, term_name, support)
    return levels, terms, locus_count


def _chance_by_level(levels: dict[int, dict[int, frozenset[str]]], level: int):
    """An exact chance function for one level: P(a random OTHER annotated node's set intersects S).

    Memoised on the claim set, because a catalogue holds only ~118 distinct COG category values —
    so the inverted-index union is paid once per distinct set rather than once per node.
    """
    holders: dict[str, set[int]] = defaultdict(set)
    for locus_id, folded in levels.items():
        for value in folded[level]:
            holders[value].add(locus_id)
    pool = {locus_id for locus_id, folded in levels.items() if folded[level]}

    @cache
    def chance(claim: frozenset[str]) -> float:
        if not pool:
            return 0.0
        matching: set[int] = set()
        for value in claim:
            matching |= holders.get(value, set())
        #: self-exclusion — the comparator is a random OTHER node, which is what the page claims
        return (len(matching) - 1) / (len(pool) - 1) if len(pool) > 1 else 0.0

    return chance


def tally_cells(
    levels: dict[int, dict[int, frozenset[str]]],
    edges,
    annotation_kind: AnnotationKind,
) -> dict[tuple[str, int], Cell]:
    """Tally every (tier, level) cell over `edges` of `(locus_id, neighbour_locus_id, cosine)`.

    ⭐ Takes plain data rather than a `Session` **so there is one implementation of the measurement,
    not two**: the API reaches it through `compute_calibration` and
    `scripts/measure_cog_function_inference.py` — which is what the write-up's numbers come from —
    reaches it directly. A second copy for the script is exactly how a page and the document that
    describes it come to disagree.
    """
    ladder = sorted(LEVEL_LABEL[annotation_kind], reverse=True)
    chance_functions = {level: _chance_by_level(levels, level) for level in ladder}
    tally: dict[tuple[str, int], list[float]] = defaultdict(lambda: [0, 0, 0.0])
    for locus_id, neighbour_id, cosine in edges:
        mine = levels.get(locus_id)
        theirs = levels.get(neighbour_id)
        if mine is None or theirs is None:
            continue
        tier = tier_for(cosine)
        if tier == NOT_CALLED:
            continue
        for level in ladder:
            if not (mine[level] and theirs[level]):
                continue
            row = tally[(tier, level)]
            row[0] += 1
            row[1] += bool(mine[level] & theirs[level])
            row[2] += chance_functions[level](mine[level])
    return {
        key: Cell(tier=key[0], level=key[1], level_label=LEVEL_LABEL[annotation_kind][key[1]],
                  pairs=int(row[0]), agreeing=int(row[1]), chance=row[2] / row[0])
        for key, row in tally.items() if row[0]
    }


def compute_calibration(
    session: Session,
    *,
    pangenome_id: int,
    annotation_kind: AnnotationKind,
    representation: EmbeddingRepresentation = EmbeddingRepresentation.ESM,
) -> Calibration:
    """Measure how often a transfer at each (tier, level) agrees, over every stored neighbour edge.

    ⭐ **Every rank, not rank 1 only.** 65 % of assignments are supplied by ranks 2-5 — that is the
    population a panel serves — and at fixed cosine the rank effect was <= 0.3 pp in the >= 0.99
    tiers and **flipped sign between species** in the weak ones, so splitting by rank would fit
    noise. It is not cosmetic: on rank-1 pairs alone COG 0.97-0.98 reads 80.9 % / 84.2 % and clears
    an 0.80 bar in both species; pooled it reads 78.1 % / 79.9 % and clears in neither.
    """
    levels, _terms, locus_count = _claims(session, pangenome_id, annotation_kind)
    edges = session.execute(
        select(LocusNearestLocus.locus_id, LocusNearestLocus.neighbour_locus_id,
               LocusNearestLocus.cross_similarity)
        .join(Locus, Locus.locus_id == LocusNearestLocus.locus_id)
        .where(
            Locus.pangenome_id == pangenome_id,
            LocusNearestLocus.representation == representation,
            LocusNearestLocus.cross_similarity.is_not(None),
        )
    ).all()
    cells = tally_cells(levels, edges, annotation_kind)
    return Calibration(
        annotation_kind=annotation_kind,
        representation=representation,
        cells=cells,
        annotated_locus_count=len(levels),
        locus_count=locus_count,
    )


def quotable_level(annotation_kind: AnnotationKind, tier: str,
                   donor_levels: dict[int, frozenset[str]]) -> int | None:
    """The tier's level, **or the deepest the donor actually states, whichever is shallower**.

    ⛔ Not a refinement — without it the page promises a depth the donor cannot supply. Of 84
    *E. coli* nodes with a >= 0.99 EC neighbour only **43** can be given a four-field code; 25 fall
    to L3, 15 to L2 and 1 to L1 because the donor's own EC is dash-padded. The fallbacks are well
    calibrated at the depth they land on (97-99.5 %), so this costs no accuracy — it only stops the
    card claiming a precision that was never there.
    """
    permitted = QUOTED_LEVEL[annotation_kind].get(tier)
    if permitted is None:
        return None
    for level in sorted(LEVEL_LABEL[annotation_kind], reverse=True):
        if level <= permitted and donor_levels[level]:
            return level
    return None
