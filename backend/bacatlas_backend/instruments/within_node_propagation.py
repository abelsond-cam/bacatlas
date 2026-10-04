"""How often a node's own call covers the members that do not carry it — the propagation base rate.

⭐ **This is the OTHER mechanism, and the page had it backwards.** `annotation_transfer.py` measures
**neighbour transfer** — a node with no call of its own borrowing one. This module measures **internal
propagation** — a node whose call is carried by *some* of its genes applying to *all* of them. That is
where the measured lift actually lives, and until now the page described it as a defect: for
*K. pneumoniae* `gumC` (100 genes, COG on 1) it said the members *"cannot be checked against each other
— a single annotated gene agrees with itself by construction"*. True, and the wrong thing to say.
David, 2026-10-04: *"the whole point of our method is effectively the carrying by one is enough."*

⭐ **The claim rests on a measured base rate, NOT on a similarity threshold, and that was measured
rather than assumed.** Within-node unanimity is 99.8-100 % at **every** band of **both**
representations: ESM's within-node similarity is saturated (kp p50 **1.000**) so it cannot discriminate
at all, and Bacformer's spreads (p10 0.743 · p50 0.975) while unanimity stays flat across it. A page
that quoted a floor would imply a gate that is doing no work. (David, 2026-10-04: *base rate only.*)

⛔ **The reference class is `locus.prevalence_band`, and the reason is a population mismatch.** Nodes
where unanimity *can* be checked have a median size of **100 genes**; nodes resting on one annotated
gene have a median of **1**. Most of the latter are RARE-band singletons where nothing propagates at
all, so a single catalogue-wide rate would quote the behaviour of 100-gene core nodes at a 1-gene one,
or the reverse. Banding fixes that with a column that is already published and already on the locus.

⛔⛔ **RARE has ZERO checkable nodes in both catalogues** (and holds most one-gene nodes: 367 kp / 410
ecoli on COG). Its rate is `None`, never 0.0 — *not measurable*, which is a different statement from
*not reliable*. `column_types.measurement()` keeps the same distinction in the schema for the same
reason.

⚠ **Unanimity is tested through `annotation_transfer`'s OWN folders** (`ec_levels` / `cog_levels` /
`single_rung`), so the rung this measures cannot drift from the rung the neighbour ladder quotes. The
test is the one `tally_cells` applies to a pair of nodes, generalised to the N annotated genes of one
node: non-empty intersection at the deepest rung every one of them states.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from sqlalchemy import text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    GENE_VALUE_COLUMN,
    LEVEL_LABEL,
    cog_levels,
    ec_levels,
    single_rung,
)
from bacatlas_backend.models.enumerations import AnnotationKind, PrevalenceBand

#: Below this many checkable nodes a band reports its count and **no rate**. The same stance and the
#: same number as `annotation_transfer.MIN_PAIRS`: a rate off a handful of nodes is not a rate.
MIN_NODES = 30


@dataclass(frozen=True)
class PropagationRate:
    """How often a band's nodes agree internally — the rate a propagated call is quoted at.

    `one_gene` is the population this rate is *applied to* and `checkable` the one it is *measured on*;
    they are reported together because they are different populations and the gap between them is the
    caveat. ⛔ A band with `checkable == 0` can still have a large `one_gene` — that is RARE.
    """

    annotation_kind: AnnotationKind
    #: ⚠ The enum's **NAME** (`CORE`), because that is what Postgres stores and what the group-by
    #: key has to match. `as_json` emits the `.value` instead — see there.
    prevalence_band: str
    checkable: int
    unanimous: int
    one_gene: int

    @property
    def rate(self) -> float | None:
        """`None` below `MIN_NODES` — ⛔ never 0.0, which would read as *measured and bad*."""
        return self.unanimous / self.checkable if self.checkable >= MIN_NODES else None

    @property
    def interval(self) -> tuple[float, float] | None:
        """95 % Wilson — the same idiom as `annotation_transfer.Cell.interval`, for the same reason.

        At the 99.8 % rates these bands reach, a normal interval runs past 1.0 and reads as nonsense.
        """
        if self.checkable < MIN_NODES:
            return None
        proportion = self.unanimous / self.checkable
        divisor = 1 + 1.96**2 / self.checkable
        centre = (proportion + 1.96**2 / (2 * self.checkable)) / divisor
        half = (
            1.96
            * ((proportion * (1 - proportion) / self.checkable + 1.96**2 / (4 * self.checkable**2)) ** 0.5)
            / divisor
        )
        return (max(0.0, centre - half), min(1.0, centre + half))

    def as_json(self) -> dict:
        """What the page reads — every rate accompanied by the node count it rests on.

        ⛔ **`prevalence_band` goes out as the enum's VALUE (`core`), not its name.** Every other
        route serves this field that way (`locus_detail_service:593`, `locus_search_service:147`,
        `audit_residual_service:105`) and the page reads it through one `bandLabel`; a second shape
        for the same field would need a second mapping, and the first reader to miss that gets
        `CORE` rendered raw beside `soft core`.
        """
        interval = self.interval
        return {
            "annotation_kind": self.annotation_kind.value,
            "prevalence_band": PrevalenceBand[self.prevalence_band].value,
            "checkable_node_count": self.checkable,
            "unanimous_node_count": self.unanimous,
            "one_gene_node_count": self.one_gene,
            "rate": self.rate,
            "interval_low": interval[0] if interval else None,
            "interval_high": interval[1] if interval else None,
            "min_nodes": MIN_NODES,
        }


def _fold_gene(annotation_kind: AnnotationKind, value) -> dict[int, frozenset[str]]:
    """One gene's stored value, folded onto the vocabulary's rungs by the SHARED folders.

    ⚠ `cog_levels` takes the locus-level category set for its L1 rung, which a single gene does not
    have — so COG is folded at L2 only here, which is the rung the page displays anyway.
    """
    if annotation_kind is AnnotationKind.EC_NUMBER:
        return ec_levels(",".join(value or ()))
    if annotation_kind is AnnotationKind.COG_ORTHOGROUP:
        return cog_levels(value, None)
    return single_rung(list(value or ()) if isinstance(value, list) else [value])


def unanimous_at_deepest_shared_rung(annotation_kind: AnnotationKind, values: list) -> bool | None:
    """Do these annotated genes agree? `None` where no rung all of them state exists.

    ⭐ **The decision the whole base rate is built from, kept pure so it can be tested without a
    database.** It is `tally_cells`' test — non-empty intersection at a shared rung — generalised from
    two nodes to the N annotated genes of one node.

    ⛔ **The DEEPEST rung every gene states, not the shallowest.** Two genes carrying `1.1.1.1` and
    `1.1.1.2` share `1.1.1` at L3 and differ at L4; taking the shallowest shared rung would call that
    agreement, which is the opposite answer. A shallower rung is a *fallback for a value that does not
    state the deeper one* — exactly what `quotable_level` does on the transfer side — and never a way
    to turn a disagreement into agreement.

    ⛔ **`None` is counted in NEITHER column.** Where no rung is stated by all of them there is nothing
    to compare, and scoring that as agreement is how a measurement quietly flatters itself.
    """
    folded = [_fold_gene(annotation_kind, value) for value in values]
    rung = next(
        (level for level in sorted(LEVEL_LABEL[annotation_kind], reverse=True) if all(one[level] for one in folded)),
        None,
    )
    if rung is None:
        return None
    shared = folded[0][rung]
    for one in folded[1:]:
        shared = shared & one[rung]
    return bool(shared)


def compute_propagation(
    session: Session, *, pangenome_id: int, annotation_kind: AnnotationKind
) -> dict[str, PropagationRate]:
    """Per prevalence band, how often a node's annotated genes agree — one pass over the catalogue.

    ⛔ **Counted from the gene rows**, like `_annotated_counts`, and for the same reason: GO has three
    `*_annotated_member_count` columns (one per namespace) and summing them double-counts every gene
    annotated in more than one.
    """
    column = GENE_VALUE_COLUMN[annotation_kind]
    rows = session.execute(
        text(f"""
            select l.locus_id, l.prevalence_band, f.{column} as value
              from locus l
              join gene_locus_membership m on m.locus_id = l.locus_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where l.pangenome_id = :pangenome_id and f.{column} is not null
        """),
        {"pangenome_id": pangenome_id},
    ).all()

    genes: dict[int, list] = defaultdict(list)
    band_of: dict[int, str] = {}
    for locus_id, band, value in rows:
        band_of[locus_id] = band.name if hasattr(band, "name") else str(band)
        genes[locus_id].append(value)

    tally: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])  # checkable, unanimous, one_gene
    for locus_id, values in genes.items():
        counts = tally[band_of[locus_id]]
        if len(values) == 1:
            counts[2] += 1
            continue
        agreed = unanimous_at_deepest_shared_rung(annotation_kind, values)
        if agreed is None:
            # ⛔ neither unanimous nor discordant — no rung all of them state, so nothing to compare
            continue
        counts[0] += 1
        counts[1] += agreed

    return {
        band: PropagationRate(
            annotation_kind=annotation_kind,
            prevalence_band=band,
            checkable=counts[0],
            unanimous=counts[1],
            one_gene=counts[2],
        )
        for band, counts in tally.items()
    }
