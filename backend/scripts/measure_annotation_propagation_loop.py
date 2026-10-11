r"""How far can annotation propagate? — intrapolation, interpolation, extrapolation, and the loop.

David, 2026-10-11: *"Our actual aim will then be a looping process where we use inter → extra →
intra etc to see how many we get to."* This script answers that, and reports what each mechanism
adds **on top of** the others rather than in isolation.

The three mechanisms, in David's terms:

====================  =========================================================  ======================
term                  mechanism                                                  instrument
====================  =========================================================  ======================
**intrapolation**     within ONE node, from its own annotated member genes       `within_node_propagation`
**interpolation**     node → node, SAME species, by ESM/Bacformer similarity     `annotation_transfer`
**extrapolation**     node → node, ACROSS species                                the cross-species map
====================  =========================================================  ======================

⛔⛔ **TWO MEASUREMENT RULES THIS WHOLE EXERCISE STANDS OR FALLS ON.**

1. **Precision is scored against GROUND TRUTH only, never against a propagated label.** Round 2's
   donors contain round 1's mistakes, so a loop scored against its own output measures nothing —
   and it looks excellent doing it. `State.truth_of` is the only door to a truth value and it
   RAISES on a node whose label was propagated; `test_annotation_propagation_loop.py` pins that.
2. **A masked held-out set, not leave-one-node-out.** The established idiom — hide one node's label
   and ask whether its donor agrees — cannot evaluate a *loop*: every other node still carries
   truth, so the loop never meets its own errors. So a fixed share of annotated nodes is masked up
   front, the whole loop runs on the remainder, and only the masked nodes are scored.

   ⚠ **The bias runs BOTH ways and does not cancel. Every table says so.**
   Masking makes the donor graph *sparser* than reality, so measured coverage is **conservative**;
   but precision is measured on nodes that really do carry a label, standing in for nodes that do
   not, and unannotated nodes sit in sparser neighbourhoods by construction — which is
   **optimistic**. Neither is quantifiable from inside this data.

⛔ **Thresholds are selected on the OTHER catalogue**, never on the set they are scored on.
Interpolation chooses on ecoli and scores on kp, and back. Extrapolation is one-directional
(ecoli → kp), so there is no other catalogue: its masked set is split in half, the bar is chosen on
one half and reported on the other, and the output says which number is which.

⭐ **The claim sets, the folds and the chance function are IMPORTED** from
`instruments/annotation_transfer.py` and `instruments/gene_ontology_fold.py` — the same code the API
serves the Function tab from. Nothing here re-implements a fold, a claim or a rung.

⛔⛔ **Nothing here goes near the web page**, and nothing measured here is RECORDED until David has
seen it (`~/.claude/CLAUDE.md`: gather the data, then STOP).

Run from `backend/` with both `*-nuna4` catalogues loaded::

    export BACATLAS_DATABASE_URL="postgresql+psycopg://$USER@localhost:5432/bacatlas_dev"
    python scripts/measure_annotation_propagation_loop.py --stage all \
        --cross-esm ../../nuna/data/proc/analysis/cross_species/ecoli_to_kp_..._cross_species_neighbours.tsv \
        --cross-bacformer ../../nuna/data/proc/analysis/cross_species/ecoli_to_kp_..._BACFORMER_....tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy import create_engine, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import _chance_by_rung, _claims, go_rung
from bacatlas_backend.models.enumerations import (
    GENE_ONTOLOGY_NAMESPACE_NAMES,
    AnnotationKind,
    EmbeddingRepresentation,
)
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_similarity import LocusNearestLocus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.services.reference_name_service import gene_ontology_fold

COG = AnnotationKind.COG_ORTHOGROUP
GO = AnnotationKind.GENE_ONTOLOGY_SLIM

#: The four claims the loop runs on, each a `(key, label, kind, rung)`. ⛔ GO's three facets are
#: SEPARATE claims and never pooled: at 0.94-0.96 the pooled rate averages a molecular-function
#: claim worth 60 % with a cellular-component one worth 84 %, which is how ~1,400 suggestions per
#: catalogue came to rest on the weak one. EC and KEGG are left out — 31.9 % and 8.2 % gene
#: coverage at source, so neither has the reach to carry a propagation argument.
CLAIMS: tuple[tuple[str, str, AnnotationKind, tuple[str | None, int]], ...] = (
    ("cog", "COG orthogroup", COG, (None, 2)),
    ("mf", "molecular function", GO, (GENE_ONTOLOGY_NAMESPACE_NAMES[0], 1)),
    ("bp", "biological process", GO, (GENE_ONTOLOGY_NAMESPACE_NAMES[1], 1)),
    ("cc", "cellular component", GO, (GENE_ONTOLOGY_NAMESPACE_NAMES[2], 1)),
)
CLAIM_LABEL = {key: label for key, label, _kind, _rung in CLAIMS}

#: ⭐ The interpolation policies the loop STARTS from — measured 2026-10-06/11 on the closed,
#: per-namespace GO claim, chosen on one catalogue and verified on the other, every one clearing
#: 90 % in BOTH species. `(esm_bar, bacformer_bar)` under the **concur** rule.
#: ⚠ Concurrence between the two geometries is doing the work a similarity threshold was doing:
#: once both maps must agree, Bacformer's own bar barely matters, which is why two of these are
#: `None`. ⛔ The `either` (union) rule generalises badly across species (83.3 % held out on BP) and
#: is not offered here.
INTERPOLATION_POLICY: dict[str, tuple[float | None, float | None]] = {
    "cog": (0.96, None),
    "mf": (0.97, None),
    "bp": (None, None),
    "cc": (None, None),
}

#: The share of annotated nodes held out, and the seed that fixes which. ⚠ One mask per
#: (catalogue, claim), reused by every mechanism, or the mechanisms are not comparable.
MASK_FRACTION = 0.20
MASK_SEED = 0
#: A policy that calls 12 nodes at 100 % is not a policy.
MIN_CALLED = 100
#: The bar every mechanism is held to.
TARGET = 0.90
#: Loop rounds before we stop regardless.
MAX_ROUNDS = 5
#: ⛔ A mechanism is retired only on a rate measured over at least this many held-out nodes. The
#: same stance and the same number as `annotation_transfer.MIN_PAIRS`: a rate off a handful is not a
#: rate, and retiring a mechanism on noise would understate the loop.
MIN_SCORED_TO_RETIRE = 30

CATALOGUES = ("ecoli-nuna4", "kp-nuna4")
REPS = ("ESM", "BACFORMER")
#: Cross-species bars swept in S5. ⚠ Cross-species cosines are a THIRD distribution: the ESM map's
#: rank-1 pairs sit at ~0.976 and Bacformer's at ~0.35, so neither within-species ladder transfers.
CROSS_BARS = {
    "ESM": (None, 0.90, 0.95, 0.96, 0.97, 0.975, 0.98),
    "BACFORMER": (None, 0.20, 0.30, 0.35, 0.40, 0.50, 0.60),
}


# ======================================================================== the problem, off the database
@dataclass
class Problem:
    """One (catalogue, claim) instance — the truth, the node sizes, and the two neighbour walks.

    ⚠ `truth` is the node's own label, taken from the rank-0 `locus_annotation_entry` rows, which
    are already a modal vote over the node's member genes. It is the same claim the Function tab's
    ladder is calibrated on, so the policies above mean here what they mean there.
    """

    catalogue: str
    claim_key: str
    truth: dict[int, frozenset[str]]
    #: Every node of the catalogue, annotated or not → its member gene count.
    gene_count: dict[int, int]
    #: Node → member genes carrying THIS claim. ⛔ Zero for a node with no truth label.
    genes_carrying: dict[int, int]
    #: rep → node → `[(rank, neighbour, cosine)]`, ascending rank.
    walk: dict[str, dict[int, list[tuple[int, int, float]]]]
    node_label: dict[int, str]
    #: UniRef50 family → the nodes holding a gene of it (S7 only).
    family_nodes: dict[str, set[int]]
    node_families: dict[int, set[str]]

    @property
    def node_total(self) -> int:
        """Every node of the catalogue — the denominator of node coverage."""
        return len(self.gene_count)

    @property
    def gene_total(self) -> int:
        """Every member gene of the catalogue — the denominator of gene coverage."""
        return sum(self.gene_count.values())


@dataclass
class State:
    """Which nodes carry a label right now, where each came from, and what the truth was.

    ⛔⛔ **`truth_of` is the ONLY door to a truth value and it RAISES on a propagated node.** That is
    measurement rule 1 expressed as code rather than as a comment: a loop scored against its own
    output measures nothing, and the mistake is invisible in the output.
    """

    label: dict[int, frozenset[str]]
    #: node → `"truth"`, or `"<mechanism>:<round>"` for anything the loop added.
    source: dict[int, str]
    #: The held-out nodes, with their hidden truth. The only nodes precision may be scored on.
    masked: dict[int, frozenset[str]]
    #: ⭐⭐ **Every certainty signal the call had available, per node.** David, 2026-10-11:
    #: *"a similarity of 0.8 with esm of 0.99 must be better than bacformer 0.5 and esm 0.99."*
    #: One cosine cannot express that, so the record carries four things:
    #:
    #: * `esm_cosine` — the ESM cosine to the donor actually used;
    #: * `bacformer_cosine` — the cosine to the donor BACFORMER independently chose, which under the
    #:   concurrence rule is a different node and measures how strongly the second geometry backed a
    #:   concurring claim;
    #: * `pair_bacformer_cosine` — the Bacformer cosine of **the same donor pair**, which is what
    #:   David described. ⛔ It is `None` wherever that pair is absent from the Bacformer shortlist,
    #:   and the share where it exists is REPORTED rather than quietly dropped: within a species the
    #:   stored walk is only five deep, so the same pair is often simply not there.
    donor: dict[int, dict] = field(default_factory=dict)

    def truth_of(self, node: int) -> frozenset[str]:
        """The held-out truth for a masked node. ⛔ RAISES for anything else."""
        if node not in self.masked:
            raise AssertionError(
                f"node {node} is not in the masked held-out set, so it has no scoreable truth — "
                "scoring a propagated label against itself is measurement rule 1"
            )
        return self.masked[node]

    def copy(self) -> State:
        """A detached copy, so a mechanism can be measured alone without disturbing the loop."""
        return State(dict(self.label), dict(self.source), self.masked, dict(self.donor))


def _pangenome_id(session: Session, catalogue: str) -> int:
    """The catalogue's surrogate key."""
    return session.execute(
        select(Pangenome.pangenome_id).where(Pangenome.catalogue_key == catalogue)
    ).scalar_one()


def _gene_carriage(session: Session, pangenome_id: int, fold) -> tuple[dict[str, dict[int, int]], dict[str, int]]:
    """Per claim, per node, how many member genes carry it — plus the stored-column gap.

    ⛔ **Counted from the GENE rows, and the GO facets are folded**, because
    `locus.go_annotated_member_count_*` counts genes carrying a *term* in that namespace, which is
    not the same as genes with a non-empty *claim* there: a gene stating only the namespace root
    (452 ecoli loci state nothing but `cellular_component`) is counted by the column and carries no
    comparable claim. The gap is returned and REPORTED rather than assumed away.
    """
    carrying: dict[str, dict[int, int]] = {key: defaultdict(int) for key, *_ in CLAIMS}
    for locus_id, count in session.execute(
        text(
            """
            select l.locus_id, l.cog_annotated_member_count
              from locus l where l.pangenome_id = :pid and l.cog_annotated_member_count > 0
            """
        ),
        {"pid": pangenome_id},
    ).all():
        carrying["cog"][locus_id] = count

    facet_key = {GENE_ONTOLOGY_NAMESPACE_NAMES[index]: key for index, key in enumerate(("mf", "bp", "cc"))}
    for locus_id, terms in session.execute(
        text(
            """
            select m.locus_id, f.gene_ontology_terms
              from gene_locus_membership m
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where m.pangenome_id = :pid and f.gene_ontology_terms is not null
            """
        ),
        {"pid": pangenome_id},
    ).all():
        for (facet, _depth), classes in go_rung([t for t in (terms or ()) if t], fold).items():
            if classes:
                carrying[facet_key[facet]][locus_id] += 1

    # ⛔ The three sums are ALIASED. Postgres names every un-aliased `sum(...)` column `sum`, so a
    # mapping built from the row would hold one key and silently answer for whichever came last.
    stored = session.execute(
        text(
            """
            select sum(go_annotated_member_count_molecular_function) as mf,
                   sum(go_annotated_member_count_biological_process) as bp,
                   sum(go_annotated_member_count_cellular_component) as cc
              from locus where pangenome_id = :pid
            """
        ),
        {"pid": pangenome_id},
    ).one()
    gap = {
        key: int(getattr(stored, key) or 0) - sum(carrying[key].values()) for key in ("mf", "bp", "cc")
    }
    return {key: dict(value) for key, value in carrying.items()}, gap


def load_problems(session: Session, catalogue: str) -> tuple[dict[str, Problem], dict[str, int]]:
    """Every claim's problem instance for one catalogue, in one pass over the heavy tables."""
    pangenome_id = _pangenome_id(session, catalogue)
    fold = gene_ontology_fold(session)

    gene_count: dict[int, int] = {}
    node_label: dict[int, str] = {}
    for locus_id, label, members in session.execute(
        select(Locus.locus_id, Locus.node_label, Locus.member_gene_count).where(
            Locus.pangenome_id == pangenome_id
        )
    ).all():
        gene_count[locus_id] = members
        node_label[locus_id] = label

    walk: dict[str, dict[int, list[tuple[int, int, float]]]] = {}
    for rep in (EmbeddingRepresentation.ESM, EmbeddingRepresentation.BACFORMER):
        one: dict[int, list[tuple[int, int, float]]] = defaultdict(list)
        for node, rank, neighbour, cosine in session.execute(
            select(
                LocusNearestLocus.locus_id,
                LocusNearestLocus.rank,
                LocusNearestLocus.neighbour_locus_id,
                LocusNearestLocus.cross_similarity,
            )
            .join(Locus, Locus.locus_id == LocusNearestLocus.locus_id)
            .where(
                Locus.pangenome_id == pangenome_id,
                LocusNearestLocus.representation == rep,
                LocusNearestLocus.cross_similarity.is_not(None),
            )
            .order_by(LocusNearestLocus.locus_id, LocusNearestLocus.rank)
        ).all():
            one[node].append((rank, neighbour, cosine))
        walk[rep.name] = dict(one)

    family_nodes: dict[str, set[int]] = defaultdict(set)
    node_families: dict[int, set[str]] = defaultdict(set)
    for locus_id, accession in session.execute(
        text(
            """
            select distinct m.locus_id, f.uniref50_accession
              from gene_locus_membership m
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where m.pangenome_id = :pid and f.uniref50_accession is not null
            """
        ),
        {"pid": pangenome_id},
    ).all():
        family_nodes[accession].add(locus_id)
        node_families[locus_id].add(accession)

    carrying, gap = _gene_carriage(session, pangenome_id, fold)

    problems: dict[str, Problem] = {}
    for key, _label, kind, rung in CLAIMS:
        folded, _terms, _total = _claims(session, pangenome_id, kind, fold)
        truth = {node: value[rung] for node, value in folded.items() if value.get(rung)}
        problems[key] = Problem(
            catalogue=catalogue,
            claim_key=key,
            truth=truth,
            gene_count=gene_count,
            genes_carrying=carrying[key],
            walk=walk,
            node_label=node_label,
            family_nodes=dict(family_nodes),
            node_families=dict(node_families),
        )
    return problems, gap


# ============================================================================= the mask, and the census
def masked_state(problem: Problem, *, fraction: float = MASK_FRACTION, seed: int = MASK_SEED) -> State:
    """Hold out a random share of the annotated nodes — the only nodes precision is scored on."""
    annotated = sorted(problem.truth)
    chosen = set(random.Random(seed).sample(annotated, round(len(annotated) * fraction)))
    return State(
        label={node: value for node, value in problem.truth.items() if node not in chosen},
        source={node: "truth" for node in problem.truth if node not in chosen},
        masked={node: problem.truth[node] for node in chosen},
    )


def census(problem: Problem, state: State) -> dict[str, float | int]:
    """⭐ David's FOUR numbers, every time — nodes (n, %) and genes (n, %).

    Two gene numbers, because intrapolation is the difference between them and nothing else is:
    `genes_carrying` counts the member genes that hold the claim themselves **anywhere in the
    catalogue**, `genes_covered` counts every gene of every labelled node plus the carrying genes of
    the nodes that have no label — which is what the node's label means once intra has run.
    ⛔ Node coverage is IDENTICAL in both: intra cannot annotate a node that had nothing.
    """
    nodes = len(state.label)
    # ⛔⛔ **A gene is annotated if IT carries the claim, OR its node holds a label.** Summing
    # `genes_carrying` over LABELLED nodes only — which this did until 2026-10-11 — silently drops
    # every gene whose node has no node-level claim in this facet. The node claim is the rank-0
    # modal vote, so a node can hold genes stating a molecular function and still have no rank-0
    # molecular-function row. Measured: that excluded 25,751 ecoli molecular-function genes from the
    # baseline (5.26 points), and 1,157 of them sat inside nodes the pass later labelled, so the
    # pass was credited with genes the annotation source had already annotated.
    carrying = sum(problem.genes_carrying.values())
    covered = sum(
        problem.gene_count[node] if node in state.label else problem.genes_carrying.get(node, 0)
        for node in problem.gene_count
    )
    return {
        "nodes": nodes,
        "node_share": nodes / problem.node_total,
        "genes_carrying": carrying,
        "genes_carrying_share": carrying / problem.gene_total,
        "genes_covered": covered,
        "genes_covered_share": covered / problem.gene_total,
    }


def precision(state: State, assigned: dict[int, frozenset[str]]) -> tuple[int, int, float | None]:
    """`(called, scoreable, precision)` — scoreable counts ONLY masked nodes.

    ⛔ A node that was never annotated has no truth and is counted in `called` but never in
    `scoreable`. That asymmetry is the point: coverage is what we gain, precision is what we can
    check, and they are measured on different populations.
    """
    scoreable = right = 0
    for node, claim in assigned.items():
        if node not in state.masked:
            continue
        scoreable += 1
        right += bool(claim & state.truth_of(node))
    return len(assigned), scoreable, (right / scoreable if scoreable else None)


# ================================================================================== the three mechanisms
def nearest_labelled(problem: Problem, state: State, node: int, rep: str) -> tuple[frozenset[str], float, int] | None:
    """The nearest neighbour of `node` that carries a label right now, in one representation.

    ⛔ **The walk is only FIVE deep** (`locus_nearest_locus` stores `rank <= 5`), and that is the
    binding constraint on everything below, not the thresholds. A node whose five nearest
    neighbours are all unlabelled can never receive a label however many rounds run, and the loop's
    gain after round 1 is exactly the nodes whose top five *became* labelled. ⚠ The donor is the
    nearest **labelled** node, not the nearest node — so its cosine may come from rank 4, which is
    what the live service does and why the tier is read off the donor rather than off rank 1.
    """
    for _rank, neighbour, cosine in problem.walk[rep].get(node, ()):
        if neighbour != node and neighbour in state.label:
            return state.label[neighbour], cosine, neighbour
    return None


def _pair_cosine(
    walk: dict[int, list[tuple[int, int, float]]], node: int, neighbour: int
) -> float | None:
    """The cosine of this exact pair in `walk`, or `None` if the pair is not in the stored shortlist.

    ⛔ `None` means *not stored*, which is a different statement from *not similar* — the within-
    species walk keeps five neighbours, so most pairs are simply absent. Treating absence as a low
    cosine would manufacture a correlation between certainty and shortlist depth.
    """
    for _rank, candidate, cosine in walk.get(node, ()):
        if candidate == neighbour:
            return cosine
    return None


def interpolate(
    problem: Problem, state: State, policy: tuple[float | None, float | None], *, rule: str = "concur"
) -> dict[int, frozenset[str]]:
    """Unlabelled nodes → the claim their two nearest labelled neighbours CONCUR on.

    ⚠ **`concur` is a consensus test, not two independent predictions**, and the difference matters
    for how its precision should be read. A COG claim is a single accession, so *"the two claims
    intersect"* is *"the two claims are equal"*; `ESM|concur` and `BAC|concur` are therefore
    identical columns by construction. What concurrence buys is the agreement of two *geometries*
    about the same answer, and 19-28 % of nodes get the very same donor from both maps, where it
    buys nothing at all.
    """
    esm_bar, bacformer_bar = policy
    assigned: dict[int, frozenset[str]] = {}
    for node in problem.gene_count:
        if node in state.label:
            continue
        esm = nearest_labelled(problem, state, node, "ESM")
        bacformer = nearest_labelled(problem, state, node, "BACFORMER")
        esm_ok = esm is not None and (esm_bar is None or esm[1] >= esm_bar)
        bacformer_ok = bacformer is not None and (bacformer_bar is None or bacformer[1] >= bacformer_bar)
        if rule == "esm":
            chosen = esm if esm_ok else None
        elif rule == "bacformer":
            chosen = bacformer if bacformer_ok else None
        else:
            chosen = esm if (esm_ok and bacformer_ok and (esm[0] & bacformer[0])) else None
        if chosen is not None:
            assigned[node] = chosen[0]
            # ⚠ The same-pair Bacformer cosine exists only if the ESM donor also sits in this node's
            # Bacformer top-5. `_pair_cosine` returns None otherwise and the stratifier reports how
            # often that happened.
            state.donor[node] = {
                "rule": rule,
                "donor": chosen[2],
                "esm_cosine": esm[1] if esm else None,
                "bacformer_cosine": bacformer[1] if bacformer else None,
                "pair_bacformer_cosine": _pair_cosine(problem.walk["BACFORMER"], node, chosen[2]),
            }
    return assigned


def intrapolate(problem: Problem, state: State) -> tuple[int, int]:
    """`(nodes gained, genes gained)` — and the first is ALWAYS zero, by construction.

    ⭐ Intrapolation is a **gene**-coverage mechanism, not a node-coverage one: it spreads a label a
    node already has to the members that do not carry it, measured at 99.8-100 % in every
    prevalence band of both representations. It cannot annotate a node that had nothing, so the
    report prints a zero rather than a flat line and says why.
    """
    genes = sum(
        problem.gene_count[node] - problem.genes_carrying.get(node, 0) for node in state.label
    )
    return 0, genes


def read_cross_map(
    path: str, *, donor: str, recipient: str, depth: int | None = None
) -> dict[str, list[tuple[str, float, int]]]:
    """The cross-species map as `recipient node label → [(donor node label, cosine, rank)]`.

    ⛔⛔ **THE DIRECTION IS NOT IN THE HEADER.** A map in either direction carries both `*_node`
    columns, so a header check reads an `ecoli → kp` map as a `kp → ecoli` one without complaint and
    every rank then means the opposite of what it says. The direction is taken from the provenance
    JSON beside the TSV, which records it explicitly; the header is a second, independent check that
    catches a truncated or older map.

    ⛔ Read with `csv`, so every node label arrives as `str` and stays one. `locus.node_label` is
    TEXT and labels *look* numeric (`"1098"`); a reader that inferred dtypes would join nothing at
    all while looking perfectly reasonable.
    """
    sidecar = Path(path).with_suffix(".json")
    if not sidecar.exists():
        raise SystemExit(f"FATAL: no provenance beside the map — expected {sidecar}")
    provenance = json.loads(sidecar.read_text())
    stated = (provenance.get("donor"), provenance.get("recipient"))
    if stated != (donor, recipient):
        raise SystemExit(
            f"FATAL: {sidecar.name} says donor={stated[0]} recipient={stated[1]}; "
            f"expected {donor} → {recipient}"
        )
    #: ⛔⛔ **The shortlist depth is the confound in every extra-vs-inter comparison.** The
    #: within-species walk stores FIVE neighbours; this map carries a mean of 23.3 (ESM) and 35.6
    #: (Bacformer) donors per recipient node, capped at 30 by the run that built it. So extrapolation
    #: looking further than interpolation is a property of how the two maps were BUILT, and `depth`
    #: exists so the two can be compared at equal depth rather than argued about.
    out: dict[str, list[tuple[str, float, int]]] = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for column in (f"{recipient}_node", f"{donor}_node", "cross", "rank_cross"):
            if column not in (reader.fieldnames or ()):
                raise SystemExit(f"FATAL: {path} has no column {column!r} — columns are {reader.fieldnames}")
        for row in reader:
            out[row[f"{recipient}_node"]].append(
                (row[f"{donor}_node"], float(row["cross"]), int(row["rank_cross"]))
            )
    ranked = {node: sorted(donors, key=lambda item: item[2]) for node, donors in out.items()}
    if depth is not None:
        ranked = {node: donors[:depth] for node, donors in ranked.items()}
    return ranked


def extrapolate(
    recipient: Problem,
    recipient_state: State,
    donor: Problem,
    donor_state: State,
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    bars: dict[str, float | None],
    *,
    rule: str = "concur",
) -> dict[int, frozenset[str]]:
    """Unlabelled recipient nodes → a claim carried across the species boundary.

    ⚠ The donor's labels are read from the **donor's own working state**, so a donor node whose
    label was itself interpolated can pass it on — which is the loop, and which is also how one
    wrong interpolation crosses a species boundary. Provenance records it.
    """
    donor_id_of = {label: node for node, label in donor.node_label.items()}
    assigned: dict[int, frozenset[str]] = {}
    for node in recipient.gene_count:
        if node in recipient_state.label:
            continue
        found: dict[str, tuple[frozenset[str], float, int]] = {}
        for rep in REPS:
            for donor_label, cosine, _rank in maps[rep].get(recipient.node_label[node], ()):
                donor_node = donor_id_of.get(donor_label)
                if donor_node is not None and donor_node in donor_state.label:
                    found[rep] = (donor_state.label[donor_node], cosine, donor_node)
                    break
        ok = {
            rep: rep in found and (bars[rep] is None or found[rep][1] >= bars[rep]) for rep in REPS
        }
        if rule == "esm":
            chosen = found.get("ESM") if ok["ESM"] else None
        elif rule == "bacformer":
            chosen = found.get("BACFORMER") if ok["BACFORMER"] else None
        else:
            chosen = (
                found["ESM"]
                if all(ok.values()) and (found["ESM"][0] & found["BACFORMER"][0])
                else None
            )
        if chosen is not None:
            assigned[node] = chosen[0]
            # ⭐ Across species the same-pair lookup almost always succeeds: the maps carry a mean of
            # 23.3 (ESM) and 35.6 (Bacformer) donors per node against the walk's five.
            donor_label = donor.node_label[chosen[2]]
            pair = next(
                (
                    cosine
                    for label, cosine, _rank in maps["BACFORMER"].get(
                        recipient.node_label[node], ()
                    )
                    if label == donor_label
                ),
                None,
            )
            recipient_state.donor[node] = {
                "rule": f"extra:{rule}",
                "donor": chosen[2],
                "esm_cosine": found["ESM"][1] if "ESM" in found else None,
                "bacformer_cosine": found["BACFORMER"][1] if "BACFORMER" in found else None,
                "pair_bacformer_cosine": pair,
            }
    return assigned


def uniref_propagate(problem: Problem, state: State) -> dict[int, frozenset[str]]:
    """Unlabelled nodes → the claim of a labelled node sharing a UniRef50 family.

    ⭐ David's benchmark arm, run in two positions — before the embedding steps and after them. The
    marginal gain in the second position is how much UniRef50 knows that ESM and Bacformer did not
    already find.

    ⚠ **Its error mode is its own**: a family spanning two nodes may span two functions, which is
    precisely the paralogue split the clustering exists to make. Where the labelled nodes of a
    family disagree, the intersection is taken and an empty intersection assigns nothing — the same
    stance as `unanimous_at_deepest_shared_rung`, where no comparable rung votes neither way.
    """
    assigned: dict[int, frozenset[str]] = {}
    for node in problem.gene_count:
        if node in state.label:
            continue
        shared: frozenset[str] | None = None
        witness: int | None = None
        for family in problem.node_families.get(node, ()):
            for other in problem.family_nodes.get(family, ()):
                if other == node or other not in state.label:
                    continue
                shared = state.label[other] if shared is None else shared & state.label[other]
                witness = other
        if shared:
            assigned[node] = shared
            state.donor[node] = {
                "rule": "uniref50",
                "donor": witness,
                "esm_cosine": None,
                "bacformer_cosine": None,
                "pair_bacformer_cosine": None,
            }
    return assigned


def precision_by_gene(problem: Problem, state: State, assigned: dict[int, frozenset[str]]) -> float | None:
    """Precision WEIGHTED BY NODE SIZE — where the amplification shows.

    ⭐ Open question 1, made visible rather than argued about. Once a node has a label, intra
    spreads it to every member, so one wrong interpolation on a 100-gene core node costs 100 genes
    and one on a singleton costs 1. Node-level precision cannot see that difference; this can.
    """
    weight = right = 0
    for node, claim in assigned.items():
        if node not in state.masked:
            continue
        weight += problem.gene_count[node]
        right += problem.gene_count[node] * bool(claim & state.truth_of(node))
    return right / weight if weight else None


def apply(state: State, assigned: dict[int, frozenset[str]], mechanism: str, round_number: int) -> None:
    """Commit a round's labels, each stamped with the mechanism and round that produced it."""
    for node, claim in assigned.items():
        state.label[node] = claim
        state.source[node] = f"{mechanism}:{round_number}"


def _pct(value: float | None) -> str:
    """A rate, or `—` where there was nothing to measure. ⛔ Never 0.0 for *not measured*."""
    return "      —" if value is None else f"{value:>7.1%}"


# ============================================================================================ reporting
def report_baseline(problems: dict[str, dict[str, Problem]], gaps: dict[str, dict[str, int]]) -> None:
    """S2 — where every claim starts, in nodes and in genes, before anything propagates."""
    print(f"\n{'=' * 118}\nS2 — THE BASELINE CENSUS\n{'=' * 118}")
    print("  ⚠ 'genes carrying' counts member genes holding the claim themselves; 'genes covered' counts")
    print("    every gene of every annotated node — the two differ by exactly what intrapolation does.\n")
    header = f"  {'catalogue':<12}{'claim':<22}{'nodes':>8}  {'of all':>7}"
    print(header + f"{'genes carry':>13}{'of all':>8}{'genes cover':>13}{'of all':>8}")
    for catalogue in CATALOGUES:
        for key, label, _kind, _rung in CLAIMS:
            problem = problems[catalogue][key]
            whole = State(dict(problem.truth), {}, {})
            numbers = census(problem, whole)
            print(
                f"  {catalogue:<12}{label:<22}{numbers['nodes']:>8,}  {numbers['node_share']:>7.1%}"
                f"{numbers['genes_carrying']:>13,}{numbers['genes_carrying_share']:>8.1%}"
                f"{numbers['genes_covered']:>13,}{numbers['genes_covered_share']:>8.1%}"
            )
        problem = problems[catalogue]["cog"]
        print(f"    {catalogue}: {problem.node_total:,} nodes · {problem.gene_total:,} member genes")
    print("\n  ⚠ Gene-carriage gap against `locus.go_annotated_member_count_*` — genes the stored column")
    print("    counts that state nothing comparable (a namespace root only):")
    for catalogue in CATALOGUES:
        print(f"    {catalogue:<12}" + " · ".join(f"{key} {gaps[catalogue][key]:+,}" for key in ("mf", "bp", "cc")))


def report_chance(problems: dict[str, dict[str, Problem]], states: dict[str, dict[str, State]]) -> None:
    """The null every agreement rate is read against — from `annotation_transfer`'s own function.

    ⛔ **A statistical decision, not a convenience** (`annotation_transfer._chance_by_rung`). Within a
    species the comparator is *"a random OTHER node of this catalogue"*, so the focal node is removed
    from both numerator and denominator. Across species the pool is the DONOR catalogue and the focal
    node was never in it, so subtracting one would deflate chance and inflate every lift read against
    it. Both are printed, because a 91 % agreement means something different at a 40 % null than at a
    2 % one.

    ⚠ David, 2026-10-06: *"the absolute agreement is what we mostly want to know"* — so this is
    context for the tables, not the headline.
    """
    print(f"\n{'=' * 118}\nS2c — THE NULL: what a random annotated donor would agree with anyway\n{'=' * 118}")
    print("  Averaged over each catalogue's masked held-out truths. 'within' self-excludes; 'cross' does not.\n")
    print(f"  {'claim':<22}{'ecoli within':>14}{'kp within':>12}{'ecoli→kp cross':>16}")
    for key, label, _kind, rung in CLAIMS:
        shown = []
        for catalogue in CATALOGUES:
            state = states[catalogue][key]
            pool = {node: {rung: claim} for node, claim in state.label.items()}
            chance = _chance_by_rung(pool, rung)
            shown.append(sum(chance(truth) for truth in state.masked.values()) / max(1, len(state.masked)))
        donor = {
            node: {rung: claim} for node, claim in states["ecoli-nuna4"][key].label.items()
        }
        cross = _chance_by_rung(donor, rung, self_excluded=False)
        kp_masked = states["kp-nuna4"][key].masked
        across = sum(cross(truth) for truth in kp_masked.values()) / max(1, len(kp_masked))
        print(f"  {label:<22}{shown[0]:>13.1%} {shown[1]:>11.1%} {across:>15.1%}")


def report_intra(problems: dict[str, dict[str, Problem]]) -> None:
    """S3 — intrapolation alone. Node coverage does not move, and that is by construction."""
    print(f"\n{'=' * 118}\nS3 — INTRAPOLATION ALONE (within one node, from its own annotated members)\n{'=' * 118}")
    print("  ⛔ Nodes gained is ZERO for every row: intra cannot annotate a node that had nothing. It is a")
    print("     GENE-coverage mechanism, measured at 99.8-100 % in every prevalence band of both maps.\n")
    print(f"  {'catalogue':<12}{'claim':<22}{'nodes +':>9}{'genes +':>12}{'gene cover':>12}{'→':>4}{'after':>9}")
    for catalogue in CATALOGUES:
        for key, label, _kind, _rung in CLAIMS:
            problem = problems[catalogue][key]
            whole = State(dict(problem.truth), {}, {})
            before = census(problem, whole)
            nodes, genes = intrapolate(problem, whole)
            print(
                f"  {catalogue:<12}{label:<22}{nodes:>9,}{genes:>12,}"
                f"{before['genes_carrying_share']:>12.1%}{'→':>4}{before['genes_covered_share']:>9.1%}"
            )


ESM_BARS = (None, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
BAC_BARS = (None, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)


def report_inter(problems: dict[str, dict[str, Problem]], states: dict[str, dict[str, State]]) -> None:
    """S4 — interpolation alone, on the masked set, with the policy RE-CHOSEN on the other catalogue."""
    print(f"\n{'=' * 118}\nS4 — INTERPOLATION ALONE (node → node, same species)\n{'=' * 118}")
    print("  Precision is measured on the masked held-out nodes ONLY. Coverage counts every node called.")
    print("  ⚠ conservative in coverage (the donor graph is sparser than reality) and optimistic in")
    print("    precision (annotated nodes stand in for unannotated ones, which sit in sparser places).\n")
    print(
        f"  {'catalogue':<12}{'claim':<22}{'policy':<26}{'called':>8}"
        f"{'scored':>8}{'prec':>8}{'by gene':>9}{'node cover':>12}"
    )
    for catalogue in CATALOGUES:
        for key, label, _kind, _rung in CLAIMS:
            problem, state = problems[catalogue][key], states[catalogue][key].copy()
            policy = INTERPOLATION_POLICY[key]
            assigned = interpolate(problem, state, policy)
            called, scored, rate = precision(state, assigned)
            by_gene = precision_by_gene(problem, state, assigned)
            after = (len(state.label) + called) / problem.node_total
            esm_bar, bacformer_bar = policy
            shown = f"concur esm≥{esm_bar or 'any'} bac≥{bacformer_bar or 'any'}"
            print(
                f"  {catalogue:<12}{label:<22}{shown:<26}{called:>8,}{scored:>8,}"
                f"{_pct(rate)}{_pct(by_gene)}{after:>12.1%}"
            )

    print("\n  ⭐ The policy RE-CHOSEN from scratch, on one catalogue's masked set, scored on the other's.")
    print("     ⛔ This is the number that counts; a grid of ~60 cells will always hold one that clears 90 %.")
    print(f"  {'claim':<22}{'chosen on':<12}{'policy':<30}{'own':>8}{'HELD OUT':>10}{'called':>9}")
    for key, label, _kind, _rung in CLAIMS:
        grids: dict[str, dict[tuple, tuple[int, float | None]]] = {}
        for catalogue in CATALOGUES:
            problem, base = problems[catalogue][key], states[catalogue][key]
            grid: dict[tuple, tuple[int, float | None]] = {}
            for rule in ("esm", "bacformer", "concur"):
                for esm_bar in ESM_BARS if rule != "bacformer" else (None,):
                    for bacformer_bar in BAC_BARS if rule != "esm" else (None,):
                        state = base.copy()
                        assigned = interpolate(problem, state, (esm_bar, bacformer_bar), rule=rule)
                        _called, scored, rate = precision(state, assigned)
                        grid[(rule, esm_bar, bacformer_bar)] = (scored, rate)
            grids[catalogue] = grid
        for chosen_on, held_out in (CATALOGUES, tuple(reversed(CATALOGUES))):
            viable = [
                (policy, scored, rate)
                for policy, (scored, rate) in grids[chosen_on].items()
                if rate is not None and rate >= TARGET and scored >= MIN_CALLED
            ]
            if not viable:
                best = max(grids[chosen_on].items(), key=lambda item: item[1][1] or 0)
                print(
                    f"  {label:<22}{chosen_on[:5]:<12}⛔ nothing clears {TARGET:.0%} "
                    f"with ≥{MIN_CALLED} scored; best {best[1][1] or 0:.1%}"
                )
                continue
            policy, scored, rate = max(viable, key=lambda item: item[1])
            held_scored, held_rate = grids[held_out][policy]
            rule, esm_bar, bacformer_bar = policy
            shown = f"{rule} esm≥{esm_bar or 'any'} bac≥{bacformer_bar or 'any'}"
            print(
                f"  {label:<22}{chosen_on[:5]:<12}{shown:<30}{_pct(rate)}{_pct(held_rate)}{held_scored:>9,}"
            )


def score_on(state: State, assigned: dict[int, frozenset[str]], subset: set[int]) -> tuple[int, float | None]:
    """`(scored, precision)` restricted to `subset` — the half-split selection's scorer."""
    scored = right = 0
    for node, claim in assigned.items():
        if node not in subset:
            continue
        scored += 1
        right += bool(claim & state.truth_of(node))
    return scored, (right / scored if scored else None)


def report_extra(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
) -> None:
    """S5 — extrapolation alone, ecoli → kp, scored THREE ways: ESM, Bacformer, and both."""
    print(f"\n{'=' * 118}\nS5 — EXTRAPOLATION ALONE (ecoli → kp, one-directional)\n{'=' * 118}")
    print("  ⛔ ONE-DIRECTIONAL. The only map on disk runs ecoli → kp, so kp can receive and ecoli cannot.")
    print("  ⛔ There is no OTHER catalogue to choose a bar on, so the masked kp nodes are split in half:")
    print("     the bar is chosen on half A and the number that counts is half B. Both are printed.")
    print("  ⚠ Cross-species cosines are a third distribution — the ESM map's rank-1 pairs sit near 0.976")
    print("    and Bacformer's near 0.35 — so neither within-species ladder transfers.\n")
    print(
        f"  {'claim':<22}{'rule':<11}{'bar chosen on A':<22}{'A prec':>8}"
        f"{'B prec':>8}{'B scored':>10}{'called':>9}{'node cover':>12}"
    )
    for key, label, _kind, _rung in CLAIMS:
        recipient, donor = problems["kp-nuna4"][key], problems["ecoli-nuna4"][key]
        recipient_base, donor_base = states["kp-nuna4"][key], states["ecoli-nuna4"][key]
        halves = sorted(recipient_base.masked)
        half_a = set(halves[::2])
        half_b = set(halves[1::2])
        for rule in ("esm", "bacformer", "concur"):
            grid = {}
            for esm_bar in CROSS_BARS["ESM"] if rule != "bacformer" else (None,):
                for bacformer_bar in CROSS_BARS["BACFORMER"] if rule != "esm" else (None,):
                    state = recipient_base.copy()
                    assigned = extrapolate(
                        recipient, state, donor, donor_base, maps,
                        {"ESM": esm_bar, "BACFORMER": bacformer_bar}, rule=rule,
                    )
                    a_scored, a_rate = score_on(state, assigned, half_a)
                    b_scored, b_rate = score_on(state, assigned, half_b)
                    grid[(esm_bar, bacformer_bar)] = (a_scored, a_rate, b_scored, b_rate, len(assigned))
            viable = [
                (bars, values) for bars, values in grid.items()
                if values[1] is not None and values[1] >= TARGET and values[0] >= MIN_CALLED / 2
            ]
            if not viable:
                best_bars, best = max(grid.items(), key=lambda item: item[1][1] or 0)
                shown = f"esm≥{best_bars[0] or 'any'} bac≥{best_bars[1] or 'any'}"
                print(
                    f"  {label:<22}{rule:<11}⛔ none clears {TARGET:.0%} on A; best {shown:<11}"
                    f"{_pct(best[1])}{_pct(best[3])}{best[2]:>10,}{best[4]:>9,}"
                )
                continue
            bars, (a_scored, a_rate, b_scored, b_rate, called) = max(viable, key=lambda item: item[1][4])
            shown = f"esm≥{bars[0] or 'any'} bac≥{bars[1] or 'any'}"
            cover = (len(recipient_base.label) + called) / recipient.node_total
            print(
                f"  {label:<22}{rule:<11}{shown:<22}{_pct(a_rate)}{_pct(b_rate)}{b_scored:>10,}"
                f"{called:>9,}{cover:>12.1%}"
            )


def report_cross_reach(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
) -> None:
    """⛔ The coverage a cross-vs-within comparison must assert BEFORE it reports a difference.

    `~/.claude/CLAUDE.md`: *a comparison must assert its own coverage before it reports a
    difference.* Extrapolation reaching more nodes than interpolation is only a statement about the
    mechanism if both are looking at the same number of candidates AND the candidate pools are
    equally well annotated. This prints the second half — how often a donor shortlist of the same
    depth contains an annotated node at all, within kp and from ecoli.
    """
    print("\n  ⛔ REACH, side by side — can a donor shortlist of this depth find an annotated node at all?")
    print(
        f"    {'claim':<22}{'kp unlabelled':>14}{'within kp':>11}{'from ecoli':>12}"
        f"{'within %':>10}{'cross %':>9}"
    )
    donor_id_of = {label: node for node, label in problems["ecoli-nuna4"]["cog"].node_label.items()}
    for key, label, _kind, _rung in CLAIMS:
        recipient, state = problems["kp-nuna4"][key], states["kp-nuna4"][key]
        donor_state = states["ecoli-nuna4"][key]
        unlabelled = [node for node in recipient.gene_count if node not in state.label]
        within = sum(nearest_labelled(recipient, state, node, "ESM") is not None for node in unlabelled)
        across = 0
        for node in unlabelled:
            for donor_label, _cosine, _rank in maps["ESM"].get(recipient.node_label[node], ()):
                donor_node = donor_id_of.get(donor_label)
                if donor_node is not None and donor_node in donor_state.label:
                    across += 1
                    break
        print(
            f"    {label:<22}{len(unlabelled):>14,}{within:>11,}{across:>12,}"
            f"{within / len(unlabelled):>10.1%}{across / len(unlabelled):>9.1%}"
        )
    ecoli, kp = problems["ecoli-nuna4"], problems["kp-nuna4"]
    print("    ⚠ The donor pools are not equally annotated, and that is a candidate explanation:")
    for key, label, _kind, _rung in CLAIMS:
        print(
            f"      {label:<22}ecoli {len(ecoli[key].truth) / ecoli[key].node_total:>5.1%} annotated"
            f"  ·  kp {len(kp[key].truth) / kp[key].node_total:>5.1%}"
        )


#: ⭐⭐ **Certainty is TWO-DIMENSIONAL, and that was measured rather than assumed.**
#:
#: David, 2026-10-11: *"The bacformer certainty must help give accuracy at least … a similarity of
#: 0.8 with esm of 0.99 must be better than bacformer 0.5 and esm 0.99."* He is right, and an earlier
#: version of this file dismissed Bacformer on a bad argument — that David's ESM band edges leave
#: every Bacformer bin empty, which is true and irrelevant. The question is whether Bacformer orders
#: accuracy **at fixed ESM**, which needs no shared scale at all. Measured over 4,620 held-out calls
#: carrying both cosines, within each ESM band, lowest Bacformer third → highest:
#:
#:   ESM >= 0.99     95.8 % → 96.8 %   (+1.0)   — saturated; Bacformer adds nothing
#:   ESM 0.95-0.99   76.8 % → 94.5 %  (+17.7)   — ⭐ the band where 45 % of all calls live
#:   ESM 0.90-0.95   61.9 % → 86.7 %  (+24.8)
#:   ESM < 0.90      42.9 % → 75.0 %  (+32.1)   — n=22, directional only
#:
#: ⛔ **The bar is 0.60 and it was chosen on one catalogue and scored on the other**, which is the
#: one thing the earlier interpolation policies did not survive. Within ESM 0.95-0.99 it reads
#: 92.7 % (ecoli) / 94.6 % (kp held out) above the bar against 85.6 % / 79.8 % below it, and every
#: bar from 0.40 to 0.75 tracks between the two catalogues within ~2 points.
BACFORMER_BAR = 0.60

#: The ESM-only bands the conditional test above is tabulated over. ⚠ NOT a certainty scheme on
#: their own — that is the point of `CERTAINTY_TIERS`.
ESM_BANDS: tuple[tuple[str, float, float], ...] = (
    (">= 0.99", 0.99, 1.01),
    ("0.95-0.99", 0.95, 0.99),
    ("0.90-0.95", 0.90, 0.95),
    ("< 0.90", -1.01, 0.90),
)

#: The six tiers, strongest first. ⭐ **D ≈ C**, which is the finding in one line: a high Bacformer
#: cosine buys back roughly one whole ESM band.
CERTAINTY_TIERS: tuple[str, ...] = (
    "A · ESM ≥ 0.99",
    "B · ESM 0.95-0.99 · Bac ≥ 0.60",
    "C · ESM 0.95-0.99 · Bac < 0.60",
    "D · ESM 0.90-0.95 · Bac ≥ 0.60",
    "E · ESM 0.90-0.95 · Bac < 0.60",
    "F · ESM < 0.90",
)
#: Which tiers clear the 90 % bar on held-out truth, with the interval to back it:
#: A 96.7 % (95.9-97.4) and B 94.0 % (92.3-95.4). C 85.0 %, D 84.6 %, E 79.6 %, F 73.9 % do not.
CALLABLE_TIERS: frozenset[str] = frozenset(CERTAINTY_TIERS[:2])


def bacformer_certainty(signals: dict) -> float | None:
    """The Bacformer cosine a call's certainty is read from — same pair first, own donor as fallback.

    ⛔ **The same pair is the better axis and it is missing a quarter of the time.** It separates
    accuracy by +17.7 points inside ESM 0.95-0.99 against the own-donor axis's +13.0, so it is
    preferred — but the within-species walk stores only five neighbours, so the ESM donor is often
    absent from the Bacformer shortlist. Falling back keeps every call in a tier (the fallback covers
    24 % of calls at 88.8 % accuracy against the same-pair 90.1 %) rather than dropping a quarter of
    them into an unmeasured bin, which would silently make every tier look better than it is.
    """
    pair = signals.get("pair_bacformer_cosine")
    return pair if pair is not None else signals.get("bacformer_cosine")


def tier_of(signals: dict) -> str:
    """Which certainty tier a call falls in, from both cosines."""
    esm = signals.get("esm_cosine")
    if esm is None or esm < 0.90:
        return CERTAINTY_TIERS[5]
    bacformer = bacformer_certainty(signals)
    strong = bacformer is not None and bacformer >= BACFORMER_BAR
    if esm >= 0.99:
        return CERTAINTY_TIERS[0]
    if esm >= 0.95:
        return CERTAINTY_TIERS[1] if strong else CERTAINTY_TIERS[2]
    return CERTAINTY_TIERS[3] if strong else CERTAINTY_TIERS[4]


def _stratify(
    problem: Problem, state: State, assigned: dict[int, frozenset[str]]
) -> list[dict]:
    """One row per certainty band: what it gained, and how right it was.

    ⭐ David's Q4: *"some will be extremely high and others lower"*. A policy reported as one rate is
    an average over bands that differ, and the average is what the page currently quotes.

    ⛔ A band's precision is `None` below `MIN_SCORED_TO_RETIRE` held-out nodes — never 0.0, which
    would read as *measured and bad* rather than *not measured*. Its gain is still counted, because
    the nodes are gained whether or not the band happens to hold enough masked nodes to check.
    """
    grouped: dict[str, list[int]] = {name: [] for name in CERTAINTY_TIERS}
    for node in assigned:
        grouped[tier_of(state.donor.get(node, {}))].append(node)
    rows = []
    for name in CERTAINTY_TIERS:
        nodes = grouped[name]
        subset = {node: assigned[node] for node in nodes}
        _called, scored, rate = precision(state, subset)
        rows.append(
            {
                "band": name,
                "callable": name in CALLABLE_TIERS,
                "nodes": len(nodes),
                "genes": sum(problem.gene_count[node] for node in nodes),
                "scored": scored,
                "precision": rate if scored >= MIN_SCORED_TO_RETIRE else None,
                "precision_unguarded": rate,
                # ⛔ **The SAME guard, and it was missing here.** `precision` was suppressed below
                # 30 held-out nodes while the gene-weighted rate beside it was not, so the table
                # printed `—   0.0%` off a single node — an honest dash beside a measured-looking
                # zero, in the one column its own preamble forbids. 18 of 72 strata rows were
                # affected; two rested on one node each.
                "precision_by_gene": (
                    precision_by_gene(problem, state, subset)
                    if scored >= MIN_SCORED_TO_RETIRE
                    else None
                ),
                "precision_by_gene_unguarded": precision_by_gene(problem, state, subset),
            }
        )
    return rows


def measure_one_pass(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    cross_bars: dict[str, float | None],
    *,
    mask_fraction: float,
) -> dict:
    """⭐ THE deliverable: one pass — inter, then extra, with intra applied continuously.

    ⛔⛔ **TWO RUNS, because coverage and precision cannot honestly come from the same one.**

    * **the coverage run** holds NO mask. Every true label is available as a donor, which is the real
      deployment, and the nodes it labels are the genuinely unannotated ones. Its node and gene
      counts are the achievable coverage.
    * **the precision run** holds the mask. It is the only run with held-out truth, so it is the only
      one that can state a rate.

    ⛔ **Reading coverage off the masked run would understate it, and badly.** With 20 % of ecoli's
    COG labels withheld, the pass recovered 224 of 999 masked nodes and the other 775 — which are
    large, ~68 genes each — simply vanish from the census, so gene coverage read 58.8 % against a
    69.5 % baseline. A figure drawn from that would show one pass LOSING coverage.

    ⚠ So every stratum carries its gain from one run and its rate from the other, and says so. The
    two runs label different node sets, so a band's `nodes`/`genes` and its `precision` are not
    numerator and denominator of each other and must never be divided.

    ecoli gets interpolation only; there is no donor species for it.
    """
    record: dict = {
        # ⛔ **The fraction, not just the seed.** It changes every rate in this record — kp
        # round-1 extrapolation reads 93.1 % at 0.20 and 84.9 % at 0.50 — and a record that names
        # only the seed cannot be told apart from a run at another fraction.
        "mask_fraction": mask_fraction,
        "mask_seed": MASK_SEED,
        "certainty_bands": list(CERTAINTY_TIERS),
        "callable_bands": sorted(CALLABLE_TIERS),
        "certainty_axis": (
            "ESM cosine of the donor used, split by the Bacformer cosine of the same pair "
            f"(own-donor cosine where the shortlist lacks it) at {BACFORMER_BAR}"
        ),
        "bacformer_bar": BACFORMER_BAR,
        "gain_from": "the UNMASKED coverage run — every true label available as a donor",
        "rates_from": "the MASKED precision run — the only run with held-out truth",
        # ⚠ The two runs also differ in DONOR QUALITY, not only in which nodes they call, because
        # ecoli is interpolated before kp is extrapolated into in both. Giving the precision run the
        # coverage run's ecoli donors raises kp extrapolation precision by 6.5 / 5.7 / 3.4 / 4.2
        # points (COG / MF / BP / CC) at mask 0.50 and 1.5 / 1.2 / 1.1 / 1.3 at 0.20. The published
        # rate is therefore CONSERVATIVE against the pool its gain was measured under.
        "donor_pool_note": (
            "extrapolation's rate is measured against a donor pool poorer than the one its gain is "
            "measured under, so it understates by ~1.1-1.5 points at mask 0.20 and ~3.4-6.5 at 0.50"
        ),
        "catalogues": {},
    }
    for catalogue in CATALOGUES:
        problem_any = problems[catalogue]["cog"]
        record["catalogues"][catalogue] = {
            "node_total": problem_any.node_total,
            "gene_total": problem_any.gene_total,
            "claims": {},
        }

    def one_pass(key: str, masked: bool) -> dict[str, tuple[Problem, State, dict[str, dict]]]:
        """Interpolate both species, then extrapolate into kp. Returns per-catalogue state."""
        working = {
            catalogue: (
                states[catalogue][key].copy()
                if masked
                else State(dict(problems[catalogue][key].truth), {}, {})
            )
            for catalogue in CATALOGUES
        }
        out: dict[str, tuple[Problem, State, dict[str, dict]]] = {}
        for catalogue in CATALOGUES:
            problem, state = problems[catalogue][key], working[catalogue]
            assigned = interpolate(problem, state, INTERPOLATION_POLICY[key])
            _called, scored, rate = precision(state, assigned)
            out[catalogue] = (
                problem,
                state,
                {
                    "inter": {
                        "called": len(assigned),
                        "scored": scored,
                        "precision": rate,
                        "precision_by_gene": precision_by_gene(problem, state, assigned),
                        "strata": _stratify(problem, state, assigned),
                    }
                },
            )
            apply(state, assigned, "inter", 1)
        recipient, donor = problems["kp-nuna4"][key], problems["ecoli-nuna4"][key]
        state = working["kp-nuna4"]
        assigned = extrapolate(recipient, state, donor, working["ecoli-nuna4"], maps, cross_bars)
        _called, scored, rate = precision(state, assigned)
        out["kp-nuna4"][2]["extra"] = {
            "called": len(assigned),
            "scored": scored,
            "precision": rate,
            "precision_by_gene": precision_by_gene(recipient, state, assigned),
            "strata": _stratify(recipient, state, assigned),
        }
        apply(state, assigned, "extra", 1)
        return out

    for key, label, _kind, _rung in CLAIMS:
        coverage = one_pass(key, masked=False)
        scored_run = one_pass(key, masked=True)
        for catalogue in CATALOGUES:
            problem, state, gains = coverage[catalogue]
            _p, _s, rates = scored_run[catalogue]
            whole = State(dict(problem.truth), {}, {})
            mechanisms: dict[str, dict] = {}
            for mechanism, gained in gains.items():
                by_band = {row["band"]: row for row in rates[mechanism]["strata"]}
                mechanisms[mechanism] = {
                    "called": gained["called"],
                    "precision": rates[mechanism]["precision"],
                    "precision_by_gene": rates[mechanism]["precision_by_gene"],
                    "scored": rates[mechanism]["scored"],
                    "strata": [
                        {
                            "band": row["band"],
                            "nodes": row["nodes"],
                            "genes": row["genes"],
                            "scored": by_band[row["band"]]["scored"],
                            "precision": by_band[row["band"]]["precision"],
                            "precision_by_gene": by_band[row["band"]]["precision_by_gene"],
                        }
                        for row in gained["strata"]
                    ],
                }
            record["catalogues"][catalogue]["claims"][key] = {
                "label": label,
                "baseline": census(problem, whole),
                "after_one_pass": census(problem, state),
                "masked": len(states[catalogue][key].masked),
                "mechanisms": mechanisms,
            }
    return record


def scoreable_calls(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    cross_bars: dict[str, float | None],
) -> list[dict]:
    """⭐⭐ One row per HELD-OUT call, with every certainty signal it had and whether it was right.

    This is the thing to emit rather than another aggregate. David, 2026-10-11: *"We want a way of
    estimating the actual accuracy"* — which needs the joint distribution of the two cosines against
    correctness, not a table already collapsed over one of them. Written once; every stratification,
    every band edge and any calibration is then a few seconds of work over the CSV instead of
    another full run.

    ⛔ Only masked nodes appear. A call on a node that was never annotated has no truth and cannot
    be a row here, which is measurement rule 1 again.
    """
    rows: list[dict] = []
    for key, label, _kind, _rung in CLAIMS:
        working = {catalogue: states[catalogue][key].copy() for catalogue in CATALOGUES}
        passes: list[tuple[str, str, Problem, State, dict[int, frozenset[str]]]] = []
        for catalogue in CATALOGUES:
            problem, state = problems[catalogue][key], working[catalogue]
            assigned = interpolate(problem, state, INTERPOLATION_POLICY[key])
            passes.append((catalogue, "inter", problem, state, assigned))
        for _catalogue, _mechanism, _problem, state, assigned in list(passes):
            apply(state, assigned, "inter", 1)
        recipient = problems["kp-nuna4"][key]
        state = working["kp-nuna4"]
        assigned = extrapolate(
            recipient, state, problems["ecoli-nuna4"][key], working["ecoli-nuna4"], maps, cross_bars
        )
        passes.append(("kp-nuna4", "extra", recipient, state, assigned))

        for catalogue, mechanism, problem, state, assigned in passes:
            for node, claim in assigned.items():
                if node not in state.masked:
                    continue
                signals = state.donor.get(node, {})
                rows.append(
                    {
                        "catalogue": catalogue,
                        "claim": key,
                        "claim_label": label,
                        "mechanism": mechanism,
                        "node": node,
                        "node_label": problem.node_label[node],
                        "gene_count": problem.gene_count[node],
                        "esm_cosine": signals.get("esm_cosine"),
                        "bacformer_cosine": signals.get("bacformer_cosine"),
                        "pair_bacformer_cosine": signals.get("pair_bacformer_cosine"),
                        "donor": signals.get("donor"),
                        "correct": int(bool(claim & state.truth_of(node))),
                    }
                )
    return rows


def report_certainty_surface(rows: list[dict]) -> None:
    """⭐ Does Bacformer add accuracy AT FIXED ESM? David's question, as a 2×N table.

    ⛔ **The test has to be conditional**, because the two cosines are correlated: a marginal table
    showing *"high Bacformer is more accurate"* would be satisfied entirely by high Bacformer pairs
    also having high ESM. Within an ESM band, Bacformer either orders accuracy or it does not.

    ⚠ Pooled over both catalogues and both mechanisms, because 16 cells over a few hundred held-out
    calls leaves nothing measurable otherwise. Every cell prints its n.
    """
    print(f"\n{'=' * 118}\nS9 — DOES BACFORMER ADD CERTAINTY AT FIXED ESM?\n{'=' * 118}")
    print("  Pooled over catalogues and mechanisms. Bacformer terciles are computed WITHIN each ESM")
    print("  band, so each row's three cells hold a third of that band's calls by construction.")
    print("  ⚠ 'same pair' uses the Bacformer cosine of the SAME donor pair where the shortlist holds")
    print("    it; 'own donor' uses the cosine to the donor Bacformer itself chose, always available.\n")
    for axis_name, column in (("own donor", "bacformer_cosine"), ("same pair", "pair_bacformer_cosine")):
        usable = [row for row in rows if row[column] is not None and row["esm_cosine"] is not None]
        print(
            f"  Bacformer axis = {axis_name}: {len(usable):,} of {len(rows):,} held-out calls carry it"
            f" ({len(usable) / len(rows):.0%})"
        )
        print(
            f"    {'ESM band':<12}{'n':>7}{'overall':>9}"
            f"{'low third':>12}{'n':>6}{'mid third':>12}{'n':>6}{'top third':>12}{'n':>6}{'top - low':>11}"
        )
        for name, low, high in ESM_BANDS:
            inside = [row for row in usable if low <= (row["esm_cosine"] or -1) < high]
            if not inside:
                print(f"    {name:<12}{0:>7}      no calls")
                continue
            ordered = sorted(inside, key=lambda row: row[column])
            cut = len(ordered) // 3
            thirds = [ordered[:cut], ordered[cut : 2 * cut], ordered[2 * cut :]] if cut else [[], [], ordered]
            overall = sum(row["correct"] for row in inside) / len(inside)
            cells = [
                (sum(row["correct"] for row in third) / len(third) if third else None, len(third))
                for third in thirds
            ]
            delta = (
                cells[2][0] - cells[0][0] if cells[2][0] is not None and cells[0][0] is not None else None
            )
            print(
                f"    {name:<12}{len(inside):>7,}{overall:>9.1%}"
                + "".join(f"{_pct(rate):>12}{n:>6,}" for rate, n in cells)
                + (f"{delta:>+11.1%}" if delta is not None else f"{'—':>11}")
            )
        print()


def report_strata(record: dict) -> None:
    """S8 — the gain stratified by certainty, which is what the headline figures are drawn from."""
    print(f"\n{'=' * 118}\nS8 — ANNOTATION GAIN, BY PREDICTION, STRATIFIED BY CERTAINTY\n{'=' * 118}")
    print("  Certainty is TWO-dimensional: the ESM cosine of the donor used, split by the Bacformer")
    print(f"  cosine of the same donor pair at {BACFORMER_BAR}. ⚠ Bacformer orders accuracy BELOW ESM 0.99")
    print("  and not at or above it, where the claim is already saturated at ~96.7 %.")
    print("  ⛔ 'nodes +' and 'genes +' come from the UNMASKED coverage run; 'prec' from the MASKED one.")
    print("     They label different node sets on purpose, so never divide one by the other.")
    print(f"  ⛔ A band's rate is — below {MIN_SCORED_TO_RETIRE} held-out nodes: not measured, not zero.\n")
    for catalogue, body in record["catalogues"].items():
        for claim in body["claims"].values():
            print(f"  {catalogue} · {claim['label']}")
            print(
                f"    {'mechanism':<10}{'tier':<32}{'nodes +':>9}{'genes +':>10}{'scored':>8}"
                f"{'prec':>8}{'by gene':>9}"
            )
            for mechanism, measured in claim["mechanisms"].items():
                for row in measured["strata"]:
                    print(
                        f"    {mechanism:<10}{row['band']:<32}{row['nodes']:>9,}{row['genes']:>10,}"
                        f"{row['scored']:>8,}{_pct(row['precision'])}{_pct(row['precision_by_gene'])}"
                    )
                print(
                    f"    {mechanism:<10}{'ALL':<32}{measured['called']:>9,}"
                    f"{sum(r['genes'] for r in measured['strata']):>10,}{measured['scored']:>8,}"
                    f"{_pct(measured['precision'])}{_pct(measured['precision_by_gene'])}"
                )
            before, after = claim["baseline"], claim["after_one_pass"]
            print(
                f"    → nodes {before['nodes']:,} ({before['node_share']:.1%})"
                f" → {after['nodes']:,} ({after['node_share']:.1%})"
                f"  ·  genes {before['genes_carrying_share']:.1%} carrying"
                f" → {before['genes_covered_share']:.1%} after intra"
                f" → {after['genes_covered_share']:.1%} after one pass"
            )


def _after_inter(problem: Problem, state: State, claim_key: str) -> State:
    """A detached copy of `state` with one round of interpolation applied — the 2×2's treatment."""
    advanced = state.copy()
    apply(advanced, interpolate(problem, advanced, INTERPOLATION_POLICY[claim_key]), "inter", 1)
    return advanced


def report_interaction(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    cross_bars: dict[str, float | None],
) -> None:
    """⭐ David's Q1: how much does interpolation change what extrapolation can do?

    A 2×2, because *"interpolation first"* is TWO different interventions and they act on opposite
    sides of the transfer:

    * **donor side** — interpolating *ecoli* gives kp more donors to read from, and those donors
      carry propagated labels rather than true ones;
    * **recipient side** — interpolating *kp* takes nodes off the table before extrapolation sees
      them, so extra is left with whatever interpolation could not reach.

    ⛔ **The four cells are not comparable on `called` alone**, and that is the trap this table is
    shaped to avoid: a recipient that has already been interpolated has fewer unlabelled nodes, so
    extra's `called` falls for a reason that has nothing to do with the donors. The last column is
    therefore precision restricted to the masked nodes **every cell calls**, which is the only
    population all four can be read against.
    """
    print(f"\n{'=' * 118}\nS5b — DOES INTERPOLATION HELP EXTRAPOLATION? (the 2×2)\n{'=' * 118}")
    print("  Recipient is always kp. Donor is always ecoli. Scored on kp's masked held-out nodes.")
    print("  ⛔ 'common' is precision over the masked nodes ALL FOUR cells call — the only fair column")
    print("     when the cells call different numbers of nodes for structural reasons.\n")
    for key, label, _kind, _rung in CLAIMS:
        recipient, donor = problems["kp-nuna4"][key], problems["ecoli-nuna4"][key]
        cells: dict[tuple[bool, bool], tuple[dict[int, frozenset[str]], State]] = {}
        for donor_interpolated in (False, True):
            donor_state = (
                _after_inter(donor, states["ecoli-nuna4"][key], key)
                if donor_interpolated
                else states["ecoli-nuna4"][key].copy()
            )
            for recipient_interpolated in (False, True):
                recipient_state = (
                    _after_inter(recipient, states["kp-nuna4"][key], key)
                    if recipient_interpolated
                    else states["kp-nuna4"][key].copy()
                )
                assigned = extrapolate(
                    recipient, recipient_state, donor, donor_state, maps, cross_bars
                )
                cells[(donor_interpolated, recipient_interpolated)] = (assigned, recipient_state)
        common = set.intersection(
            *(
                {node for node in assigned if node in state.masked}
                for assigned, state in cells.values()
            )
        )
        print(f"  {label}")
        print(
            f"    {'donor ecoli':<14}{'recipient kp':<16}{'called':>8}{'scored':>8}{'prec':>8}"
            f"{'by gene':>9}{'common':>9}{'n common':>10}{'kp nodes':>10}"
        )
        for (donor_interpolated, recipient_interpolated), (assigned, state) in cells.items():
            _called, scored, rate = precision(state, assigned)
            by_gene = precision_by_gene(problem=recipient, state=state, assigned=assigned)
            restricted = {node: assigned[node] for node in common}
            _c, _s, common_rate = precision(state, restricted)
            cover = (len(state.label) + len(assigned)) / recipient.node_total
            print(
                f"    {'+ inter' if donor_interpolated else 'truth only':<14}"
                f"{'+ inter' if recipient_interpolated else 'truth only':<16}"
                f"{len(assigned):>8,}{scored:>8,}{_pct(rate)}{_pct(by_gene)}{_pct(common_rate)}"
                f"{len(common):>10,}{cover:>10.1%}"
            )


def report_concordance(
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    cross_bars: dict[str, float | None],
) -> None:
    """⭐ David's Q2: where interpolation and extrapolation both call a kp node, do they agree?

    ⭐⭐ **This is the one accuracy measure here that needs no held-out truth at all**, which makes it
    the only one that could be run on the nodes we actually want to annotate rather than on annotated
    nodes standing in for them. Two mechanisms reading DIFFERENT donor pools — one inside kp, one
    inside ecoli — landing on the same claim is evidence about the claim.

    ⚠ **And it is genuinely two predictions, unlike the ESM/Bacformer concurrence test.** That one
    compared two geometries over the SAME donor pool, and because a COG claim is a single accession
    *"the claims intersect"* collapsed to *"the claims are equal"*; 19-28 % of nodes even got the
    identical donor from both maps. Here the donors are in different species and cannot coincide.
    ⚠ They are still correlated through the representation: both read ESM cosines.

    The payoff column is the last pair — precision WHEN THEY AGREE against precision when they do
    not. If agreement separates them, it is a filter that costs reach and buys precision without
    needing a threshold on either side.
    """
    print(f"\n{'=' * 118}\nS5c — DO INTERPOLATION AND EXTRAPOLATION AGREE? (kp, both from TRUTH donors)\n{'=' * 118}")
    print("  Both mechanisms run independently on the same kp state, each reading only TRUE labels.")
    print("  ⭐ Needs no held-out truth to compute, so it is the one measure that transfers to the")
    print("     nodes we actually want to annotate. Precision columns still use the masked set.\n")
    print(
        f"  {'claim':<22}{'inter':>8}{'extra':>8}{'both':>8}{'agree':>8}{'rate':>8}"
        f"{'prec|agree':>12}{'n':>6}{'prec|differ':>13}{'n':>6}"
    )
    rows: list[tuple[str, dict]] = []
    for key, label, _kind, _rung in CLAIMS:
        recipient, donor = problems["kp-nuna4"][key], problems["ecoli-nuna4"][key]
        state = states["kp-nuna4"][key]
        inter_calls = interpolate(recipient, state.copy(), INTERPOLATION_POLICY[key])
        extra_calls = extrapolate(
            recipient, state.copy(), donor, states["ecoli-nuna4"][key], maps, cross_bars
        )
        both = set(inter_calls) & set(extra_calls)
        agreeing = {node for node in both if inter_calls[node] & extra_calls[node]}
        scored_agree = {node: inter_calls[node] for node in agreeing if node in state.masked}
        scored_differ = {
            node: inter_calls[node] for node in (both - agreeing) if node in state.masked
        }
        _c, n_agree, rate_agree = precision(state, scored_agree)
        _c, n_differ, rate_differ = precision(state, scored_differ)
        print(
            f"  {label:<22}{len(inter_calls):>8,}{len(extra_calls):>8,}{len(both):>8,}"
            f"{len(agreeing):>8,}{len(agreeing) / len(both) if both else 0:>8.1%}"
            f"{_pct(rate_agree)}{n_agree:>6,}{_pct(rate_differ)}{n_differ:>6,}"
        )
        rows.append((label, {"inter": inter_calls, "extra": extra_calls, "both": both}))

    print("\n  ⚠ What each mechanism reaches ALONE, and what agreement costs in reach:")
    print(
        f"  {'claim':<22}{'inter only':>12}{'extra only':>12}{'both':>8}{'agreeing':>10}"
        f"{'union':>8}{'agree/union':>13}"
    )
    for label, parts in rows:
        inter_only = len(set(parts["inter"]) - parts["both"])
        extra_only = len(set(parts["extra"]) - parts["both"])
        union = len(set(parts["inter"]) | set(parts["extra"]))
        agreeing = sum(
            1 for node in parts["both"] if parts["inter"][node] & parts["extra"][node]
        )
        print(
            f"  {label:<22}{inter_only:>12,}{extra_only:>12,}{len(parts['both']):>8,}"
            f"{agreeing:>10,}{union:>8,}{agreeing / union if union else 0:>13.1%}"
        )


def run_loop(
    claim_key: str,
    problems: dict[str, dict[str, Problem]],
    states: dict[str, dict[str, State]],
    maps: dict[str, dict[str, list[tuple[str, float, int]]]],
    cross_bars: dict[str, float | None],
    *,
    pre_uniref: bool = False,
    post_uniref: bool = False,
) -> list[dict]:
    """Inter → extra → (intra, continuously) → inter → … until it stops adding or stops being right.

    ⭐ **Intrapolation is NOT a step, and that is a finding rather than an omission.** It cannot
    annotate a node that had nothing, so it never changes the donor pool — it converts a node label
    into gene coverage, which the `genes_covered` column already carries at every round. The loop
    is therefore inter → extra → inter → extra, with intra applied continuously and for free.

    ⚠ **ecoli runs inter → inter only.** There is no donor species for it: the map is one-directional.
    """
    working = {catalogue: states[catalogue][claim_key].copy() for catalogue in CATALOGUES}
    rows: list[dict] = []
    #: (catalogue, mechanism) pairs whose last measured round fell below the bar — they stop running.
    retired: set[tuple[str, str]] = set()

    def record(catalogue: str, mechanism: str, round_number: int, assigned: dict[int, frozenset[str]]) -> int:
        problem, state = problems[catalogue][claim_key], working[catalogue]
        called, scored, rate = precision(state, assigned)
        by_gene = precision_by_gene(problem, state, assigned)
        apply(state, assigned, mechanism, round_number)
        numbers = census(problem, state)
        rows.append(
            {
                "catalogue": catalogue,
                "round": round_number,
                "mechanism": mechanism,
                "new": called,
                "scored": scored,
                "precision": rate,
                "precision_by_gene": by_gene,
                **numbers,
            }
        )
        return called

    if pre_uniref:
        for catalogue in CATALOGUES:
            record(catalogue, "uniref50", 0, uniref_propagate(problems[catalogue][claim_key], working[catalogue]))

    for round_number in range(1, MAX_ROUNDS + 1):
        added = 0
        for catalogue in CATALOGUES:
            if (catalogue, "inter") in retired:
                continue
            added += record(
                catalogue, "inter", round_number,
                interpolate(problems[catalogue][claim_key], working[catalogue], INTERPOLATION_POLICY[claim_key]),
            )
        if ("kp-nuna4", "extra") not in retired:
            added += record(
                "kp-nuna4", "extra", round_number,
                extrapolate(
                    problems["kp-nuna4"][claim_key], working["kp-nuna4"],
                    problems["ecoli-nuna4"][claim_key], working["ecoli-nuna4"],
                    maps, cross_bars,
                ),
            )
        # ⛔ A mechanism is retired INDEPENDENTLY, and only on a rate that is measurable. Stopping
        # the whole loop because one mechanism degraded would discard the rounds another was still
        # earning; retiring on a rate off four held-out nodes would stop it on noise.
        for row in list(rows):
            if row["round"] != round_number or row["precision"] is None:
                continue
            if row["scored"] >= MIN_SCORED_TO_RETIRE and row["precision"] < TARGET:
                retired.add((row["catalogue"], row["mechanism"]))
                rows.append(
                    {
                        **row,
                        "mechanism": f"⛔ {row['mechanism']} retired",
                        "new": 0,
                        "scored": 0,
                        "precision": None,
                        "precision_by_gene": None,
                    }
                )
        if added == 0:
            break
        if len(retired) >= 3:
            rows.append(
                {
                    "catalogue": "—",
                    "round": round_number,
                    "mechanism": "⛔ STOPPED: every mechanism fell below the bar and retired",
                    "new": 0,
                    "scored": 0,
                    "precision": None,
                    "precision_by_gene": None,
                    "nodes": 0,
                    "node_share": 0.0,
                    "genes_carrying": 0,
                    "genes_carrying_share": 0.0,
                    "genes_covered": 0,
                    "genes_covered_share": 0.0,
                }
            )
            break

    if post_uniref:
        # ⭐ The marginal gain here is the whole point of the arm: what UniRef50 knows that ESM and
        # Bacformer did not already find.
        for catalogue in CATALOGUES:
            record(
                catalogue, "uniref50", MAX_ROUNDS + 1,
                uniref_propagate(problems[catalogue][claim_key], working[catalogue]),
            )
    return rows


def print_rounds(title: str, rows: list[dict], problems: dict[str, dict[str, Problem]], claim_key: str) -> None:
    """One table: what each round added, how right it was, and where coverage stood after it."""
    print(f"\n  {title}")
    print(
        f"    {'cat':<7}{'rd':>3}{'mechanism':<12}{'new':>8}{'scored':>8}{'prec':>8}"
        f"{'by gene':>9}{'nodes':>9}{'of all':>8}{'genes':>11}{'of all':>8}"
    )
    for row in rows:
        if row["new"] == 0 and row["mechanism"].startswith("⛔"):
            print(f"    {row['mechanism']}")
            continue
        print(
            f"    {row['catalogue'][:5]:<7}{row['round']:>3}{row['mechanism']:<12}{row['new']:>8,}{row['scored']:>8,}"
            f"{_pct(row['precision'])}{_pct(row['precision_by_gene'])}"
            f"{row['nodes']:>9,}{row['node_share']:>8.1%}{row['genes_covered']:>11,}{row['genes_covered_share']:>8.1%}"
        )
    for catalogue in CATALOGUES:
        problem = problems[catalogue][claim_key]
        mine = [row for row in rows if row["catalogue"] == catalogue]
        if not mine:
            continue
        final = mine[-1]
        truth_share = len(problem.truth) / problem.node_total
        print(
            f"    {catalogue[:5]} END  nodes {final['nodes']:,} ({final['node_share']:.1%}) \u00b7 "
            f"genes {final['genes_covered']:,} ({final['genes_covered_share']:.1%})"
            f"  |  true baseline {len(problem.truth):,} nodes ({truth_share:.1%}), of which "
            f"{len(problem.truth) - mine[0]['nodes'] + mine[0]['new']:,} were masked out of the donor pool"
        )


def report_reach(problems: dict[str, dict[str, Problem]], states: dict[str, dict[str, State]]) -> None:
    """⛔ The binding constraint, measured before any policy is applied: can the walk REACH a donor?

    `locus_nearest_locus` stores five neighbours. A node whose five are all unlabelled cannot be
    called by any threshold, in any round — so this is the ceiling every number below sits under,
    and it is a property of the stored map rather than of the policy.
    """
    print(f"\n{'=' * 118}\nS2b — THE CEILING: how many unlabelled nodes can REACH a labelled one at all\n{'=' * 118}")
    print("  ⛔ The walk is five deep. This is the ceiling on interpolation, and it is not a threshold.\n")
    print(
        f"  {'catalogue':<12}{'claim':<22}{'unlabelled':>11}{'ESM reach':>11}"
        f"{'BAC reach':>11}{'both':>8}{'neither':>9}"
    )
    for catalogue in CATALOGUES:
        for key, label, _kind, _rung in CLAIMS:
            problem, state = problems[catalogue][key], states[catalogue][key]
            unlabelled = [node for node in problem.gene_count if node not in state.label]
            esm = sum(nearest_labelled(problem, state, node, "ESM") is not None for node in unlabelled)
            bacformer = sum(nearest_labelled(problem, state, node, "BACFORMER") is not None for node in unlabelled)
            both = sum(
                nearest_labelled(problem, state, node, "ESM") is not None
                and nearest_labelled(problem, state, node, "BACFORMER") is not None
                for node in unlabelled
            )
            print(
                f"  {catalogue:<12}{label:<22}{len(unlabelled):>11,}{esm:>11,}{bacformer:>11,}{both:>8,}"
                f"{len(unlabelled) - max(esm, bacformer):>9,}"
            )


def main() -> None:
    """Run the stages David asked for, in his order, and print nothing that is not measured here."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        default="all",
        choices=(
            "all",
            "baseline",
            "intra",
            "inter",
            "extra",
            "interaction",
            "concordance",
            "loop",
            "uniref",
            "strata",
        ),
    )
    parser.add_argument("--cross-esm", required=False, help="ecoli→kp ESM cross-species map (TSV, with its .json)")
    parser.add_argument("--cross-bacformer", required=False, help="the same map built on Bacformer")
    parser.add_argument(
        "--mask-fraction",
        type=float,
        default=MASK_FRACTION,
        help=(
            "⭐ the share of annotated nodes held out. 0.20 is exhausted by round 1 — nearly every "
            "masked node is labelled there, so rounds 2+ have almost no held-out truth left to be "
            "scored against. A LARGER mask is the only way to measure a later round, and it costs "
            "donor density, so the two are reported side by side rather than one being chosen."
        ),
    )
    parser.add_argument(
        "--cross-depth",
        type=int,
        default=None,
        help="keep only this many cross-species donors per node — 5 matches the within-species walk",
    )
    parser.add_argument("--cross-esm-bar", type=float, default=None)
    parser.add_argument(
        "--json-out",
        default=None,
        help="write the one-pass record here, so a figure reads measured numbers rather than typed ones",
    )
    parser.add_argument(
        "--calls-out",
        default=None,
        help="write one row per held-out call with both cosines and whether it was right",
    )
    parser.add_argument("--cross-bacformer-bar", type=float, default=None)
    arguments = parser.parse_args()

    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        raise SystemExit("FATAL: set BACATLAS_DATABASE_URL")
    engine = create_engine(url, future=True)

    problems: dict[str, dict[str, Problem]] = {}
    gaps: dict[str, dict[str, int]] = {}
    with Session(engine) as session:
        for catalogue in CATALOGUES:
            problems[catalogue], gaps[catalogue] = load_problems(session, catalogue)
    states = {
        catalogue: {
            key: masked_state(problems[catalogue][key], fraction=arguments.mask_fraction)
            for key, *_ in CLAIMS
        }
        for catalogue in CATALOGUES
    }

    print(f"\n{'=' * 118}")
    print("HOW FAR CAN ANNOTATION PROPAGATE? — intrapolation, interpolation, extrapolation, looped")
    print(f"{'=' * 118}")
    print(
        f"  held out: {arguments.mask_fraction:.0%} of annotated nodes, seed {MASK_SEED} — "
        "the ONLY nodes precision is scored on"
    )
    print("  ⚠ coverage is CONSERVATIVE (a sparser donor graph than reality) and precision is OPTIMISTIC")
    print("    (annotated nodes standing in for unannotated ones, which sit in sparser neighbourhoods).")
    for catalogue in CATALOGUES:
        counts = {key: len(states[catalogue][key].masked) for key, *_ in CLAIMS}
        print(f"  {catalogue:<12}masked: " + " · ".join(f"{CLAIM_LABEL[k]} {v:,}" for k, v in counts.items()))

    if arguments.stage in ("all", "baseline"):
        report_baseline(problems, gaps)
        report_reach(problems, states)
        report_chance(problems, states)
    if arguments.stage in ("all", "intra"):
        report_intra(problems)
    if arguments.stage in ("all", "inter"):
        report_inter(problems, states)

    maps: dict[str, dict[str, list[tuple[str, float, int]]]] = {}
    if arguments.cross_esm and arguments.cross_bacformer:
        maps = {
            "ESM": read_cross_map(
                arguments.cross_esm, donor="ecoli", recipient="kp", depth=arguments.cross_depth
            ),
            "BACFORMER": read_cross_map(
                arguments.cross_bacformer, donor="ecoli", recipient="kp", depth=arguments.cross_depth
            ),
        }
        depths = [len(donors) for donors in maps["ESM"].values()]
        bacformer_depths = [len(donors) for donors in maps["BACFORMER"].values()]
        print(
            f"\n  cross-species shortlist depth — ESM mean {sum(depths) / len(depths):.1f}"
            f" · Bacformer mean {sum(bacformer_depths) / len(bacformer_depths):.1f}"
            f" donors per kp node{' (CAPPED at ' + str(arguments.cross_depth) + ')' if arguments.cross_depth else ''}"
        )
        print("  ⛔ The within-species walk stores FIVE. Extrapolation looking further is a property of")
        print("     how the maps were BUILT, not of the mechanism — run with --cross-depth 5 to compare.")
    elif arguments.stage in (
        "all", "extra", "interaction", "concordance", "loop", "uniref", "strata",
    ):
        raise SystemExit("FATAL: --cross-esm and --cross-bacformer are required for the extra/loop stages")

    cross_bars = {"ESM": arguments.cross_esm_bar, "BACFORMER": arguments.cross_bacformer_bar}
    if arguments.stage in ("all", "extra"):
        report_extra(problems, states, maps)
        report_cross_reach(problems, states, maps)
    if arguments.stage in ("all", "strata"):
        record = measure_one_pass(
            problems, states, maps, cross_bars, mask_fraction=arguments.mask_fraction
        )
        report_strata(record)
        calls = scoreable_calls(problems, states, maps, cross_bars)
        report_certainty_surface(calls)
        if arguments.calls_out:
            with open(arguments.calls_out, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(calls[0]))
                writer.writeheader()
                writer.writerows(calls)
            print(f"  ⭐ {len(calls):,} held-out calls written to {arguments.calls_out}")
        if arguments.json_out:
            # ⛔ The figure reads THIS, never a number typed into a plot script. `~/.claude/CLAUDE.md`:
            # a number is of record only if it names the artifact it was read from.
            Path(arguments.json_out).write_text(json.dumps(record, indent=2))
            print(f"\n  ⭐ one-pass record written to {arguments.json_out}")
    if arguments.stage in ("all", "interaction"):
        report_interaction(problems, states, maps, cross_bars)
    if arguments.stage in ("all", "concordance"):
        report_concordance(problems, states, maps, cross_bars)

    if arguments.stage in ("all", "loop"):
        print(f"\n{'=' * 118}\nS6 — THE LOOP: inter → extra → inter → … (intra applied continuously)\n{'=' * 118}")
        print(f"  cross-species bars held at ESM≥{cross_bars['ESM']} · Bacformer≥{cross_bars['BACFORMER']}")
        print("  ⛔ 'prec' is this row's NEW labels against held-out truth. Error compounding shows there")
        print("     and hides in a cumulative figure, so the cumulative one is not printed.")
        for key, label, _kind, _rung in CLAIMS:
            rows = run_loop(key, problems, states, maps, cross_bars)
            print_rounds(label, rows, problems, key)

    if arguments.stage in ("all", "uniref"):
        print(f"\n{'=' * 118}\nS7 — THE UNIREF50 BENCHMARK, IN BOTH POSITIONS\n{'=' * 118}")
        print("  ⭐ David's design. Arm A runs UniRef50 FIRST, Arm B runs it LAST. Arm B's marginal gain is")
        print("     how much UniRef50 knows that ESM and Bacformer did not already find.")
        print("  ⚠ Its error mode is its own: a family spanning two nodes may span two functions.")
        for key, label, _kind, _rung in CLAIMS:
            arm_a = run_loop(key, problems, states, maps, cross_bars, pre_uniref=True)
            arm_b = run_loop(key, problems, states, maps, cross_bars, post_uniref=True)
            print_rounds(f"{label} — ARM A (UniRef50 first)", arm_a, problems, key)
            print_rounds(f"{label} — ARM B (UniRef50 last)", arm_b, problems, key)

    print(f"\n{'=' * 118}")
    print("⛔ NOTHING HERE IS RECORDED. Gather the data, then STOP — these numbers go to David first.")
    print(f"{'=' * 118}\n")


if __name__ == "__main__":
    main()
