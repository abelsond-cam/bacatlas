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
    #: The donor that supplied each propagated label, for the worked examples.
    donor: dict[int, tuple[str, int, float | None]] = field(default_factory=dict)

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
    `genes_carrying` counts only the member genes that hold the claim themselves, `genes_covered`
    counts every gene of every labelled node — which is what the node's label means once intra has
    run. ⛔ Node coverage is IDENTICAL in both: intra cannot annotate a node that had nothing.
    """
    nodes = len(state.label)
    carrying = sum(problem.genes_carrying.get(node, 0) for node in state.label)
    covered = sum(problem.gene_count[node] for node in state.label)
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
            state.donor[node] = (rule, chosen[2], chosen[1])
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
            recipient_state.donor[node] = (f"extra:{rule}", chosen[2], chosen[1])
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
            state.donor[node] = ("uniref50", witness, None)
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
        choices=("all", "baseline", "intra", "inter", "extra", "loop", "uniref"),
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
    elif arguments.stage in ("all", "extra", "loop", "uniref"):
        raise SystemExit("FATAL: --cross-esm and --cross-bacformer are required for the extra/loop stages")

    cross_bars = {"ESM": arguments.cross_esm_bar, "BACFORMER": arguments.cross_bacformer_bar}
    if arguments.stage in ("all", "extra"):
        report_extra(problems, states, maps)
        report_cross_reach(problems, states, maps)

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
