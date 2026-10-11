"""The propagation loop's harness — the gates without which its numbers would mean nothing.

⛔ The first test here is the one that matters most. A loop scored against its own output measures
nothing and looks excellent doing it, so `State.truth_of` must REFUSE a node whose label the loop
supplied. Everything else in the harness can be wrong and be noticed; that one cannot.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from measure_annotation_propagation_loop import (  # noqa: E402
    CLAIMS,
    MASK_FRACTION,
    TARGET,
    Problem,
    State,
    _stratify,
    apply,
    census,
    interpolate,
    intrapolate,
    masked_state,
    precision,
    precision_by_gene,
    uniref_propagate,
)

DATABASE_URL = os.environ.get("BACATLAS_DATABASE_URL")
needs_database = pytest.mark.skipif(not DATABASE_URL, reason="BACATLAS_DATABASE_URL is not set")


def toy(**overrides) -> Problem:
    """Five nodes in a line: 1 and 2 annotated, 3-5 not. ESM and Bacformer walks are given."""
    defaults = dict(
        catalogue="toy",
        claim_key="cog",
        truth={1: frozenset({"A"}), 2: frozenset({"A"}), 3: frozenset({"B"})},
        gene_count={1: 100, 2: 10, 3: 1, 4: 1, 5: 1},
        genes_carrying={1: 40, 2: 5, 3: 1},
        walk={
            "ESM": {3: [(1, 1, 0.99)], 4: [(1, 3, 0.99)], 5: [(1, 4, 0.99)]},
            "BACFORMER": {3: [(1, 1, 0.80)], 4: [(1, 3, 0.80)], 5: [(1, 4, 0.80)]},
        },
        node_label={node: str(node) for node in (1, 2, 3, 4, 5)},
        family_nodes={},
        node_families={},
    )
    defaults.update(overrides)
    return Problem(**defaults)


# ============================================ ⛔ measurement rule 1: truth is truth, never our output
def test_truth_of_refuses_a_node_that_was_never_masked():
    """⛔ THE gate. Scoring a propagated label against itself is the one error that voids everything."""
    state = State(label={}, source={}, masked={7: frozenset({"A"})})
    assert state.truth_of(7) == frozenset({"A"})
    with pytest.raises(AssertionError, match="masked held-out set"):
        state.truth_of(8)


def test_a_propagated_label_is_never_scoreable_even_after_it_is_applied():
    """A node the loop labelled is in `label` but not in `masked`, so it cannot be scored."""
    state = State(label={}, source={}, masked={})
    apply(state, {9: frozenset({"A"})}, "inter", 1)
    assert state.label[9] == frozenset({"A"})
    assert state.source[9] == "inter:1"
    with pytest.raises(AssertionError):
        state.truth_of(9)


def test_precision_counts_an_unscoreable_call_in_called_and_not_in_scored():
    """⛔ Coverage and precision are measured on DIFFERENT populations, and the gap is the point."""
    state = State(label={}, source={}, masked={1: frozenset({"A"})})
    called, scored, rate = precision(state, {1: frozenset({"A"}), 2: frozenset({"Z"})})
    assert (called, scored, rate) == (2, 1, 1.0)


def test_precision_by_gene_exposes_the_amplification_node_precision_hides():
    """One wrong call on a 100-gene node and one right call on a singleton: 50 % by node, 1 % by gene."""
    problem = toy()
    state = State(label={}, source={}, masked={1: frozenset({"A"}), 3: frozenset({"B"})})
    assigned = {1: frozenset({"WRONG"}), 3: frozenset({"B"})}
    _called, _scored, by_node = precision(state, assigned)
    by_gene = precision_by_gene(problem, state, assigned)
    assert by_node == 0.5
    assert by_gene == pytest.approx(1 / 101)
    assert by_gene < by_node, "node-level precision cannot see that the wrong call covered 100 genes"


# ======================================================================= the census, and what intra is
def test_intrapolation_gains_genes_and_never_a_node():
    """⛔ By construction, and the report prints the zero rather than a flat line."""
    problem = toy()
    state = State(label=dict(problem.truth), source={}, masked={})
    nodes, genes = intrapolate(problem, state)
    assert nodes == 0
    assert genes == (100 - 40) + (10 - 5) + (1 - 1)


def test_census_separates_genes_carrying_from_genes_covered():
    """The two gene numbers differ by exactly what intrapolation does — that is why both are printed."""
    problem = toy()
    numbers = census(problem, State(label=dict(problem.truth), source={}, masked={}))
    assert numbers["nodes"] == 3
    assert numbers["genes_carrying"] == 46
    assert numbers["genes_covered"] == 111
    assert numbers["genes_covered_share"] == pytest.approx(111 / 113)


def test_a_gene_carrying_the_claim_counts_even_where_its_node_holds_no_label():
    """⛔ The defect an audit found on 2026-10-11, and the toy above could not have caught it.

    A node's label is the rank-0 modal vote over its members, so a node can hold genes stating a
    molecular function and still carry no rank-0 molecular-function row. Summing `genes_carrying`
    over LABELLED nodes only dropped every such gene from the baseline — 25,751 ecoli
    molecular-function genes, 5.26 points — and 1,157 of them sat inside nodes the pass later
    labelled, so the pass was credited with genes the annotation source had already annotated.
    """
    # Node 4 carries the claim on 1 of its 1 genes but has NO node-level label.
    problem = toy(
        truth={1: frozenset({"A"})},
        gene_count={1: 100, 4: 8},
        genes_carrying={1: 40, 4: 3},
    )
    numbers = census(problem, State(label=dict(problem.truth), source={}, masked={}))
    assert numbers["genes_carrying"] == 43, "node 4's three annotated genes are annotated at source"
    assert numbers["genes_covered"] == 103, "node 1 contributes all 100; node 4 only its three"
    assert numbers["nodes"] == 1, "node 4 still holds no node-level claim"


def test_the_pass_is_not_credited_with_genes_the_source_already_annotated():
    """The other half of the same defect: labelling node 4 must gain 5 genes, not 8."""
    problem = toy(
        truth={1: frozenset({"A"})},
        gene_count={1: 100, 4: 8},
        genes_carrying={1: 40, 4: 3},
    )
    state = State(label=dict(problem.truth), source={}, masked={})
    before = census(problem, state)
    apply(state, {4: frozenset({"A"})}, "inter", 1)
    after = census(problem, state)
    assert after["genes_covered"] - before["genes_covered"] == 5


def test_a_strata_rate_below_the_floor_is_none_in_BOTH_columns():
    """⛔ `precision_by_gene` escaped the guard and printed 0.0 % off one node beside an honest dash."""
    problem = toy(truth={1: frozenset({"A"}), 3: frozenset({"B"})})
    state = State(label={1: frozenset({"A"})}, source={}, masked={3: frozenset({"B"})})
    assigned = interpolate(problem, state, (None, None))
    rows = _stratify(problem, state, assigned)
    called = [row for row in rows if row["nodes"]]
    assert called, "node 3 should have been called"
    for row in called:
        assert row["scored"] < 30
        assert row["precision"] is None
        assert row["precision_by_gene"] is None, "the guard must cover the gene-weighted rate too"
        assert row["precision_unguarded"] == 0.0, "the unguarded value is kept, just not published"


def test_the_strata_gene_gain_sums_to_the_census_gene_gain():
    """⛔ The cross-figure identity. Two figures drawing the same gain must agree on its size.

    An audit found the census correction had reached `census` and not `_stratify`, so the
    certainty figure re-credited the pass with 1,157-6,672 genes per panel that the other figure
    had stopped claiming — both COG panels agreeing exactly, which is what identified it as the
    facet-splitting bug rather than a new one. This asserts the two can never drift again.
    """
    problem = toy(
        truth={1: frozenset({"A"})},
        gene_count={1: 100, 3: 1, 4: 8},
        genes_carrying={1: 40, 4: 3},
    )
    state = State(label=dict(problem.truth), source={}, masked={})
    before = census(problem, state)
    assigned = interpolate(problem, state, (None, None))
    rows = _stratify(problem, state, assigned)
    apply(state, assigned, "inter", 1)
    after = census(problem, state)
    assert sum(row["genes"] for row in rows) == after["genes_covered"] - before["genes_covered"]


def test_mask_removes_the_label_from_the_working_state_and_keeps_the_truth():
    """A masked node is genuinely unlabelled for the run, and its truth is reachable only via `masked`."""
    problem = toy(truth={node: frozenset({"A"}) for node in range(1, 11)})
    state = masked_state(problem, fraction=0.5, seed=1)
    assert len(state.masked) == 5
    assert len(state.label) == 5
    assert not (set(state.masked) & set(state.label))
    assert all(state.source[node] == "truth" for node in state.label)


def test_the_mask_is_the_same_set_for_the_same_seed_and_different_for_another():
    """Every mechanism must be measured on ONE mask, or the mechanisms are not comparable."""
    problem = toy(truth={node: frozenset({"A"}) for node in range(1, 51)})
    assert set(masked_state(problem, seed=0).masked) == set(masked_state(problem, seed=0).masked)
    assert set(masked_state(problem, seed=0).masked) != set(masked_state(problem, seed=7).masked)


# ============================================================================= reaching a donor at all
def test_interpolation_cannot_reach_past_the_stored_walk():
    """⛔ The ceiling. Node 5's only neighbour is 4, which is unlabelled, so no bar can call it."""
    problem = toy()
    state = State(label={1: frozenset({"A"})}, source={1: "truth"}, masked={})
    assigned = interpolate(problem, state, (None, None))
    assert 3 in assigned, "node 3's neighbour 1 is labelled in both maps"
    assert 4 not in assigned, "node 4's only neighbour is 3, which is not labelled YET"
    assert 5 not in assigned


def test_a_second_round_reaches_what_the_first_could_not():
    """That is the loop: round 2's gain is the nodes whose five neighbours BECAME labelled."""
    problem = toy()
    state = State(label={1: frozenset({"A"})}, source={1: "truth"}, masked={})
    apply(state, interpolate(problem, state, (None, None)), "inter", 1)
    second = interpolate(problem, state, (None, None))
    assert set(second) == {4}


def test_a_round_that_adds_nothing_is_empty_so_the_loop_can_terminate():
    """The loop's stop condition is a measurement, not a round counter."""
    problem = toy()
    state = State(label=dict(problem.truth), source={}, masked={})
    apply(state, interpolate(problem, state, (None, None)), "inter", 1)
    apply(state, interpolate(problem, state, (None, None)), "inter", 2)
    assert interpolate(problem, state, (None, None)) == {}


def test_concurrence_refuses_where_the_two_maps_carry_different_claims():
    """⚠ `concur` is a consensus test: two geometries agreeing on one answer, not two predictions."""
    problem = toy(
        walk={
            "ESM": {3: [(1, 1, 0.99)]},
            "BACFORMER": {3: [(1, 2, 0.80)]},
        },
        truth={1: frozenset({"A"}), 2: frozenset({"DIFFERENT"})},
    )
    state = State(label=dict(problem.truth), source={}, masked={})
    assert interpolate(problem, state, (None, None), rule="concur") == {}
    assert interpolate(problem, state, (None, None), rule="esm") == {3: frozenset({"A"})}


def test_a_bar_refuses_a_donor_below_it():
    """The bar is read off the DONOR's cosine, which may come from rank 4 rather than rank 1."""
    problem = toy(walk={"ESM": {3: [(1, 1, 0.95)]}, "BACFORMER": {3: [(1, 1, 0.80)]}})
    state = State(label={1: frozenset({"A"})}, source={}, masked={})
    assert interpolate(problem, state, (0.96, None)) == {}
    assert interpolate(problem, state, (0.94, None)) == {3: frozenset({"A"})}


# ======================================================================== error compounding, explicitly
def test_a_round_one_error_is_inherited_by_round_two_and_the_round_rate_shows_it():
    """⭐ The synthetic case the plan demands: round 2's precision must be able to FALL.

    Node 3's donor is node 1, which carries `A`; the held-out truth for 3 is `B`, so round 1 is
    wrong. Node 4's donor is node 3, so round 2 inherits `A` against a held-out truth of `B` — and
    the per-round rate drops while a cumulative rate would hide it behind round 1's successes.
    """
    problem = toy(truth={1: frozenset({"A"}), 3: frozenset({"B"}), 4: frozenset({"B"})})
    state = State(
        label={1: frozenset({"A"})},
        source={1: "truth"},
        masked={3: frozenset({"B"}), 4: frozenset({"B"})},
    )
    first = interpolate(problem, state, (None, None))
    _called, scored_one, rate_one = precision(state, first)
    apply(state, first, "inter", 1)
    second = interpolate(problem, state, (None, None))
    _called, scored_two, rate_two = precision(state, second)
    assert (scored_one, rate_one) == (1, 0.0)
    assert set(second) == {4}
    assert (scored_two, rate_two) == (1, 0.0), "round 2 inherited round 1's wrong label"
    apply(state, second, "inter", 2)
    assert state.label[4] == frozenset({"A"}), "the error crossed a round, which is the risk being measured"


def test_uniref_assigns_nothing_where_a_family_s_labelled_nodes_disagree():
    """⚠ A family spanning two nodes may span two functions — the paralogue split, re-merged."""
    problem = toy(
        family_nodes={"U1": {1, 2, 3}},
        node_families={1: {"U1"}, 2: {"U1"}, 3: {"U1"}},
        truth={1: frozenset({"A"}), 2: frozenset({"B"})},
    )
    state = State(label={1: frozenset({"A"}), 2: frozenset({"B"})}, source={}, masked={})
    assert uniref_propagate(problem, state) == {}

    agreeing = toy(
        family_nodes={"U1": {1, 2, 3}},
        node_families={1: {"U1"}, 2: {"U1"}, 3: {"U1"}},
        truth={1: frozenset({"A"}), 2: frozenset({"A", "C"})},
    )
    state = State(label={1: frozenset({"A"}), 2: frozenset({"A", "C"})}, source={}, masked={})
    assert uniref_propagate(agreeing, state) == {3: frozenset({"A"})}


# ============================================================ the census identity, against SQL not itself
@needs_database
@pytest.mark.parametrize("catalogue", ["ecoli-nuna4", "kp-nuna4"])
def test_node_coverage_equals_nodes_with_at_least_one_annotated_gene(catalogue):
    """⛔ The oracle is SQL, independent of the harness — not the harness compared with itself.

    COG only: its `locus.cog_annotated_member_count` column was verified gene-by-gene (0 mismatches
    out of 17,531 and 15,670), so it is an independent statement of the same fact. The GO columns
    are NOT, which is why the harness folds the gene rows and reports the gap instead.
    """
    from measure_annotation_propagation_loop import load_problems  # noqa: PLC0415

    engine = create_engine(DATABASE_URL, future=True)
    with Session(engine) as session:
        problems, _gap = load_problems(session, catalogue)
        expected = session.execute(
            text(
                """
                select count(*) from locus
                 where pangenome_id = (select pangenome_id from pangenome where catalogue_key = :key)
                   and cog_annotated_member_count > 0
                """
            ),
            {"key": catalogue},
        ).scalar_one()
    assert len(problems["cog"].truth) == expected


@needs_database
@pytest.mark.parametrize("catalogue", ["ecoli-nuna4", "kp-nuna4"])
def test_the_masked_nodes_are_absent_from_the_donor_pool(catalogue):
    """⛔ A masked node that could still donate would make the held-out set worthless."""
    from measure_annotation_propagation_loop import load_problems  # noqa: PLC0415

    engine = create_engine(DATABASE_URL, future=True)
    with Session(engine) as session:
        problems, _gap = load_problems(session, catalogue)
    for key, _label, _kind, _rung in CLAIMS:
        problem = problems[key]
        state = masked_state(problem)
        assert len(state.masked) == round(len(problem.truth) * MASK_FRACTION)
        assert not (set(state.masked) & set(state.label))
        assigned = interpolate(problem, state, (None, None))
        for node in assigned:
            assert node not in state.label, "a node that already had a label was re-called"


def test_the_target_is_the_bar_the_plan_states():
    """A constant, pinned, because the whole exercise is reported against it."""
    assert TARGET == 0.90
