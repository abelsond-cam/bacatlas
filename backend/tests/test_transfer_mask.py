"""The held-out mask — ⛔ the one error that would make every number in this strand meaningless.

Design: `nuna/docs/model_evaluation/annotation_transfer_design.md` §2.3b step 1. `L_open` is the only
label any transfer step may read; `L_true` is the scorer's alone. These tests are written **before** the
module, because a leak here does not fail — it produces a better number.
"""

from __future__ import annotations

import pytest

from bacatlas_backend.transfer import mask as mask_module
from bacatlas_backend.transfer.mask import HeldOutMask, LabelLeak

# 20 carrying nodes (ids deliberately non-contiguous and unsorted) and 5 that carry nothing.
CARRYING = {
    node: frozenset({f"C{node}"})
    for node in (907, 3, 55, 2, 18, 400, 7, 91, 1, 33, 64, 12, 5, 77, 21, 8, 150, 42, 6, 29)
}
UNANNOTATED = (999, 998, 997, 996, 995)
ALL_NODES = tuple(CARRYING) + UNANNOTATED


def build(fraction: float = 0.2, seed: int = 0, claims=None) -> HeldOutMask:
    return HeldOutMask.build(claims if claims is not None else CARRYING, fraction=fraction, seed=seed)


# ====================================================================== what the mask IS
def test_the_mask_is_the_asked_fraction_of_the_CARRYING_nodes():
    m = build()
    assert len(m.held_out) == 4, "20 % of 20 carrying nodes"
    assert m.n_carrying == 20
    assert m.n_open == 16


def test_the_mask_is_drawn_ONLY_from_nodes_that_carry_the_claim():
    """⛔ Masking an unannotated node holds back nothing and inflates the held-out denominator with
    rows that have no truth to score against."""
    m = build()
    assert m.held_out <= set(CARRYING)
    assert not (m.held_out & set(UNANNOTATED))


def test_same_seed_same_mask_and_a_different_seed_moves_it():
    assert build(seed=0).held_out == build(seed=0).held_out
    assert build(seed=1).held_out != build(seed=0).held_out


def test_the_mask_does_not_depend_on_the_order_the_claims_ARRIVED_in():
    """⛔ The subtlest leak in the file. If the draw followed dict insertion order, two runs over the
    same catalogue with the same seed would hold out different nodes whenever the query returned rows
    in a different order — and every accuracy figure would be unreproducible while looking fine."""
    shuffled = {node: CARRYING[node] for node in sorted(CARRYING, reverse=True)}
    assert build(claims=shuffled).held_out == build().held_out


# =============================================================== L_open: the only readable label
def test_open_label_is_the_true_claim_for_a_carrying_unmasked_node():
    m = build()
    for node in set(CARRYING) - m.held_out:
        assert m.open_label(node) == CARRYING[node]


def test_a_masked_node_is_INDISTINGUISHABLE_from_an_unannotated_one():
    """⭐ The whole point of the mask: the model must not be able to tell "held back" from "never
    knew". Both are the empty set, and `open_label` is the only thing any step may ask."""
    m = build()
    for node in m.held_out:
        assert m.open_label(node) == frozenset()
    for node in UNANNOTATED:
        assert m.open_label(node) == frozenset()


def test_open_label_of_an_unknown_node_is_empty_not_an_error():
    """A node absent from the claim map is unannotated for this claim, which is the population the
    whole exercise exists to annotate — absence is the subject, never an error."""
    assert build().open_label(123_456) == frozenset()


# ======================================================= L_true: the scorer's, and nobody else's
def test_held_out_truth_returns_the_real_claim_for_a_held_out_node():
    m = build()
    for node in m.held_out:
        assert m.held_out_truth(node) == CARRYING[node]


def test_held_out_truth_REFUSES_a_node_the_model_could_already_see():
    """⛔ Scoring on a node whose label was never hidden is leakage, and it reads as a high accuracy."""
    m = build()
    open_node = next(iter(set(CARRYING) - m.held_out))
    with pytest.raises(LabelLeak, match="not held out"):
        m.held_out_truth(open_node)


def test_held_out_truth_REFUSES_an_unannotated_node():
    """There is no truth to return, so returning the empty set would silently score a prediction
    against "nothing" and count it wrong."""
    with pytest.raises(LabelLeak, match="carries no claim"):
        build().held_out_truth(UNANNOTATED[0])


def test_held_out_truth_REFUSES_A_NODE_THAT_HAS_BEEN_PROPAGATED_ONTO():
    """⛔⛔ The design's named rule. Once a predicted label has been written onto a node, its "truth"
    is a prediction — scoring against it compares the model with itself and cannot fail."""
    m = build()
    node = next(iter(m.held_out))
    m.note_propagated(node)
    with pytest.raises(LabelLeak, match="propagated"):
        m.held_out_truth(node)


def test_noting_a_propagation_does_not_change_any_other_node():
    m = build()
    a, b = sorted(m.held_out)[:2]
    m.note_propagated(a)
    assert m.held_out_truth(b) == CARRYING[b]


def test_a_propagated_node_still_reports_an_empty_open_label():
    """⚠ `note_propagated` records that a prediction was made; it does not make the prediction
    readable as a label. A step that could read its own output back as evidence would compound."""
    m = build()
    node = next(iter(m.held_out))
    m.note_propagated(node)
    assert m.open_label(node) == frozenset()


# ======================================================================== per-claim independence
def test_each_claim_gets_its_own_mask():
    """⚠ A node held out for COG may be open for molecular function. Each claim is its own model with
    its own carrying population, so one shared node list would hold out the wrong nodes for most of
    them — and would couple thirteen independent measurements through one draw."""
    other = {node: frozenset({f"G{node}"}) for node in list(CARRYING)[:10]}
    assert build().held_out != HeldOutMask.build(other, fraction=0.2, seed=0).held_out


# =============================================================================== the small-claim edge
def test_a_claim_too_small_to_hold_anything_back_is_refused_not_rounded_to_zero():
    """⛔ `int(3 * 0.2) == 0`. A mask of nothing yields an empty held-out set, and every downstream
    rate becomes 0/0 — which prints as a dash or a NaN rather than as "this claim was not scoreable"."""
    with pytest.raises(ValueError, match="too few"):
        HeldOutMask.build({1: frozenset({"a"}), 2: frozenset({"b"})}, fraction=0.2, seed=0)


def test_the_default_fraction_and_seed_are_the_design_s():
    assert mask_module.HELD_OUT_FRACTION == 0.2
    assert mask_module.SEED == 0
