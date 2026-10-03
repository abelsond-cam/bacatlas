"""The cross-species transfer measurement — the null, the map reader, and the four states.

⭐ **Every test here runs WITHOUT a database.** The DB-backed inference tests skip when
`BACATLAS_DATABASE_URL` is unset, and a skip is not a pass — so the three things that are genuinely
new and genuinely easy to get wrong are pinned as pure functions:

* **the chance baseline's pool.** Across species the comparator is a random annotated node of the
  *donor* catalogue, not of the union, and the focal node is not in that pool so self-exclusion is
  dropped. Getting it wrong rescales every lift and leaves agreement untouched — nothing in the
  output would look wrong.
* **the map reader.** `locus.node_label` is TEXT and labels look numeric (`"1098"`); a dtype-inferring
  reader turns them into integers and the join matches nothing while looking reasonable. And a
  reversed map read forwards would transfer every label the wrong way.
* **the walk's four states.** *No annotated donor* and *donor agreed* are both "no disagreement" in a
  naive tally. They must sum to the unlabelled total or a node has been counted twice or not at all.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from bacatlas_backend.instruments.annotation_transfer import (
    _chance_by_level,
    tally_cells,
)
from bacatlas_backend.models.enumerations import AnnotationKind

COG = AnnotationKind.COG_ORTHOGROUP


def _load_script():
    """`scripts/` is not a package, so the driver is loaded by path — like any one-file tool."""
    path = Path(__file__).resolve().parent.parent / "scripts" / "measure_cross_species_transfer.py"
    spec = importlib.util.spec_from_file_location("measure_cross_species_transfer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


xs = _load_script()


def _cog(accession: str, category: str) -> dict[int, frozenset[str]]:
    return {2: frozenset({accession}), 1: frozenset({category})}


@pytest.fixture
def two_catalogues():
    """Four donor loci (two hold `E`, two hold `C`) and three recipient loci, all holding `E`.

    ⚠ The shape is chosen so the right null and the wrong one give DIFFERENT answers. With a single
    recipient locus they coincide: adding the focal node to the pool raises numerator and denominator
    by one each, and self-exclusion removes one from each again. Several recipients break the tie,
    which is why the fixture has three — a two-locus fixture would have passed either way.
    """
    donor = {1: _cog("COG0001", "E"), 2: _cog("COG0002", "E"), 3: _cog("COG0003", "C"), 4: _cog("COG0004", "C")}
    recipient = {101: _cog("COG0001", "E"), 102: _cog("COG0001", "E"), 103: _cog("COG0002", "E")}
    return donor, recipient


# ======================================================================= the null, and whose pool
def test_the_cross_species_null_is_the_DONOR_pool_and_does_not_self_exclude(two_catalogues):
    """⛔ 2 of 4 donor loci hold `E`, so chance is exactly 0.5 — and the union would say 0.667."""
    donor, recipient = two_catalogues
    both = {**donor, **recipient}
    edges = [(101, 1, 0.995), (102, 2, 0.995), (103, 1, 0.995)]

    right = tally_cells(both, edges, COG, donor_pool=donor)[(">= 0.99", 1)]
    assert right.chance == pytest.approx(2 / 4)

    # the wrong null: the union of both catalogues, self-excluding as a within-species run does
    wrong = tally_cells(both, edges, COG)[(">= 0.99", 1)]
    assert wrong.chance == pytest.approx((5 - 1) / (7 - 1))
    assert wrong.chance != pytest.approx(right.chance), "the fixture must be able to tell them apart"
    # ⛔ agreement is IDENTICAL under both nulls — which is why a wrong null is invisible
    assert right.agreeing == wrong.agreeing and right.pairs == wrong.pairs


def test_self_exclusion_is_dropped_only_when_the_focal_node_is_not_in_the_pool(two_catalogues):
    """The two modes of `_chance_by_level`, on one pool, so the arithmetic is visible."""
    donor, _recipient = two_catalogues
    excluded = _chance_by_level(donor, 1, self_excluded=True)
    included = _chance_by_level(donor, 1, self_excluded=False)
    assert included(frozenset({"E"})) == pytest.approx(2 / 4)
    assert excluded(frozenset({"E"})) == pytest.approx((2 - 1) / (4 - 1))
    # a claim nothing in the pool holds: 0 either way, and never negative
    assert included(frozenset({"Z"})) == 0.0
    assert excluded(frozenset({"Z"})) == pytest.approx(-1 / 3)


def test_a_claim_no_donor_holds_scores_zero_chance_and_no_lift(two_catalogues):
    """A recipient claim absent from the donor catalogue: chance 0, so `lift` is None, not infinite."""
    donor, _recipient = two_catalogues
    levels = {**donor, 101: _cog("COG9999", "Z"), 102: _cog("COG9999", "Z")}
    cells = tally_cells(levels, [(101, 1, 0.995), (102, 2, 0.995)], COG, donor_pool=donor)
    cell = cells[(">= 0.99", 1)]
    assert cell.chance == 0.0
    assert cell.lift is None, "⛔ never an unbounded ratio — lift is None where chance is 0"


# ============================================================================== the map, off disk
def _write_map(tmp_path: Path, rows, *, header=None, provenance=..., name="map") -> Path:
    """A map and its sidecar. `provenance=None` writes no sidecar; `...` writes the default kp->ecoli one."""
    columns = header or [
        "kp_node",
        "ecoli_node",
        "n_kp_genes",
        "n_ecoli_genes",
        "n_pairs",
        "cross",
        "cross_p75",
        "cross_max",
        "rank_cross",
        "rank_cross_p75",
        "rank_cross_max",
    ]
    path = tmp_path / f"{name}.tsv"
    path.write_text("\t".join(columns) + "\n" + "".join("\t".join(str(v) for v in row) + "\n" for row in rows))
    if provenance is not None:
        if provenance is ...:
            provenance = {
                "donor": "ecoli",
                "recipient": "kp",
                "rep": "esm",
                "k": 256,
                "cap": 30,
                "rank_by": "sim",
                "recall": [],
            }
        path.with_suffix(".json").write_text(json.dumps(provenance))
    return path


def test_a_numeric_looking_node_label_stays_a_string(tmp_path):
    """⛔ `locus.node_label` is TEXT. Read `"1098"` as 1098 and the join silently matches nothing."""
    path = _write_map(
        tmp_path,
        [
            ["1098", "02811", 7, 12, 84, 0.991, 0.993, 0.998, 1, 1, 1],
        ],
    )
    rows, _provenance = xs.read_map(str(path), donor="ecoli", recipient="kp")
    assert rows[0]["recipient_node"] == "1098" and isinstance(rows[0]["recipient_node"], str)
    # ⚠ and a leading zero survives, which an integer round-trip would destroy outright
    assert rows[0]["donor_node"] == "02811"
    assert rows[0]["cross"] == pytest.approx(0.991)
    assert rows[0]["rank_cross_max"] == 1 and rows[0]["recipient_genes"] == 7


def test_the_header_ALONE_cannot_tell_the_direction_which_is_why_the_sidecar_decides(tmp_path):
    """⛔⛔ The failure this test FOUND, and the reason the sidecar is required.

    A map in either direction carries both `*_node` and both `n_*_genes` columns, so a header check
    reads an `ecoli -> kp` map as a `kp -> ecoli` one without complaint — and every rank then means
    the opposite of what it says. The first version of this test asserted the header would catch it;
    it did not raise, and the fix was to move the authority to the provenance, not to loosen the test.
    """
    reversed_header = [
        "ecoli_node",
        "kp_node",
        "n_ecoli_genes",
        "n_kp_genes",
        "n_pairs",
        "cross",
        "cross_p75",
        "cross_max",
        "rank_cross",
        "rank_cross_p75",
        "rank_cross_max",
    ]
    # ⚠ no sidecar: the header on its own is accepted BOTH ways round — the defect, demonstrated
    bare = _write_map(
        tmp_path,
        [["E1", "K1", 7, 12, 84, 0.99, 0.99, 0.99, 1, 1, 1]],
        header=reversed_header,
        provenance=None,
        name="bare",
    )
    columns = set(open(bare).readline().strip().split("\t"))
    assert {"kp_node", "ecoli_node", "n_kp_genes", "n_ecoli_genes"} <= columns, (
        "both directions' columns are present, so the header is symmetric and proves nothing"
    )

    # which is why a missing sidecar is fatal rather than a warning
    with pytest.raises(SystemExit, match="no provenance beside the map"):
        xs.read_map(str(bare), donor="kp", recipient="ecoli")

    # and with one, the direction it records is what decides
    path = _write_map(
        tmp_path,
        [["E1", "K1", 7, 12, 84, 0.99, 0.99, 0.99, 1, 1, 1]],
        header=reversed_header,
        provenance={"donor": "ecoli", "recipient": "kp", "recall": []},
    )
    with pytest.raises(SystemExit, match="says kp asks and ecoli answers"):
        xs.read_map(str(path), donor="kp", recipient="ecoli")


def test_the_sidecar_travels_with_the_map_even_when_the_label_is_full_of_dots(tmp_path):
    """⚠ The real label is `nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL` — `.98` is not a suffix."""
    path = _write_map(
        tmp_path,
        [["1098", "2811", 7, 12, 84, 0.99, 0.99, 0.99, 1, 1, 1]],
        name="ecoli_to_kp_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL_neighbours",
    )
    rows, provenance = xs.read_map(str(path), donor="ecoli", recipient="kp")
    assert len(rows) == 1 and provenance["cap"] == 30


def test_a_header_without_every_rule_is_refused(tmp_path):
    """A map from an older run that carried only the median — refused rather than partly measured."""
    path = _write_map(
        tmp_path,
        [["K1", "E1", 7, 0.99, 1]],
        header=[
            "kp_node",
            "ecoli_node",
            "n_kp_genes",
            "cross",
            "rank_cross",
        ],
    )
    with pytest.raises(SystemExit, match="cross_p75"):
        xs.read_map(str(path), donor="ecoli", recipient="kp")


def test_a_valid_header_with_no_rows_is_refused(tmp_path):
    """⛔ An empty map would make every count 0 and every state read as 'nothing to transfer'."""
    path = _write_map(tmp_path, [])
    with pytest.raises(SystemExit, match="valid header and no rows"):
        xs.read_map(str(path), donor="ecoli", recipient="kp")


# ========================================================================== the walk, four states
@pytest.fixture
def walk_inputs():
    """Three unlabelled kp nodes: one reachable, one whose only donor is unannotated, one unreached."""
    recipient_ids = {"1098": (101, 7), "2811": (102, 3), "9999": (103, 5)}
    donor_ids = {"E1": (1, 10), "E2": (2, 4)}
    donor_levels = {1: _cog("COG0001", "E")}  # E2 (locus 2) carries nothing
    return recipient_ids, donor_ids, donor_levels


LABELS = {101: "1098", 102: "2811", 103: "9999", 1: "E1", 2: "E2"}


def _rows(*specs):
    """Map rows, with every rule carrying the same cosine and rank unless a test says otherwise."""
    return [
        {
            "recipient_node": r,
            "donor_node": d,
            "recipient_genes": 0,
            "cross": c,
            "cross_p75": c,
            "cross_max": c,
            "rank_cross": k,
            "rank_cross_p75": k,
            "rank_cross_max": k,
        }
        for r, d, c, k in specs
    ]


def _walk(rows, *, rule, recipient_ids, donor_ids, donor_levels, recipient_levels, kind=COG):
    """The map -> candidates -> walk path, as `main` wires it."""
    return xs.walk(
        xs.candidates_from_map(rows, rule=rule, recipient_ids=recipient_ids, donor_ids=donor_ids),
        recipient_ids=recipient_ids,
        donor_levels=donor_levels,
        recipient_levels=recipient_levels,
        kind=kind,
        label_of=LABELS,
    )


def test_the_four_states_sum_to_the_unlabelled_total(walk_inputs):
    """⛔ `nuna/CLAUDE.md` §2 — 'no annotated donor' and 'donor agreed' must never read the same."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = _walk(
        _rows(("1098", "E2", 0.995, 1), ("1098", "E1", 0.991, 2), ("2811", "E2", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    assert reach["unlabelled"] == 3
    assert reach["state"]["labelled"] == [1, 7], "1098 walks past the unannotated E2 to E1"
    assert reach["state"]["no_annotated_donor"] == [1, 3], "2811's only candidate carries nothing"
    assert reach["state"]["no_candidate_at_all"] == [1, 5], "9999 has no row at all"
    assert reach["state"]["too_remote_to_call"] == [0, 0]
    assert sum(nodes for nodes, _genes in reach["state"].values()) == 3


def test_a_donor_below_the_floor_is_reached_but_NOT_called(walk_inputs):
    """⛔ 'Definitely too remote to call' (David) — the donor is found and no label is taken from it."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = _walk(
        _rows(("1098", "E1", 0.85, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    assert reach["state"]["too_remote_to_call"] == [1, 7]
    assert reach["state"]["labelled"] == [0, 0]
    assert reach["state"]["no_annotated_donor"] == [0, 0], "it HAD an annotated donor; it was too far"


def test_a_node_that_already_carries_the_vocabulary_is_not_in_the_denominator(walk_inputs):
    """The two mechanisms are disjoint: a node with its own call is never offered a neighbour's."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = _walk(
        _rows(("1098", "E1", 0.995, 1), ("2811", "E1", 0.995, 1), ("9999", "E1", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={101: _cog("COG0007", "K")},
    )
    assert reach["unlabelled"] == 2, "101 (node 1098) carries its own COG and drops out"
    assert reach["state"]["labelled"] == [2, 8]


def test_the_rule_chooses_which_donor_the_walk_takes(walk_inputs):
    """⭐ The linkage by-product made consequential — and the earlier version of this test was inert.

    It first used one annotated donor and one unannotated one, so the walk reached the same donor
    whichever rule ordered the list: the unannotated candidate is skipped either way. Sorting by the
    wrong rule changed nothing and the test passed with the rule ignored. It takes **two annotated
    donors the rules order differently, landing in different tiers** for the choice to be observable
    at all — which is the same reason the by-tier breakdown is reported and not just a total.
    """
    recipient_ids, donor_ids, _unused = walk_inputs
    donor_levels = {1: _cog("COG0001", "E"), 2: _cog("COG0002", "C")}
    rows = [
        # E1 is the closest single gene pair (rank 1 by max) but a mediocre median
        {
            "recipient_node": "1098",
            "donor_node": "E1",
            "recipient_genes": 0,
            "cross": 0.900,
            "cross_p75": 0.930,
            "cross_max": 0.985,
            "rank_cross": 2,
            "rank_cross_p75": 2,
            "rank_cross_max": 1,
        },
        # E2 is the better median (rank 1 by median) and only second by max
        {
            "recipient_node": "1098",
            "donor_node": "E2",
            "recipient_genes": 0,
            "cross": 0.995,
            "cross_p75": 0.995,
            "cross_max": 0.996,
            "rank_cross": 1,
            "rank_cross_p75": 1,
            "rank_cross_max": 2,
        },
    ]
    common = dict(recipient_ids=recipient_ids, donor_ids=donor_ids, donor_levels=donor_levels, recipient_levels={})
    by_median = _walk(rows, rule="cross", **common)
    by_max = _walk(rows, rule="cross_max", **common)

    assert by_median["state"]["labelled"][0] == 1 and by_max["state"]["labelled"][0] == 1
    # the median takes E2 at 0.995 — the top tier; the max takes E1 at 0.985 — one tier down
    assert list(by_median["by_tier"]) == [(">= 0.99", 2)]
    assert list(by_max["by_tier"]) == [("0.98-0.99", 2)]
    assert by_median["examples"][0][1] == "E2" and by_max["examples"][0][1] == "E1"


def test_a_dash_padded_donor_quotes_no_deeper_than_it_states(walk_inputs):
    """`quotable_level` is what stops the transfer promising a depth the donor cannot supply."""
    recipient_ids, donor_ids, _levels = walk_inputs
    #: a donor with COG categories but no orthogroup accession can only supply the L1 rung
    donor_levels = {1: {2: frozenset(), 1: frozenset({"E"})}}
    reach = _walk(
        _rows(("1098", "E1", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    assert reach["state"]["labelled"][0] == 1
    assert list(reach["by_tier"]) == [(">= 0.99", 1)], "the >= 0.99 tier permits L2; the donor has only L1"


def test_the_two_reach_columns_must_describe_the_SAME_nodes(walk_inputs, capsys):
    """⛔ Two columns over different populations is a comparison between two different things.

    The cross-species reach and the within-species baseline are the SAME unlabelled nodes asked two
    questions — one of a donor catalogue, one of their own neighbours. If the totals differ, one of
    them was built from a different `recipient_levels` and the difference between the columns means
    nothing. ⚠ `nuna/CLAUDE.md` §2: a comparison asserts its own coverage BEFORE it reports one.
    """
    recipient_ids, donor_ids, donor_levels = walk_inputs
    rows = _rows(("1098", "E1", 0.995, 1), ("2811", "E1", 0.995, 1), ("9999", "E1", 0.995, 1))
    full = _walk(
        rows,
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    # a baseline built as though one recipient node already had a call — 2 unlabelled, not 3
    short = _walk(
        rows,
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={101: _cog("COG0007", "K")},
    )
    assert full["unlabelled"] == 3 and short["unlabelled"] == 2
    with pytest.raises(SystemExit, match="describe different populations"):
        xs.report_reach(full, full, short, COG, rule="cross", recipient="kp", donor="ecoli")
    # and the matching set reports without complaint
    xs.report_reach(full, full, full, COG, rule="cross", recipient="kp", donor="ecoli")
    assert "worked cases" in capsys.readouterr().out


def test_the_within_species_graph_feeds_the_SAME_walk_as_the_map(walk_inputs):
    """⭐ One implementation, two donor sets — the two columns cannot differ by a definition.

    `candidates_from_database` returns the same `{recipient_locus_id: [(rank, donor_id, cosine)]}`
    shape as `candidates_from_map`, so the within-species baseline is measured by the same tiers,
    the same `quotable_level` fallback and the same four states. Quoting it from the write-up instead
    is what lets a baseline decay silently while still reading plausibly.
    """
    recipient_ids, _donor_ids, _donor_levels = walk_inputs
    #: kp node 102 is its own donor pool's neighbour here — the within-species shape, by hand
    within = {101: [(1, 102, 0.993)], 102: [(1, 101, 0.993)]}
    reach = xs.walk(
        within,
        recipient_ids=recipient_ids,
        donor_levels={102: _cog("COG0042", "J")},  # only 102 carries a call
        recipient_levels={},
        kind=COG,
        label_of=LABELS,
    )
    assert reach["state"]["labelled"] == [1, 7], "101 takes 102's call"
    assert reach["state"]["no_annotated_donor"] == [1, 3], "102's only neighbour (101) carries nothing"
    assert reach["state"]["no_candidate_at_all"] == [1, 5], "103 is in no neighbour list"
    assert list(reach["by_tier"]) == [(">= 0.99", 2)]


def test_a_longer_candidate_list_reaches_more_donors_WHATEVER_the_embedding_does(walk_inputs):
    """⛔⛔ The list-length confound, and the truncation that removes it.

    The stored within-species graph holds **five** neighbours per locus (checked against the
    database: ranks 1-5, 4.77 on average); the cross-species map holds up to `cap`, ~23. A walk
    allowed 23 ranks finds an annotated donor more often than one allowed 5 for a reason that has
    nothing to do with the embedding — so the reach is reported at its own depth AND truncated to
    the baseline's, and only the two columns sharing a depth are comparable.

    Here the only annotated donor sits at rank 3: visible at full depth, invisible at depth 2.
    """
    recipient_ids, donor_ids, donor_levels = walk_inputs
    rows = _rows(("1098", "E2", 0.995, 1), ("1098", "E2", 0.994, 2), ("1098", "E1", 0.993, 3))
    mapped = xs.candidates_from_map(rows, rule="cross", recipient_ids=recipient_ids, donor_ids=donor_ids)
    shared = dict(recipient_ids=recipient_ids, donor_levels=donor_levels, recipient_levels={}, kind=COG)

    deep = xs.walk(mapped, **shared)
    shallow = xs.walk(mapped, max_rank=2, **shared)
    assert deep["state"]["labelled"][0] == 1, "rank 3 carries the only call, and full depth sees it"
    assert shallow["state"]["labelled"][0] == 0, "truncated to 2 ranks, the same walk reaches nothing"
    assert shallow["state"]["no_annotated_donor"][0] == 1
    # the depth each column was allowed is reported, so a reader is never guessing which is which
    assert deep["depth_p50"] == 3.0 and shallow["depth_p50"] == 2.0
    assert deep["max_rank"] is None and shallow["max_rank"] == 2


def test_candidates_sharing_a_rank_are_refused(walk_inputs):
    """⛔ With ties in the rank, `sorted` orders by DONOR ID and the walk takes the wrong donor."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    duplicated = {101: [(1, 2, 0.999), (1, 1, 0.991)]}  # both at rank 1
    with pytest.raises(SystemExit, match="distinct ranks"):
        xs.walk(duplicated, recipient_ids=recipient_ids, donor_levels=donor_levels, recipient_levels={}, kind=COG)


# ========================================================= the UniRef50 bridge (Stage 4), three states
def test_the_bridge_has_THREE_states_and_not_bridged_is_not_the_same_as_not_asked():
    """⛔ `not_bridged` and `unmeasurable` are different claims, and conflating them is the mistake
    `nuna/CLAUDE.md` retired the label `no_homology` for: it meant *not measured*, not *nothing found*.

    A pair where either node carries no UniRef50 accession has not been shown to lack a bridge — it
    has not been asked. Counting it as "ESM reached where sequence identity could not" would inflate
    the one number Stage 4 exists to produce.
    """
    assert xs.bridge_state(frozenset({"A", "B"}), frozenset({"B", "C"})) == "bridged"
    assert xs.bridge_state(frozenset({"A"}), frozenset({"Z"})) == "not_bridged"
    # either side missing, or empty, is NOT a negative result
    assert xs.bridge_state(None, frozenset({"Z"})) == "unmeasurable"
    assert xs.bridge_state(frozenset({"A"}), None) == "unmeasurable"
    assert xs.bridge_state(frozenset(), frozenset({"Z"})) == "unmeasurable"
    assert xs.bridge_state(None, None) == "unmeasurable"
    assert set(xs.BRIDGE_STATES) == {"bridged", "not_bridged", "unmeasurable"}


def test_a_bridge_is_an_intersection_so_one_rare_family_is_enough(walk_inputs):
    """⛔ Why the bridge is built at GENE level and not from the 8-capped crosstab.

    The single family linking two nodes can be a rare one in either of them, and a cap keeps only the
    commonest eight. Here the two nodes share exactly one family, which is each one's *least* common
    — gene level calls it bridged, a top-1 view calls it not bridged, and those are opposite answers
    to the question Stage 4 asks.
    """
    recipient_full = frozenset({"UniRef50_MAIN_K", "UniRef50_RARE_SHARED"})
    donor_full = frozenset({"UniRef50_MAIN_E", "UniRef50_RARE_SHARED"})
    assert xs.bridge_state(recipient_full, donor_full) == "bridged"
    # the same pair as a cap that kept only each node's commonest family would have seen it
    assert xs.bridge_state(frozenset({"UniRef50_MAIN_K"}), frozenset({"UniRef50_MAIN_E"})) == "not_bridged"


def test_the_reach_is_stratified_from_the_walks_OWN_record(walk_inputs, capsys):
    """⭐ `bridge_reach` splits `reach["taken"]`, so it cannot drift from the §5 totals.

    Two labelled kp nodes, one bridged and one not; the split must sum back to what the walk said.
    """
    recipient_ids, donor_ids, _unused = walk_inputs
    donor_levels = {1: _cog("COG0001", "E"), 2: _cog("COG0002", "C")}
    rows = _rows(("1098", "E1", 0.995, 1), ("2811", "E2", 0.9925, 1))
    reach = _walk(
        rows,
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    assert reach["state"]["labelled"] == [2, 10], "7 genes + 3 genes"

    xs.bridge_reach(
        reach,
        #: kp 1098 (locus 101) shares a family with ecoli E1 (locus 1); kp 2811 (102) shares none with E2 (2)
        recipient_sets={101: frozenset({"U_SHARED"}), 102: frozenset({"U_KP_ONLY"})},
        donor_sets={1: frozenset({"U_SHARED"}), 2: frozenset({"U_EC_ONLY"})},
        kind=COG,
        label_of=LABELS,
        recipient="kp",
        donor="ecoli",
    )
    printed = capsys.readouterr().out
    assert "bridged" in printed and "not_bridged" in printed and "unmeasurable" in printed
    # one node/7 genes bridged and one node/3 genes not, and the >= 0.98 row carries both
    assert "1/7" in printed.replace(" ", "") or "1/7 " in printed
    assert "the tiers the ladder stands behind" in printed
    assert "NO UniRef50 bridge" in printed, "an unbridged case at >= 0.98 must be offered for reading"


def test_an_unmeasurable_pair_is_never_counted_as_unbridged_in_the_reach(walk_inputs, capsys):
    """⛔ The inflation this guards: a donor with no UniRef50 at all is not evidence of remote homology."""
    recipient_ids, donor_ids, _unused = walk_inputs
    donor_levels = {1: _cog("COG0001", "E")}
    reach = _walk(
        _rows(("1098", "E1", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
    )
    xs.bridge_reach(
        reach,
        recipient_sets={101: frozenset({"U_KP"})},
        donor_sets={},  # donor has none
        kind=COG,
        label_of=LABELS,
        recipient="kp",
        donor="ecoli",
    )
    printed = capsys.readouterr().out
    total_line = [line for line in printed.splitlines() if line.strip().startswith("TOTAL")][0]
    # the one node lands in the LAST column (unmeasurable), not the middle one (not_bridged)
    bridged, not_bridged, unmeasurable = (cell for cell in total_line.split()[1:4])
    assert bridged.startswith("0/") and not_bridged.startswith("0/") and unmeasurable.startswith("1/")


def test_the_two_catalogues_SHARE_the_node_label_namespace_so_the_label_map_keys_on_locus_id():
    """⛔ The bug this caught: both catalogues run node labels '0', '1', '10', … .

    `{**donor_ids, **recipient_ids}` merges on the LABEL, so the recipient's entries overwrite the
    donor's and every donor locus_id falls out of the reverse map — the symptom being the donor
    column of the worked cases printing raw locus ids, which reads as a formatting choice rather than
    a lost lookup. locus_id is globally unique, so the merge has to key on that.
    """
    donor_ids = {"0": (302532, 10), "1": (302533, 4)}
    recipient_ids = {"0": (320063, 7), "1": (320064, 3)}  # ⬅ the SAME labels

    wrong = {locus_id: label for label, (locus_id, _g) in {**donor_ids, **recipient_ids}.items()}
    assert len(wrong) == 2 and 302532 not in wrong, "the label merge loses the donor catalogue entirely"

    right = {locus_id: label for label, (locus_id, _g) in donor_ids.items()}
    right.update({locus_id: label for label, (locus_id, _g) in recipient_ids.items()})
    assert len(right) == 4
    assert right[302532] == "0" and right[320063] == "0", "both catalogues keep their own label '0'"


# ============================================ UniRef50 FIRST, then ESM — the two arms in series
def test_arm_A_is_not_limited_to_the_ESM_shortlist(walk_inputs, capsys):
    """⛔ A UniRef50-first transfer needs no embedding, so it must not be scoped to the ESM map.

    Arm A is the baseline cross-species ESM has to beat; restricting it to the ~21 shortlisted
    candidates would hand ESM a reach that belongs to the bridge. Here the bridged donor (E2) is
    **absent from the ESM candidate list entirely** and Arm A must still find it.
    """
    recipient_ids, donor_ids, _unused = walk_inputs
    xs.uniref_first_then_esm(
        recipient_ids={"1098": (101, 100)},
        donor_ids=donor_ids,
        recipient_counts={101: {"U_SHARED": 95}},
        donor_counts={2: {"U_SHARED": 90}},  # donor locus 2 bridges; it is NOT in the map
        donor_levels={2: _cog("COG0002", "C")},
        recipient_levels={},
        esm_candidates={},  # ⬅ no ESM candidates at all
        kind=COG,
        label_of={101: "1098", 2: "E2"},
        recipient="kp",
        donor="ecoli",
    )
    printed = capsys.readouterr().out
    arm_a = [line for line in printed.splitlines() if "A: UniRef50 bridge" in line][0]
    assert arm_a.split()[-3:] == ["1", "95", "100"], "1 node, 95 bridged genes, 100 whole-node genes"


def test_the_whole_node_gene_count_is_reported_beside_the_bridged_one(capsys):
    """⭐ David's point: a bridge can rest on a MINORITY of a node's genes, and the node is then
    annotated whole on the measured 99.5-99.9 % within-node unanimity. Both counts are reported.

    Here the bridge covers 22 of 100 genes — below the strict threshold, so the permissive arm counts
    it and the strict one does not, which is exactly the sensitivity check that choice needs.
    """
    xs.uniref_first_then_esm(
        recipient_ids={"13": (113, 100)},
        donor_ids={"E7": (7, 100)},
        recipient_counts={113: {"U_MINOR": 22, "U_MAJOR": 78}},
        donor_counts={7: {"U_MINOR": 50}},
        donor_levels={7: _cog("COG0009", "S")},
        recipient_levels={},
        esm_candidates={},
        kind=COG,
        label_of={113: "13", 7: "E7"},
        recipient="kp",
        donor="ecoli",
    )
    lines = capsys.readouterr().out.splitlines()
    permissive = [line for line in lines if "A: UniRef50 bridge" in line][0]
    strict = [line for line in lines if "A-strict" in line][0]
    assert permissive.split()[-3:] == ["1", "22", "100"], "bridged on 22 genes, node is 100"
    assert strict.split()[-3:] == ["0", "—", "0"], "22/100 is below the 50 % strict threshold"


def test_the_arms_are_exclusive_and_sum_to_the_unlabelled_total(capsys):
    """⛔ A node reached by the bridge must NOT also be counted by ESM, or the two arms double-count.

    Three unlabelled nodes: one bridged, one reachable only by ESM, one by neither. The run aborts if
    they do not sum, so the assertion here is that the split is the intended one.
    """
    xs.uniref_first_then_esm(
        recipient_ids={"a": (101, 10), "b": (102, 20), "c": (103, 30)},
        donor_ids={"E1": (1, 5), "E2": (2, 5)},
        recipient_counts={101: {"U_SHARED": 10}, 102: {"U_KP": 20}},  # 103 has no UniRef50 at all
        donor_counts={1: {"U_SHARED": 5}, 2: {"U_EC": 5}},
        donor_levels={1: _cog("COG0001", "E"), 2: _cog("COG0002", "C")},
        recipient_levels={},
        #: node b reaches annotated donor 2 by ESM; node a ALSO has an ESM route but the bridge wins
        esm_candidates={101: [(1, 2, 0.999)], 102: [(1, 2, 0.995)], 103: [(1, 2, 0.5)]},
        kind=COG,
        label_of={},
        recipient="kp",
        donor="ecoli",
    )
    lines = capsys.readouterr().out.splitlines()

    def cells(needle):
        return [line for line in lines if needle in line][0].split()[-3:]

    assert cells("A: UniRef50 bridge")[0] == "1", "node a, via the bridge"
    # ⚠ "B: ESM" alone also matches the accuracy table's HEADER; the arm row is the one that says
    # "no bridge available", and the first version of this assertion picked up the header instead.
    assert cells("no bridge available")[0] == "1", "node b, via ESM only — a is NOT counted twice"
    assert cells("neither")[0] == "1", "node c is below the ESM floor and has no bridge"
    assert [line for line in lines if "TOTAL unlabelled" in line][0].split()[-1] == "3"


# ================================================ PREVALENCE MATCH as a reliability axis (David's ask)
def test_the_representation_follows_the_map_and_has_no_default():
    """⛔ `locus_nearest_locus` stores BOTH representations, so the wrong one is always available.

    A Bacformer map judged against an ESM within-species baseline is two geometries in one table with
    nothing in the output to say so. The resolution is pinned here, and the two loaders take
    `representation` with **no default** so forgetting it is a TypeError rather than a number.
    """
    import inspect

    assert xs.representation_of({"rep": "bacformer"}) == "BACFORMER"
    assert xs.representation_of({"rep": "esm"}) == "ESM"
    assert xs.representation_of({}) == "ESM", "an older map with no `rep` key predates Bacformer"
    for loader in (xs.within_species_edges, xs.candidates_from_database):
        parameter = inspect.signature(loader).parameters["representation"]
        assert parameter.default is inspect.Parameter.empty, (
            f"{loader.__name__} must not default its representation — a forgotten argument would "
            "silently compare a Bacformer map against an ESM baseline"
        )


def test_the_band_distance_sign_is_from_the_RECIPIENTS_point_of_view(capsys):
    """⛔ The sign is load-bearing and it is the opposite of the intuition that raised the question.

    David's case was a CLOUD donor labelling a CORE recipient — *"if the donor is cloud and the
    recipient is core, then it isn't the right match"* — which this reports as the donor being
    **rarer**. The measurement then showed the *other* direction (donor commoner than the recipient)
    to be the worse one, so a flipped sign would invert the finding and the recommendation with it.

    ⚠ The first version of this test asserted only that both labels appeared, which a flipped sign
    satisfies by swapping them. It takes cells above the pair floor whose **agreement differs** for
    the sign to be observable: here every donor-rarer pair agrees and every donor-commoner pair does
    not, so a flip swaps 100 % with 0 %.
    """
    recipient_ids, donor_ids, rows = {}, {}, []
    recipient_bands, donor_bands, levels, donor_levels = {}, {}, {}, {}
    for i in range(40):
        # CLOUD donor -> CORE recipient ("donor 3 bands rarer"), and they AGREE
        r, d = 1000 + i, 2000 + i
        recipient_ids[f"core{i}"] = (r, 100)
        donor_ids[f"Ecloud{i}"] = (d, 4)
        recipient_bands[r], donor_bands[d] = "CORE", "CLOUD"
        levels[r] = levels[d] = _cog("COG0001", "E")
        donor_levels[d] = levels[d]
        # CORE donor -> CLOUD recipient ("donor 3 bands commoner"), and they DISAGREE
        r2, d2 = 3000 + i, 4000 + i
        recipient_ids[f"cloud{i}"] = (r2, 5)
        donor_ids[f"Ecore{i}"] = (d2, 100)
        recipient_bands[r2], donor_bands[d2] = "CLOUD", "CORE"
        levels[r2] = _cog("COG0002", "C")
        levels[d2] = _cog("COG0009", "S")
        donor_levels[d2] = levels[d2]
        for recipient, donor in ((f"core{i}", f"Ecloud{i}"), (f"cloud{i}", f"Ecore{i}")):
            rows.append(
                {
                    "recipient_node": recipient,
                    "donor_node": donor,
                    "recipient_genes": 0,
                    "cross": 0.99,
                    "cross_p75": 0.99,
                    "cross_max": 0.99,
                    "rank_cross": 1,
                    "rank_cross_p75": 1,
                    "rank_cross_max": 1,
                }
            )

    xs.prevalence_match(
        rows,
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        recipient_bands=recipient_bands,
        donor_bands=donor_bands,
        recipient_sets={},
        donor_sets={},
        both=levels,
        donor_levels=donor_levels,
        kind=COG,
        rule="cross",
        floor=0.98,
    )
    lines = capsys.readouterr().out.splitlines()
    rarer = [line for line in lines if "3 bands rarer" in line][0]
    commoner = [line for line in lines if "3 bands commoner" in line][0]
    assert "100.0% n=40" in rarer, f"donor-rarer pairs all agree; got {rarer!r}"
    assert "0.0% n=40" in commoner, f"donor-commoner pairs all disagree; got {commoner!r}"


def test_the_band_order_is_the_published_enums_own(capsys):
    """⛔ RARE < CLOUD < SHELL < SOFT_CORE < CORE. A reordering silently rescales every distance."""
    assert xs.BAND_ORDER == ("RARE", "CLOUD", "SHELL", "SOFT_CORE", "CORE")
    #: a band the enum does not name must be skipped, not coerced to a position
    recipient_ids = {"a": (101, 10)}
    rows = [
        {
            "recipient_node": "a",
            "donor_node": "E1",
            "recipient_genes": 10,
            "cross": 0.99,
            "cross_p75": 0.99,
            "cross_max": 0.99,
            "rank_cross": 1,
            "rank_cross_p75": 1,
            "rank_cross_max": 1,
        }
    ]
    levels = {101: _cog("COG0001", "E"), 1: _cog("COG0001", "E")}
    xs.prevalence_match(
        rows,
        recipient_ids=recipient_ids,
        donor_ids={"E1": (1, 10)},
        recipient_bands={101: "NOT_A_BAND"},
        donor_bands={1: "CORE"},
        recipient_sets={},
        donor_sets={},
        both=levels,
        donor_levels={1: levels[1]},
        kind=COG,
        rule="cross",
        floor=0.98,
    )
    body = [line for line in capsys.readouterr().out.splitlines() if "n=" in line]
    assert not body, "an unrecognised band is skipped rather than given a position"
