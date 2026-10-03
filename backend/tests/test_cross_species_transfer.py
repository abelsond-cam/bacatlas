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
    """Three unlabelled kp nodes: one reachable, one whose only donor is unannotated, one off the map."""
    recipient_ids = {"1098": (101, 7), "2811": (102, 3), "9999": (103, 5)}
    donor_ids = {"E1": (1, 10), "E2": (2, 4)}
    donor_levels = {1: _cog("COG0001", "E")}  # E2 (locus 2) carries nothing
    return recipient_ids, donor_ids, donor_levels


def _rows(*specs):
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


def test_the_four_states_sum_to_the_unlabelled_total(walk_inputs):
    """⛔ `nuna/CLAUDE.md` §2 — 'no annotated donor' and 'donor agreed' must never read the same."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = xs.walk(
        _rows(("1098", "E2", 0.995, 1), ("1098", "E1", 0.991, 2), ("2811", "E2", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
        kind=COG,
    )
    assert reach["unlabelled"] == 3
    assert reach["state"]["labelled"] == [1, 7], "1098 walks past the unannotated E2 to E1"
    assert reach["state"]["no_annotated_donor"] == [1, 3], "2811's only candidate carries nothing"
    assert reach["state"]["absent_from_map"] == [1, 5], "9999 has no row at all"
    assert reach["state"]["too_remote_to_call"] == [0, 0]
    assert sum(nodes for nodes, _genes in reach["state"].values()) == 3


def test_a_donor_below_the_floor_is_reached_but_NOT_called(walk_inputs):
    """⛔ 'Definitely too remote to call' (David) — the donor is found and no label is taken from it."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = xs.walk(
        _rows(("1098", "E1", 0.85, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
        kind=COG,
    )
    assert reach["state"]["too_remote_to_call"] == [1, 7]
    assert reach["state"]["labelled"] == [0, 0]
    assert reach["state"]["no_annotated_donor"] == [0, 0], "it HAD an annotated donor; it was too far"


def test_a_node_that_already_carries_the_vocabulary_is_not_in_the_denominator(walk_inputs):
    """The two mechanisms are disjoint: a node with its own call is never offered a neighbour's."""
    recipient_ids, donor_ids, donor_levels = walk_inputs
    reach = xs.walk(
        _rows(("1098", "E1", 0.995, 1), ("2811", "E1", 0.995, 1), ("9999", "E1", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={101: _cog("COG0007", "K")},
        kind=COG,
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
    common = dict(
        recipient_ids=recipient_ids, donor_ids=donor_ids, donor_levels=donor_levels, recipient_levels={}, kind=COG
    )
    by_median = xs.walk(rows, rule="cross", **common)
    by_max = xs.walk(rows, rule="cross_max", **common)

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
    reach = xs.walk(
        _rows(("1098", "E1", 0.995, 1)),
        rule="cross",
        recipient_ids=recipient_ids,
        donor_ids=donor_ids,
        donor_levels=donor_levels,
        recipient_levels={},
        kind=COG,
    )
    assert reach["state"]["labelled"][0] == 1
    assert list(reach["by_tier"]) == [(">= 0.99", 1)], "the >= 0.99 tier permits L2; the donor has only L1"
