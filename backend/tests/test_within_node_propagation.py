"""The propagation base rate — the number the Function tab now rests its main claim on.

⭐ The page used to call a one-gene call unverifiable. It now says the call covers the whole syntelogue
and quotes a rate, so that rate has to be right, and it has to be right for the **band it is applied
to** rather than for the catalogue as a whole: nodes where unanimity can be checked have a median size
of 100 genes and nodes resting on one annotated gene have a median of 1.

What is pinned here, and why each one is a failure that would otherwise be silent:

* ⛔ `RARE` returns `None`, never 0.0 — it has essentially no checkable nodes and holds most of the
  one-gene population, so a 0.0 would print "0 % agree" on the very nodes the page is quietest about;
* ⛔ a band under the floor reports its count and no rate, like `annotation_transfer.Cell`;
* ⚠ the folding goes through the **shared** folders, so a dash-padded or comma-joined EC value behaves
  here exactly as it does in the neighbour ladder.

⚠ The raw-SQL oracle over the gene rows went with the per-gene tables in D8; the stored tally is
checked against the stored verdicts in `tests/test_locus_annotation_rollups.py`. One database check
remains here: RARE is the band with nothing to measure, read from `pangenome_propagation_rate`.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.within_node_propagation import (
    MIN_NODES,
    PropagationRate,
    _fold_gene,
    unanimous_at_deepest_shared_rung,
)
from bacatlas_backend.models.enumerations import AnnotationKind

COG = AnnotationKind.COG_ORTHOGROUP
EC = AnnotationKind.EC_NUMBER


def rate_of(**overrides) -> PropagationRate:
    fields = {
        "annotation_kind": COG,
        "prevalence_band": "CORE",
        "checkable": 1000,
        "unanimous": 998,
        "one_gene": 40,
    }
    return PropagationRate(**{**fields, **overrides})


# ============================================================ the floor, and what NULL has to mean
def test_a_band_with_no_checkable_nodes_reports_None_and_never_zero():
    """⛔ RARE is this case, and it holds MOST one-gene nodes — 367 kp / 410 ecoli on COG.

    A 0.0 would print "0 % of these agree" on exactly the nodes the page should be quietest about,
    and it would be read as a measurement. `None` says *not measurable*, which is what is true.
    """
    empty = rate_of(prevalence_band="RARE", checkable=0, unanimous=0, one_gene=367)
    assert empty.rate is None
    assert empty.interval is None
    assert empty.as_json()["rate"] is None
    # ⛔ and the one-gene population is still reported, so the silence is visible rather than absent
    assert empty.as_json()["one_gene_node_count"] == 367
    assert empty.as_json()["checkable_node_count"] == 0


def test_a_band_under_the_floor_reports_its_count_and_no_rate():
    """The same stance as `annotation_transfer.Cell.agreement`, and the same number."""
    assert rate_of(checkable=MIN_NODES - 1, unanimous=MIN_NODES - 1).rate is None
    at_floor = rate_of(checkable=MIN_NODES, unanimous=MIN_NODES)
    assert at_floor.rate == 1.0, "the floor is inclusive, as MIN_PAIRS is"
    assert at_floor.interval is not None


def test_the_interval_stays_inside_zero_and_one_at_the_rates_these_bands_reach():
    """⚠ Why Wilson and not a normal interval: every real band here is above 99 %.

    A normal interval at 276/276 runs past 1.0 and prints as nonsense on the page.
    """
    perfect = rate_of(checkable=276, unanimous=276)
    low, high = perfect.interval
    assert high <= 1.0 and low > 0.98
    # the measured kp CORE cell, as a regression on the arithmetic itself
    core = rate_of(checkable=3253, unanimous=3248)
    assert core.rate == pytest.approx(0.99846, abs=1e-5)
    assert 0.996 < core.interval[0] < core.rate < core.interval[1] < 1.0


# =================================================== the folding goes through the SHARED folders
def test_an_EC_value_folds_exactly_as_the_neighbour_ladder_folds_it():
    """⚠ EC values are comma-joined SETS, dash-padded to an unknown depth — `gumC` is `2.7.10.-`.

    If this module folded them itself, the rung the page quotes for propagation could differ from the
    rung it quotes for transfer, on the same locus, in the same card.
    """
    # ⚠ A rung is `(facet, depth)` since 2026-10-06. EC has one axis, so its facet is None.
    full = _fold_gene(EC, ["1.1.1.1"])
    assert full[(None, 4)] == frozenset({"1.1.1.1"})
    assert full[(None, 1)] == frozenset({"1"})
    # dash-padded: states L3 and above, and NOTHING at L4 — incomparable, not a disagreement
    padded = _fold_gene(EC, ["2.7.10.-"])
    assert padded[(None, 4)] == frozenset()
    assert padded[(None, 3)] == frozenset({"2.7.10"})
    # a gene carrying a SET of codes keeps both
    both = _fold_gene(EC, ["1.6.5.9", "7.1.1.-"])
    assert both[(None, 1)] == frozenset({"1", "7"})


def test_a_COG_gene_folds_at_the_rung_the_page_shows():
    """⚠ `cog_levels`' L1 rung is the LOCUS's category set, which one gene does not have."""
    folded = _fold_gene(COG, "COG3206")
    assert folded[(None, 2)] == frozenset({"COG3206"})
    assert folded[(None, 1)] == frozenset(), "a single gene carries no locus-level category set"


# ================================================= the stored tally, where a database is loaded
@pytest.fixture(scope="module")
def session():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the oracle needs a loaded database")
    with Session(create_engine(url, future=True)) as open_session:
        yield open_session


def test_RARE_really_is_the_band_with_nothing_to_measure(session):
    """⛔ The design rests on this, so it is asserted rather than assumed: RARE is near-empty.

    If a catalogue ever arrived where RARE had thousands of checkable nodes, the three-case split on
    the page would be wrong and this test is where that shows up.
    """
    for pangenome_id in (1, 2):
        row = session.execute(
            text("select checkable, unanimous, one_gene from pangenome_propagation_rate "
                 "where pangenome_id = :pangenome_id and annotation_kind = :kind and prevalence_band = 'RARE'"),
            {"pangenome_id": pangenome_id, "kind": COG.name},
        ).one_or_none()
        assert row is not None, "RARE should appear — it holds most of the one-gene population"
        rare = rate_of(prevalence_band="RARE", checkable=row.checkable, unanimous=row.unanimous, one_gene=row.one_gene)
        assert rare.checkable < MIN_NODES, f"RARE now has {rare.checkable} checkable nodes"
        assert rare.rate is None
        assert rare.one_gene > 100, f"RARE holds only {rare.one_gene} one-gene nodes"


# ============================================ the rung the decision is taken at, without a database
def test_the_deepest_shared_rung_is_used_and_NOT_the_shallowest():
    """⛔ The two choices give OPPOSITE answers, which is why this is pinned.

    `1.1.1.1` and `1.1.1.2` share `1.1.1` at L3 and differ at L4. Taking the deepest rung both state
    calls it a disagreement, which it is; taking the shallowest would call it agreement and silently
    inflate every EC rate on the page.
    """
    assert unanimous_at_deepest_shared_rung(EC, [["1.1.1.1"], ["1.1.1.2"]]) is False
    assert unanimous_at_deepest_shared_rung(EC, [["1.1.1.1"], ["1.1.1.1"]]) is True
    # ⚠ and a genuine fallback still works: both are dash-padded, so L3 is the deepest either states
    assert unanimous_at_deepest_shared_rung(EC, [["2.7.10.-"], ["2.7.10.-"]]) is True
    assert unanimous_at_deepest_shared_rung(EC, [["2.7.10.-"], ["2.7.11.-"]]) is False


def test_no_shared_rung_is_counted_in_NEITHER_column():
    """⛔ Scoring "nothing to compare" as agreement is how a measurement flatters itself.

    One gene states L4 and the other only L1, so no rung is stated by both.
    """
    assert unanimous_at_deepest_shared_rung(EC, [["1.1.1.1"], ["2.-.-.-"]]) is False, (
        "L1 IS shared here — both state a class, and they differ"
    )
    # a value that folds to nothing at all cannot be compared with anything
    assert unanimous_at_deepest_shared_rung(EC, [["1.1.1.1"], []]) is None


def test_a_single_annotated_gene_is_never_asked_whether_it_agrees_with_itself():
    """⚠ The one-gene case is counted separately upstream and never reaches the decision.

    Asked directly it would answer True by construction — which is the claim this whole module exists
    to stop the page making.
    """
    assert unanimous_at_deepest_shared_rung(COG, ["COG3206"]) is True, (
        "true by construction — which is why `propagation_verdict` routes it to `one_gene` instead"
    )
