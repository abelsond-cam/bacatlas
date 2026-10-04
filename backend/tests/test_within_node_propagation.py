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
  here exactly as it does in the neighbour ladder;
* ⭐ and against a **raw-SQL oracle sharing none of the instrument's code**, which is the only thing
  that can catch the instrument and its author being wrong the same way.

⚠ Only the oracle needs a database; everything else runs without one, because a skip is not a pass.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    ANNOTATED_COUNT_PREDICATE,
    GENE_VALUE_COLUMN,
    SUPPORTED_KINDS,
)
from bacatlas_backend.instruments.within_node_propagation import (
    MIN_NODES,
    PropagationRate,
    _fold_gene,
    compute_propagation,
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
    full = _fold_gene(EC, ["1.1.1.1"])
    assert full[4] == frozenset({"1.1.1.1"})
    assert full[1] == frozenset({"1"})
    # dash-padded: states L3 and above, and NOTHING at L4 — incomparable, not a disagreement
    padded = _fold_gene(EC, ["2.7.10.-"])
    assert padded[4] == frozenset()
    assert padded[3] == frozenset({"2.7.10"})
    # a gene carrying a SET of codes keeps both
    both = _fold_gene(EC, ["1.6.5.9", "7.1.1.-"])
    assert both[1] == frozenset({"1", "7"})


def test_a_COG_gene_folds_at_the_rung_the_page_shows():
    """⚠ `cog_levels`' L1 rung is the LOCUS's category set, which one gene does not have."""
    folded = _fold_gene(COG, "COG3206")
    assert folded[2] == frozenset({"COG3206"})
    assert folded[1] == frozenset(), "a single gene carries no locus-level category set"


def test_every_supported_vocabulary_has_a_gene_column_and_a_derived_predicate():
    """⛔ The mapping that replaced four hand-written copies — a missing kind would skip silently."""
    assert set(GENE_VALUE_COLUMN) == set(SUPPORTED_KINDS)
    assert set(ANNOTATED_COUNT_PREDICATE) == set(SUPPORTED_KINDS)
    for kind, column in GENE_VALUE_COLUMN.items():
        assert ANNOTATED_COUNT_PREDICATE[kind] == f"f.{column} is not null", (
            "the predicate must be DERIVED from the column, or the two can disagree again"
        )


# ======================================================= the oracle, in SQL that shares no code
@pytest.fixture(scope="module")
def session():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the oracle needs a loaded database")
    with Session(create_engine(url, future=True)) as open_session:
        yield open_session


@pytest.mark.parametrize("pangenome_id", [1, 2])
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_every_vocabulary_is_measured_on_ITS_OWN_column(session, pangenome_id, kind):
    """⛔⛔ The defect this instrument replaced, kept as a gate.

    `measure_cog_function_inference.gene_level_values` selected `f.cog_id` and `f.ec_numbers` and
    nothing else, so its GO and KEGG rows measured **COG** agreement over a GO- or KEGG-filtered
    subset of genes. It went unnoticed for a year of edits because the numbers looked entirely
    reasonable — COG agreement is also ~99.5 %, so a wrong measurement of the right shape reads as
    a right one.

    ⭐ **The tell is the POPULATION, not the rate.** A gene carrying a GO term and no COG folds to
    nothing, so its node states no rung at all, falls out of `checkable` AND `one_gene`, and the two
    stop summing to the number of nodes carrying the vocabulary. That identity is exact on both
    catalogues for all four vocabularies, and it is also this module's coverage statement: every
    node the rate is quoted for is in exactly one of the two columns.
    """
    rates = compute_propagation(session, pangenome_id=pangenome_id, annotation_kind=kind)
    carrying = session.execute(
        text(f"""
            select count(distinct m.locus_id)
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where l.pangenome_id = :pangenome_id and {ANNOTATED_COUNT_PREDICATE[kind]}
        """),
        {"pangenome_id": pangenome_id},
    ).scalar_one()
    assert carrying > 0, "non-vacuity: this catalogue does carry this vocabulary"
    assert sum(one.checkable + one.one_gene for one in rates.values()) == carrying


@pytest.mark.parametrize("pangenome_id", [1, 2])
def test_the_COG_rates_match_a_recount_done_entirely_in_SQL(session, pangenome_id):
    """⭐ The anti-vacuity oracle: `count(distinct cog_id) = 1`, which shares none of the folding.

    COG is single-valued per gene, so unanimity reduces to one distinct id — a statement SQL can make
    on its own. ⛔ A frozen oracle cannot see a drift it shares, so this one is a different algorithm,
    not a copied number.
    """
    oracle = {
        row.band: (row.checkable, row.unanimous, row.one_gene)
        for row in session.execute(
            text("""
                with per as (
                  select l.locus_id,
                         l.prevalence_band as band,
                         count(*) as annotated,
                         count(distinct f.cog_id) as distinct_calls
                    from locus l
                    join gene_locus_membership m on m.locus_id = l.locus_id
                    join gene_functional_annotation f
                         on f.genome_id = m.genome_id and f.flat_index = m.flat_index
                   where l.pangenome_id = :pangenome_id and f.cog_id is not null
                   group by 1, 2)
                select band,
                       count(*) filter (where annotated >= 2) as checkable,
                       count(*) filter (where annotated >= 2 and distinct_calls = 1) as unanimous,
                       count(*) filter (where annotated = 1) as one_gene
                  from per group by 1
            """),
            {"pangenome_id": pangenome_id},
        ).all()
    }
    assert sum(counts[0] for counts in oracle.values()) > 0, (
        "the oracle recounted ZERO checkable nodes — the oracle itself is vacuous"
    )

    measured = compute_propagation(session, pangenome_id=pangenome_id, annotation_kind=COG)
    for band, (checkable, unanimous, one_gene) in oracle.items():
        name = band.name if hasattr(band, "name") else str(band)
        assert name in measured, f"the instrument reported no {name} band"
        got = measured[name]
        assert (got.checkable, got.unanimous, got.one_gene) == (checkable, unanimous, one_gene), (
            f"{name}: instrument says {(got.checkable, got.unanimous, got.one_gene)}, "
            f"SQL says {(checkable, unanimous, one_gene)}"
        )


def test_RARE_really_is_the_band_with_nothing_to_measure(session):
    """⛔ The design rests on this, so it is asserted rather than assumed: RARE is near-empty.

    If a catalogue ever arrived where RARE had thousands of checkable nodes, the three-case split on
    the page would be wrong and this test is where that shows up.
    """
    for pangenome_id in (1, 2):
        rates = compute_propagation(session, pangenome_id=pangenome_id, annotation_kind=COG)
        rare = rates.get("RARE")
        assert rare is not None, "RARE should appear — it holds most of the one-gene population"
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
        "true by construction — which is why `compute_propagation` routes it to `one_gene` instead"
    )
