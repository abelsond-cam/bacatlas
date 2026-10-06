"""The Function tab's inference block — the ladder, the walk, and the rate it quotes.

⚠ Runs against `BACATLAS_DATABASE_URL` with both catalogues loaded and published, like
`test_api_endpoints.py`. The measurement is a property of the catalogue, so a fixture cannot stand in
for it: the whole point of this module is that the rate the page shows is **recomputed from the
loaded data** rather than copied out of
`nuna/docs/model_evaluation/function_inference.md`.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, func, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.instruments.annotation_transfer import (
    CALLING_FLOOR,
    MIN_PAIRS,
    NOT_CALLED,
    QUOTED_LEVEL,
    TIERS,
    Calibration,
    Cell,
    callable_level,
    cog_levels,
    ec_levels,
    quotable_level,
    tier_for,
)
from bacatlas_backend.instruments.within_node_propagation import MIN_NODES
from bacatlas_backend.models.enumerations import AnnotationKind, EmbeddingRepresentation, PrevalenceBand
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.services.function_inference_service import (
    _PROPAGATION,
    calibration_for,
    clear_calibration_cache,
    clear_propagation_cache,
    propagation_for,
)


@pytest.fixture(scope="module")
def database_url():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the inference tests need a loaded database")
    return url


@pytest.fixture(scope="module")
def application(database_url):
    with Session(create_engine(database_url, future=True)) as session:
        published = session.execute(
            select(func.count())
            .select_from(PathogenSpecies)
            .where(PathogenSpecies.default_pangenome_id.is_not(None))
        ).scalar_one()
    if published < 2:
        pytest.skip(f"only {published} species are published in {database_url}")
    return create_application(Configuration(database_url=database_url))


@pytest.fixture(scope="module")
def client(application):
    return application.test_client()


@pytest.fixture(scope="module")
def session(database_url):
    with Session(create_engine(database_url, future=True)) as open_session:
        yield open_session


# ── the ladder is a DECISION, so it is pinned ───────────────────────────────────────────────────
def test_the_ladder_is_the_one_David_chose_and_a_change_to_it_must_be_deliberate():
    """⛔ David, 2026-10-03 and 2026-10-06 — nuna `PROJECT_STATE.md` §6. Not a tuning parameter.

    Each boundary was argued from a measurement: nothing below 0.90 appears at all ("definitely too
    remote to call"); EC quotes its full code only in the top tier, dropping a level at 0.98-0.99
    where the full code is 84 % and the sub-subclass 95 %.

    ⭐ **The bottom tier split at 0.94 on 2026-10-06, and this test is what made it deliberate** — it
    failed the moment `TIERS` changed, which is the whole reason it is written as an equality against
    a literal rather than as a property. It had been one band, 0.90-0.96, because 0.95-0.96 and
    0.90-0.95 were statistically indistinguishable; that was true of the GO claim as it was then
    compared, and stopped being true once a GO claim was closed upward through the slim. Measured on
    the closed claim, GO steps at 0.94 in BOTH species — ecoli 73.2 % → 80.2 %, kp 72.7 % → 80.0 %
    — so `0.94-0.96` reads 82.2 % / 84.0 % and clears the calling floor while `0.90-0.94` reads
    71.2 % / 72.6 % and does not.

    ⚠ `TIERS` is shared by all four vocabularies, so the split re-cut every bottom cell. Checked
    rather than assumed: COG tops out at 35.8 % / 37.9 %, EC at 40.0 % / 43.5 %, KEGG at 25.9 % /
    34.5 %, so both halves go on refusing and the split costs them nothing.
    """
    assert [name for name, _, _ in TIERS] == [
        ">= 0.99", "0.98-0.99", "0.97-0.98", "0.96-0.97", "0.94-0.96", "0.90-0.94",
    ]
    assert [low for _, low, _ in TIERS] == [0.99, 0.98, 0.97, 0.96, 0.94, 0.90]
    assert tier_for(0.8999) == NOT_CALLED, "below 0.90 is never called"
    assert tier_for(None) == NOT_CALLED, "a missing cosine is not a tier"
    # ⚠ The new boundary itself, both sides of it — half-open, so 0.94 is the STRONGER tier.
    assert tier_for(0.94) == "0.94-0.96"
    assert tier_for(0.9399) == "0.90-0.94"
    assert QUOTED_LEVEL[AnnotationKind.EC_NUMBER] == {
        ">= 0.99": 4, "0.98-0.99": 3, "0.97-0.98": 1, "0.96-0.97": 1, "0.94-0.96": 1, "0.90-0.94": 1,
    }
    assert QUOTED_LEVEL[AnnotationKind.COG_ORTHOGROUP] == {
        ">= 0.99": 2, "0.98-0.99": 2, "0.97-0.98": 1, "0.96-0.97": 1, "0.94-0.96": 1, "0.90-0.94": 1,
    }


def test_an_EC_value_is_a_SET_of_codes_each_resolved_to_its_own_depth():
    """The two shapes that make `split_part` on the string return a code neither side holds.

    Measured in the published catalogues: 1,034 of 9,398 values are comma-joined lists and 2,257
    carry `-` placeholders.
    """
    plain = ec_levels("2.7.10.1")
    assert plain[4] == {"2.7.10.1"} and plain[1] == {"2"}

    padded = ec_levels("3.1.-.-")
    assert padded[2] == {"3.1"}
    assert padded[3] == frozenset() and padded[4] == frozenset(), (
        "a dash means the SOURCE does not state this level — incomparable, not a disagreement"
    )

    joined = ec_levels("1.6.5.9,7.1.1.-")
    assert joined[4] == {"1.6.5.9"}, "only the complete code contributes at level 4"
    assert joined[3] == {"1.6.5", "7.1.1"} and joined[1] == {"1", "7"}
    assert "9,7" not in joined[4], "⛔ the `split_part` bug: a code neither side holds"


def test_a_tier_quotes_no_deeper_than_the_DONOR_actually_states():
    """⛔ Of 84 ecoli nodes with a >= 0.99 EC neighbour only 43 can be given a four-field code."""
    complete = ec_levels("2.7.10.1")
    assert quotable_level(AnnotationKind.EC_NUMBER, ">= 0.99", complete) == 4

    padded = ec_levels("2.7.-.-")
    assert quotable_level(AnnotationKind.EC_NUMBER, ">= 0.99", padded) == 2, (
        "the tier permits 4, the donor states 2, so 2 is quoted — not 4"
    )
    assert quotable_level(AnnotationKind.EC_NUMBER, "0.98-0.99", padded) == 2

    assert quotable_level(AnnotationKind.EC_NUMBER, NOT_CALLED, complete) is None
    assert quotable_level(AnnotationKind.COG_ORTHOGROUP, ">= 0.99",
                          cog_levels(None, ["M"])) == 1, (
        "a donor with categories but no orthogroup still has a category to give"
    )
    assert quotable_level(AnnotationKind.COG_ORTHOGROUP, ">= 0.99",
                          cog_levels(None, None)) is None


# ── ⛔ the floor: whether to speak at all, which is NOT how deep ─────────────────────────────────
def _ladder(kind, rates: dict[tuple[str, int], float | None]) -> Calibration:
    """A hand-built ladder. ⚠ `None` means a cell under `MIN_PAIRS` — a count and no rate."""
    cells = {}
    for (tier, level), rate in rates.items():
        pairs = 1_000 if rate is not None else MIN_PAIRS - 1
        agreeing = round((rate or 0) * pairs)
        cells[(tier, level)] = Cell(tier=tier, level=level, level_label="x", pairs=pairs,
                                    agreeing=agreeing, chance=0.01)
    return Calibration(annotation_kind=kind, representation=EmbeddingRepresentation.ESM, cells=cells,
                       annotated_locus_count=1, locus_count=1)


def test_the_floor_is_a_SEPARATE_decision_from_the_depth():
    """⛔⛔ David, 2026-10-05. `QUOTED_LEVEL` says how deep a tier may be read; it never said whether
    the call was worth making, and nothing else did either — so `gumC` was offered a KEGG orthology
    from a 0.929 neighbour where agreement is **24.3 %**.

    ⭐ The floor is applied as one more rung of the SAME fallback, not as a gate beside it: a level
    nobody can be quoted at accurately is a level the donor cannot usefully supply.
    """
    assert CALLING_FLOOR == 0.80
    complete = ec_levels("2.7.10.1")

    # every rung of the >= 0.99 EC cell clears the floor, so the full code still goes out
    rich = _ladder(AnnotationKind.EC_NUMBER, {(">= 0.99", level): 0.99 for level in (1, 2, 3, 4)})
    assert callable_level(AnnotationKind.EC_NUMBER, ">= 0.99", complete, rich) == 4
    assert quotable_level(AnnotationKind.EC_NUMBER, ">= 0.99", complete) == 4, "unchanged by the floor"

    # ⭐ the deepest rung misses it and a shallower one clears: quote the shallower one
    mixed = _ladder(AnnotationKind.EC_NUMBER,
                    {(">= 0.99", 4): 0.60, (">= 0.99", 3): 0.97, (">= 0.99", 2): 0.99, (">= 0.99", 1): 0.99})
    assert callable_level(AnnotationKind.EC_NUMBER, ">= 0.99", complete, mixed) == 3

    # ⛔ COG category at 0.97-0.98 is 78.1 % ecoli / 79.9 % kp with NO shallower rung — say nothing
    thin = _ladder(AnnotationKind.COG_ORTHOGROUP, {("0.97-0.98", 1): 0.781})
    assert callable_level(AnnotationKind.COG_ORTHOGROUP, "0.97-0.98", cog_levels("COG1", ["J"]), thin) is None

    # ⛔ a cell with no measured rate cannot clear a floor. The card used to print the value anyway.
    unmeasured = _ladder(AnnotationKind.KEGG_ORTHOLOGY, {("0.90-0.94", 1): None})
    assert callable_level(AnnotationKind.KEGG_ORTHOLOGY, "0.90-0.94",
                          {1: frozenset({"K01991"})}, unmeasured) is None
    # and a cell that is simply absent
    assert callable_level(AnnotationKind.KEGG_ORTHOLOGY, "0.96-0.97",
                          {1: frozenset({"K01991"})}, unmeasured) is None

    assert callable_level(AnnotationKind.EC_NUMBER, NOT_CALLED, complete, rich) is None


def test_gumC_is_no_longer_offered_a_KEGG_it_would_be_wrong_about(client):
    """⭐ The worked case for the floor: `wza` is a rank-5 neighbour at **0.9293** — close, carrying
    a KEGG — and kp KEGG agreement at that tier is **24.3 % over 481 pairs**.

    ⛔ The donor is still NAMED and the rate still carried, because "too remote to call" would be
    false here and the reader is owed the number that refused it.
    """
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
    candidate = by_kind["kegg_orthology"]["candidate"]

    assert candidate is not None, "the donor is still found and still reported"
    assert candidate["donor"]["display_name"] == "wza" and candidate["rank"] == 5
    assert candidate["cosine"] > 0.92, "non-vacuity: this donor is CLOSE, not remote"
    assert candidate["tier"] != NOT_CALLED
    assert candidate["level"] is None and candidate["value"] is None, "⛔ nothing is suggested"
    assert candidate["calibration"]["agreement"] < CALLING_FLOOR, (
        "and the measured rate that refused it travels with the refusal"
    )
    assert payload["calibration"]["kegg_orthology"]["calling_floor"] == CALLING_FLOOR


def test_the_floor_is_not_vacuous_a_STRONG_neighbour_still_suggests(client, session):
    """⛔ A floor that suppressed everything would pass every test above and serve a blank page."""
    label = session.execute(
        text("""
            select focal.node_label
              from locus focal
              join pangenome p on p.pangenome_id = focal.pangenome_id
              join locus_nearest_locus n
                   on n.locus_id = focal.locus_id and n.representation = 'ESM'
              join locus donor on donor.locus_id = n.neighbour_locus_id
             where p.catalogue_key = 'kp-nuna4'
               and focal.cog_distinct_id_count = 0
               and donor.cog_distinct_id_count > 0
               and n.cross_similarity >= 0.99
               and n.rank = 1
             order by focal.node_label limit 1
        """)
    ).scalar()
    assert label is not None, "the catalogue has >= 0.99 COG donors, or this test is vacuous"
    payload = client.get(f"/api/v1/species/kp/loci/{label}/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
    candidate = by_kind["cog_orthogroup"]["candidate"]
    assert candidate["level"] is not None and candidate["value"], (
        "a >= 0.99 COG neighbour agrees 99.4 % of the time and is still quoted"
    )


# ── the calibration is MEASURED, and an independent oracle says so ──────────────────────────────
def test_the_COG_calibration_matches_a_recount_done_entirely_in_SQL(session):
    """⭐ The anti-vacuity gate: an oracle that shares none of the instrument's code.

    COG level 2 is the one rung whose agreement is plain string equality, so SQL can recount it
    without reimplementing the level folding. If `compute_calibration` silently stopped counting —
    returning empty cells, or every pair as agreeing — this is what notices.
    """
    clear_calibration_cache()
    calibration = calibration_for(
        session, pangenome_id=1, annotation_kind=AnnotationKind.COG_ORTHOGROUP
    )
    for tier, low, high in TIERS:
        recounted = session.execute(
            text("""
                select count(*) as pairs,
                       count(*) filter (where mine.term_value = theirs.term_value) as agreeing
                  from locus_nearest_locus n
                  join locus focal on focal.locus_id = n.locus_id
                  join locus_annotation_entry mine
                       on mine.locus_id = n.locus_id
                      and mine.annotation_kind = 'COG_ORTHOGROUP' and mine.rank_within_locus = 0
                  join locus_annotation_entry theirs
                       on theirs.locus_id = n.neighbour_locus_id
                      and theirs.annotation_kind = 'COG_ORTHOGROUP' and theirs.rank_within_locus = 0
                 where focal.pangenome_id = 1 and n.representation = 'ESM'
                   and n.cross_similarity >= :low and n.cross_similarity < :high
            """),
            {"low": low, "high": high},
        ).one()
        cell = calibration.cell(tier, 2)
        assert cell is not None, f"{tier} has no level-2 cell at all"
        assert recounted.pairs > 0, f"{tier} recounted ZERO pairs — the oracle itself is vacuous"
        assert (cell.pairs, cell.agreeing) == (recounted.pairs, recounted.agreeing), tier
        assert 0.0 < cell.agreement < 1.0 or cell.agreement in (1.0,), tier


def test_chance_is_EXACT_and_not_a_sampled_null(session):
    """`chance` must be the share of OTHER annotated nodes whose set intersects the focal node's.

    ⚠ Checked on COG level 2, where a claim set holds exactly one value, so the quantity has a
    closed form SQL can state: `(holders - 1) / (pool - 1)`, averaged over the cell's pairs. A
    sampled null would land near this and not on it, which is precisely why the sampled version was
    replaced — the page and the write-up disagreed in the third digit.
    """
    clear_calibration_cache()
    calibration = calibration_for(
        session, pangenome_id=1, annotation_kind=AnnotationKind.COG_ORTHOGROUP
    )
    tier, low, high = TIERS[0]
    expected = session.execute(
        text("""
            with holders as (
                select a.term_value, count(*) as n
                  from locus_annotation_entry a join locus l on l.locus_id = a.locus_id
                 where l.pangenome_id = 1 and a.annotation_kind = 'COG_ORTHOGROUP'
                   and a.rank_within_locus = 0
                 group by 1),
            pool as (select sum(n) as n from holders)
            select avg((holders.n - 1)::float / (pool.n - 1)) as chance
              from locus_nearest_locus n
              join locus focal on focal.locus_id = n.locus_id
              join locus_annotation_entry mine
                   on mine.locus_id = n.locus_id
                  and mine.annotation_kind = 'COG_ORTHOGROUP' and mine.rank_within_locus = 0
              join locus_annotation_entry theirs
                   on theirs.locus_id = n.neighbour_locus_id
                  and theirs.annotation_kind = 'COG_ORTHOGROUP' and theirs.rank_within_locus = 0
              join holders on holders.term_value = mine.term_value
              cross join pool
             where focal.pangenome_id = 1 and n.representation = 'ESM'
               and n.cross_similarity >= :low and n.cross_similarity < :high
        """),
        {"low": low, "high": high},
    ).scalar_one()
    cell = calibration.cell(tier, 2)
    assert cell.chance == pytest.approx(expected, rel=1e-9)
    assert cell.lift > 100, "the top tier is hundreds of times better than chance, not marginally"


def test_a_cell_below_the_pair_floor_reports_its_n_and_NO_rate(session):
    """⛔ `agreement` is `None`, never a number, below `MIN_PAIRS` — a rate off 7 pairs is not a rate."""
    clear_calibration_cache()
    for pangenome_id in (1, 2):
        for kind in (AnnotationKind.COG_ORTHOGROUP, AnnotationKind.EC_NUMBER):
            calibration = calibration_for(session, pangenome_id=pangenome_id, annotation_kind=kind)
            assert calibration.cells, f"{pangenome_id}/{kind} produced no cells at all"
            for cell in calibration.cells.values():
                assert cell.pairs > 0
                if cell.pairs < MIN_PAIRS:
                    assert cell.agreement is None and cell.interval is None and cell.lift is None
                else:
                    assert cell.agreement is not None and cell.interval is not None
                    low, high = cell.interval
                    assert low <= cell.agreement <= high


def test_the_two_catalogues_do_not_share_a_cached_ladder(session):
    """⚠ The cache is keyed per catalogue, so one species' rate can never answer for the other."""
    clear_calibration_cache()
    ecoli = calibration_for(session, pangenome_id=1, annotation_kind=AnnotationKind.COG_ORTHOGROUP)
    kp = calibration_for(session, pangenome_id=2, annotation_kind=AnnotationKind.COG_ORTHOGROUP)
    assert ecoli.annotated_locus_count != kp.annotated_locus_count
    assert ecoli.cell(">= 0.99", 2).pairs != kp.cell(">= 0.99", 2).pairs


# ── the endpoint, on the cases that teach the caveats ───────────────────────────────────────────
def test_gumC_says_its_ONE_COG_call_PROPAGATES_to_the_whole_syntelogue(client, session):
    """⭐ The worked case, and the one the page had BACKWARDS: kp node 722, 100 members, COG3206 on
    **1** of them, EC on 16.

    ⛔ The two vocabularies still must not read alike on the *check*: one vote had sixteen voters,
    the other one, and `checkable` is what keeps that apart. ⭐ **But the CLAIM is the same in both
    cases** — the call covers all 100 member genes either way, at the measured rate for the node's
    prevalence band. David, 2026-10-04: *"the whole point of our method is effectively the carrying
    by one is enough."* This test previously asserted the opposite copy ("cannot be checked"), which
    is why it is rewritten here rather than deleted.
    """
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}

    cog = by_kind["cog_orthogroup"]["own"]
    assert cog["term"] == "COG3206"
    assert (cog["gene_count"], cog["annotated_gene_count"], cog["member_gene_count"]) == (1, 1, 100)
    assert cog["checkable"] is False, "one voter is still not a vote, and the page still says so"

    # ⛔ COVERAGE BEFORE THE RATE. The two populations must account for every kp CORE node that
    # carries a COG at all — `checkable` is what the rate is MEASURED on, `one_gene` what it is
    # APPLIED to, and a rate quoted without them would be a number from an unnamed denominator.
    band = cog["propagation"]
    assert band["prevalence_band"] == "core", "gumC is a 100-genome core node, not a singleton"
    carrying = session.execute(
        text("""
            select count(*)
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = 'kp-nuna4' and l.prevalence_band = 'CORE'
               and l.cog_distinct_id_count > 0
        """)
    ).scalar_one()
    assert band["checkable_node_count"] + band["one_gene_node_count"] == carrying, (
        "every CORE node carrying a COG is in exactly one of the two columns"
    )
    assert band["checkable_node_count"] >= MIN_NODES, "and there are enough of them to quote a rate"
    assert band["rate"] is not None and band["rate"] > 0.99
    assert band["interval_low"] <= band["rate"] <= band["interval_high"] <= 1.0, (
        "⛔ a Wilson interval, because a normal one runs past 1.0 at these rates"
    )

    ec = by_kind["ec_number"]["own"]
    assert ec["term"] == "2.7.10.-"
    assert ec["annotated_gene_count"] == 16 and ec["member_gene_count"] == 100
    assert ec["checkable"] is True
    assert ec_levels(ec["term"])[4] == frozenset(), (
        "and its EC states only three levels however similar a reader's neighbour is"
    )
    # ⭐ The sixteen-voter vocabulary carries the SAME kind of statement as the one-voter one —
    # only its first clause differs. That is the A1 decision, and it is the whole fix.
    assert ec["propagation"]["prevalence_band"] == "core"
    assert ec["propagation"]["rate"] is not None


# ── the TWO mechanisms, and the populations they serve ──────────────────────────────────────────
@pytest.mark.parametrize("catalogue,pangenome_id", [("ecoli-nuna4", 1), ("kp-nuna4", 2)])
def test_the_two_mechanisms_PARTITION_the_catalogue_with_nothing_left_over(
    client, session, catalogue, pangenome_id
):
    """⛔ Step 1 (its own genes) and step 2 (a neighbour) are DISJOINT and together exhaust the
    catalogue — the coverage statement the two-mechanism claim rests on.

    David, 2026-10-04: *"1/ Impute within the node… THEN 2/ Impute between nodes."* The service runs
    the walk only where `own is None`, so a node is served by exactly one mechanism; if these
    populations did not sum, some nodes would be served by neither and nothing would say so.
    """
    counted = session.execute(
        text("""
            select count(*) as loci,
                   count(*) filter (where carries) as step_one,
                   count(*) filter (where not carries) as step_two
              from (select exists (
                       select 1 from locus_annotation_entry e
                        where e.locus_id = l.locus_id
                          and e.annotation_kind = 'COG_ORTHOGROUP'
                          and e.rank_within_locus = 0) as carries
                      from locus l where l.pangenome_id = :pangenome_id) as node
        """),
        {"pangenome_id": pangenome_id},
    ).one()
    assert counted.step_one > 0 and counted.step_two > 0, (
        "non-vacuity: both populations are non-empty in this catalogue"
    )
    assert counted.step_one + counted.step_two == counted.loci

    # ⭐ And the SERVICE agrees with that split, which is the half SQL cannot assert: one node from
    # each population, served by exactly one mechanism and never by both or neither.
    for carries in (True, False):
        label = session.execute(
            text(f"""
                select l.node_label from locus l
                 where l.pangenome_id = :pangenome_id
                   and {'' if carries else 'not'} exists (
                       select 1 from locus_annotation_entry e
                        where e.locus_id = l.locus_id
                          and e.annotation_kind = 'COG_ORTHOGROUP'
                          and e.rank_within_locus = 0)
                 order by l.member_gene_count desc limit 1
            """),
            {"pangenome_id": pangenome_id},
        ).scalar_one()
        species = "ecoli" if catalogue.startswith("ecoli") else "kp"
        payload = client.get(f"/api/v1/species/{species}/loci/{label}/function").get_json()
        served = next(
            entry for entry in payload["inference"]["vocabularies"]
            if entry["annotation_kind"] == "cog_orthogroup"
        )
        assert (served["own"] is not None) is carries
        assert (served["walk"] == []) is carries, (
            "⛔ a node with its own call is never walked, and one without it always is"
        )


def test_a_donor_CARRIES_a_vocabulary_by_the_ANY_GENE_rule_not_by_a_medoid(client, session):
    """⭐ A node whose COG rests on ONE of its hundred genes is a valid donor, and is used.

    ⛔ **A test rather than a comment.** The medoid lookup is gone from `services/`, `api/` and
    `annotation_transfer.py` — every remaining mention records its removal — but a removal recorded
    only in prose is a removal that can come back. If a donor had to carry the call at its medoid, a
    1-of-100 donor would be skipped and the walk would run on to a further neighbour, which is a
    silent change in what the page suggests rather than an error.
    """
    row = session.execute(
        text("""
            select focal.node_label as focal, donor.node_label as donor,
                   donor.cog_annotated_member_count as annotated,
                   donor.member_gene_count as members
              from locus focal
              join pangenome p on p.pangenome_id = focal.pangenome_id
              join locus_nearest_locus n
                   on n.locus_id = focal.locus_id and n.representation = 'ESM'
              join locus donor on donor.locus_id = n.neighbour_locus_id
             where p.catalogue_key = 'kp-nuna4'
               and focal.cog_distinct_id_count = 0
               and donor.cog_distinct_id_count > 0
               and donor.cog_annotated_member_count = 1
               and donor.member_gene_count >= 10
               and not exists (
                     select 1 from locus_nearest_locus earlier
                       join locus e on e.locus_id = earlier.neighbour_locus_id
                      where earlier.locus_id = focal.locus_id
                        and earlier.representation = 'ESM'
                        and earlier.rank < n.rank
                        and e.cog_distinct_id_count > 0)
             order by focal.node_label limit 1
        """)
    ).first()
    assert row is not None, "the catalogue has one-gene donors, or this test is vacuous"
    assert row.annotated == 1 and row.members >= 10

    payload = client.get(f"/api/v1/species/kp/loci/{row.focal}/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
    candidate = by_kind["cog_orthogroup"]["candidate"]
    assert candidate is not None, "⛔ a 1-of-100 donor was skipped — the any-gene rule is not in force"
    assert candidate["donor"]["node_label"] == row.donor
    assert candidate["donor"]["annotated_gene_count"] == 1
    assert candidate["donor"]["checkable"] is False, "and it is still marked as unchecked"
    assert candidate["donor"]["propagation"]["rate"] is not None, (
        "⭐ which is exactly why it is a valid donor: its call covers its node at a measured rate"
    )


def test_the_band_the_propagation_names_is_the_SAME_STRING_the_locus_view_serves(client):
    """⛔ One field, one shape, across two routes.

    The instrument keys its tally on the enum NAME, because that is what Postgres stores; every
    route that has ever served `prevalence_band` serves the enum VALUE, and the page reads all of
    them through one `bandLabel`. A second shape for the same field is invisible on the wire and
    renders `CORE` raw beside `soft core` on the first page that forgets.
    """
    locus = client.get("/api/v1/species/kp/loci/722").get_json()["locus"]["prevalence_band"]
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    stated = {
        entry["own"]["propagation"]["prevalence_band"]
        for entry in payload["inference"]["vocabularies"]
        if entry["own"] is not None
    }
    assert stated, "non-vacuity: this node does carry calls to propagate"
    assert stated == {locus}, "the Function tab and the locus view name the band identically"


def test_the_band_rate_does_NOT_depend_on_how_many_genes_carried_the_call(client, session):
    """⭐ One node resting on ONE COG and one resting on many quote the IDENTICAL band cell.

    ⛔ The rate answers "is the modal call right for the members that do not carry it", and that
    question does not change with the number that do. A rate that moved between the two cases would
    mean the page was quoting a property of the node rather than of its reference class.
    """
    crowded = session.execute(
        text("""
            select l.node_label
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = 'kp-nuna4' and l.prevalence_band = 'CORE'
               and l.cog_annotated_member_count >= 50
             order by l.node_label limit 1
        """)
    ).scalar_one()

    def cog_own(label):
        payload = client.get(f"/api/v1/species/kp/loci/{label}/function").get_json()
        by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
        return by_kind["cog_orthogroup"]["own"]

    one_gene, many = cog_own("722"), cog_own(crowded)
    assert one_gene["annotated_gene_count"] == 1 and many["annotated_gene_count"] >= 50, (
        "non-vacuity: the two nodes really are the two cases"
    )
    assert one_gene["checkable"] is False and many["checkable"] is True
    assert one_gene["propagation"] == many["propagation"], "the same reference class, the same rate"


def test_a_RARE_node_is_told_the_rate_is_NOT_MEASURABLE_and_never_zero(client, session):
    """⛔⛔ The distinction the whole three-state design exists for.

    RARE holds most of the one-gene nodes (367 of them on kp COG) and has **zero** checkable ones,
    because a singleton node has one gene to check. `0.0` there would read as *measured, and the
    label propagates nowhere* — the opposite of *there is no reference class*, which is the truth.
    It is the mistake `no_homology` was retired for, in a new place.
    """
    label = session.execute(
        text("""
            select l.node_label
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = 'kp-nuna4' and l.prevalence_band = 'RARE'
               and l.cog_distinct_id_count > 0
             order by l.node_label limit 1
        """)
    ).scalar_one()
    payload = client.get(f"/api/v1/species/kp/loci/{label}/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
    band = by_kind["cog_orthogroup"]["own"]["propagation"]

    assert band["prevalence_band"] == "rare"
    assert band["checkable_node_count"] == 0, "a singleton node has one gene to check"
    assert band["rate"] is None, "⛔ not measurable"
    assert band["rate"] != 0.0, "⛔ and certainly not measured as zero"
    assert band["interval_low"] is None and band["interval_high"] is None
    assert band["one_gene_node_count"] > 0, "non-vacuity: the band this applies to is populated"


def test_the_DONOR_card_quotes_the_DONORS_band_and_not_the_recipients(client, session):
    """⛔ A transferred call propagates inside the DONOR's node, so it is the donor's band that
    governs — and the two bands differ often enough that taking the recipient's would be wrong
    rather than merely imprecise.

    Chosen by a query that says what it exercises, like the walk test below, rather than by a pinned
    label that could quietly stop exercising it.
    """
    row = session.execute(
        text("""
            select focal.node_label as focal, lower(focal.prevalence_band::text) as focal_band,
                   donor.node_label as donor, lower(donor.prevalence_band::text) as donor_band
              from locus focal
              join pangenome p on p.pangenome_id = focal.pangenome_id
              join locus_nearest_locus n
                   on n.locus_id = focal.locus_id and n.representation = 'ESM'
              join locus donor on donor.locus_id = n.neighbour_locus_id
             where p.catalogue_key = 'kp-nuna4'
               and focal.cog_distinct_id_count = 0
               and donor.cog_distinct_id_count > 0
               and donor.prevalence_band <> focal.prevalence_band
               and not exists (
                     select 1 from locus_nearest_locus earlier
                       join locus e on e.locus_id = earlier.neighbour_locus_id
                      where earlier.locus_id = focal.locus_id
                        and earlier.representation = 'ESM'
                        and earlier.rank < n.rank
                        and e.cog_distinct_id_count > 0)
             order by focal.node_label limit 1
        """)
    ).first()
    assert row is not None, "the catalogue has focal/donor band mismatches, or this test is vacuous"

    payload = client.get(f"/api/v1/species/kp/loci/{row.focal}/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}
    candidate = by_kind["cog_orthogroup"]["candidate"]
    assert candidate["donor"]["node_label"] == row.donor
    assert candidate["donor"]["propagation"]["prevalence_band"] == row.donor_band
    assert row.donor_band != row.focal_band, "non-vacuity: the two bands really do differ here"


def test_the_propagation_cache_is_per_catalogue_like_the_ladder_above(session):
    """⛔ The trap `clear_calibration_cache` exists for, in a second cache: one catalogue's base rate
    answering for another's would be invisible on the page and wrong in both.
    """
    clear_propagation_cache()
    kinds = {"pangenome_id": 2, "annotation_kind": AnnotationKind.COG_ORTHOGROUP}
    kp = propagation_for(session, prevalence_band=PrevalenceBand.CORE, **kinds)
    ecoli = propagation_for(
        session, pangenome_id=1, annotation_kind=AnnotationKind.COG_ORTHOGROUP,
        prevalence_band=PrevalenceBand.CORE,
    )
    assert kp.checkable != ecoli.checkable, "two catalogues, two measurements"
    assert kp.checkable >= MIN_NODES and ecoli.checkable >= MIN_NODES


def test_a_band_the_catalogue_never_measured_arrives_as_NOT_MEASURABLE_not_MISSING():
    """⛔ Missing key and not-measurable must reach the page as the SAME shape.

    ⚠ Seeds the private cache deliberately: the published catalogues populate all five bands, so the
    absent-band branch has no live case and would otherwise go untested — and it is the branch that
    decides whether the page has to tell *no such band* from *too few nodes*. With the cache seeded
    the session is never touched, which is why `None` is a safe argument here.
    """
    clear_propagation_cache()
    _PROPAGATION[(-1, AnnotationKind.KEGG_ORTHOLOGY)] = {}
    cell = propagation_for(
        None, pangenome_id=-1, annotation_kind=AnnotationKind.KEGG_ORTHOLOGY,
        prevalence_band=PrevalenceBand.SOFT_CORE,
    )
    clear_propagation_cache()
    assert cell.prevalence_band == "SOFT_CORE", "the dataclass keeps the Postgres NAME"
    assert cell.as_json()["prevalence_band"] == "soft_core", "and the wire carries the value"
    assert (cell.checkable, cell.unanimous, cell.one_gene) == (0, 0, 0)
    assert cell.rate is None and cell.interval is None
    assert cell.as_json()["rate"] is None


def test_a_node_with_its_own_call_is_never_offered_a_neighbours(client):
    """⛔ A suggestion beside a real annotation would compete with it."""
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    for entry in payload["inference"]["vocabularies"]:
        if entry["own"] is not None:
            assert entry["candidate"] is None and entry["walk"] == []


def test_the_walk_shows_the_ranks_that_carried_NOTHING(client, session):
    """⚠ For ecoli COG only 1,679 of 4,143 assignments come from rank 1, so this is the common case.

    The node is chosen by a query that says what it exercises — the FIRST annotated neighbour sits
    past rank 1 — rather than by a pinned label that could quietly stop exercising it.
    """
    label = session.execute(
        text("""
            select focal.node_label
              from locus focal
              join pangenome p on p.pangenome_id = focal.pangenome_id
              join locus_nearest_locus n
                   on n.locus_id = focal.locus_id and n.representation = 'ESM'
              join locus donor on donor.locus_id = n.neighbour_locus_id
             where p.catalogue_key = 'ecoli-nuna4'
               and focal.cog_distinct_id_count = 0
               and donor.cog_distinct_id_count > 0
               and n.cross_similarity >= 0.98
               and not exists (
                     select 1 from locus_nearest_locus earlier
                       join locus e on e.locus_id = earlier.neighbour_locus_id
                      where earlier.locus_id = focal.locus_id
                        and earlier.representation = 'ESM'
                        and earlier.rank < n.rank
                        and e.cog_distinct_id_count > 0)
               and n.rank >= 3
             limit 1
        """)
    ).scalar()
    assert label is not None, "no ecoli node has its first COG donor at rank 3+ — query is vacuous"

    entry = next(
        row for row in client.get(f"/api/v1/species/ecoli/loci/{label}/function").get_json()
        ["inference"]["vocabularies"] if row["annotation_kind"] == "cog_orthogroup"
    )
    candidate = entry["candidate"]
    assert candidate is not None and candidate["rank"] >= 3
    assert [step["rank"] for step in entry["walk"]] == sorted(
        step["rank"] for step in entry["walk"]
    ), "the walk is reported in rank order, outward"
    earlier = [step for step in entry["walk"] if step["rank"] < candidate["rank"]]
    assert earlier and not any(step["carries_annotation"] for step in earlier), (
        "every rank before the donor must be shown, and shown as empty"
    )


def test_every_candidate_carries_the_measured_rate_for_the_depth_it_quotes(client):
    """⭐ The rate is the thing the reader judges the suggestion by, so it travels WITH it."""
    payload = client.get("/api/v1/species/ecoli/loci/17173/function").get_json()
    cells = {(cell["tier"], cell["level"]): cell
             for cell in payload["calibration"]["cog_orthogroup"]["cells"]}
    assert cells, "the response ships no calibration at all"
    for entry in payload["inference"]["vocabularies"]:
        candidate = entry["candidate"]
        if candidate is None or candidate["level"] is None:
            continue
        quoted = candidate["calibration"]
        assert quoted is not None
        assert quoted["pairs"] >= MIN_PAIRS
        assert quoted["tier"] == candidate["tier"] and quoted["level"] == candidate["level"]
        if entry["annotation_kind"] == "cog_orthogroup":
            assert quoted == cells[(candidate["tier"], candidate["level"])], (
                "the rate beside the suggestion is the SAME cell the shipped ladder holds"
            )
        assert candidate["value"], "a quoted level must name what it asserts"
