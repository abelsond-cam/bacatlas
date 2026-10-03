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
    MIN_PAIRS,
    NOT_CALLED,
    QUOTED_LEVEL,
    TIERS,
    cog_levels,
    ec_levels,
    quotable_level,
    tier_for,
)
from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.services.function_inference_service import (
    calibration_for,
    clear_calibration_cache,
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
    """⛔ David, 2026-10-03 — nuna `PROJECT_STATE.md` §6. Not a tuning parameter.

    Three separate decisions live in these two constants, and each was argued from a measurement:
    below 0.96 is ONE tier because 0.95-0.96 and 0.90-0.95 were indistinguishable; nothing below
    0.90 appears at all ("definitely too remote to call"); and EC quotes its full code only in the
    top tier, dropping a level at 0.98-0.99 where the full code is 84 % and the sub-subclass 95 %.
    """
    assert [name for name, _, _ in TIERS] == [
        ">= 0.99", "0.98-0.99", "0.97-0.98", "0.96-0.97", "0.90-0.96",
    ]
    assert [low for _, low, _ in TIERS] == [0.99, 0.98, 0.97, 0.96, 0.90]
    assert tier_for(0.8999) == NOT_CALLED, "below 0.90 is never called"
    assert tier_for(None) == NOT_CALLED, "a missing cosine is not a tier"
    assert QUOTED_LEVEL[AnnotationKind.EC_NUMBER] == {
        ">= 0.99": 4, "0.98-0.99": 3, "0.97-0.98": 1, "0.96-0.97": 1, "0.90-0.96": 1,
    }
    assert QUOTED_LEVEL[AnnotationKind.COG_ORTHOGROUP] == {
        ">= 0.99": 2, "0.98-0.99": 2, "0.97-0.98": 1, "0.96-0.97": 1, "0.90-0.96": 1,
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
def test_gumC_says_its_COG_rests_on_ONE_gene_and_cannot_be_checked(client):
    """⭐ The worked case: kp node 722, 100 members, COG3206 on **1** of them, EC on 16.

    ⛔ The two vocabularies must not read alike on this card. One rests on a vote with sixteen
    voters; the other on a vote with one, which is unanimous by construction.
    """
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    by_kind = {entry["annotation_kind"]: entry for entry in payload["inference"]["vocabularies"]}

    cog = by_kind["cog_orthogroup"]["own"]
    assert cog["term"] == "COG3206"
    assert (cog["gene_count"], cog["annotated_gene_count"], cog["member_gene_count"]) == (1, 1, 100)
    assert cog["checkable"] is False

    ec = by_kind["ec_number"]["own"]
    assert ec["term"] == "2.7.10.-"
    assert ec["annotated_gene_count"] == 16 and ec["member_gene_count"] == 100
    assert ec["checkable"] is True
    assert ec_levels(ec["term"])[4] == frozenset(), (
        "and its EC states only three levels however similar a reader's neighbour is"
    )


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
