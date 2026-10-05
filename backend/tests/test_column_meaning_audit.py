"""The adversarial pass, as a standing test — a column whose meaning and contents differ.

⛔ **This is not "find a bug".** Every defect this pass exists to catch passed every test in the
suite, because every test compared a value to a value: `cog_cat` held `CP` in a column that fits it,
a Pfam architecture was truncated at 256 characters into a valid-looking different architecture, and
`uniref50_impurity` was NULL on all 17,531 loci — which is *exactly* its documented meaning
(*"NULL below 5 labelled members"*) and was in fact "never loaded".

⭐ **So the instrument measures contents and puts them beside the declaration**, and the tests below
pin the answers the pass established. Each one is a question the census raised and a verification
that settled it — never a tolerance, and never a snapshot that could be re-baselined.

⚠ These run against `BACATLAS_DATABASE_URL` with both catalogues loaded. A census over an empty
table says nothing, and says so.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

import bacatlas_backend.models  # noqa: F401  (registers every table on the metadata)
from bacatlas_backend.database import Base
from bacatlas_backend.instruments.column_census import (
    WIDTH_HEADROOM_WARNING,
    ask_the_questions,
    run_column_census,
)
from bacatlas_backend.models.intergenic_gap import IntergenicGap
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import PangenomeEvaluation


@pytest.fixture(scope="module")
def loaded_session():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the audit runs against a loaded database")
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        loci = session.execute(select(func.count()).select_from(Locus)).scalar_one()
        if loci < 30_000:
            pytest.skip(f"only {loci:,} loci are loaded; the audit needs both published catalogues")
        yield session


@pytest.fixture(scope="module")
def census(loaded_session):
    observations = run_column_census(loaded_session)
    return {(row.table, row.column): row for row in observations}


# ── coverage, stated ───────────────────────────────────────────────────────────────────────────
def test_the_census_examines_every_column_and_NAMES_what_it_could_not(census, loaded_session):
    """⛔ A census that omitted its empty tables would report full coverage over half the schema.

    ⭐ **The oracle is the registered metadata, not a count.** This read `== 29` and had been red
    since the 31st table landed — so the one test whose job is to notice a new table had stopped
    noticing, and two tables were added through a permanently-failing assertion before a third and
    fourth (`enzyme_class`, `kegg_orthology`) arrived on 2026-10-05. A hand-maintained integer is
    the wrong instrument for "did the census cover everything": it goes stale silently, and the
    staleness looks exactly like the failure it exists to report. `Base.metadata` already knows
    every table, so the census is now compared against it by NAME, which says which table is
    missing instead of by how many.
    """
    observations = list(census.values())
    assert len(observations) >= 320, f"examined {len(observations)} columns"
    examined = {row.table for row in observations}
    # ⚠ `alembic_version` is Alembic's, not ours — it carries no meaning to audit and is not a
    # mapped model, so it is absent from the metadata and must be absent here too.
    registered = set(Base.metadata.tables) - {"alembic_version"}
    assert examined == registered, (
        f"not censused: {sorted(registered - examined)}; "
        f"censused but not registered: {sorted(examined - registered)}"
    )
    unexamined = sorted({row.table for row in observations if row.is_empty})
    # ⚠ One table is legitimately empty: the roster build report belongs to the Klebsiella lineage.
    assert unexamined == ["genome_collection_build_report"]


# ── the pass's own findings, pinned ────────────────────────────────────────────────────────────
def test_uniref50_impurity_is_LOADED_and_its_NULL_means_what_it_says(loaded_session):
    """⛔⛔ **The pass's sharpest find.** This column was 100 % NULL, and NULL is its DOCUMENTED
    value for a locus below five labelled members — so `WHERE uniref50_impurity IS NULL` returned
    every locus and read as a fact about the biology. It was never loaded at all."""
    measured, total = loaded_session.execute(
        select(func.count(Locus.uniref50_impurity), func.count()).select_from(Locus)
    ).one()
    assert 0 < measured < total, f"{measured:,} of {total:,} loci carry a measured impurity"
    # And where it is NULL, the documented reason must hold: fewer than five labelled members.
    contradictions = loaded_session.execute(
        select(func.count())
        .select_from(Locus)
        .where(Locus.uniref50_impurity.is_(None), Locus.uniref50_labelled_member_count >= 5)
    ).scalar_one()
    assert contradictions == 0, (
        f"{contradictions:,} loci have ≥5 labelled members and no impurity — NULL then means "
        "something other than what the column says it means"
    )


def test_the_other_cluster_table_columns_are_loaded_too(loaded_session):
    """Each was declared, documented, present in the source, and 100 % NULL.

    ⚠ `embed_within_over_nearest` left this list on 2026-09-24. It was the ratio of a locus's two
    MEDOID distances, so it went with them; `medoid_genome_id` and `medoid_flat_index` stayed,
    because they NAME the representative member rather than measuring anything from it.
    """
    filled = loaded_session.execute(
        select(
            func.count(Locus.uniref50_coverage),
            func.count(Locus.medoid_genome_id),
            func.count(Locus.medoid_flat_index),
        ).select_from(Locus)
    ).one()
    total = loaded_session.execute(select(func.count()).select_from(Locus)).scalar_one()
    assert all(value == total for value in filled), dict(zip(
        ("uniref50_coverage", "medoid_genome_id", "medoid_flat_index"),
        filled, strict=True,
    ))


def test_seqid_coverage_is_NULL_because_the_AUDIT_skipped_it_and_that_is_recorded(loaded_session):
    """⚠ The same value carrying a different fact. Before the ingest read the column, NULL meant
    *we never looked*; now it means *the audit ran with `--skip-seqid-to-medoid`* — and the
    pangenome's `omitted_sections` says so, which is what makes the two distinguishable."""
    from bacatlas_backend.models.pangenome import Pangenome

    measured = loaded_session.execute(
        select(func.count(Locus.seqid_coverage)).select_from(Locus)
    ).scalar_one()
    assert measured == 0
    omissions = loaded_session.execute(select(Pangenome.omitted_sections)).scalars().all()
    assert all("seqid_to_medoid" in (record or {}) for record in omissions), omissions


def test_an_audit_metric_name_is_never_truncated_into_a_collision(loaded_session, census):
    """⛔ `metric_name` is part of the UNIQUE key, so a truncation is not a lost suffix — it is two
    different metrics becoming one, and the loser refused with a constraint name for a reason."""
    observation = census[("pangenome_evaluation", "metric_name")]
    assert observation.declared_text_length >= 128
    assert observation.maximum_text_length < observation.declared_text_length
    longest = loaded_session.execute(
        select(PangenomeEvaluation.metric_name)
        .order_by(func.length(PangenomeEvaluation.metric_name).desc())
        .limit(1)
    ).scalar_one()
    assert longest == "family_split_across_clusters_n_split_genes_singleton"
    assert len(longest) == 52


def test_the_git_sha_recorded_NAMES_A_NUNA_COMMIT_and_not_an_ingest_tree_one(loaded_session):
    """⛔ A live bug for the length of one commit: `_git_sha()` runs in the cwd, so the ingest wrote
    BacAtlas's HEAD into a column that says *"the version of the registry this row was read from"*.
    A real short sha, the right length and shape, and about the wrong repository.

    ⚠ **The assertion is that the sha RESOLVES in nuna, not that it equals nuna's HEAD.** The stored
    value records which nuna the ingest read, so it is *supposed* to go stale the moment anyone
    commits to nuna — and it did, within an hour of being written. A test pinned to HEAD would fail
    on every unrelated commit, and a test that fails for a reason the reader does not care about is
    a test that gets deleted rather than fixed. Resolving the sha in each tree distinguishes the two
    repositories, which is the thing the bug got wrong, and is stable as both of them move.
    """
    import subprocess
    from pathlib import Path

    import nuna

    from bacatlas_backend.models.pangenome import Pangenome

    stored = {value for value in loaded_session.execute(select(Pangenome.git_sha)).scalars() if value}
    if not stored:
        pytest.skip("no git sha was recorded — nuna may not be a git checkout here")

    def resolves_in(tree: Path, sha: str) -> bool:
        try:
            return subprocess.run(
                ["git", "cat-file", "-e", f"{sha}^{{commit}}"], cwd=tree, capture_output=True
            ).returncode == 0
        except OSError:
            return False

    nuna_tree = Path(nuna.__file__).resolve().parent
    ingest_tree = Path(__file__).resolve().parents[1]
    for sha in stored:
        assert resolves_in(nuna_tree, sha), f"{sha} is not a commit in the nuna checkout"
        assert not resolves_in(ingest_tree, sha), (
            f"{sha} resolves in the INGEST tree — which is the bug this test exists for"
        )


# ── the invariants the census questions, answered ──────────────────────────────────────────────
def test_a_negative_gap_length_is_an_OVERLAP_and_the_signed_convention_is_intact(loaded_session):
    """⛔ 18.8 % of adjacent pairs overlap, at a median of −4 bases. Clamping at zero made *abuts
    exactly* and *overlaps by 190 bases* the same number."""
    overlapping, total = loaded_session.execute(
        select(
            func.count().filter(IntergenicGap.median_signed_length_nt < 0),
            func.count(),
        ).select_from(IntergenicGap)
    ).one()
    assert 0.1 < overlapping / total < 0.3, f"{overlapping:,} of {total:,} adjacencies overlap"
    deepest = loaded_session.execute(
        select(func.min(IntergenicGap.median_signed_length_nt))
    ).scalar_one()
    assert deepest < -100, "no deep overlap survives — a clamp may have been reintroduced"


def test_a_gap_variance_of_zero_is_a_MEASURED_zero_and_the_majority_case(loaded_session):
    """⛔ White on the track has to mean *identical in every genome*, never *small*."""
    zeros, total = loaded_session.execute(
        select(
            func.count().filter(IntergenicGap.length_variance_score == 0.0),
            func.count(),
        ).select_from(IntergenicGap)
    ).one()
    assert 0.8 < zeros / total < 0.9, f"{zeros:,} of {total:,} gaps vary not at all"
    nulls = loaded_session.execute(
        select(func.count()).select_from(IntergenicGap).where(
            IntergenicGap.length_variance_score.is_(None)
        )
    ).scalar_one()
    # ⚠ Read densely from `gap_table`, so NULL would mean *this run did not measure variance* — a
    # state that does not arise here. The payload's sparsity was a JSON encoding artifact.
    assert nulls == 0


def test_interest_score_negatives_are_a_sentinel_AND_a_real_term_and_the_two_differ(loaded_session):
    """⚠ *"Negative means barred"* is wrong: `(major − 1) × 2e6` makes loci negative without barring."""
    sentinel = loaded_session.execute(
        select(func.count()).select_from(Locus).where(Locus.interest_score == -1e12)
    ).scalar_one()
    negative = loaded_session.execute(
        select(func.count()).select_from(Locus).where(Locus.interest_score < 0)
    ).scalar_one()
    assert sentinel > 0
    assert negative > sentinel, "every negative score is the sentinel — the second term is untested"


# ── the instrument itself ──────────────────────────────────────────────────────────────────────
def test_a_fixed_width_value_that_fills_its_column_is_NOT_flagged_as_near_overflow():
    """⚠ A sha256 is 64 characters by construction. Flagging it drowns the real width findings."""
    from bacatlas_backend.instruments.column_census import ColumnObservation

    fixed = ColumnObservation(
        table="t", column="sha", declared_type="VARCHAR(64)", nullable=True, row_count=10,
        declared_text_length=64, minimum_text_length=64, maximum_text_length=64,
    )
    varying = ColumnObservation(
        table="t", column="name", declared_type="VARCHAR(64)", nullable=True, row_count=10,
        declared_text_length=64, minimum_text_length=4, maximum_text_length=63,
    )
    assert not [q for q in ask_the_questions(fixed) if q.startswith("WIDTH")]
    assert [q for q in ask_the_questions(varying) if q.startswith("WIDTH")]
    assert 0 < WIDTH_HEADROOM_WARNING < 1


def test_a_column_in_an_EMPTY_table_reports_that_it_was_not_examined():
    """⛔ *Not looked at* and *no problem found* must never be the same output."""
    from bacatlas_backend.instruments.column_census import ColumnObservation

    empty = ColumnObservation(table="t", column="c", declared_type="TEXT", nullable=True)
    assert empty.is_empty
    assert ask_the_questions(empty) == []
    assert "NOT EXAMINED" in empty.render()
