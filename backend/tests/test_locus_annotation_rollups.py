"""The locus annotation rollups against the gene rows they replace — ⭐ the free oracle of this stage.

⛔ **This file exists because of WHEN it can run.** `gene_functional_annotation` and
`gene_locus_membership` are one row per gene and are being deleted; the rollups on `locus` and in
`pangenome_propagation_rate` are what replace them. For exactly one stage both representations are
live, and that is the only window in which the new columns can be checked against the measurement
they are replacing rather than against themselves. Afterwards the gene-grain source is gone and the
only available check is that ingest agreed with itself.

**Three representations, and each pair catches what the other cannot:**

1. `locus.propagation_verdict_*` — written at ingest from nuna's annotation frame;
2. `pangenome_propagation_rate` — tallied from (1) in the same pass, so it proves the WRITE PATH:
   a column projected to the wrong position, an enum coerced to its value instead of its name, a
   partial `COPY`;
3. `within_node_propagation.compute_propagation` — reading `gene_functional_annotation` at request
   time, so it proves the DECISION: two readers of one parquet reaching the same value.

⭐ And for COG a **fourth**, sharing no code with any of them: `count(distinct cog_id) = 1` in pure
SQL. COG is single-valued per gene, so its five-state verdict is something SQL can decide on its own
— the same stance as `test_within_node_propagation`'s oracle, at per-locus grain instead of per-band.

⚠ **Every comparison asserts its own coverage before it reports agreement.** A rollup test that
`continue`s past a NULL column reports "0 differ" while having compared nothing, and the columns are
nullable by design — NULL means *this ingest generation did not compute it*. So each test counts what
it examined and fails on zero.
"""

from __future__ import annotations

import os
from collections import defaultdict

import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.derive_locus_propagation import PROPAGATION_VOCABULARIES
from bacatlas_backend.instruments.annotation_transfer import (
    ANNOTATED_COUNT_PREDICATE,
    GENE_VALUE_COLUMN,
    SUPPORTED_KINDS,
)
from bacatlas_backend.instruments.within_node_propagation import (
    VERDICT_TALLY,
    compute_propagation,
    propagation_verdict,
)
from bacatlas_backend.models.enumerations import (
    GENE_ONTOLOGY_NAMESPACE_NAMES,
    AnnotationKind,
    PropagationVerdict,
)
from bacatlas_backend.services.reference_name_service import gene_ontology_fold

#: The two catalogues every per-locus comparison runs over. ⚠ Both, never one: the ecoli and kp
#: cohorts differ in COG coverage by a factor the rollups are sensitive to, and a bug in a
#: vocabulary one of them barely carries would hide in the other.
CATALOGUES = [1, 2]

#: ⛔ The suffix each vocabulary's `locus` column uses, from the ingest's own mapping rather than
#: re-spelled here. `AnnotationKind.COG_ORTHOGROUP` is not `cog` anywhere in the schema.
SUFFIX_OF_KIND = {kind: suffix for suffix, kind in PROPAGATION_VOCABULARIES.items()}


def verdict_of(value) -> PropagationVerdict | None:
    """A `propagation_verdict_*` cell from a raw `text()` query → the enum member, or None.

    ⛔ **A `text()` statement carries no type information, so psycopg hands back the Postgres enum
    LABEL as a `str`** — and the labels are the enum NAMES (`UNANIMOUS`), not its values
    (`unanimous`), because that is what `staging_table_loader.coerce_value` writes. Comparing a str
    to a `PropagationVerdict` is silently False for every row, which in a difference loop reads as
    *everything differs* and in a membership test as *nothing is present*. Both have happened here.
    """
    if value is None or isinstance(value, PropagationVerdict):
        return value
    return PropagationVerdict[str(value)]


def band_of(value) -> str:
    """A `prevalence_band` cell → its NAME, whichever shape the driver returned."""
    return value.name if hasattr(value, "name") else str(value)


@pytest.fixture(scope="module")
def session():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — these oracles need a loaded database")
    with Session(create_engine(url, future=True)) as open_session:
        yield open_session


@pytest.fixture(scope="module")
def fold(session):
    """⛔ The SAME fold the service runs with, and never `None`.

    `test_within_node_propagation`'s oracle calls `compute_propagation` with no fold, where
    `go_rung` pools all three GO facets into the first. That is a configuration production never
    uses, and on it a GO node can never abstain — so the identity it asserts holds there and does
    not hold here. Every comparison in this file is against the folded path.
    """
    loaded = gene_ontology_fold(session)
    assert loaded.slim_of, "the GO reference is not loaded — a folded comparison would prove nothing"
    return loaded


def _verdicts_from_the_gene_rows(session, pangenome_id, kind, fold) -> dict[int, PropagationVerdict]:
    """`{locus_id: verdict}` recomputed from `gene_functional_annotation`, for every locus.

    ⚠ Includes a locus with NO annotated gene as `NO_COVERAGE`, which the instrument's own query
    cannot produce — it filters the column to non-null, so a silent locus simply has no row. The
    stored column has one for every locus, so the comparison has to as well or it would skip the
    majority of the catalogue (12,538 of 17,531 ecoli loci on COG).
    """
    column = GENE_VALUE_COLUMN[kind]
    grouped: dict[int, list] = defaultdict(list)
    rows = session.execute(
        text(f"""
            select l.locus_id, f.{column} as value
              from locus l
              join gene_locus_membership m on m.locus_id = l.locus_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where l.pangenome_id = :pangenome_id and f.{column} is not null
        """),
        {"pangenome_id": pangenome_id},
    ).all()
    for locus_id, value in rows:
        grouped[locus_id].append(value)
    every_locus = [
        row[0]
        for row in session.execute(
            text("select locus_id from locus where pangenome_id = :pangenome_id"),
            {"pangenome_id": pangenome_id},
        ).all()
    ]
    return {
        locus_id: propagation_verdict(kind, grouped.get(locus_id, []), fold) for locus_id in every_locus
    }


# ═══════════════════════════════════════════════ 1 · the verdict, locus by locus, against the genes
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_the_stored_verdict_matches_the_one_recomputed_from_the_GENE_rows(session, fold, pangenome_id, kind):
    """⭐⭐ The comparison this whole stage is sequenced around, at per-locus grain.

    Band totals can agree while individual loci disagree in compensating directions — two loci
    swapping `unanimous` and `discordant` within a band is invisible to every tally. This is the test
    that cannot be fooled that way, and it is only runnable while the gene rows exist.
    """
    suffix = SUFFIX_OF_KIND[kind]
    stored = {
        row.locus_id: verdict_of(row.verdict)
        for row in session.execute(
            text(f"select locus_id, propagation_verdict_{suffix} as verdict from locus "
                 "where pangenome_id = :pangenome_id"),
            {"pangenome_id": pangenome_id},
        ).all()
    }
    recomputed = _verdicts_from_the_gene_rows(session, pangenome_id, kind, fold)

    # ⛔ Coverage FIRST: the column is nullable and a catalogue loaded before it existed holds NULL
    # for every locus, on which a difference loop would report perfect agreement.
    assert stored, f"pangenome {pangenome_id} has no locus rows at all"
    assert set(stored) == set(recomputed), "the two sides do not describe the same loci"
    written = sum(1 for verdict in stored.values() if verdict is not None)
    assert written == len(stored), (
        f"{len(stored) - written:,} of {len(stored):,} loci hold NULL for propagation_verdict_{suffix}"
        " — re-ingest this catalogue; a NULL column compares equal to nothing and reports no difference"
    )

    differences = [
        (locus_id, stored[locus_id], recomputed[locus_id])
        for locus_id in stored
        if stored[locus_id] is not recomputed[locus_id]
    ]
    assert not differences, (
        f"{len(differences):,} of {len(stored):,} loci disagree, e.g. {differences[:5]} — the ingest "
        "frame and the gene rows are two readers of one parquet and they have diverged"
    )
    # ⚠ Non-vacuity: at least two DIFFERENT verdicts must occur, or a column stuck on one value
    # would pass every assertion above.
    assert len({verdict for verdict in stored.values()}) >= 2, (
        f"every locus has the same propagation_verdict_{suffix} — the comparison proves nothing"
    )


# ═══════════════════════════════════════════════════ 2 · the tally, against the verdicts and the genes
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_the_stored_tally_is_the_GROUP_BY_of_the_stored_verdicts_and_matches_the_genes(
    session, fold, pangenome_id, kind
):
    """Three numbers per band from three places, and all three must agree.

    ⛔ `checkable` is **`unanimous + discordant`**, not *two or more annotated genes*. The two
    coincide for COG only because a COG id always states its rung, so `agreed` is never `None`
    there; for GO, 81 ecoli loci state no shared facet and belong in neither column. Deriving the
    tally as `annotated >= 2` would move them into the denominator and drop the published rate with
    nothing turning red — hence `VERDICT_TALLY` rather than a hand-written predicate.
    """
    suffix = SUFFIX_OF_KIND[kind]
    from_verdicts: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    rows = session.execute(
        text(f"select prevalence_band, propagation_verdict_{suffix} as verdict from locus "
             "where pangenome_id = :pangenome_id"),
        {"pangenome_id": pangenome_id},
    ).all()
    assert rows, f"pangenome {pangenome_id} has no locus rows"
    for band, raw in rows:
        verdict = verdict_of(raw)
        assert verdict is not None, f"a NULL propagation_verdict_{suffix} cannot be tallied"
        counts = from_verdicts[band_of(band)]
        for index, increment in enumerate(VERDICT_TALLY[verdict]):
            counts[index] += increment

    stored = {
        band_of(row.prevalence_band): (row.checkable, row.unanimous, row.one_gene)
        for row in session.execute(
            text("select prevalence_band, checkable, unanimous, one_gene "
                 "from pangenome_propagation_rate "
                 "where pangenome_id = :pangenome_id and annotation_kind = :kind"),
            {"pangenome_id": pangenome_id, "kind": kind.name},
        ).all()
    }
    measured = compute_propagation(session, pangenome_id=pangenome_id, annotation_kind=kind, fold=fold)

    assert stored, f"no pangenome_propagation_rate rows for pangenome {pangenome_id} {kind.value}"
    assert set(stored) == set(from_verdicts), (
        "the stored tally and the verdict columns cover different bands: "
        f"only stored {sorted(set(stored) - set(from_verdicts))}, "
        f"only derived {sorted(set(from_verdicts) - set(stored))}"
    )
    for band, counts in sorted(from_verdicts.items()):
        assert stored[band] == tuple(counts), f"{kind.value} {band}: stored {stored[band]} vs verdicts {counts}"
        from_genes = measured.get(band)
        from_genes = (
            (from_genes.checkable, from_genes.unanimous, from_genes.one_gene) if from_genes else (0, 0, 0)
        )
        assert stored[band] == from_genes, f"{kind.value} {band}: stored {stored[band]} vs genes {from_genes}"
    assert sum(sum(counts) for counts in from_verdicts.values()) > 0, (
        f"{kind.value} tallies to zero in every band — the comparison is vacuous"
    )


# ═══════════════════════════════════════ 3 · COG, against SQL that shares none of the folding
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
def test_the_COG_verdict_matches_a_five_state_recount_done_entirely_in_SQL(session, pangenome_id):
    """⭐ The independent algorithm. A frozen oracle cannot see a drift it shares.

    COG is single-valued per gene, so each of the five states reduces to a count SQL can take on its
    own: no annotated gene, exactly one, two or more agreeing, two or more disagreeing. ⚠ There is no
    `not_comparable` leg because a non-null `cog_id` always states its rung — which is itself worth
    asserting, since its appearance would mean `_fold_gene` had changed under COG.
    """
    oracle = {
        row.locus_id: row.verdict
        for row in session.execute(
            text("""
                with per as (
                  select l.locus_id,
                         count(f.cog_id) as annotated,
                         count(distinct f.cog_id) as distinct_calls
                    from locus l
                    join gene_locus_membership m on m.locus_id = l.locus_id
                    left join gene_functional_annotation f
                         on f.genome_id = m.genome_id and f.flat_index = m.flat_index
                   where l.pangenome_id = :pangenome_id
                   group by 1)
                select locus_id,
                       case when annotated = 0 then 'NO_COVERAGE'
                            when annotated = 1 then 'ONE_GENE'
                            when distinct_calls = 1 then 'UNANIMOUS'
                            else 'DISCORDANT' end as verdict
                  from per
            """),
            {"pangenome_id": pangenome_id},
        ).all()
    }
    stored = {
        row.locus_id: verdict_of(row.verdict)
        for row in session.execute(
            text("select locus_id, propagation_verdict_cog as verdict from locus "
                 "where pangenome_id = :pangenome_id"),
            {"pangenome_id": pangenome_id},
        ).all()
    }
    assert oracle and set(oracle) == set(stored), "the oracle and the column describe different loci"
    assert PropagationVerdict.NOT_COMPARABLE not in set(stored.values()), (
        "a COG locus reported `not_comparable`, which cannot happen while a non-null cog_id states "
        "its own rung — `_fold_gene`'s COG branch has changed and this oracle no longer applies"
    )
    differences = [
        (locus_id, stored[locus_id], oracle[locus_id])
        for locus_id in oracle
        if stored[locus_id] is None or stored[locus_id].name != oracle[locus_id]
    ]
    assert not differences, f"{len(differences):,} loci differ, e.g. {differences[:5]}"
    # ⚠ Non-vacuity, stated per state: all four must actually occur, or the case that does not is
    # untested and the test would survive its own deletion.
    seen = {verdict for verdict in oracle.values()}
    assert seen == {"NO_COVERAGE", "ONE_GENE", "UNANIMOUS", "DISCORDANT"}, f"only {sorted(seen)} occur"


# ═══════════════════════════════════════════════════ 4 · the annotated counts, against the gene rows
@pytest.mark.parametrize("pangenome_id", [1, 2, 4, 5])
def test_every_annotated_member_count_matches_the_gene_rows_locus_by_locus(session, pangenome_id):
    """⭐ Including `go_annotated_member_count_any`, which is why the per-request join can go.

    ⛔ The three GO namespace columns CANNOT be summed to get it, and this test measures the factor
    rather than asserting the warning: on ecoli their sum is **2.30x** the real count of genes
    carrying a GO term. `annotation_transfer.py:231-234` is the comment; this is the number.
    """
    row = session.execute(
        text("""
            with recount as (
              select m.locus_id,
                     count(*) filter (where f.cog_id is not null)             as cog,
                     count(*) filter (where f.ec_numbers is not null)          as ec,
                     count(*) filter (where f.kegg_orthology_id is not null)   as kegg,
                     count(*) filter (where f.gene_ontology_terms is not null) as go_any
                from gene_locus_membership m
                left join gene_functional_annotation f
                       on f.genome_id = m.genome_id and f.flat_index = m.flat_index
               where m.pangenome_id = :pangenome_id
               group by 1
            )
            select count(*) as loci,
                   count(*) filter (where l.cog_annotated_member_count  is distinct from coalesce(r.cog,0))  as cog_bad,
                   count(*) filter (where l.ec_annotated_member_count   is distinct from coalesce(r.ec,0))   as ec_bad,
                   count(*) filter (where l.kegg_annotated_member_count is distinct from coalesce(r.kegg,0))
                       as kegg_bad,
                   count(*) filter (where l.go_annotated_member_count_any is null) as go_any_null,
                   count(*) filter (where l.go_annotated_member_count_any is distinct from coalesce(r.go_any,0))
                       as go_any_bad,
                   sum(coalesce(r.go_any,0)) as go_any_genes,
                   sum(coalesce(l.go_annotated_member_count_molecular_function,0)
                      +coalesce(l.go_annotated_member_count_biological_process,0)
                      +coalesce(l.go_annotated_member_count_cellular_component,0)) as three_column_sum
              from locus l left join recount r on r.locus_id = l.locus_id
             where l.pangenome_id = :pangenome_id
        """),
        {"pangenome_id": pangenome_id},
    ).one()
    assert row.loci > 0, f"pangenome {pangenome_id} has no loci"
    assert row.go_any_null == 0, (
        f"{row.go_any_null:,} of {row.loci:,} loci hold NULL for go_annotated_member_count_any — "
        "re-ingest; a NULL is not a mismatch and would not be counted below"
    )
    assert (row.cog_bad, row.ec_bad, row.kegg_bad, row.go_any_bad) == (0, 0, 0, 0), (
        f"mismatching loci: cog {row.cog_bad}, ec {row.ec_bad}, kegg {row.kegg_bad}, go_any {row.go_any_bad}"
    )
    assert row.go_any_genes > 0, "no gene in this catalogue carries a GO term — the test is vacuous"
    # ⛔ The reason `go_annotated_member_count_any` exists, asserted as an inequality rather than
    # quoted: the three namespace columns genuinely over-count, by a wide margin.
    # ⚠ `sum()` comes back as `Decimal`, so the factor is written as integer arithmetic rather
    # than `* 1.5`, which raises on a Decimal and would have failed as a TypeError rather than a
    # measurement.
    assert 2 * int(row.three_column_sum) > 3 * int(row.go_any_genes), (
        f"the three GO namespace columns sum to {row.three_column_sum:,} against "
        f"{row.go_any_genes:,} genes carrying a term — if that factor has fallen to ~1, GO is no "
        "longer multi-namespace and the `_any` column's justification needs re-reading"
    )


# ═════════════════════════════════════════════════════ 5 · the claim counts, and what bounds them
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
def test_a_folded_GO_claim_is_never_more_common_than_a_GO_term(session, pangenome_id):
    """⛔ `go_claiming_member_count_*` <= `go_annotated_member_count_*`, per locus, per namespace.

    A gene asserting only the namespace root carries a term and states no claim, so the claim count
    is bounded above by the annotated count and the gap is the root-only population. ⚠ The gap is
    asserted NON-ZERO somewhere: equal columns would mean the fold is doing nothing, which is what a
    missing GO reference looks like, and that would pass an inequality on its own.
    """
    for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES:
        row = session.execute(
            text(f"""
                select count(*) as loci,
                       count(*) filter (where go_claiming_member_count_{namespace} is null) as unwritten,
                       count(*) filter (where go_claiming_member_count_{namespace}
                                            > go_annotated_member_count_{namespace}) as above,
                       sum(go_annotated_member_count_{namespace} - go_claiming_member_count_{namespace}) as gap
                  from locus where pangenome_id = :pangenome_id
            """),
            {"pangenome_id": pangenome_id},
        ).one()
        assert row.loci > 0 and row.unwritten == 0, (
            f"{row.unwritten:,} of {row.loci:,} loci hold NULL for "
            f"go_claiming_member_count_{namespace} — re-ingest, or the fold was absent"
        )
        assert row.above == 0, (
            f"{row.above:,} loci claim more {namespace} than they carry, which is not possible: a "
            "claim is folded FROM the terms"
        )
        assert row.gap > 0, (
            f"no {namespace} gap anywhere — every gene's terms survive the fold, which is what an "
            "unloaded GO reference looks like"
        )


# ════════════════════════════════════════════ 6 · the folded identity the unfolded one cannot state
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_the_populations_account_for_every_node_carrying_the_vocabulary(session, fold, pangenome_id, kind):
    """⭐ The coverage identity, repaired for the folded path — and `not_comparable` is the repair.

    `test_within_node_propagation` asserts `checkable + one_gene == nodes carrying`, which is exact
    only because it calls `compute_propagation` with NO fold: unfolded, GO pools its three facets and
    a node can never abstain. On the folded path production runs, **81 ecoli and 33 kp loci** carry
    GO on two or more genes with no facet every one of them states, and they are in neither column.
    The stored `not_comparable` count is exactly that shortfall, which makes the identity exact
    again — and is a use for the column that no tally provides.
    """
    suffix = SUFFIX_OF_KIND[kind]
    rates = compute_propagation(session, pangenome_id=pangenome_id, annotation_kind=kind, fold=fold)
    counted = sum(one.checkable + one.one_gene for one in rates.values())
    not_comparable = session.execute(
        text(f"select count(*) from locus where pangenome_id = :pangenome_id "
             f"and propagation_verdict_{suffix} = 'NOT_COMPARABLE'"),
        {"pangenome_id": pangenome_id},
    ).scalar_one()
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
    assert counted + not_comparable == carrying, (
        f"{kind.value}: checkable+one_gene {counted:,} plus not_comparable {not_comparable:,} is "
        f"{counted + not_comparable:,}, against {carrying:,} loci carrying the vocabulary"
    )


def test_GO_really_does_abstain_somewhere_so_the_identity_above_is_not_vacuous(session):
    """⛔ Without this, test 6 is `x + 0 == x` and would pass with the column always NULL or 0.

    The whole point of distinguishing `not_comparable` from `no_coverage` is that the first occurs.
    """
    abstaining = session.execute(
        text("select pangenome_id, count(*) from locus "
             "where propagation_verdict_go = 'NOT_COMPARABLE' group by 1 order by 1")
    ).all()
    assert abstaining, (
        "no locus in any catalogue reports `not_comparable` for GO — either the column is unwritten "
        "or the fold is absent, and the five-state vocabulary is then four states with a spare"
    )
    assert all(count > 0 for _pangenome_id, count in abstaining)


# ══════════════════════════════════════════════════════ 7 · the fold a stored verdict is read against
def test_every_catalogue_records_the_fold_its_GO_columns_were_computed_under(session, fold):
    """⛔ A stored GO verdict with no recorded fold cannot be known to mean what a reader assumes.

    The verdict is a comparison of sets closed upward through the slim; change the slim and the
    verdict moves. Before it was stored, `compute_propagation` re-folded on every cold process, so a
    vendored-table update propagated with no re-ingest. Stored, it does not — so the version is
    written beside it and this is where a drift surfaces.
    """
    rows = session.execute(
        text("select pangenome_id, annotation_fold_version from pangenome order by 1")
    ).all()
    assert rows, "no catalogues at all"
    live = fold.version()
    stale = [(pangenome_id, version) for pangenome_id, version in rows if version != live]
    assert not stale, (
        f"the loaded fold is {live} and these catalogues were built under another: {stale}. Their "
        "propagation_verdict_go and go_claiming_member_count_* describe a different reference; "
        "re-ingest them or do not quote those columns"
    )


def test_the_fold_version_moves_when_the_fold_does_and_not_when_a_NAME_does(fold):
    """⚠ Three maps decide a verdict and the fourth is display only.

    Folding `name_of` into the digest would make a corrected GO *label* invalidate every stored
    verdict in the database and force a re-ingest that changes no number. ⛔ And an EMPTY fold must
    have its own version rather than none: with no reference `go_rung` pools the three facets, which
    is a different comparison, so *"no reference loaded"* has to be as detectable as a different one.
    """
    from dataclasses import replace

    from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold

    assert replace(fold, name_of={}).version() == fold.version(), "a display name moved the version"
    assert replace(fold, name_of={"GO:0000001": "renamed"}).version() == fold.version()

    moved = dict(fold.slim_of)
    moved["GO:0000001"] = frozenset({"GO:0016020"})
    assert replace(fold, slim_of=moved).version() != fold.version(), "a changed fold did NOT move it"
    assert GeneOntologyFold().version() != fold.version()
    assert GeneOntologyFold().version().startswith("go:0:0:"), "an empty fold must say it is empty"


# ═════════════════════════════════════════════════════════════════ 8 · the table's own shape
def test_the_rate_table_holds_one_row_per_vocabulary_and_band_with_no_rate_computed_on_too_few(session):
    """⚠ The floor lives in Python, so the stored numbers must not be read as a rate directly.

    `PropagationRate.rate` returns `None` below `MIN_NODES = 30` and the table has no `rate` column,
    so anything computing `unanimous / checkable` off these rows re-introduces the 0.0 that property
    refuses. This test pins what the table DOES guarantee — one row per band present, per
    vocabulary, `unanimous <= checkable` — and names the bands that sit under the floor so their
    presence is deliberate rather than noticed later.
    """
    rows = session.execute(
        text("""
            select p.pangenome_id, count(*) as rows,
                   count(*) filter (where r.unanimous > r.checkable) as impossible,
                   count(*) filter (where r.checkable < 30 and r.one_gene > 0) as under_the_floor
              from pangenome p join pangenome_propagation_rate r on r.pangenome_id = p.pangenome_id
             group by 1 order by 1
        """)
    ).all()
    assert rows, "no pangenome has any propagation rate rows"
    for pangenome_id, count, impossible, under_the_floor in rows:
        bands = session.execute(
            text("select count(distinct prevalence_band) from locus where pangenome_id = :pangenome_id"),
            {"pangenome_id": pangenome_id},
        ).scalar_one()
        assert count == len(PROPAGATION_VOCABULARIES) * bands, (
            f"pangenome {pangenome_id}: {count} rate rows for {bands} bands x "
            f"{len(PROPAGATION_VOCABULARIES)} vocabularies — a missing row reaches the page as a "
            "missing key where it expects 'not measurable'"
        )
        assert impossible == 0, f"pangenome {pangenome_id}: {impossible} rows claim unanimous > checkable"
        assert under_the_floor > 0, (
            f"pangenome {pangenome_id} has no band under MIN_NODES with a one-gene population — "
            "RARE has been that band in every catalogue measured, and its absence would mean the "
            "three-case split on the page is now describing something else"
        )


def test_the_four_vocabularies_are_the_ones_the_page_asks_for(session):
    """⛔ The rate table's `annotation_kind` must cover exactly `SUPPORTED_KINDS` and no more.

    `AnnotationKind` carries seven values; `gene_symbol`, `protein_product` and `pfam_architecture`
    have no rung ladder and no propagation question, so a row for one of them would be a measurement
    of nothing wearing a real column name.
    """
    stored = {
        row[0] if isinstance(row[0], AnnotationKind) else AnnotationKind[str(row[0])]
        for row in session.execute(text("select distinct annotation_kind from pangenome_propagation_rate")).all()
    }
    assert stored == set(SUPPORTED_KINDS), (
        f"stored kinds {sorted(k.value for k in stored)} against the page's "
        f"{sorted(k.value for k in SUPPORTED_KINDS)}"
    )
    assert set(PROPAGATION_VOCABULARIES.values()) == set(SUPPORTED_KINDS), (
        "the ingest writes a different set of vocabularies than the page reads"
    )
