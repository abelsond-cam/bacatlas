"""The locus annotation rollups — the stored columns against each other and against their bounds.

⛔ **The gene-row halves are gone.** While `gene_functional_annotation` and `gene_locus_membership`
existed, this file checked every stored verdict and count against a recount over them; those tables
were dropped in D8 and that oracle with them. What remains checks what can still be checked:

1. `pangenome_propagation_rate` is the `GROUP BY` of `locus.propagation_verdict_*` — the WRITE PATH:
   a column projected to the wrong position, an enum coerced to its value instead of its name, a
   partial `COPY`;
2. `checkable + one_gene + not_comparable` equals the loci whose stored annotated count is non-zero
   — the verdicts against the counts, two ingest outputs computed by different code;
3. the GO claim counts are bounded by the term counts, `_any` is written, and the fold is recorded.

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
    LOCUS_ANNOTATED_COUNT_COLUMN,
    SUPPORTED_KINDS,
)
from bacatlas_backend.instruments.within_node_propagation import VERDICT_TALLY
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
    """⛔ The SAME fold the service runs with, and never `None` — the stored GO columns are folded."""
    loaded = gene_ontology_fold(session)
    assert loaded.slim_of, "the GO reference is not loaded — a folded comparison would prove nothing"
    return loaded


# ═══════════════════════════════════════════════════════════ 1 · the tally, against the verdicts
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_the_stored_tally_is_the_GROUP_BY_of_the_stored_verdicts(session, pangenome_id, kind):
    """Three numbers per band, stored and re-derived from the stored verdicts, and they must agree.

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
    assert stored, f"no pangenome_propagation_rate rows for pangenome {pangenome_id} {kind.value}"
    assert set(stored) == set(from_verdicts), (
        "the stored tally and the verdict columns cover different bands: "
        f"only stored {sorted(set(stored) - set(from_verdicts))}, "
        f"only derived {sorted(set(from_verdicts) - set(stored))}"
    )
    for band, counts in sorted(from_verdicts.items()):
        assert stored[band] == tuple(counts), f"{kind.value} {band}: stored {stored[band]} vs verdicts {counts}"
    assert sum(sum(counts) for counts in from_verdicts.values()) > 0, (
        f"{kind.value} tallies to zero in every band — the comparison is vacuous"
    )


# ════════════════════════════════════════════ 2 · `go_annotated_member_count_any`, and why it exists
@pytest.mark.parametrize("pangenome_id", [1, 2, 4, 5])
def test_go_any_is_written_and_the_three_namespace_columns_overcount_it(session, pangenome_id):
    """⛔ `function_inference_service._annotated_counts` substitutes 0 for a NULL here, and this is the
    gate that makes that safe: no catalogue may hold a NULL `go_annotated_member_count_any`.

    The three GO namespace columns CANNOT be summed to get it, and the factor is asserted rather than
    quoted (ecoli measured **2.30x**; `annotation_transfer.LOCUS_ANNOTATED_COUNT_COLUMN`).
    """
    row = session.execute(
        text("""
            select count(*) as loci,
                   count(*) filter (where go_annotated_member_count_any is null) as go_any_null,
                   sum(coalesce(go_annotated_member_count_any, 0)) as go_any_genes,
                   sum(coalesce(go_annotated_member_count_molecular_function,0)
                      +coalesce(go_annotated_member_count_biological_process,0)
                      +coalesce(go_annotated_member_count_cellular_component,0)) as three_column_sum
              from locus where pangenome_id = :pangenome_id
        """),
        {"pangenome_id": pangenome_id},
    ).one()
    assert row.loci > 0, f"pangenome {pangenome_id} has no loci"
    assert row.go_any_null == 0, (
        f"{row.go_any_null:,} of {row.loci:,} loci hold NULL for go_annotated_member_count_any — re-ingest"
    )
    assert row.go_any_genes > 0, "no gene in this catalogue carries a GO term — the test is vacuous"
    # ⚠ `sum()` comes back as `Decimal`, so the factor is integer arithmetic rather than `* 1.5`.
    assert 2 * int(row.three_column_sum) > 3 * int(row.go_any_genes), (
        f"the three GO namespace columns sum to {row.three_column_sum:,} against "
        f"{row.go_any_genes:,} genes carrying a term — if that factor has fallen to ~1, GO is no "
        "longer multi-namespace and the `_any` column's justification needs re-reading"
    )


# ═════════════════════════════════════════════════════ 3 · the claim counts, and what bounds them
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


# ══════════════════════════════════════════ 4 · the verdicts against the annotated counts, per locus
@pytest.mark.parametrize("pangenome_id", CATALOGUES)
@pytest.mark.parametrize("kind", SUPPORTED_KINDS, ids=lambda kind: kind.value)
def test_the_populations_account_for_every_node_carrying_the_vocabulary(session, pangenome_id, kind):
    """⭐ `checkable + one_gene + not_comparable == loci carrying the vocabulary`, on the folded path.

    The verdicts and the annotated counts are computed by different code from the same frame, so a
    locus whose verdict says *no_coverage* while its count says it carries the vocabulary (or the
    reverse) breaks the identity. On the folded path **81 ecoli and 33 kp loci** carry GO on two or
    more genes with no facet every one of them states — they are `not_comparable`, which is why that
    term is in the sum.
    """
    suffix = SUFFIX_OF_KIND[kind]
    tallied = session.execute(
        text("select coalesce(sum(checkable + one_gene), 0) from pangenome_propagation_rate "
             "where pangenome_id = :pangenome_id and annotation_kind = :kind"),
        {"pangenome_id": pangenome_id, "kind": kind.name},
    ).scalar_one()
    not_comparable = session.execute(
        text(f"select count(*) from locus where pangenome_id = :pangenome_id "
             f"and propagation_verdict_{suffix} = 'NOT_COMPARABLE'"),
        {"pangenome_id": pangenome_id},
    ).scalar_one()
    carrying = session.execute(
        text(f"select count(*) from locus where pangenome_id = :pangenome_id "
             f"and {LOCUS_ANNOTATED_COUNT_COLUMN[kind]} > 0"),
        {"pangenome_id": pangenome_id},
    ).scalar_one()
    assert carrying > 0, "non-vacuity: this catalogue does carry this vocabulary"
    assert tallied + not_comparable == carrying, (
        f"{kind.value}: checkable+one_gene {tallied:,} plus not_comparable {not_comparable:,} is "
        f"{tallied + not_comparable:,}, against {carrying:,} loci whose "
        f"{LOCUS_ANNOTATED_COUNT_COLUMN[kind]} is non-zero"
    )


def test_GO_really_does_abstain_somewhere_so_the_identity_above_is_not_vacuous(session):
    """⛔ Without this, test 4 is `x + 0 == x` and would pass with the column always NULL or 0.

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


# ══════════════════════════════════════════════════════ 5 · the fold a stored verdict is read against
def test_every_catalogue_records_the_fold_its_GO_columns_were_computed_under(session, fold):
    """⛔ A stored GO verdict with no recorded fold cannot be known to mean what a reader assumes.

    The verdict is a comparison of sets closed upward through the slim; change the slim and the
    verdict moves. Before it was stored, the request-time instrument re-folded on every cold process, so a
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


# ═════════════════════════════════════════════════════════════════ 6 · the table's own shape
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
