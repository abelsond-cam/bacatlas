"""What each HTTP ROUTE costs — measured at the route, not at the service behind it.

⛔⛔ **Why this file exists.** Every cost test before it measured a SERVICE (`load_locus_detail`,
`search_loci`, `load_gene_sequences`), and every one passed, while every route that is not the shell
was resolving its species through `load_species_catalogue` — the whole shell, including a GROUP BY
over every locus of the species — to learn one integer. A locus click cost ~17 statements against a
budget of 9, and search paid a whole-catalogue aggregation on every keystroke. Nothing could see it,
because the thing measured was not the thing served.

So each budget here is the service's own logical minimum PLUS ONE statement to resolve the species,
derived rather than recorded, and one property is asserted for every route at once: **only the shell
may aggregate over the catalogue.**
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle
from bacatlas_backend.services.function_inference_service import (
    clear_calibration_cache,
    clear_propagation_cache,
)

#: Resolving `species_key` → its published pangenome: one statement, on every route.
SPECIES_RESOLUTION = 1

#: `test_api_endpoints.test_a_locus_view_issues_ONE_statement_PER_TABLE_and_no_more` derives this:
#: NINE tables plus one neighbour-resolution statement. ⚠ It was 9 until 2026-09-24, when the medoid
#: geometry's single table became two — the per-locus similarities and the ranked nearest loci, which
#: have different cardinalities and so are deliberately not joined into one statement.
LOCUS_VIEW_MINIMUM = 10

#: The Function tab on a WARM process: the locus view's own statements plus the species resolution,
#: the function block and the inference walk. Measured on both published catalogues.
FUNCTION_TAB_WARM = 14

#: ⭐ The Function tab on a COLD process — **measured 2026-10-26: ecoli 28, kp 26**, and the larger
#: is the budget. This number is new because it is now worth bounding: every one of these statements
#: is a bounded read. The cold tab used to issue four catalogue-wide scans of
#: `gene_functional_annotation` on top, at 0.22-1.94 s each, and a statement budget said nothing
#: useful about that — the cost was in the rows, not the count.
FUNCTION_TAB_COLD = 28

#: ⚠ A genome PROJECTED onto the ecoli catalogue, if one is loaded. The projection is an optional
#: addition to a published catalogue, so these two tests skip rather than fail where none exists —
#: a database without a projection is a valid database.
PROJECTED_SAMPLE_ID = os.environ.get("BACATLAS_PROJECTED_SAMPLE_ID", "SAMN05374479")


@pytest.fixture(scope="module")
def application():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — route costs are measured on a loaded database")
    create_engine(url, future=True)
    return create_application(
        Configuration(
            database_url=url,
            artifact_roots={"gff": __import__("pathlib").Path(os.environ["BACATLAS_ROOT_GFF"])}
            if os.environ.get("BACATLAS_ROOT_GFF")
            else {},
        )
    )


def _measure(application, path, **query):
    client = application.test_client()
    with SqlCostOracle(application.extensions["bacatlas_database"].engine) as report:
        response = client.get(path, query_string=query)
    assert response.status_code == 200, (path, response.status_code, response.get_data(as_text=True)[:300])
    return report


def _catalogue_aggregations(report):
    """Statements that aggregate over `locus` — the census's signature."""
    return [
        record.summary()
        for record in report.statements
        if "GROUP BY" in record.sql.upper() and "FROM locus" in record.sql
    ]


def test_a_locus_CLICK_costs_the_locus_view_plus_one(application):
    report = _measure(application, "/api/v1/species/ecoli/loci/2811")
    report.assert_at_most(LOCUS_VIEW_MINIMUM + SPECIES_RESOLUTION, what="one locus click, at the route")


def test_an_ANCHORED_locus_click_adds_exactly_the_genome_lookup(application):
    report = _measure(application, "/api/v1/species/ecoli/loci/2811", anchor="SAMEA103923484")
    report.assert_at_most(
        LOCUS_VIEW_MINIMUM + SPECIES_RESOLUTION + 1, what="one anchored locus click, at the route"
    )


def test_a_PROJECTED_locus_click_adds_the_genome_lookup_and_its_placements(application):
    """⚠ Two statements more than an unanchored click, and they are DIFFERENT statements from the
    anchor's: a projected genome is in no `member_genome_ids`, so it is resolved from
    `projected_genome` and its genes at this locus are read from `projected_gene_placement`. Both
    are index scans on `(pangenome_id, …)` — the point of the budget is that neither grows with the
    catalogue."""
    report = _measure(
        application, "/api/v1/species/ecoli/loci/2811", projected=PROJECTED_SAMPLE_ID
    )
    report.assert_at_most(
        LOCUS_VIEW_MINIMUM + SPECIES_RESOLUTION + 2, what="one projected locus click, at the route"
    )


def test_the_projected_genome_LIST_is_two_statements(application):
    """One to resolve the species, one to list. It is O(projected genomes), never O(catalogue)."""
    report = _measure(application, "/api/v1/species/ecoli/projected-genomes")
    report.assert_at_most(1 + SPECIES_RESOLUTION, what="the projected-genome list, at the route")


def test_a_search_KEYSTROKE_is_two_statements(application):
    report = _measure(application, "/api/v1/species/ecoli/search", q="ligase")
    report.assert_at_most(1 + SPECIES_RESOLUTION, what="one search keystroke, at the route")


def test_the_residual_lists_are_two_statements(application):
    report = _measure(application, "/api/v1/species/kp/audit/residual-loci")
    report.assert_at_most(1 + SPECIES_RESOLUTION, what="the footer's residual lists, at the route")


@pytest.mark.parametrize(
    "query",
    [{"limit": 1}, {"limit": 1000}, {"q": "samea22", "limit": 1000}, {"q": "NOT-A-GENOME"}],
)
def test_an_anchor_picker_KEYSTROKE_is_two_statements_whatever_the_limit(application, query):
    """⛔ The per-genome counts are READ. Aggregated here they would be ~412 M rows per keystroke.

    One statement for the genomes — their matched count rides on every row via a window, so it does
    not grow with the limit, the collection or the query — plus the species resolution.
    """
    report = _measure(application, "/api/v1/species/ecoli/genomes", **query)
    report.assert_at_most(1 + SPECIES_RESOLUTION, what="one anchor-picker keystroke, at the route")
    # ⛔ And no aggregation over a membership table at all: the counts come from their own table.
    assert not [
        record.summary()
        for record in report.statements
        if "GROUP BY" in record.sql.upper()
        or "gene_locus_membership" in record.sql
        or "unnest" in record.sql.lower()
    ]


@pytest.mark.parametrize(
    "path",
    [
        "/api/v1/species/ecoli/loci/2811",
        "/api/v1/species/ecoli/loci/2811/function",
        "/api/v1/species/ecoli/loci/2811/arrangements",
        "/api/v1/species/ecoli/search?q=ligase",
        "/api/v1/species/kp/audit/residual-loci",
        "/api/v1/species/ecoli/genomes?q=samea22",
        "/api/v1/species/ecoli/projected-genomes",
    ],
)
def test_ONLY_the_shell_aggregates_over_the_catalogue(application, path):
    """⛔ The property the whole file exists for, over every per-locus and per-keystroke route."""
    route, _, query = path.partition("?")
    report = _measure(application, route, **dict(pair.split("=") for pair in query.split("&") if pair))
    assert _catalogue_aggregations(report) == []


def test_and_the_shell_DOES_aggregate_so_the_property_above_is_not_vacuous(application):
    report = _measure(application, "/api/v1/species/ecoli")
    assert len(_catalogue_aggregations(report)) == 1


def _gene_grain_catalogue_scans(report):
    """Statements reading `gene_functional_annotation` across a whole catalogue.

    ⚠ It still excludes `locus_id = any(:locus_ids)`, which used to be `_annotated_counts` reading
    six loci rather than a catalogue. That reader has since moved to `locus.*_annotated_member_count`
    and issues no such statement at all, so the exclusion now matches nothing on these routes — it is
    kept because the detector is about SCANS and a future bounded reader of the gene rows should not
    trip it.
    """
    return [
        record.summary()
        for record in report.statements
        if "gene_functional_annotation" in record.sql and "locus_id = any(" not in record.sql
    ]


def _gene_grain_reads(report):
    """ANY statement touching either one-row-per-gene table, bounded or not."""
    return [
        record.summary()
        for record in report.statements
        if "gene_functional_annotation" in record.sql or "gene_locus_membership" in record.sql
    ]


@pytest.mark.parametrize("species,label", [("ecoli", "2811"), ("kp", "722")])
def test_the_FUNCTION_tab_pays_its_catalogue_SCANS_once_per_process_and_never_again(
    application, species, label
):
    """⭐⭐ **The Function tab now reads NO gene row at all — cold or warm — and that is the change.**

    This test used to assert the opposite half and say so: *"the guarantee is therefore not that the
    scan is cheap but that it happens ONCE"*. The scan was `within_node_propagation`'s propagation
    base rate over every annotated gene of the catalogue, measured at **0.22-1.94 s** per (catalogue,
    vocabulary), so a cold tab paid 4.2-5.2 s wall. It is now a primary-key read of at most five rows
    from `pangenome_propagation_rate`, and `_annotated_counts` reads four stored `locus` columns
    instead of joining six loci to their gene rows. ⛔ So the old non-vacuity assertion —
    `assert _gene_grain_catalogue_scans(cold)` — is now FALSE, and deleting it would have left
    `== []` passing forever over a detector that can no longer match anything.

    ⚠ **The non-vacuity therefore moved rather than going away**: the detector is exercised directly
    against `compute_propagation`, which still scans, so it is proved able to fire in the same run
    that proves the route does not fire it.

    ⚠ Still clears both caches, so the cold half is really cold whatever order the suite runs in.
    """
    clear_calibration_cache()
    clear_propagation_cache()
    cold = _measure(application, f"/api/v1/species/{species}/loci/{label}/function")
    assert _gene_grain_reads(cold) == [], (
        f"a COLD Function tab still reads a one-row-per-gene table: {_gene_grain_reads(cold)}"
    )

    warm = _measure(application, f"/api/v1/species/{species}/loci/{label}/function")
    assert _gene_grain_reads(warm) == []
    warm.assert_at_most(FUNCTION_TAB_WARM, what="a warm Function tab, at the route")
    # ⚠ And the cold tab is no longer dramatically more expensive than the warm one: the only thing
    # the caches now save is a handful of bounded reads, not a catalogue scan.
    cold.assert_at_most(FUNCTION_TAB_COLD, what="a COLD Function tab, at the route")


def test_the_gene_grain_scan_detector_can_still_FIRE_so_the_test_above_is_not_vacuous(application):
    """⛔ A detector that cannot match anything reports an empty list for every route.

    `compute_propagation` is the scan the Function tab used to pay. It still exists — it is the
    oracle `tests/test_locus_annotation_rollups.py` compares the stored rollups against — so it can
    be run deliberately here to show that `_gene_grain_reads` recognises it. ⚠ If this ever goes
    red because `compute_propagation` has gone, the assertions above stop meaning anything and need
    a different witness, not a deletion.
    """
    from sqlalchemy.orm import Session

    from bacatlas_backend.instruments.within_node_propagation import compute_propagation
    from bacatlas_backend.models.enumerations import AnnotationKind

    engine = application.extensions["bacatlas_database"].engine
    with SqlCostOracle(engine) as report, Session(engine) as session:
        compute_propagation(session, pangenome_id=1, annotation_kind=AnnotationKind.COG_ORTHOGROUP)
    assert _gene_grain_reads(report), "the detector did not fire on a statement that does scan"
    assert _gene_grain_catalogue_scans(report), "the catalogue-scan detector did not fire either"
