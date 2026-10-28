"""The run layer: how a step is named and ordered, and the graded headline.

⛔ **Two whole subjects left this file with the model axis (2026-10-27)**, and the loss is named
rather than quietly absorbed:

* **the exclusivity decode** — four tests on the `-excl`-is-a-prefix-of-`-exclLOGP` hazard. Nothing
  in this service reads a run_id token any more; the id is stored verbatim and nuna's own suite
  owns that vocabulary. ⚠ The structurally identical trap in the CATALOGUE KEYS is still live and
  is pinned in `test_catalogue_addressing.py`.
* **the step cross-check** — three tests on comparing each manifest step with `nuna_model_step`,
  including the ρ comparison that would otherwise fire on every multi-step load. The registry is
  gone, so there is nothing to compare against; `pangenome_step` records what the run actually did.

What remains is about a *reading* rather than a value: a step named by what it consumed, steps
ordered deepest-first, and an audit summary that belongs to another model. Each has a test that
fails if it is got wrong in the plausible direction.
"""

import pytest

from bacatlas_backend.ingest.ingest_pangenome_run import (
    KNOWN_INPUT_KEYS,
    PangenomeIngestError,
    _step_name_from_manifest,
    audit_evaluation_rows,
    build_step_rows,
    guard_audit_describes_this_model,
    guard_input_keys_have_not_drifted,
)
from bacatlas_backend.models.enumerations import EvaluationKind

nuna_pipeline = pytest.importorskip(
    "nuna.tl.cluster.nuna_pipeline", reason="nuna is the ingest extra, not a serving dependency"
)
run_manifest = pytest.importorskip("nuna.tl.cluster.run_manifest")
MODELS = nuna_pipeline.MODELS

ECOLI_RUN = "ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0"
STANDARD_RUN = ECOLI_RUN.replace("-excl_", "-exclLOGP_")


# ── naming a realised step ─────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    ("manifest", "expected"),
    [
        ({"pipeline": "density_potts", "node_mode": "mmseq"}, "step2"),
        ({"pipeline": "density_potts", "node_mode": "precluster"}, "step4"),
        ({"pipeline": "grouped_rho_cpm"}, "step3b"),
        ({"pipeline": "grouped_rho_cpm", "atoms_assign": "/x/y.tsv"}, "step3c"),
    ],
)
def test_a_step_is_named_by_what_it_CONSUMED_not_by_its_filename(manifest, expected):
    assert _step_name_from_manifest(manifest) == expected


def test_steps_are_ordered_earliest_first_which_is_DEEPEST_first():
    """`_provenance_from_chain`'s rule: deeper in the chain is earlier in the pipeline."""
    chain = [
        {"pipeline": "density_potts", "node_mode": "precluster", "_depth": 0, "_assign": "/a/4.tsv"},
        {"pipeline": "grouped_rho_cpm", "_depth": 1, "_assign": "/a/3b.tsv"},
        {"pipeline": "density_potts", "node_mode": "mmseq", "_depth": 2, "_assign": "/a/2.tsv"},
    ]
    rows = build_step_rows(chain)
    assert [row["step_name"] for row in rows] == ["step2", "step3b", "step4"]
    assert [row["step_ordinal"] for row in rows] == [2, 3, 4]


# ── the graded headline ────────────────────────────────────────────────────────────────────────
def test_every_summary_key_becomes_a_row_and_a_None_becomes_a_row_that_SAYS_so():
    rows = audit_evaluation_rows(
        {"n_clusters_total": 17531, "over_merge_gene_rate": 0.000343, "label": "ecoli_nuna4",
         "single_family_divergent_n_clusters": None},
        __import__("pathlib").Path("/a/summary.json"),
    )
    by_name = {row["metric_name"]: row for row in rows}
    assert len(rows) == 4, "a dropped key makes 'not measured' and 'not read' the same absence"
    assert by_name["n_clusters_total"]["numeric_value"] == 17531.0
    assert by_name["label"]["numeric_value"] is None
    assert by_name["label"]["detail"] == "ecoli_nuna4"
    assert by_name["single_family_divergent_n_clusters"]["detail"] == "not measured by this audit run"
    assert all(row["evaluation_kind"] is EvaluationKind.ACCESSORY_AUDIT for row in rows)


def test_an_audit_summary_for_another_model_is_a_HARD_error(ecoli_artifacts):
    with pytest.raises(PangenomeIngestError, match="grades"):
        guard_audit_describes_this_model({"label": "kp_nuna4_something"}, ecoli_artifacts.model_label)


def test_the_real_ecoli_summary_passes_its_own_guard_and_yields_rows(ecoli_artifacts):
    import json

    summary = json.loads(ecoli_artifacts.audit_summary.read_text())
    guard_audit_describes_this_model(summary, ecoli_artifacts.model_label)
    rows = audit_evaluation_rows(summary, ecoli_artifacts.audit_summary)
    assert len(rows) == len(summary) >= 70
    assert {row["metric_name"] for row in rows} >= {"n_clusters_total", "synteny_only_n_clusters"}


# ── the constraint guard ───────────────────────────────────────────────────────────────────────
def test_the_input_key_guard_matches_nuna_today_and_names_the_drift_if_it_stops():
    guard_input_keys_have_not_drifted(run_manifest.INPUT_KEYS)
    assert tuple(run_manifest.INPUT_KEYS) == KNOWN_INPUT_KEYS
    with pytest.raises(PangenomeIngestError, match="CHECK admits only"):
        guard_input_keys_have_not_drifted((*KNOWN_INPUT_KEYS, "step3c_assign"))


# ── reconciliation: a zero needs its denominator ───────────────────────────────────────────────
def test_a_run_with_no_PUBLISHED_TRIPLE_is_not_reconcilable_and_says_so():
    """⛔ `reconcile_pangenome` returns no differences for an unknown run_id — having compared nothing.

    Two runs have a checked-in triple. Every cohort loaded from the browser's cohort axis onward
    has none, so the empty list is about to become the NORMAL result — and `main` printed
    *"reconciled: every count the loader reported is what the tables hold"* on the strength of it.
    `pangenome_is_reconcilable` is that zero's denominator: it is what lets the caller tell
    *nothing differed* from *nothing was compared*.
    """
    from pathlib import Path

    from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
    from bacatlas_backend.ingest.ingest_command_line import pangenome_is_reconcilable
    from bacatlas_backend.ingest.published_catalogues import PUBLISHED_CATALOGUES

    published = PUBLISHED_CATALOGUES[0]
    known = CatalogueArtifacts(
        data_root=Path("/nowhere"),
        set_key=published.set_key,
        model_label=published.model_label,
        run_id=published.run_id,
        species_key=published.species_key,
    )
    assert pangenome_is_reconcilable(known) is True

    # ⚠ A plausible FUTURE run of the same species and the same model — only the run differs. It is
    # the shape every cohort load takes, not a nonsense string, so the test fails for the reason it
    # names rather than because the input was obviously unlike a real one.
    unknown = CatalogueArtifacts(
        data_root=Path("/nowhere"),
        set_key=published.set_key,
        model_label=published.model_label,
        run_id=published.run_id + "_regenerated",
        species_key=published.species_key,
    )
    assert pangenome_is_reconcilable(unknown) is False
