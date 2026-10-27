"""The run layer: exclusivity decoding, the step cross-check, and the graded headline.

The three things that can go wrong silently here are all about a *reading* rather than a value: the
`-excl` prefix hazard, a ρ comparison that fires on every load, and an audit summary that belongs to
another model. Each has a test that would fail if it were got wrong in the plausible direction.
"""

import pytest

from bacatlas_backend.ingest.ingest_nuna_model_registry import expand_model_steps
from bacatlas_backend.ingest.ingest_pangenome_run import (
    KNOWN_INPUT_KEYS,
    PangenomeIngestError,
    _step_name_from_manifest,
    audit_evaluation_rows,
    build_step_rows,
    cross_check_steps,
    guard_audit_describes_this_model,
    guard_input_keys_have_not_drifted,
    resolve_exclusivity_form,
)
from bacatlas_backend.models.enumerations import (
    EvaluationKind,
    ExclusivityForm,
    ExclusivityFormSource,
    RhoRule,
)
from bacatlas_backend.models.nuna_model import NunaModelStep

nuna_pipeline = pytest.importorskip(
    "nuna.tl.cluster.nuna_pipeline", reason="nuna is the ingest extra, not a serving dependency"
)
run_manifest = pytest.importorskip("nuna.tl.cluster.run_manifest")
MODELS = nuna_pipeline.MODELS

ECOLI_RUN = "ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0"
STANDARD_RUN = ECOLI_RUN.replace("-excl_", "-exclLOGP_")


# ── the prefix hazard ──────────────────────────────────────────────────────────────────────────
def test_the_published_runs_decode_as_DAMPED_because_their_token_is_bare_excl():
    form, source, token = resolve_exclusivity_form(ECOLI_RUN, None)
    assert form is ExclusivityForm.DAMPED_EXCLUSION
    assert token == "-excl"
    assert source is ExclusivityFormSource.RUN_ID_TOKEN


def test_exclLOGP_decodes_as_the_STANDARD_form_and_keeps_its_whole_token():
    """⛔ `-excl` is a PREFIX of `-exclLOGP`. A shortest-first walk returns the wrong form AND a
    truncated token, and the row would then satisfy the CHECK against the wrong pair."""
    form, _, token = resolve_exclusivity_form(STANDARD_RUN, None)
    assert form is ExclusivityForm.EXCLUSION
    assert token == "-exclLOGP"


def test_a_run_id_with_no_exclusivity_token_is_NONE_and_not_a_NULL():
    """The step-2/step-3b parent case. `form_from_run_id` returns None; a NOT NULL column cannot."""
    form, _, token = resolve_exclusivity_form("ecoli_esm_clever_exploded_mmseq0.98_k100_res0.98_seed0", None)
    assert form is ExclusivityForm.NONE
    assert token is None


def test_a_manifest_that_records_a_form_BEATS_the_token():
    """A token is a guess about a default; a manifest is a record."""
    form, source, token = resolve_exclusivity_form(ECOLI_RUN, {"exclusivity_form": "logp"})
    assert form is ExclusivityForm.EXCLUSION, "the legacy alias `logp` canonicalises to exclusion"
    assert source is ExclusivityFormSource.RUN_MANIFEST
    assert token == "-excl", "the literal token is still recorded, so the disagreement stays visible"


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


# ── the cross-check ────────────────────────────────────────────────────────────────────────────
def _registry_steps(model_key):
    return [
        NunaModelStep(
            step_ordinal=row.step_ordinal,
            step_name=row.step_name,
            gamma=row.gamma,
            rho_rule=row.rho_rule,
            rho_ceiling=row.rho_ceiling,
            uses_exclusivity=row.uses_exclusivity,
        )
        for row in expand_model_steps(MODELS[model_key])
    ]


def test_step_2_does_NOT_report_a_conflict_even_though_its_strings_differ():
    """⛔ The manifest says `ceiling`@1e9 and the registry says `off`@1e9. Comparing the literals
    would fire on every multi-step model on every ingest — and a check that always fires is one
    nobody reads."""
    rows = build_step_rows(
        [
            {
                "pipeline": "density_potts", "node_mode": "mmseq", "_depth": 0, "_assign": "/a/2.tsv",
                "gamma": 0.98, "rho_rule": "ceiling", "rho_ceiling": 1e9, "exclusivity": False,
            }
        ]
    )
    assert cross_check_steps(rows, _registry_steps("nuna5")) == []


def test_a_gamma_that_actually_differs_IS_reported():
    rows = build_step_rows(
        [
            {
                "pipeline": "density_potts", "node_mode": "precluster", "_depth": 0, "_assign": "/a/4.tsv",
                "gamma": 0.1, "rho_rule": "ceiling", "rho_ceiling": 1.0, "exclusivity": True,
            }
        ]
    )
    messages = cross_check_steps(rows, _registry_steps("nuna5"))
    assert len(messages) == 1
    assert "γ=0.1" in messages[0] and "0.5" in messages[0]


def test_a_rho_rule_that_actually_differs_IS_reported():
    rows = build_step_rows(
        [
            {
                "pipeline": "grouped_rho_cpm", "_depth": 0, "_assign": "/a/3b.tsv",
                "gamma": 0.5, "rho_rule": "ceiling", "rho_ceiling": 1.0, "exclusivity": False,
            }
        ]
    )
    messages = cross_check_steps(rows, _registry_steps("nuna5"))
    assert len(messages) == 1 and "ρ rule" in messages[0]


def test_the_published_ecoli_manifest_agrees_with_nuna4_on_every_step_it_describes(ecoli_artifacts):
    """The real file, not a fixture — and it is a ONE-step chain locally, which the test says."""
    chain = run_manifest.load_run_chain(ecoli_artifacts.assignment)
    assert len(chain) == 1, "precluster_assign is an absolute CSD3 path, so the chain does not walk here"
    rows = build_step_rows(chain)
    assert [row["step_name"] for row in rows] == ["step4"]
    assert rows[0]["rho_rule"] is RhoRule.CEILING
    assert cross_check_steps(rows, _registry_steps("nuna4")) == []


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
