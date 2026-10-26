"""`nuna_pipeline.MODELS` → `nuna_model` + `nuna_model_step`.

⛔ **The registry is the executable definition of what a model IS.** `CLAUDE.md` says so outright:
where a doc and the registry disagree, the registry wins, because it is the thing that actually ran.
So these rows are read from `ModelSpec` and from the argument lists `run_nuna_chain` builds — never
from a table in a document, and never parsed out of a run_id.

⭐ **Every parameter that varies by step lands on the STEP.** γ, the ρ rail, the representation, the
node mode and whether the step carries the exclusivity weight are all per-step; `nuna5` has five
steps and no two share a γ/ρ pair. A schema that put any of them on the model would assert something
false about every multi-step chain — which is the whole reason this table exists.

⚠ **Two things here are DERIVED, and each says so in its row.**

1. **`rho_rule = OFF` at step 2 is a reading of the ceiling, not a flag.** `nuna_pipeline` passes
   step 2 `--rho-ceiling 1e9` and no `--rho-rule` at all, so `density_potts`' own default (`ceiling`)
   is what the run manifest records. `CLAUDE.md`'s table and `run_manifest.rho_phrase` both call that
   state *ρ off* — `rho_phrase` returns `"ρ off"` for any ceiling ≥ 1e6 — so that is the model's
   semantic and it is what the registry row carries. The manifest's literal `"ceiling"` is what
   `pangenome_step` carries. ⛔ **They are not a disagreement**, and the cross-check must normalise
   through :func:`rho_rule_is_effectively_off` rather than comparing the two strings, or every
   multi-step model would report a spurious step-2 conflict on every ingest.
2. **`edge_weight_expression` follows the exclusivity form**, because the two forms enter the
   expression at different points: `damped_exclusion` multiplies *after* `mᵢmⱼ` exactly as published,
   while `f` caps at 1 and so must hit `s` *before* `mᵢmⱼ` or γ means something different at every
   node size.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import (
    EmbeddingRepresentation,
    ExclusivityForm,
    ExclusivityFormSource,
    RhoRule,
)
from bacatlas_backend.models.nuna_model import NunaModel, NunaModelStep

#: A ρ ceiling at or above this is not a constraint — it is how a deliberate over-merge step says
#: *no rail at all*. The threshold is `run_manifest.rho_phrase`'s own, quoted rather than re-chosen.
RHO_OFF_CEILING = 1e6

#: Step 2's ceiling, as `nuna_pipeline` passes it: the literal string `"1e9"`.
#: ⛔ `nuna_model_step.rho_ceiling` is `Numeric(18, 4)` **because of this value** — the original
#: `Numeric(12, 4)` allows 8 digits left of the point and this needs 10, so ingesting the step 2 of
#: every multi-step model raised `NumericValueOutOfRange`.
STEP2_RHO_CEILING = 1e9

#: `-c` is a pipeline constant and is never encoded in a run_id — `_dedup_str` says so in as many
#: words: *"coverage = pipeline default, not in run_id"*.
DEDUP_COVERAGE = 0.95

#: `accessory_audit_run._GRAPH_CONSTRUCTION`'s vocabulary, for the one method these chains use.
MERGE_ALGORITHM = "CPM (leidenalg CPMVertexPartition)"


def rho_rule_is_effectively_off(rho_rule: str | None, rho_ceiling: float | None) -> bool:
    """Whether a (rule, ceiling) pair means *ρ off*, by `run_manifest.rho_phrase`'s own rule.

    ⛔ The single place this judgement is made. Step 2 records `rho_rule = "ceiling"` in its manifest
    and `rho_ceiling = 1e9`, and that pair **is** the OFF state — comparing the strings alone reports
    a conflict on every multi-step model, every time it is ingested.
    """
    if rho_ceiling is not None and float(rho_ceiling) >= RHO_OFF_CEILING:
        return True
    return rho_rule is not None and str(rho_rule).lower() == RhoRule.OFF.value


@dataclass(frozen=True)
class StepRow:
    """One expanded step of a chain — the shape `nuna_model_step` stores."""

    step_ordinal: int
    step_name: str
    representation: EmbeddingRepresentation | None
    gamma: float | None
    rho_rule: RhoRule | None
    rho_ceiling: float | None
    rho_exit: bool | None
    rho_continue: str | None
    node_mode: str | None
    uses_exclusivity: bool
    edge_weight_expression: str | None
    node_sizes: str | None


#: ⛔⛔ **nuna's canonical hypergeometric form → the value `nuna_model.exclusivity_form` STORES.**
#: The two vocabularies are different and the difference is a date: nuna renamed the field
#: `ModelSpec.exclusivity_form` → `hypergeometric_form` and its values `exclusion` →
#: `hypergeometric`, `damped_exclusion` → `damped` on **2026-10-23**, recording the old spellings in
#: `hypergeometric.HYPERGEOMETRIC_FORM_ALIASES` as legacy. This schema still stores the legacy
#: spellings, and it must: `nuna_model` is on the never-rename list, its values are in published
#: quotes, and `pangenome.run_id_exclusivity_token` cross-checks against them.
#:
#: ⚠ So this dict is the whole of the translation, and it is one line per form on purpose. When
#: nuna renames a third time, `spec.hypergeometric_form` raises `AttributeError` on the first model
#: of the first ingest — which is exactly how this breakage was found, 12-22 h after it landed, and
#: is far better than a tolerant `getattr` that would have written a wrong form for every model.
EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC: dict[str, ExclusivityForm] = {
    "hypergeometric": ExclusivityForm.EXCLUSION,
    "damped": ExclusivityForm.DAMPED_EXCLUSION,
}


class ModelNotRecordable(RuntimeError):
    """A `ModelSpec` carries a parameter `nuna_model_step` has no column for.

    ⛔ **Raised rather than absorbed, because the alternative is a row that describes a different
    model.** `nuna_model_step` claims to be the executable definition of a chain; a step row that
    says *ρ ceiling 1.0* for a pass that actually ran with the **ρ-free floor** describes a stricter
    rail than the model used — and the floor is what decides whether an edge exists at all
    (`nuna/CLAUDE.md`: *"the floor decides whether an edge exists, and disjointness decides what it
    carries"*). An omission here is a falsehood about the model's identity, which is the one thing
    this table exists not to be.
    """


#: ⛔ Model parameters nuna's `ModelSpec` carries that this schema has **no column for**. Each is a
#: `(field, what it changes)` pair, and a spec setting any of them is refused by name.
#:
#: ⚠ Both are the **ρ-free floor** (nuna added them 2026-10-24): the rail refuses a merge that would
#: take ρ above 1, but for a small node ρ = 1 by near-certainty and carries no information, so the
#: floor lifts the rail when both endpoints are small enough that spanning distinct genomes is
#: unsurprising. It is cohort-adaptive — m < 24 at 100 genomes, < 285 at 13,573 — so it is not a
#: constant that could be folded into `rho_ceiling`. ⛔ **Recording it is a migration and a judgement
#: about what a model IS, not a mapping entry**, so it is surfaced rather than invented here.
UNRECORDABLE_SPEC_FIELDS: dict[str, str] = {
    "rho_free_p": "the ρ-free floor on the Bacformer CONTEXT pass",
    "rho_free_p_merge": "the ρ-free floor on the MERGE pass",
}


def unrecordable_parameters(spec) -> dict[str, str]:
    """Which `UNRECORDABLE_SPEC_FIELDS` this spec actually sets, as `{field: what it changes}`."""
    return {
        field: description
        for field, description in UNRECORDABLE_SPEC_FIELDS.items()
        if getattr(spec, field, None) is not None
    }


def exclusivity_form_of(spec) -> ExclusivityForm:
    """A `ModelSpec`'s hypergeometric form, as the value this schema stores for it.

    ⚠ Normalised through **nuna's own** `hypergeometric.canonical_form` rather than against a copy of
    its alias table, so a spec carrying a legacy spelling resolves the same way here as it does in
    the run that produced it. A second alias table is how two spellings of one model get into one
    record.
    """
    from nuna.tl.cluster.hypergeometric import canonical_form  # noqa: PLC0415 — nuna is ingest-only

    canonical = canonical_form(spec.hypergeometric_form)
    if canonical not in EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC:
        raise ValueError(
            f"nuna names a hypergeometric form {canonical!r} that this schema has no value for; "
            f"`nuna_model.exclusivity_form` stores {sorted(EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC)} and "
            "adding one is a migration, not a mapping entry"
        )
    return EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC[canonical]


def _edge_weight(uses_exclusivity: bool, form: ExclusivityForm) -> str:
    """The edge weight a step actually computes, spelled as the method record spells it."""
    if not uses_exclusivity:
        return "s·mᵢmⱼ"
    return "f·mᵢmⱼ" if form is ExclusivityForm.EXCLUSION else "s·mᵢmⱼ·(1−P₀)"


def expand_model_steps(spec) -> list[StepRow]:
    """A `ModelSpec` → one row per step, following `run_nuna_chain`'s own argument lists.

    ⚠ Read against `nuna_pipeline.run_nuna_chain` and not against any table in a document. The steps
    it builds are, in order: the MMseqs dedup prerequisite (not a CPM pass, hence a row with no γ and
    no rail), the ESM over-merge, the Bacformer context split, an optional within-group exclusivity
    pass, and then either one global merge or the split 4b → 4c → (4d).
    """
    form = exclusivity_form_of(spec)
    rows = [
        StepRow(
            step_ordinal=1,
            step_name="step1_dedup",
            representation=None,
            # ⛔ Not a CPM pass at all: MMseqs `close_seq` is a PREREQUISITE that runs in a different
            # environment and has no manifest. So every CPM parameter is NULL here, and that is the
            # honest shape rather than a zero.
            gamma=None,
            rho_rule=None,
            rho_ceiling=None,
            rho_exit=None,
            rho_continue=None,
            node_mode=None,
            uses_exclusivity=False,
            edge_weight_expression=None,
            node_sizes=None,
        ),
    ]
    # ⛔ **`has_esm` is nuna's OWN predicate, not `gamma2 is not None` re-derived here.** Both models
    # of record (`base`, `hypergeo`) have no ESM step, and the consequence is not just a missing row:
    # the ESM pass is what gives the Bacformer pass a homology GROUP to run within, so without it the
    # Bacformer pass runs GLOBALLY over the mmseqs alleles — `--node-mode mmseq`, not `precluster`.
    # A row that kept `precluster` would name a pass that never ran.
    if spec.has_esm:
        rows.append(
            StepRow(
                step_ordinal=len(rows) + 1,
                step_name="step2",
                representation=EmbeddingRepresentation.ESM,
                gamma=float(spec.gamma2),
            # See the module docstring: OFF is the model's semantic for a 1e9 ceiling, and the
            # manifest's literal "ceiling" is not a disagreement with it.
                rho_rule=RhoRule.OFF,
                rho_ceiling=STEP2_RHO_CEILING,
                rho_exit=None,
                rho_continue=None,
                node_mode="mmseq",
                uses_exclusivity=False,
                edge_weight_expression=_edge_weight(False, form),
                node_sizes="gene_count",
            )
        )
    rows.append(
        StepRow(
            step_ordinal=len(rows) + 1,
            step_name="step3b",
            representation=EmbeddingRepresentation.BACFORMER,
            gamma=float(spec.gamma3b),
            rho_rule=RhoRule(spec.rho_rule_3b),
            rho_ceiling=1.0,
            rho_exit=True,
            rho_continue=None,
            node_mode="precluster" if spec.has_esm else "mmseq",
            # ⛔ NO exclusivity here, and it is a precondition rather than a preference:
            # `pairwise_max` deliberately admits colliding pairs, and the exclusivity weight reads
            # only the endpoints' genome counts — meaningful only where every edge is already
            # genome-disjoint. `grouped_rho_cpm` refuses the combination outright.
            uses_exclusivity=False,
            edge_weight_expression=_edge_weight(False, form),
            node_sizes="gene_count",
        )
    )
    ordinal = len(rows) + 1
    if spec.has_3c:
        rows.append(
            StepRow(
                step_ordinal=ordinal,
                step_name="step3c",
                representation=EmbeddingRepresentation.BACFORMER,
                gamma=float(spec.gamma3c),
                rho_rule=RhoRule.CEILING,
                rho_ceiling=1.0,
                rho_exit=True,
                rho_continue=spec.rho_continue,
                node_mode="precluster",
                uses_exclusivity=True,
                edge_weight_expression=_edge_weight(True, form),
                node_sizes="gene_count",
            )
        )
        ordinal += 1

    if not spec.has_split_step4:
        if spec.gamma4 is None:
            # ⚠ The nuna2 shape, which `base` is: there is no merge pass at all and the Bacformer
            # context pass above is the final partition. A `step4` row here would invent one.
            return rows
        rows.append(
            StepRow(
                step_ordinal=ordinal,
                step_name="step4",
                representation=EmbeddingRepresentation.BACFORMER,
                gamma=float(spec.gamma4),
                rho_rule=RhoRule.CEILING,
                rho_ceiling=1.0,
                rho_exit=True,
                rho_continue=spec.rho_continue,
                node_mode="precluster",
                uses_exclusivity=True,
                edge_weight_expression=_edge_weight(True, form),
                node_sizes="gene_count",
            )
        )
        return rows

    # The split step 4. ⚠ 4b runs with NO exclusivity, which is what lets it take `pairwise_max`:
    # the ceiling rail is a precondition FOR the weight, so with no weight there is nothing to
    # protect. `spec.gamma4` is ignored entirely when the split is set — asserted, not assumed.
    split = [
        ("step4b", spec.gamma4b, RhoRule(spec.rho_rule_4b), False),
        ("step4c", spec.gamma4c, RhoRule.CEILING, True),
    ]
    if spec.gamma4d is not None:
        split.append(("step4d", spec.gamma4d, RhoRule.CEILING, True))
    for name, gamma, rule, excl in split:
        rows.append(
            StepRow(
                step_ordinal=ordinal,
                step_name=name,
                representation=EmbeddingRepresentation.BACFORMER,
                gamma=float(gamma),
                rho_rule=rule,
                rho_ceiling=1.0,
                rho_exit=True,
                # ⚠ NULL where the pass carries no exclusivity weight, because `_merge_pass` passes
                # `--rho-continue` only inside its `if excl` block. Storing the spec's value on 4b
                # would claim a parameter that pass never received.
                rho_continue=spec.rho_continue if excl else None,
                node_mode="precluster",
                uses_exclusivity=excl,
                edge_weight_expression=_edge_weight(excl, form),
                node_sizes="gene_count",
            )
        )
        ordinal += 1
    return rows


def model_key_for_label(label: str, models: dict) -> str:
    """The registry key whose `ModelSpec.label` is `label`, or a failure naming what is known.

    ⚠ `label` here is the bare `ModelSpec.label`, **not** the audit's `--model-label`, which prefixes
    it with the set token (`ecoli_nuna4_…`). Callers holding the latter strip the prefix first — see
    :func:`model_key_for_audit_label` — because two sets share one model and the registry does not
    know about sets at all.
    """
    for key, spec in models.items():
        if spec.label == label:
            return key
    raise KeyError(
        f"no MODELS entry has label {label!r}. Known labels: {sorted(s.label for s in models.values())}. "
        "A model is identified by its registry entry and never by a filename, so this is raised "
        "rather than invented."
    )


def model_key_for_audit_label(audit_label: str, set_key: str, models: dict) -> str:
    """`ecoli_nuna4_…` + `ecoli` → the registry key. The set token is a prefix, never part of it."""
    prefix = f"{set_key}_"
    if not audit_label.startswith(prefix):
        raise ValueError(
            f"audit model label {audit_label!r} does not begin with the set token {prefix!r}; "
            "every artifact of this catalogue is addressed by one or the other."
        )
    return model_key_for_label(audit_label[len(prefix) :], models)


def ingest_nuna_model(
    session: Session, model_key: str, *, models: dict, registry_git_sha: str | None = None
) -> int:
    """Ensure one `nuna_model` and its steps exist; return the `nuna_model_id`.

    Idempotent by `model_key`, and **it rewrites the steps** when the row already exists. That is the
    right blast radius: the registry is code, so a step whose γ changed upstream must not survive in
    a table that claims to be its executable definition — while a model already referenced by an
    ingested pangenome keeps its id, so nothing downstream is orphaned.
    """
    spec = models[model_key]
    unrecordable = unrecordable_parameters(spec)
    if unrecordable:
        raise ModelNotRecordable(
            f"`{model_key}` ({spec.label}) sets "
            + ", ".join(f"{field}={getattr(spec, field)!r} — {what}" for field, what in unrecordable.items())
            + ", and `nuna_model_step` has no column for it. Recording a step row without it would "
            "claim a rail the pass did not run under, so the row is refused until the column exists."
        )
    existing = session.execute(
        select(NunaModel).where(NunaModel.model_key == model_key)
    ).scalar_one_or_none()

    form = exclusivity_form_of(spec)
    steps = expand_model_steps(spec)

    if existing is None:
        existing = NunaModel(model_key=model_key)
        session.add(existing)
    existing.label = spec.label
    existing.dedup_identity = float(spec.dedup_id)
    existing.dedup_coverage = DEDUP_COVERAGE
    existing.knn_k = int(spec.knn_k)
    existing.random_seed = int(spec.seed)
    existing.exclusivity_form = form
    existing.merge_algorithm = MERGE_ALGORITHM
    existing.graph_construction = "density(s) · ρ≤1 · explode(mᵢmⱼ)"
    existing.step_count = len(steps)
    existing.registry_source_git_sha = registry_git_sha
    session.flush()

    existing.steps.clear()
    session.flush()
    for row in steps:
        session.add(
            NunaModelStep(
                nuna_model_id=existing.nuna_model_id,
                step_ordinal=row.step_ordinal,
                step_name=row.step_name,
                representation=row.representation,
                gamma=row.gamma,
                rho_rule=row.rho_rule,
                rho_ceiling=row.rho_ceiling,
                rho_exit=row.rho_exit,
                rho_continue=row.rho_continue,
                node_mode=row.node_mode,
                uses_exclusivity=row.uses_exclusivity,
                edge_weight_expression=row.edge_weight_expression,
                node_sizes=row.node_sizes,
                # ⚠ These rows come from the REGISTRY, so they say so. `pangenome_step` reads the
                # same shape from the manifests, and the two are cross-checked rather than merged.
                provenance_source=ExclusivityFormSource.MODEL_REGISTRY,
            )
        )
    session.flush()
    return existing.nuna_model_id


def ingest_model_registry(
    session: Session, *, models: dict, registry_git_sha: str | None = None
) -> dict[str, int]:
    """Every RECORDABLE entry in the registry → `{model_key: nuna_model_id}`.

    All of them, not just the published one: the registry is *the vocabulary of the whole system*,
    and a graded run that was never published still needs a model to point at.

    ⛔ **A model this schema cannot record is SKIPPED and NAMED, not fatal — and the blast radius is
    the point.** `nuna_model_step` has no column for the ρ-free floor, so three registry entries
    cannot be written faithfully (see :data:`UNRECORDABLE_SPEC_FIELDS`). Aborting here made the
    *whole pangenome layer* unloadable because of a parameter on a model the catalogue being
    ingested does not use — which is how a registry drift became a total ingest outage.

    ⭐ The failure boundary moves to where it belongs: `load_pangenome_layer` resolves its OWN model
    with `scalar_one()`, so loading a catalogue whose model is unrecordable still fails, loudly, on
    that model. Skipping one you do not need does not.

    ⚠ Every refusal is printed on every run. A missing row in the authority table is not something a
    reader should have to infer from a `KeyError` three screens later.
    """
    recorded, refused = {}, []
    for key in models:
        try:
            recorded[key] = ingest_nuna_model(
                session, key, models=models, registry_git_sha=registry_git_sha
            )
        except ModelNotRecordable as refusal:
            refused.append(str(refusal))
    for refusal in refused:
        print(f"  ⚠ model NOT recorded: {refusal}")
    if refused:
        print(
            f"  ⚠ {len(recorded)} of {len(models)} registry models recorded; "
            f"{len(refused)} refused for a parameter this schema has no column for"
        )
    return recorded
