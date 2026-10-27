"""One run of one model over one collection — `pangenome`, its steps, its lineage and its grades.

⛔ **Read from the run manifests, never parsed from a run_id.** Run ids emit tokens *non-default
only*, so an absent token means "whatever the default was when that run was made" and is not
recoverable from the id at all — two models differing only in a ρ rule have produced byte-identical
provenance. `runs/{run_id}.json` is the only thing that can describe a run, and it is stored verbatim
alongside the columns read out of it.

⭐ **The chain walks locally through `mirrored_manifest_chain`, and that is what makes the
provenance real.** A manifest names its parent by the ABSOLUTE path of the machine that wrote it, so
a bare `load_run_chain` reaches exactly one manifest here and records a four-step model as a
one-step one — silently, a short chain being what a genuine one-step model produces too. Measured:
without the adapter the footer came out as 7 rows from 1 manifest, of which the only pipeline row
read *"step 1 · dedup: not recorded in any manifest"*; with it, 9 rows from 3 manifests, **byte-
identical to the published page's footer** on both species.

⚠ A parent run that has no `pangenome` row of its own is still a `pangenome_input_edge` with
`parent_run_id` set and `parent_pangenome_id` NULL — a gap **admitted**, which is the shape that
column exists for.

⛔ **`rho_rule` is compared through `rho_rule_is_effectively_off`, never as a string.** Step 2 runs
with `--rho-ceiling 1e9` and no `--rho-rule`, so its manifest records `density_potts`' default
`"ceiling"` while the model means *ρ off*. Comparing the two literals reports a conflict on every
multi-step model on every ingest, and a cross-check that always fires is a cross-check nobody reads.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_nuna_model_registry import rho_rule_is_effectively_off
from bacatlas_backend.ingest.mirrored_manifest_chain import manifests_resolved_against
from bacatlas_backend.models.enumerations import (
    EmbeddingRepresentation,
    EvaluationKind,
    ExclusivityForm,
    ExclusivityFormSource,
    RhoRule,
)
from bacatlas_backend.models.nuna_model import NunaModelStep
from bacatlas_backend.models.pangenome import (
    Pangenome,
    PangenomeEvaluation,
    PangenomeInputEdge,
    PangenomeStep,
)

#: `run_manifest.INPUT_KEYS`, restated so :func:`guard_input_keys_have_not_drifted` can compare.
#: ⛔ `pangenome_input_edge` carries a `CHECK (input_key IN (…))` over exactly these three. If nuna
#: adds a fourth, every lineage edge of the new shape is refused by Postgres with a constraint name
#: and no explanation — so the drift is caught here, where it can say what actually happened.
KNOWN_INPUT_KEYS = ("precluster_assign", "atoms_assign", "step2_assign")


class PangenomeIngestError(RuntimeError):
    """A run that cannot be recorded without inventing part of its identity."""


@dataclass
class PangenomeRunReport:
    """What was recorded about the run, and every gap that was admitted rather than filled."""

    run_id: str
    pangenome_id: int
    steps_written: int
    input_edges_written: int
    unresolved_parent_run_ids: list[str]
    evaluation_rows_written: int
    exclusivity_form: str
    exclusivity_form_source: str
    step_disagreements: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def render(self) -> str:
        """The report as the CLI prints it — what was written, then what was not."""
        lines = [
            f"pangenome {self.run_id}",
            f"  id                    {self.pangenome_id}",
            f"  steps (from manifest) {self.steps_written}",
            f"  lineage edges         {self.input_edges_written}"
            + (f" ({len(self.unresolved_parent_run_ids)} parent runs not ingested)"
               if self.unresolved_parent_run_ids else ""),
            f"  evaluation rows       {self.evaluation_rows_written}",
            f"  exclusivity           {self.exclusivity_form} (from {self.exclusivity_form_source})",
        ]
        for disagreement in self.step_disagreements:
            lines.append(f"  ⚠ STEP DISAGREEMENT {disagreement}")
        for note in self.notes:
            lines.append(f"  NOTE {note}")
        return "\n".join(lines)


# ── exclusivity ────────────────────────────────────────────────────────────────────────────────
def resolve_exclusivity_form(
    run_id: str, manifest: dict | None
) -> tuple[ExclusivityForm, ExclusivityFormSource, str | None]:
    """`(form, how it was determined, the literal token)` for one run.

    ⛔ `hypergeometric.form_from_run_id` is the only permitted decoder of the id, and the token is
    found by walking `HYPERGEOMETRIC_TOKENS` **longest first** — `-excl` is a PREFIX of `-exclLOGP`,
    and the three ways to get that wrong (`in`, a half-`replace`, an `endswith` chain) were all live
    bugs on one day. The manifest wins where it recorded a form, because a token is a guess about a
    default and a manifest is a record.

    ⛔⛔ **nuna renamed this whole vocabulary on 2026-10-23 and the names here did not follow.** The
    module `nuna.tl.cluster.exclusivity` became `hypergeometric`, `EXCLUSIVITY_TOKENS` became
    `HYPERGEOMETRIC_TOKENS`, and — the part that is not a rename — `canonical_form` and
    `form_from_run_id` now return `hypergeometric` / `damped` where they returned `exclusion` /
    `damped_exclusion`. Those two strings are what `ExclusivityForm` STORES, so the translation is
    explicit and shared with the registry loader rather than repeated: see
    `ingest_nuna_model_registry.EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC`. ⚠ The run_id TOKENS are
    unchanged and must stay readable — `-excl` and `-exclLOGP` name the two published catalogues.

    ⚠ Returns `ExclusivityForm.NONE` where nothing names a form. That is the step-2/step-3b parent
    case: those runs carry no exclusivity at all, their ids carry no token, and `form_from_run_id`
    returns `None` for them — which a NOT NULL column cannot store.
    """
    from nuna.tl.cluster.hypergeometric import (  # noqa: PLC0415 — nuna is an ingest-only dependency
        HYPERGEOMETRIC_TOKENS,
        canonical_form,
        form_from_run_id,
    )

    from bacatlas_backend.ingest.ingest_nuna_model_registry import (  # noqa: PLC0415 — one mapping, no cycle
        EXCLUSIVITY_FORM_OF_HYPERGEOMETRIC as STORED_FORM,
    )

    token = next((tok for tok, _ in HYPERGEOMETRIC_TOKENS if tok in run_id), None)

    manifest_form = (manifest or {}).get("exclusivity_form")
    if manifest_form:
        return STORED_FORM[canonical_form(str(manifest_form))], ExclusivityFormSource.RUN_MANIFEST, token

    decoded = form_from_run_id(run_id)
    if decoded is not None:
        return STORED_FORM[decoded], ExclusivityFormSource.RUN_ID_TOKEN, token
    return ExclusivityForm.NONE, ExclusivityFormSource.RUN_ID_TOKEN, token


# ── the steps a run actually took ──────────────────────────────────────────────────────────────
def _step_name_from_manifest(manifest: dict) -> str:
    """A short name for a realised step, from what it consumed — never from its filename.

    ⚠ Deliberately coarse. `run_manifest.stage_name` produces the prose the footer prints; this is
    the join key against `nuna_model_step.step_name`, and the honest mapping is by what the run
    read: a `precluster` merge is a step-4-family pass, a `mmseq` one is step 2.
    """
    pipeline = manifest.get("pipeline")
    node_mode = manifest.get("node_mode")
    if pipeline == "grouped_rho_cpm":
        return "step3c" if manifest.get("atoms_assign") else "step3b"
    if pipeline == "density_potts":
        return "step4" if node_mode == "precluster" else "step2"
    return str(pipeline or "unknown")[:32]


def _rho_rule(value) -> RhoRule | None:
    if value is None:
        return None
    try:
        return RhoRule(str(value))
    except ValueError:
        return None


def build_step_rows(chain: list[dict]) -> list[dict]:
    """The chain's manifests → `pangenome_step` field dicts, earliest step first.

    ⚠ Ordered by `_depth` DESCENDING, because deeper in the chain is *earlier* in the pipeline —
    `_provenance_from_chain`'s own rule. Ordinals are then 2, 3, 4… so they line up with the
    registry's, whose step 1 is the MMseqs dedup that has no manifest at all.
    """
    stages = sorted(chain, key=lambda meta: -meta.get("_depth", 0))
    rows = []
    for ordinal, manifest in enumerate(stages, start=2):
        rep = manifest.get("rep")
        rows.append(
            {
                "step_ordinal": ordinal,
                "step_name": _step_name_from_manifest(manifest),
                "representation": EmbeddingRepresentation(rep) if rep in ("esm", "bacformer") else None,
                "gamma": manifest.get("gamma"),
                "rho_rule": _rho_rule(manifest.get("rho_rule")),
                "rho_ceiling": manifest.get("rho_ceiling"),
                "node_mode": manifest.get("node_mode"),
                "uses_exclusivity": manifest.get("exclusivity"),
                "step_assignment_path": manifest.get("_assign"),
                "gene_count": None,
                "node_count": manifest.get("n_graph_nodes"),
                "wall_seconds": manifest.get("wall_seconds"),
                # ⚠ NULL where the manifest recorded `n/a`. Not measured is not zero.
                "peak_rss_gigabytes": manifest.get("peak_rss_gb"),
                "_manifest": manifest,
            }
        )
    return rows


def cross_check_steps(step_rows: list[dict], registry_steps: list[NunaModelStep]) -> list[str]:
    """Where the run's manifests and the registry disagree, as one line each.

    ⚠ **Recorded, not resolved.** That is how a run whose manifest contradicts the registry gets
    caught rather than silently normalised to whichever the loader happened to prefer. An empty list
    means the two agree on every step both describe.

    ⛔ ρ is compared through the OFF normalisation, not as a string — see the module docstring.
    """
    by_name = {step.step_name: step for step in registry_steps}
    messages = []
    for row in step_rows:
        expected = by_name.get(row["step_name"])
        if expected is None:
            messages.append(f"{row['step_name']}: the registry model has no step of that name")
            continue
        if expected.gamma is not None and row["gamma"] is not None:
            if abs(float(expected.gamma) - float(row["gamma"])) > 1e-9:
                messages.append(
                    f"{row['step_name']}: manifest γ={row['gamma']} but the registry says {expected.gamma}"
                )
        run_off = rho_rule_is_effectively_off(
            row["rho_rule"].value if row["rho_rule"] else None, row["rho_ceiling"]
        )
        registry_off = rho_rule_is_effectively_off(
            expected.rho_rule.value if expected.rho_rule else None, expected.rho_ceiling
        )
        if run_off != registry_off:
            messages.append(
                f"{row['step_name']}: manifest ρ={row['rho_rule']}@{row['rho_ceiling']} "
                f"but the registry says ρ={expected.rho_rule}@{expected.rho_ceiling}"
            )
        elif not run_off and row["rho_rule"] is not None and expected.rho_rule is not None:
            if row["rho_rule"] is not expected.rho_rule:
                messages.append(
                    f"{row['step_name']}: manifest ρ rule {row['rho_rule'].value} "
                    f"but the registry says {expected.rho_rule.value}"
                )
        if expected.uses_exclusivity is not None and row["uses_exclusivity"] is not None:
            if bool(expected.uses_exclusivity) != bool(row["uses_exclusivity"]):
                messages.append(
                    f"{row['step_name']}: manifest exclusivity={row['uses_exclusivity']} "
                    f"but the registry says {expected.uses_exclusivity}"
                )
    return messages


# ── the graded headline ────────────────────────────────────────────────────────────────────────
def audit_evaluation_rows(audit_summary: dict, source_path: Path) -> list[dict]:
    """Every audit-summary key → one `pangenome_evaluation` row, copied verbatim.

    ⛔ **Never recomputed, never rounded, never renamed** — `AUDIT_HEADLINE_KEYS` says exactly that
    of the subset the page quotes, and the reason generalises: a page (or an API) quoting a lookalike
    it derived itself is the failure this table exists to make impossible.

    ⚠ A `None` in the summary is *not measured* and is stored as a row with a NULL value and a
    `detail` saying so — dropping the key would make "the audit did not measure this" and "this
    ingest did not read it" the same absence.
    """
    rows = []
    for key, value in sorted(audit_summary.items()):
        numeric, detail = None, None
        if isinstance(value, bool):
            detail = str(value)
        elif isinstance(value, (int, float)):
            numeric = float(value)
        elif value is None:
            detail = "not measured by this audit run"
        else:
            detail = str(value)
        rows.append(
            {
                "evaluation_kind": EvaluationKind.ACCESSORY_AUDIT,
                # ⛔ NOT truncated: a slice here would let two long keys collide on the unique
                # constraint, and the loser would be refused with a constraint name for a reason.
                "metric_name": key,
                "numeric_value": numeric,
                "verdict": None,
                "detail": detail,
                "source_artifact_path": str(source_path),
            }
        )
    return rows


# ── the run itself ─────────────────────────────────────────────────────────────────────────────
def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def catalogue_key_for(species_key: str, model_key: str, override: str | None = None) -> str:
    """`{species_key}-{model_key}`, or `override` — the string a reader addresses this catalogue by.

    ⛔ **The separator is a HYPHEN, and the guard below is the whole reason.** Catalogue keys must be
    PREFIX-FREE: nuna's payload lookup globs `locus_browser_<key>_*.json` and takes the newest match,
    so a key that is a prefix of another key's payload name silently steals it — `ecoli` matches a
    `locus_browser_ecoli_nuna5_*.json` and, being the later export, wins. The live page then renders
    a different model's catalogue with nothing in any output saying so. No species token contains a
    hyphen, so `ecoli-nuna5` cannot be captured by `ecoli` the way `ecoli_nuna5` can.

    This is the same shape as `-excl` being a prefix of `-exclLOGP` in the run_id tokens, which four
    readers got wrong in one day. Refuse it here rather than write a comment asking people not to.
    """
    key = override or f"{species_key}-{model_key}"
    if not key or any(c.isspace() for c in key) or "/" in key:
        raise PangenomeIngestError(f"catalogue key {key!r} is not URL-safe; it goes in a path segment.")

    # ⛔ NO UNDERSCORE, ANYWHERE IN THE KEY. This is what makes prefix-freeness STRUCTURAL rather
    # than a comparison against whatever keys happen to exist today. nuna's payload lookup globs
    # `locus_browser_<key>_*.json`, so key A captures key B's payloads exactly when B starts with
    # `A + "_"` — which requires an underscore inside B. Forbid the character and no key can ever be
    # another's prefix, whatever is added later and in whatever order.
    #
    # ⚠ An earlier version of this guard only refused `{species_key}_…`. That checked one key against
    # one species key, never against the SET, and so still admitted `ecoli-nuna5` beside
    # `ecoli-nuna5_damped` — the precise silent repoint the hyphen was chosen to prevent.
    if "_" in key:
        raise PangenomeIngestError(
            f"catalogue key {key!r} contains an underscore. nuna's payload lookup globs "
            f"`locus_browser_<key>_*.json` and takes the newest match, so a key with an underscore "
            f"can be captured by a shorter one and silently serve another model's catalogue. Use "
            f"hyphens: {key.replace('_', '-')!r}"
        )

    # ⛔ AND IT MUST CONTAIN A HYPHEN, because that is the invariant the RESOLVER reads. `_resolve_
    # pangenome` treats a key with a hyphen as a catalogue and one without as a species, so a
    # hyphenless key is either unreachable (`/catalogues/sensitive` → "no species 'sensitive'") or,
    # worse, silently resolves to the species DEFAULT instead of the catalogue that carries it
    # (`--catalogue-key ecoli`). A guard that asserts a property the resolver depends on, without
    # enforcing it, is the same shape as a run_id token nothing checks.
    if "-" not in key:
        raise PangenomeIngestError(
            f"catalogue key {key!r} has no hyphen. The API tells a catalogue key from a species key "
            f"by the hyphen, so this one would resolve as a species — returning the species' default "
            f"catalogue, or nothing at all, rather than this one. Use `{species_key}-{key}`."
        )
    return key


def ingest_pangenome_run(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    pathogen_species_id: int,
    genome_collection_id: int,
    nuna_model_id: int | None,
    registry_steps: list[NunaModelStep],
    genome_count: int,
    gene_count: int,
    locus_count: int,
    catalogue_key: str,
    ingest_generation: int = 1,
) -> tuple[int, PangenomeRunReport]:
    """Record the run and everything that describes it; return `(pangenome_id, report)`.

    ⚠ `locus_count` is passed rather than counted here, because the loci do not exist yet when this
    runs — and it is reconciled against `COUNT(*)` once they do, by the catalogue loader.
    """
    from nuna.tl.cluster import run_manifest
    from nuna.tl.probe.accessory_audit_run import model_provenance, parse_model_id

    guard_input_keys_have_not_drifted(run_manifest.INPUT_KEYS)
    run_id = artifacts.run_id
    assignment = artifacts.assignment
    # ⭐ The whole chain, resolved into THIS machine's mirror — see `mirrored_manifest_chain`. A
    # manifest names its parent by the absolute path of the machine that wrote it, so without this
    # the walk stops at one manifest and a four-step model is recorded as a one-step one, silently.
    with manifests_resolved_against(artifacts.processed_root):
        manifest = run_manifest.read_manifest(assignment)
        chain = run_manifest.load_run_chain(assignment)
        provenance = [
            list(pair)
            for pair in model_provenance(run_id, parse_model_id(assignment), assign=assignment)
        ]
    if not chain:
        raise PangenomeIngestError(
            f"no run manifest is reachable from {assignment}. A run_id cannot describe a multi-step "
            "model — tokens are emitted non-default-only — so ingesting without one would record a "
            "model this run may not have been."
        )

    form, form_source, token = resolve_exclusivity_form(run_id, manifest)

    existing = session.execute(
        select(Pangenome).where(
            Pangenome.run_id == run_id, Pangenome.ingest_generation == ingest_generation
        )
    ).scalar_one_or_none()
    if existing is None:
        existing = Pangenome(run_id=run_id, ingest_generation=ingest_generation)
        session.add(existing)

    existing.catalogue_key = catalogue_key
    # ⭐ The model, as one free-text string the footer publishes. `artifacts.model_label` is nuna's
    # own `--model-label` — the same string `guard_audit_describes_this_model` checks the accessory
    # audit against below, so a row and the audit that graded it cannot name different models.
    existing.model_label = artifacts.model_label
    existing.pathogen_species_id = pathogen_species_id
    existing.genome_collection_id = genome_collection_id
    existing.nuna_model_id = nuna_model_id
    existing.exclusivity_form = form
    existing.exclusivity_form_source = form_source
    existing.run_id_exclusivity_token = token
    existing.genome_count = genome_count
    existing.gene_count = gene_count
    existing.locus_count = locus_count
    existing.assignment_file_path = str(assignment)
    existing.assignment_sha256 = _sha256(assignment)
    existing.run_manifest_json = manifest
    existing.chain_manifest_json = {
        "depth": len(chain),
        "runs": [
            {"assign": entry.get("_assign"), "depth": entry.get("_depth"), "run_id": entry.get("run_id")}
            for entry in chain
        ],
    }
    existing.provenance_rows = provenance
    # ⛔ WHICH VERSION of the science package this row was read from. Left NULL it loses exactly
    # the fact the column exists to carry, and a registry is code: `nuna5` meant a different chain
    # before 2026-08-24, so a row that cannot say which nuna it was read from cannot be checked.
    existing.git_sha = nuna_git_sha()
    existing.ingested_at = datetime.now(UTC)
    session.flush()
    pangenome_id = existing.pangenome_id

    # ── steps ──────────────────────────────────────────────────────────────────────────────────
    step_rows = build_step_rows(chain)
    disagreements = cross_check_steps(step_rows, registry_steps)
    session.query(PangenomeStep).filter(PangenomeStep.pangenome_id == pangenome_id).delete(
        synchronize_session=False
    )
    for row in step_rows:
        stage_manifest = row.pop("_manifest")
        session.add(
            PangenomeStep(
                pangenome_id=pangenome_id,
                stage_name=run_manifest.stage_name(stage_manifest)[:128],
                stage_detail=run_manifest.stage_detail(stage_manifest),
                **row,
            )
        )

    # ── lineage ────────────────────────────────────────────────────────────────────────────────
    session.query(PangenomeInputEdge).filter(
        PangenomeInputEdge.child_pangenome_id == pangenome_id
    ).delete(synchronize_session=False)
    edges, unresolved = 0, []
    # ⛔ ONE EDGE PER (input_key, parent), NOT ONE PER CHAIN STAGE — the unique constraint is
    # `(child_pangenome_id, input_key, parent_run_id)` and a longer chain makes collisions ordinary.
    # In a five-step model both step 3c and step 4 consume step 2's assignment, so the naive loop
    # emitted `step2_assign -> <same step2 run>` twice and the whole ingest died on a UniqueViolation
    # part-way through. It never fired on nuna4 only because that chain is one stage shorter.
    #
    # ⚠ Deduplicating loses nothing: the per-stage detail is already kept verbatim in
    # `chain_manifest_json`. This table answers "what was this catalogue built from", which is a
    # property of the pair, not of each stage that read it.
    #
    # The SHALLOWEST depth wins, because depth counts backwards along the chain (deeper is earlier in
    # the pipeline) — so the smallest value is the nearest stage that consumed this parent, which is
    # the honest answer to "how directly did this catalogue depend on it".
    nearest: dict[tuple[str, str], int | None] = {}
    for entry in chain:
        for key in run_manifest.INPUT_KEYS:
            parent_path = entry.get(key)
            if not parent_path:
                continue
            pair = (key, Path(parent_path).stem)
            depth = entry.get("_depth")
            seen = nearest.get(pair, "absent")
            if seen == "absent" or (depth is not None and (seen is None or depth < seen)):
                nearest[pair] = depth

    for (key, parent_run_id), depth in nearest.items():
        parent_id = session.execute(
            select(Pangenome.pangenome_id).where(Pangenome.run_id == parent_run_id)
        ).scalar_one_or_none()
        if parent_id is None:
            unresolved.append(parent_run_id)
        session.add(
            PangenomeInputEdge(
                child_pangenome_id=pangenome_id,
                parent_run_id=parent_run_id[:256],
                parent_pangenome_id=parent_id,
                input_key=key,
                depth=depth,
            )
        )
        edges += 1

    # ── the graded headline ────────────────────────────────────────────────────────────────────
    session.query(PangenomeEvaluation).filter(
        PangenomeEvaluation.pangenome_id == pangenome_id
    ).delete(synchronize_session=False)
    summary = json.loads(artifacts.audit_summary.read_text())
    guard_audit_describes_this_model(summary, artifacts.model_label)
    evaluation_rows = audit_evaluation_rows(summary, artifacts.audit_summary)
    for row in evaluation_rows:
        session.add(PangenomeEvaluation(pangenome_id=pangenome_id, **row))
    session.flush()

    report = PangenomeRunReport(
        run_id=run_id,
        pangenome_id=pangenome_id,
        steps_written=len(step_rows),
        input_edges_written=edges,
        unresolved_parent_run_ids=unresolved,
        evaluation_rows_written=len(evaluation_rows),
        exclusivity_form=form.value,
        exclusivity_form_source=form_source.value,
        step_disagreements=disagreements,
    )
    if len(chain) == 1:
        report.notes.append(
            "the manifest chain is one deep: `precluster_assign` is an absolute CSD3 path, so the "
            "parent runs are recorded as lineage edges with no pangenome behind them"
        )
    # ⛔ No `copy_structure_gate` rows: the gate ledger lives on CSD3 and is not in the local mirror.
    # Said out loud, because an absent evaluation kind and an unmeasured one look identical in a
    # table that only holds what it was given.
    report.notes.append(
        "no copy_structure_gate rows: $P/analysis/ledgers/step4/gate_ledger.tsv is not in the local "
        "mirror. G1-G4 never reach a page, so this blocks nothing but the publish decision"
    )
    return pangenome_id, report


def nuna_git_sha() -> str | None:
    """The **nuna** repository's commit — resolved from nuna's own directory, not the cwd.

    ⛔⛔ **`export_payload._git_sha()` cannot be reused here, and reusing it was a live bug for the
    length of one commit.** It runs `git rev-parse` with no `cwd`, which is correct where it lives
    (the exporter runs inside the nuna checkout) and wrong from here: the ingest runs in the
    *bacatlas* tree, so it returned **bacatlas's** HEAD and wrote it into a column whose docstring
    says *"the version of the registry this row was read from"*. The value was a real short sha, the
    right length and the right shape — and about the wrong repository.
    ⚠ Anchored on `nuna.__file__` rather than on a configured path, so it follows the package that
    was actually imported and cannot drift from it.
    """
    import subprocess
    from pathlib import Path

    import nuna

    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=Path(nuna.__file__).resolve().parent,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def guard_input_keys_have_not_drifted(input_keys) -> None:
    """Refuse to load lineage if `run_manifest.INPUT_KEYS` no longer matches the CHECK constraint.

    Caught here rather than at the INSERT, where it surfaces as a constraint violation naming a
    constraint and not the fact that the pipeline grew a new kind of parent edge.
    """
    if tuple(input_keys) != KNOWN_INPUT_KEYS:
        raise PangenomeIngestError(
            f"nuna's run_manifest.INPUT_KEYS is now {tuple(input_keys)!r} but "
            f"pangenome_input_edge's CHECK admits only {KNOWN_INPUT_KEYS!r}. Add a migration "
            "widening the constraint before loading a run whose chain uses the new key."
        )


def guard_audit_describes_this_model(audit_summary: dict, model_label: str) -> None:
    """Refuse an audit summary that grades a different model — `_check_audit_describes_this_model`'s rule.

    ⛔ A hard error, not a warning. The audit supplies the collapse tier, the Pfam verdict and the
    graded headline; a summary belonging to another model would populate all three with numbers that
    are individually plausible and jointly describe nothing that ran.
    """
    labelled = audit_summary.get("label")
    if labelled and str(labelled) != model_label:
        raise PangenomeIngestError(
            f"the audit summary grades {labelled!r} but this catalogue is {model_label!r}. Its "
            "collapse tiers, Pfam verdicts and headline would all be another model's."
        )
