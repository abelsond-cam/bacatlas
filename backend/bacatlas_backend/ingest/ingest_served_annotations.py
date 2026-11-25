"""Load one catalogue's annotation-transfer predictions into `served_locus_annotation`.

The file, `{catalogue}_served_annotations.tsv`, is nuna's `summarise_transfer_gains.py` output: one row
per (node, claim) that GAINS a prediction, already cut at the version's serve bar, with the donor it was
copied from and both medoid cosines (`nuna/docs/model_evaluation/annotation_transfer_gains.md` §8).

This stage is `ingest --stage served-annotations`, and it is deliberately NOT part of `all`: a donor may
sit in ANOTHER catalogue (on the probe, 23,556 of Kp's 25,475 predictions come from *E. coli* donors), so
both catalogues must already be loaded, and a catalogue build must not wait on another catalogue's load.

⭐ **The catalogue's state is always stated, never inferred** (`ServedAnnotationState`). The stage takes
exactly one of three explicit instructions: load THIS file (`--served-annotations PATH`), record that the
model produced none (`--no-served-annotations`), or clear what is loaded (`--clear-served-annotations`).
⛔ A path nobody gave is never read as "none produced" (review of 9cfa481: an unmounted data root, or a
set key the producer never writes, used to delete a catalogue's predictions and say none were made).

⛔ **Refused, by name and with every problem at once, rather than loaded wrong:**

* no `annotation_model_version` column, an empty value, two values, or a version this build does not
  register — a prediction whose model cannot be named cannot be shown with its accuracy;
* no `catalogue_run_id` column, or one that is not this catalogue's `run_id` — ⭐ the IDENTITY check. The
  gene-count check below is only a size-histogram check: labels are assigned near size order, so a
  relabelled run of the same shape passes it (shifting every label by one still matches 97 % of rows);
* a row below the version's serve cut, an accuracy or cosine out of range, a repeated (node, claim), a
  claim the version does not predict;
* a COG class that is not one category letter, or is `S` (*function unknown* — the producer strips it);
* a GO class that is not a SLIM class (a namespace root, or a raw term), is in another namespace than its
  claim's, or a set that is not one slim class closed upward through the reference's own slim;
* a node label that does not resolve, or resolves to a locus of a different SIZE than the file says;
* a donor catalogue that is not loaded — ⚠ the whole file, because the producer kept ONE donor per node
  (the higher-probability direction), so the rows from an absent catalogue cannot be dropped and the rest
  served: the losing calls for those nodes are not in the file.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from sqlalchemy import select, text
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.staging_table_loader import copy_rows
from bacatlas_backend.instruments.gene_ontology_fold import NAMESPACE_ROOTS, GeneOntologyFold
from bacatlas_backend.models.enumerations import ServedAnnotationState
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.reference_vocabulary import GeneOntologyTerm
from bacatlas_backend.models.served_locus_annotation import (
    CLAIM_ADDRESS,
    COSINE_BOUND,
    KNOWN_ANNOTATION_MODEL_VERSIONS,
    AnnotationModelVersion,
    ServedLocusAnnotation,
)

#: ⛔ The column that names the model. The probe files of CSD3 job 39455187 (written 2026-10-28) predate the
#: annotation model's versioning (2026-11-14) and do NOT carry it; nuna's `summarise_transfer_gains.py` must
#: write it. Until then a probe file is refused, by design.
VERSION_COLUMN = "annotation_model_version"

#: ⛔ The column that names the RUN the file's node labels belong to — the assignment the medoids were built
#: on — which must equal the catalogue's `run_id`. Requested of nuna's producer (B5) on 2026-11-14; the
#: probe files predate it and are refused as written.
CATALOGUE_RUN_COLUMN = "catalogue_run_id"

#: ⛔ The run each row's DONOR label belongs to — equal to `catalogue_run_id` for a same-catalogue donor.
#: Every distinct value must be a loaded catalogue's `run_id`; the producer takes both from the medoid
#: sidecars' `assignment` path, extension stripped, which is exactly what `pangenome.run_id` holds.
DONOR_RUN_COLUMN = "donor_run_id"

#: Every column this stage reads. ⚠ The producer writes more (`rank_esm`, `rank_bac`, `direction`,
#: `both_directions`, `directions_agree`); they are diagnostics, and `served_locus_annotation`'s module
#: docstring says why each is not kept.
REQUIRED_COLUMNS = (
    "recipient_node",
    "claim",
    "predicted_classes",
    "n_predicted_classes",
    "p_serve",
    "donor_node",
    "esm_medoid_sim",
    "bac_medoid_sim",
    "n_recipient_genes",
    "n_donor_genes",
    "donor_catalogue",
    VERSION_COLUMN,
    CATALOGUE_RUN_COLUMN,
    DONOR_RUN_COLUMN,
)

#: ⛔ The producer's set separator (`nuna.tl.annotation_transfer.derive_node_claims.CLASS_SEP`).
CLASS_SEP = "|"

#: COG category `S` is *function unknown*; the producer drops it from every claim (`COG_UNKNOWN`).
COG_UNKNOWN = "S"

#: How many examples a refusal names per problem — enough to find the pattern, short enough to read.
EXAMPLES = 5


class ServedAnnotationsRefused(ValueError):
    """The file cannot be loaded as it stands. The message lists EVERY problem found, each with its count."""


class ServedDonorsBlockReload(RuntimeError):
    """Another catalogue's predictions take donors from the catalogue being reloaded — refused, by name."""


@dataclass
class ServedAnnotationReport:
    """What the stage did, in counts that can be reconciled against the file and the table."""

    catalogue_key: str
    state: ServedAnnotationState
    path: Path | None = None
    version: str | None = None
    sha256: str | None = None
    rows_in_file: int = 0
    rows_written: int = 0
    rows_replaced: int = 0
    by_claim: dict[str, int] = field(default_factory=dict)
    by_donor_catalogue: dict[str, int] = field(default_factory=dict)

    def render(self) -> str:
        """The report as the lines the CLI prints."""
        if self.state is ServedAnnotationState.NONE_PRODUCED:
            return (
                f"served annotations for {self.catalogue_key}: recorded NONE PRODUCED, as instructed "
                f"(--no-served-annotations); {self.rows_replaced:,} earlier rows removed"
            )
        if self.state is ServedAnnotationState.NOT_LOADED:
            return (
                f"served annotations for {self.catalogue_key}: CLEARED — now not loaded; "
                f"{self.rows_replaced:,} rows removed"
            )
        claims = " · ".join(f"{claim} {count:,}" for claim, count in sorted(self.by_claim.items()))
        donors = " · ".join(f"{key} {count:,}" for key, count in sorted(self.by_donor_catalogue.items()))
        return "\n".join(
            [
                f"served annotations for {self.catalogue_key}: {self.version} from {self.path}",
                f"  file sha256       {self.sha256}",
                f"  rows in file      {self.rows_in_file:,}",
                f"  rows written      {self.rows_written:,}  ({claims})",
                f"  donor catalogue   {donors}",
                f"  earlier rows      {self.rows_replaced:,} replaced",
            ]
        )


def _examples(values) -> str:
    listed = list(values)
    shown = ", ".join(str(value) for value in listed[:EXAMPLES])
    return f"{shown}{' …' if len(listed) > EXAMPLES else ''}"


def _one_constant(frame: pd.DataFrame, column: str, what: str) -> tuple[str | None, list[str]]:
    """The single value a constant column holds, or the problems that stop it being one."""
    values = frame[column].str.strip()
    absent = int((values == "").sum())
    named = sorted(set(values[values != ""]))
    problems = []
    if absent:
        problems.append(f"{absent:,} row(s) have an empty {column}")
    if len(named) > 1:
        problems.append(f"the file names {len(named)} {what}s {named}; one file is one {what}")
    return (named[0] if len(named) == 1 else None), problems


def read_served_annotations(path: Path) -> tuple[pd.DataFrame, AnnotationModelVersion, str]:
    """Read the file and name its model and its run — or refuse. Returns `(frame, version, run_id)`.

    ⛔ Node labels are read as TEXT. They look numeric and are not: a label read as an integer joins to
    `locus.node_label` as a different string, or not at all (`a-numeric-looking-text-key-coerces-silently`).
    """
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "recipient_node": str,
            "donor_node": str,
            "claim": str,
            "predicted_classes": str,
            "donor_catalogue": str,
            VERSION_COLUMN: str,
            CATALOGUE_RUN_COLUMN: str,
            DONOR_RUN_COLUMN: str,
        },
        keep_default_na=False,
        na_values={"p_serve": [""], "esm_medoid_sim": [""], "bac_medoid_sim": [""]},
    )
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        hints = []
        if VERSION_COLUMN in missing:
            hints.append(
                "⛔ No `annotation_model_version` column, so the model that made these predictions cannot be "
                "named — and a prediction shown with an accuracy must say whose accuracy it is. nuna's "
                "`summarise_transfer_gains.py` must write it (one constant column, e.g. `at-v0.2`); the probe "
                "files of CSD3 job 39455187 predate the annotation model's versioning and carry none."
            )
        if CATALOGUE_RUN_COLUMN in missing or DONOR_RUN_COLUMN in missing:
            hints.append(
                "⛔ No `catalogue_run_id` / `donor_run_id`, so nothing says which run the node labels belong to — "
                "and a label resolves in any run. nuna's producer writes both (the assignment run each medoid "
                "file was built on) from CSD3 job 40993415 on; the job 39455187 files carry neither."
            )
        raise ServedAnnotationsRefused(f"{path}: missing column(s) {missing}. " + " ".join(hints))
    problems = []
    if not len(frame):
        problems.append("the file has no rows, so it names no annotation model version and no run")
    tag, version_problems = _one_constant(frame, VERSION_COLUMN, "version")
    run_id, run_problems = _one_constant(frame, CATALOGUE_RUN_COLUMN, "run")
    problems += version_problems + run_problems
    if tag is not None and tag not in KNOWN_ANNOTATION_MODEL_VERSIONS:
        problems.append(
            f"version {tag!r} is not registered in KNOWN_ANNOTATION_MODEL_VERSIONS "
            f"(known: {sorted(KNOWN_ANNOTATION_MODEL_VERSIONS)}) — a new version is a code change that states "
            "its cut and its claims, never something a file introduces by spelling a tag"
        )
    if problems:
        raise ServedAnnotationsRefused(f"{path}: " + "; ".join(problems))
    return frame, KNOWN_ANNOTATION_MODEL_VERSIONS[tag], run_id


def file_sha256(path: Path) -> str:
    """The file's sha256 — its provenance, and the load fingerprint the API's validators carry."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gene_ontology_slim_classes(session: Session) -> frozenset[str]:
    """The slim's own classes, less the three namespace roots — the only classes a GO prediction may name.

    ⚠ Read from `gene_ontology_term` (a slim class is the one term whose fold is itself), NOT from the fold:
    `GeneOntologyFold.slim_of` maps a raw term that reaches no slim class to itself too, so there the two
    look the same.
    """
    rows = session.execute(
        select(GeneOntologyTerm.go_id, GeneOntologyTerm.slim_go_ids)
    ).all()
    return frozenset(go_id for go_id, slim in rows if slim is not None and list(slim) == [go_id]) - NAMESPACE_ROOTS


def order_lineage(classes: frozenset[str], fold: GeneOntologyFold) -> list[str] | None:
    """A closed GO lineage, most specific class first — or `None` if the set is not one.

    The producer's GO claim is the donor's voted slim class closed upward (`derive_node_claims.go_claims`),
    so exactly one class has every other as an ancestor, and the set is exactly that class plus its slim
    ancestors in the reference. `None` where either fails, which the caller refuses.

    ⚠ The ancestors follow in an order where each class precedes its own ancestors — sorted by how many
    slim ancestors each has, descending. That is a valid ORDERING (an ancestor's ancestor set is a strict
    subset of its descendant's) and nothing more: ⛔ it is NOT a depth, which in a DAG is the longest path,
    not a closure size (`closure-size-is-not-depth`). Nothing here reads it as one.
    """
    leaves = [go_class for go_class in classes if classes - {go_class} <= fold.ancestors.get(go_class, frozenset())]
    if len(leaves) != 1:
        return None
    leaf = leaves[0]
    if classes != {leaf} | fold.ancestors.get(leaf, frozenset()):
        return None
    above = sorted(classes - {leaf}, key=lambda go_class: (-len(fold.ancestors.get(go_class, ())), go_class))
    return [leaf, *above]


def _resolve_loci(session: Session, pangenome_id: int) -> pd.DataFrame:
    """`node_label → (locus_id, member_gene_count)` for one pangenome — one statement."""
    rows = session.execute(
        select(Locus.node_label, Locus.locus_id, Locus.member_gene_count).where(Locus.pangenome_id == pangenome_id)
    ).all()
    return pd.DataFrame(rows, columns=["node_label", "locus_id", "member_gene_count"]).set_index("node_label")


def newest_pangenome(session: Session, run_id: str) -> Pangenome | None:
    """The newest ingest generation of a run — the one every other stage addresses by `run_id`."""
    return session.execute(
        select(Pangenome).where(Pangenome.run_id == run_id).order_by(Pangenome.ingest_generation.desc()).limit(1)
    ).scalar_one_or_none()


def _delete_existing(session: Session, pangenome_id: int) -> int:
    """Remove this catalogue's predictions, whatever version wrote them — a load replaces, never appends."""
    result = session.execute(
        text(
            "DELETE FROM served_locus_annotation s USING locus l "
            "WHERE s.locus_id = l.locus_id AND l.pangenome_id = :pangenome_id"
        ),
        {"pangenome_id": pangenome_id},
    )
    return int(result.rowcount or 0)


def _set_state(
    pangenome: Pangenome, state: ServedAnnotationState, version: str | None = None, sha256: str | None = None
) -> None:
    """The three columns move together — the pangenome's CHECKs refuse any other combination."""
    pangenome.served_annotation_state = state
    pangenome.annotation_model_version = version
    pangenome.served_annotation_sha256 = sha256


def record_none_produced(session: Session, *, pangenome: Pangenome) -> ServedAnnotationReport:
    """⭐ The operator states the model produced NONE for this catalogue (`--no-served-annotations`)."""
    replaced = _delete_existing(session, pangenome.pangenome_id)
    _set_state(pangenome, ServedAnnotationState.NONE_PRODUCED)
    session.flush()
    return ServedAnnotationReport(
        catalogue_key=pangenome.catalogue_key, state=ServedAnnotationState.NONE_PRODUCED, rows_replaced=replaced
    )


def clear_served_annotations(session: Session, *, pangenome: Pangenome) -> ServedAnnotationReport:
    """Remove this catalogue's predictions and say they are NOT LOADED (`--clear-served-annotations`).

    ⭐ The command a blocked reload of ANOTHER catalogue names: it removes the rows that hold donors there.
    """
    replaced = _delete_existing(session, pangenome.pangenome_id)
    _set_state(pangenome, ServedAnnotationState.NOT_LOADED)
    session.flush()
    return ServedAnnotationReport(
        catalogue_key=pangenome.catalogue_key, state=ServedAnnotationState.NOT_LOADED, rows_replaced=replaced
    )


def reset_served_state(session: Session, pangenome_id: int) -> None:
    """⛔ A catalogue reload's half of the bargain: its rows cascade away with its loci, so it is NOT LOADED.

    Called by `_delete_pangenome_layer`. Without it a reloaded catalogue kept `at-v0.1` over zero rows and
    every locus read "declined" (review of 9cfa481).
    """
    session.execute(
        text(
            "UPDATE pangenome SET served_annotation_state = 'NOT_LOADED', annotation_model_version = NULL, "
            "served_annotation_sha256 = NULL WHERE pangenome_id = :pangenome_id"
        ),
        {"pangenome_id": pangenome_id},
    )


def clear_command(pangenome: Pangenome) -> str:
    """The exact command that clears a catalogue's predictions — for a refusal to name."""
    set_key = pangenome.model_label.split("_", 1)[0]
    return (
        "python -m bacatlas_backend.ingest --stage served-annotations --clear-served-annotations "
        f"--set-key {set_key} --model-label {pangenome.model_label} --run-id {pangenome.run_id}"
    )


def refuse_reload_while_donating(session: Session, pangenome_id: int) -> None:
    """⛔ Refuse to reload a catalogue while ANOTHER catalogue's predictions take donors from it.

    The donor key is deferred to commit, so a reload holding both ends passes; but a row in another
    catalogue naming a donor here would fail at COMMIT — after the frame build, which is hours at 71k. So
    this runs FIRST, and names each such catalogue, its row count and the command that clears it.
    """
    holders = session.execute(
        text(
            "SELECT p.pangenome_id, count(*) FROM served_locus_annotation s "
            "JOIN locus d ON d.locus_id = s.donor_locus_id "
            "JOIN locus r ON r.locus_id = s.locus_id "
            "JOIN pangenome p ON p.pangenome_id = r.pangenome_id "
            "WHERE d.pangenome_id = :pangenome_id AND r.pangenome_id <> :pangenome_id "
            "GROUP BY p.pangenome_id ORDER BY p.pangenome_id"
        ),
        {"pangenome_id": pangenome_id},
    ).all()
    if not holders:
        return
    lines = []
    for holder_id, count in holders:
        holder = session.get(Pangenome, holder_id)
        lines.append(f"{holder.catalogue_key}: {count:,} prediction(s) — clear with: {clear_command(holder)}")
    raise ServedDonorsBlockReload(
        "this catalogue donates to annotation-transfer predictions held by another catalogue, and reloading it "
        "would delete those donors. Nothing was built. Clear them first, reload, then load them again:\n  "
        + "\n  ".join(lines)
    )


def load_served_annotations(
    session: Session,
    *,
    pangenome: Pangenome,
    set_key: str,
    path: Path,
    fold: GeneOntologyFold,
    slim_classes: frozenset[str],
) -> ServedAnnotationReport:
    """Validate the whole file, resolve both ends of every row, then replace this catalogue's predictions.

    ``set_key`` is the dataset token the producer writes as `donor_catalogue` for a donor in THIS catalogue
    (`ecoli`, `kp`); it is checked against the run columns, never used to find a run. ⭐ Each row's donor is
    resolved through its OWN `donor_run_id` — the producer's statement of which run the donor label belongs
    to — so no operator flag supplies it, and none can disagree with the file. ``slim_classes`` is
    `gene_ontology_slim_classes(session)` — the only GO classes a prediction may name.
    """
    frame, version, catalogue_run_id = read_served_annotations(path)
    problems: list[str] = []
    if catalogue_run_id != pangenome.run_id:
        problems.append(
            f"⛔ the file's {CATALOGUE_RUN_COLUMN} is {catalogue_run_id!r}, and this catalogue is run "
            f"{pangenome.run_id!r}: its node labels belong to another run, where every one resolves to a "
            "different locus"
        )

    # ── the rows themselves ───────────────────────────────────────────────────────────────────────
    unknown_claims = sorted(set(frame["claim"]) - version.claims)
    if unknown_claims:
        problems.append(f"claim(s) {unknown_claims} are not predicted by {version.tag} ({sorted(version.claims)})")
    p = pd.to_numeric(frame["p_serve"], errors="coerce")
    bad_p = frame.loc[~(p >= version.serve_cut) | ~(p <= 1.0), ["recipient_node", "claim", "p_serve"]]
    if len(bad_p):
        problems.append(
            f"{len(bad_p):,} row(s) have p_serve outside [{version.serve_cut}, 1] — {version.tag} serves only at "
            f"p ≥ {version.serve_cut}: {_examples(bad_p.itertuples(index=False, name=None))}"
        )
    for column in ("esm_medoid_sim", "bac_medoid_sim"):
        values = pd.to_numeric(frame[column], errors="coerce")
        bad = frame.loc[~(values.abs() <= COSINE_BOUND), "recipient_node"]
        if len(bad):
            problems.append(f"{len(bad):,} row(s) have a {column} that is not a cosine: nodes {_examples(bad)}")
    repeated = frame.loc[frame.duplicated(["recipient_node", "claim"], keep=False), ["recipient_node", "claim"]]
    if len(repeated):
        problems.append(
            f"{len(repeated):,} row(s) repeat a (node, claim) — one prediction per claim: "
            f"{_examples(repeated.drop_duplicates().itertuples(index=False, name=None))}"
        )

    # ── the class sets: split once per DISTINCT string, which is a few hundred, not one per row ──
    needs_reference = bool(set(frame["claim"]) - {"cog"})
    if needs_reference and not fold.namespace_of:
        raise ServedAnnotationsRefused(
            f"{path}: the GO reference is not loaded (`gene_ontology_term` is empty), so no GO class can be "
            "checked against its namespace or ordered within its lineage. Run `ingest --stage reference` first."
        )
    ordered: dict[tuple[str, str], list[str] | str] = {}
    for claim, classes in frame[["claim", "predicted_classes"]].drop_duplicates().itertuples(index=False):
        members = [member.strip() for member in str(classes).split(CLASS_SEP)]
        if not classes or any(not member for member in members) or len(set(members)) != len(members):
            ordered[(claim, classes)] = "an empty or repeated class"
            continue
        if claim == "cog":
            if any(len(member) != 1 or not member.isalpha() or not member.isupper() for member in members):
                ordered[(claim, classes)] = "a COG class that is not one category letter"
            elif COG_UNKNOWN in members:
                ordered[(claim, classes)] = "COG `S` (function unknown), which is not a prediction"
            else:
                ordered[(claim, classes)] = sorted(members)
            continue
        if claim not in CLAIM_ADDRESS:
            continue  # already refused above, as an unknown claim
        _kind, namespace = CLAIM_ADDRESS[claim]
        if not set(members) <= slim_classes:
            ordered[(claim, classes)] = "a GO class that is not a slim class (a namespace root, or a raw GO term)"
            continue
        if any(fold.namespace_of.get(member) != namespace for member in members):
            ordered[(claim, classes)] = f"a GO class the reference does not place in namespace {namespace}"
            continue
        lineage = order_lineage(frozenset(members), fold)
        ordered[(claim, classes)] = (
            lineage if lineage is not None else "not one slim class closed upward through the reference's slim"
        )
    refused_sets = {key: reason for key, reason in ordered.items() if isinstance(reason, str)}
    for reason, count in Counter(refused_sets.values()).items():
        sets = [f"{claim}:{classes}" for (claim, classes), why in refused_sets.items() if why == reason]
        problems.append(f"{count:,} distinct class set(s) hold {reason}: {_examples(sets)}")
    counted = frame["predicted_classes"].map(lambda value: len(str(value).split(CLASS_SEP)))
    miscounted = frame.loc[counted != pd.to_numeric(frame["n_predicted_classes"], errors="coerce"), "recipient_node"]
    if len(miscounted):
        problems.append(
            f"{len(miscounted):,} row(s) have n_predicted_classes ≠ the classes listed: {_examples(miscounted)}"
        )

    # ── both ends of every row, and the guard that they are THIS run ─────────────────────────────
    recipients = _resolve_loci(session, pangenome.pangenome_id)
    donor_runs = frame[DONOR_RUN_COLUMN].str.strip()
    if (donor_runs == "").any():
        problems.append(f"{int((donor_runs == '').sum()):,} row(s) have an empty {DONOR_RUN_COLUMN}")
    # ⛔ The token and the run must tell one story: a same-catalogue donor is THIS run, another catalogue's
    # is not, and one token names one run. Either contradiction means the columns were assembled wrongly.
    here = frame["donor_catalogue"] == set_key
    contradicts = (here & (donor_runs != catalogue_run_id)) | (~here & (donor_runs == catalogue_run_id))
    if contradicts.any():
        pairs = frame.loc[contradicts, ["donor_catalogue", DONOR_RUN_COLUMN]].drop_duplicates()
        problems.append(
            f"{int(contradicts.sum()):,} row(s) have a donor_catalogue that contradicts their {DONOR_RUN_COLUMN} "
            f"(a donor in {set_key!r} must be run {catalogue_run_id!r}, and only it): "
            f"{_examples(pairs.itertuples(index=False, name=None))}"
        )
    runs_per_token = frame.assign(_run=donor_runs).groupby("donor_catalogue")["_run"].nunique()
    for token in runs_per_token[runs_per_token > 1].index:
        problems.append(f"donor catalogue {token!r} names {int(runs_per_token[token])} different donor runs")
    donor_pangenomes: dict[str, int] = {}
    for run_id in sorted(set(donor_runs) - {""}):
        rows_there = int((donor_runs == run_id).sum())
        donor_pangenome = pangenome if run_id == pangenome.run_id else newest_pangenome(session, run_id)
        if donor_pangenome is None:
            tokens = sorted(set(frame.loc[donor_runs == run_id, "donor_catalogue"]))
            problems.append(
                f"{rows_there:,} of {len(frame):,} row(s) take their donor from catalogue {tokens} (run {run_id}), "
                "which is NOT LOADED in this database. ⚠ The file cannot be loaded without it: the producer kept "
                "one donor per node, so those rows cannot be dropped and the rest served"
            )
            continue
        donor_pangenomes[run_id] = donor_pangenome.pangenome_id
    donor_loci = {run_id: _resolve_loci(session, pid) for run_id, pid in donor_pangenomes.items()}

    recipient_id = frame["recipient_node"].map(recipients["locus_id"])
    recipient_size = frame["recipient_node"].map(recipients["member_gene_count"])
    unresolved = frame.loc[recipient_id.isna(), "recipient_node"]
    if len(unresolved):
        problems.append(
            f"{len(unresolved):,} recipient node(s) are not loci of {pangenome.catalogue_key}: {_examples(unresolved)}"
        )
    wrong_size = frame.loc[
        recipient_id.notna() & (recipient_size != pd.to_numeric(frame["n_recipient_genes"], errors="coerce")),
        ["recipient_node", "n_recipient_genes"],
    ]
    if len(wrong_size):
        problems.append(
            f"⛔ {len(wrong_size):,} recipient node(s) resolve to a locus of a DIFFERENT size than the file says "
            f"— the file was built on another run: {_examples(wrong_size.itertuples(index=False, name=None))}"
        )
    donor_id = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    donor_size = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    for run_id, loci in donor_loci.items():
        rows = donor_runs == run_id
        donor_id[rows] = frame.loc[rows, "donor_node"].map(loci["locus_id"]).astype("Int64")
        donor_size[rows] = frame.loc[rows, "donor_node"].map(loci["member_gene_count"]).astype("Int64")
    resolvable = donor_runs.isin(list(donor_loci))
    unresolved = frame.loc[resolvable & donor_id.isna(), ["donor_catalogue", "donor_node"]]
    if len(unresolved):
        problems.append(
            f"{len(unresolved):,} donor node(s) are not loci of their catalogue: "
            f"{_examples(unresolved.itertuples(index=False, name=None))}"
        )
    wrong_size = frame.loc[
        donor_id.notna() & (donor_size != pd.to_numeric(frame["n_donor_genes"], errors="coerce")),
        ["donor_catalogue", "donor_node", "n_donor_genes"],
    ]
    if len(wrong_size):
        problems.append(
            f"⛔ {len(wrong_size):,} donor node(s) resolve to a locus of a DIFFERENT size than the file says — "
            f"the donor catalogue is another run: {_examples(wrong_size.itertuples(index=False, name=None))}"
        )

    if problems:
        raise ServedAnnotationsRefused(
            f"{path}: {len(problems)} problem(s), nothing loaded:\n  " + "\n  ".join(problems)
        )

    # ── write: replace this catalogue's predictions whole ─────────────────────────────────────────
    sha256 = file_sha256(path)
    replaced = _delete_existing(session, pangenome.pangenome_id)

    def rows():
        for index, row in enumerate(frame.itertuples(index=False)):
            kind, namespace = CLAIM_ADDRESS[row.claim]
            yield (
                int(recipient_id.iat[index]),
                kind.name,
                namespace,
                ordered[(row.claim, row.predicted_classes)],
                float(row.p_serve),
                int(donor_id.iat[index]),
                float(row.esm_medoid_sim),
                float(row.bac_medoid_sim),
                version.tag,
            )

    written = copy_rows(
        session,
        ServedLocusAnnotation.__table__,
        (
            "locus_id",
            "annotation_kind",
            "gene_ontology_namespace",
            "term_values",
            "estimated_accuracy",
            "donor_locus_id",
            "esm_medoid_similarity",
            "bacformer_medoid_similarity",
            "annotation_model_version",
        ),
        rows(),
    )
    if written != len(frame):
        raise ServedAnnotationsRefused(f"{path}: wrote {written:,} rows of {len(frame):,}")
    _set_state(pangenome, ServedAnnotationState.LOADED, version=version.tag, sha256=sha256)
    session.flush()
    return ServedAnnotationReport(
        catalogue_key=pangenome.catalogue_key,
        state=ServedAnnotationState.LOADED,
        path=path,
        version=version.tag,
        sha256=sha256,
        rows_in_file=len(frame),
        rows_written=written,
        rows_replaced=replaced,
        by_claim=dict(Counter(frame["claim"])),
        by_donor_catalogue=dict(Counter(frame["donor_catalogue"])),
    )
