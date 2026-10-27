"""A pangenome = one `nuna_model` applied to one `genome_collection`. It is what owns loci."""

from __future__ import annotations

from sqlalchemy import (
    BigInteger,
    Boolean,
    CheckConstraint,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    SmallInteger,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY, JSONB
from sqlalchemy.orm import Mapped, mapped_column, relationship

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import (
    EmbeddingRepresentation,
    EvaluationKind,
    ExclusivityForm,
    ExclusivityFormSource,
    GateVerdict,
    RhoRule,
)


class Pangenome(Base):
    """One run of one model over one collection."""

    __tablename__ = "pangenome"
    __table_args__ = (
        UniqueConstraint("run_id", "ingest_generation"),
        # `uq_pangenome__catalogue_key` — 27 chars, well inside Postgres' 63-byte identifier limit.
        # A name over that is truncated with a hash while Alembic compares against the name it ASKED
        # for, so the table reports as drifted from its own migration forever.
        UniqueConstraint("catalogue_key"),
        # ⛔ THE `-excl` PREFIX HAZARD, made unrepresentable. `-excl` is a prefix of `-exclLOGP`,
        # and four readers got this wrong silently in one day. A row that disagrees with its own
        # token cannot exist, so no query downstream has to be careful.
        CheckConstraint(
            "run_id_exclusivity_token IS NULL"
            " OR (run_id_exclusivity_token = '-exclLOGP' AND exclusivity_form = 'EXCLUSION')"
            " OR (run_id_exclusivity_token = '-excl' AND exclusivity_form = 'DAMPED_EXCLUSION')",
            name="exclusivity_token_agrees",
        ),
    )

    pangenome_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)

    #: The on-disk assignment stem.
    run_id: Mapped[str] = mapped_column(String(256), nullable=False)

    #: ⭐ **How a reader addresses this catalogue** — short, URL-safe and stable: `ecoli-nuna4`,
    #: `ecoli-nuna5`, `kp-nuna5`. It is the same key nuna's exporter writes as `--dset` into
    #: `meta.dset`, so one vocabulary spans the static pages and this service.
    #:
    #: ⛔ **The keys must be PREFIX-FREE, and the separator is a HYPHEN for that reason.** nuna's
    #: payload lookup globs `locus_browser_<key>_*.json` and takes the newest match, so a key that is
    #: a prefix of another's payload name silently steals it — `ecoli` would match a
    #: `locus_browser_ecoli_nuna5_*.json` and, being the later export, win. No species token contains
    #: a hyphen, so `ecoli-nuna5` is safe where `ecoli_nuna5` is not. Same trap as `-excl` being a
    #: prefix of `-exclLOGP`, in a second vocabulary.
    #:
    #: ⚠ A species key alone cannot serve here: one species holds several catalogues.
    catalogue_key: Mapped[str] = mapped_column(String(64), nullable=False)

    #: ⭐ A re-ingest of the same run_id bumps this rather than mutating rows in place. 889k loci
    #: with 49M arrangements cannot be diffed row-by-row cheaply, and a half-applied upsert leaves a
    #: catalogue that is partly two models — the worst possible state. Build alongside, verify, flip
    #: the species pointer, then drop the old generation.
    ingest_generation: Mapped[int] = mapped_column(Integer, nullable=False, default=1)

    pathogen_species_id: Mapped[int] = mapped_column(
        ForeignKey("pathogen_species.pathogen_species_id"), nullable=False
    )
    genome_collection_id: Mapped[int] = mapped_column(
        ForeignKey("genome_collection.genome_collection_id"), nullable=False
    )
    nuna_model_id: Mapped[int | None] = mapped_column(ForeignKey("nuna_model.nuna_model_id"), nullable=True)

    exclusivity_form: Mapped[ExclusivityForm] = mapped_column(nullable=False)
    exclusivity_form_source: Mapped[ExclusivityFormSource] = mapped_column(nullable=False)
    #: The literal token found in the run_id, or NULL. Stored so the decision is auditable without
    #: re-parsing — and so the CHECK above has something to check against.
    run_id_exclusivity_token: Mapped[str | None] = mapped_column(String(16), nullable=True)

    #: ⭐ **The model, in full, as one free-text string** — `CPM_gamma0.5_hypergeo_gamma0.6`,
    #: `nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL`. David, 2026-10-27: the footer publishes it
    #: and *"not explained in any more detail than that. Just allows cross reference."*
    #:
    #: ⛔ **This replaces a whole axis, and that is the point.** The model was 32 rows in
    #: `nuna_model` + 129 in `nuna_model_step`, an enum, a registry ingest and nine decoded
    #: provenance rows — all to print a footer nobody picked a catalogue by. The cohort is the axis
    #: now; the model is a label beside `run_id`, which is what actually addresses the artifacts.
    #:
    #: ⚠ It is nuna's own `export_payload --model-label`, already stored as the
    #: `pangenome_evaluation` row with `metric_name='label'` — which is where `meta.model_label` has
    #: always been served from, NOT from `nuna_model.label`. Promoted to a column here so the footer
    #: does not reach through an evaluation row to find the name of its own model.
    model_label: Mapped[str] = mapped_column(String(256), nullable=False)

    genome_count: Mapped[int] = mapped_column(Integer, nullable=False)
    gene_count: Mapped[int] = mapped_column(BigInteger, nullable=False)
    locus_count: Mapped[int] = mapped_column(Integer, nullable=False)

    assignment_file_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    #: ⚠ The assignment TSV is byte-identical-pinned by a test upstream. A different sha for the
    #: same run_id means something changed that should not have, so ingest refuses by default.
    assignment_sha256: Mapped[str | None] = mapped_column(String(64), nullable=True)

    #: `runs/{run_id}.json` VERBATIM. Kept because tokens are emitted non-default-only, so the id
    #: cannot describe the run and the manifest is the only thing that can.
    run_manifest_json: Mapped[dict | None] = mapped_column(JSONB, nullable=True)
    chain_manifest_json: Mapped[dict | None] = mapped_column(JSONB, nullable=True)
    #: `model_provenance`'s `[field, value]` rows, as the footer prints them.
    provenance_rows: Mapped[list | None] = mapped_column(JSONB, nullable=True)
    #: `{section: reason}` for any input that was absent. The page prints these rather than showing
    #: an empty panel, so an absent section is never mistaken for a measured zero.
    omitted_sections: Mapped[dict | None] = mapped_column(JSONB, nullable=True)

    #: ⛔⛔ **Which GO slim produced this catalogue's folded columns — so a change to it is
    #: DETECTABLE rather than silent.** `locus.propagation_verdict_go` and the three
    #: `locus.go_claiming_member_count_*` are computed by closing each gene's terms upward through
    #: the slim (`instruments/annotation_transfer.go_rung`). Before they were stored, a vendored-table
    #: update moved every one of those verdicts on the next cold process with no re-ingest; stored,
    #: it does not, and nothing in the schema would say so.
    #:
    #: ⚠ **The fold changes the ANSWER, not merely its precision.** With no fold loaded `go_rung`
    #: pools all three facets into the first one, which is a different comparison, not a coarser one:
    #: the module's own example is `gumC`, where 16 genes saying *plasma membrane* and one saying
    #: *membrane* read as *"classes differ"* until each side was closed through the slim. So a NULL
    #: here means *the fold was never recorded*, and a stored verdict with no recorded fold cannot be
    #: trusted to mean what a reader assumes.
    annotation_fold_version: Mapped[str | None] = mapped_column(String(64), nullable=True)

    #: ⭐ **The locus the page opens on, and the four chips beside it.** Stored, because neither is
    #: derivable from `locus.interest_score`: `render_page.landing` FILTERS the ranking for the
    #: first locus whose modal Bakta symbol is `fimA` (`LANDING_SYMBOL`) — *"a stranger arriving
    #: cold gets a locus they may recognise rather than whichever one scored highest that week"* —
    #: and falls back to the ranking's own first only where a species has none. On the published
    #: catalogues that preference actually fires: ecoli opens on 2811 and kp on 1098, neither of
    #: which is its `examples[0]`.
    #: ⚠ They were a RENDER-TIME payload mutation, which is how every site page deployed
    #: 2026-08-20..2026-08-25 came to open on locus 0 with no reference names — a caller serialised
    #: the payload before `render` had mutated it, and nothing failed. As columns written by ingest
    #: and checked by the verification gate, that failure is not expressible.
    landing_locus_id: Mapped[int | None] = mapped_column(BigInteger, nullable=True)
    example_locus_ids: Mapped[list[int] | None] = mapped_column(ARRAY(BigInteger), nullable=True)

    git_sha: Mapped[str | None] = mapped_column(String(40), nullable=True)
    built_at: Mapped[str | None] = mapped_column(String(32), nullable=True)
    ingested_at: Mapped[object | None] = mapped_column(DateTime(timezone=True), nullable=True)
    #: ⛔ **Visibility, NOT "the one we serve".** True means *offer this catalogue in the picker*; a
    #: species may have several true at once. Which one answers a bare species request is
    #: `pathogen_species.default_pangenome_id`. The two were the same fact only while each species
    #: had exactly one catalogue.
    is_published: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    steps: Mapped[list[PangenomeStep]] = relationship(
        back_populates="pangenome", order_by="PangenomeStep.step_ordinal", cascade="all, delete-orphan"
    )

    def __repr__(self) -> str:
        return f"<Pangenome {self.run_id} gen{self.ingest_generation}>"


class PangenomeStep(Base):
    """What a step of THIS run actually did — the same shape as the model's step, plus its cost.

    Cross-checked against `nuna_model_step` at ingest. ⚠ A disagreement is **recorded, not
    resolved**: that is how a run whose manifest contradicts the registry gets caught rather than
    silently normalised to one of them.
    """

    __tablename__ = "pangenome_step"
    __table_args__ = (UniqueConstraint("pangenome_id", "step_ordinal"),)

    pangenome_step_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )

    step_ordinal: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    step_name: Mapped[str] = mapped_column(String(32), nullable=False)
    representation: Mapped[EmbeddingRepresentation | None] = mapped_column(nullable=True)
    gamma: Mapped[float | None] = mapped_column(Float, nullable=True)
    rho_rule: Mapped[RhoRule | None] = mapped_column(nullable=True)
    rho_ceiling: Mapped[float | None] = mapped_column(Float, nullable=True)
    node_mode: Mapped[str | None] = mapped_column(String(32), nullable=True)
    uses_exclusivity: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    step_assignment_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    gene_count: Mapped[int | None] = mapped_column(BigInteger, nullable=True)
    node_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    wall_seconds: Mapped[float | None] = mapped_column(Float, nullable=True)
    #: ⚠ NULL where the chain manifest recorded `n/a`. Not measured is not zero.
    peak_rss_gigabytes: Mapped[float | None] = mapped_column(Float, nullable=True)

    #: The manifest's own strings, verbatim.
    stage_name: Mapped[str | None] = mapped_column(String(128), nullable=True)
    stage_detail: Mapped[str | None] = mapped_column(Text, nullable=True)

    pangenome: Mapped[Pangenome] = relationship(back_populates="steps")


class PangenomeInputEdge(Base):
    """run → parent-run lineage, from `run_manifest.INPUT_KEYS`.

    ⚠ `parent_pangenome_id` is NULL where the parent run was never ingested. That is a gap
    **admitted**, not invented: `parent_run_id` still records what the manifest named.
    """

    __tablename__ = "pangenome_input_edge"
    __table_args__ = (
        UniqueConstraint("child_pangenome_id", "input_key", "parent_run_id"),
        CheckConstraint(
            "input_key IN ('precluster_assign', 'atoms_assign', 'step2_assign')",
            name="input_key_is_known",
        ),
    )

    pangenome_input_edge_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    child_pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )
    parent_run_id: Mapped[str] = mapped_column(String(256), nullable=False)
    parent_pangenome_id: Mapped[int | None] = mapped_column(
        ForeignKey("pangenome.pangenome_id"), nullable=True
    )
    input_key: Mapped[str] = mapped_column(String(32), nullable=False)
    depth: Mapped[int | None] = mapped_column(SmallInteger, nullable=True)


class PangenomeEvaluation(Base):
    """How a pangenome graded — gate, accessory audit, and (reserved) merge evidence.

    ⭐ One table with an `evaluation_kind`, so the spurious-merge assessment lands as a third kind
    rather than a fourth table.

    ⛔ `verdict` is three-valued and NOT NULL. `gates.py` returns `passed=None` for a metric a run
    could not produce, and *"an unmeasured gate must not masquerade as a pass or a fail"*. As a
    nullable boolean, `bool_and` over SQL NULLs silently converts *incomplete* into *not counted*.
    ⚠ This matters for choosing which pangenome to publish and never reaches a page: the browser
    reads the audit headline, not the gate — G1–G4 appear nowhere in `app.js`.
    """

    __tablename__ = "pangenome_evaluation"
    __table_args__ = (
        UniqueConstraint("pangenome_id", "evaluation_kind", "metric_name"),
        Index("ix_pangenome_evaluation__kind_verdict", "evaluation_kind", "verdict"),
        CheckConstraint(
            "numeric_value IS NULL OR numeric_value <> 'NaN'::double precision",
            name="numeric_value_is_not_nan",
        ),
    )

    pangenome_evaluation_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )
    evaluation_kind: Mapped[EvaluationKind] = mapped_column(nullable=False)
    #: `G1`…`G4`, `GATE`, `synteny_only_gene_rate`, `single_copy_core`, …
    #: ⛔ **`String(128)`, and the ingest does NOT truncate.** The audit's keys are compositional
    #: (`{family}_{measure}_{qualifier}_{sub}`) and the longest today is 52 characters
    #: (`family_split_across_clusters_n_split_genes_singleton`) against a declared 64 — 81 % of the
    #: width. That is close enough that a longer key is a matter of one more qualifier, and the
    #: failure would not be a truncation alone: `metric_name` is part of the UNIQUE constraint, so
    #: two long keys sharing a 64-character prefix would collide and one of them would be refused
    #: — or, with an earlier `[:64]` slice, silently overwrite the other.
    metric_name: Mapped[str] = mapped_column(String(128), nullable=False)
    numeric_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    verdict: Mapped[GateVerdict | None] = mapped_column(nullable=True)
    detail: Mapped[str | None] = mapped_column(Text, nullable=True)
    source_artifact_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
