"""A genome collection is the roster a pangenome was built over — the "which genomes" question.

⚠ Two incompatible roster lineages exist and they do not share a schema. Naming the lineage on the
collection is what stops a loader assuming either one.
"""

from __future__ import annotations

from sqlalchemy import (
    Boolean,
    ForeignKey,
    Index,
    Integer,
    String,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import RosterLineage


class GenomeCollection(Base):
    """`probe_ecoli`, `kp_HI_s0`, `union_roster`, `blood_faeces`, `kpsc_all`."""

    __tablename__ = "genome_collection"
    __table_args__ = (UniqueConstraint("pathogen_species_id", "collection_key"),)

    genome_collection_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    pathogen_species_id: Mapped[int] = mapped_column(
        ForeignKey("pathogen_species.pathogen_species_id"), nullable=False
    )
    collection_key: Mapped[str] = mapped_column(String(64), nullable=False)

    #: ⭐ **What the picker SHOWS for this cohort** — "E. coli probe", "AMR panel", "Kp". The first
    #: display text this table has ever carried, and the reason it now has a read path at all: the
    #: browser's axis is the cohort, not the model (David, 2026-10-27).
    #:
    #: ⛔ **It must NOT contain the genome count.** The page renders `"{label} — {count:,} genomes"`
    #: from this and :attr:`genome_count` beside it, which is *"passed, never inferred"*. A label
    #: that spelled its own count could disagree with the roster it names, and nothing would say so
    #: — the example labels David gave ("blood/faeces 12,150 genomes") are the rendered form, not
    #: the stored one.
    collection_label: Mapped[str] = mapped_column(String(128), nullable=False)

    roster_lineage: Mapped[RosterLineage] = mapped_column(nullable=False)

    #: ⭐ N — the prevalence denominator AND the exclusivity denominator. **Passed, never inferred**
    #: from the rows present, which is `profile_from_pairs`' own rule: a genome that contributed no
    #: gene to a locus still counts against it.
    genome_count: Mapped[int] = mapped_column(Integer, nullable=False)

    roster_source_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    roster_source_sha256: Mapped[str | None] = mapped_column(String(64), nullable=True)

    #: Nested rosters — `kp_HI_anchor` ⊂ `kp_HI_s0` ⊂ `union_roster`. This is what makes the
    #: measured genome→locus curve unconfounded, so the nesting is a fact worth storing.
    parent_genome_collection_id: Mapped[int | None] = mapped_column(
        ForeignKey("genome_collection.genome_collection_id"), nullable=True
    )


class GenomeCollectionMembership(Base):
    """Which genomes are in a collection, and **at which ordinal**."""

    __tablename__ = "genome_collection_membership"
    __table_args__ = (
        UniqueConstraint(
            "genome_collection_id", "collection_genome_ordinal"
        ),
        Index("ix_genome_collection_membership__genome_id", "genome_id"),
    )

    genome_collection_id: Mapped[int] = mapped_column(
        ForeignKey("genome_collection.genome_collection_id", ondelete="CASCADE"), primary_key=True
    )
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id"), primary_key=True)

    #: ⭐ The position in the sorted accession list that `meta.genomes` publishes and `arr.gid`
    #: indexes into. **Stored, never re-derived by sorting** — it is a contract with already
    #: published payloads, and re-deriving it is exactly the implicit ordering this project punishes.
    collection_genome_ordinal: Mapped[int] = mapped_column(Integer, nullable=False)

    is_anchor: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    assembly_file_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    annotation_file_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    #: `long_read` | `short_read`. ⚠ Never mixed within a genome — `resolve_file_pair`'s rule.
    file_path_source: Mapped[str | None] = mapped_column(String(16), nullable=True)

    #: What the roster ASKED for, which may differ from `genome.sample_id` only when the join went
    #: through `bare_accession`. Stored with `matched_by` so a bare-accession match is visible.
    requested_sample_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    matched_by: Mapped[str | None] = mapped_column(String(32), nullable=True)


class GenomeCollectionBuildReport(Base):
    """The four ways a genome fails to reach a run, each counted and named.

    ⛔ **This is not decoration.** "13,602 samples silently becoming 13,171 runnable" is the failure
    the roster builder exists to prevent, and a database that stored only the survivors would
    reintroduce it exactly. ⚠ Every count is NULLABLE and **NULL means the store was not given**,
    not zero — a coverage of 0 and a coverage never measured are different findings.
    """

    __tablename__ = "genome_collection_build_report"

    genome_collection_id: Mapped[int] = mapped_column(
        ForeignKey("genome_collection.genome_collection_id", ondelete="CASCADE"), primary_key=True
    )

    requested_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    matched_exact_key_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    matched_bare_accession_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    not_in_metadata_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    no_usable_file_pair_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    path_present_file_absent_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    no_embedding_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    roster_count: Mapped[int | None] = mapped_column(Integer, nullable=True)

    coverage_meta_parquet: Mapped[int | None] = mapped_column(Integer, nullable=True)
    coverage_esm_embedding: Mapped[int | None] = mapped_column(Integer, nullable=True)
    coverage_bacformer_embedding: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: A few examples, so a count is actionable rather than just alarming.
    not_in_metadata_examples: Mapped[list[str] | None] = mapped_column(ARRAY(String), nullable=True)
