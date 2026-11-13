"""A genome the model never clustered, placed on a published pangenome's loci.

⛔ **This layer only ever ADDS.** A projected genome is not a member of the collection, contributes no
gene to `gene_locus_membership`, and changes no count anywhere else: prevalence, bands, the census,
the marginals, arrangement genome counts and the genome picker all stay the modelled 100. Two tables,
written by `--stage projection`, and nothing else in the schema is touched — which is what keeps the
publish gate, the column audit and the parity suites meaningful while this feature exists.

⛔ **A placement is not what the model would have done.** The model is a CPM chain (ESM step 2 →
Bacformer 3b → step 4); this is one nearest-neighbour hop in Bacformer space. `rule_label` carries that
sentence into the database so no reader of these rows has to be told separately, and the API and the
page repeat it. Method and limitations: nuna `docs/locus_projection.md`.

⚠ **There is no novelty status** (David, 2026-09-22). Every gene is placed, `locus_id` is NOT NULL, and
a distant placement is visible through `nearest_cosine` and the agreement rather than hidden behind a
threshold. A schema with a nullable locus would invite one.
"""

from __future__ import annotations

from sqlalchemy import (
    BigInteger,
    Boolean,
    CheckConstraint,
    DateTime,
    ForeignKey,
    ForeignKeyConstraint,
    Index,
    Integer,
    SmallInteger,
    String,
    Text,
    func,
)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.column_types import measurement, nan_guards


class ProjectedGenome(Base):
    """One genome placed on one pangenome, with the rule it was placed by and what came of it."""

    __tablename__ = "projected_genome"
    __table_args__ = (
        Index("ix_projected_genome__pangenome_id", "pangenome_id"),
        CheckConstraint("placed_gene_count >= 0", name="placed_gene_count_is_not_negative"),
        # ⚠ The two are recorded separately because they are different decisions: how many neighbours
        # were SEARCHED (and stored, so a larger check can be tried without re-searching) and how many
        # were used to CHECK the placement. Reading one as the other would misreport the evidence.
        CheckConstraint(
            "neighbours_reported <= neighbours_searched",
            name="neighbours_reported_within_neighbours_searched",
        ),
        *nan_guards("nearest_cosine_median", "nearest_cosine_fifth_percentile", "nearest_cosine_minimum"),
    )

    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), primary_key=True
    )
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)

    #: The sentence a reader is owed, stored beside the rows rather than only in a doc.
    rule_label: Mapped[str] = mapped_column(Text, nullable=False)
    #: k, as searched (100) — every one is kept on the cluster, so 5/25/50/100 can be tried later.
    neighbours_searched: Mapped[int] = mapped_column(Integer, nullable=False)
    #: n, as used for the agreement (10 today).
    neighbours_reported: Mapped[int] = mapped_column(Integer, nullable=False)
    representation: Mapped[str] = mapped_column(String(32), nullable=False, default="bacformer")

    #: Provenance. ⚠ `assignment_sha256` must equal `pangenome.assignment_sha256` — a projection
    #: computed against a different assignment would name loci that mean something else, and the
    #: ingest refuses it by name rather than letting the labels resolve by accident.
    source_file_path: Mapped[str] = mapped_column(Text, nullable=False)
    source_sha256: Mapped[str | None] = mapped_column(String(64), nullable=True)
    assignment_sha256: Mapped[str | None] = mapped_column(String(64), nullable=True)
    nuna_git_sha: Mapped[str | None] = mapped_column(String(64), nullable=True)
    ingested_at: Mapped[object] = mapped_column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    #: ⛔ Four counts that are NOT the same number, and merging any two makes the page claim something
    #: it cannot support: how many genes the genome has; how many were placed; how many have no
    #: Bacformer vector at all (past its 6,000-protein cut) and so no placement; and how many are
    #: alone on their contig, which have a placement but NO ±5 window — an absence of observation, not
    #: a neighbourhood that matched nothing.
    gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    placed_gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    gene_without_vector_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    gene_without_window_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    contested_gene_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    distinct_locus_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    multi_copy_locus_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    window_matched_gene_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    #: The distribution, not a mean: with no novelty threshold this is the only thing that makes a
    #: distant placement visible on the page.
    nearest_cosine_median: Mapped[float | None] = measurement()
    nearest_cosine_fifth_percentile: Mapped[float | None] = measurement()
    nearest_cosine_minimum: Mapped[float | None] = measurement()

    def __repr__(self) -> str:
        return f"<ProjectedGenome pangenome={self.pangenome_id} genome={self.genome_id}>"


class ProjectedGenePlacement(Base):
    """One gene of a projected genome, on the locus of its nearest modelled gene.

    ⚠ **`locus_id` is NOT NULL.** Every gene is placed (David: no threshold, no novel status), so
    there is no row here that means "nowhere". A gene with no Bacformer vector has no row at all and is
    counted in `projected_genome.gene_without_vector_count` instead — absent, not null.
    """

    __tablename__ = "projected_gene_placement"
    __table_args__ = (
        # ⛔ No foreign key onto the gene's coordinates any more (migration `9e8f7a6b5c4d`, 2026-11-13).
        # The gene layer is one row per genome — `genome_genes.gene_spans`, a block indexed by
        # `flat_index` — so there is no `(genome_id, flat_index)` row for a constraint to point at.
        # The guarantee it gave (every placed gene is one of THIS genome's genes) now lives in
        # `ingest_projected_genomes._check_the_genes_are_this_genomes`, which refuses a `flat_index`
        # outside `[0, genome_genes.gene_count)`. A check in a loader is not a constraint in a schema —
        # `storing_the_gene_layer.md` §6 item 3 records that cost as chosen, not overlooked.
        ForeignKeyConstraint(
            ["pangenome_id", "genome_id"],
            ["projected_genome.pangenome_id", "projected_genome.genome_id"],
            ondelete="CASCADE",
            name="fk_projected_gene_placement__projected_genome",
        ),
        # ⭐ The hot path: "what has this projected genome placed at this locus", resolved in one
        # index scan while the arrangements load.
        Index("ix_projected_gene_placement__pangenome_id__locus_id", "pangenome_id", "locus_id"),
        Index("ix_projected_gene_placement__locus_id__genome_id", "locus_id", "genome_id"),
        CheckConstraint(
            "agreeing_neighbour_count >= 0",
            name="agreeing_neighbour_count_is_not_negative",
        ),
        *nan_guards("nearest_cosine", "placed_summed_cosine", "runner_up_summed_cosine"),
    )

    pangenome_id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    #: ⚠ Integer, matching `genome.genome_id`; `pangenome_id` is BigInteger, matching its own.
    genome_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    flat_index: Mapped[int] = mapped_column(Integer, primary_key=True)

    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)
    nearest_cosine: Mapped[float | None] = measurement()

    #: ⛔ **Read `agreeing_neighbour_count` against `available_neighbour_count`, never alone.** A locus
    #: with *m* modelled genes can supply at most min(n, m) of the n checkers, so the raw count is
    #: bounded by the locus's size and not only by the evidence. Measured on the first ten genomes:
    #: contested is 0.39 % at loci with ≥ 10 modelled genes and 42.6 % below that, 96.3 % at
    #: singletons. The denominator is stored so a serialiser cannot forget it.
    agreeing_neighbour_count: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    available_neighbour_count: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    placed_summed_cosine: Mapped[float | None] = measurement()

    #: The runner-up is ABSENT, not zero, when all n checkers sit in the placed locus — a real and
    #: good case that a 0.0 would misreport as a beaten rival.
    runner_up_locus_id: Mapped[int | None] = mapped_column(
        ForeignKey("locus.locus_id", ondelete="SET NULL"), nullable=True
    )
    runner_up_summed_cosine: Mapped[float | None] = measurement()
    is_contested: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    #: ρ > 1 is ordinary in this model, so copies are counted, not resolved — "copy 2 of 3 here".
    copy_ordinal: Mapped[int] = mapped_column(SmallInteger, nullable=False, default=1)
    copies_at_locus: Mapped[int] = mapped_column(SmallInteger, nullable=False, default=1)

    #: The gene's own ±5 neighbourhood, packed exactly as `locus_arrangement.neighbour_slot_codes` is
    #: (`catalogue_ordinal * 2 + same_strand`, −1 at a contig end) — built by the model's own window
    #: code, so the two are comparable by equality.
    #: ⚠ NULL means the gene is **alone on its contig**: no window exists at all. That is not the same
    #: as a window matching no arrangement, which is `matched_locus_arrangement_id IS NULL` with a
    #: vector present.
    neighbour_slot_codes: Mapped[list[int] | None] = mapped_column(ARRAY(Integer), nullable=True)
    matched_locus_arrangement_id: Mapped[int | None] = mapped_column(
        ForeignKey("locus_arrangement.locus_arrangement_id", ondelete="SET NULL"), nullable=True
    )

    def __repr__(self) -> str:
        return f"<ProjectedGenePlacement genome={self.genome_id} gene={self.flat_index} locus={self.locus_id}>"
