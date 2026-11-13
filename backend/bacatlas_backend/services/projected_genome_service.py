"""Genomes placed on a pangenome they were never part of — read side.

⛔ **These rows are drawn BESIDE the model, never inside it.** Nothing here reaches a count the page
already quotes: prevalence, bands, the census, the marginals, arrangement genome counts and the
genome picker are all the modelled 100, whatever is projected. Every function in this module reads
`projected_*` and nothing else.

⚠ **A projected genome is not an anchor genome, and the two must never share a code path.** An
anchor is a member — its genes are in `genome_locus_assignment` and its id is inside
`locus_arrangement.member_genome_ids`. A projected genome is in neither, so a lookup that searched
the collection would find it absent and report *"this genome has no gene at this locus"*, which is a
claim about a genome the pangenome says nothing about.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.projected_genome import ProjectedGenePlacement, ProjectedGenome


@dataclass(frozen=True)
class ProjectedGenomeRow:
    """One projected genome as the picker and the summary need it."""

    sample_id: str
    genome_id: int
    rule_label: str
    neighbours_searched: int
    neighbours_reported: int
    gene_count: int
    placed_gene_count: int
    gene_without_vector_count: int
    gene_without_window_count: int
    contested_gene_count: int
    distinct_locus_count: int
    multi_copy_locus_count: int
    window_matched_gene_count: int
    nearest_cosine_median: float | None
    nearest_cosine_fifth_percentile: float | None
    nearest_cosine_minimum: float | None
    nuna_git_sha: str | None


def list_projected_genomes(session: Session, pangenome_id: int) -> list[ProjectedGenomeRow]:
    """Every genome projected onto this pangenome, with its counts. **One statement.**"""
    rows = session.execute(
        select(ProjectedGenome, Genome.sample_id)
        .join(Genome, Genome.genome_id == ProjectedGenome.genome_id)
        .where(ProjectedGenome.pangenome_id == pangenome_id)
        .order_by(Genome.sample_id)
    ).all()
    return [
        ProjectedGenomeRow(
            sample_id=str(sample_id),
            genome_id=int(projected.genome_id),
            rule_label=projected.rule_label,
            neighbours_searched=projected.neighbours_searched,
            neighbours_reported=projected.neighbours_reported,
            gene_count=projected.gene_count,
            placed_gene_count=projected.placed_gene_count,
            gene_without_vector_count=projected.gene_without_vector_count,
            gene_without_window_count=projected.gene_without_window_count,
            contested_gene_count=projected.contested_gene_count,
            distinct_locus_count=projected.distinct_locus_count,
            multi_copy_locus_count=projected.multi_copy_locus_count,
            window_matched_gene_count=projected.window_matched_gene_count,
            nearest_cosine_median=projected.nearest_cosine_median,
            nearest_cosine_fifth_percentile=projected.nearest_cosine_fifth_percentile,
            nearest_cosine_minimum=projected.nearest_cosine_minimum,
            nuna_git_sha=projected.nuna_git_sha,
        )
        for projected, sample_id in rows
    ]


def resolve_projected_genome(session: Session, pangenome_id: int, sample_id: str) -> int | None:
    """A genome PROJECTED onto this pangenome, by sample id — or `None`. One statement.

    ⛔ Scoped to `projected_genome`, not to the species: the database holds 22 *E. coli* and 58 Kp
    genomes that are in neither the collection nor any projection, and resolving one of those would
    produce a page describing a genome nothing was ever computed for.
    """
    return session.execute(
        select(ProjectedGenome.genome_id)
        .join(Genome, Genome.genome_id == ProjectedGenome.genome_id)
        .where(ProjectedGenome.pangenome_id == pangenome_id, Genome.sample_id == sample_id)
    ).scalar_one_or_none()


def load_placements_at_locus(
    session: Session, *, pangenome_id: int, genome_id: int, locus_id: int
) -> list[ProjectedGenePlacement]:
    """Every gene of a projected genome placed at one locus, in copy order. **One statement.**

    ⚠ A **list**, for the same reason `anchor_arrangement_ranks` is: several genes of one genome can
    land on one locus, as the model's own loci have ρ > 1, and returning the first would silently
    drop the rest.
    """
    return list(
        session.execute(
            select(ProjectedGenePlacement)
            .where(
                ProjectedGenePlacement.pangenome_id == pangenome_id,
                ProjectedGenePlacement.genome_id == genome_id,
                ProjectedGenePlacement.locus_id == locus_id,
            )
            .order_by(ProjectedGenePlacement.copy_ordinal)
        ).scalars()
    )
