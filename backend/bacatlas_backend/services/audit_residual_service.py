"""The audit's residual loci — two lists of what the 2026-11-26 audit measured, each locus one click away.

⭐ **What they list (David, 2026-11-26).** The footer used to list *every locus grouped on context alone*
and *every Pfam conflict*. At 71,435 genomes both were decided by stray genes — kp-71k locus 4008, Wzi on
69,825 of 69,884 genes, sat in both — and David ruled them audit bugs: the page never says "no homology",
"over-merge", "grouped on context alone" or "conflicting". The two lists now name facts, not verdicts:

1. **Loci grouped by neither MMseqs nor a shared Pfam architecture** — the ladder's last rung, `esm_only`.
   The locus card quotes the ESM number for each.
2. **Loci where at most 95 % of genes share one Pfam architecture** — `pfam_architecture_shared` false,
   among non-singletons with at least one annotated gene. ⚠ A singleton is excluded because it has no
   members to share anything; a locus with no annotated gene because "no gene carries a domain" is a
   coverage fact the card states, not an architecture that fails to be shared.

⛔ **Two lists, never one, and they overlap.** An `esm_only` locus is by construction one whose Pfam
architecture is not shared, unless it carries no Pfam at all. A locus in both appears in both.

⭐ **The predicates live HERE, once.** The shell counts both lists inside its census statement
(`species_catalogue_service`) and this module lists them, from the same two clauses — so a footer summary
can never promise a count its own list does not hold.

⚠ **Fetched when the reader opens the list, not with the species.** At the probe scale both lists are tens
of loci; at the 80,000-genome design target they could be thousands, and a closed `<details>` element is
the wrong place to spend a first page load.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import and_, or_, select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import EmbeddingRepresentation
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_similarity import LocusSimilarity

#: The 2026-11-26 ladder's last rung: MMseqs groups the locus at no identity, and its genes share no
#: Pfam architecture on more than 95 % of them.
ESM_ONLY_TIER = "esm_only"


def esm_only_clause():
    """List 1 — the `esm_only` tier."""
    return Locus.collapse_tier == ESM_ONLY_TIER


def pfam_architecture_not_shared_clause():
    """List 2 — the audit's `pfam_shared` false, on a locus with members to share and Pfam to share.

    ⛔ `IS FALSE`, never `NOT … IS TRUE`: a catalogue audited before 2026-11-26 holds NULL here, which is
    *not measured*, and must list nothing rather than every locus.
    """
    return and_(
        Locus.pfam_architecture_shared.is_(False),
        Locus.member_gene_count > 1,
        Locus.pfam_annotated_member_count >= 1,
    )


@dataclass(frozen=True)
class ResidualLocus:
    """One row of either list — enough evidence beside the name to judge it before clicking."""

    label: str
    display_name: str
    prevalence_band: str
    gene_count: int
    uniref50_family_count: int | None
    pfam_architecture_count: int | None
    #: ⭐ The two facts the lists are about, so a row states its own reason: the genes on the commonest
    #: Pfam domain set (over `gene_count`, ALL the locus's genes) and the ESM number the card quotes.
    pfam_dominant_domain_gene_count: int | None
    esm_member_similarity_p05: float | None
    syntenic_a5: float | None
    #: ⛔ Similarities now, not the `1 − d` distances the medoid card stored. A residual row shows
    #: "ESM 0.98/0.61" either way, but the numbers behind it are a median over every within-locus
    #: gene pair and the highest such median against another locus — not a member's distance to one
    #: gene. The client no longer converts; there is nothing to convert.
    esm_within_similarity: float | None
    esm_nearest_similarity: float | None


@dataclass(frozen=True)
class AuditResiduals:
    """Both lists, each in catalogue order."""

    esm_only: list[ResidualLocus]
    pfam_architecture_not_shared: list[ResidualLocus]


def load_audit_residuals(session: Session, *, pangenome_id: int) -> AuditResiduals:
    """Both residual lists, in ONE statement — the two sets overlap, so each row is read once."""
    esm_only = esm_only_clause().label("is_esm_only")
    not_shared = pfam_architecture_not_shared_clause().label("is_pfam_architecture_not_shared")
    rows = session.execute(
        select(
            Locus.node_label,
            Locus.display_name,
            Locus.prevalence_band,
            Locus.member_gene_count,
            Locus.uniref50_family_count,
            Locus.pfam_architecture_count,
            Locus.pfam_dominant_domain_gene_count,
            Locus.esm_member_similarity_p05,
            Locus.syntenic_a5,
            LocusSimilarity.within_similarity,
            LocusSimilarity.nearest_similarity,
            esm_only,
            not_shared,
        )
        # ⚠ An OUTER join: a locus with no similarity row would otherwise drop out of its own list.
        .outerjoin(
            LocusSimilarity,
            and_(
                LocusSimilarity.locus_id == Locus.locus_id,
                LocusSimilarity.representation == EmbeddingRepresentation.ESM,
            ),
        )
        .where(
            Locus.pangenome_id == pangenome_id,
            or_(esm_only_clause(), pfam_architecture_not_shared_clause()),
        )
        .order_by(Locus.catalogue_ordinal)
    ).all()

    first: list[ResidualLocus] = []
    second: list[ResidualLocus] = []
    for row in rows:
        entry = ResidualLocus(
            label=row.node_label,
            display_name=row.display_name,
            prevalence_band=row.prevalence_band.value,
            gene_count=row.member_gene_count,
            uniref50_family_count=row.uniref50_family_count,
            # ⚠ `-1`/NULL means Pfam could not judge the locus — not "no architectures".
            pfam_architecture_count=(
                row.pfam_architecture_count
                if row.pfam_architecture_count is not None and row.pfam_architecture_count >= 0
                else None
            ),
            pfam_dominant_domain_gene_count=row.pfam_dominant_domain_gene_count,
            esm_member_similarity_p05=row.esm_member_similarity_p05,
            syntenic_a5=row.syntenic_a5,
            esm_within_similarity=row.within_similarity,
            esm_nearest_similarity=row.nearest_similarity,
        )
        # ⛔ Each list is decided by the database's OWN evaluation of its clause, read back as a column —
        # never by re-testing the row in Python, which would be a second spelling of the predicate.
        if row.is_esm_only:
            first.append(entry)
        if row.is_pfam_architecture_not_shared:
            second.append(entry)
    return AuditResiduals(esm_only=first, pfam_architecture_not_shared=second)
