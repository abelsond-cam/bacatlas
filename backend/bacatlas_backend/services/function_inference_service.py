"""The Function tab's inference block — what this node's function rests on, and how sure that is.

Two mechanisms, reported separately because they are not equally trustworthy and the page must never
let them read alike (`nuna/docs/model_evaluation/function_inference.md`):

* **its own genes.** The node's call is already a modal vote over its member genes, and `checkable`
  says whether that vote had more than one voter. ⛔ A node with exactly ONE annotated gene is
  unanimous *by construction*; 518 of 4,993 *E. coli* COG nodes are in that state, `gumC` among them
  (COG3206 on 1 of 100 genes), and a card that showed it like any other would be claiming a check
  that never happened.
* **a neighbour.** Walk the stored ESM neighbours outward until one carries the vocabulary, then
  quote it only as deep as its similarity earns — and carry the **measured** agreement rate for that
  depth so the reader judges the suggestion instead of trusting it.

⚠ **The walk is reported, not just its result.** David, 2026-10-03: *"If they don't have any
annotation, then look to next neighbour etc. Then note the similarity."* A card that printed only the
supplying neighbour would hide that ranks 1-4 were empty — and for *E. coli* COG only 1,679 of 4,143
assignments come from rank 1, so that is the common case rather than the exception.
"""

from __future__ import annotations

from collections import defaultdict

from sqlalchemy import select, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    NOT_CALLED,
    SUPPORTED_KINDS,
    Calibration,
    cog_levels,
    compute_calibration,
    ec_levels,
    quotable_level,
    single_rung,
    tier_for,
)
from bacatlas_backend.models.enumerations import AnnotationKind, EmbeddingRepresentation
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry
from bacatlas_backend.models.locus_similarity import LocusNearestLocus

#: ⚠ A process-lifetime cache, legitimate only because this service is READ-ONLY and a catalogue
#: cannot change without a re-ingest and a restart — the same ground the `ETag: "{pangenome_id}"` on
#: every response already stands on. Measured at **0.19-0.42 s** per (catalogue, vocabulary) on the
#: published catalogues, so the first request after a restart pays that once.
#: ⛔ Tests must call `clear_calibration_cache()` between catalogues, or one fixture's ladder answers
#: another's — the trap `gene_sequence_service.clear_parsed_genome_cache` exists for.
_CALIBRATIONS: dict[tuple[int, AnnotationKind, EmbeddingRepresentation], Calibration] = {}

#: ⛔ **Counted from the gene rows, not read from `locus.*_annotated_member_count`.** Those columns
#: exist for COG, EC and KEGG but NOT for GO: GO has three, one per namespace, and summing them
#: double counts every gene annotated in more than one. One bounded statement over at most six loci
#: (the focal node and its five neighbours) gives all four axes on the same definition.
#: ⚠ Verified against the three columns that do exist: 0 mismatching nodes out of 17,531 and 15,670.
ANNOTATED_COUNT_PREDICATE = {
    AnnotationKind.COG_ORTHOGROUP: "f.cog_id is not null",
    AnnotationKind.GENE_ONTOLOGY_SLIM: "f.gene_ontology_terms is not null",
    AnnotationKind.EC_NUMBER: "f.ec_numbers is not null",
    AnnotationKind.KEGG_ORTHOLOGY: "f.kegg_orthology_id is not null",
}


def _annotated_counts(session: Session, locus_ids: set[int]) -> dict[int, dict[AnnotationKind, int]]:
    """Per locus, how many member genes carry each vocabulary — the denominator of `checkable`."""
    if not locus_ids:
        return {}
    columns = ", ".join(
        f"count(*) filter (where {predicate}) as {kind.name.lower()}"
        for kind, predicate in ANNOTATED_COUNT_PREDICATE.items()
    )
    rows = session.execute(
        text(f"""
            select m.locus_id, {columns}
              from gene_locus_membership m
              left join gene_functional_annotation f
                     on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where m.locus_id = any(:locus_ids)
             group by m.locus_id
        """),
        {"locus_ids": list(locus_ids)},
    ).all()
    return {
        row.locus_id: {kind: getattr(row, kind.name.lower()) for kind in ANNOTATED_COUNT_PREDICATE}
        for row in rows
    }


def clear_calibration_cache() -> None:
    """Drop every cached ladder. Call between catalogues in tests."""
    _CALIBRATIONS.clear()


def calibration_for(
    session: Session,
    *,
    pangenome_id: int,
    annotation_kind: AnnotationKind,
    representation: EmbeddingRepresentation = EmbeddingRepresentation.ESM,
) -> Calibration:
    """The agreement ladder for this catalogue, computed once per process."""
    key = (pangenome_id, annotation_kind, representation)
    if key not in _CALIBRATIONS:
        _CALIBRATIONS[key] = compute_calibration(
            session,
            pangenome_id=pangenome_id,
            annotation_kind=annotation_kind,
            representation=representation,
        )
    return _CALIBRATIONS[key]


def _fold(annotation_kind: AnnotationKind, terms: list[str], categories: list[str] | None):
    """A node's stated terms folded onto the ladder for its vocabulary.

    ⚠ `terms` is a LIST because rank 0 is not unique per locus: `top_go_slim` ranks within
    (locus, namespace), so a GO-bearing node states up to three rank-0 terms and reading one would
    discard the rest.
    """
    if annotation_kind is AnnotationKind.EC_NUMBER:
        return ec_levels(terms[0]) if terms else None
    if annotation_kind is AnnotationKind.COG_ORTHOGROUP:
        return cog_levels(terms[0] if terms else None, categories)
    return single_rung(terms)


def _top_entries(session: Session, locus_ids: set[int]) -> dict[tuple[int, AnnotationKind], list]:
    """Every supported vocabulary's rank-0 entries for every locus, in one statement.

    ⛔ A LIST per (locus, kind), not a single row — see `_fold`.
    """
    if not locus_ids:
        return {}
    rows = session.execute(
        select(
            LocusAnnotationEntry.locus_id,
            LocusAnnotationEntry.annotation_kind,
            LocusAnnotationEntry.term_value,
            LocusAnnotationEntry.term_name,
            LocusAnnotationEntry.member_gene_count,
        ).where(
            LocusAnnotationEntry.locus_id.in_(locus_ids),
            LocusAnnotationEntry.annotation_kind.in_(SUPPORTED_KINDS),
            LocusAnnotationEntry.rank_within_locus == 0,
        )
    ).all()
    grouped: dict[tuple[int, AnnotationKind], list] = defaultdict(list)
    for locus_id, kind, term, name, support in rows:
        grouped[(locus_id, kind)].append((term, name, support))
    return grouped


def load_inference(
    session: Session,
    *,
    pangenome_id: int,
    locus: Locus,
    representation: EmbeddingRepresentation = EmbeddingRepresentation.ESM,
) -> dict:
    """What can be said about this node's function, by which mechanism, and at what measured rate."""
    neighbours = session.execute(
        select(
            LocusNearestLocus.rank,
            LocusNearestLocus.cross_similarity,
            Locus.locus_id,
            Locus.node_label,
            Locus.catalogue_ordinal,
            Locus.display_name,
            Locus.member_gene_count,
            Locus.modal_cog_categories,
        )
        .join(Locus, Locus.locus_id == LocusNearestLocus.neighbour_locus_id)
        .where(
            LocusNearestLocus.locus_id == locus.locus_id,
            LocusNearestLocus.representation == representation,
        )
        .order_by(LocusNearestLocus.rank)
    ).all()

    touched = {locus.locus_id, *(row.locus_id for row in neighbours)}
    entries = _top_entries(session, touched)
    annotated_counts = _annotated_counts(session, touched)
    vocabularies = []
    for annotation_kind in SUPPORTED_KINDS:
        said = entries.get((locus.locus_id, annotation_kind), [])
        own_term, own_name, own_support = said[0] if said else (None, None, None)
        annotated = annotated_counts.get(locus.locus_id, {}).get(annotation_kind, 0)
        own = None
        if own_term is not None:
            own = {
                "term": own_term,
                "name": own_name,
                "gene_count": own_support,
                "annotated_gene_count": annotated,
                "member_gene_count": locus.member_gene_count,
                # ⛔ The whole point of the field: one voter is not a vote.
                "checkable": annotated >= 2,
            }
        vocabularies.append(
            {
                "annotation_kind": annotation_kind.value,
                "own": own,
                # ⚠ Only ever offered where the node has NO call of its own — a suggestion beside a
                # real annotation would compete with it.
                **(
                    _walk(annotation_kind, neighbours, entries, annotated_counts,
                          session, pangenome_id, representation)
                    if own is None
                    else {"walk": [], "candidate": None}
                ),
            }
        )
    return {"representation": representation.value, "vocabularies": vocabularies}


def _walk(
    annotation_kind: AnnotationKind,
    neighbours,
    entries,
    annotated_counts,
    session: Session,
    pangenome_id: int,
    representation: EmbeddingRepresentation,
) -> dict:
    """Outward along the ranks until one neighbour carries the vocabulary; report every step."""
    walk, candidate = [], None
    for row in neighbours:
        said = entries.get((row.locus_id, annotation_kind), [])
        term, name, support = said[0] if said else (None, None, None)
        tier = tier_for(row.cross_similarity)
        walk.append(
            {
                "rank": row.rank,
                "cosine": row.cross_similarity,
                "tier": tier,
                "node_label": row.node_label,
                "catalogue_ordinal": row.catalogue_ordinal,
                "display_name": row.display_name,
                "carries_annotation": term is not None,
            }
        )
        if candidate is not None or term is None:
            continue
        #: ⛔ The first donor found ends the walk even when its tier is `NOT_CALLED`: a nearer
        #: neighbour that carries nothing cannot be skipped over in favour of a FURTHER one that
        #: does, or the cosine the rate is read from would no longer be the walk's own.
        donor_annotated = annotated_counts.get(row.locus_id, {}).get(annotation_kind, 0)
        folded = _fold(annotation_kind, [t for t, _n, _s in said], row.modal_cog_categories)
        level = quotable_level(annotation_kind, tier, folded) if (folded and tier != NOT_CALLED) else None
        cell = (
            calibration_for(
                session,
                pangenome_id=pangenome_id,
                annotation_kind=annotation_kind,
                representation=representation,
            ).cell(tier, level)
            if level is not None
            else None
        )
        candidate = {
            "rank": row.rank,
            "cosine": row.cross_similarity,
            "tier": tier,
            "donor": {
                "node_label": row.node_label,
                "catalogue_ordinal": row.catalogue_ordinal,
                "display_name": row.display_name,
                "term": term,
                "name": name,
                "gene_count": support,
                "annotated_gene_count": donor_annotated,
                "member_gene_count": row.member_gene_count,
                "checkable": donor_annotated >= 2,
            },
            # `null` where the tier is too remote to call, or the donor states nothing this deep.
            "level": level,
            "value": sorted(folded[level]) if level is not None else None,
            "calibration": cell.as_json() if cell is not None else None,
        }
    return {"walk": walk, "candidate": candidate}
