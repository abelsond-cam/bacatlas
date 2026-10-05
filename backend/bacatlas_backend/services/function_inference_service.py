"""The Function tab's inference block — what this node's function rests on, and how sure that is.

Two mechanisms, reported separately because they are not equally trustworthy and the page must never
let them read alike (`nuna/docs/model_evaluation/function_inference.md`):

* **its own genes.** The node's call is already a modal vote over its member genes, and `checkable`
  says whether that vote had more than one voter. ⛔ A node with exactly ONE annotated gene is
  unanimous *by construction*; 518 of 4,993 *E. coli* COG nodes are in that state, `gumC` among them
  (COG3206 on 1 of 100 genes), and a card that showed it like any other would be claiming a check
  that never happened. ⭐ **But that is a missing check, NOT a missing claim** — the call still
  covers the node's other 99 genes, and `propagation` is the measured rate at which it is right to.
  (David, 2026-10-04: *"the whole point of our method is effectively the carrying by one is
  enough."*)
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
    ANNOTATED_COUNT_PREDICATE,
    NOT_CALLED,
    SUPPORTED_KINDS,
    Calibration,
    callable_level,
    cog_levels,
    compute_calibration,
    ec_levels,
    quotable_level,
    single_rung,
    tier_for,
)
from bacatlas_backend.instruments.within_node_propagation import (
    PropagationRate,
    compute_propagation,
)
from bacatlas_backend.models.enumerations import AnnotationKind, EmbeddingRepresentation, PrevalenceBand
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

#: ⚠ The same stance as `_CALIBRATIONS` and the same justification, but a **heavier** scan: it reads
#: every annotated GENE rather than every stored edge, measured at **0.22-1.94 s** per (catalogue,
#: vocabulary) on the published catalogues — so a cold *E. coli* Function tab pays ~5.8 s across the
#: four, once per process. ⛔ Keyed WITHOUT a representation, which is a finding and not an omission:
#: within-node unanimity is 99.8-100 % at every band of BOTH representations (ESM's within-node
#: similarity is saturated at p50 1.000 and cannot discriminate at all), so no similarity enters this
#: number. See `within_node_propagation`'s header.
_PROPAGATION: dict[tuple[int, AnnotationKind], dict[str, PropagationRate]] = {}

def _annotated_counts(session: Session, locus_ids: set[int]) -> dict[int, dict[AnnotationKind, int]]:
    """Per locus, how many member genes carry each vocabulary — the denominator of `checkable`.

    One bounded statement over at most six loci (the focal node and its five neighbours) gives all
    four axes on the same definition. ⚠ `ANNOTATED_COUNT_PREDICATE` lives in
    `instruments/annotation_transfer.py` (moved 2026-10-04) because `within_node_propagation` needs
    the same mapping and an instrument cannot import a service; it is derived there from
    `GENE_VALUE_COLUMN`, so the column and the predicate are one fact in one place.
    """
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
    return {row.locus_id: {kind: getattr(row, kind.name.lower()) for kind in ANNOTATED_COUNT_PREDICATE} for row in rows}


def clear_calibration_cache() -> None:
    """Drop every cached ladder. Call between catalogues in tests."""
    _CALIBRATIONS.clear()


def clear_propagation_cache() -> None:
    """Drop every cached base rate. ⛔ Call between catalogues in tests, like its twin above."""
    _PROPAGATION.clear()


def propagation_for(
    session: Session,
    *,
    pangenome_id: int,
    annotation_kind: AnnotationKind,
    prevalence_band: PrevalenceBand,
) -> PropagationRate:
    """How often a call in this band covers the node's unannotated members — one band of one pass.

    ⛔ **A band absent from the tally is returned as zeros, never dropped.** A band with no annotated
    node at all and a band with 29 checkable ones must reach the page as the SAME shape carrying the
    same `rate: None`, or the page has to distinguish *missing key* from *not measurable* — and three
    states collapsing into two is the mistake `no_homology` was retired for.
    """
    key = (pangenome_id, annotation_kind)
    if key not in _PROPAGATION:
        _PROPAGATION[key] = compute_propagation(
            session, pangenome_id=pangenome_id, annotation_kind=annotation_kind
        )
    measured = _PROPAGATION[key].get(prevalence_band.name)
    if measured is not None:
        return measured
    return PropagationRate(
        annotation_kind=annotation_kind,
        prevalence_band=prevalence_band.name,
        checkable=0,
        unanimous=0,
        one_gene=0,
    )


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
        # ⛔ **ORDER BY, because `said[0]` is load-bearing.** `_fold` and `_walk` take the FIRST
        # entry of each (locus, kind) list, and a GO-bearing locus has THREE rank-0 rows — one per
        # namespace. With no ordering the term the card printed was whichever namespace Postgres
        # happened to return first, which is arbitrary and can change between runs. Namespace order
        # (molecular function, biological process, cellular component) then term, so it is at least
        # the same answer every time. ⚠ Stage 2 replaces the choice with all three, faceted; until
        # then determinism is the fix, not the final answer.
        .order_by(
            LocusAnnotationEntry.locus_id,
            LocusAnnotationEntry.annotation_kind,
            LocusAnnotationEntry.gene_ontology_namespace.nulls_first(),
            LocusAnnotationEntry.term_value,
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
            # ⚠ Read for the DONOR's propagation rate: a donor resting on one annotated gene is a
            # valid donor, and the reason is its own band's base rate, not the focal node's.
            Locus.prevalence_band,
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
                # ⭐ And the whole point of THIS one: whether or not the vote could be checked, the
                # call covers the node's unannotated members, at the measured rate for its band.
                # Reported on every node that carries a call, not only the one-gene ones — the claim
                # and the number are the same in both cases and only the first clause differs.
                "propagation": propagation_for(
                    session,
                    pangenome_id=pangenome_id,
                    annotation_kind=annotation_kind,
                    prevalence_band=locus.prevalence_band,
                ).as_json(),
            }
        vocabularies.append(
            {
                "annotation_kind": annotation_kind.value,
                "own": own,
                # ⚠ Only ever offered where the node has NO call of its own — a suggestion beside a
                # real annotation would compete with it.
                **(
                    _walk(annotation_kind, neighbours, entries, annotated_counts, session, pangenome_id, representation)
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
        ladder = calibration_for(
            session,
            pangenome_id=pangenome_id,
            annotation_kind=annotation_kind,
            representation=representation,
        )
        #: ⭐ **Two levels, and the difference between them is the whole floor.** `quoted` is how
        #: deep this tier and this donor COULD be read; `level` is how deep it may actually be
        #: SUGGESTED, which also requires the measured rate at that depth to clear `CALLING_FLOOR`.
        #: ⛔ Where they differ the card must still name the donor and quote the rate — otherwise a
        #: refusal at 24.3 % agreement reads as "too remote to call", which it is not: the neighbour
        #: is at 0.929 and perfectly close. It is the *measurement* that refuses, and the reader is
        #: owed the number that did it.
        measurable = folded and tier != NOT_CALLED
        quoted = quotable_level(annotation_kind, tier, folded) if measurable else None
        level = callable_level(annotation_kind, tier, folded, ladder) if measurable else None
        cell = ladder.cell(tier, level if level is not None else quoted) if quoted is not None else None
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
                # ⛔ The DONOR's band, not the recipient's — what is being claimed here is that the
                # donor's call covers the donor's node, which is a statement about that node alone.
                "propagation": propagation_for(
                    session,
                    pangenome_id=pangenome_id,
                    annotation_kind=annotation_kind,
                    prevalence_band=row.prevalence_band,
                ).as_json(),
            },
            # `null` where the tier is too remote to call, or the donor states nothing this deep.
            "level": level,
            "value": sorted(folded[level]) if level is not None else None,
            "calibration": cell.as_json() if cell is not None else None,
        }
    return {"walk": walk, "candidate": candidate}
