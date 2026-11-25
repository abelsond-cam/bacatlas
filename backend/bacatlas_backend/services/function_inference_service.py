"""The Function tab's own-call block — what this node's function rests on, and how sure that is.

* **its own genes.** The node's call is already a modal vote over its member genes, and `checkable`
  says whether that vote had more than one voter. ⛔ A node with exactly ONE annotated gene is
  unanimous *by construction*; 518 of 4,993 *E. coli* COG nodes are in that state, `gumC` among them
  (COG3206 on 1 of 100 genes), and a card that showed it like any other would be claiming a check
  that never happened. ⭐ **But that is a missing check, NOT a missing claim** — the call still
  covers the node's other 99 genes, and `propagation` is the measured rate at which it is right to.
  (David, 2026-10-04: *"the whole point of our method is effectively the carrying by one is
  enough."*)

⛔ **The second mechanism that lived here — walk the five stored ESM neighbours until one carries the
vocabulary, and quote it as deep as a measured agreement ladder allowed — is RETIRED** (David,
2026-11-14). annotation-transfer v0.1 replaces it: `served_annotation_service` serves the model's
predictions, each with its own calibrated accuracy. The ladder's MEASUREMENT survives in
`instruments/annotation_transfer.py`, which `within_node_propagation`, `transfer/claims` and the
propagation ingest read; only the per-request walk, its process-cached calibration and the route's
`calibration` block went.
"""

from __future__ import annotations

from collections import defaultdict

from sqlalchemy import select, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import LOCUS_ANNOTATED_COUNT_COLUMN, SUPPORTED_KINDS
from bacatlas_backend.instruments.within_node_propagation import PropagationRate
from bacatlas_backend.models.enumerations import AnnotationKind, PrevalenceBand
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry
from bacatlas_backend.models.pangenome_propagation_rate import PangenomePropagationRate
from bacatlas_backend.services.reference_name_service import reference_names

#: ⚠ A process-lifetime cache, legitimate only because this service is READ-ONLY and a catalogue
#: cannot change without a re-ingest and a restart — the ground the `ETag: "{pangenome_id}"` on every
#: response already stands on. It holds at most five stored rows per (catalogue, vocabulary)
#: (`_stored_propagation`). ⛔ Keyed WITHOUT a representation, which is a finding and not an omission:
#: within-node unanimity is 99.8-100 % at every band of BOTH representations (ESM's within-node
#: similarity is saturated at p50 1.000 and cannot discriminate at all), so no similarity enters this
#: number. See `within_node_propagation`'s header. ⛔ Tests call `clear_propagation_cache()` between
#: catalogues, or one fixture's rates answer another's.
_PROPAGATION: dict[tuple[int, AnnotationKind], dict[str, PropagationRate]] = {}


def _annotated_counts(session: Session, locus_ids: set[int]) -> dict[int, dict[AnnotationKind, int]]:
    """Per locus, how many member genes carry each vocabulary — the denominator of `checkable`.

    ⭐ **Four stored columns, read from the locus row** (`LOCUS_ANNOTATED_COUNT_COLUMN`, which
    records the gene-grain recount they were verified against: 0 mismatching loci across all four
    catalogues). GO's `_any` column is what made the fourth axis expressible.
    """
    if not locus_ids:
        return {}
    columns = ", ".join(
        f"{column} as {kind.name.lower()}" for kind, column in LOCUS_ANNOTATED_COUNT_COLUMN.items()
    )
    rows = session.execute(
        text(f"select locus_id, {columns} from locus where locus_id = any(:locus_ids)"),
        {"locus_ids": list(locus_ids)},
    ).all()
    # ⚠ `0` where the column is NULL — a catalogue ingested before `go_annotated_member_count_any`
    # existed holds NULL, and `checkable` is a comparison against it. ⛔ That is a MEASURED-ZERO
    # substitution and it is the one place in this file where NULL and 0 are deliberately conflated:
    # the alternative is a `None` reaching `annotated >= 2` on the page. The gate that makes it safe
    # is `test_locus_annotation_rollups`, which refuses a catalogue holding any NULL here.
    return {
        row.locus_id: {
            kind: getattr(row, kind.name.lower()) or 0 for kind in LOCUS_ANNOTATED_COUNT_COLUMN
        }
        for row in rows
    }


def clear_propagation_cache() -> None:
    """Drop every cached base rate. ⛔ Call between catalogues in tests."""
    _PROPAGATION.clear()


def _stored_propagation(
    session: Session, *, pangenome_id: int, annotation_kind: AnnotationKind
) -> dict[str, PropagationRate]:
    """One catalogue's band tally for one vocabulary, read from `pangenome_propagation_rate`.

    ⭐⭐ **This is what the stage bought.** It used to be `within_node_propagation`'s
    `compute_propagation`, a scan over **every annotated gene row of the catalogue** — measured at
    0.22-1.94 s per (catalogue, vocabulary), so a cold *E. coli* Function tab paid ~5.8 s across the
    four, and at the 80,000-genome design target the same scan is over ~412 M rows. It is now a
    primary-key range read of at most five rows.

    ⛔ **It still returns `PropagationRate` objects and not the stored numbers.** The table has no
    `rate` column, by design: `PropagationRate.rate` returns `None` below `MIN_NODES = 30`, and
    anything computing `unanimous / checkable` off these rows re-introduces the 0.0 that property
    refuses — on RARE, which has zero checkable nodes in every catalogue measured. The floor, the
    Wilson interval and the band's NAME→VALUE conversion all live on that dataclass, so the read
    goes through it.

    ⚠ The cache above keeps the same shape and the same justification it always had, but now holds
    twenty rows rather than a measurement: the catalogue cannot change without a re-ingest and a
    restart, which is the ground the `ETag: "{pangenome_id}"` already stands on.
    """
    rows = session.execute(
        select(
            PangenomePropagationRate.prevalence_band,
            PangenomePropagationRate.checkable,
            PangenomePropagationRate.unanimous,
            PangenomePropagationRate.one_gene,
        ).where(
            PangenomePropagationRate.pangenome_id == pangenome_id,
            PangenomePropagationRate.annotation_kind == annotation_kind,
        )
    ).all()
    return {
        band.name: PropagationRate(
            annotation_kind=annotation_kind,
            # ⚠ The enum's NAME, because that is the key `propagation_for` looks up and what
            # `compute_propagation` used to return. `as_json` emits the VALUE; see there.
            prevalence_band=band.name,
            checkable=checkable,
            unanimous=unanimous,
            one_gene=one_gene,
        )
        for band, checkable, unanimous, one_gene in rows
    }


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
        _PROPAGATION[key] = _stored_propagation(session, pangenome_id=pangenome_id, annotation_kind=annotation_kind)
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


def _top_entries(session: Session, locus_ids: set[int]) -> dict[tuple[int, AnnotationKind], list]:
    """Every supported vocabulary's rank-0 entries for every locus, in one statement.

    ⛔ A LIST per (locus, kind), not a single row: `top_go_slim` ranks within (locus, namespace), so a
    GO-bearing node states up to three rank-0 terms and reading one would discard the rest.
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
        )
        .where(
            LocusAnnotationEntry.locus_id.in_(locus_ids),
            LocusAnnotationEntry.annotation_kind.in_(SUPPORTED_KINDS),
            LocusAnnotationEntry.rank_within_locus == 0,
        )
        # ⛔ **ORDER BY, because `said[0]` is load-bearing.** `load_inference` takes the FIRST
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
    # ⭐ The EC and KEGG names are filled HERE rather than in each consumer. `term_name` still wins
    # where it exists; see `serialise_annotation_entry` for why that ordering is the guarantee.
    names = reference_names(session)
    grouped: dict[tuple[int, AnnotationKind], list] = defaultdict(list)
    for locus_id, kind, term, name, support in rows:
        grouped[(locus_id, kind)].append((term, name or names.name_for(kind, term), support))
    return grouped


def load_inference(session: Session, *, pangenome_id: int, locus: Locus) -> dict:
    """What this node's own genes say, per vocabulary, and the measured rate at which it propagates."""
    entries = _top_entries(session, {locus.locus_id})
    annotated_counts = _annotated_counts(session, {locus.locus_id})
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
        vocabularies.append({"annotation_kind": annotation_kind.value, "own": own})
    return {"vocabularies": vocabularies}
