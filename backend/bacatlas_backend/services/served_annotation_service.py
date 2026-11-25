"""The Function tab's predictions — what annotation-transfer suggests for a locus, each with its own accuracy.

⭐ **This replaced the 5-nearest-loci ladder** (David, 2026-11-14): annotation-transfer (v0.1; v0.2 from
2026-11-23) copies the best donor's classes, single donor, and serves a prediction only at a calibrated
accuracy of 0.80 or more. The page shows each prediction's OWN number and no bands (David, 2026-10-28), and
the version's label (`annotation_model_label`) beside them.

⛔ **A donor in another catalogue is linked only where that catalogue is PUBLISHED on this site**
(`donor.catalogue_published` = the donor pangenome's `is_published`, which means *offer this one in the
picker*). The ingest guarantees the donor's catalogue is in this database, not that it is offered: a
catalogue loaded only to donate, or staged and not yet shown, would otherwise be one click away from a
reader. The page names such a donor — catalogue, locus and gene — without a link.

⛔ **Five answers, keyed by `state`, and none may read like another** (review of 9cfa481 found three of
them sharing one sentence):

* `loaded`, with calls — the predictions, one per claim;
* `loaded`, no calls — the model DECLINED this locus (no prediction reached the bar for it);
* `none_produced` — an operator stated the model produced none for this catalogue: a finding;
* `not_loaded` — nobody loaded this catalogue's predictions (or a reload removed them): a statement about
  this database, never about the model;
* `unknown_model` — loaded from a version this build does not register. ⛔ FAIL CLOSED: no calls, because a
  prediction must be shown with the bar it cleared and this build cannot state that bar.

⚠ **One statement**, keyed by the locus — the unique index leads with `locus_id` — plus the GO fold the
names come from, which is cached per process (`reference_name_service.gene_ontology_fold`). Nothing here
aggregates over the catalogue: the version is read off the pangenome row the route already holds.
"""

from __future__ import annotations

from hashlib import blake2b

from sqlalchemy import select
from sqlalchemy.orm import Session, aliased

from bacatlas_backend.models.enumerations import (
    ServedAnnotationKind,
    ServedAnnotationState,
    gene_ontology_namespace_name,
)
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.models.served_locus_annotation import (
    KNOWN_ANNOTATION_MODEL_VERSIONS,
    ServedLocusAnnotation,
    annotation_model_label,
)
from bacatlas_backend.services.reference_name_service import gene_ontology_fold


def load_predictions(session: Session, *, pangenome: Pangenome, locus: Locus) -> dict:
    """The `predictions` block of the Function tab for one locus."""
    state = pangenome.served_annotation_state
    version = pangenome.annotation_model_version
    known = KNOWN_ANNOTATION_MODEL_VERSIONS.get(version) if version is not None else None
    block = {
        "state": state.value,
        "annotation_model_version": version,
        "annotation_model_label": annotation_model_label(version),
        "serve_cut": known.serve_cut if known is not None else None,
        "calls": [],
    }
    if state is not ServedAnnotationState.LOADED:
        return block
    if known is None:
        # ⛔ Fail closed: rows this build cannot put a bar beside are not shown at all.
        return {**block, "state": "unknown_model"}

    donor = aliased(Locus)
    donor_pangenome = aliased(Pangenome)
    donor_species = aliased(PathogenSpecies)
    rows = session.execute(
        select(
            ServedLocusAnnotation.annotation_kind,
            ServedLocusAnnotation.gene_ontology_namespace,
            ServedLocusAnnotation.term_values,
            ServedLocusAnnotation.estimated_accuracy,
            ServedLocusAnnotation.esm_medoid_similarity,
            ServedLocusAnnotation.bacformer_medoid_similarity,
            ServedLocusAnnotation.annotation_model_version,
            donor.node_label,
            donor.display_name,
            donor.member_gene_count,
            donor.pangenome_id,
            donor_pangenome.catalogue_key,
            donor_pangenome.is_published,
            donor_species.scientific_name,
        )
        .join(donor, donor.locus_id == ServedLocusAnnotation.donor_locus_id)
        .join(donor_pangenome, donor_pangenome.pangenome_id == donor.pangenome_id)
        .join(donor_species, donor_species.pathogen_species_id == donor_pangenome.pathogen_species_id)
        .where(ServedLocusAnnotation.locus_id == locus.locus_id)
        # COG first, then GO in the namespace order the cards use (molecular function, process, component).
        .order_by(ServedLocusAnnotation.annotation_kind, ServedLocusAnnotation.gene_ontology_namespace)
    ).all()

    fold = gene_ontology_fold(session) if any(
        row.annotation_kind is ServedAnnotationKind.GENE_ONTOLOGY_SLIM for row in rows
    ) else None
    calls = []
    for row in rows:
        is_go = row.annotation_kind is ServedAnnotationKind.GENE_ONTOLOGY_SLIM
        calls.append(
            {
                "annotation_kind": row.annotation_kind.value,
                "gene_ontology_namespace": gene_ontology_namespace_name(row.gene_ontology_namespace),
                # ⭐ Most specific first, as stored. GO is named here from the reference; ⚠ a COG category
                # letter is named on the page (`functionVocabulary.COG_CATEGORIES`), as the COG card
                # already names `modal_cog_categories`, so `name` is null for COG by design.
                "terms": [
                    {"value": value, "name": fold.name_of.get(value) if is_go and fold is not None else None}
                    for value in row.term_values
                ],
                # ⛔ The row's OWN calibrated number — the page shows it, never a band.
                "estimated_accuracy": row.estimated_accuracy,
                "esm_medoid_similarity": row.esm_medoid_similarity,
                "bacformer_medoid_similarity": row.bacformer_medoid_similarity,
                "annotation_model_version": row.annotation_model_version,
                "donor": {
                    "node_label": row.node_label,
                    "display_name": row.display_name,
                    "gene_count": row.member_gene_count,
                    # ⭐ The donor may live in ANOTHER catalogue (E. coli → Kp); the page then links across
                    # catalogues by this key, because locus labels are catalogue-private.
                    "catalogue_key": row.catalogue_key,
                    "species_name": row.scientific_name,
                    "same_catalogue": row.pangenome_id == pangenome.pangenome_id,
                    # ⛔ Whether a LINK to the donor may be offered: its catalogue is published on this site.
                    # Read from the donor's OWN pangenome row, so it is part of this body that neither
                    # the recipient's id nor its load fingerprint tracks — `predictions_data_version`.
                    "catalogue_published": bool(row.is_published),
                },
            }
        )
    return {**block, "calls": calls}


def predictions_data_version(pangenome: Pangenome, predictions: dict) -> str:
    """What the Function tab's validator keys on for its `predictions` block.

    ⭐ The catalogue's load fingerprint (state, version, the file's sha256), and the published state of every
    donor's catalogue. ⛔ Both change the body without moving the recipient's pangenome id: a reload of a
    corrected file (review of 9cfa481), and publishing or hiding a donor's catalogue, which turns a named donor
    into a link or back. A tag that tracks only the first would replay a stale link state
    (`a-validator-must-track-shape-not-just-data`).

    ⚠ Every donor, the same-catalogue ones included: the page walks to those whatever their flag says, but
    the flag is still in the body, and a validator tracks the body, not what one client makes of it. A locus
    with no call carries no flag, so its version is the fingerprint alone, exactly as before.
    """
    fingerprint = pangenome.served_annotation_fingerprint()
    donors = sorted(
        {(call["donor"]["catalogue_key"], call["donor"]["catalogue_published"]) for call in predictions["calls"]}
    )
    if not donors:
        return fingerprint
    published = ",".join(f"{key}={int(flag)}" for key, flag in donors)
    return f"{fingerprint}-d{blake2b(published.encode(), digest_size=4).hexdigest()}"
