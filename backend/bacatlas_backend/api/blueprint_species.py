"""`/api/v1/species…` — the catalogue, one locus, its arrangements, its function tab, and search.

⭐ **Every response is immutable for a build**, because writes are ours alone and offline. So each
carries `ETag: "{pangenome_id}"` and a long `Cache-Control`, and the cache is a pure function of its
key with nothing that can invalidate it — the single largest payoff of the read-only scope.

⛔ **`locus_label` is the only locus identifier that appears in a URL.** It is the existing page's
hash, so every link a reader has saved keeps working, and it is durable across a re-export while the
surrogate id and the catalogue ordinal are not. The integer id appears only inside responses.

⚠ **A failed request and an empty result must never render the same.** The rule generalises from
`app.js:4604` — *"a sequence panel that fails silently is one a reader will read as 'this genome has
nothing here', which is a different claim and a false one"* — so a 404 says which of the species, the
pangenome or the locus was missing, and never returns an empty body with a 200.
"""

from __future__ import annotations

from flask import Blueprint, current_app, jsonify, request
from sqlalchemy import select

from bacatlas_backend.instruments.annotation_transfer import SUPPORTED_KINDS
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.serialisers.locus_serialiser import (
    serialise_annotation_entry,
    serialise_arrangement,
    serialise_gene_sequence,
    serialise_locus_detail,
)
from bacatlas_backend.serialisers.projected_genome_serialiser import serialise_projected_genome
from bacatlas_backend.services.audit_residual_service import load_audit_residuals
from bacatlas_backend.services.function_inference_service import calibration_for, load_inference
from bacatlas_backend.services.gene_sequence_service import (
    SequenceUnavailable,
    load_gene_sequences,
)
from bacatlas_backend.services.genome_listing_service import clamp_genome_limit, list_genomes
from bacatlas_backend.services.locus_detail_service import (
    LocusNotFound,
    load_arrangement_page,
    load_function_block,
    load_locus_detail,
)
from bacatlas_backend.services.locus_search_service import DEFAULT_RESULT_LIMIT, search_loci
from bacatlas_backend.services.projected_genome_service import (
    list_projected_genomes,
    resolve_projected_genome,
)
from bacatlas_backend.services.species_catalogue_service import (
    SpeciesNotPublished,
    list_catalogues,
    list_published_species,
    load_catalogue,
    load_species_catalogue,
    resolve_catalogue,
    resolve_published_pangenome,
)

species_blueprint = Blueprint("species", __name__)

#: One day. Long, because a response cannot go stale without a deploy that changes the ETag.
CACHE_CONTROL = "public, max-age=86400"

#: Arrangements per page past the display cut.
ARRANGEMENT_PAGE_LIMIT = 50


def _immutable(payload, pangenome_id: int | None):
    """Attach the caching headers every response of a built catalogue is entitled to."""
    response = jsonify(payload)
    response.headers["Cache-Control"] = CACHE_CONTROL
    if pangenome_id is not None:
        response.headers["ETag"] = f'"{pangenome_id}"'
    return response


def _session():
    return current_app.extensions["bacatlas_database"].session()


def _not_found(message: str):
    """⚠ A named 404, never an empty 200 — see the module docstring."""
    return jsonify({"error": "not_found", "detail": message}), 404


@species_blueprint.get("/species")
def get_species():
    """Every species and the pangenome it serves. Replaces `published.tsv`."""
    with _session() as session:
        rows = list_published_species(session)
        return _immutable(
            {
                "species": [
                    {
                        "key": species.species_key,
                        "scientific_name": species.scientific_name,
                        "ncbi_taxonomy_id": species.ncbi_taxonomy_id,
                        # ⚠ A species with nothing published is LISTED, with `null`. Filtering it
                        # out would read as though the species did not exist.
                        "published": pangenome is not None,
                        "genome_count": pangenome.genome_count if pangenome else None,
                        "gene_count": pangenome.gene_count if pangenome else None,
                        "locus_count": pangenome.locus_count if pangenome else None,
                        # ⭐ The species' DEFAULT catalogue — what a bare species request resolves
                        # to, and the key a page should address once it has followed one. The full
                        # menu is `GET /catalogues`.
                        # ⚠ This replaced a field called `model_label` that was a hardcoded `None`.
                        # A name that describes something other than its content is worse than an
                        # absent field, and nothing read it.
                        "catalogue_key": pangenome.catalogue_key if pangenome else None,
                    }
                    for species, pangenome in rows
                ]
            },
            None,
        )


@species_blueprint.get("/catalogues")
def get_catalogues():
    """Every catalogue on offer — what the model picker is built from.

    ⭐ A species appears once per catalogue it holds, so E. coli on nuna4 and nuna5 is two rows.
    `is_default` marks the one a bare species request resolves to; `key` is what a URL should carry
    to pin this clustering.

    ⚠ Listed from `pangenome.is_published`, which means *offer this one*, NOT
    `default_pangenome_id`. A catalogue that is loaded but not offered is reachable by key and absent
    from the menu — which is how a model is staged before anyone is shown it.
    """
    with _session() as session:
        rows = list_catalogues(session)
        return _immutable(
            {
                "catalogues": [
                    {
                        "key": pangenome.catalogue_key,
                        "species": {
                            "key": species.species_key,
                            "scientific_name": species.scientific_name,
                        },
                        "model": None
                        if model is None
                        else {
                            "key": model.model_key,
                            "label": model.label,
                            "step_count": model.step_count,
                            "exclusivity_form": model.exclusivity_form.value,
                        },
                        "is_default": species.default_pangenome_id == pangenome.pangenome_id,
                        "genome_count": pangenome.genome_count,
                        "gene_count": pangenome.gene_count,
                        "locus_count": pangenome.locus_count,
                        "run_id": pangenome.run_id,
                    }
                    for pangenome, species, model in rows
                ]
            },
            None,
        )


@species_blueprint.get("/catalogues/<species_key>")
@species_blueprint.get("/species/<species_key>")
def get_species_catalogue(species_key: str):
    """The census, the model provenance and the landing locus — everything before a locus."""
    with _session() as session:
        try:
            # A hyphen means a catalogue key and pins the model; a bare species follows its default.
            # Same rule as `_resolve_pangenome`, and the two must not drift apart.
            catalogue = (
                load_catalogue(session, species_key)
                if "-" in species_key
                else load_species_catalogue(session, species_key)
            )
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        pangenome = catalogue.pangenome
        return _immutable(
            {
                "species": {
                    "key": catalogue.species.species_key,
                    "scientific_name": catalogue.species.scientific_name,
                },
                "pangenome": {
                    # ⭐ The key this catalogue is addressed by. Without it a page that followed a
                    # bare species cannot name what it was given, so it cannot pin the address to
                    # the clustering the reader is actually looking at — and it cannot be
                    # reconstructed from (species, model) either, because that pair is not unique
                    # and `catalogue_key_for` admits deliberate keys like `ecoli-sensitive`.
                    "catalogue_key": pangenome.catalogue_key,
                    "run_id": pangenome.run_id,
                    "genome_count": pangenome.genome_count,
                    "gene_count": pangenome.gene_count,
                    "locus_count": pangenome.locus_count,
                    # ⛔ `built_at` is NULL on all four loaded catalogues — nothing in `ingest`
                    # writes it — so the footer printed "built —" and the page could not answer
                    # "is this current?", which is the whole complaint it exists to settle.
                    # `ingested_at` IS populated (2026-09-25 for the nuna4 pair) and is the honest
                    # answer to "how old is the catalogue THIS stack is serving": when it was loaded.
                    # Both ship; the footer prefers whichever it has.
                    "built_at": pangenome.built_at,
                    "ingested_at": (
                        pangenome.ingested_at.isoformat() if pangenome.ingested_at else None
                    ),
                    "git_sha": pangenome.git_sha,
                    "exclusivity_form": pangenome.exclusivity_form.value,
                    # ⚠ An absent section is NAMED. The page prints these rather than showing an
                    # empty panel, so an omission is never mistaken for a measured zero.
                    "omitted_sections": pangenome.omitted_sections or {},
                },
                "model": (
                    {
                        "key": catalogue.model.model_key,
                        "label": catalogue.model.label,
                        "exclusivity_form": catalogue.model.exclusivity_form.value,
                        "knn_k": catalogue.model.knn_k,
                        "step_count": catalogue.model.step_count,
                    }
                    if catalogue.model
                    else None
                ),
                # ⭐ One row per pipeline stage, as the footer prints them — read from the run's own
                # manifests, never decoded from the run_id.
                "steps": [
                    {
                        "ordinal": step.step_ordinal,
                        "name": step.step_name,
                        "representation": step.representation.value if step.representation else None,
                        "gamma": step.gamma,
                        "rho_rule": step.rho_rule.value if step.rho_rule else None,
                        "rho_ceiling": step.rho_ceiling,
                        "uses_exclusivity": step.uses_exclusivity,
                        "stage_name": step.stage_name,
                        "stage_detail": step.stage_detail,
                    }
                    for step in catalogue.steps
                ],
                "provenance_rows": pangenome.provenance_rows or [],
                "prevalence_census": catalogue.prevalence_census,
                # ⭐ The same bands counted by GENE. Divided by `genome_count` this is the page's
                # "per genome" line, and its three parts sum to the whole by construction.
                "prevalence_gene_census": catalogue.prevalence_gene_census,
                "audit_headline": catalogue.audit_headline,
                # ⛔ `map_projections` was REPLACED, not renamed. It carried a UMAP over every
                # locus's MEDOID, that fit's extent and scale, a whole-catalogue sprite descriptor,
                # and a null baseline sampled over random pairs of MEDOIDS. All of it described a
                # construction that reduced a locus to one member.
                "similarity_baselines": [
                    {
                        "representation": baseline.representation.value,
                        "form": baseline.similarity_form,
                        # ⛔ The random GENE-pair floor, and NOT the retired medoid null: on the
                        # published ecoli catalogue those sit at 0.0587 and 0.0651, close enough to
                        # look interchangeable and not be. Without it a cosine has no scale — ESM's
                        # floor is ~0.742 and Bacformer's ~0.059, so the same 0.41 reads oppositely.
                        "floor_median": baseline.floor_median,
                        # ⭐ Its spread, so the card can draw a box rather than a bare tick: a median
                        # alone cannot say whether 0.41 sits far outside random or inside its
                        # shoulder. `null` where the run's audit JSON was not beside its CSV.
                        "floor_p25": baseline.floor_p25,
                        "floor_p75": baseline.floor_p75,
                        "floor_p99": baseline.floor_p99,
                        # ⭐ The other half of "p12 of 12,104 loci".
                        "measurable_locus_count": baseline.measurable_locus_count,
                        "neighbour_knn_k": baseline.neighbour_knn_k,
                    }
                    for baseline in catalogue.similarity_baselines
                ],
                "landing_locus": catalogue.landing_locus_label,
                "example_loci": catalogue.example_locus_labels,
                "example_locus_rows": catalogue.example_locus_rows,
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/audit/residual-loci")
@species_blueprint.get("/species/<species_key>/audit/residual-loci")
def get_audit_residual_loci(species_key: str):
    """Every locus grouped on context alone, and every Pfam conflict — two lists, each clickable.

    Fetched when the footer's list is opened. ⚠ The two lists are NOT disjoint: a locus can be both,
    and it then appears in both, exactly as the published footer listed it.
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        residuals = load_audit_residuals(session, pangenome_id=pangenome.pangenome_id)

        def rows(entries):
            return [
                {
                    "label": entry.label,
                    "display_name": entry.display_name,
                    "prevalence_band": entry.prevalence_band,
                    "gene_count": entry.gene_count,
                    "uniref50_family_count": entry.uniref50_family_count,
                    "pfam_architecture_count": entry.pfam_architecture_count,
                    "syntenic_a5": entry.syntenic_a5,
                    # ⚠ DISTANCES, as stored. The footer has always shown similarities; the client
                    # converts, exactly as it does for the locus card, so there is one conversion.
                    # ⛔ SIMILARITIES now, not distances: the client used to convert with
                    # `1 − d` and there is nothing left to convert. The numbers behind them are a
                    # median over every within-locus gene pair and the highest such median against
                    # another locus, not a member's distance to one gene.
                    "esm_within_similarity": entry.esm_within_similarity,
                    "esm_nearest_similarity": entry.esm_nearest_similarity,
                }
                for entry in entries
            ]

        return _immutable(
            {
                "grouped_on_context_alone": rows(residuals.grouped_on_context_alone),
                "pfam_conflicts": rows(residuals.pfam_conflicts),
            },
            pangenome.pangenome_id,
        )


def _resolve_catalogue_genome(session, pangenome, sample_id: str) -> int | None:
    """A genome of THIS pangenome's collection, by sample id — or `None`. One statement.

    ⛔⛔ **The collection, not the species.** The database holds 122 ecoli and 158 kp genomes against
    100-genome pangenomes, so a species-scoped lookup still resolved 22 + 58 genomes this catalogue
    never modelled. Anchored, one read as *"anchored — this genome has no gene at this locus"*; on the
    Sequence tab as `genes: []`, *"has no gene here"* — a claim about a genome this pangenome says
    nothing about at all. (First unscoped, then species-scoped in `acb8ea4`, which still left these.)
    """
    from bacatlas_backend.models.genome_collection import GenomeCollectionMembership

    return session.execute(
        select(Genome.genome_id)
        .join(GenomeCollectionMembership, GenomeCollectionMembership.genome_id == Genome.genome_id)
        .where(
            Genome.sample_id == sample_id,
            GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id,
        )
    ).scalar_one_or_none()


def _resolve_pangenome(session, species_key: str):
    """The catalogue this request addresses — a catalogue key if it is one, else a species' default.

    ⭐ **Both forms share one path segment, and they cannot collide.** `ecoli-nuna5` names a
    catalogue exactly; `ecoli` names a species and gets whatever that species currently serves.
    Catalogue keys always contain a hyphen and species keys never do — which is the same property
    that makes the keys prefix-free, doing a second job here.

    ⚠ The difference is not cosmetic: a request naming a catalogue is PINNED and keeps meaning the
    same clustering after the default moves, which is the whole point of putting the model in the
    URL. A request naming a species deliberately follows the default.

    ⛔ One statement either way. See `resolve_published_pangenome` for what this used to cost every
    route.
    """
    if "-" in species_key:
        return resolve_catalogue(session, species_key)
    return resolve_published_pangenome(session, species_key)


# ⛔ `GET /species/<key>/map/<rep>/scatter.png` was DELETED on 2026-09-24 with the whole-catalogue
# sprite it served. It was the only endpoint here that did not return JSON.
#
# ⚠ Worth keeping if a catalogue-wide picture ever returns: its ETag was the CONTENT digest and not the
# pangenome id, because an image is cached hard and for a long time by browsers and proxies we do not
# control, so a re-render under the same id had to be able to invalidate it.


@species_blueprint.get("/catalogues/<species_key>/loci/<path:locus_label>")
@species_blueprint.get("/species/<species_key>/loci/<path:locus_label>")
def get_locus(species_key: str, locus_label: str):
    """⭐ The hot path. One round trip, and the popover is then offline.

    `anchor` names a genome by `sample_id`; given one, the arrangement that genome carries is
    included even if it sits past the display cap.
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))

        anchor_genome_id = None
        anchor = request.args.get("anchor")
        projected = request.args.get("projected")
        # ⛔ One genome is anchored at a time, and the two are DIFFERENT relations to this
        # pangenome — a member and a genome placed beside it. Serving both would have to pick one
        # silently, so it is a named 400 instead.
        if anchor and projected:
            return (
                jsonify(
                    error="anchor and projected are mutually exclusive: a genome is either one of "
                    "this pangenome's modelled genomes or one placed onto it afterwards, and the "
                    "page anchors to one genome at a time"
                ),
                400,
            )
        if anchor:
            anchor_genome_id = _resolve_catalogue_genome(session, pangenome, anchor)
            if anchor_genome_id is None:
                return _not_found(f"no genome {anchor!r} in the {species_key!r} catalogue")

        projected_genome_id = None
        if projected:
            projected_genome_id = resolve_projected_genome(session, pangenome.pangenome_id, projected)
            # ⚠ A genome that EXISTS but was never projected is a 404 that says which, rather than
            # an empty 200 — the database holds 22 ecoli and 58 kp genomes in neither the collection
            # nor any projection, and a blank answer for one reads as "it has nothing here".
            if projected_genome_id is None:
                return _not_found(
                    f"no genome {projected!r} placed on the {species_key!r} catalogue — it is not one "
                    "of the genomes projected onto this pangenome"
                )

        try:
            detail = load_locus_detail(
                session,
                pangenome_id=pangenome.pangenome_id,
                node_label=locus_label,
                anchor_genome_id=anchor_genome_id,
                projected_genome_id=projected_genome_id,
            )
        except LocusNotFound as error:
            return _not_found(str(error))

        payload = serialise_locus_detail(detail)
        payload["resolved_neighbour_count"] = detail.resolved_neighbour_count
        return _immutable(payload, pangenome.pangenome_id)


@species_blueprint.get("/catalogues/<species_key>/loci/<path:locus_label>/arrangements")
@species_blueprint.get("/species/<species_key>/loci/<path:locus_label>/arrangements")
def get_locus_arrangements(species_key: str, locus_label: str):
    """Arrangements past the display cut — the full scroller, paged."""
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        try:
            detail = load_locus_detail(
                session,
                pangenome_id=pangenome.pangenome_id,
                node_label=locus_label,
                arrangement_limit=0,
            )
        except LocusNotFound as error:
            return _not_found(str(error))
        offset = max(0, request.args.get("offset", default=0, type=int))
        page = load_arrangement_page(
            session, locus_id=detail.locus.locus_id, offset=offset, limit=ARRANGEMENT_PAGE_LIMIT
        )
        return _immutable(
            {
                "arrangements": [
                    serialise_arrangement(arrangement, detail.neighbour_display_rows)
                    for arrangement in page
                ],
                "offset": offset,
                "total": detail.locus.total_arrangement_count,
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/loci/<path:locus_label>/function")
@species_blueprint.get("/species/<species_key>/loci/<path:locus_label>/function")
def get_locus_function(species_key: str, locus_label: str):
    """The Function tab — fetched on tab open, not on every walk.

    ⭐ Carries its own **calibration**: the measured agreement rate for every (similarity tier,
    depth) the inference ladder can quote, recomputed from the loaded catalogue rather than copied
    out of the write-up. ~3.6 kB of cells on a cold-path response, so the reader can see the whole
    ladder instead of one rate presented without its neighbours.

    ⛔ **It no longer carries `go_verdicts`, and that is a retirement rather than an omission**
    (David, 2026-10-05). The verdict was `worst_relation` — Pfam's domain-architecture comparator —
    applied to GO slim classes with no ontology: `goslim_metagenomics` keeps `membrane`, `plasma
    membrane` and `outer membrane` as siblings, and the three namespace ROOTS besides, so two genes
    making compatible statements at different depths scored `disjoint`. 11 of the 43 `disjoint`
    verdicts in the two catalogues involve a root — a term that is an ancestor of everything.
    ⚠ The columns, the enum and the ingest stay: `instruments/payload_serialiser.py` writes the
    verdict into the STATIC payload the published pages read, and `test_catalogue_parity` pins the
    two to each other. This route simply stops serving it.

    ⭐ And its **propagation** base rate, which rides INSIDE `inference` rather than beside the
    calibration — nested on `own` and on `candidate.donor`, because it qualifies one node's call
    (how often that call covers the members not carrying it) rather than the ladder as a whole. It
    needs no top-level block for the same reason: a catalogue-wide table of band rates would be a
    second copy of the one cell the page actually quotes.
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        try:
            detail = load_locus_detail(
                session,
                pangenome_id=pangenome.pangenome_id,
                node_label=locus_label,
                arrangement_limit=0,
            )
        except LocusNotFound as error:
            return _not_found(str(error))
        locus = detail.locus
        grouped = load_function_block(session, locus_id=locus.locus_id)
        return _immutable(
            {
                # ⛔ Separate from `annotations`, never merged into it: one is what this node's genes
                # SAY, the other is what a neighbour SUGGESTS, and a reader cannot be left to tell
                # them apart by context.
                "inference": load_inference(
                    session, pangenome_id=pangenome.pangenome_id, locus=locus
                ),
                "calibration": {
                    kind.value: calibration_for(
                        session, pangenome_id=pangenome.pangenome_id, annotation_kind=kind
                    ).as_json()
                    for kind in SUPPORTED_KINDS
                },
                "annotations": {
                    kind: [serialise_annotation_entry(entry) for entry in entries]
                    for kind, entries in grouped.items()
                },
                # ⛔ Coverage BEFORE every verdict, and against the LOCUS size — a share against the
                # annotated subset would read 100 % where one gene in forty carries a label.
                "coverage": {
                    "gene_count": locus.member_gene_count,
                    "cog_annotated_gene_count": locus.cog_annotated_member_count,
                    # ⛔ A COUNT, never a relation: COG ids are single-valued per gene, so two in one
                    # locus is an ordinary consequence of grouping above family level.
                    "cog_distinct_id_count": locus.cog_distinct_id_count,
                    "modal_cog_categories": locus.modal_cog_categories,
                    "ec_annotated_gene_count": locus.ec_annotated_member_count,
                    "kegg_annotated_gene_count": locus.kegg_annotated_member_count,
                    "go_annotated_gene_count": {
                        "molecular_function": locus.go_annotated_member_count_molecular_function,
                        "biological_process": locus.go_annotated_member_count_biological_process,
                        "cellular_component": locus.go_annotated_member_count_cellular_component,
                    },
                },
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/genomes/<sample_id>/loci/<path:locus_label>/sequence")
@species_blueprint.get("/species/<species_key>/genomes/<sample_id>/loci/<path:locus_label>/sequence")
def get_gene_sequence(species_key: str, sample_id: str, locus_label: str):
    """⭐ The Sequence tab — a gene's DNA, its flanks and its protein, sliced from the original GFF.

    ⚠ **~5 kB out, where the published page pulled 1.3 MB in.** It had to: GitHub Pages applies
    `Range` to the compressed stream, so the browser fetched whole `.nseq` files and decoded DNA
    itself. A server slices the file Bakta wrote, and the custom format is retired rather than ported.

    ⛔ **Three different answers, and they must never look alike.** A genome with no gene at this
    locus is an ANSWER (`genes: []`); a gene whose bases could not be read is a FAILURE with a named
    reason; a locus or genome that does not exist is a 404. The published page's own rule, at
    `app.js:4604`: *"a sequence panel that fails silently is one a reader will read as 'this genome
    has nothing here', which is a different claim and a false one."*
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))

        genome_id = _resolve_catalogue_genome(session, pangenome, sample_id)
        if genome_id is None:
            return _not_found(f"no genome {sample_id!r} in the {species_key!r} catalogue")
        locus_id = session.execute(
            select(Locus.locus_id).where(
                Locus.pangenome_id == pangenome.pangenome_id, Locus.node_label == locus_label
            )
        ).scalar_one_or_none()
        if locus_id is None:
            return _not_found(f"no locus {locus_label!r} in {species_key}")

        configuration = current_app.config["BACATLAS"]
        try:
            gff_root = configuration.artifact_roots["gff"]
        except KeyError:
            # ⚠ A configuration failure, named as one. "Not configured" and "the file is missing"
            # are indistinguishable in a 404, and only one of them is ours to fix.
            return (
                jsonify(
                    {
                        "error": "sequence_unavailable",
                        "detail": "this server has no annotation store configured, so no sequence "
                        "can be read (set BACATLAS_ROOT_GFF)",
                    }
                ),
                503,
            )

        try:
            rows = load_gene_sequences(
                session,
                pangenome_id=pangenome.pangenome_id,
                node_label=locus_label,
                sample_id=sample_id,
                gff_root=gff_root,
            )
        except SequenceUnavailable as error:
            # ⛔ 503, not 404 and not an empty list: the gene is there and we could not read it.
            return jsonify({"error": "sequence_unavailable", "detail": str(error)}), 503

        return _immutable(
            {
                "genome": {"sample_id": sample_id},
                "locus": {"label": locus_label},
                # ⭐ A LIST, because rho > 1 puts one genome in a locus twice — and an empty one is
                # the answer "this genome has no gene here", not a failure.
                "genes": [serialise_gene_sequence(row) for row in rows],
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/search")
@species_blueprint.get("/species/<species_key>/search")
def get_search(species_key: str):
    """Substring search, with the same semantics the page has today. Replaces the 3.0 MB `HAY`."""
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        limit = min(request.args.get("limit", default=DEFAULT_RESULT_LIMIT, type=int), 100)
        result = search_loci(
            session,
            pangenome_id=pangenome.pangenome_id,
            query=request.args.get("q", default="", type=str),
            limit=limit,
        )
        return _immutable(
            {
                "query": result.query,
                # ⚠ Which mode answered. A 1–2 character query searches less of the haystack, and a
                # reader is entitled to know that rather than to conclude the catalogue is empty.
                "mode": result.mode,
                "truncated": result.truncated,
                "hits": [
                    {
                        "label": hit.node_label,
                        "display_name": hit.display_name,
                        "best_product": hit.best_product,
                        "gene_count": hit.member_gene_count,
                        "genome_count": hit.member_genome_count,
                        "prevalence_band": hit.prevalence_band,
                        "rank_band": hit.rank_band,
                    }
                    for hit in result.hits
                ],
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/genomes")
@species_blueprint.get("/species/<species_key>/genomes")
def get_genomes(species_key: str):
    """The anchor control's list — `meta.genomes` + `anchorSearch` + `GENOME_N`, as one query.

    `q` is a case-insensitive literal substring of the accession (empty matches every genome) and
    `limit` defaults to 100, held to [1, 1000]. ⭐ Two statements whatever either asks for: one to
    resolve the species, one for the genomes — their counts are READ from what ingest stored, never
    aggregated here.

    ⚠ **The ETag stays `"{pangenome_id}"` although `q` and `limit` change the body.** An entity tag
    distinguishes representations of ONE resource (RFC 9110 §8.8.3), and the query string is part of
    the target URI — so `?q=SAMEA2` and `?q=SAMEA3` are two resources, each cached and revalidated
    under its own URI, and for any one URI the body changes only when the pangenome does. Folding `q`
    into the tag would add nothing a cache uses; a cache keyed on the path alone would be broken by
    `max-age` long before any ETag was consulted. `/search` relies on the same reading.
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))
        listing = list_genomes(
            session,
            pangenome,
            query=request.args.get("q", default="", type=str),
            limit=clamp_genome_limit(request.args.get("limit", type=int)),
        )
        return _immutable(
            {
                "species_key": species_key,
                "query": listing.query,
                "genome_count": listing.genome_count,
                # ⛔ Before the limit. A list cut at the limit must not read as the whole answer.
                "matched_genome_count": listing.matched_genome_count,
                "truncated": listing.truncated,
                "genomes": [
                    {
                        "sample_id": genome.sample_id,
                        "collection_genome_ordinal": genome.collection_genome_ordinal,
                        # ⛔ TWO counts, and they are different facts: loci where the genome has a
                        # gene, and loci where it sits in an arrangement — the second is what the
                        # published dropdown printed. They differ for every probe genome.
                        "locus_count": genome.locus_count,
                        "arrangement_locus_count": genome.arrangement_locus_count,
                    }
                    for genome in listing.genomes
                ],
            },
            pangenome.pangenome_id,
        )


@species_blueprint.get("/catalogues/<species_key>/projected-genomes")
@species_blueprint.get("/species/<species_key>/projected-genomes")
def get_projected_genomes(species_key: str):
    """⭐ Genomes placed on this catalogue after the model was built — never members of it.

    ⛔ **Every count on this page stays the modelled 100.** These rows are an addition drawn beside
    the model: prevalence, bands, the census, the marginals, arrangement genome counts and the
    genome picker are unchanged by anything here, and the `no-count-change` test asserts it.

    ⚠ The response carries the RULE and both neighbour numbers, because a placement means nothing
    without them: the nearest modelled gene decided the locus, `neighbours_reported` of the nearest
    checked it, and `neighbours_searched` were kept so a wider check can be tried without searching
    again.
    """
    with _session() as session:
        try:
            pangenome = _resolve_pangenome(session, species_key)
        except SpeciesNotPublished as error:
            return _not_found(str(error))

        rows = list_projected_genomes(session, pangenome.pangenome_id)
        payload = {
            "species_key": species_key,
            "genomes": [serialise_projected_genome(row) for row in rows],
            # ⛔ The sentence the page must repeat wherever these genomes appear.
            "caveat": (
                "placed after the model was built: each gene sits on the locus of its nearest "
                "modelled gene, which approximates what the model would have done and is not the "
                "model's own clustering"
            ),
        }
        return _immutable(payload, pangenome.pangenome_id)
