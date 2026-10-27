"""The species list and one species' census — everything the page needs before a locus.

⭐ **`default_pangenome_id` decides what a species serves when the URL names no model.** Ingest
builds a new generation alongside the live one and never touches it; `publish_pangenome` flips it in
its own transaction. That is what makes a re-ingest safe on a running service: the new catalogue is
complete and verified before anything points at it, and a rollback is the same one-row update
backwards.

⛔ **It is a DEFAULT, not the only catalogue a species has.** A species may hold several at once —
nuna4 beside nuna5, or a sensitive model beside a less-sensitive one — and each is addressed
directly by `pangenome.catalogue_key` (`ecoli-nuna5`). `resolve_catalogue` is that path;
`resolve_published_pangenome` is the bare-species path and returns the default. Keep them distinct:
collapsing them is how a link that pins a model quietly stops pinning it.

⚠ **The census is read, not counted per request.** `pangenome.genome_count` / `gene_count` /
`locus_count` are written by ingest and reconciled against the checked-in published triple, so the
number the page prints is the number the loader wrote and not a `COUNT(*)` that could disagree with
it after a partial write.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import (
    EvaluationKind,
    PrevalenceBand,
)
from bacatlas_backend.models.genome_collection import GenomeCollection
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_similarity import PangenomeSimilarityBaseline
from bacatlas_backend.models.pangenome import Pangenome, PangenomeEvaluation
from bacatlas_backend.models.pathogen_species import PathogenSpecies

#: The audit headline keys the footer prints, in `export_payload.AUDIT_HEADLINE_KEYS` order.
#: ⛔ A whitelist, and copied verbatim — *"the page quoting a lookalike it derived itself is the
#: failure mode this exists to make impossible"*.
AUDIT_HEADLINE_KEYS = (
    "n_clusters_total", "n_genes_total", "n_singleton_clusters",
    "synteny_only_n_clusters", "synteny_only_n_genes", "synteny_only_gene_rate",
    "no_homology_n_clusters", "no_homology_n_genes", "no_homology_gene_rate",
    "over_merge_gene_rate", "over_merge_gene_rate_num",
    "split_gene_rate", "split_gene_rate_excl_singletons",
    "family_split_across_clusters_n_families",
    "pfam_conflict_n_clusters", "pfam_conflict_n_genes", "pfam_conflict_gene_rate",
    "pfam_judgeable_n_clusters", "pfam_no_coverage_n_clusters",
    "n_clusters_esm_rescued", "esm_rescued_gene_num",
)  # fmt: skip


class SpeciesNotPublished(LookupError):
    """The species exists, or does not, and either way it has nothing to serve."""


class CatalogueNotFound(SpeciesNotPublished):
    """No catalogue carries this key.

    ⚠ Subclasses `SpeciesNotPublished` deliberately, so every route that already turns that into a
    404 keeps doing so. It is a distinct class only so a caller CAN tell "no such catalogue" from
    "this species serves nothing" when it wants to.
    """


@dataclass
class SpeciesCatalogue:
    """One species' published pangenome, and everything the shell renders before a locus."""

    species: PathogenSpecies
    pangenome: Pangenome
    #: ⭐ **THE COHORT — the browser's axis since 2026-10-27, replacing `model`.** The picker offers
    #: cohorts and the model is one free-text label in the footer (David: *"not explained in any
    #: more detail than that. Just allows cross reference."*).
    #:
    #: ⛔ **Required, not Optional**, even though the old `model` was nullable. `pangenome.
    #: genome_collection_id` is NOT NULL with a foreign key, so a catalogue without a cohort cannot
    #: exist — and an `| None` here would invite an `if catalogue.cohort:` branch whose false arm is
    #: unreachable and therefore untestable, which is the shape `a-fixture-can-collapse-the-
    #: distinction-it-tests` warns about. The nullability belonged to `nuna_model_id`, not here.
    cohort: GenomeCollection
    prevalence_census: dict = field(default_factory=dict)
    #: ⭐ The SAME catalogue partitioned by GENE — member genes summed per band. The page prints both
    #: (`app.js::renderPangenome`) because they answer different questions: a third of the loci are
    #: singletons, while a typical genome carries about fifty singleton genes. Given only the locus
    #: count a reader reasonably concludes the wrong thing about what is in a genome.
    prevalence_gene_census: dict = field(default_factory=dict)
    audit_headline: dict = field(default_factory=dict)
    similarity_baselines: list = field(default_factory=list)
    #: ⚠ Descriptors ONLY — never the bytes. The species response is JSON and the sprite is a
    #: megabyte of PNG; they travel by different routes on purpose.
    landing_locus_label: str | None = None
    example_locus_labels: list = field(default_factory=list)
    #: The example loci as the chips draw them — label, name and the UniRef50 families the chip
    #: quotes (`render_page._examples`). A label alone would have to be resolved per chip.
    example_locus_rows: list = field(default_factory=list)


def list_published_species(session: Session) -> list[tuple[PathogenSpecies, Pangenome | None]]:
    """Every species and the pangenome it serves — `published.tsv`, as a query.

    ⚠ A species with no published pangenome is RETURNED with `None`, not filtered out. The picker
    lists the others and says why they are unavailable; a silently short list would read as though
    the species did not exist.
    """
    rows = session.execute(
        select(PathogenSpecies, Pangenome)
        .outerjoin(Pangenome, Pangenome.pangenome_id == PathogenSpecies.default_pangenome_id)
        .order_by(PathogenSpecies.species_key)
    ).all()
    return [(species, pangenome) for species, pangenome in rows]


def resolve_published_pangenome(session: Session, species_key: str) -> Pangenome:
    """The pangenome a species serves — ONE statement, for every route that is not the shell.

    ⛔⛔ **Not `load_species_catalogue(...).pangenome`.** That is the whole shell — the band census
    (a GROUP BY over every locus of the species), the audit headline, the projections, the sprite
    descriptors and the example chips — and every locus navigation, search keystroke and sequence
    view paid for it just to learn one integer. The cost oracle measured the SERVICES and so could
    not see it: at 17,531 loci it is a few milliseconds per click, at the 889,160-locus design target
    it is an aggregation over the whole catalogue on the hot path.
    """
    row = session.execute(
        select(PathogenSpecies.default_pangenome_id, Pangenome)
        .outerjoin(Pangenome, Pangenome.pangenome_id == PathogenSpecies.default_pangenome_id)
        .where(PathogenSpecies.species_key == species_key)
    ).one_or_none()
    if row is None:
        raise SpeciesNotPublished(f"no species {species_key!r}")
    if row.Pangenome is None:
        raise SpeciesNotPublished(
            f"{species_key!r} has no published pangenome. ⚠ That is a distinct state from 'no such "
            "species', and the picker says so rather than omitting it."
        )
    return row.Pangenome


def resolve_catalogue(session: Session, catalogue_key: str) -> Pangenome:
    """The pangenome a catalogue key names — ONE statement, and no species indirection.

    ⛔ This is the path that makes a link PIN a model. `resolve_published_pangenome` answers "what
    does this species serve by default"; this answers "give me exactly this catalogue", which is what
    a URL carrying `ecoli-nuna5` is asking for. A reader who shared that link and a reader who opens
    it must see the same clustering even after the default moves.
    """
    pangenome = session.execute(
        select(Pangenome).where(Pangenome.catalogue_key == catalogue_key)
    ).scalar_one_or_none()
    if pangenome is None:
        raise CatalogueNotFound(
            f"no catalogue {catalogue_key!r}. ⚠ Keys are `{{species}}-{{model}}` with a HYPHEN — "
            "`ecoli-nuna5`, never `ecoli_nuna5`."
        )
    return pangenome


@dataclass(frozen=True)
class OfferedCatalogue:
    """One row of the picker: the catalogue, the species it belongs to, and the COHORT it is over.

    ⛔ **A dataclass and not a tuple, deliberately, because the change it carries is silent
    otherwise.** `list_catalogues` used to return `(Pangenome, PathogenSpecies, NunaModel | None)`,
    and every one of its eleven call sites destructured it as `for p, _s, _m in …` while reading
    only the first two. Swapping the third element's TYPE would therefore have changed the picker's
    axis from the model to the cohort with every test still green — the same shape as a diff loop
    that `continue`s past a renamed column and reports "0 differ". A dataclass is not iterable, so
    each of those sites had to be read.
    """

    pangenome: Pangenome
    species: PathogenSpecies
    cohort: GenomeCollection


def list_catalogues(session: Session) -> list[OfferedCatalogue]:
    """Every catalogue offered in the picker, with the species and the COHORT that identify it.

    ⚠ Filtered on `pangenome.is_published`, which means *offer this one* — NOT
    `default_pangenome_id`, which is one per species. A catalogue that is loaded but not offered is
    reachable by its key and simply absent from the menu, which is how a catalogue is staged before
    it is shown to anyone.

    ⛔ **The collection join is INNER, where the model join was OUTER, and that is a statement not
    an oversight.** `pangenome.genome_collection_id` is NOT NULL with a foreign key, so the two are
    provably equivalent on any row that can exist; OUTER would copy a nullability that belonged to
    `nuna_model_id` and would quietly serve a picker option with no cohort name. ⚠ INNER's own
    failure mode is a catalogue silently missing from the menu, which is why
    `test_the_picker_offers_EVERY_published_catalogue` counts the rows against `is_published`
    directly rather than trusting the join.
    """
    rows = session.execute(
        select(Pangenome, PathogenSpecies, GenomeCollection)
        .join(PathogenSpecies, PathogenSpecies.pathogen_species_id == Pangenome.pathogen_species_id)
        .join(
            GenomeCollection,
            GenomeCollection.genome_collection_id == Pangenome.genome_collection_id,
        )
        .where(Pangenome.is_published.is_(True))
        .order_by(PathogenSpecies.species_key, Pangenome.catalogue_key)
    ).all()
    return [OfferedCatalogue(pangenome, species, cohort) for pangenome, species, cohort in rows]


def load_catalogue(session: Session, catalogue_key: str) -> SpeciesCatalogue:
    """The whole shell for ONE named catalogue — the model-pinned twin of `load_species_catalogue`."""
    pangenome = resolve_catalogue(session, catalogue_key)
    species = session.get(PathogenSpecies, pangenome.pathogen_species_id)
    return _load_catalogue_shell(session, species, pangenome)


def load_species_catalogue(session: Session, species_key: str) -> SpeciesCatalogue:
    """The whole species shell in a fixed number of statements, none of them per locus.

    Serves the species' DEFAULT catalogue. `load_catalogue` serves a named one.
    """
    species = session.execute(
        select(PathogenSpecies).where(PathogenSpecies.species_key == species_key)
    ).scalar_one_or_none()
    if species is None:
        raise SpeciesNotPublished(f"no species {species_key!r}")
    if species.default_pangenome_id is None:
        raise SpeciesNotPublished(
            f"{species_key!r} has no published pangenome. ⚠ That is a distinct state from 'no such "
            "species', and the picker says so rather than omitting it."
        )
    pangenome = session.get(Pangenome, species.default_pangenome_id)
    return _load_catalogue_shell(session, species, pangenome)


def _load_catalogue_shell(
    session: Session, species: PathogenSpecies, pangenome: Pangenome
) -> SpeciesCatalogue:
    """The shell itself, once the catalogue is decided. Both entry points share it exactly."""
    # ⭐ One `session.get` on the primary key, replacing TWO statements: the `nuna_model` fetch and
    # the `pangenome_step` ordering. The shell is a statement cheaper than it was, not dearer.
    # ⚠ `pangenome_step` is still WRITTEN — it is what the run actually did, with wall time and peak
    # RSS — it is simply no longer served, because no component ever read `catalogue.steps`.
    cohort = session.get(GenomeCollection, pangenome.genome_collection_id)
    catalogue = SpeciesCatalogue(species=species, pangenome=pangenome, cohort=cohort)

    # ⚠ The census IS counted, because it is a distribution rather than a total and no column holds
    # it. One grouped statement over an index, not one per band.
    # ⚠ The gene sums ride in the SAME grouped statement: a second one would be a second read of the
    # same rows, and the cost oracle counts statements.
    counts = session.execute(
        select(Locus.prevalence_band, func.count(), func.sum(Locus.member_gene_count))
        .where(Locus.pangenome_id == pangenome.pangenome_id)
        .group_by(Locus.prevalence_band)
    ).all()
    catalogue.prevalence_census = {band.value: count for band, count, _ in counts}
    catalogue.prevalence_gene_census = {band.value: int(genes or 0) for band, _, genes in counts}
    for band in PrevalenceBand:
        # ⛔ A band with no loci is `0`, not absent: the page prints every band and an absent key
        # would render as a gap rather than as the measured zero it is.
        catalogue.prevalence_census.setdefault(band.value, 0)
        catalogue.prevalence_gene_census.setdefault(band.value, 0)

    catalogue.audit_headline = {
        row.metric_name: row.numeric_value if row.numeric_value is not None else row.detail
        for row in session.execute(
            select(PangenomeEvaluation).where(
                PangenomeEvaluation.pangenome_id == pangenome.pangenome_id,
                PangenomeEvaluation.evaluation_kind == EvaluationKind.ACCESSORY_AUDIT,
                PangenomeEvaluation.metric_name.in_(AUDIT_HEADLINE_KEYS),
            )
        ).scalars()
    }

    catalogue.similarity_baselines = list(
        session.execute(
            select(PangenomeSimilarityBaseline).where(
                PangenomeSimilarityBaseline.pangenome_id == pangenome.pangenome_id
            )
        ).scalars()
    )

    wanted = [pangenome.landing_locus_id, *(pangenome.example_locus_ids or [])]
    rows = {
        row.locus_id: row
        for row in session.execute(
            select(
                Locus.locus_id,
                Locus.node_label,
                Locus.display_name,
                Locus.uniref50_family_count,
                Locus.uniref50_major_family_count,
            ).where(Locus.locus_id.in_([value for value in wanted if value is not None]))
        ).all()
    }
    landing = rows.get(pangenome.landing_locus_id)
    catalogue.landing_locus_label = landing.node_label if landing else None
    examples = [rows[locus_id] for locus_id in (pangenome.example_locus_ids or []) if locus_id in rows]
    catalogue.example_locus_labels = [row.node_label for row in examples]
    catalogue.example_locus_rows = [
        {
            "label": row.node_label,
            "display_name": row.display_name,
            # ⚠ The MAJOR family count, as the chip always quoted (`nodes.n_u50_major` over `n_u50`):
            # a locus with one real family and a scatter of one-member families is not the point
            # the chip is making. Falls back to the full count where no major count was measured.
            "uniref50_family_count": (
                row.uniref50_major_family_count
                if row.uniref50_major_family_count is not None
                else row.uniref50_family_count
            ),
        }
        for row in examples
    ]
    return catalogue


# ⛔ `load_scatter_sprite_descriptors` and `load_scatter_sprite` were DELETED on 2026-09-24 with the
# whole-catalogue sprite they served: a pre-rendered PNG of every locus's MEDOID, which went with the
# medoid geometry it was a picture of. The route, its ETag handling and the restore drill's check on it
# went in the same commit.
#
# ⚠ Worth keeping if a catalogue-wide picture ever returns: the reason those two were separate functions
# is that the descriptor query NAMED ITS COLUMNS ONE BY ONE. `select(Model)` would have pulled a megabyte
# of PNG into the species response for every page load, in order to send a dozen numbers.
