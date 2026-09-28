"""The pangenome layer, written — loci, lists, both neighbourhood views, gaps and geometry.

⛔ **One pangenome is the unit of work.** Every table here is scoped to `pangenome_id` and replaced
wholesale, because 889k loci with 49M arrangements cannot be diffed row-by-row cheaply and a
half-applied upsert leaves a catalogue that is partly two models — the worst possible state.

⭐ **The derived columns are computed here and never at request time.** `display_name`, `best_product`,
`search_text`, the separation percentiles and `interest_score` are the fan-out fix: rendering one
locus touches a median 15–19 other loci because a neighbour's *name* is a transitive read through its
Pfam and product lists. Materialised, that becomes `WHERE locus_id = ANY($1)` over ~20 ids.

⚠ **`landing_locus_id` and `example_locus_ids` are set after the loci exist**, in the same
transaction. They were a render-time payload mutation, which is how every page deployed
2026-08-20..2026-08-25 opened on locus 0 with no reference names — a caller serialised before
`render` had mutated, and nothing failed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from sqlalchemy import select, update
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.catalogue_frames import (
    OFFSETS,
    SIMILARITY_COLUMNS,
    CatalogueFrames,
    build_catalogue_frames,
    load_pfam_name_table,
)
from bacatlas_backend.ingest.derive_locus_display import best_product, display_name, search_text
from bacatlas_backend.ingest.derive_locus_ranking import (
    build_interest_inputs,
    interest_score,
    landing_index,
    pfam_concordance,
    ranking,
    similarity_index,
)
from bacatlas_backend.ingest.ingest_genome_locus_counts import write_genome_locus_counts
from bacatlas_backend.ingest.staging_table_loader import copy_rows
from bacatlas_backend.models.enumerations import (
    AnnotationKind,
    EmbeddingRepresentation,
    GeneOntologyAgreementVerdict,
    PrevalenceBand,
)
from bacatlas_backend.models.gene import GeneLocusMembership
from bacatlas_backend.models.intergenic_gap import IntergenicGap, IntergenicGapFeature
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry, LocusUnirefFamilyCrosstab
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant
from bacatlas_backend.models.locus_similarity import (
    LocusNearestLocus,
    LocusSimilarity,
    PangenomeSimilarityBaseline,
)
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pangenome_genome_locus_count import PangenomeGenomeLocusCount
from bacatlas_backend.models.projected_genome import ProjectedGenome

#: `catalogue_map._COS` — the factor the 15 upper-triangle cosines are stored at.
COSINE_SCALE_FACTOR = 10_000

#: `catalogue_map._Q` — the half-extent the quantised map positions occupy, leaving the int16 ends
#: unused so a rounding step can never wrap a point to the opposite edge of the map.
QUANTISATION_HALF_EXTENT = 32_500


@dataclass
class CatalogueLoadReport:
    """What the pangenome layer wrote, per table, plus the coverage it examined."""

    pangenome_id: int
    loci: int = 0
    annotation_entries: int = 0
    uniref_crosstab_rows: int = 0
    arrangements: int = 0
    offset_occupants: int = 0
    gene_memberships: int = 0
    genes_without_a_window: int = 0
    genome_locus_counts: int = 0
    intergenic_gaps: int = 0
    intergenic_gap_features: int = 0
    similarity_baselines: int = 0
    similarity_rows: int = 0
    nearest_locus_rows: int = 0
    landing_locus_label: str | None = None
    example_locus_labels: list[str] = field(default_factory=list)
    display_name_sources: dict = field(default_factory=dict)
    separation_measurable: dict = field(default_factory=dict)
    coverage: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def render(self) -> str:
        """The report as the CLI prints it — every table, then what it was derived from."""
        lines = [
            f"pangenome layer (pangenome_id={self.pangenome_id})",
            f"  loci                    {self.loci:,}",
            f"  annotation entries      {self.annotation_entries:,}",
            f"  uniref50 cross-tab      {self.uniref_crosstab_rows:,}",
            f"  arrangements            {self.arrangements:,}",
            f"  offset occupants        {self.offset_occupants:,}",
            f"  gene memberships        {self.gene_memberships:,}"
            + (f"  ({self.genes_without_a_window:,} with no window)" if self.genes_without_a_window else ""),
            f"  genome locus counts     {self.genome_locus_counts:,} genomes",
            f"  intergenic gaps         {self.intergenic_gaps:,} "
            f"({self.intergenic_gap_features:,} named features)",
            f"  similarity baselines    {self.similarity_baselines} "
            f"({self.similarity_rows:,} locus rows, {self.nearest_locus_rows:,} nearest-locus rows)",
            f"  landing locus           {self.landing_locus_label}",
            f"  example loci            {', '.join(self.example_locus_labels)}",
            f"  display names           {self.display_name_sources}",
            f"  separation measurable   {self.separation_measurable}",
        ]
        for key, value in sorted(self.coverage.items()):
            lines.append(f"  read {key:<22s} {value}")
        for note in self.notes:
            lines.append(f"  NOTE {note}")
        return "\n".join(lines)


def _clean(value):
    """A pandas cell → a Python scalar or None, without a NaN reaching a `measurement()` column."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    item = getattr(value, "item", None)
    if item is not None and type(value).__module__.startswith("numpy"):
        return _clean(item())
    return value


def _int_or_none(value):
    value = _clean(value)
    return None if value is None else int(value)


def _int_or_zero(value):
    """A count that is zero when nothing was found, because the looking DID happen."""
    value = _clean(value)
    return 0 if value is None else int(value)


# ── loci ───────────────────────────────────────────────────────────────────────────────────────
def _derive_locus_columns(frames: CatalogueFrames, pfam_names: dict, policy: dict) -> dict:
    """Every column that exists only because the page computes it, computed once.

    ⚠ The annotation lists are grouped ONCE into per-locus dicts and reused by all four derivations.
    Re-scanning the long frame per locus would be quadratic over 17.5k loci and 72k rows.
    """
    loci = frames.loci
    nodes = [str(value) for value in loci["node"]]
    entries = frames.annotation_entries

    def rows_by_node(kind: str) -> dict[str, list[tuple[str, int]]]:
        if entries is None or not len(entries):
            return {}
        subset = entries[entries["annotation_kind"] == kind].sort_values(
            ["node", "rank"], kind="mergesort"
        )
        grouped: dict[str, list[tuple[str, int]]] = {}
        for node, term, count in zip(
            subset["node"], subset["term"], subset["gene_count"], strict=True
        ):
            grouped.setdefault(str(node), []).append((str(term), int(count)))
        return grouped

    architectures = rows_by_node("pfam_architecture")
    products = rows_by_node("protein_product")
    families: dict[str, list[str]] = {}
    crosstab = frames.uniref_crosstab
    if crosstab is not None and len(crosstab):
        ordered = crosstab.sort_values(["node", "rank"], kind="mergesort")
        for node, accession in zip(ordered["node"], ordered["uniref50"], strict=True):
            families.setdefault(str(node), []).append(str(accession))

    symbols = [_clean(value) for value in loci["bakta_gene_symbol"]]
    names, sources, accessions, best_products, haystacks = [], [], [], [], []
    for node, symbol in zip(nodes, symbols, strict=True):
        resolved = display_name(
            node, symbol, architectures.get(node, []), products.get(node, []), pfam_names
        )
        names.append(resolved.name)
        sources.append(resolved.source)
        accessions.append(resolved.source_accession)
        best_products.append(best_product(symbol, products.get(node, [])))
        haystacks.append(search_text(node, symbol, products.get(node, []), families.get(node, [])))

    concordance = pfam_concordance(
        [[count for _, count in architectures.get(node, [])] for node in nodes],
        [_int_or_none(value) for value in loci["pfam_annotated_member_count"]],
        [int(value) for value in loci["member_gene_count"]],
    )
    scores = [
        interest_score(row, policy)
        for row in build_interest_inputs(
            syntenic_a5=[_clean(value) for value in loci["syntenic_a5"]],
            uniref50_major_family_count=[int(value) for value in loci["uniref50_major_family_count"]],
            collapse_tier=[_clean(value) for value in loci["collapse_tier"]],
            resolved_threshold=[_clean(value) for value in loci["resolved_threshold"]],
            prevalence_band_index=[
                PREVALENCE_BAND_ORDER.index(str(value)) for value in loci["prevalence_band"]
            ],
            member_gene_count=[int(value) for value in loci["member_gene_count"]],
            concordance=concordance,
            pfam_concordance_class=[_clean(value) for value in loci["pfam_concordance_class"]],
        )
    ]
    # ⛔ Three midranks per representation now, not one. The three views are NOT three pictures of
    # one thing: only 18–34 % of flagged loci are flagged by all three (PROJECT_STATE Layer 6), so a
    # single percentile would be a rank for one view and a decoration on the other two.
    similarities = {
        representation: similarity_index(
            *(
                [_clean(v) for v in block["loci"][f"{stem}_{block['form']}"]]
                for stem in SIMILARITY_COLUMNS.values()
            )
        )
        for representation, block in frames.similarity.items()
    }
    return {
        "display_name": names,
        "display_name_source": sources,
        "display_name_source_accession": accessions,
        "best_product": best_products,
        "search_text": haystacks,
        "interest_score": scores,
        "similarities": similarities,
    }


#: `export_payload.BAND_ORDER` — the index `_interest` reads, so it must be that order and not the
#: enum's declaration order.
PREVALENCE_BAND_ORDER = ["core", "soft_core", "shell", "cloud", "rare"]

LOCUS_COLUMNS = (
    "pangenome_id",
    "pathogen_species_id",
    "node_label",
    "catalogue_ordinal",
    "member_gene_count",
    "member_genome_count",
    "prevalence_band",
    "display_name",
    "display_name_source",
    "display_name_source_accession",
    "best_product",
    "bakta_gene_symbol",
    "named_member_count",
    "uniref50_family_count",
    "uniref50_major_family_count",
    "uniref50_labelled_member_count",
    "pfam_annotated_member_count",
    "pfam_architecture_count",
    "pfam_concordance_class",
    "syntenic_a5",
    "collapse_tier",
    "collapse_bucket",
    "resolved_threshold",
    "seqid_coverage",
    "uniref50_impurity",
    "uniref50_coverage",
    "medoid_genome_id",
    "medoid_flat_index",
    "cog_annotated_member_count",
    "cog_distinct_id_count",
    "modal_cog_categories",
    "ec_annotated_member_count",
    "kegg_annotated_member_count",
    "go_annotated_member_count_molecular_function",
    "go_annotated_member_count_biological_process",
    "go_annotated_member_count_cellular_component",
    "go_verdict_molecular_function",
    "go_verdict_biological_process",
    "go_verdict_cellular_component",
    "median_gene_length_nt",
    "gene_length_interquartile_range_nt",
    "context_observed_member_counts",
    "total_arrangement_count",
    "arrangement_member_gene_count",
    "arrangement_member_genome_count",
    "interest_score",
    "search_text",
)


def _locus_rows(
    frames: CatalogueFrames,
    derived: dict,
    pangenome_id: int,
    species_id: int,
    genome_id_by_sample: dict,
):
    """One tuple per locus, in `LOCUS_COLUMNS` order.

    ⚠ `medoid_sample_id` is resolved to a `genome_id` HERE and not stored as an accession: the
    medoid is one real gene, and a column holding its genome as text would be a second, unjoinable
    spelling of a key the schema already has.
    """
    loci = frames.loci
    for index in range(len(loci)):
        row = loci.iloc[index]
        yield (
            pangenome_id,
            species_id,
            str(row["node"]),
            int(row["catalogue_ordinal"]),
            int(row["member_gene_count"]),
            int(row["member_genome_count"]),
            PrevalenceBand(str(row["prevalence_band"])),
            derived["display_name"][index],
            derived["display_name_source"][index],
            derived["display_name_source_accession"][index],
            derived["best_product"][index],
            _clean(row["bakta_gene_symbol"]),
            int(row["named_member_count"]),
            _int_or_none(row["uniref50_family_count"]),
            _int_or_none(row["uniref50_major_family_count"]),
            _int_or_none(row["uniref50_labelled_member_count"]),
            _int_or_none(row["pfam_annotated_member_count"]),
            _int_or_none(row["pfam_architecture_count"]),
            _clean(row["pfam_concordance_class"]),
            _clean(row["syntenic_a5"]),
            _clean(row["collapse_tier"]),
            _clean(row["collapse_bucket"]),
            _clean(row["resolved_threshold"]),
            _clean(row["seqid_coverage"]),
            _clean(row["uniref50_impurity"]),
            _clean(row["uniref50_coverage"]),
            genome_id_by_sample.get(str(_clean(row["medoid_sample_id"]) or "")),
            _int_or_none(row["medoid_flat_index"]),
            _int_or_none(row["cog_annotated_member_count"]),
            _int_or_none(row["cog_distinct_id_count"]),
            _cog_categories(row["modal_cog_category"]),
            _int_or_none(row["ec_annotated_member_count"]),
            _int_or_none(row["kegg_annotated_member_count"]),
            _int_or_none(row["go_annotated_member_count_molecular_function"]),
            _int_or_none(row["go_annotated_member_count_biological_process"]),
            _int_or_none(row["go_annotated_member_count_cellular_component"]),
            _verdict(row["go_verdict_molecular_function"]),
            _verdict(row["go_verdict_biological_process"]),
            _verdict(row["go_verdict_cellular_component"]),
            _int_or_none(row["median_gene_length_nt"]),
            _int_or_none(row["gene_length_interquartile_range_nt"]),
            list(row["context_observed_member_counts"]),
            int(row["total_arrangement_count"]),
            int(row["arrangement_member_gene_count"]),
            int(row["arrangement_member_genome_count"]),
            derived["interest_score"][index],
            derived["search_text"][index],
        )


def _cog_categories(value) -> list[str] | None:
    """`"EHJQ"` → `["E", "H", "J", "Q"]` — the SET Bakta concatenated, taken apart again.

    ⛔ Splitting is not cosmetic. The value fits a scalar column, so storing it whole fails nothing
    and makes `WHERE modal_cog_category = 'C'` miss every multi-category locus — 95 of 118 distinct
    values on the published *E. coli* catalogue. ⚠ An empty string yields `None`, not `[]`: no COG
    category is *not annotated*, which is different from *annotated with nothing*.
    """
    value = _clean(value)
    if value is None:
        return None
    letters = [character for character in str(value).strip() if character.isalpha()]
    return letters or None


def _verdict(value) -> GeneOntologyAgreementVerdict | None:
    """A GO agreement verdict, or None.

    ⛔ `no_coverage` is a VALUE, never a NULL — fewer than two annotated members is neither
    agreement nor disagreement, and must not be counted as either.
    """
    value = _clean(value)
    return None if value is None else GeneOntologyAgreementVerdict(str(value))


# ── the load ───────────────────────────────────────────────────────────────────────────────────
def _delete_pangenome_layer(session: Session, pangenome_id: int) -> None:
    """Every table scoped to this pangenome, in FK order. The blast radius is one catalogue."""
    for model, column in (
        (PangenomeGenomeLocusCount, PangenomeGenomeLocusCount.pangenome_id),
        (GeneLocusMembership, GeneLocusMembership.pangenome_id),
        (IntergenicGapFeature, None),
        (IntergenicGap, IntergenicGap.pangenome_id),
        (LocusOffsetOccupant, LocusOffsetOccupant.pangenome_id),
        (LocusArrangement, LocusArrangement.pangenome_id),
        (PangenomeSimilarityBaseline, PangenomeSimilarityBaseline.pangenome_id),
    ):
        if column is None:
            session.query(IntergenicGapFeature).filter(
                IntergenicGapFeature.intergenic_gap_id.in_(
                    select(IntergenicGap.intergenic_gap_id).where(
                        IntergenicGap.pangenome_id == pangenome_id
                    )
                )
            ).delete(synchronize_session=False)
            continue
        session.query(model).filter(column == pangenome_id).delete(synchronize_session=False)
    # ⛔ **The projection goes too, and it is announced.** `projected_gene_placement.locus_id`
    # cascades from `locus`, so re-ingesting a catalogue would silently empty the placements while
    # leaving `projected_genome` rows whose counts describe nothing — a summary saying "4,676 genes
    # placed on 4,652 loci" over a table with no rows in it. Deleting both together and SAYING so is
    # the only version of this a reader can act on: the placements are nuna artifacts and are
    # re-loadable with `--stage projection`, but nothing else would tell you they had gone.
    projected = session.query(ProjectedGenome).filter(
        ProjectedGenome.pangenome_id == pangenome_id
    ).delete(synchronize_session=False)
    if projected:
        print(
            f"  ⚠ {projected} projected genome(s) deleted with this catalogue's loci — "
            "re-run `--stage projection` to place them again"
        )

    # `locus` last: annotation entries, the cross-tab and the geometry cascade from it.
    session.query(Locus).filter(Locus.pangenome_id == pangenome_id).delete(synchronize_session=False)
    session.flush()


def ingest_locus_catalogue(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    pangenome_id: int,
    pathogen_species_id: int,
    genome_id_by_ordinal: list[int],
    frames: CatalogueFrames | None = None,
) -> CatalogueLoadReport:
    """Load the whole pangenome layer for one catalogue. One transaction, one blast radius."""
    frames = frames if frames is not None else build_catalogue_frames(artifacts)
    policy = _policy()
    derived = _derive_locus_columns(frames, load_pfam_name_table(), policy)

    _delete_pangenome_layer(session, pangenome_id)
    report = CatalogueLoadReport(pangenome_id=pangenome_id, coverage=dict(frames.coverage))

    report.loci = copy_rows(
        session,
        Locus.__table__,
        LOCUS_COLUMNS,
        _locus_rows(
            frames, derived, pangenome_id, pathogen_species_id,
            _genome_ids_by_sample(session, frames.samples),
        ),
    )
    session.flush()

    locus_id_by_label = {
        label: locus_id
        for label, locus_id in session.execute(
            select(Locus.node_label, Locus.locus_id).where(Locus.pangenome_id == pangenome_id)
        ).all()
    }
    locus_id_by_ordinal = [None] * len(frames.loci)
    for label, ordinal in zip(frames.loci["node"], frames.loci["catalogue_ordinal"], strict=True):
        locus_id_by_ordinal[int(ordinal)] = locus_id_by_label[str(label)]

    report.annotation_entries = _load_annotation_entries(session, frames, locus_id_by_label)
    report.uniref_crosstab_rows = _load_uniref_crosstab(session, frames, locus_id_by_label)
    arrangement_ids = _load_arrangements(
        session, frames, pangenome_id, locus_id_by_label, genome_id_by_ordinal
    )
    report.arrangements = len(arrangement_ids)
    report.offset_occupants = _load_offset_occupants(
        session, frames, pangenome_id, locus_id_by_label, locus_id_by_ordinal
    )
    memberships, without_window = _load_gene_memberships(
        session, frames, pangenome_id, locus_id_by_label, arrangement_ids
    )
    report.gene_memberships, report.genes_without_a_window = memberships, without_window
    # ⭐ AFTER both the arrangements and the memberships: it counts the rows just written, so the
    # picker's "N loci" and the locus view cannot describe two different catalogues.
    report.genome_locus_counts = write_genome_locus_counts(session, pangenome_id)
    report.intergenic_gaps, report.intergenic_gap_features = _load_intergenic_gaps(
        session, frames, pangenome_id, locus_id_by_label
    )
    (
        report.similarity_baselines,
        report.similarity_rows,
        report.nearest_locus_rows,
    ) = _load_cluster_similarity(session, frames, pangenome_id, locus_id_by_label, derived)

    # ── the landing locus and the example chips ────────────────────────────────────────────────
    order = ranking(derived["interest_score"])
    symbols = [_clean(value) for value in frames.loci["bakta_gene_symbol"]]
    landing = landing_index(order, symbols)
    examples = order[:4]
    session.execute(
        update(Pangenome)
        .where(Pangenome.pangenome_id == pangenome_id)
        .values(
            landing_locus_id=locus_id_by_ordinal[landing],
            example_locus_ids=[locus_id_by_ordinal[index] for index in examples],
            locus_count=report.loci,
            omitted_sections=frames.omitted or None,
        )
    )
    session.flush()

    report.landing_locus_label = str(frames.loci["node"].iloc[landing])
    report.example_locus_labels = [str(frames.loci["node"].iloc[index]) for index in examples]
    counts: dict[str, int] = {}
    for source in derived["display_name_source"]:
        counts[source] = counts.get(source, 0) + 1
    report.display_name_sources = counts
    report.separation_measurable = {
        key: index.measurable_count for key, index in derived["similarities"].items()
    }
    for section, reason in frames.omitted.items():
        report.notes.append(f"omitted {section}: {reason}")
    return report


def _policy() -> dict:
    """`export_payload.POLICY` — the failure/rescue tiers and the contested Pfam verdict."""
    from nuna.tl.locus_browser.export_payload import POLICY

    return POLICY


# ── the per-table loaders ──────────────────────────────────────────────────────────────────────
def _load_annotation_entries(session, frames, locus_id_by_label) -> int:
    entries = frames.annotation_entries
    if entries is None or not len(entries):
        return 0
    columns = (
        "locus_id",
        "annotation_kind",
        "rank_within_locus",
        "term_value",
        "term_name",
        "member_gene_count",
        "gene_ontology_namespace",
    )
    names = _reference_names()

    def rows():
        for node, kind, rank, term, count, namespace in zip(
            entries["node"],
            entries["annotation_kind"],
            entries["rank"],
            entries["term"],
            entries["gene_count"],
            entries["namespace"],
            strict=True,
        ):
            kind_value = str(kind)
            yield (
                locus_id_by_label[str(node)],
                AnnotationKind(kind_value),
                int(rank),
                str(term),
                names.get(kind_value, {}).get(str(term)),
                int(count),
                _int_or_none(namespace),
            )

    return copy_rows(session, LocusAnnotationEntry.__table__, columns, rows())


def _reference_names() -> dict:
    """The vendored name tables, joined into the ROW that uses them.

    ⛔ Joined in, not shipped as a lookup table — which is what makes the render-time `attach_*`
    mutation (five days of pages printing accessions with no names) structurally impossible.
    ⚠ **KEGG has no entry and never will**: its terms permit linking freely but not redistributing
    its content, and ~880 KO descriptions in a served page is redistribution.
    """
    from nuna.tl.locus_browser.render_page import ASSETS
    from nuna.tl.locus_browser.vendor_cog import cog_reference
    from nuna.tl.locus_browser.vendor_go import go_reference

    cog = cog_reference(ASSETS / "cog_names.tsv.gz")
    go = go_reference(ASSETS / "go_names.tsv.gz")
    return {
        "cog_orthogroup": {key: str(value[0]) for key, value in cog.items()},
        "gene_ontology_slim": {key: str(value[0]) for key, value in go.items()},
    }


def _load_uniref_crosstab(session, frames, locus_id_by_label) -> int:
    crosstab = frames.uniref_crosstab
    if crosstab is None or not len(crosstab):
        return 0
    columns = (
        "locus_id",
        "rank_within_locus",
        "uniref50_accession",
        "member_gene_count",
        "modal_bakta_product",
        "modal_pfam_architecture",
        "pfam_annotated_member_count",
        "modal_bakta_gene_symbol",
        "distinct_real_symbol_count",
        "named_member_count",
    )

    def rows():
        for row in crosstab.itertuples(index=False):
            yield (
                locus_id_by_label[str(row.node)],
                int(row.rank),
                str(row.uniref50)[:64],
                int(row.n),
                _clean(getattr(row, "product", None)),
                _clean(getattr(row, "arch", None)),
                # ⚠ `fillna(0)`, matching the payload: the family's genes WERE looked at and none
                # carried a domain, which is a measured zero and not an absence.
                _int_or_zero(getattr(row, "n_pfam", None)),
                _clean(getattr(row, "gene", None)),
                _int_or_zero(getattr(row, "n_sym", None)),
                # ⚠ A measured zero for a family no gene of which is named — the whole point of
                # the column. Same `fillna(0)` reasoning as `n_pfam` two lines up.
                _int_or_zero(getattr(row, "n_named", None)),
            )

    return copy_rows(session, LocusUnirefFamilyCrosstab.__table__, columns, rows())


def _load_arrangements(
    session, frames, pangenome_id, locus_id_by_label, genome_id_by_ordinal
) -> dict:
    """The joint view. Returns `{(node, rank): locus_arrangement_id}` for the per-gene link.

    ⛔ `member_genome_ids` is resolved from the collection ORDINAL to `genome_id` here, so nothing
    downstream ever holds a positional index that a re-export would renumber.
    """
    arrangements = frames.arrangements
    if arrangements is None or not len(arrangements):
        return {}
    slot_columns = [f"s{index}" for index in range(len(OFFSETS))]
    columns = (
        "locus_id",
        "pangenome_id",
        "rank_within_locus",
        "member_gene_count",
        "member_genome_count",
        "is_recorded_reverse_complement",
        "neighbour_slot_codes",
        "member_genome_ids",
    )

    def rows():
        for row in arrangements.itertuples(index=False):
            yield (
                locus_id_by_label[str(row.node)],
                pangenome_id,
                int(row.rank),
                int(row.genes),
                int(row.genomes),
                bool(row.flip),
                [int(getattr(row, name)) for name in slot_columns],
                [int(genome_id_by_ordinal[int(value)]) for value in row.gset],
            )

    copy_rows(session, LocusArrangement.__table__, columns, rows())
    session.flush()
    return {
        (label, rank): arrangement_id
        for label, rank, arrangement_id in session.execute(
            select(
                Locus.node_label,
                LocusArrangement.rank_within_locus,
                LocusArrangement.locus_arrangement_id,
            )
            .join(Locus, Locus.locus_id == LocusArrangement.locus_id)
            .where(LocusArrangement.pangenome_id == pangenome_id)
        ).all()
    }


def _load_offset_occupants(
    session, frames, pangenome_id, locus_id_by_label, locus_id_by_ordinal
) -> int:
    """The marginal view. ⚠ `neigh_pos` is a payload ORDINAL and is resolved to a `locus_id` here."""
    occupants = frames.offset_occupants
    if occupants is None or not len(occupants):
        return 0
    columns = (
        "locus_id",
        "pangenome_id",
        "signed_offset",
        "rank_within_offset",
        "neighbour_locus_id",
        "member_gene_count",
        "same_strand_member_count",
    )

    def rows():
        for row in occupants.itertuples(index=False):
            yield (
                locus_id_by_label[str(row.node)],
                pangenome_id,
                int(row.d),
                int(row.rank),
                locus_id_by_ordinal[int(row.neigh_pos)],
                int(row.cnt),
                int(row.same),
            )

    return copy_rows(session, LocusOffsetOccupant.__table__, columns, rows())


def _load_gene_memberships(
    session, frames, pangenome_id, locus_id_by_label, arrangement_ids
) -> tuple[int, int]:
    """One row per gene, with the arrangement it realises — `(written, without a window)`.

    ⚠ **A gene with no window gets a row with a NULL arrangement, not no row.** *"This genome has no
    gene at this locus"* and *"this genome has a gene here with no recorded neighbourhood"* are
    different sentences on the page, and only a row can carry the second.
    """
    memberships = frames.gene_memberships
    if memberships is None or not len(memberships):
        return 0, 0
    genome_id_by_sample = _genome_ids_by_sample(session, frames.samples)
    missing = [sample for sample in frames.samples if sample not in genome_id_by_sample]
    if missing:
        raise KeyError(
            f"{len(missing)} of {len(frames.samples)} roster genomes have no `genome` row "
            f"({missing[:5]}). Every membership is keyed by genome_id, so a missing one would "
            "silently drop that genome's genes from every locus it is in."
        )

    columns = ("pangenome_id", "genome_id", "flat_index", "locus_id", "locus_arrangement_id")
    without_window = 0

    def rows():
        nonlocal without_window
        for row in memberships.itertuples(index=False):
            rank = row.arrangement_rank
            missing = rank is None or (isinstance(rank, float) and math.isnan(rank))
            if missing:
                without_window += 1
            yield (
                pangenome_id,
                genome_id_by_sample[str(row.sample_id)],
                int(row.flat_index),
                locus_id_by_label[str(row.node)],
                None if missing else arrangement_ids.get((str(row.node), int(rank))),
            )

    written = copy_rows(session, GeneLocusMembership.__table__, columns, rows())
    return written, without_window


def _genome_ids_by_sample(session, samples: list[str]) -> dict:
    from bacatlas_backend.models.genome import Genome

    return {
        sample_id: genome_id
        for sample_id, genome_id in session.execute(
            select(Genome.sample_id, Genome.genome_id).where(Genome.sample_id.in_(samples))
        ).all()
    }


def _load_intergenic_gaps(session, frames, pangenome_id, locus_id_by_label) -> tuple[int, int]:
    """The gaps, keyed on the canonical pair — which is by **node LABEL**, not by any integer.

    ⛔ `gene_adjacencies` sorts the pair with `left <= right` on the node ids, and node ids are TEXT.
    The published page keys its lookup by payload INDEX instead, which is why 7,379 of 22,838 *E.
    coli* gaps cannot be found on it. `a`/`b` here are the labels, in the label order, so the API
    resolves a pair by sorting two labels and the miss cannot recur.

    ⚠ **The variance and the mode are stored DENSELY**, from `gap_table` itself. The payload's sparse
    `vi`/`vd`/`vmd` triple existed only because a JSON array cannot be sparse, and there absence from
    `vi` meant *agreement* (86.1 % of gaps) while absence of `vmd` meant the mode was simply not
    carried. Read from the source both are present for every gap, so `length_variance_score` is 0.0
    where every genome agrees and `modal_length_nt` is never NULL for want of an index entry.
    """
    gaps = frames.gaps
    if gaps is None or not len(gaps):
        return 0, 0
    columns = (
        "pangenome_id",
        "flanking_locus_id_a",
        "flanking_locus_id_b",
        "observed_genome_count",
        "median_signed_length_nt",
        "quartile1_signed_length_nt",
        "quartile3_signed_length_nt",
        "minimum_signed_length_nt",
        "maximum_signed_length_nt",
        "length_variance_score",
        "modal_length_nt",
        "distinct_named_feature_count",
    )
    ordered = list(gaps.itertuples(index=False))
    for row in ordered:
        if str(row.a) > str(row.b):
            raise ValueError(
                f"gap ({row.a}, {row.b}) is not in label order. `gene_adjacencies` canonicalises on "
                "the node ids as TEXT, and every reader downstream inherits that ordering."
            )

    def rows():
        for row in ordered:
            yield (
                pangenome_id,
                locus_id_by_label[str(row.a)],
                locus_id_by_label[str(row.b)],
                int(row.n),
                int(row.nt),
                int(row.q1),
                int(row.q3),
                int(row.mn),
                int(row.mx),
                _clean(row.v),
                _int_or_none(row.md),
                int(row.n_feat),
            )

    written = copy_rows(session, IntergenicGap.__table__, columns, rows())
    session.flush()

    gap_ids = {
        (a, b): gap_id
        for a, b, gap_id in session.execute(
            select(
                IntergenicGap.flanking_locus_id_a,
                IntergenicGap.flanking_locus_id_b,
                IntergenicGap.intergenic_gap_id,
            ).where(IntergenicGap.pangenome_id == pangenome_id)
        ).all()
    }
    feature_columns = (
        "intergenic_gap_id",
        "rank_within_gap",
        "feature_label",
        "feature_type",
        "observed_genome_count",
    )

    def feature_rows():
        for row in ordered:
            key = (locus_id_by_label[str(row.a)], locus_id_by_label[str(row.b)])
            gap_id = gap_ids[key]
            for rank, (label, feature_type, count) in enumerate(
                zip(row.labels, row.types, row.counts, strict=True)
            ):
                yield (gap_id, rank, str(label)[:256], str(feature_type)[:64], int(count))

    features = copy_rows(session, IntergenicGapFeature.__table__, feature_columns, feature_rows())
    return written, features


def _load_cluster_similarity(session, frames, pangenome_id, locus_id_by_label, derived) -> tuple[int, int, int]:
    """The per-representation baseline, the per-locus similarities, and the ranked nearest loci.

    ⛔ **What this replaced on 2026-09-24.** ``_load_map_geometry`` wrote a UMAP over every locus's
    MEDOID, that medoid's 6×6 local geometry, and the whole-catalogue dust rendered from the same
    quantised arrays. All three described a construction that reduced a locus to one member, and they
    went together because the picture was a picture *of* the numbers.

    ⚠ **Rows, not an ordinal array, for the nearest loci.** The retired ``nearest_locus_ordinals``
    packed five catalogue ordinals with ``-1`` for absent, and every reader had to resolve them; a
    ``-1`` read as an index names locus 0, which is real and looks entirely plausible. A neighbour
    outside the catalogue now simply has no row, and each row carries its own similarity.
    """
    similarity = frames.similarity
    if not similarity:
        return 0, 0, 0

    baselines = 0
    similarity_rows = 0
    nearest_rows = 0
    for representation, block in similarity.items():
        index = derived["similarities"][representation]
        floor = block["floor"]

        # ⛔ The denominator the card prints — *"of 12,104 loci"* — is asserted against the artifact's
        # own count rather than trusted. They are two derivations of one number: the `.meta` counts
        # loci with a within-locus pair, and `similarity_index` counts loci with a rankable value.
        declared = floor["measurable_locus_count"]
        if declared is not None and declared != index.measurable_count:
            raise ValueError(
                f"{representation}: the artifact declares {declared:,} measurable loci and the "
                f"ranking found {index.measurable_count:,}. The card prints one of these under a "
                "percentile computed from the other, so they cannot differ."
            )
        session.add(
            PangenomeSimilarityBaseline(
                pangenome_id=pangenome_id,
                representation=EmbeddingRepresentation(representation),
                similarity_form=block["form"],
                floor_median=floor["median"],
                floor_p25=floor["p25"],
                floor_p75=floor["p75"],
                floor_p99=floor["p99"],
                measurable_locus_count=index.measurable_count,
                neighbour_knn_k=floor["knn_k"],
                source_csv_path=block["source"],
            )
        )
        baselines += 1

        loci = block["loci"]
        form = block["form"]
        labels = [str(value) for value in loci["node"]]
        column_of = {name: f"{stem}_{form}" for name, stem in SIMILARITY_COLUMNS.items()}
        values = {name: [_clean(v) for v in loci[column]] for name, column in column_of.items()}

        similarity_columns = (
            "locus_id",
            "representation",
            *SIMILARITY_COLUMNS,
            "separation_percentile",
            "weak_margin_percentile",
            "own_fraction_percentile",
        )

        def similarity_rows_for(labels=labels, values=values, index=index, representation=representation):
            for position, label in enumerate(labels):
                yield (
                    locus_id_by_label[label],
                    EmbeddingRepresentation(representation),
                    *(values[name][position] for name in SIMILARITY_COLUMNS),
                    index.separation_percentile[position],
                    index.weak_margin_percentile[position],
                    index.own_fraction_percentile[position],
                )

        similarity_rows += copy_rows(
            session, LocusSimilarity.__table__, similarity_columns, similarity_rows_for()
        )

        # ── the five nearest, ranked ───────────────────────────────────────────────────────────
        nearest = block["nearest"]
        nearest_columns = ("locus_id", "representation", "rank", "neighbour_locus_id", "cross_similarity")
        head_by_label: dict[str, float] = {}

        def nearest_rows_for(nearest=nearest, representation=representation, head=head_by_label):
            for label, rank, neighbour, cross in zip(
                nearest["node"], nearest["rank"], nearest["neighbour"], nearest["cross"], strict=True
            ):
                label, neighbour = str(label), str(neighbour)
                # ⚠ A neighbour outside this catalogue is DROPPED, not stored as a sentinel. The
                # artifact is computed over the model's own loci, but a hand-trimmed catalogue may
                # hold fewer, and an unresolvable id has no honest placeholder in a foreign key.
                if label not in locus_id_by_label or neighbour not in locus_id_by_label:
                    continue
                if int(rank) == 1:
                    head[label] = float(cross)
                yield (
                    locus_id_by_label[label],
                    EmbeddingRepresentation(representation),
                    int(rank),
                    locus_id_by_label[neighbour],
                    _clean(cross),
                )

        nearest_rows += copy_rows(
            session, LocusNearestLocus.__table__, nearest_columns, nearest_rows_for()
        )

        # ⛔ `nearest_similarity` IS rank 1 of the list, stored twice so the card's headline renders
        # without reading the list. Asserted rather than hoped: the two files are written together by
        # one run, so a disagreement can only mean they came from two — and that is the exact failure
        # the card being replaced had no way to notice.
        disagreed = [
            label
            for position, label in enumerate(labels)
            if label in head_by_label
            and values["nearest_similarity"][position] is not None
            and abs(head_by_label[label] - values["nearest_similarity"][position]) > 1e-9
        ]
        if disagreed:
            raise ValueError(
                f"{representation}: {len(disagreed):,} loci (e.g. {disagreed[:3]}) whose "
                "`near` column disagrees with rank 1 of their own nearest-loci list — the two "
                "artifacts are not two views of one run"
            )
    return baselines, similarity_rows, nearest_rows
