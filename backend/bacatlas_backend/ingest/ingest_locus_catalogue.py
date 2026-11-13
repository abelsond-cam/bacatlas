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
from pathlib import Path

from sqlalchemy import ARRAY, Integer, String, any_, literal, select, update
from sqlalchemy.orm import Session

from bacatlas_backend.gene_blocks import (
    NO_ARRANGEMENT_RANK,
    UNASSIGNED_LOCUS_ORDINAL,
    GeneBlockError,
    encode_int32_vector,
)
from bacatlas_backend.gene_blocks.nseq_gene_table import loci_from_bytes
from bacatlas_backend.ingest.artifact_locator import (
    DEFAULT_READ_WORKERS,
    CatalogueArtifacts,
    listed_files,
    read_ahead_in_order,
)
from bacatlas_backend.ingest.catalogue_frames import (
    OFFSETS,
    SIMILARITY_COLUMNS,
    CatalogueFrameError,
    CatalogueFrames,
    build_catalogue_frames,
    load_pfam_name_table,
)
from bacatlas_backend.ingest.derive_locus_display import best_product, display_name, search_text
from bacatlas_backend.ingest.derive_locus_propagation import PROPAGATION_VOCABULARIES
from bacatlas_backend.ingest.derive_locus_ranking import (
    build_interest_inputs,
    interest_score,
    landing_index,
    pfam_concordance,
    ranking,
    similarity_index,
)
from bacatlas_backend.ingest.ingest_genome_locus_counts import write_genome_locus_counts
from bacatlas_backend.ingest.ingest_propagation_rates import write_propagation_rates
from bacatlas_backend.ingest.phase_report import phase
from bacatlas_backend.ingest.staging_table_loader import copy_rows
from bacatlas_backend.models.enumerations import (
    GENE_ONTOLOGY_NAMESPACE_NAMES,
    AnnotationKind,
    EmbeddingRepresentation,
    GeneOntologyAgreementVerdict,
    PrevalenceBand,
    PropagationVerdict,
)
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
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
from bacatlas_backend.models.pangenome_propagation_rate import PangenomePropagationRate
from bacatlas_backend.models.projected_genome import ProjectedGenome
from bacatlas_backend.services.reference_name_service import gene_ontology_fold

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
    #: Membership rows in the frame — every one is encoded into its genome's assignment block.
    gene_memberships: int = 0
    genes_without_a_window: int = 0
    #: One `genome_locus_assignment` row per collection genome, written from the memberships frame.
    genome_locus_assignments: int = 0
    #: ⭐ The `.loci` oracle, as three numbers and never a percentage: how many genomes' encoded
    #: ordinals were compared byte-for-byte against a `.loci` nuna published, how many had no such
    #: file, and where the files were looked for (`None` = the oracle was not looked at at all).
    loci_oracle_checked: int = 0
    loci_oracle_absent: int = 0
    loci_oracle_root: str | None = None
    genome_locus_counts: int = 0
    propagation_rates: int = 0
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
            f"  genome assignments      {self.genome_locus_assignments:,} genomes  "
            + self._loci_oracle_line(),
            f"  genome locus counts     {self.genome_locus_counts:,} genomes",
            # ⚠ The denominator is printed with it: a vocabulary that was never measured writes no
            # rows at all, and a bare count gives a reader nothing to notice that against.
            f"  propagation rates       {self.propagation_rates:,} rows of "
            f"{len(PROPAGATION_VOCABULARIES) * len(PREVALENCE_BAND_ORDER)} (4 vocabularies x 5 bands)",
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

    def _loci_oracle_line(self) -> str:
        """*Not looked at* and *looked at, none found* are different findings, and both are said."""
        if self.loci_oracle_root is None:
            return "(.loci oracle NOT looked at: the artifacts name no nseq_root)"
        return (
            f"({self.loci_oracle_checked:,} checked byte-for-byte against a .loci, "
            f"{self.loci_oracle_absent:,} with no .loci under {self.loci_oracle_root})"
        )


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
    "present_genome_ids",
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
    "go_annotated_member_count_any",
    "go_claiming_member_count_molecular_function",
    "go_claiming_member_count_biological_process",
    "go_claiming_member_count_cellular_component",
    "propagation_verdict_cog",
    "propagation_verdict_ec",
    "propagation_verdict_kegg",
    "propagation_verdict_go",
    "median_gene_length_nt",
    "gene_length_interquartile_range_nt",
    "context_observed_member_counts",
    "total_arrangement_count",
    "arrangement_member_gene_count",
    "arrangement_member_genome_count",
    "interest_score",
    "search_text",
)


def _present_genome_ids_by_node(frames: CatalogueFrames, genome_id_by_sample: dict) -> dict:
    """`{node label: sorted distinct genome_ids with a member gene there}`, built in ONE pass.

    ⭐ From the memberships frame — the same rows `genome_locus_assignment` is written from — so the
    two cannot disagree about which genomes carry a locus. ⛔ Vectorised: the frame is 374 M rows at
    71,435 genomes, so this factorises node and sample once and takes the unique `(node, genome_id)`
    pairs with numpy; nothing here is per row in Python. A per-locus scan inside `_locus_rows` would
    be quadratic on top.
    """
    import numpy as np

    memberships = frames.gene_memberships
    if memberships is None:
        raise CatalogueFrameError(
            "`present_genome_ids` cannot be built: the frames carry no gene memberships, and the "
            "column is NOT NULL — a locus with no stated genomes would be a locus nothing can find."
        )
    if not len(memberships):
        return {}
    node_codes, node_uniques = _factorize(memberships["node"])
    label_of_code = [str(value) for value in node_uniques]
    genome_ids_of_code = _genome_ids_of_samples(memberships["sample_id"], genome_id_by_sample)
    sample_codes = genome_ids_of_code["codes"]
    genome_id_of_code = genome_ids_of_code["genome_ids"]
    # One int64 per membership: `node_code * (max genome_id + 1) + genome_id`, so `np.unique` sorts
    # by node first and by genome_id within a node — exactly the order the column wants.
    genome_ids = genome_id_of_code[sample_codes]
    base = int(genome_id_of_code.max()) + 1
    pairs = np.unique(node_codes.astype(np.int64) * base + genome_ids.astype(np.int64))
    pair_nodes = pairs // base
    pair_genomes = pairs % base
    boundaries = np.flatnonzero(np.diff(pair_nodes)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [len(pairs)]))
    return {
        label_of_code[int(pair_nodes[start])]: [int(value) for value in pair_genomes[start:end]]
        for start, end in zip(starts, ends, strict=True)
    }


def _factorize(values):
    """`(codes, uniques)` for a frame column of any dtype — pandas' C-speed hashing, never a loop."""
    import pandas as pd

    codes, uniques = pd.factorize(values, use_na_sentinel=True)
    if (codes < 0).any():
        raise CatalogueFrameError(
            f"{int((codes < 0).sum()):,} membership rows carry a null in a key column "
            "(`node` or `sample_id`); a gene with no locus or no genome has no row to be in."
        )
    return codes, list(uniques)


def _genome_ids_of_samples(samples, genome_id_by_sample: dict) -> dict:
    """Factorise `sample_id` and resolve each DISTINCT sample to its `genome_id` once.

    ⛔ Refuses by name a sample that has no `genome` row rather than dropping its genes: every
    membership is keyed by genome_id, so a missing one would silently empty that genome from every
    locus it is in.
    """
    import numpy as np

    codes, uniques = _factorize(samples)
    missing = [str(sample) for sample in uniques if str(sample) not in genome_id_by_sample]
    if missing:
        raise KeyError(
            f"{len(missing)} membership sample(s) have no `genome` row ({missing[:5]}). "
            "Every membership is keyed by genome_id, so a missing one would silently drop that "
            "genome's genes from every locus it is in."
        )
    return {
        "codes": codes,
        "samples": [str(sample) for sample in uniques],
        "genome_ids": np.asarray([genome_id_by_sample[str(sample)] for sample in uniques], dtype=np.int64),
    }


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

    ⛔ `present_genome_ids` is asserted against `member_genome_count` on EVERY locus and the first
    disagreement refuses the catalogue by node. The count comes from nuna's export and the array
    from the memberships frame; they are two derivations of one fact, and a locus whose array is
    shorter than its count would be a locus some of whose genomes nothing can find.
    """
    loci = frames.loci
    present_by_node = _present_genome_ids_by_node(frames, genome_id_by_sample)
    # ⛔ **Built BY NAME and projected through `LOCUS_COLUMNS`, never as a positional tuple.**
    # `copy_rows` validates that every name exists on the table (`staging_table_loader.py:125`) and
    # nothing validates arity or order, so a value inserted at the wrong index shifts every later
    # column one place left. Half of these are nullable integers, so the shift loads cleanly and
    # produces rollups that are plausible and wrong — `go_claiming_member_count_molecular_function`
    # and `_biological_process` would simply swap. The dict below cannot express that error, and the
    # key-set assertion is checked ONCE rather than per locus.
    for index in range(len(loci)):
        row = loci.iloc[index]
        label = str(row["node"])
        member_genome_count = int(row["member_genome_count"])
        present_genome_ids = present_by_node.get(label, [])
        if len(present_genome_ids) != member_genome_count:
            raise CatalogueFrameError(
                f"locus {label!r}: the memberships name {len(present_genome_ids):,} distinct "
                f"genome(s) but member_genome_count says {member_genome_count:,}. "
                "`present_genome_ids` is what answers 'which genomes carry this locus', so the two "
                "cannot differ."
            )
        values = {
            "pangenome_id": pangenome_id,
            "pathogen_species_id": species_id,
            "node_label": label,
            "catalogue_ordinal": int(row["catalogue_ordinal"]),
            "member_gene_count": int(row["member_gene_count"]),
            "member_genome_count": member_genome_count,
            "present_genome_ids": present_genome_ids,
            "prevalence_band": PrevalenceBand(str(row["prevalence_band"])),
            "display_name": derived["display_name"][index],
            "display_name_source": derived["display_name_source"][index],
            "display_name_source_accession": derived["display_name_source_accession"][index],
            "best_product": derived["best_product"][index],
            "bakta_gene_symbol": _clean(row["bakta_gene_symbol"]),
            "named_member_count": int(row["named_member_count"]),
            "uniref50_family_count": _int_or_none(row["uniref50_family_count"]),
            "uniref50_major_family_count": _int_or_none(row["uniref50_major_family_count"]),
            "uniref50_labelled_member_count": _int_or_none(row["uniref50_labelled_member_count"]),
            "pfam_annotated_member_count": _int_or_none(row["pfam_annotated_member_count"]),
            "pfam_architecture_count": _int_or_none(row["pfam_architecture_count"]),
            "pfam_concordance_class": _clean(row["pfam_concordance_class"]),
            "syntenic_a5": _clean(row["syntenic_a5"]),
            "collapse_tier": _clean(row["collapse_tier"]),
            "collapse_bucket": _clean(row["collapse_bucket"]),
            "resolved_threshold": _clean(row["resolved_threshold"]),
            "seqid_coverage": _clean(row["seqid_coverage"]),
            "uniref50_impurity": _clean(row["uniref50_impurity"]),
            "uniref50_coverage": _clean(row["uniref50_coverage"]),
            "medoid_genome_id": genome_id_by_sample.get(str(_clean(row["medoid_sample_id"]) or "")),
            "medoid_flat_index": _int_or_none(row["medoid_flat_index"]),
            "cog_annotated_member_count": _int_or_none(row["cog_annotated_member_count"]),
            "cog_distinct_id_count": _int_or_none(row["cog_distinct_id_count"]),
            "modal_cog_categories": _cog_categories(row["modal_cog_category"]),
            "ec_annotated_member_count": _int_or_none(row["ec_annotated_member_count"]),
            "kegg_annotated_member_count": _int_or_none(row["kegg_annotated_member_count"]),
            "go_annotated_member_count_any": _int_or_none(row["go_annotated_member_count_any"]),
            "propagation_verdict_cog": _propagation(row["propagation_verdict_cog"]),
            "propagation_verdict_ec": _propagation(row["propagation_verdict_ec"]),
            "propagation_verdict_kegg": _propagation(row["propagation_verdict_kegg"]),
            "propagation_verdict_go": _propagation(row["propagation_verdict_go"]),
            "median_gene_length_nt": _int_or_none(row["median_gene_length_nt"]),
            "gene_length_interquartile_range_nt": _int_or_none(row["gene_length_interquartile_range_nt"]),
            "context_observed_member_counts": list(row["context_observed_member_counts"]),
            "total_arrangement_count": int(row["total_arrangement_count"]),
            "arrangement_member_gene_count": int(row["arrangement_member_gene_count"]),
            "arrangement_member_genome_count": int(row["arrangement_member_genome_count"]),
            "interest_score": derived["interest_score"][index],
            "search_text": derived["search_text"][index],
        }
        # ⚠ Nine columns are `<prefix>_{molecular_function,biological_process,cellular_component}`
        # and they are written by LOOP rather than by hand. Three near-identical lines per prefix is
        # where a copy-paste files a biological process under cellular component — both nullable
        # integers, so Postgres accepts it and the rollup is plausible and wrong.
        for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES:
            values[f"go_annotated_member_count_{namespace}"] = _int_or_none(
                row[f"go_annotated_member_count_{namespace}"]
            )
            values[f"go_claiming_member_count_{namespace}"] = _int_or_none(
                row[f"go_claiming_member_count_{namespace}"]
            )
            values[f"go_verdict_{namespace}"] = _verdict(row[f"go_verdict_{namespace}"])
        if index == 0 and set(values) != set(LOCUS_COLUMNS):
            raise CatalogueFrameError(
                "the locus row builder and LOCUS_COLUMNS disagree: "
                f"only in the row {sorted(set(values) - set(LOCUS_COLUMNS))}, "
                f"only in the column list {sorted(set(LOCUS_COLUMNS) - set(values))}"
            )
        yield tuple(values[name] for name in LOCUS_COLUMNS)


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


def _propagation(value) -> PropagationVerdict | None:
    """A propagation verdict, or None where the vocabulary was never measured.

    ⛔ **It REFUSES anything else rather than returning None.** `_verdict` one function down coerces
    with `GeneOntologyAgreementVerdict(str(value))`, which raises on an unknown string — the right
    behaviour, and the reason a near-copy here that swallowed the unknown case would be dangerous:
    `locus.propagation_verdict_*` is nullable, so a coercer that returned None on everything it did
    not recognise would write NULL for all 17,531 loci and no test would fail. NULL has a meaning
    here (*this ingest never computed it*) and it must be reachable only on purpose.
    """
    if value is None or value != value:  # NaN compares unequal to itself
        return None
    if isinstance(value, PropagationVerdict):
        return value
    raise CatalogueFrameError(
        f"a propagation verdict reached the loader as {type(value).__name__} {value!r}; "
        "`derive_locus_propagation` yields PropagationVerdict or None and nothing else"
    )


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
        (PangenomePropagationRate, PangenomePropagationRate.pangenome_id),
        # References pangenome and genome only, so its place in the order is free.
        (GenomeLocusAssignment, GenomeLocusAssignment.pangenome_id),
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
    allow_missing_strand: bool = False,
    genomes_without_gene_table: frozenset[str] | None = None,
    read_workers: int = DEFAULT_READ_WORKERS,
) -> CatalogueLoadReport:
    """Load the whole pangenome layer for one catalogue. One transaction, one blast radius.

    ⭐ ``genomes_without_gene_table``, when given, is the NAMED set of collection genomes that hold no
    `genome_genes` row (`--genomes-without-gene-table`): the genomes the genome stage actually loaded
    without one must be exactly this set, or the load is refused naming the difference both ways.
    Absent, a genome without a gene table is still accepted wherever the genome stage recorded WHY
    (`genome.gene_table_omitted_reason`), and refused by name wherever it did not.
    """
    # ⛔ The fold comes from `gene_ontology_term`, which is what `function_inference_service` reads,
    # and NOT from nuna's `slim_of(go_reference())` that `catalogue_frames` already loads for
    # `go_verdicts` — that mapping does no upward closure and strips no root. Two folds for one
    # question is how a stored verdict and a served one disagree with nothing to say so.
    annotation_fold = gene_ontology_fold(session)
    frames = frames if frames is not None else build_catalogue_frames(
        artifacts, annotation_fold=annotation_fold, allow_missing_strand=allow_missing_strand
    )
    policy = _policy()
    derived = _derive_locus_columns(frames, load_pfam_name_table(), policy)
    phase("load: derived locus columns")

    _delete_pangenome_layer(session, pangenome_id)
    report = CatalogueLoadReport(pangenome_id=pangenome_id, coverage=dict(frames.coverage))

    # ⛔ Resolved ONCE and handed to both the locus rows and the memberships: `present_genome_ids`
    # and `genome_locus_assignment` are two readings of one frame, and two lookups could differ.
    genome_id_by_sample = _genome_ids_by_sample(session, frames.samples)
    missing_genomes = [sample for sample in frames.samples if sample not in genome_id_by_sample]
    if missing_genomes:
        raise KeyError(
            f"{len(missing_genomes)} of {len(frames.samples)} roster genomes have no `genome` row "
            f"({missing_genomes[:5]}). Every membership is keyed by genome_id, so a missing one would "
            "silently drop that genome's genes from every locus it is in."
        )
    # ⛔ Before a single locus row: the gene count of every collection genome, from its block or — for
    # a genome NAMED as holding none — from `genome.coding_gene_count`. Refuses by name on a gap.
    gene_counts = _gene_counts_by_genome(
        session, [genome_id_by_sample[sample] for sample in frames.samples], frames.samples
    )
    if genomes_without_gene_table is not None:
        _refuse_unless_named_set_matches(genomes_without_gene_table, gene_counts.omitted_reasons)

    report.loci = copy_rows(
        session,
        Locus.__table__,
        LOCUS_COLUMNS,
        _locus_rows(frames, derived, pangenome_id, pathogen_species_id, genome_id_by_sample),
    )
    session.flush()
    phase("load: locus rows")

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
    phase("load: annotation entries + crosstab")
    report.arrangements = _load_arrangements(
        session, frames, pangenome_id, locus_id_by_label, genome_id_by_ordinal
    )
    phase("load: arrangements")
    report.offset_occupants = _load_offset_occupants(
        session, frames, pangenome_id, locus_id_by_label, locus_id_by_ordinal
    )
    phase("load: offset occupants")
    # ⚠ `nseq_root` is being added to `CatalogueArtifacts` by the genome stage; until it lands the
    # attribute may not exist, and an artifact set without it means the `.loci` oracle is NOT looked
    # at — which the report says in those words rather than printing "0 of 0".
    loci_root = getattr(artifacts, "nseq_root", None)
    written = _load_gene_memberships(
        session,
        frames,
        pangenome_id,
        genome_id_by_sample,
        loci_root=Path(loci_root).expanduser() if loci_root is not None else None,
        gene_counts=gene_counts,
        read_workers=read_workers,
    )
    report.gene_memberships = written.memberships
    report.genes_without_a_window = written.without_window
    report.genome_locus_assignments = written.assignments
    report.loci_oracle_checked = written.loci_checked
    report.loci_oracle_absent = written.loci_absent
    report.loci_oracle_root = written.loci_root
    phase("load: genome locus assignments")
    # ⭐ AFTER both the arrangements and the assignments: it counts the rows just written, so the
    # picker's "N loci" and the locus view cannot describe two different catalogues.
    report.genome_locus_counts = write_genome_locus_counts(session, pangenome_id)
    phase("load: genome locus counts")
    report.propagation_rates = write_propagation_rates(
        session,
        pangenome_id=pangenome_id,
        bands=[str(band) for band in frames.loci["prevalence_band"]],
        verdicts={
            suffix: list(frames.loci[f"propagation_verdict_{suffix}"])
            for suffix in PROPAGATION_VOCABULARIES
        },
    )
    report.intergenic_gaps, report.intergenic_gap_features = _load_intergenic_gaps(
        session, frames, pangenome_id, locus_id_by_label
    )
    (
        report.similarity_baselines,
        report.similarity_rows,
        report.nearest_locus_rows,
    ) = _load_cluster_similarity(session, frames, pangenome_id, locus_id_by_label, derived)
    phase("load: gaps + similarity")

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
            omitted_sections=_omitted_sections(frames.omitted, gene_counts, len(frames.samples)),
            # ⚠ Names the fold the GO-folded columns were computed under, and ONLY those: the three
            # `go_annotated_member_count_*` beside them come from nuna's own vendored table read at
            # frame-build time, so this string cannot certify their vintage. It governs
            # `propagation_verdict_go` and `go_claiming_member_count_*`.
            annotation_fold_version=annotation_fold.version(),
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
    for section, reason in (_omitted_sections(frames.omitted, gene_counts, len(frames.samples)) or {}).items():
        report.notes.append(f"omitted {section}: {reason}")
    return report


def _omitted_sections(omitted: dict, gene_counts: GeneCounts, roster_size: int) -> dict | None:
    """`frames.omitted`, plus `gene_table` when any collection genome holds no gene table.

    ⭐ The footer renders every value as *"Not shown on this build: …"*, so the sentence names what is
    missing (the Sequence tab), for how many of the roster, why, and which genomes — the count alone
    would not let a reader find the one they anchored.
    """
    sections = dict(omitted)
    named = gene_counts.omitted_reasons
    if named:
        reasons = sorted(set(named.values()))
        sections["gene_table"] = (
            f"the gene table (Sequence tab) of {len(named):,} of {roster_size:,} genomes — "
            + "; ".join(reasons)
            + f" ({', '.join(sorted(named))}). They are clustered, assigned and counted like every "
            "other genome"
        )
    return sections or None


def _refuse_unless_named_set_matches(named: frozenset[str], omitted_reasons: dict[str, str]) -> None:
    """⛔ The genomes loaded without a gene table must be EXACTLY the named set — both differences named."""
    actual = set(omitted_reasons)
    unexpected, missing = sorted(actual - named), sorted(named - actual)
    if unexpected or missing:
        raise KeyError(
            "--genomes-without-gene-table does not match the database: "
            + "; ".join(
                part
                for part in (
                    f"{len(unexpected)} collection genome(s) hold no gene table but are not on the list: "
                    f"{unexpected[:10]}" if unexpected else "",
                    f"{len(missing)} listed genome(s) are not collection genomes without a gene table "
                    f"(they hold one, or are not in this collection): {missing[:10]}" if missing else "",
                )
                if part
            )
        )


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
) -> int:
    """The joint view. Returns how many arrangements were written.

    ⚠ It used to read every row back as `{(node, rank): locus_arrangement_id}` for a caller that took
    only `len()` of it — the per-gene consumer of the ids went with `gene_locus_membership` on
    2026-10-26. At 71,435 genomes that was tens of millions of rows fetched into a dict (~372 B and
    1.4 µs each, measured) held through every later phase, to learn a number `copy_rows` already
    returns. ⚠ The two are the same number: `(locus_id, rank_within_locus)` and
    `(pangenome_id, node_label)` are both unique and the pangenome's layer is deleted first, so the
    read-back held exactly the rows this COPY wrote.

    ⛔ `member_genome_ids` is resolved from the collection ORDINAL to `genome_id` here, so nothing
    downstream ever holds a positional index that a re-export would renumber.
    """
    arrangements = frames.arrangements
    if arrangements is None or not len(arrangements):
        return 0
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

    written = copy_rows(session, LocusArrangement.__table__, columns, rows())
    session.flush()
    return written


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


@dataclass
class GeneMembershipsWritten:
    """What `_load_gene_memberships` wrote and what it checked, each as a count a reader can size."""

    #: Membership rows read from the frame, all of them encoded into the assignment blocks.
    memberships: int = 0
    without_window: int = 0
    assignments: int = 0
    loci_checked: int = 0
    loci_absent: int = 0
    loci_root: str | None = None


@dataclass(frozen=True)
class GeneCounts:
    """Every collection genome's gene count, and which of them hold no gene table and why."""

    #: `{genome_id: gene_count}` — `genome_genes.gene_count`, or `genome.coding_gene_count` for a
    #: genome named as holding no gene table.
    by_genome: dict
    #: `{sample_id: genome.gene_table_omitted_reason}` for exactly the genomes with no block.
    omitted_reasons: dict


ASSIGNMENT_COLUMNS = (
    "pangenome_id",
    "genome_id",
    "gene_count",
    "placed_gene_count",
    "distinct_locus_count",
    "gene_locus_ordinals",
    "gene_arrangement_ranks",
)


def _load_gene_memberships(
    session,
    frames,
    pangenome_id,
    genome_id_by_sample,
    *,
    loci_root: Path | None = None,
    gene_counts: GeneCounts | None = None,
    read_workers: int = DEFAULT_READ_WORKERS,
) -> GeneMembershipsWritten:
    """The per-genome assignment blocks (`genome_locus_assignment`), from the memberships frame.

    ⚠ **A gene with no window is still placed.** *"This genome has no gene at this locus"* and
    *"this genome has a gene here with no recorded neighbourhood"* are different sentences on the
    page: the second is the gene's ordinal with `NO_ARRANGEMENT_RANK` as its rank, and it is counted
    here from the frame's own `arrangement_rank`.

    ⭐ Every genome of the collection gets an assignment row — one with no membership at all gets
    `gene_count` entries of `-1` and two zero counts, because the picker lists every genome and a
    missing row is a genome it cannot offer. ⛔ `gene_count` is read from `genome_genes`, never
    inferred from the memberships: the block is indexed by `flat_index` and a length guessed from
    the highest index seen would be short by exactly the unassigned tail. ⚠ A genome NAMED as holding
    no gene table (`genome.gene_table_omitted_reason`) takes `genome.coding_gene_count` — the `_meta`
    parquet's row count, checked as `0..n-1` at the genome stage — and is written like any other;
    it has no `.loci`, so it falls under the oracle's "absent" count.

    ⛔ **The catalogue ordinal has ONE source**, `frames.loci["catalogue_ordinal"]` — the same column
    `_locus_rows` writes to `locus.catalogue_ordinal` — so a `.loci` written by nuna and the block
    written here index the same numbering by construction.

    ⭐ **The `.loci` oracle.** Where `loci_root` is given and `<sample>.loci` exists, the encoded
    ordinals are compared byte-for-byte against nuna's file and the first mismatch refuses the
    catalogue naming the sample and the index. The counts of checked and absent genomes travel in
    the return value: a report that said only "checked" could mean 0 of 71,435.
    ⚠ Which `.loci` exist is ONE scan of `loci_root` (`listed_files`), never a cold `is_file()` per
    genome; the files are read by ``read_workers`` threads ahead of this generator and compared here,
    in collection order, byte for byte as before. A read that fails raises when its genome is reached.

    ⚠ **Vectorised, because the frame is 374 M rows at 71,435 genomes.** Node and sample are
    factorised once (C-speed hashing), the rows are sorted by genome once (an integer argsort), and
    each genome's block is built by numpy fancy indexing on its slice. The only Python loop is one
    iteration per GENOME.
    """
    import numpy as np

    memberships = frames.gene_memberships
    if memberships is None:
        return GeneMembershipsWritten()
    genome_ids_in_collection = [genome_id_by_sample[sample] for sample in frames.samples]
    if gene_counts is None:
        gene_counts = _gene_counts_by_genome(session, genome_ids_in_collection, frames.samples)
    gene_count_by_genome = gene_counts.by_genome

    # ── the frame as four aligned arrays ──────────────────────────────────────────────────────
    ordinals = _catalogue_ordinals_of(memberships["node"], frames.loci)
    resolved = _genome_ids_of_samples(memberships["sample_id"], genome_id_by_sample)
    sample_codes, sample_of_code = resolved["codes"], resolved["samples"]
    flat_indexes = np.asarray(memberships["flat_index"].to_numpy(), dtype=np.int64)
    ranks, without_window = _arrangement_ranks_of(memberships["arrangement_rank"])

    written = GeneMembershipsWritten(
        memberships=len(flat_indexes),
        without_window=without_window,
        loci_root=None if loci_root is None else str(loci_root),
    )

    # ── the per-genome blocks: one slice per genome after one sort ────────────────────────────
    # ⛔ Building the blocks is where a flat_index past the end, a repeated one, or a `.loci` that
    # disagrees is refused BY NAME.
    order = np.argsort(sample_codes, kind="stable")
    sorted_codes = sample_codes[order]
    starts = np.searchsorted(sorted_codes, np.arange(len(sample_of_code)), side="left")
    ends = np.searchsorted(sorted_codes, np.arange(len(sample_of_code)), side="right")
    slice_of_sample = {
        sample_of_code[code]: (int(starts[code]), int(ends[code])) for code in range(len(sample_of_code))
    }

    loci_present = frozenset() if loci_root is None else listed_files(loci_root, ".loci")

    def read_oracle(pair):
        """The `.loci`'s bytes (or the error reading them raised); `None` where there is no file."""
        sample = pair[0]
        if sample not in loci_present:
            return None
        try:
            return True, (loci_root / f"{sample}.loci").read_bytes()
        except Exception as error:  # noqa: BLE001 — re-raised where the comparison reads it
            return False, error

    def assignment_rows():
        in_order = read_ahead_in_order(
            zip(frames.samples, genome_ids_in_collection, strict=True),
            read_oracle,
            read_workers if loci_root is not None else 0,
        )
        for (sample, genome_id), oracle_read in in_order:
            gene_count = gene_count_by_genome[genome_id]
            start, end = slice_of_sample.get(sample, (0, 0))
            rows = order[start:end]
            block = _assignment_block(
                sample,
                gene_count,
                flat_indexes[rows],
                ordinals[rows],
                ranks[rows],
            )
            if loci_root is not None:
                if oracle_read is not None:
                    _refuse_unless_loci_agree(
                        loci_root / f"{sample}.loci", sample, gene_count, block["ordinals"], read=oracle_read
                    )
                    written.loci_checked += 1
                else:
                    written.loci_absent += 1
            yield (
                pangenome_id,
                genome_id,
                gene_count,
                block["placed_gene_count"],
                block["distinct_locus_count"],
                encode_int32_vector(block["ordinals"], name=f"{sample} gene_locus_ordinals"),
                encode_int32_vector(block["ranks"], name=f"{sample} gene_arrangement_ranks"),
            )

    written.assignments = copy_rows(
        session, GenomeLocusAssignment.__table__, ASSIGNMENT_COLUMNS, assignment_rows()
    )

    return written


def _gene_counts_by_genome(session, genome_ids: list, samples: list) -> GeneCounts:
    """Every collection genome's gene count — from its block, or for a NAMED genome from `genome`.

    ⛔ The genome stage writes `genome_genes`; a genome it has not reached has no block length, and
    inferring one from the memberships would silently shorten the block by its unassigned tail. So a
    genome with no block is accepted ONLY where the genome stage recorded why
    (`gene_table_omitted_reason` non-NULL), and then its count is `genome.coding_gene_count` — the
    `_meta` parquet's own row count. Every other gap is refused by name.

    ⚠ `= ANY(array)`, never `IN (…)`: at 71,435 genomes an `IN` list is 71,435 bind parameters and
    Postgres' wire protocol carries at most 65,535 — measured, psycopg refuses before sending.
    """
    rows = session.execute(
        select(
            Genome.genome_id,
            GenomeGenes.gene_count,
            Genome.gene_table_omitted_reason,
            Genome.coding_gene_count,
        )
        .outerjoin(GenomeGenes, GenomeGenes.genome_id == Genome.genome_id)
        .where(Genome.genome_id == any_(literal(list(genome_ids), ARRAY(Integer))))
    ).all()
    found = {genome_id: (block_count, reason, coding) for genome_id, block_count, reason, coding in rows}
    by_genome: dict = {}
    omitted_reasons: dict = {}
    missing: list[str] = []
    for sample, genome_id in zip(samples, genome_ids, strict=True):
        block_count, reason, coding = found.get(genome_id, (None, None, None))
        if block_count is not None:
            by_genome[genome_id] = block_count
        elif reason is not None and coding is not None:
            by_genome[genome_id] = coding
            omitted_reasons[sample] = reason
        else:
            missing.append(sample)
    if missing:
        raise KeyError(
            f"{len(missing)} of {len(samples)} collection genomes have no `genome_genes` row and are not "
            f"named as holding none ({missing[:5]}). `genome_locus_assignment.gene_count` must equal "
            "the block's, so the genome stage has to have written the block — or recorded why there is "
            "none (`--genomes-without-gene-table`) — before the pangenome stage can index it."
        )
    return GeneCounts(by_genome=by_genome, omitted_reasons=omitted_reasons)


def _catalogue_ordinals_of(nodes, loci):
    """Each membership's catalogue ordinal, from `frames.loci` — the ONE numbering `_locus_rows` writes."""
    import numpy as np

    ordinal_by_label = {
        str(label): int(ordinal)
        for label, ordinal in zip(loci["node"], loci["catalogue_ordinal"], strict=True)
    }
    codes, uniques = _factorize(nodes)
    unknown = [str(label) for label in uniques if str(label) not in ordinal_by_label]
    if unknown:
        raise CatalogueFrameError(
            f"{len(unknown):,} membership node(s) are not in the locus frame ({unknown[:5]}); "
            "a gene assigned to a locus the catalogue does not hold has no ordinal to carry."
        )
    ordinal_of_code = np.asarray([ordinal_by_label[str(label)] for label in uniques], dtype=np.int64)
    return ordinal_of_code[codes]


def _arrangement_ranks_of(ranks):
    """`arrangement_rank` as int32 with `NO_ARRANGEMENT_RANK` where the frame holds NaN/None.

    Returns `(ranks, how many had none)`. ⚠ A rank that is neither absent nor a non-negative whole
    number is refused: the column is a position in `locus_arrangement`, and `int()` of 1.5 is 1.
    """
    import numpy as np
    import pandas as pd

    values = pd.to_numeric(ranks, errors="raise").to_numpy(dtype="float64", na_value=np.nan)
    absent = np.isnan(values)
    present = values[~absent]
    if len(present) and ((present < 0) | (present != np.floor(present))).any():
        bad = present[(present < 0) | (present != np.floor(present))][:5]
        raise CatalogueFrameError(
            f"arrangement_rank holds values that are not non-negative whole numbers ({bad.tolist()}); "
            "a rank is a position in `locus_arrangement` and cannot be rounded into one."
        )
    out = np.where(absent, NO_ARRANGEMENT_RANK, values).astype(np.int32)
    return out, int(absent.sum())


def _assignment_block(sample: str, gene_count: int, flat_indexes, ordinals, ranks) -> dict:
    """One genome's two int32 vectors and its two counts, from its slice of the memberships.

    ⛔ Refuses by name a `flat_index` outside `[0, gene_count)` — which would either raise an index
    error a long way from the genome or, negative, silently write the LAST gene — and a repeated one,
    which would let one membership overwrite another with nothing to say so.
    """
    import numpy as np

    if len(flat_indexes):
        lowest, highest = int(flat_indexes.min()), int(flat_indexes.max())
        if lowest < 0 or highest >= gene_count:
            offender = lowest if lowest < 0 else highest
            raise CatalogueFrameError(
                f"{sample}: membership flat_index {offender} is outside [0, {gene_count}) — "
                f"`genome_genes.gene_count` says this genome has {gene_count:,} genes, so either the "
                "genome stage and the catalogue describe different annotations or the index is not "
                "a flat_index."
            )
        distinct_indexes, counts = np.unique(flat_indexes, return_counts=True)
        if (counts > 1).any():
            repeated = int(distinct_indexes[counts > 1][0])
            raise CatalogueFrameError(
                f"{sample}: flat_index {repeated} has {int(counts[counts > 1][0])} membership rows; a gene "
                "is in exactly one locus, and a second row would silently overwrite the first."
            )
    ordinal_vector = np.full(gene_count, UNASSIGNED_LOCUS_ORDINAL, dtype=np.int32)
    rank_vector = np.full(gene_count, NO_ARRANGEMENT_RANK, dtype=np.int32)
    ordinal_vector[flat_indexes] = ordinals
    rank_vector[flat_indexes] = ranks
    placed = ordinal_vector != UNASSIGNED_LOCUS_ORDINAL
    return {
        "ordinals": ordinal_vector,
        "ranks": rank_vector,
        "placed_gene_count": int(placed.sum()),
        "distinct_locus_count": int(len(np.unique(ordinal_vector[placed]))),
    }


def _refuse_unless_loci_agree(path: Path, sample: str, gene_count: int, ordinals, *, read=None) -> None:
    """Compare the block against nuna's `.loci` for this genome and name the first difference.

    ⛔ Byte-for-byte via `read_loci`'s own decode (`loci_from_bytes`), which already refuses a file
    whose length is not `gene_count * 4` — that refusal is re-raised naming the sample, because the
    file's own name says nothing about which `genome_genes` row it was measured against. ``read`` is
    `(True, bytes)` or `(False, error)` from a read done ahead; absent, the file is read here.
    """
    import numpy as np

    try:
        if read is None:
            raw = path.read_bytes()
        else:
            succeeded, raw = read
            if not succeeded:
                raise raw
        published = loci_from_bytes(raw, path, gene_count=gene_count)
    except GeneBlockError as error:
        raise CatalogueFrameError(
            f"{sample}: the published {path.name} does not hold {gene_count:,} genes "
            f"(`genome_genes.gene_count`): {error}"
        ) from error
    # `read_loci` returns a native-endian `array("i")`; numpy reads its buffer without a copy.
    expected = np.frombuffer(published, dtype=np.int32)
    differing = np.flatnonzero(expected != ordinals)
    if len(differing):
        first = int(differing[0])
        raise CatalogueFrameError(
            f"{sample}: gene_locus_ordinals disagrees with the published {path.name} at "
            f"{len(differing):,} of {gene_count:,} genes, first at flat_index {first} "
            f"(ours {int(ordinals[first])}, the file's {int(expected[first])}). The block and the "
            "file are two encodings of one assignment, so this catalogue was not built from the "
            "frame the file was."
        )


def _genome_ids_by_sample(session, samples: list[str]) -> dict:
    # ⚠ `= ANY(array)`, never `IN (…)`: the bind-parameter limit, as in `_gene_counts_by_genome`.
    return {
        sample_id: genome_id
        for sample_id, genome_id in session.execute(
            select(Genome.sample_id, Genome.genome_id).where(
                Genome.sample_id == any_(literal(list(samples), ARRAY(String)))
            )
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
