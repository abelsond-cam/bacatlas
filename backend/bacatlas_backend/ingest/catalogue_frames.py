"""One catalogue, computed from its source artifacts — the pandas half, with no database in it.

⛔ **This calls nuna's own functions and reimplements none of them.** `node_order`, `top_counts`,
`family_pfam`, `family_modal_symbols`, `family_modal_products`, `oriented_windows`, `arrangements`,
`neighbour_counts`, `near_synteny_agreement`, `gap_table` and even the two rounding helpers
(`_floats`) are the science; a second implementation of any of them would be a second
thing to keep in step, and the two would agree until the day they did not. What this module does is
**assemble** their outputs into rows rather than into the payload's columnar, string-interned,
run-length-encoded blocks.

⭐ That split is what makes the parity check mean something. Ingest and the export share the science
and differ in the assembly, so a comparison between the database and the published catalogue tests
the assembly — which is the part that is new.

⚠ **Read the SOURCE artifacts, never the published payload.** The payload is the acceptance oracle;
loading it would make the test circular. The one thing taken from outside the source tree is the
vendored Pfam/GO/COG name tables, which are public reference metadata and not a property of any model.

⛔ **The audit's own verdicts are READ, never re-derived** — `class_clan`, `n_arch`, the collapse
tier, the resolved threshold and the two embedding pairs. The page once re-derived the Pfam verdict
by a different rule and disagreed with the report on 2 of 22,624 loci: *"a page that quotes the
report must not be able to contradict it, ever."*
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from bacatlas_backend.ingest.allele_variant_symbols import fold_symbol_column
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.derive_locus_propagation import (
    PROPAGATION_VOCABULARIES,
    go_claiming_counts,
    propagation_verdicts,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import scalar_or_none, split_set_value
from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold

#: The signed offsets, in display order. ⛔ `0` is absent — it is the focal gene.
#: Asserted against `export_payload.OFFSETS` at load, because `locus_offset_occupant` has a CHECK
#: over this vocabulary and `context_observed_member_counts` is an array in exactly this order.
OFFSETS = (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)

#: `render_page`'s `nodes.a5` / `nodes.resolved` rounding, which `interest_score` must see.
SYNTENY_DECIMALS = 4
RESOLVED_DECIMALS = 3

#: The seven annotation vocabularies, as `(AnnotationKind value, payload list key, value column)`.
#: ⚠ `u50` is deliberately absent: it has five extra columns and gets its own cross-tab table.
ANNOTATION_SOURCES = (
    ("gene_symbol", "sym", "gene"),
    ("protein_product", "prod", "product"),
    ("pfam_architecture", "pfam", "arch"),
    ("cog_orthogroup", "cog", "cog_id"),
    ("ec_number", "ec", "ec"),
    ("kegg_orthology", "kegg", "kegg_ko"),
)


class CatalogueFrameError(RuntimeError):
    """A catalogue that cannot be assembled without inventing part of it."""


@dataclass
class CatalogueFrames:
    """Every table of one catalogue, as pandas frames in payload locus order."""

    #: One row per locus, ordered by `catalogue_ordinal`, carrying every `locus` column.
    loci: object
    #: `(node, annotation_kind, rank, term, gene_count, namespace)` — long, one row per term.
    annotation_entries: object
    #: `(node, rank, uniref50, n, product, arch, n_pfam, gene, n_sym)` — the cross-tab.
    uniref_crosstab: object
    #: `(node, rank, genes, genomes, flip, s0..s9, gset)` — the JOINT view.
    arrangements: object
    #: `(node, d, neigh_pos, cnt, same, obs, rank)` — the MARGINAL view.
    offset_occupants: object
    #: One row per member gene: `(sample_id, flat_index, node, arrangement_rank)`.
    gene_memberships: object
    #: `(a, b, …)` keyed by node LABEL, canonicalised the way `gene_adjacencies` canonicalises.
    gaps: object
    #: `{representation: {"info", "form", "source", "loci", "nearest", "floor"}}` — set-to-set similarity.
    similarity: dict

    #: The genome vocabulary `arr.gid` indexes into, and the ordinal every membership is a position in.
    samples: list = field(default_factory=list)
    n_genes: int = 0
    audit_summary: dict = field(default_factory=dict)
    #: Sections the source tree could not supply, so an absent block is never read as a measured zero.
    omitted: dict = field(default_factory=dict)
    #: ⭐ What was actually examined, so a report states its own coverage rather than saying "done".
    coverage: dict = field(default_factory=dict)


# ── loading ────────────────────────────────────────────────────────────────────────────────────
class MissingStrand(CatalogueFrameError):
    """A cohort whose strand sidecars are short — refused, because the damage is invisible."""


def _gate_strand(processed, samples: list[str], omitted: dict, *, allow_missing: bool) -> None:
    """⛔ Refuse a cohort whose `{BS}_strand.parquet` set is short, naming what is absent.

    `load_gene_coords` fills an absent sidecar with `+` for that genome's **whole gene set** and
    only *prints* a counter. Every track arrow on that genome then points forward, its windows are
    never reversed, and `syntenic_A5` — a mean over the cohort — is pulled toward 0.5 with no error
    anywhere. ⚠ **Partial coverage is worse than none**: one neighbourhood yields two different
    `s0..s9` vectors purely on which genomes happened to have a file.

    ⭐ The rule and both its checks are nuna's — `cohort_a5_run` already answers this for the same
    statistic — so this is the second CALLER of one gate, not a second gate. `strand_balance` is
    here because a sidecar that exists and reads `+` everywhere passes `strand_coverage` and then
    produces exactly the same uniform map.
    """
    from nuna.tl.probe.cohort_a5_run import strand_coverage

    n_present, missing = strand_coverage(processed, samples)
    if not missing:
        return
    head = ", ".join(missing[:5]) + (" …" if len(missing) > 5 else "")
    message = (
        f"{len(missing):,} of {len(samples):,} genomes have no `_strand.parquet` ({head}). "
        f"{n_present:,} present. load_gene_coords fills those with '+' for the genome's whole "
        "gene set, which un-reverses its windows, points every track arrow forward and pulls "
        "syntenic_A5 toward 0.5 with no error."
    )
    if not allow_missing:
        raise MissingStrand(message + " Refusing — pass --allow-missing-strand to proceed.")
    omitted["strand"] = message + " --allow-missing-strand given; PROCEEDING."
    print(f"⚠ {omitted['strand']}", flush=True)


def _gate_strand_balance(coords, omitted: dict, *, allow_missing: bool) -> None:
    """⚠ The SECOND silent-failure path: sidecars that all exist and all read `+`.

    Genes are roughly evenly split between strands in a real bacterial genome, so a cohort outside
    `cohort_a5_run.STRAND_BALANCE_BAND` is not a cohort — it is a map that was never read. This
    passes `strand_coverage` by construction, which is why it is a separate check and not a
    branch of it.
    """
    from nuna.tl.probe.cohort_a5_run import strand_balance

    verdict = strand_balance(coords)
    if not verdict["strand_degenerate"]:
        return
    message = (
        f"strand is degenerate: {verdict['share_forward']:.4f} of genes read '+' "
        f"({verdict['strand_counts']}). Every sidecar may be present and still carry one strand, "
        "which produces exactly the uniform map an absent sidecar does."
    )
    if not allow_missing:
        raise MissingStrand(message + " Refusing — pass --allow-missing-strand to proceed.")
    omitted["strand_balance"] = message + " --allow-missing-strand given; PROCEEDING."
    print(f"⚠ {omitted['strand_balance']}", flush=True)


def _load_source_frames(artifacts: CatalogueArtifacts, *, allow_missing_strand: bool = False) -> dict:
    """Every source frame `export_payload.main` loads, by the same calls and in the same order."""
    import numpy
    import pandas
    from nuna.eval.assignments import gene_universe, load_assignments
    from nuna.eval.reference_annotation import (
        load_bakta_dbxrefs,
        load_bakta_noncoding,
        load_bakta_products,
        load_gene_coords,
    )

    processed = artifacts.processed_root
    # ⚠ `close_seq=` and not the derived `{set}_close_seq.parquet`: the derived name is a naming
    # CONVENTION a big cohort does not keep, and the locator is the one place that knows which
    # file this catalogue's universe actually is.
    universe = gene_universe(processed, artifacts.set_key, close_seq=artifacts.gene_universe)
    samples = sorted(str(value) for value in universe["sample_id"].unique())

    omitted: dict[str, str] = {}
    _gate_strand(processed, samples, omitted, allow_missing=allow_missing_strand)
    assign = load_assignments(artifacts.assignment, universe=universe)[
        ["sample_id", "flat_index", "node"]
    ]
    products = load_bakta_products(processed, samples)
    coords = load_gene_coords(processed, samples)
    _gate_strand_balance(coords, omitted, allow_missing=allow_missing_strand)
    noncoding = load_bakta_noncoding(processed, samples)

    dbxrefs = load_bakta_dbxrefs(processed, samples)
    if not len(dbxrefs):
        # ⚠ A payload without COG is *a page without a function tab*, not a broken export — the same
        # judgement `build_payload` makes. Recorded so the absence is never read as a measured zero.
        omitted["uniref50"] = "no {BS}_dbxref.parquet under embeddings/meta/ — allele families unavailable"
        dbxrefs = None

    pfam = None
    if artifacts.pfam_architectures.exists():
        pfam = pandas.read_parquet(artifacts.pfam_architectures)
        pfam["sample_id"] = pfam["sample_id"].astype(str)
    else:
        omitted["pfam"] = f"{artifacts.pfam_architectures.name} not found — domain architectures unavailable"

    return {
        "numpy": numpy,
        "pandas": pandas,
        "samples": samples,
        "n_genes": int(len(universe)),
        "assign": assign,
        "products": products,
        "coords": coords,
        "noncoding": noncoding,
        "dbxrefs": dbxrefs,
        "pfam": pfam,
        "omitted": omitted,
    }


def _load_audit_frames(artifacts: CatalogueArtifacts, omitted: dict) -> dict:
    """The audit's own outputs — read, and recorded with WHICH artifact each came from.

    ⚠ The published *E. coli* export did **not** read `{label}_cluster_table.parquet` (it did not
    exist on CSD3 that day), so its Bacformer pair was recomputed from the full embedding matrix.
    Measured against the shipped catalogue afterwards the two agree to 5.0e-5 — the payload's own
    4-significant-figure rounding — on all 12,104 measurable loci. So the cluster table is a sound
    source, and it is a *different artifact* than the page used: which one supplied each column is
    recorded rather than left for a reader to assume.
    """
    import pandas
    from nuna.tl.locus_browser.export_payload import load_audit_evidence

    evidence, source = load_audit_evidence(artifacts.cluster_table, artifacts.homology_waterfall)
    if evidence is None:
        raise CatalogueFrameError(
            f"neither {artifacts.cluster_table.name} nor {artifacts.homology_waterfall.name} exists. "
            "The collapse tier, the resolved threshold and both embedding pairs come from one of "
            "them, and a catalogue without those is not the one the pages serve."
        )

    concordance = None
    if artifacts.pfam_concordance.exists():
        concordance = pandas.read_csv(artifacts.pfam_concordance, sep="\t")
        concordance["node"] = concordance["node"].astype(str)
        concordance = concordance.drop_duplicates("node").set_index("node")
    else:
        omitted["pfam_concordance"] = (
            f"{artifacts.pfam_concordance.name} not found — the Pfam verdict would fall back to a "
            "derivation over the top-5 architectures, which can differ from the audit's over all"
        )

    return {
        "evidence": evidence.set_index("node"),
        "evidence_source": source,
        "concordance": concordance,
        "audit_summary": json.loads(artifacts.audit_summary.read_text()),
    }


def load_pfam_name_table() -> dict:
    """The vendored Pfam table — `render_page._pfam_lookup`, the only reader of that file.

    ⛔ **`display_name`'s second fall-through needs this and the exported payload does not carry
    it.** `attach_pfam_names` runs at RENDER time, so an ingest reading only the export would demote
    370 *E. coli* loci to `locus <label>` — which reads as a data gap rather than a missing table.
    """
    from nuna.tl.locus_browser.render_page import _pfam_lookup

    return _pfam_lookup()


# ── the wide annotation frame ──────────────────────────────────────────────────────────────────
def _annotation_frame(assign, products, dbxrefs, pfam, numpy):
    """`build_payload`'s one wide frame — every per-locus reduction below is a groupby on it."""
    from nuna.eval.reference_annotation import DBXREF_VALUE_COLS
    from nuna.tl.locus_browser.export_payload import real_gene_names

    function_columns = [column for column in DBXREF_VALUE_COLS if column != "sample_id"]
    names = real_gene_names(products)
    # ⭐ The ONE place a catalogue symbol is decided, so the fold has one home: the name vote, the
    # locus symbol list and the cross-tab's modal symbol all read this column.
    # `gene.bakta_gene_symbol` is loaded by the genome layer from its own source and keeps the RAW
    # allele — see `allele_variant_symbols`, which also says why `named_member_count` cannot move.
    names = names.assign(gene=fold_symbol_column(names["gene"]))
    frame = assign.merge(names, on=["sample_id", "flat_index"], how="left").merge(
        products[["sample_id", "flat_index", "product"]], on=["sample_id", "flat_index"], how="left"
    )
    if dbxrefs is not None:
        # ⚠ `reindex` and not a column select: a frame from before a column existed yields NaN for
        # it rather than a KeyError, which is what keeps an older mirror loadable.
        frame = frame.merge(
            dbxrefs.reindex(columns=["sample_id", "flat_index", *function_columns]),
            on=["sample_id", "flat_index"],
            how="left",
        )
    else:
        for column in function_columns:
            frame[column] = numpy.nan
    if pfam is not None:
        frame = frame.merge(
            pfam[["sample_id", "flat_index", "arch"]], on=["sample_id", "flat_index"], how="left"
        )
    else:
        frame["arch"] = numpy.nan
    return frame


def _family_named_symbol_counts(annotation, pandas):
    """Per (locus, UniRef50 family), how many of ITS genes carry any gene symbol — the denominator.

    ⭐ **The symbol column’s missing coverage number.** `family_modal_symbols` gives the family’s
    modal symbol and how many DISTINCT symbols it holds, and neither can say how many of its genes
    were named at all. A family of 9 genes of which **2** carry `rfbX` reports `n_sym = 1` — one
    distinct symbol — which the card draws with no marker and a reader reads as nine genes agreeing.
    Measured on kp locus 3992: its 9-gene family names 2 genes, and the locus (97 genes, 7 named) is
    called `mviN` on a vote of 5 to 2.

    ⚠ This is the same distinction `pfam_annotated_member_count` exists for one column to its right
    — *missing evidence is not different evidence* — applied to the column that never had it. **`0`
    is a MEASURED ZERO**: the family’s genes were looked at and none was named.
    """
    named = annotation.dropna(subset=["uniref50", "gene"])
    if not len(named):
        return pandas.DataFrame(columns=["node", "uniref50", "n_named"])
    return (
        named.groupby(["node", "uniref50"], observed=True)
        .size()
        .rename("n_named")
        .reset_index()
    )


def _values_by_node(annotation, column, convert, node_keys):
    """Node → one value per gene that CARRIES `column`, in the shape the database stores it.

    ⛔ **Genes that do not carry it are DROPPED, not carried as `None`.**
    `within_node_propagation.propagation_verdict` counts this list's length, and that length is the
    `one_gene` population the whole propagation claim is about — a `None` per unannotated gene would
    turn every node into a many-gene node and empty the `one_gene` column.

    ⛔⛔ **It asserts its own coverage, because the failure is silent and total.** Node labels look
    numeric and are not: read back through numpy a label can arrive as `2811` or `2811.0` while
    `_locus_rows` writes `str(row["node"])`, and a key that does not match makes every locus
    `no_coverage` — a full column of plausible measured zeros, no error anywhere. So the grouped gene
    count is compared against the column's own non-null count and every key against the node set.
    """
    grouped: dict[str, list] = {}
    grouped_genes = 0
    for node, value in zip(annotation["node"].to_numpy(), annotation[column].to_numpy(), strict=True):
        converted = convert(value)
        if converted is None:
            continue
        grouped.setdefault(str(node), []).append(converted)
        grouped_genes += 1
    carried = int(annotation[column].notna().sum())
    if grouped_genes != carried:
        raise CatalogueFrameError(
            f"`{column}` grouped {grouped_genes:,} genes but {carried:,} rows are non-null — the "
            "parser and the frame disagree about what counts as carried"
        )
    unknown = set(grouped) - set(node_keys)
    if unknown:
        raise CatalogueFrameError(
            f"`{column}` grouped {len(unknown):,} node keys that are in no locus row, e.g. "
            f"{sorted(unknown)[:3]} — the node label was coerced on one side and not the other"
        )
    return grouped


# ── the catalogue ──────────────────────────────────────────────────────────────────────────────
def build_catalogue_frames(
    artifacts: CatalogueArtifacts,
    *,
    annotation_fold: GeneOntologyFold | None = None,
    allow_missing_strand: bool = False,
) -> CatalogueFrames:
    """Assemble every table of one catalogue. The only entry point; everything else is a helper.

    ⭐ **`annotation_fold` arrives as plain dicts and not as a session**, so this stays testable without a
    database while the GO-folded rollups are computed from the SAME reference the request path reads
    — `reference_name_service.gene_ontology_fold`, over the `gene_ontology_term` rows.

    ⛔⛔ **It is NOT nuna's `go_fold`, and the name says so after that cost a real defect.** This
    function also builds `go_fold = slim_of(go_reference())` for `go_verdicts`, a plain dict that
    does no upward closure and strips no namespace root. Written as a parameter called `go_fold` it
    was SHADOWED by that local assignment, so the propagation verdict was computed with the unclosed
    mapping — the `gumC` defect `instruments/gene_ontology_fold.py` exists to fix, under a column
    whose name promises the opposite. Nothing but a type mismatch made it fail rather than compute a
    plausible wrong answer. Two folds answer one question here; they must not share a name.
    Absent a fold, the GO-folded columns are **NULL**.
    """
    import numpy
    import pandas
    from nuna.eval.metrics.oversubscription import node_membership_profile
    from nuna.eval.metrics.syntology import NEAR_SHELLS, near_synteny_agreement
    from nuna.eval.prevalence import categorise
    from nuna.tl.locus_browser.export_payload import (
        MAJOR_SHARE,
        TOP_COG,
        TOP_EC,
        TOP_FAMILIES,
        TOP_GO,
        TOP_KEGG,
        TOP_NEIGHBOURS,
        TOP_PFAM,
        TOP_PRODUCTS,
        TOP_SYMBOLS,
        _floats,
        family_modal_products,
        family_modal_symbols,
        family_pfam,
        neighbour_counts,
        node_order,
        oriented_windows,
        top_counts,
    )
    from nuna.tl.locus_browser.export_payload import (
        OFFSETS as NUNA_OFFSETS,
    )
    from nuna.tl.locus_browser.export_payload import (
        arrangements as build_arrangements,
    )
    from nuna.tl.locus_browser.functional import NS_CODE, cog_summary, go_verdicts, top_go_slim
    from nuna.tl.locus_browser.intergenic import gap_table
    from nuna.tl.locus_browser.vendor_go import go_reference, slim_of

    if tuple(NUNA_OFFSETS) != OFFSETS:
        raise CatalogueFrameError(
            f"nuna's OFFSETS is now {tuple(NUNA_OFFSETS)!r} but this schema stores "
            f"`context_observed_member_counts` as an array in the order {OFFSETS!r} and has a CHECK "
            "over that vocabulary. A reordering would silently put every count under a different "
            "position, so it is refused rather than absorbed."
        )

    source = _load_source_frames(artifacts, allow_missing_strand=allow_missing_strand)
    omitted = source["omitted"]
    audit = _load_audit_frames(artifacts, omitted)
    samples, assign = source["samples"], source["assign"]
    n_genomes = len(samples)
    sample_position = {sample: index for index, sample in enumerate(samples)}

    # ── the locus roster, in the payload's canonical order ─────────────────────────────────────
    profile = node_membership_profile(assign, n_genomes)
    profile["band"] = categorise(profile["genomes_present"].to_numpy(), n_genomes)
    profile = node_order(profile)
    nodes = profile["node"].tolist()
    node_position = {node: index for index, node in enumerate(nodes)}
    n_loci = len(nodes)

    annotation = _annotation_frame(
        assign, source["products"], source["dbxrefs"], source["pfam"], numpy
    )

    # ── per-locus reductions ───────────────────────────────────────────────────────────────────
    named_counts = annotation.assign(_named=annotation["gene"].notna()).groupby("node")["_named"].sum()
    distinct_families = annotation.groupby("node")["uniref50"].nunique()
    labelled = annotation.dropna(subset=["uniref50"])
    labelled_counts = labelled.groupby("node").size()
    family_counts = (
        labelled.groupby(["node", "uniref50"], observed=True).size().rename("n").reset_index()
    )
    family_counts["_total"] = family_counts["node"].map(labelled_counts)
    # ⛔ "Major" = holding ≥ MAJOR_SHARE of the LABELLED members — the audit's own bar. The distinct
    # count above includes single-gene stragglers, so it overstates how many families a locus spans:
    # the card claims with the major count and lists with the distinct one.
    major_families = (
        family_counts[family_counts["n"] >= MAJOR_SHARE * family_counts["_total"]]
        .groupby("node")
        .size()
    )

    top_symbols = top_counts(annotation, "gene", TOP_SYMBOLS)
    top_products = top_counts(annotation, "product", TOP_PRODUCTS)
    top_families = top_counts(annotation, "uniref50", TOP_FAMILIES)
    if len(top_families):
        top_families = (
            top_families.merge(family_modal_products(annotation), on=["node", "uniref50"], how="left")
            .merge(family_pfam(annotation), on=["node", "uniref50"], how="left")
            .merge(family_modal_symbols(annotation), on=["node", "uniref50"], how="left")
            .merge(_family_named_symbol_counts(annotation, pandas), on=["node", "uniref50"], how="left")
        )
        # ⚠ A family no gene of which is named has no row in the join, and its denominator is a
        # MEASURED ZERO — not an absence. Same reason `n_pfam` is written `fillna(0)`.
        top_families["n_named"] = top_families["n_named"].fillna(0).astype("int64")
    top_architectures = top_counts(annotation, "arch", TOP_PFAM)
    # ⛔ How many members carry ANY domain. `top_counts` drops nulls, so without this denominator the
    # page cannot tell a locus whose members DISAGREE about architecture from one where most members
    # were never annotated — and the second read as the first is the "14/99 have different Pfam
    # structures" misreading this column exists to prevent.
    pfam_annotated = annotation.groupby("node")["arch"].count()

    go_fold = slim_of(go_reference())
    go_namespaces = {term: value[1] for term, value in go_reference().items()}
    cog = cog_summary(annotation).set_index("node")
    go_verdict_frame = go_verdicts(annotation, go_fold, go_namespaces).set_index("node")
    top_cog = top_counts(annotation, "cog_id", TOP_COG)
    top_ec = top_counts(annotation, "ec", TOP_EC)
    top_kegg = top_counts(annotation, "kegg_ko", TOP_KEGG)
    top_go = top_go_slim(annotation, go_fold, go_namespaces, TOP_GO)
    ec_annotated = annotation.groupby("node")["ec"].count()
    kegg_annotated = annotation.groupby("node")["kegg_ko"].count()

    # ── the ±5 neighbourhood, built ONCE (10 rows per gene) ────────────────────────────────────
    assign_coords = source["coords"].merge(assign, on=["sample_id", "flat_index"], how="inner")
    # ⚠ 1-based INCLUSIVE, hence the +1: these are the GFF's own coordinates, and getting it wrong is
    # a one-base error in every block width — invisible on screen and wrong in the data.
    gene_lengths = (
        assign_coords["end"].astype("int64") - assign_coords["start"].astype("int64") + 1
    ).rename("nt")
    by_node = pandas.DataFrame(
        {"node": assign_coords["node"].to_numpy(), "nt": gene_lengths.to_numpy()}
    ).groupby("node")["nt"]
    median_length = by_node.median()
    length_iqr = by_node.quantile(0.75) - by_node.quantile(0.25)

    gaps = gap_table(assign_coords, source["noncoding"]) if len(assign_coords) else None
    if gaps is not None and len(gaps):
        gaps = gaps[gaps["a"].isin(node_position) & gaps["b"].isin(node_position)]
    windows = oriented_windows(assign_coords, node_position, radius=NEAR_SHELLS)
    occupants = neighbour_counts(windows, top=TOP_NEIGHBOURS)
    slot_of = {offset: index for index, offset in enumerate(OFFSETS)}
    # ⛔ `top=None` keeps EVERY arrangement — the page's own default. At 4, a quarter of loci were
    # truncated and the genome a reader anchors to often had no row at all.
    arrangement_frame = build_arrangements(
        windows, slot_of, top=None, sample_pos=sample_position
    )
    # A5 groups on locus identity alone (the audit's definition), so the strand bit must not split
    # its groups — hence the rename rather than the packed `neigh_pos`.
    synteny = near_synteny_agreement(
        windows.rename(columns={"neigh_pos": "neigh_id"}), radius=NEAR_SHELLS
    ).set_index("node")["syntenic_A5"]
    gene_memberships = _gene_memberships(
        assign, windows, arrangement_frame, slot_of, numpy, pandas
    )
    del windows

    # ── the audit's own columns ────────────────────────────────────────────────────────────────
    evidence = audit["evidence"]
    if "syntenic_A5" in evidence.columns:
        # The audit computed A5 over this same window; prefer its value verbatim and fill the rest
        # from ours. On the published ecoli catalogue the audit covers 12,020 and the windows 4,889
        # more, so the union is what the page actually shipped.
        synteny = evidence["syntenic_A5"].dropna().combine_first(synteny)

    def evidence_column(name):
        if name not in evidence.columns:
            return pandas.Series(numpy.nan, index=nodes)
        return evidence[name].reindex(nodes)

    concordance = audit["concordance"]

    def concordance_column(name, default):
        if concordance is None or name not in concordance.columns:
            return pandas.Series(default, index=nodes)
        return concordance[name].reindex(nodes)

    modal_symbol = (
        top_symbols[top_symbols["rank"] == 0].set_index("node")["gene"]
        if len(top_symbols)
        else pandas.Series(dtype=object)
    )

    loci = pandas.DataFrame({"node": nodes})
    loci["catalogue_ordinal"] = range(n_loci)
    loci["member_gene_count"] = profile["size"].to_numpy()
    loci["member_genome_count"] = profile["genomes_present"].to_numpy()
    loci["prevalence_band"] = profile["band"].to_numpy()
    loci["bakta_gene_symbol"] = [modal_symbol.get(node) for node in nodes]
    loci["named_member_count"] = named_counts.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["uniref50_family_count"] = distinct_families.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["uniref50_major_family_count"] = major_families.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["uniref50_labelled_member_count"] = labelled_counts.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["pfam_annotated_member_count"] = pfam_annotated.reindex(nodes).fillna(0).astype("int64").to_numpy()
    # ⛔⛔ **SIX COLUMNS THE CLUSTER TABLE CARRIES AND THIS INGEST ONCE DID NOT READ.** Each is
    # documented on `locus` in detail and each was 100 % NULL, which for `uniref50_impurity` is
    # indistinguishable from its DOCUMENTED meaning — *"NULL below 5 labelled members"*. A reader
    # querying `WHERE uniref50_impurity IS NULL` would have concluded that all 17,531 loci have
    # fewer than five labelled members, when in fact **8,576 carry a measured impurity** and the
    # rest were simply never loaded. That is the exact shape of "a column whose meaning and its
    # contents differ", and the payload cannot catch it because the payload does not carry them.
    loci["uniref50_impurity"] = _floats(evidence_column("u50_impurity"), nd=6)
    loci["uniref50_coverage"] = _floats(evidence_column("u50_coverage"), nd=6)
    loci["seqid_coverage"] = _floats(evidence_column("seqid_coverage"), nd=6)
    if "seqid_coverage" not in evidence.columns or not evidence["seqid_coverage"].notna().any():
        # ⛔ NULL now means *the audit did not measure this*, and it says so. Before this ingest read
        # the column at all, NULL meant *we never looked* — the same value carrying a different
        # fact, which is exactly the confusion an omission record exists to prevent. nuna's own
        # comment on the source: *"seqid lens not measured — NaN (not 0, which would read as
        # 'divergent')"*.
        omitted["seqid_to_medoid"] = (
            "the audit ran with --skip-seqid-to-medoid, so member-vs-medoid identity was never "
            "measured — seqid_coverage is NULL for every locus, and that is a fact about the RUN"
        )
    #: ⚠ The medoid is a GENE — one real member, named by its genome and its flat_index — and not a
    #: centroid. ⛔ The four `*_medoid_distance` columns and `embed_within_over_nearest` went on
    #: 2026-09-24 with the card that showed them; these two stay because they NAME the representative
    #: member, which is an ordinary fact about a locus rather than a geometry measured from it.
    loci["medoid_sample_id"] = evidence_column("medoid_sample_id").to_numpy()
    loci["medoid_flat_index"] = [
        None if not numpy.isfinite(float(value)) else int(value)
        for value in evidence_column("medoid_flat_index")
    ]
    loci["syntenic_a5"] = _floats(synteny.reindex(nodes), nd=SYNTENY_DECIMALS)
    loci["collapse_tier"] = evidence_column("tier").to_numpy()
    loci["collapse_bucket"] = evidence_column("bucket").to_numpy()
    loci["resolved_threshold"] = _floats(evidence_column("resolved_threshold"), nd=RESOLVED_DECIMALS)
    loci["pfam_concordance_class"] = concordance_column("class_clan", None).to_numpy()
    loci["pfam_architecture_count"] = [
        None if not numpy.isfinite(float(value)) else int(value)
        for value in concordance_column("n_arch", numpy.nan)
    ]
    loci["cog_annotated_member_count"] = (
        cog["n_cog"].reindex(nodes).fillna(0).astype("int64").to_numpy() if len(cog) else 0
    )
    loci["cog_distinct_id_count"] = (
        cog["cog_n_distinct"].reindex(nodes).fillna(0).astype("int64").to_numpy() if len(cog) else 0
    )
    loci["modal_cog_category"] = (
        cog["cog_cat"].reindex(nodes).to_numpy() if len(cog) else None
    )
    loci["ec_annotated_member_count"] = ec_annotated.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["kegg_annotated_member_count"] = kegg_annotated.reindex(nodes).fillna(0).astype("int64").to_numpy()
    for namespace, code in NS_CODE.items():
        loci[f"go_verdict_{namespace}"] = go_verdict_frame[f"go_{code}"].reindex(nodes).to_numpy()
        loci[f"go_annotated_member_count_{namespace}"] = (
            go_verdict_frame[f"n_go_{code}"].reindex(nodes).fillna(0).astype("int64").to_numpy()
        )

    # ── propagation: does each node's own call agree with itself? ──────────────────────────────
    # ⛔ Read from the SAME frame every rollup above reads, through the SAME two parsers the genome
    # layer writes `gene_functional_annotation` with. Two readers of one parquet reaching one value
    # is the requirement; `tests/test_locus_annotation_rollups.py` is what says whether they do.
    annotation_measured = source["dbxrefs"] is not None
    node_keys = [str(node) for node in nodes]
    # ⚠ `ec` and `go` are comma-joined in the frame and ARRAYS in the database, so they are split to
    # the stored shape; `cog_id` and `kegg_ko` are scalars in both. ⛔ `go` must be split and `ec`
    # need not be — `ec_levels` re-joins its list with commas, so the round trip is a no-op there —
    # but both go through the one splitter, because "this one happens not to need it" is how the
    # other one ends up not getting it either.
    carried_values = {
        "cog": _values_by_node(annotation, "cog_id", scalar_or_none, node_keys),
        "ec": _values_by_node(annotation, "ec", split_set_value, node_keys),
        "kegg": _values_by_node(annotation, "kegg_ko", scalar_or_none, node_keys),
        "go": _values_by_node(annotation, "go", split_set_value, node_keys),
    }
    for suffix, annotation_kind in PROPAGATION_VOCABULARIES.items():
        loci[f"propagation_verdict_{suffix}"] = propagation_verdicts(
            node_keys,
            carried_values[suffix],
            annotation_kind=annotation_kind,
            fold=annotation_fold,
            # ⚠ GO alone is withheld without a fold. The other three use none, so a deployment with
            # no GO reference must not lose them.
            measured=annotation_measured and (suffix != "go" or annotation_fold is not None),
        )
    # ⛔ Genes carrying ANY GO term — NOT the sum of the three namespace counts above, which
    # double-counts every gene annotated in more than one. See `models/locus.py`, where the
    # measured factor is recorded.
    loci["go_annotated_member_count_any"] = (
        annotation.groupby("node")["go"].count().reindex(nodes).fillna(0).astype("int64").to_numpy()
        if annotation_measured
        else None
    )
    claiming = go_claiming_counts(
        node_keys, carried_values["go"], fold=annotation_fold if annotation_measured else None
    )
    for namespace, counts in claiming.items():
        loci[f"go_claiming_member_count_{namespace}"] = counts
    loci["median_gene_length_nt"] = median_length.reindex(nodes).fillna(0).astype("int64").to_numpy()
    loci["gene_length_interquartile_range_nt"] = (
        length_iqr.reindex(nodes).fillna(0).astype("int64").to_numpy()
    )
    loci["context_observed_member_counts"] = _observed_counts(occupants, node_position, n_loci, numpy)
    loci["total_arrangement_count"] = _total_arrangements(arrangement_frame, node_position, n_loci, numpy)
    loci["arrangement_member_gene_count"] = _arrangement_member_total(
        arrangement_frame, node_position, n_loci, numpy
    )
    loci["arrangement_member_genome_count"] = _arrangement_member_genome_total(
        arrangement_frame, node_position, n_loci, numpy, pandas
    )

    annotation_entries = _annotation_entries(
        {
            "sym": top_symbols,
            "prod": top_products,
            "pfam": top_architectures,
            "cog": top_cog,
            "ec": top_ec,
            "kegg": top_kegg,
        },
        top_go,
        node_position,
        pandas,
    )

    return CatalogueFrames(
        loci=loci,
        annotation_entries=annotation_entries,
        uniref_crosstab=_uniref_crosstab(top_families, node_position),
        arrangements=arrangement_frame,
        offset_occupants=occupants,
        gene_memberships=gene_memberships,
        gaps=gaps,
        similarity=_cluster_similarity(artifacts, node_position, n_loci, pandas),
        samples=samples,
        n_genes=source["n_genes"],
        audit_summary=audit["audit_summary"],
        omitted=omitted,
        coverage={
            "loci": n_loci,
            "genomes": n_genomes,
            "genes_assigned": int(len(assign)),
            "genes_with_coordinates": int(len(assign_coords)),
            "arrangements": int(len(arrangement_frame)),
            "offset_occupant_rows": int(len(occupants)),
            "gaps": int(len(gaps)) if gaps is not None else 0,
            "audit_evidence_source": audit["evidence_source"],
            "audit_evidence_loci": int(len(evidence)),
        },
    )


# ── the pieces ─────────────────────────────────────────────────────────────────────────────────
def _observed_counts(occupants, node_position, n_loci, numpy):
    """`ctx.obs` per locus — member genes for which each offset EXISTS, in `OFFSETS` order.

    ⛔ Counted **before** the top-N cut, which is what lets the page say "and 14 others" honestly;
    read against `member_gene_count` it exposes contig-edge truncation, which a short-read assembly
    has a great deal of.
    """
    slot_of = {offset: index for index, offset in enumerate(OFFSETS)}
    observed = numpy.zeros((n_loci, len(OFFSETS)), dtype=numpy.int64)
    if len(occupants):
        first = occupants[occupants["rank"] == 0]
        rows = first["node"].map(node_position).to_numpy(dtype=int)
        columns = first["d"].map(slot_of).to_numpy(dtype=int)
        observed[rows, columns] = first["obs"].to_numpy(dtype=numpy.int64)
    return [row.tolist() for row in observed]


def _total_arrangements(arrangement_frame, node_position, n_loci, numpy):
    """`arr.tot` — arrangements in TOTAL, never moved by any display cap."""
    totals = numpy.zeros(n_loci, dtype=numpy.int64)
    if len(arrangement_frame):
        first = arrangement_frame[arrangement_frame["rank"] == 0]
        totals[first["node"].map(node_position).to_numpy(dtype=int)] = first["total"].to_numpy(
            dtype=numpy.int64
        )
    return totals.tolist()


def _arrangement_member_genome_total(arrangement_frame, node_position, n_loci, numpy, pandas):
    """DISTINCT genomes reaching SOME arrangement — the UNION over all of them, per locus.

    ⛔ **A union, never a sum.** `sum(genomes)` over a locus's arrangements double-counts every
    genome at ρ > 1, which sits in two of them — so a locus where every genome is accounted for
    would report *more* genomes than the locus has, and the API's "is the membership complete"
    question would answer yes at loci where it is no. There is deliberately no uniqueness
    constraint on (locus, genome) for exactly this reason.

    ⚠ `gset` is the arrangement's membership as collection ordinals, written by `build_arrangements`
    in the same pass as the row it rides on, so the two cannot describe different rows.

    Exploded rather than looped: this is O(genes), 486,717 entries at 100 genomes, and the same
    order as the gene-membership pass beside it.
    """
    totals = numpy.zeros(n_loci, dtype=numpy.int64)
    if not len(arrangement_frame) or "gset" not in arrangement_frame.columns:
        return totals.tolist()
    pairs = arrangement_frame[["node", "gset"]].explode("gset").dropna()
    if not len(pairs):
        return totals.tolist()
    counted = pairs.drop_duplicates().groupby("node").size()
    positions = counted.index.map(node_position)
    keep = positions.notna() if hasattr(positions, "notna") else None
    if keep is not None and not keep.all():
        counted = counted[keep]
        positions = positions[keep]
    totals[positions.to_numpy(dtype=int)] = counted.to_numpy(dtype=numpy.int64)
    del pairs
    return totals.tolist()


def _arrangement_member_total(arrangement_frame, node_position, n_loci, numpy):
    """Member genes sitting in SOME arrangement, summed over ALL of them.

    ⛔ Over every arrangement, not the exported few — this is the denominator that lets the API tell
    *"in an arrangement past the display cap"* apart from *"no coordinates for the gene, so no
    window"*. Summed here rather than aggregated per request so the locus view stays at one
    statement per table.
    """
    totals = numpy.zeros(n_loci, dtype=numpy.int64)
    if len(arrangement_frame):
        summed = arrangement_frame.groupby("node")["genes"].sum()
        positions = summed.index.map(node_position)
        keep = positions.notna() if hasattr(positions, "notna") else None
        if keep is not None and not keep.all():
            summed = summed[keep]
            positions = positions[keep]
        totals[positions.to_numpy(dtype=int)] = summed.to_numpy(dtype=numpy.int64)
    return totals.tolist()


def _annotation_entries(top_frames: dict, top_go, node_position, pandas):
    """The six single-column vocabularies plus GO, as one long frame in rank order."""
    kind_by_key = {key: kind for kind, key, _ in ANNOTATION_SOURCES}
    column_by_key = {key: column for _, key, column in ANNOTATION_SOURCES}
    pieces = []
    for key, frame in top_frames.items():
        if frame is None or not len(frame):
            continue
        piece = frame[["node", column_by_key[key], "n", "rank"]].rename(
            columns={column_by_key[key]: "term", "n": "gene_count"}
        )
        piece = piece[piece["node"].isin(node_position)].copy()
        piece["annotation_kind"] = kind_by_key[key]
        piece["namespace"] = None
        pieces.append(piece)
    if top_go is not None and len(top_go):
        piece = top_go[["node", "slim", "n", "rank", "ns"]].rename(
            columns={"slim": "term", "n": "gene_count", "ns": "namespace"}
        )
        piece = piece[piece["node"].isin(node_position)].copy()
        piece["annotation_kind"] = "gene_ontology_slim"
        pieces.append(piece)
    if not pieces:
        return pandas.DataFrame(
            columns=["node", "term", "gene_count", "rank", "annotation_kind", "namespace"]
        )
    return pandas.concat(pieces, ignore_index=True)


def _uniref_crosstab(top_families, node_position):
    """The UniRef50 cross-tab, filtered to loci this catalogue actually has."""
    if top_families is None or not len(top_families):
        return top_families
    return top_families[top_families["node"].isin(node_position)].reset_index(drop=True)


def _gene_memberships(assign, windows, arrangement_frame, slot_of, numpy, pandas):
    """One row per ASSIGNED gene, with the arrangement rank it realises where it has one.

    ⭐ Written in the same pass as `member_genome_ids`, so the array and the rows cannot disagree
    about which genomes carry an arrangement.

    ⛔⛔ **The base is the ASSIGNMENT, not the windows — exactly one row per gene, always.**
    `flank_windows` emits nothing at all for a gene with no neighbour, so a gene alone on its contig
    never reaches the window table. Measured on the published *E. coli* catalogue: **2,352 of
    489,146 genes (0.48 %)**, which is entirely ordinary in draft assemblies whose mean contig holds
    ~14 genes. Building this frame from the windows alone left `gene_locus_membership` — a table
    whose docstring says *"exactly one row per gene"* — short by exactly that many, with nothing
    anywhere reporting it.

    ⚠ Those genes get a row with a **NULL arrangement rank**, which is a different fact from *"this
    genome has no gene at this locus"*. The page says the two in different words, and only a row can
    carry the second.
    """
    base = assign[["sample_id", "flat_index", "node"]].copy()
    base["sample_id"] = base["sample_id"].astype(str)
    if not len(windows) or not len(arrangement_frame):
        return base.assign(arrangement_rank=None)[
            ["node", "sample_id", "flat_index", "arrangement_rank"]
        ]

    columns = [f"s{index}" for index in range(len(slot_of))]
    gene_columns = ["node", "sample_id", "flat_index"]
    codes, unique = pandas.MultiIndex.from_frame(windows[gene_columns]).factorize()
    unique = unique.set_names(gene_columns)
    slots = windows["d"].map(slot_of)
    usable = slots.notna().to_numpy()
    matrix = numpy.full((len(unique), len(slot_of)), -1, dtype=numpy.int64)
    matrix[codes[usable], slots[usable].to_numpy(dtype=int)] = (
        windows["neigh_pos"].to_numpy(dtype=numpy.int64)[usable] * 2
        + windows["same_strand"].to_numpy().astype(numpy.int64)[usable]
    )
    genes = pandas.DataFrame(matrix, columns=columns)
    genes["node"] = unique.get_level_values("node")
    genes["sample_id"] = unique.get_level_values("sample_id")
    genes["flat_index"] = unique.get_level_values("flat_index")
    ranks = arrangement_frame[["node", *columns, "rank"]].rename(columns={"rank": "arrangement_rank"})
    windowed = genes.merge(ranks, on=["node", *columns], how="left")[
        ["node", "sample_id", "flat_index", "arrangement_rank"]
    ]
    windowed["sample_id"] = windowed["sample_id"].astype(str)
    # ⛔ LEFT join from the assignment, so every assigned gene has a row and the ones with no window
    # arrive with a NULL rank rather than being absent.
    return base.merge(
        windowed[["sample_id", "flat_index", "arrangement_rank"]],
        on=["sample_id", "flat_index"],
        how="left",
    )[["node", "sample_id", "flat_index", "arrangement_rank"]]


#: The five per-locus quantities, by the stem they carry in the CSV (each suffixed with the FORM).
#: ⛔ `sep` and `weak_margin` are in the CSV and deliberately NOT read: each is the difference of two
#: columns that ARE, the API subtracts, and a separately-rounded difference drifting from the two
#: rounded numbers printed beside it is the class of quiet disagreement this card was rebuilt to end.
SIMILARITY_COLUMNS = {
    "within_similarity": "within",
    "nearest_similarity": "near",
    "weak_own_similarity": "weak_intra",
    "weak_other_similarity": "weak_inter",
    "own_neighbour_fraction": "nn_own_frac",
}


def _cluster_similarity(artifacts: CatalogueArtifacts, node_position, n_loci: int, pandas) -> dict:
    """Per representation: one row per locus, the ranked nearest loci, and the gene-pair baseline.

    ⛔ **What this replaced on 2026-09-24.** ``_map_geometry`` read a UMAP over every locus's MEDOID
    plus that medoid's 6×6 local geometry; ``_null_baselines`` read the distribution of random MEDOID
    pairs. Both described a construction that reduced a locus to one member. These numbers are
    medians over every gene pair *inside* a locus, and over every gene pair *spanning* two loci.

    ⛔ **The nearest-loci CSV is addressed relative to the similarity CSV**, ``_nearest_sibling``'s own
    rule, so one run's medians can never be paired with another run's neighbour list.
    """
    from bacatlas_backend.ingest.artifact_locator import REPRESENTATIONS

    out: dict[str, dict] = {}
    locus_sets: dict[str, set] = {}
    for representation in REPRESENTATIONS:
        csv_path = artifacts.cluster_similarity(representation)
        metadata_path = artifacts.cluster_similarity_metadata(representation)
        info = dict(
            line.split("=", 1) for line in metadata_path.read_text().splitlines() if "=" in line
        )
        if not info.get("rep"):
            raise CatalogueFrameError(
                f"{metadata_path.name} names no `rep=` — with two representations on one card an "
                "unlabelled artifact would be attached to whichever came first"
            )
        declared = info.get("n_loci")
        if declared and int(declared) != n_loci:
            raise CatalogueFrameError(
                f"{csv_path.name} was computed over {int(declared):,} loci but this catalogue has "
                f"{n_loci:,} — it is a different model. Re-run the similarity for this one."
            )
        # ⛔ Which form the card shows is a DECISION of record (raw — David, 2026-09-24), so a
        # two-form diagnostic run is refused rather than silently halved.
        forms = [form for form in (info.get("forms") or "").split(",") if form]
        if len(forms) != 1:
            raise CatalogueFrameError(
                f"{metadata_path.name} names forms={info.get('forms')!r} — the card shows one. "
                "Ingest the run you mean (`build_cluster_similarity --forms raw`)."
            )
        form = forms[0]

        loci = pandas.read_csv(csv_path, dtype={"node": str})
        unknown = set(loci["node"]) - set(node_position)
        if unknown:
            raise CatalogueFrameError(
                f"{csv_path.name} carries {len(unknown):,} loci this catalogue does not have (e.g. "
                f"{sorted(unknown)[:3]}) — it was built from a different model or assignment."
            )
        missing = [
            f"{stem}_{form}" for stem in SIMILARITY_COLUMNS.values() if f"{stem}_{form}" not in loci
        ]
        if missing:
            raise CatalogueFrameError(
                f"{csv_path.name} has no {', '.join(missing)} — it is not a form={form} artifact"
            )

        locus_sets[representation] = set(loci["node"])
        first = next(iter(locus_sets))
        if locus_sets[representation] != locus_sets[first]:
            difference = len(locus_sets[representation] ^ locus_sets[first])
            raise CatalogueFrameError(
                f"the {representation} and {first} similarity runs disagree about {difference:,} "
                "loci — they are not two views of one catalogue"
            )

        nearest = pandas.read_csv(
            artifacts.cluster_nearest(representation), dtype={"node": str, "neighbour": str}
        )
        nearest = nearest[nearest["form"] == form]

        out[representation] = {
            "info": info,
            "form": form,
            "source": str(csv_path),
            "loci": loci,
            "nearest": nearest,
            "floor": _similarity_floor(artifacts, representation, form, info),
        }
    return out


def _similarity_floor(artifacts: CatalogueArtifacts, representation: str, form: str, info: dict) -> dict:
    """The random GENE-pair baseline: a median from the `.meta`, its quartiles from the audit JSON.

    ⛔ **Not the retired medoid null.** On the published *E. coli* catalogue the gene-pair floor is
    **0.0587** and the medoid one was **0.0651** — close enough to look interchangeable and not be.

    ⚠ The quartiles are optional by design: without the audit JSON the card draws a bare tick instead
    of a p25–p75 box, which is a degradation rather than a failure.
    """
    import json

    floor = {
        "median": float(info[f"floor_{form}"]) if info.get(f"floor_{form}") else None,
        "p25": None,
        "p75": None,
        "p99": None,
        "measurable_locus_count": int(info["n_measurable"]) if info.get("n_measurable") else None,
        "knn_k": int(info["k"]) if info.get("k") else None,
    }
    audit_path = artifacts.cluster_similarity_audit(representation)
    if audit_path.exists():
        forms = json.loads(audit_path.read_text()).get("forms") or {}
        block = (forms.get(form) or {}).get("floor") or {}
        for key, source in (("p25", "floor_p25"), ("p75", "floor_p75"), ("p99", "floor_p99")):
            if source in block:
                floor[key] = float(block[source])
        if floor["median"] is None and "floor" in block:
            floor["median"] = float(block["floor"])
    return floor
