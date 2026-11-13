"""The propagation and GO-claim columns, from distinct values, against the per-gene route they replaced.

`_values_by_node` used to hold one Python value per carrying gene, for all four vocabularies at once —
~250 B a gene retained and ~340 B at peak, ~94 / 128 GB at the 71,435-genome cohort's 374 M genes —
and `propagation_verdicts` / `go_claiming_counts` then folded every one of those genes again. It now
holds each node's DISTINCT values with how many genes carry each. The decision
(`within_node_propagation.propagation_verdict`) and the claim (`go_rung`) are unchanged and still the
only doors; nothing a locus row stores may move. So the old route is kept HERE, verbatim, as the
oracle, and the two are compared on frames built to reach every state of the verdict.
"""

from __future__ import annotations

import random

import numpy as np
import pandas as pd
import pytest

from bacatlas_backend.ingest import catalogue_frames as cf
from bacatlas_backend.ingest.catalogue_frames import CatalogueFrameError, _propagation_columns
from bacatlas_backend.ingest.derive_locus_propagation import PROPAGATION_VOCABULARIES
from bacatlas_backend.ingest.ingest_genome_gene_table import scalar_or_none, split_set_value
from bacatlas_backend.instruments.annotation_transfer import go_rung
from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold
from bacatlas_backend.instruments.within_node_propagation import propagation_verdict
from bacatlas_backend.models.enumerations import GENE_ONTOLOGY_NAMESPACE_NAMES, PropagationVerdict


# ───────────────────────────────────────────────────────── the oracle: the per-gene originals, verbatim
def values_by_node_per_gene(annotation, column, convert, node_keys):
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


def propagation_verdicts_per_gene(nodes, values_by_node, *, annotation_kind, fold=None, measured=True):
    if not measured:
        return [None] * len(nodes)
    return [propagation_verdict(annotation_kind, values_by_node.get(node, []), fold) for node in nodes]


def go_claiming_counts_per_gene(nodes, terms_by_node, *, fold):
    if fold is None:
        return {namespace: [None] * len(nodes) for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES}
    counts = {namespace: [] for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES}
    for node in nodes:
        per_namespace = dict.fromkeys(GENE_ONTOLOGY_NAMESPACE_NAMES, 0)
        for terms in terms_by_node.get(node, []):
            for (facet, _depth), classes in go_rung([term for term in terms if term], fold).items():
                if classes:
                    per_namespace[facet] += 1
        for namespace, count in per_namespace.items():
            counts[namespace].append(count)
    return counts


def propagation_columns_per_gene(annotation, node_keys, annotation_fold, annotation_measured):
    """`build_catalogue_frames`' block as it stood: all four vocabularies grouped first, then judged."""
    carried_values = {
        "cog": values_by_node_per_gene(annotation, "cog_id", scalar_or_none, node_keys),
        "ec": values_by_node_per_gene(annotation, "ec", split_set_value, node_keys),
        "kegg": values_by_node_per_gene(annotation, "kegg_ko", scalar_or_none, node_keys),
        "go": values_by_node_per_gene(annotation, "go", split_set_value, node_keys),
    }
    verdicts = {}
    for suffix, annotation_kind in PROPAGATION_VOCABULARIES.items():
        verdicts[f"propagation_verdict_{suffix}"] = propagation_verdicts_per_gene(
            node_keys,
            carried_values[suffix],
            annotation_kind=annotation_kind,
            fold=annotation_fold,
            measured=annotation_measured and (suffix != "go" or annotation_fold is not None),
        )
    claiming = go_claiming_counts_per_gene(
        node_keys, carried_values["go"], fold=annotation_fold if annotation_measured else None
    )
    return verdicts, claiming


# ─────────────────────────────────────────────────────────────────────────────────── the reference
MEMBRANE, PLASMA, OUTER = "GO:0016020", "GO:0005886", "GO:0019867"
CYTOPLASM, INTRACELLULAR, KINASE, PROCESS = "GO:0005737", "GO:0005622", "GO:0016301", "GO:0008152"
CC_ROOT, MF_ROOT = "GO:0005575", "GO:0003674"
FOLD = GeneOntologyFold(
    slim_of={
        MEMBRANE: frozenset({MEMBRANE}), PLASMA: frozenset({PLASMA}), OUTER: frozenset({OUTER}),
        CYTOPLASM: frozenset({CYTOPLASM}), INTRACELLULAR: frozenset({INTRACELLULAR}),
        KINASE: frozenset({KINASE}), PROCESS: frozenset({PROCESS}),
        CC_ROOT: frozenset({CC_ROOT}), MF_ROOT: frozenset({MF_ROOT}),
        "GO:0031226": frozenset({PLASMA}),  # a RAW term folding onto a class
    },
    ancestors={PLASMA: frozenset({MEMBRANE}), OUTER: frozenset({MEMBRANE}), CYTOPLASM: frozenset({INTRACELLULAR})},
    namespace_of={
        MEMBRANE: 2, PLASMA: 2, OUTER: 2, CYTOPLASM: 2, INTRACELLULAR: 2, KINASE: 0, PROCESS: 1,
        CC_ROOT: 2, MF_ROOT: 0, "GO:0031226": 2,
    },
)  # fmt: skip
# GO:0099999 is in no table — obsolete, or newer than the reference: no namespace, so no claim.

#: One locus per state, named for it: (cog_id, ec, kegg_ko, go) per gene.
EDGES: dict[str, list[tuple]] = {
    "no_coverage": [(None, None, None, None)] * 3,
    "one_gene": [("COG1", "1.1.1.1", "K1", PLASMA), (None, None, None, None)],
    "unanimous_many": [("COG1", "1.1.1.1", "K1", f"{PLASMA},{KINASE}")] * 5,
    "discordant": [
        ("COG1", "1.1.1.1", "K1", KINASE),
        ("COG2", "1.1.1.2", "K2", PROCESS),
        ("COG2", "1.1.1.2", "K2", PROCESS),
    ],
    # EC agrees only at a shallow level: L4 is not stated by both, so L3 is compared — and agrees
    "shallow_ec": [(None, "1.1.1.-", None, None), (None, "1.1.1.7", None, None)],
    "ec_set_values": [(None, "1.6.5.9,7.1.1.-", None, None), (None, "7.1.1.2", None, None)],
    "ec_not_comparable": [(None, "-.-.-.-", None, None), (None, "1.1.1.1", None, None)],
    # GO agrees only through the closure: plasma membrane under membrane (the gumC case)
    "shallow_go": [(None, None, None, PLASMA), (None, None, None, MEMBRANE), (None, None, None, "GO:0031226")],
    "go_root_only": [(None, None, None, CC_ROOT), (None, None, None, f"{CC_ROOT},{MF_ROOT}")],
    "go_root_and_term": [(None, None, None, f"{CC_ROOT},{OUTER}"), (None, None, None, OUTER)],
    "go_two_namespaces": [(None, None, None, f"{KINASE},{PLASMA}"), (None, None, None, f"{KINASE},{CYTOPLASM}")],
    "go_unknown_term": [(None, None, None, "GO:0099999"), (None, None, None, "GO:0099999")],
    "go_duplicates": [(None, None, None, f"{KINASE},{KINASE}"), (None, None, None, KINASE)],
    "go_one_annotated": [(None, None, None, PROCESS), (None, None, None, None), (None, None, None, None)],
    "go_siblings": [(None, None, None, MEMBRANE), (None, None, None, CYTOPLASM)],
    # numeric-looking labels: the coercion the coverage assertion exists for
    "2811": [("COG9", "2.7.1.1", "K9", KINASE)] * 2,
    "07": [("COG9", None, None, None)],
}


def _edge_frame() -> tuple[pd.DataFrame, list[str]]:
    rows = [(node, *gene) for node, genes in EDGES.items() for gene in genes]
    random.Random(0).shuffle(rows)  # a locus's genes are not contiguous in the real frame
    frame = pd.DataFrame(rows, columns=["node", "cog_id", "ec", "kegg_ko", "go"])
    keys = list(EDGES) + ["never_seen"]  # a locus row with no gene in this frame at all
    return frame, keys


def _random_frame(seed: int, n_genes: int = 500, n_nodes: int = 30) -> tuple[pd.DataFrame, list[str]]:
    rng = random.Random(seed)
    go_vocabulary = [MEMBRANE, PLASMA, OUTER, CYTOPLASM, INTRACELLULAR, KINASE, PROCESS, CC_ROOT, MF_ROOT,
                     "GO:0031226", "GO:0099999"]  # fmt: skip
    ec_vocabulary = ["1.1.1.1", "1.1.1.2", "1.1.1.-", "1.1.-.-", "2.7.1.1", "-.-.-.-", "7.1.1.2"]

    def maybe(make, p=0.35):
        return None if rng.random() < p else make()

    rows = []
    for _ in range(n_genes):
        rows.append((
            f"{rng.randrange(n_nodes)}",
            maybe(lambda: rng.choice(["COG1", "COG2", "COG3"])),
            maybe(lambda: ",".join(rng.sample(ec_vocabulary, rng.randint(1, 2)))),
            maybe(lambda: rng.choice(["K1", "K2"]), 0.6),
            maybe(lambda: ",".join(rng.choice(go_vocabulary) for _ in range(rng.randint(1, 3)))),
        ))  # fmt: skip
    frame = pd.DataFrame(rows, columns=["node", "cog_id", "ec", "kegg_ko", "go"])
    return frame, [str(n) for n in range(n_nodes)]


def _same(got, want, n_nodes):
    """Coverage first — the same columns, every list as long as the locus roster — then exact values."""
    for g, w in zip(got, want, strict=True):
        assert list(g) == list(w) and len(w) > 0
        for column in w:
            assert len(g[column]) == len(w[column]) == n_nodes, column
            for i, (a, b) in enumerate(zip(g[column], w[column], strict=True)):
                assert type(a) is type(b) and a == b, (column, i, a, b)


FRAMES = {"edges": _edge_frame, **{f"random{seed}": (lambda seed=seed: _random_frame(seed)) for seed in range(25)}}


@pytest.mark.parametrize("frame", FRAMES)
@pytest.mark.parametrize("fold", ["closed", "empty", "none"])
@pytest.mark.parametrize("measured", [True, False])
def test_propagation_columns_are_bit_exact_against_the_per_gene_route(frame, fold, measured):
    annotation, node_keys = FRAMES[frame]()
    folds = {"closed": FOLD, "empty": GeneOntologyFold(), "none": None}
    want = propagation_columns_per_gene(annotation, node_keys, folds[fold], measured)
    got = _propagation_columns(annotation, node_keys, folds[fold], measured)
    _same(got, want, len(node_keys))


def test_the_edge_frame_reaches_every_state_it_was_built_for():
    """Coverage of the oracle itself — a fixture that never produced a state could not test it."""
    annotation, node_keys = _edge_frame()
    verdicts, claiming = propagation_columns_per_gene(annotation, node_keys, FOLD, True)
    seen = {v for column in verdicts.values() for v in column}
    assert seen == set(PropagationVerdict), seen
    by_node = dict(zip(node_keys, verdicts["propagation_verdict_go"], strict=True))
    assert by_node["shallow_go"] is PropagationVerdict.UNANIMOUS, "the closure must meet plasma/membrane"
    assert by_node["go_root_only"] is PropagationVerdict.NOT_COMPARABLE
    assert by_node["go_one_annotated"] is PropagationVerdict.ONE_GENE
    ec_by_node = dict(zip(node_keys, verdicts["propagation_verdict_ec"], strict=True))
    assert ec_by_node["shallow_ec"] is PropagationVerdict.UNANIMOUS
    assert any(claiming["cellular_component"]) and any(claiming["molecular_function"])


def test_integer_node_labels_stringify_exactly_as_before():
    annotation, node_keys = _random_frame(3)
    annotation["node"] = annotation["node"].astype("int64")
    _same(
        _propagation_columns(annotation, node_keys, FOLD, True),
        propagation_columns_per_gene(annotation, node_keys, FOLD, True),
        len(node_keys),
    )


def test_a_dbxref_frame_that_was_never_measured_is_all_nan():
    """`_annotation_frame` fills every function column with NaN when there is no dbxref parquet."""
    annotation, node_keys = _random_frame(5)
    for column in ("cog_id", "ec", "kegg_ko", "go"):
        annotation[column] = np.nan
    _same(
        _propagation_columns(annotation, node_keys, FOLD, False),
        propagation_columns_per_gene(annotation, node_keys, FOLD, False),
        len(node_keys),
    )


@pytest.mark.parametrize(
    "spoil, message",
    [
        (lambda f: f.assign(ec=f["ec"].where(f.index != 4, "")), "parser and the frame disagree"),
        (lambda f: f.assign(go=f["go"].where(f.index != 6, "  ")), "parser and the frame disagree"),
        (lambda f: f.assign(node=f["node"].where(f.index != 2, "999.0")), "node label was coerced"),
    ],
)
def test_the_coverage_assertions_still_refuse_with_the_same_words(spoil, message):
    """⛔ The silent-and-total failures `_values_by_node` exists to refuse: same exception, same text."""
    annotation, node_keys = _random_frame(9)
    annotation = spoil(annotation)
    with pytest.raises(CatalogueFrameError) as old:
        propagation_columns_per_gene(annotation, node_keys, FOLD, True)
    with pytest.raises(CatalogueFrameError) as new:
        _propagation_columns(annotation, node_keys, FOLD, True)
    assert message in str(old.value) and str(new.value) == str(old.value)


def test_values_by_node_holds_distinct_values_whose_counts_sum_to_the_carriers():
    annotation, node_keys = _random_frame(1)
    grouped = cf._values_by_node(annotation, "go", split_set_value, node_keys)
    per_gene = values_by_node_per_gene(annotation, "go", split_set_value, node_keys)
    assert set(grouped) == set(per_gene)
    for node, entries in grouped.items():
        assert sum(genes for _value, genes in entries) == len(per_gene[node])
        assert sorted(map(tuple, (v for v, _ in entries))) == sorted(set(map(tuple, per_gene[node])))
