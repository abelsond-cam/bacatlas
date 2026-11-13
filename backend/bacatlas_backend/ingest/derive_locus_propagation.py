"""Per locus, whether its own annotated genes agree — the rollup, computed once at ingest.

⭐ **The decision is NOT here.** It is `instruments/within_node_propagation.propagation_verdict`.
This module is only the loop over a catalogue's nodes
and the two awkward details that loop has to get right: the shape each vocabulary's value is held
in, and what to write when a vocabulary was never measured at all.

The input is nuna's annotation frame. (While the per-gene tables existed, the stored verdicts were
checked locus by locus against a recount over them; that oracle went with the tables in D8.)

⚠ **A missing vocabulary must reach the database as NULL and never as `no_coverage`.** When a
catalogue has no `{BS}_dbxref.parquet`, `catalogue_frames._annotation_frame` fills every function
column with NaN, every node then carries no value, and a verdict derived without care says
`no_coverage` for all 17,531 loci — *we looked and found nothing* — when the truth is *we never
looked*. That is the exact failure recorded at `ingest/catalogue_frames.py` for the six columns that
were 100 % NULL while their documented meaning was "not measured", and `measured=False` is how this
module refuses to express it.
"""

from __future__ import annotations

from bacatlas_backend.instruments.annotation_transfer import go_rung
from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold
from bacatlas_backend.instruments.within_node_propagation import propagation_verdict
from bacatlas_backend.models.enumerations import (
    GENE_ONTOLOGY_NAMESPACE_NAMES,
    AnnotationKind,
    PropagationVerdict,
)

#: ⛔ The four vocabularies a propagation verdict exists for, and the column suffix each is stored
#: under. ⚠ The suffixes are the house abbreviations the neighbouring `*_annotated_member_count`
#: columns already use, NOT `AnnotationKind`'s own names — `cog_orthogroup_annotated_member_count`
#: is not a column anywhere and inventing a second spelling here is how one half of a response comes
#: to read a different column from the other.
PROPAGATION_VOCABULARIES: dict[str, AnnotationKind] = {
    "cog": AnnotationKind.COG_ORTHOGROUP,
    "ec": AnnotationKind.EC_NUMBER,
    "kegg": AnnotationKind.KEGG_ORTHOLOGY,
    "go": AnnotationKind.GENE_ONTOLOGY_SLIM,
}

#: ⚠ Only the GO verdict and the GO claim counts are folded. `cog`, `ec` and `kegg` reach
#: `propagation_verdict` unfolded and are unaffected by whether a reference is loaded, so they are
#: written whatever the fold is — and must be, or a deployment with no GO reference would lose three
#: measurements that do not depend on it.
FOLDED_VOCABULARIES = frozenset({"go"})


def _decision_values(entries: list[tuple]) -> list:
    """A node's distinct `(value, genes)` entries → the values list `propagation_verdict` is asked about.

    ⭐ **Exact, not a sample.** `propagation_verdict` reads its list's LENGTH only as none / one / more,
    and its CONTENT only through an all-genes test and an intersection over the folded values, which no
    repeat can change. So the distinct values give the verdict the full per-gene list gives, provided
    the length still says how many genes there were: a node whose many genes all carry ONE value is
    asked about that value twice, and stays a many-gene node rather than becoming `one_gene`.
    """
    values = [value for value, _genes in entries]
    if len(values) == 1 and entries[0][1] >= 2:
        return values * 2
    return values


def propagation_verdicts(
    nodes: list[str],
    values_by_node: dict[str, list[tuple]],
    *,
    annotation_kind: AnnotationKind,
    fold: GeneOntologyFold | None = None,
    measured: bool = True,
) -> list[PropagationVerdict | None]:
    """One verdict per node, in `nodes` order. `None` for every node when `measured` is False.

    `values_by_node` holds, per node, the DISTINCT values its carrying genes hold with how many genes
    hold each — `catalogue_frames._values_by_node`'s shape — so a node absent from it has no annotated
    gene and gets `no_coverage`, which is a measured statement rather than an absence.
    """
    if not measured:
        return [None] * len(nodes)
    return [
        propagation_verdict(annotation_kind, _decision_values(values_by_node.get(node, [])), fold)
        for node in nodes
    ]


def go_claiming_counts(
    nodes: list[str],
    terms_by_node: dict[str, list[tuple[list[str], int]]],
    *,
    fold: GeneOntologyFold | None,
) -> dict[str, list[int | None]]:
    """Per namespace, how many of each node's genes state a NON-EMPTY folded claim there.

    ⛔ **A different population from `go_annotated_member_count_*`, and always the smaller one.** That
    column counts a gene carrying a term in the namespace; this counts a gene whose terms survive the
    fold, so a gene asserting only the namespace root is in the first and not the second. The
    measured gaps are on `models/locus.py`, beside the columns.

    ⚠ **No fold means no count — NULL, not zero.** Without a reference `go_rung` pools all three
    facets into the first (`annotation_transfer.py:308-315`), so counting there would report every
    gene's whole claim under *molecular function* and nothing under the other two: three populated,
    plausible, wrong columns. `None` says the fold was absent, which is the truth.

    `terms_by_node` is `catalogue_frames._values_by_node`'s shape: per node, each distinct term list
    with how many genes carry it. A claim is a function of the term list alone, so each distinct list
    is folded once and its facets are counted once per gene carrying it.
    """
    if fold is None:
        return {namespace: [None] * len(nodes) for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES}
    counts: dict[str, list[int | None]] = {
        namespace: [] for namespace in GENE_ONTOLOGY_NAMESPACE_NAMES
    }
    claimed: dict[tuple[str, ...], tuple[str, ...]] = {}
    for node in nodes:
        per_namespace = dict.fromkeys(GENE_ONTOLOGY_NAMESPACE_NAMES, 0)
        for terms, genes in terms_by_node.get(node, []):
            key = tuple(terms)
            facets = claimed.get(key)
            if facets is None:
                # ⚠ `go_rung` is the ONE door every GO claim passes through — the ladder, within-node
                # propagation and this count all read the same closed, root-free sets. Its key is
                # `(facet, depth)`, never a namespace index, whatever its return annotation says.
                facets = claimed[key] = tuple(
                    facet
                    for (facet, _depth), classes in go_rung([term for term in terms if term], fold).items()
                    if classes
                )
            for facet in facets:
                per_namespace[facet] += genes
        for namespace, count in per_namespace.items():
            counts[namespace].append(count)
    return counts
