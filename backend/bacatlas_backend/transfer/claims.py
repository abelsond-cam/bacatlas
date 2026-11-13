"""What a node claims — the ONE definition every mechanism shares.

⛔ **This module WRAPS the existing folders; it does not reimplement them.** `gene_ontology_fold`,
`go_rung`, `ec_levels` and `cog_levels` are correct and were hard-won — the GO closure in particular
cost a page that read *"classes differ"* for a term and its own parent. Re-deriving a claim here is
how the two halves of one answer come to disagree, which is the original defect this whole strand
exists to have fixed.

⭐ **A claim is `(vocabulary, facet) → a set of classes`, and the facet is what GO needed.** COG, EC and
KEGG have one axis. GO has three independent ones — what the protein does, what process it is in,
where it is — and pooling them meant two nodes "agreed on GO" by sharing a cellular-component class,
which is not an answer to *what does this protein do*.

⚠ **Correct means INTERSECTS, not equals.** For COG a claim is a single accession, so intersection is
equality; for GO it is a genuine overlap of closed class sets. The two are therefore not the same
measurement, which is why every model downstream is fitted per claim and why their chance baselines
differ by ~400× (COG 0.1 %, cellular component 38.5 %).
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    GENE_ONTOLOGY_FACETS,
    _claims,
)
from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.services.reference_name_service import gene_ontology_fold


@dataclass(frozen=True)
class ClaimKind:
    """One thing a node can claim, and how to get it.

    `rung` is `(facet, depth)` — the key the shared folders return. The facet is `None` for every
    one-axis vocabulary and one of `GENE_ONTOLOGY_FACETS` for GO.
    """

    key: str
    label: str
    annotation_kind: AnnotationKind
    rung: tuple[str | None, int]


#: ⛔ The four claims this work covers. EC (31.9 % of genes) and KEGG (8.2 %) are left out: neither
#: has the reach to carry a transfer argument, and KEGG's licence is unresolved for publication.
CLAIM_KINDS: tuple[ClaimKind, ...] = (
    ClaimKind("cog", "COG orthogroup", AnnotationKind.COG_ORTHOGROUP, (None, 2)),
    ClaimKind("mf", "molecular function", AnnotationKind.GENE_ONTOLOGY_SLIM, (GENE_ONTOLOGY_FACETS[0], 1)),
    ClaimKind("bp", "biological process", AnnotationKind.GENE_ONTOLOGY_SLIM, (GENE_ONTOLOGY_FACETS[1], 1)),
    ClaimKind("cc", "cellular component", AnnotationKind.GENE_ONTOLOGY_SLIM, (GENE_ONTOLOGY_FACETS[2], 1)),
)
BY_KEY = {kind.key: kind for kind in CLAIM_KINDS}


def pangenome_id(session: Session, catalogue: str) -> int:
    """The catalogue's surrogate key."""
    return session.execute(
        select(Pangenome.pangenome_id).where(Pangenome.catalogue_key == catalogue)
    ).scalar_one()


@dataclass(frozen=True)
class Catalogue:
    """One loaded catalogue: every node's size, label, and claim under each vocabulary.

    ⚠ `claims[key]` holds only the nodes that STATE something at that rung. A node absent from it is
    unannotated for that claim, which is the population this whole exercise exists to annotate — so
    absence is the subject, never an error.
    """

    catalogue: str
    pangenome_id: int
    gene_count: dict[int, int]
    node_label: dict[int, str]
    claims: dict[str, dict[int, frozenset[str]]]
    genes_carrying: dict[str, dict[int, int]]

    @property
    def node_total(self) -> int:
        """Every node, annotated or not — the denominator of node coverage."""
        return len(self.gene_count)

    @property
    def gene_total(self) -> int:
        """Every member gene — the denominator of gene coverage."""
        return sum(self.gene_count.values())

    def id_of_label(self) -> dict[str, int]:
        """`node_label → locus_id`. ⛔ Node labels are TEXT and are the durable key across artifacts."""
        return {label: node for node, label in self.node_label.items()}


def _carrying(session: Session, identifier: int) -> dict[str, dict[int, int]]:
    """Per claim, per node, how many member genes state it THEMSELVES — read from stored `locus` columns.

    COG is `cog_annotated_member_count`. GO is `go_claiming_member_count_<facet>`: genes whose terms,
    passed through `go_rung` at ingest (`derive_locus_propagation.go_claiming_counts`), give a
    NON-EMPTY class set in that facet. ⛔ Not `go_annotated_member_count_*`, which counts genes
    carrying a *term* in the namespace — a gene stating only the namespace root is counted there and
    states nothing comparable (`models/locus.py`, beside both columns).

    ⚠ The GO counts are folded at INGEST, so they are meaningful against
    `pangenome.annotation_fold_version`; a NULL (no fold at ingest, or no dbxrefs measured) reads as
    no node carrying that facet.
    """
    facet_columns = {
        "mf": Locus.go_claiming_member_count_molecular_function,
        "bp": Locus.go_claiming_member_count_biological_process,
        "cc": Locus.go_claiming_member_count_cellular_component,
        "cog": Locus.cog_annotated_member_count,
    }
    rows = session.execute(
        select(Locus.locus_id, *facet_columns.values()).where(Locus.pangenome_id == identifier)
    ).all()
    carrying: dict[str, dict[int, int]] = {kind.key: {} for kind in CLAIM_KINDS}
    for locus_id, *counts in rows:
        for key, count in zip(facet_columns, counts, strict=True):
            if count:
                carrying[key][locus_id] = count
    return carrying


def load(session: Session, catalogue: str) -> Catalogue:
    """Load one catalogue's nodes and every claim they state."""
    identifier = pangenome_id(session, catalogue)
    fold = gene_ontology_fold(session)

    gene_count: dict[int, int] = {}
    node_label: dict[int, str] = {}
    for locus_id, label, members in session.execute(
        select(Locus.locus_id, Locus.node_label, Locus.member_gene_count).where(
            Locus.pangenome_id == identifier
        )
    ).all():
        gene_count[locus_id] = members
        node_label[locus_id] = label

    claims: dict[str, dict[int, frozenset[str]]] = {}
    for kind in CLAIM_KINDS:
        folded, _terms, _total = _claims(session, identifier, kind.annotation_kind, fold)
        claims[kind.key] = {
            node: value[kind.rung] for node, value in folded.items() if value.get(kind.rung)
        }

    return Catalogue(
        catalogue=catalogue,
        pangenome_id=identifier,
        gene_count=gene_count,
        node_label=node_label,
        claims=claims,
        genes_carrying=_carrying(session, identifier),
    )


def agrees(one: frozenset[str], other: frozenset[str]) -> bool:
    """⛔ The single definition of *correct* — a non-empty intersection, never equality."""
    return bool(one & other)
