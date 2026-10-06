"""Names for the two vocabularies that reach the page as bare codes.

`term_name` is populated for 100 % of the catalogue's COG and GO entries and **0 %** of its EC and
KEGG ones, because Bakta does not emit a name for either. So a locus card showed `2.7.10.-` and
`K01991` and left the reader to look them up. These two tables close that, joined at serve time
rather than written into `locus_annotation_entry`, so refreshing a reference is
`ingest --stage reference` and never a catalogue re-ingest.

⛔ **Cached per PROCESS, like `_CALIBRATIONS`** — a global vocabulary does not change between
requests, and re-reading 28,516 KEGG rows per request would be visible in the route budget.
`clear_reference_name_cache()` exists for the tests that load a reference mid-process.

⚠ **An absent table is the normal state, not a failure.** Every deployment that has not accepted the
KEGG licence position simply does not load it; both lookups then return `None` and the card falls
back to a bare linked accession, exactly as the published static pages do permanently.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold
from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.models.reference_vocabulary import (
    EnzymeClass,
    GeneOntologySlimAncestor,
    GeneOntologyTerm,
    KeggOrthology,
)

#: How a multi-code EC value's names are joined. ⚠ The same separator the card already uses between
#: distinct entries, because an EC *value* is itself a set and the two levels read alike.
EC_JOIN = " · "


@dataclass(frozen=True)
class ReferenceNames:
    """The two lookups, loaded once. Empty dicts are a valid, expected state."""

    enzyme: dict[str, str] = field(default_factory=dict)
    kegg: dict[str, str] = field(default_factory=dict)

    def enzyme_name(self, value: str) -> str | None:
        """`2.7.10.-` → *Protein-tyrosine kinases*; `1.6.5.9,7.1.1.-` → both names, joined.

        ⛔ **An EC value is a SET, not a string.** 1,032 of the catalogue's 9,398 EC entries are
        comma-joined lists and some carry whitespace after the comma, so this splits and strips
        before every lookup. Doing it in SQL with `split_part` returns a code neither side holds.

        ⚠ **All or nothing.** If any code in the list is unnamed the whole value is left unnamed,
        because a partial join (*"Nucleotidyltransferases"* for `2.7.7.-,3.1.-.-`) reads as the
        complete answer and is the half that happens to be first.
        """
        codes = [code.strip() for code in value.split(",") if code.strip()]
        if not codes:
            return None
        names = [self.enzyme.get(code) for code in codes]
        if any(name is None for name in names):
            return None
        # A value may repeat a code; keep first-seen order rather than sorting, which would state an
        # ordering the annotation does not have.
        seen: dict[str, None] = {}
        for name in names:
            assert name is not None
            seen.setdefault(name, None)
        return EC_JOIN.join(seen)

    def kegg_name(self, value: str) -> str | None:
        """`K01991` → *polysaccharide biosynthesis/export protein*.

        ⚠ The **definition**, not the symbol. `K01991`'s symbol is `wza, gfcE`, and the page already
        shows this node's own gene symbols a few hundred pixels away — repeating them says nothing,
        whereas the definition is the thing a bare accession withholds.
        """
        return self.kegg.get(value.strip()) or None

    def name_for(self, annotation_kind: AnnotationKind | str, value: str | None) -> str | None:
        """The name for one (vocabulary, value), or `None` where this is not a vocabulary we name."""
        if value is None:
            return None
        kind = annotation_kind.value if isinstance(annotation_kind, AnnotationKind) else annotation_kind
        if kind == AnnotationKind.EC_NUMBER.value:
            return self.enzyme_name(value)
        if kind == AnnotationKind.KEGG_ORTHOLOGY.value:
            return self.kegg_name(value)
        return None


_NAMES: ReferenceNames | None = None


def reference_names(session: Session) -> ReferenceNames:
    """The process-cached lookups, loading them on first use.

    ⚠ Two statements on a cold process and **zero** thereafter, which is why the Function tab's warm
    route budget is unmoved. A deployment with neither table loaded pays the two statements once and
    caches two empty dicts — it does not retry per request.
    """
    global _NAMES
    if _NAMES is None:
        enzyme = dict(session.execute(select(EnzymeClass.ec_code, EnzymeClass.name)).all())
        kegg = dict(session.execute(select(KeggOrthology.ko_id, KeggOrthology.definition)).all())
        _NAMES = ReferenceNames(enzyme=enzyme, kegg=kegg)
    return _NAMES


def clear_reference_name_cache() -> None:
    """Drop the cache — for tests that load a reference after the process has already read it."""
    global _NAMES
    _NAMES = None


_FOLD: GeneOntologyFold | None = None


def gene_ontology_fold(session: Session) -> GeneOntologyFold:
    """The GO reference the measurement compares through, cached per process.

    ⚠ **38,092 terms and 99 ancestor edges, read once.** Two statements on a cold process and zero
    thereafter, exactly like `reference_names` and `_CALIBRATIONS` — the Function tab's warm route
    budget is what makes that a requirement rather than an optimisation.

    ⛔ **An empty fold is a valid, expected state**, not a failure: a deployment that has not run
    `ingest --stage reference` gets one, and `GeneOntologyFold.claims` then returns what it was given
    minus the roots. That is the behaviour this code had before the reference existed, so a missing
    table degrades the measurement rather than emptying it.
    """
    global _FOLD
    if _FOLD is None:
        slim_of: dict[str, frozenset[str]] = {}
        namespace_of: dict[str, int] = {}
        rows = session.execute(
            select(
                GeneOntologyTerm.go_id,
                GeneOntologyTerm.namespace_index,
                GeneOntologyTerm.slim_go_ids,
                GeneOntologyTerm.name,
            )
        ).all()
        names: dict[str, str] = {}
        for go_id, namespace_index, slim_go_ids, go_name in rows:
            names[go_id] = go_name
            # ⚠ `vendor_go.slim_of`'s rule, carried over: a term reaching NO slim class maps to
            # itself. Folding it away would make an annotated gene read as unannotated, which is a
            # different and false finding from "annotated with nothing comparable".
            slim_of[go_id] = frozenset(slim_go_ids) if slim_go_ids else frozenset({go_id})
            namespace_of[go_id] = namespace_index
        ancestors: dict[str, set[str]] = {}
        for go_id, ancestor in session.execute(
            select(GeneOntologySlimAncestor.go_id, GeneOntologySlimAncestor.ancestor_go_id)
        ).all():
            ancestors.setdefault(go_id, set()).add(ancestor)
        _FOLD = GeneOntologyFold(
            slim_of=slim_of,
            ancestors={k: frozenset(v) for k, v in ancestors.items()},
            namespace_of=namespace_of,
            name_of=names,
        )
    return _FOLD


def clear_gene_ontology_fold_cache() -> None:
    """Drop the cached fold — for tests that load the reference after the process has read it."""
    global _FOLD
    _FOLD = None
