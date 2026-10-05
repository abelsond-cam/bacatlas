"""What a gene's GO terms CLAIM — closed upward through the slim, with the roots removed.

⭐⭐ **This is the correction David's fourth report exposed.** On *Klebsiella* `gumC` the page said
the members' cellular-component classes *differ*: 16 genes say `plasma membrane`, one says
`membrane`. They do not differ. One is the parent of the other, and nothing in the comparison could
see it.

**Why nothing could.** The comparator is `nuna.eval.pfam_annotation.worst_relation` — set algebra
with one escape hatch, the *clan* fold. On Pfam that is sound: clans partition families and families
never nest. GO passes no clan; its fold is `goslim_metagenomics`, and a slim is **a selection, not a
cut through the DAG**. It keeps `membrane`, `plasma membrane` and `outer membrane` as three separate
classes, and all three namespace roots besides. So a term and its own parent landed in different
buckets and were then called `disjoint` — the most condemning verdict available.

⛔ **Full ancestor closure is still rejected, and for the reason `vendor_go` gives**: under it every
pair shares its namespace root and `disjoint` becomes unreachable. What is done here is closure
restricted to the **111 slim classes with the three roots removed**, which is a different thing: two
sibling branches still share nothing. `membrane` and `cytoplasm` remain disjoint, and a test says so.

⚠ **The two measurement paths did not even agree on their INPUT, which is worth saying plainly.**
The neighbour ladder reads `locus_annotation_entry.term_value` — already a slim class, with its
namespace already stored. Within-node propagation reads
`gene_functional_annotation.gene_ontology_terms` — raw accessions, all three namespaces pooled into
one array. Both called `single_rung`, which wraps whatever it is handed, so one compared slim classes
and the other compared raw terms and the page printed both as *"agree"*. `claims` takes either, by
mapping a raw term through the slim and leaving a class that is already slim alone.

⚠ **Empty means *states nothing comparable*, which is neither agreement nor disagreement.** A gene
whose every term is a namespace root — 452 *E. coli* loci state nothing but `cellular_component` —
has an empty claim set. A caller must drop it from the comparison, never score it as agreeing.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

#: ⛔ The three namespace roots. `cellular_component` is an ancestor of every cellular-component
#: term, so two genes sharing it share nothing — it says *"this is somewhere in a cell"*. They are
#: removed from a claim set BEFORE the closure, which is what keeps `disjoint` reachable. Mirrors
#: `nuna.tl.locus_browser.vendor_go.NAMESPACE_ROOTS`; the ingest strips them from the tables too, so
#: this is the second of two guards rather than the only one.
NAMESPACE_ROOTS = frozenset({"GO:0003674", "GO:0008150", "GO:0005575"})


@dataclass(frozen=True)
class GeneOntologyFold:
    """The GO reference, as the measurement needs it. Plain dicts, so it is testable without a database.

    ⚠ **An EMPTY fold is a valid state and must behave like today.** A deployment that has not loaded
    the reference gets `{}` for all three maps; `claims` then returns the terms it was given, minus
    the roots and unclosed — which is exactly the old behaviour, so a missing reference degrades the
    measurement rather than emptying it.
    """

    #: Raw GO term → its nearest slim classes. A term already in the slim maps to itself; a term that
    #: reaches none maps to itself too (`vendor_go.slim_of`'s rule), because folding it away would
    #: make an annotated gene read as unannotated.
    slim_of: dict[str, frozenset[str]] = field(default_factory=dict)
    #: Slim class → every slim class above it, transitively. No root appears on either side.
    ancestors: dict[str, frozenset[str]] = field(default_factory=dict)
    #: Slim class → 0 molecular_function · 1 biological_process · 2 cellular_component.
    namespace_of: dict[str, int] = field(default_factory=dict)

    def claims(self, terms: Iterable[str]) -> frozenset[str]:
        """Raw or slim GO terms → the closed, root-free set of slim classes they claim.

        The order is load-bearing and is the one thing to get right here:

        1. **fold** each term onto its slim classes (identity for a class already slim);
        2. **drop the roots** — before the closure, so a root can neither be claimed nor dragged
           back in as an ancestor of something else;
        3. **close upward** within the slim, which is what lets `plasma membrane` meet `membrane`.

        ⛔ Reversing 2 and 3 would put the roots straight back: every cellular-component class has
        `GO:0005575` above it, so closing first and stripping after is not the same operation — it
        only looks like it.
        """
        folded: set[str] = set()
        for term in terms:
            if not term:
                continue
            folded.update(self.slim_of.get(term, frozenset({term})))
        folded -= NAMESPACE_ROOTS
        closed = set(folded)
        for go_class in folded:
            closed |= self.ancestors.get(go_class, frozenset())
        return frozenset(closed)

    def claims_by_namespace(self, terms: Iterable[str]) -> dict[int, frozenset[str]]:
        """The same, split into the three namespaces — only the namespaces this gene states.

        ⭐ **The split is why GO's chance baseline was 11-17 % against COG's 0.2 %.** Pooled, two
        nodes "agree on GO" if they share *any* class in *any* namespace, so a shared
        `plasma membrane` made them agree about what the protein does. A molecular function and a
        cellular component are not rival answers to one question; they are answers to two.

        ⚠ A namespace the gene says nothing in is **absent**, not empty — absent means *not
        comparable here*, and a caller must drop such a pair rather than score it either way.
        """
        split: dict[int, set[str]] = {}
        for go_class in self.claims(terms):
            namespace = self.namespace_of.get(go_class)
            if namespace is not None:
                split.setdefault(namespace, set()).add(go_class)
        return {namespace: frozenset(classes) for namespace, classes in split.items()}
