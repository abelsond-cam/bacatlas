"""Public reference tables the catalogue JOINS against — owned by nobody here, versioned by us.

⛔ **Not a pangenome's data, and not a generic key/value table.** A Pfam family has six specific
fields the page reads individually — the short name it is chipped with, the description a hover
shows, and the InterPro entry a reader following a domain actually wants — so a generic
`(vocabulary, key, value)` shape would flatten exactly the structure the page needs and turn every
render into six lookups. Other vocabularies get their own tables here when they earn one.

⛔ **KEGG is here now, and the licence is still UNRESOLVED.** It was absent on the reasoning that
its terms permit linking and not redistribution — which remains true of any *public* deployment.
David decided on 2026-10-05 that this internal prototype names KOs anyway and carries the question
on the page: *"Don't worry about KEGG licence! Just flag it on site. Will look at it later."* So
`KeggOrthology` exists, the service joins it, and **nothing may ship to BacAtlas.org until the
academic service provider licence is confirmed**. The published static pages still link out and name
nothing; see `nuna.tl.locus_browser.vendor_kegg`.

⭐ **EC is the opposite case, and we had it backwards.** ExPASy ENZYME is **CC BY 4.0** — permissive
including derivatives, on one condition, attribution — not the CC BY-ND an earlier note in the Vue
app claimed. The cost of that error was a page that said *"transferase · 4th level not stated"*
where it could have said *"Protein-tyrosine kinases"*.
"""

from __future__ import annotations

from sqlalchemy import ForeignKey, Index, SmallInteger, String, Text
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base


class PfamFamily(Base):
    """One Pfam-A family: what it is called, what it is, and where it sits in InterPro and its clan.

    Vendored from `nuna.tl.locus_browser.vendor_reference`, which is the **single reader** of
    `pfam_names.tsv.gz` — three call sites wanted it and three private gzip parsers is exactly the
    drift that lets two of them disagree about which column the clan is in.

    ⚠ **A clanless family has an empty `clan_accession`, and that is not an identity.** Only ~46 %
    of families are in a clan at all, so mapping the clanless to a shared blank would make any two
    of them look like the same superfamily — the common case, not the corner. `clan_of` in nuna maps
    a clanless family to *itself* for that reason; nothing here may collapse them.
    """

    __tablename__ = "pfam_family"
    __table_args__ = (
        # The page resolves chips by accession; a reader searching InterPro comes the other way.
        Index("ix_pfam_family__interpro_accession", "interpro_accession"),
    )

    #: `PF00126` — **version-stripped**, as `pfam_reference` strips it. An annotation carrying
    #: `PF00126.29` must be cut at the dot before it is looked up here or it silently misses.
    pfam_accession: Mapped[str] = mapped_column(String(16), primary_key=True)
    #: `Sigma70_r2` — what the chip says. Empty string where the table has none, never NULL: the
    #: vendored reader pads short rows rather than skipping them.
    short_name: Mapped[str] = mapped_column(String(128), nullable=False, default="")
    description: Mapped[str] = mapped_column(Text, nullable=False, default="")
    #: ⭐ Preferred over the Pfam entry for the link: it is the integrated record, and the page a
    #: reader following a domain actually wants. Empty where there is no integrated entry.
    interpro_accession: Mapped[str] = mapped_column(String(16), nullable=False, default="")
    interpro_name: Mapped[str] = mapped_column(Text, nullable=False, default="")
    clan_accession: Mapped[str] = mapped_column(String(16), nullable=False, default="")
    clan_name: Mapped[str] = mapped_column(String(128), nullable=False, default="")

    def __repr__(self) -> str:
        return f"<PfamFamily {self.pfam_accession} {self.short_name}>"


class EnzymeClass(Base):
    """One EC code and its name, at **every level ExPASy names** — not just the seven classes.

    Vendored from `nuna.tl.locus_browser.vendor_ec`, which is the single reader of
    `ec_names.tsv.gz`, built from two upstream files because neither covers the hierarchy alone:
    `enzclass.txt` names levels 1-3 and `enzyme.dat` names level 4.

    ⛔ **The key is a single code, and an annotation is a SET of them.** 1,032 of the catalogue's
    9,398 EC entries are comma-joined lists (`1.6.5.9,7.1.1.-`), so a caller must split on commas
    and strip each code before looking it up here. Splitting in SQL with `split_part` returns a
    code neither side holds — see the `ec-numbers-are-sets-not-strings` note.

    ⚠ **Dash-padded codes are rows like any other.** `2.7.10.-` is a real key with a real name,
    *Protein-tyrosine kinases*; 126 of the 1,410 codes in the two catalogues state fewer than four
    fields, and those are precisely the ones a bare number says least about.

    ⚠ **`is_withdrawn` rows are KEPT.** `name` may read *"Transferred entry: 1.1.1.303 and
    1.1.1.304"*, which is upstream's own text and a true, useful statement. Dropping those rows
    would render as *"this code has no name"*, which is a different and false one.
    """

    __tablename__ = "enzyme_class"

    #: `2.7.10.-`, `1.1.1.1` — dash-padded exactly as the annotation carries it, never normalised.
    ec_code: Mapped[str] = mapped_column(String(32), primary_key=True)
    #: `Protein-tyrosine kinases`. Trailing full stop removed; internal ones kept, because a
    #: transferred-entry message is made of EC codes and stripping every `.` would mangle them.
    name: Mapped[str] = mapped_column(Text, nullable=False, default="")

    def __repr__(self) -> str:
        return f"<EnzymeClass {self.ec_code} {self.name[:32]}>"


class KeggOrthology(Base):
    """One KEGG KO: the gene symbols it is known by, and what it does.

    ⛔⛔ **Serving this is subject to an unresolved licence** — see the module docstring. The table
    being present is a deployment decision, not a default: `kegg_reference()` returns `{}` where the
    vendored file is absent, and the card then falls back to a bare linked accession.

    ⚠ **A symbol is optional and sometimes plural.** `K01991` is `wza, gfcE`; 1,186 of 28,516 KOs
    have no symbol at all, and upstream uses the KO id itself as a placeholder, which the vendored
    reader reads as *no symbol* rather than inventing one. Empty string, never NULL.
    """

    __tablename__ = "kegg_orthology"

    #: `K01991`.
    ko_id: Mapped[str] = mapped_column(String(16), primary_key=True)
    #: `wza, gfcE` — a comma-joined list of gene symbols, or empty where upstream publishes none.
    symbol: Mapped[str] = mapped_column(String(256), nullable=False, default="")
    #: `polysaccharide biosynthesis/export protein`. Keeps upstream's own `[EC:...]` cross
    #: reference where it carries one: removing it is a silent edit of a reference table.
    definition: Mapped[str] = mapped_column(Text, nullable=False, default="")

    def __repr__(self) -> str:
        return f"<KeggOrthology {self.ko_id} {self.symbol or self.definition[:32]}>"


class GeneOntologyTerm(Base):
    """One GO term: what it is called, which namespace it belongs to, and the slim classes it folds onto.

    ⛔ **The 38,092-row table exists for ONE consumer**: within-node propagation, which reads
    `gene_functional_annotation.gene_ontology_terms` — raw accessions, all three namespaces pooled
    into one array. The neighbour ladder does not need it: `locus_annotation_entry.term_value` is
    already a slim class with its namespace already stored beside it.

    ⚠ `namespace_index`, not a namespace string — `enumerations.GENE_ONTOLOGY_NAMESPACE_NAMES`
    already owns the index↔name contract, and this joins straight to
    `locus_annotation_entry.gene_ontology_namespace`. A second spelling of the same fact is how two
    parts of one response come to disagree about which namespace a term is in.

    ⚠ **A term folds onto EVERY nearest slim class, not one.** `nearest_slim` keeps all of them
    because dropping the others would make two genes that share a dropped ancestor read as disjoint.
    So `slim_go_ids` is an array, and it can be empty — a term that reaches no slim class at all
    self-maps on read in nuna, and a consumer here must do the same or silently drop the gene.
    """

    __tablename__ = "gene_ontology_term"
    __table_args__ = (
        # The propagation path groups a locus's genes by namespace before comparing them.
        Index("ix_gene_ontology_term__namespace_index", "namespace_index"),
    )

    #: `GO:0005886`.
    go_id: Mapped[str] = mapped_column(String(16), primary_key=True)
    #: `plasma membrane`.
    name: Mapped[str] = mapped_column(Text, nullable=False, default="")
    #: 0 molecular_function · 1 biological_process · 2 cellular_component.
    namespace_index: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    #: Every NEAREST slim class, which is usually one and is sometimes none.
    slim_go_ids: Mapped[list[str]] = mapped_column(ARRAY(String(16)), nullable=False)

    def __repr__(self) -> str:
        return f"<GeneOntologyTerm {self.go_id} {self.name[:32]}>"


class GeneOntologySlimAncestor(Base):
    """One edge of the slim's own hierarchy: `go_id` is under `ancestor_go_id`, transitively.

    ⭐⭐ **This is the table that fixes the comparison, and the reason it has to exist at all.**
    `goslim_metagenomics` is a SELECTION, not a cut through the DAG: it holds `membrane`,
    `plasma membrane` and `outer membrane` as three separate classes. Folding onto it therefore put
    a term and its own parent in different buckets, and the set comparison that followed called them
    `disjoint` — its most condemning verdict. Klebsiella `gumC` read *"classes differ"* for 16 genes
    saying plasma membrane and one saying membrane.

    ⛔ **No row names a namespace root, by construction.** `cellular_component` is an ancestor of
    every cellular-component class, so keeping it would make every such pair agree and `disjoint`
    unreachable — which is `vendor_go`'s own standing objection to full ancestor closure, and it is
    correct. `slim_ancestry` strips the three roots on read and the ingest never writes them.

    ⚠ **TRANSITIVE, not parent-of.** `plasma membrane` names `membrane` whether that is one step up
    or three. A caller must not walk this table recursively; one hop is the whole answer.

    An edge table rather than an array column because the FKs make "an ancestor that is not a known
    class" impossible, and because the reverse direction — *which classes are under `membrane`?* — is
    a real query.
    """

    __tablename__ = "gene_ontology_slim_ancestor"
    __table_args__ = (Index("ix_gene_ontology_slim_ancestor__ancestor", "ancestor_go_id"),)

    go_id: Mapped[str] = mapped_column(String(16), ForeignKey("gene_ontology_term.go_id"), primary_key=True)
    ancestor_go_id: Mapped[str] = mapped_column(String(16), ForeignKey("gene_ontology_term.go_id"), primary_key=True)

    def __repr__(self) -> str:
        return f"<GeneOntologySlimAncestor {self.go_id} < {self.ancestor_go_id}>"
