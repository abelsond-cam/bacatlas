"""`_claims` reads only the loci that can fold to a claim — held against loci built to tell the cases apart.

⛔ **Why a probe fixture and not the loaded catalogues.** `test_function_inference`'s full-scan oracle
compares the new read with the old one on the published catalogues, and in every one of them
`modal_cog_categories` is present on exactly the loci that carry a COG rank-0 entry (kp-71k: 83,134 and
83,134). There the outer join `_claims` uses for COG gives the same rows as an inner join, and as a
filter on the categories alone — two wrong reads the oracle cannot tell from the right one. Here each
case exists on purpose: categories with no entry, an entry with no categories, neither, an empty
category list, a three-namespace GO claim, and an EC value that resolves to no level.

Runs against the probe database (`conftest.engine`), which the default suite provides.
"""

from __future__ import annotations

from bacatlas_backend.instruments.annotation_transfer import GENE_ONTOLOGY_FACETS, _claims
from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry

from .conftest import make_locus


def _entry(locus, kind, term, *, namespace=None):
    return LocusAnnotationEntry(
        locus_id=locus.locus_id,
        annotation_kind=kind,
        rank_within_locus=0,
        term_value=term,
        term_name=None,
        member_gene_count=5,
        gene_ontology_namespace=namespace,
    )


def test_claims_keeps_exactly_the_loci_that_fold_to_a_claim(session, seeded):
    loci = {
        "categories_only": make_locus(seeded, ordinal=300, label="c300", modal_cog_categories=["E"]),
        "cog_entry_only": make_locus(seeded, ordinal=301, label="c301", modal_cog_categories=None),
        "bare": make_locus(seeded, ordinal=302, label="c302", modal_cog_categories=None),
        "empty_categories": make_locus(seeded, ordinal=303, label="c303", modal_cog_categories=[]),
        "go_three_namespaces": make_locus(seeded, ordinal=304, label="c304", modal_cog_categories=None),
        "ec_unresolvable": make_locus(seeded, ordinal=305, label="c305", modal_cog_categories=None),
        "ec_resolvable": make_locus(seeded, ordinal=306, label="c306", modal_cog_categories=None),
    }
    session.add_all(loci.values())
    session.flush()
    session.add_all(
        [
            _entry(loci["cog_entry_only"], AnnotationKind.COG_ORTHOGROUP, "COG0001"),
            _entry(loci["go_three_namespaces"], AnnotationKind.GENE_ONTOLOGY_SLIM, "GO:0003674", namespace=0),
            _entry(loci["go_three_namespaces"], AnnotationKind.GENE_ONTOLOGY_SLIM, "GO:0008150", namespace=1),
            _entry(loci["go_three_namespaces"], AnnotationKind.GENE_ONTOLOGY_SLIM, "GO:0005575", namespace=2),
            _entry(loci["ec_unresolvable"], AnnotationKind.EC_NUMBER, "-.-.-.-"),
            _entry(loci["ec_resolvable"], AnnotationKind.EC_NUMBER, "2.7.7.7"),
        ]
    )
    session.flush()
    pangenome_id = seeded["pangenome"].pangenome_id
    ids = {name: locus.locus_id for name, locus in loci.items()}

    cog, _terms, cog_count = _claims(session, pangenome_id, AnnotationKind.COG_ORTHOGROUP)
    # ⭐ The two cases the loaded catalogues cannot hold apart. An inner join loses the first; a
    # filter on the categories alone loses the second.
    assert ids["categories_only"] in cog, "a locus with categories and no COG entry still states the category rung"
    assert cog[ids["categories_only"]] == {(None, 2): frozenset(), (None, 1): frozenset({"E"})}
    assert ids["cog_entry_only"] in cog, "a locus with a COG entry and no categories still states the orthogroup"
    assert cog[ids["cog_entry_only"]] == {(None, 2): frozenset({"COG0001"}), (None, 1): frozenset()}
    assert set(cog) == {ids["categories_only"], ids["cog_entry_only"]}, "a bare or empty-category locus states nothing"

    go, _terms, _count = _claims(session, pangenome_id, AnnotationKind.GENE_ONTOLOGY_SLIM)
    assert set(go) == {ids["go_three_namespaces"]}
    # Rank 0 is not unique per locus for GO: all three namespaces' rows must reach the fold (with no
    # fold, `go_rung` pools them under the first facet).
    assert go[ids["go_three_namespaces"]][(GENE_ONTOLOGY_FACETS[0], 1)] == frozenset(
        {"GO:0003674", "GO:0008150", "GO:0005575"}
    )

    ec, _terms, _count = _claims(session, pangenome_id, AnnotationKind.EC_NUMBER)
    assert set(ec) == {ids["ec_resolvable"]}, "an EC value with no resolved level is not a claim"

    kegg, _terms, kegg_count = _claims(session, pangenome_id, AnnotationKind.KEGG_ORTHOLOGY)
    assert kegg == {}
    # ⚠ The count is of EVERY locus of the catalogue, the bare ones included — and it must survive a
    # vocabulary that returns no row at all, where the scalar subquery has nothing to ride on.
    assert cog_count == kegg_count == len(loci)
