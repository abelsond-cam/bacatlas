"""Naming the two vocabularies Bakta leaves as bare codes — and the one shape that bites.

`term_name` is populated for 100 % of the catalogue's COG and GO entries and **0 %** of its EC and
KEGG ones, so these names arrive by a serve-time join. The join itself is trivial. What is not is
that **an EC value is a SET**: 1,032 of the catalogue's 9,398 EC entries are comma-joined lists, so
every lookup has to split first, and a partial answer is worse than none.

⛔ These are pure tests of `ReferenceNames` — no database. The dicts it holds are exactly what the
two tables return, and building them by hand is what lets the awkward cases be written down at all.
"""

from __future__ import annotations

import pytest

from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.services.reference_name_service import ReferenceNames

ENZYME = {
    "2.7.10.-": "Protein-tyrosine kinases",
    "2.-.-.-": "Transferases",
    "1.1.1.1": "alcohol dehydrogenase",
    "2.7.7.-": "Nucleotidyltransferases",
    "1.1.1.5": "Transferred entry: 1.1.1.303 and 1.1.1.304",
}
KEGG = {"K01991": "polysaccharide biosynthesis/export protein", "K99999": ""}


@pytest.fixture
def names() -> ReferenceNames:
    return ReferenceNames(enzyme=dict(ENZYME), kegg=dict(KEGG))


def test_the_code_David_asked_for_by_name(names):
    """⭐ `https://enzyme.expasy.org/EC/2.7.10.-` is the page he gave; `gumC` carries exactly it."""
    assert names.enzyme_name("2.7.10.-") == "Protein-tyrosine kinases"


def test_a_dash_padded_code_is_an_ORDINARY_key_not_a_fallback(names):
    # 126 of the 1,410 codes in the two catalogues state fewer than four fields. Treating them as
    # "level 4 missing, fall back to the class word" is what produced "transferase · 4th level not
    # stated" where ExPASy had a name all along.
    assert names.enzyme_name("2.-.-.-") == "Transferases"


def test_a_comma_joined_value_names_EVERY_code(names):
    # ⛔ An EC value is a set. `split_part` in SQL would return a code neither side holds.
    assert names.enzyme_name("1.1.1.1,2.7.7.-") == "alcohol dehydrogenase · Nucleotidyltransferases"


def test_whitespace_after_the_comma_is_stripped(names):
    # Measured in the live catalogue: `1.1.1.-` and `1.1.1.- ` are both present as distinct values.
    assert names.enzyme_name(" 1.1.1.1 , 2.7.7.- ").startswith("alcohol dehydrogenase")


def test_a_repeated_code_is_named_ONCE(names):
    assert names.enzyme_name("1.1.1.1,1.1.1.1") == "alcohol dehydrogenase"


def test_ONE_unnamed_code_leaves_the_WHOLE_value_unnamed(names):
    """⛔⛔ The gate that matters most here.

    A partial join reads as the complete answer, and the half it keeps is whichever happened to be
    first. *"Nucleotidyltransferases"* under `2.7.7.-,9.9.9.9` would state that this locus catalyses
    one reaction when the annotation names two.
    """
    assert names.enzyme_name("2.7.7.-,9.9.9.9") is None
    assert names.enzyme_name("9.9.9.9") is None


def test_a_withdrawn_entry_is_NAMED_with_upstream_s_own_words(names):
    # "Transferred entry: …" tells a reader where the annotation they are looking at went. Dropping
    # it would render as "this code has no name", which is a different and false statement.
    assert names.enzyme_name("1.1.1.5") == "Transferred entry: 1.1.1.303 and 1.1.1.304"


def test_the_kegg_definition_is_the_name_not_the_symbol(names):
    # The symbol (`wza, gfcE`) repeats gene names the page already shows a few hundred pixels away;
    # the definition is the thing a bare accession withholds.
    assert names.kegg_name("K01991") == "polysaccharide biosynthesis/export protein"


def test_an_empty_definition_is_NOT_a_name(names):
    # "" would render as a named term with nothing after it — worse than the bare accession.
    assert names.kegg_name("K99999") is None
    assert names.kegg_name("K00000") is None


def test_name_for_answers_only_for_the_two_vocabularies_it_owns(names):
    assert names.name_for(AnnotationKind.EC_NUMBER, "2.7.10.-") == "Protein-tyrosine kinases"
    assert names.name_for("ec_number", "2.7.10.-") == "Protein-tyrosine kinases"
    assert names.name_for(AnnotationKind.KEGG_ORTHOLOGY, "K01991") is not None
    # ⛔ COG and GO carry `term_name` from the ingest; answering for them here would put a second,
    # differently-sourced name in play for terms that already have one.
    assert names.name_for(AnnotationKind.COG_ORTHOGROUP, "COG3206") is None
    assert names.name_for(AnnotationKind.GENE_ONTOLOGY_SLIM, "GO:0016301") is None
    assert names.name_for(AnnotationKind.EC_NUMBER, None) is None


def test_an_UNLOADED_reference_names_nothing_and_raises_nothing():
    """⛔ The normal state of any deployment that has not accepted the KEGG licence position.

    The card must then fall back to a bare linked accession — exactly what the published static
    pages do permanently — rather than fail the request.
    """
    empty = ReferenceNames()
    assert empty.enzyme_name("2.7.10.-") is None
    assert empty.kegg_name("K01991") is None
    assert empty.name_for(AnnotationKind.EC_NUMBER, "2.7.10.-") is None
