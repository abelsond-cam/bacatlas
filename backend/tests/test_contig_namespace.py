"""`resolve_contig_names` — exact, then ONE `NZ_` rule proved on the whole genome, else refuse by name.

⛔ The failure this guards is silent: a bare `NZ_` strip that maps one contig name onto a header another
name already owns hands back the wrong contig's bases at correct-looking coordinates. Every refusal here
is a case where guessing would have produced plausible DNA.
"""

from __future__ import annotations

import pytest

from bacatlas_backend.sequence_source.contig_namespace import ContigNamespaceMismatch, resolve_contig_names


def test_exact_matches_need_no_rule():
    assert resolve_contig_names(["CP1.1", "CP2.1"], ["CP2.1", "CP1.1", "CP3.1"]) == {
        "CP1.1": "CP1.1",
        "CP2.1": "CP2.1",
    }


def test_NZ_is_ADDED_when_the_assembly_is_RefSeq_and_the_coordinates_are_GenBank():
    """⭐ The measured case: GFF `CP034045.1`, assembly header `NZ_CP034045.1`."""
    assert resolve_contig_names(["CP1.1", "CP2.1"], ["NZ_CP1.1", "NZ_CP2.1"]) == {
        "CP1.1": "NZ_CP1.1",
        "CP2.1": "NZ_CP2.1",
    }


def test_NZ_is_STRIPPED_when_the_namespaces_run_the_other_way():
    assert resolve_contig_names(["NZ_CP1.1"], ["CP1.1"]) == {"NZ_CP1.1": "CP1.1"}


def test_MIXED_exact_and_prefixed_names_resolve_per_contig():
    """⚠ The rule applies to the RESIDUE, so an assembly in two namespaces resolves contig by contig."""
    assert resolve_contig_names(["plasmid_A", "CP1.1"], ["plasmid_A", "NZ_CP1.1"]) == {
        "plasmid_A": "plasmid_A",
        "CP1.1": "NZ_CP1.1",
    }


def test_an_AMBIGUOUS_map_is_refused_rather_than_guessed():
    """⛔ `CP1.1` and `NZ_CP1.1` both wanted, one header: the rule would land both on it. Injectivity is
    the guard, and the refusal names the contig whose bases would otherwise have been wrong."""
    with pytest.raises(ContigNamespaceMismatch) as raised:
        resolve_contig_names(["CP1.1", "NZ_CP1.1"], ["NZ_CP1.1"])
    assert "'CP1.1'" in str(raised.value)


def test_ONE_missing_name_refuses_the_whole_genome_BY_NAME():
    """⛔ All-or-nothing: a rule that resolves four of five contigs is not accepted for the four."""
    with pytest.raises(ContigNamespaceMismatch) as raised:
        resolve_contig_names(["CP1.1", "CP2.1", "CP9.1"], ["NZ_CP1.1", "NZ_CP2.1"])
    message = str(raised.value)
    assert "'CP9.1'" in message, message
    assert "1 of 3" in message


def test_a_name_that_matches_nothing_either_way_is_named():
    with pytest.raises(ContigNamespaceMismatch, match="'contig00001'"):
        resolve_contig_names(["contig00001"], ["NZ_CP1.1"])


def test_an_empty_request_resolves_to_an_empty_map():
    assert resolve_contig_names([], ["NZ_CP1.1"]) == {}


def test_the_result_is_INJECTIVE_whenever_it_is_returned():
    """A property every accepted map must have, stated once so a new rule cannot quietly break it."""
    cases = [
        (["CP1.1", "CP2.1"], ["NZ_CP1.1", "NZ_CP2.1"]),
        (["NZ_CP1.1", "X"], ["CP1.1", "X"]),
        (["a", "b"], ["a", "b", "c"]),
    ]
    for wanted, available in cases:
        resolved = resolve_contig_names(wanted, available)
        assert len(set(resolved.values())) == len(resolved), resolved
        assert set(resolved) == set(wanted)
