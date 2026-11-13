"""The coordinate layer's contig names → an assembly's headers: exactly, or by ONE rule proved per genome.

⛔ **The namespace split is real and it is silent.** The long-read slice's gene coordinates were built from
GenBank-style GFFs whose seqids read ``CP034045.1``; the RefSeq assembly NCBI serves for the same genome
carries headers ``NZ_CP034045.1`` (measured on the live API, 2026-11-13: ``GCF_000512165.1``'s deflines read
``>NZ_CP003999.1 Klebsiella pneumoniae …``). Every contig lookup then misses, and the endpoint would say
*"contig CP034045.1 is not in the assembly"* — a sentence about the data that is really about a prefix.

⭐ **The rule is nuna's, with the same behaviour** —
``nuna.tl.locus_browser.genome_sequence.resolve_contig_names``, the function the 71,435-genome ``.nseq``
build ran on every one of the 4,806 long-read genomes (`PROJECT_STATE.md` §Layer 6, stage 5: *"L2's contig
mapping was right on all 4,806 — 0 unresolved seqids"*). Exact match first; then a prefix rule — ``NZ_``
added, else ``NZ_`` stripped — applied to the RESIDUE only, accepted only if it resolves *every* residual
name and leaves the whole map injective; otherwise refuse, naming the first unresolved contig.

⚠ **Proved per genome, never trusted.** The rule is not "strip ``NZ_``": a bare strip would happily map a
name onto a header another name already owns, and the bases would then come from the wrong contig at
correct-looking coordinates — the one failure no page can detect. Injectivity is the guard, and refusing
by name is the honest answer when it does not hold. Nuna had a second oracle (the stored protein had to
reproduce from the chosen bases); the Sequence tab has none, which is why this refuses rather than guesses.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Final

#: The one prefix the record names. ⚠ Not ``NC_``: no observed failure carried it, and a rule nobody
#: measured is a rule that can be wrong silently. Add a second only with the genome that needed it.
REFSEQ_PREFIX: Final = "NZ_"

#: In the order nuna tries them. Each is applied to the names that did NOT match exactly.
_PREFIX_RULES: Final[tuple[tuple[str, Callable[[str], str]], ...]] = (
    (f"+{REFSEQ_PREFIX}", lambda name: REFSEQ_PREFIX + name),
    (f"-{REFSEQ_PREFIX}", lambda name: name.removeprefix(REFSEQ_PREFIX)),
)


class ContigNamespaceMismatch(LookupError):
    """The assembly's headers do not account for every contig the coordinates name — permanently.

    ⛔ Not retryable: the same assembly will come back from the same accession tomorrow. The message
    names the first unresolved contig, so the failure reads as what it is (a coordinate/assembly
    mismatch on a named contig) and not as a missing file.
    """


def resolve_contig_names(wanted: Iterable[str], available: Iterable[str]) -> dict[str, str]:
    """``{wanted name: assembly header}`` for every wanted name, or raise naming the first that has none.

    Exact matches are taken first. Only the residue is then offered to the prefix rules, in order:
    ``NZ_`` added, then ``NZ_`` stripped. A rule is accepted only if it resolves the WHOLE residue and
    the full map stays injective — two wanted names landing on one header is the wrong-contig failure
    this function exists to refuse. An empty ``wanted`` resolves to an empty map.
    """
    needed = set(wanted)
    headers = set(available)
    chosen = {name: name for name in needed if name in headers}
    missing = sorted(needed - chosen.keys())
    if not missing:
        return chosen
    # What each rule could NOT place, so the refusal names the contig that is actually the problem —
    # not merely the first that failed the exact match, which the rule may have resolved fine.
    unplaced_by_rule: list[list[str]] = []
    for _rule, apply in _PREFIX_RULES:
        mapped = {name: apply(name) for name in missing}
        absent = sorted(name for name, header in mapped.items() if header not in headers)
        if absent:
            unplaced_by_rule.append(absent)
            continue
        full = chosen | mapped
        if len(set(full.values())) == len(full):
            return full
        taken: dict[str, str] = {}
        colliding: set[str] = set()
        for name, header in full.items():
            if header in taken:
                colliding.update((name, taken[header]))
            taken[header] = name
        unplaced_by_rule.append(sorted(colliding))
    unplaced = min(unplaced_by_rule, key=len)
    raise ContigNamespaceMismatch(
        f"{len(unplaced)} of {len(needed)} contig names have no unambiguous assembly header, exactly or by "
        f"the {REFSEQ_PREFIX} prefix rule; the first is {unplaced[0]!r}"
    )
