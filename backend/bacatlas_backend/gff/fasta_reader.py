"""A FASTA block → ``{seqid: bases}``, in one place because two sources now feed it.

⛔ **One rule, two callers, and that is the point of the module.** The bases for a gene come either
from the ``##FASTA`` block of the Bakta GFF we hold, or from the AllTheBacteria assembly fetched for
that genome (``sequence_source/``). Those files are the same assembly — proved byte-identical on all
280 probe genomes, 2026-10-23 — but only if they are read by the same rule, and there are exactly
two places that rule can go wrong:

* **The seqid is the FIRST whitespace-delimited token**, never the whole defline. Our GFF writes a
  bare ``>SAMN04393397.contig00001``; ATB writes
  ``>SAMN04393397.contig00001 len=296336 cov=30.8 corr=0 origname=NODE_1_… sw=shovill-spades/1.0.4``.
  Keeping the defline would key every ATB contig under a string no database row can match, and the
  failure would read as "this contig is not in the assembly" — a sentence about the data rather than
  about the parse.
* **Bases are upper-cased here, once.** So nothing downstream has to decide whether a lower-case
  soft-masked base is a base, and no comparison between the two sources can differ by case alone.

⚠ Line wrapping is NOT preserved and must not be: our GFF wraps at 60, ATB at 70, and a reader that
kept the newlines would make two identical assemblies compare unequal.
"""

from __future__ import annotations

from typing import IO


def read_fasta_block(handle: IO[str]) -> dict[str, str]:
    """Every record in ``handle``, by seqid, upper-cased and unwrapped.

    The handle may be positioned anywhere a defline is next — ``parse_genome_annotation`` hands it
    over already past ``##FASTA``, and the ATB reader hands over a whole file.
    """
    sequences: dict[str, str] = {}
    seqid: str | None = None
    chunks: list[str] = []
    for line in handle:
        if line.startswith(">"):
            if seqid is not None:
                sequences[seqid] = "".join(chunks)
            # ⛔ First token only — see the module docstring. ATB's defline carries five more fields.
            seqid = line[1:].split()[0]
            chunks = []
        elif seqid is not None:
            chunks.append(line.strip().upper())
    if seqid is not None:
        sequences[seqid] = "".join(chunks)
    return sequences
