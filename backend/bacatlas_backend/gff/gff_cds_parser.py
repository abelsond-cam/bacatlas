"""One pass over a Bakta GFF → its CDS features and its contig sequences.

⛔ **One pass, not two.** The CDS lines and the ``##FASTA`` block are in the same file, and a
1.9 MB gzipped GFF costs more to decompress twice than to hold once. The reader therefore returns
both halves together and the caller decides what to keep.

⛔ **No pandas.** ``extract_strand.parse_gff_cds`` returns a DataFrame, which is right for an
ingest job and wrong here: pandas is in the backend's ``ingest`` extra, not its serving
dependencies, and the sequence endpoint runs on a machine that has neither pandas nor ``nuna``.

Coordinates are the GFF's own: **1-based inclusive**, and ``end`` includes the stop codon.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from bacatlas_backend.gff.fasta_reader import read_fasta_block
from bacatlas_backend.gff.gff_text_reader import open_gff_text

FASTA_DIRECTIVE = "##FASTA"

#: ⛔ **The `partial` attribute values Bakta uses for a 5'-partial CDS**, matching
#: `nuna.tl.locus_browser.genome_sequence._PARTIAL_5P`. `10` and `11` are the two-bit form whose
#: FIRST bit is the 5' end; `1` and `true` are the plain form.
#: ⚠ Not one of the 120,792 CDS across 25 probe GFFs carries a `partial=` attribute at all, so on
#: THIS cohort `phase != 0` alone is indistinguishable from the full rule. That is precisely why the
#: full rule is written here: the two would diverge silently on the first cohort where Bakta emits
#: it, and the divergence changes TRANSLATION — the initiator is promoted to `M` only when the CDS
#: is not 5'-partial.
FIVE_PRIME_PARTIAL_VALUES = frozenset({"true", "1", "10", "11"})

#: ⛔ **A `pseudo` CDS is skipped, exactly as the extractor skipped it.** The extractor that wrote
#: `protein_sequence` dropped them, so a feature kept here would have no protein to check against —
#: and, more importantly, `flat_index` is a running counter over the CDS the extractor KEPT. Admit
#: one it dropped and every index after it names a different gene. Measured: 204 pseudo CDS over 25
#: probe genomes.
PSEUDO_VALUES = frozenset({"true", "1"})


@dataclass(frozen=True)
class CodingFeature:
    """One CDS line. ``phase`` and ``is_five_prime_partial`` are what make translation reproducible."""

    seqid: str
    start_position: int          # 1-based inclusive
    end_position: int            # inclusive, INCLUDES the stop codon
    strand: str                  # '+' or '-'
    phase: int                   # 0/1/2; non-zero means the CDS does not begin at a start codon
    locus_tag: str | None
    is_five_prime_partial: bool


@dataclass(frozen=True)
class ParsedGenomeAnnotation:
    """A whole GFF: its CDS features in file order, and its contig sequences by seqid."""

    #: ⚠ **The CDS the EXTRACTOR would keep, in file order** — `pseudo` already dropped, so this is
    #: the same set `nuna.genome_sequence.parse_gff` returns. It is still NOT the gene table: the
    #: extractor additionally drops a CDS whose contig is missing, whose translation is empty, or
    #: which carries an internal stop, and those three need the sequence. `flat_index` is the index
    #: after ALL of them, which is why the ingest reproduces the whole chain rather than zipping.
    coding_features: tuple[CodingFeature, ...]
    contig_sequences: dict[str, str]

    @property
    def carries_sequence(self) -> bool:
        """Whether the file had a ``##FASTA`` block at all.

        ⚠ Checked rather than assumed. Measured 2026-09-04: all 280 probe GFFs carry one, but the
        endpoint must fail with a named reason on one that does not, never with a KeyError.
        """
        return bool(self.contig_sequences)


def _attributes(field: str) -> dict[str, str]:
    """GFF3 column 9 → a dict, **values URL-decoded**. Malformed pairs are skipped, not guessed at.

    ⛔ **The decode is not optional and it is not cosmetic.** GFF3 percent-encodes the separators
    (`;`, `=`, `,`) inside a value, and Bakta products are full of them: measured over 40 probe
    GFFs, **12,210 CDS lines carry a `%XX` escape** — every
    `UDP-N-acetylmuramoyl-L-alanyl-D-glutamate--2%2C6-diaminopimelate ligase` is one. Left encoded,
    a product reads `2%2C6` on the page and, worse, fails to match the same product read through
    `nuna`'s reader, which has always called `unquote` (`genome_sequence._attrs`).

    ⚠ This diverged silently because the cross-check compares only the fields nuna's `parse_gff`
    RETURNS — seqid, coordinates, strand, phase, partiality — and attributes are not among them. A
    copy is only safe while something fails when it diverges; here nothing did, so the test now
    covers the attributes too.
    """
    out: dict[str, str] = {}
    for item in field.rstrip(";").split(";"):
        key, separator, value = item.partition("=")
        if separator:
            out[key.strip()] = unquote(value.strip())
    return out


def parse_genome_annotation(path: Path, *, want_sequence: bool = True) -> ParsedGenomeAnnotation:
    """Read one GFF. Set ``want_sequence=False`` to stop at ``##FASTA`` and skip the bases."""
    features: list[CodingFeature] = []
    sequences: dict[str, str] = {}

    with open_gff_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith(FASTA_DIRECTIVE):
                    break
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9 or columns[2] != "CDS":
                continue
            attributes = _attributes(columns[8])
            if attributes.get("pseudo", "").lower() in PSEUDO_VALUES:
                continue
            phase = int(columns[7]) if columns[7].isdigit() else 0
            features.append(
                CodingFeature(
                    seqid=columns[0],
                    start_position=int(columns[3]),
                    end_position=int(columns[4]),
                    strand=columns[6],
                    # A '.' phase is 0 by GFF3 convention; Bakta writes an integer, but a
                    # reader that assumes so raises on a file it should have handled.
                    phase=phase,
                    locus_tag=attributes.get("locus_tag"),
                    # ⛔ BOTH signals, matching nuna's rule. A non-zero phase implies 5'-partial,
                    # but the converse does not hold: Bakta can mark a CDS partial with phase 0.
                    is_five_prime_partial=(
                        phase != 0 or attributes.get("partial", "").lower() in FIVE_PRIME_PARTIAL_VALUES
                    ),
                )
            )
        else:
            # Loop finished without hitting ##FASTA: the file carries no sequence.
            return ParsedGenomeAnnotation(tuple(features), {})

        if want_sequence:
            sequences = read_fasta_block(handle)

    return ParsedGenomeAnnotation(tuple(features), sequences)
