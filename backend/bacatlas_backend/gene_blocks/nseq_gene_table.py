"""Reading nuna's `.nseq` and `.loci` — the external oracle for the packed blocks.

⭐ **This exists to be a differential oracle, not to serve a request.** `gene_spans` claims to be
byte-identical to the `.nseq` gene table and `gene_locus_ordinals` to the `.loci` vector. A claim
like that is worth nothing unless something reads both and compares, so this reads the files nuna
already published — written by `nuna.tl.locus_browser.genome_sequence.encode_nseq`, a different
author, a different language on the page side, and frozen in git — and hands back the raw gene-table
bytes. `tests/test_gene_blocks.py` compares them to what `packed_gene_blocks` produces from the
decoded fields, over all 200 committed genomes.

**Nothing here is on the serving path** and nothing here decodes bases. The 2-bit stream is the
page's business; the database holds no bases and this module does not change that. Only the header,
the gene table and the locus vector are read.

⚠ **`app.js` validates seven bytes of the eight-byte magic.** `decodeNseq` tests
`utf8(u8.subarray(0, 7)) !== "NUNASEQ"`, so the trailing container-version byte is never looked at —
while `nuna/schemas/genome_sequence.schema.md` says *"the page tests the version to decide whether
it can read the file."* The doc describes an intention the code does not implement. This reader
checks all eight, because a container-version bump is exactly the change that must not be read as
though nothing happened.
"""

from __future__ import annotations

import json
import struct
from array import array
from dataclasses import dataclass
from pathlib import Path

from bacatlas_backend.gene_blocks.packed_gene_blocks import (
    GENE_SPAN_STRIDE,
    INT32_VECTOR_STRIDE,
    GeneBlockError,
    check_block_length,
    decode_int32_vector,
)

#: The eight-byte file magic. The trailing byte is the CONTAINER version, bumped when the layout
#: changes and never when the contents do.
NSEQ_MAGIC = b"NUNASEQ\x01"

#: Where the `uint32` header length sits, and where the JSON header starts.
_HEADER_LENGTH_OFFSET = 8
_HEADER_OFFSET = 12


@dataclass(frozen=True)
class NseqFile:
    """One `.nseq`, read as far as the gene table — the header, the table's bytes, and nothing else."""

    path: Path
    sample: str
    contig_names: tuple[str, ...]
    contig_lengths: tuple[int, ...]
    gene_count: int
    base_count: int
    ambiguous_runs: tuple[tuple[int, int], ...]
    #: ⭐ The gene table exactly as published: `gene_count * 16` bytes, four little-endian `uint32`
    #: per gene. This is the oracle `gene_spans` is compared against, so it is kept as bytes and
    #: never as decoded fields — decoding it here would make the comparison compare one decoder
    #: against itself.
    gene_table: bytes
    packed_base_bytes: int

    def contig_base_offsets(self) -> tuple[int, ...]:
        """The prefix sum of `contig_lengths` — `app.js`'s `coff`, which `.nseq` deliberately omits."""
        offsets = [0]
        for length in self.contig_lengths:
            offsets.append(offsets[-1] + length)
        return tuple(offsets)


def read_nseq(path: str | Path, *, expect_sample: str | None = None) -> NseqFile:
    """Read one `.nseq` up to and including its gene table.

    ⭐ `expect_sample` is what makes the header's `sample` key do the job the format spec gives
    it — *"the BioSample accession, so a mis-served file names itself"*. Without it the key is
    read and never compared to anything, so a file published under the wrong name decodes
    perfectly and serves another genome's coordinates. Pass the name you asked for.

    ⚠ Every length is cross-checked against the file's own size before a slice is taken. A truncated
    download is the realistic failure, and a short gene table sliced without checking yields a table
    that is a whole number of rows and simply stops early.
    """
    path = Path(path)
    raw = path.read_bytes()
    if len(raw) < _HEADER_OFFSET:
        raise GeneBlockError(f"{path}: {len(raw)} bytes is shorter than the 12-byte container header")
    if raw[:8] != NSEQ_MAGIC:
        raise GeneBlockError(
            f"{path}: magic is {raw[:8]!r}, not {NSEQ_MAGIC!r}. ⚠ A container-version bump changes "
            "the eighth byte, and the layout with it."
        )

    (header_length,) = struct.unpack_from("<I", raw, _HEADER_LENGTH_OFFSET)
    header_end = _HEADER_OFFSET + header_length
    if header_end > len(raw):
        raise GeneBlockError(
            f"{path}: the header claims {header_length} bytes, which runs past the file's {len(raw)}"
        )
    header = json.loads(raw[_HEADER_OFFSET:header_end].decode("utf-8"))

    gene_count = int(header["n_genes"])
    table_end = header_end + gene_count * GENE_SPAN_STRIDE
    if table_end > len(raw):
        raise GeneBlockError(
            f"{path}: n_genes {gene_count} needs a gene table ending at {table_end}, past the "
            f"file's {len(raw)} bytes"
        )
    gene_table = raw[header_end:table_end]
    check_block_length(gene_table, GENE_SPAN_STRIDE, gene_count, name=f"{path.name} gene table")

    contig_names = tuple(str(name) for name in header["cname"])
    contig_lengths = tuple(int(length) for length in header["clen"])
    if len(contig_names) != len(contig_lengths):
        raise GeneBlockError(
            f"{path}: {len(contig_names)} contig names against {len(contig_lengths)} lengths — "
            "`cname` and `clen` are parallel and in `contig_idx` order"
        )
    base_count = int(header["n_bases"])
    if base_count != sum(contig_lengths):
        raise GeneBlockError(
            f"{path}: n_bases {base_count} is not the sum of clen {sum(contig_lengths)}; the "
            "contig base offsets are the prefix sum of clen, so the two disagreeing moves every base"
        )
    expected_packed = -(-base_count // 4)
    packed = len(raw) - table_end
    if packed != expected_packed:
        raise GeneBlockError(
            f"{path}: {packed} packed bytes after the gene table, but n_bases {base_count} needs "
            f"{expected_packed} at four bases a byte"
        )

    sample = str(header["sample"])
    if expect_sample is not None and sample != expect_sample:
        raise GeneBlockError(
            f"{path}: the header names sample {sample!r} but {expect_sample!r} was asked for. "
            "⛔ A file published under the wrong name decodes perfectly and serves another "
            "genome's coordinates; this is the only check that can see it."
        )

    return NseqFile(
        path=path,
        sample=sample,
        contig_names=contig_names,
        contig_lengths=contig_lengths,
        gene_count=gene_count,
        base_count=base_count,
        ambiguous_runs=tuple((int(start), int(length)) for start, length in header.get("ambiguous", [])),
        gene_table=gene_table,
        packed_base_bytes=packed,
    )


def read_loci(path: str | Path, *, gene_count: int | None = None) -> array:
    """Read one `.loci` — `n_genes` little-endian `int32`, `-1` where the gene is not in the catalogue.

    ⭐ This is the same 4 B/gene layout as `genome_locus_assignment.gene_locus_ordinals`, which is
    why that column can be compared against a file nuna published rather than only against itself.
    """
    path = Path(path)
    raw = path.read_bytes()
    return decode_int32_vector(raw, name=f"{path.name} (.loci)", expected_count=gene_count)


def loci_stride() -> int:
    """Bytes per gene in a `.loci` file — named so a test can state the layout it relies on."""
    return INT32_VECTOR_STRIDE
