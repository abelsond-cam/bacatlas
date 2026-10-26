"""Fixed-stride binary blocks — the gene layer at one row per genome, not one row per gene.

⛔ **STDLIB ONLY, and that is a hard constraint, not a preference.** These blocks are decoded on
the SERVING path, and `pyproject.toml`'s serving dependencies are flask, sqlalchemy, psycopg and
alembic. `numpy` is in the `ingest` extra and is **not installed in `backend/Dockerfile`'s image** —
the two existing `import numpy` lines in this package are both function-local, inside `ingest/`. A
module-scope numpy import here would pass every test on a laptop that has the ingest extra and fail
the Sequence tab in production, on the one line that reads a gene's coordinates.

⭐ **`gene_spans` is byte-identical to nuna's `.nseq` gene table**, not merely the same fields:
four little-endian `uint32` per gene — `start · end · contig · flags` — with row *j* being
`flat_index` *j* (`nuna/schemas/genome_sequence.schema.md`). That is the whole reason the format was
chosen. `app.js` already decodes those bytes, so its decoder and this one are differential oracles
for each other rather than two readings of one author's intent, and a symmetric encode/decode bug
cannot cancel. `nseq_gene_table.py` is what makes the comparison runnable.

**Columns in, bytes out; bytes in, O(1) accessors out.** Nothing here materialises an object per
gene in either direction. Ingest holds parquet columns and would pay 380 M dataclass constructions
at 72,801 genomes; the serving path needs one gene's fields out of ~4,900 and would pay the same
cost to read one row. So the encoders take parallel sequences and the decoders wrap the buffer —
`array("I", block)` is one C-level copy with no Python object per element, exactly what
`new Uint32Array(buf.slice(...))` is on the page side.

**Every block is `stride × gene_count` bytes and nothing else** — no count prefix, no terminator,
no padding, no alignment slack. The count lives in `genome_genes.gene_count` and the length is a
CHECK constraint, so a truncated block is a constraint violation at write time rather than a short
read at serve time. `check_block_length` is that CHECK, in Python, for the paths that never touch
Postgres.
"""

from __future__ import annotations

import sys
from array import array
from collections.abc import Iterable, Sequence

#: Bytes per gene in each block. A block is `stride * gene_count` bytes, exactly.
GENE_SPAN_STRIDE = 16
GENE_NAME_REF_STRIDE = 12
INT32_VECTOR_STRIDE = 4

#: `genome_genes.span_encoding`. **1 is nuna's `.nseq` gene-table layout** and is the only value
#: this module writes. It is stored rather than assumed so a future layout can be told from this one
#: by reading a row, not by knowing which release wrote it.
SPAN_ENCODING_NSEQ_GENE_TABLE = 1

#: The `flags` word, bit for bit as `.nseq` defines it and `app.js` reads it.
FLAG_MINUS_STRAND = 1 << 0
FLAG_FIVE_PRIME_PARTIAL = 1 << 1
FLAG_PHASE_SHIFT = 2
FLAG_PHASE_MASK = 0b11
#: ⚠ Every bit the layout defines. A set bit outside this mask is refused on the way IN rather than
#: ignored on the way out: an unknown bit means the writer knew something this reader does not, and
#: silently dropping it is how two versions of a format diverge without an error.
DEFINED_FLAG_BITS = FLAG_MINUS_STRAND | FLAG_FIVE_PRIME_PARTIAL | (FLAG_PHASE_MASK << FLAG_PHASE_SHIFT)

#: Absence, per column. ⭐ Each is a value the real domain cannot take, so no separate NULL mask is
#: needed: a term id is a positive identity column, a catalogue ordinal is a position, a rank counts
#: from 0.
ABSENT_TERM_ID = -1
ABSENT_LOCUS_TAG_NUMBER = -1
UNASSIGNED_LOCUS_ORDINAL = -1
NO_ARRANGEMENT_RANK = -1

UINT32_MAX = 0xFFFFFFFF
INT32_MIN = -0x8000_0000
INT32_MAX = 0x7FFF_FFFF

_IS_LITTLE_ENDIAN = sys.byteorder == "little"

# ⛔ C does not fix the width of `unsigned int`, so `array`'s typecodes are checked once here
# rather than trusted. On every platform this service runs on both are 4 bytes;
# if that ever stops being true the blocks would be read at the wrong stride, which is a plausible
# wrong answer rather than a failure.
if (array("I").itemsize, array("i").itemsize) != (4, 4):
    raise ImportError(
        "gene_blocks needs array('I') and array('i') to be 4 bytes; this platform reports "
        f"{array('I').itemsize} and {array('i').itemsize}. "
        "The blocks are a fixed on-disk layout, so a different word size would silently misread "
        "every row rather than fail."
    )


class GeneBlockError(ValueError):
    """A block's bytes do not say what the schema claims — named with the block and the row."""


# ── the length invariant, which is the CHECK constraint in Python ─────────────────────────────
def byte_length(block: bytes | memoryview) -> int:
    """How many BYTES a buffer holds — which is not what `len()` returns for every buffer.

    ⛔ `len()` on a `memoryview` counts **items, not bytes**. `memoryview(array("I", block))`
    reports `itemsize` 4, so `len()` is a quarter of the byte length and the stride arithmetic ends
    up measuring a different quantity than the decode does. Measured 2026-10-26: a 4,900-gene
    `gene_spans` block passed as such a view counts as **1,225** genes, and — the quiet direction —
    a 19,600-gene block is **accepted** as 4,900 and then read as 19,600. Both are exactly the
    "whole genome of plausible wrong coordinates with no error anywhere" that `check_block_length`
    exists to refuse, arriving through the check itself.
    """
    return block.nbytes if isinstance(block, memoryview) else len(block)


def block_gene_count(block: bytes | memoryview, stride: int, *, name: str) -> int:
    """How many genes `block` holds, refusing a length that is not a whole number of rows."""
    if stride <= 0:
        raise GeneBlockError(f"{name}: stride must be positive, not {stride}")
    length = byte_length(block)
    if length % stride:
        raise GeneBlockError(
            f"{name}: {length} bytes is not a whole number of {stride}-byte rows "
            f"({length / stride:.2f} genes). The schema's own CHECK is "
            f"octet_length({name}) = gene_count * {stride}."
        )
    return length // stride


def check_block_length(block: bytes | memoryview, stride: int, gene_count: int, *, name: str) -> int:
    """Refuse a block whose length disagrees with `gene_count`, as the CHECK constraint does.

    ⛔ Called by every decoder that is given a count. Without it a block one row short is read one
    gene short, and every gene after the missing one is read as the one before it — a whole genome
    of plausible wrong coordinates with no error anywhere.
    """
    found = block_gene_count(block, stride, name=name)
    if found != gene_count:
        raise GeneBlockError(
            f"{name}: {byte_length(block)} bytes is {found} genes at {stride} B/gene, but "
            f"gene_count says {gene_count}."
        )
    return found


# ── the two vector primitives every block is built from ───────────────────────────────────────
def _decode_vector(
    typecode: str,
    block: bytes | memoryview,
    *,
    name: str,
    stride: int,
    expected_count: int | None,
) -> array:
    """`block` → an `array` of `typecode`, little-endian whatever the CPU is."""
    if expected_count is None:
        block_gene_count(block, stride, name=name)
    else:
        check_block_length(block, stride, expected_count, name=name)
    values = array(typecode)
    values.frombytes(bytes(block))
    if not _IS_LITTLE_ENDIAN:
        values.byteswap()
    return values


def _encode_vector(
    typecode: str,
    values: Iterable[int],
    *,
    name: str,
    low: int,
    high: int,
    fields_per_gene: int = 1,
) -> bytes:
    """An iterable of ints → little-endian bytes, naming the row that was out of range.

    ⚠ `array.append` raises `OverflowError` on an out-of-range value, which is the right outcome but
    names neither the column nor the row. The range is therefore checked here, so a negative
    coordinate reads as *"gene 4,102's start_position is -1"* and not as a bare overflow a long way
    from the genome that caused it.
    """
    out = array(typecode)
    for index, value in enumerate(values):
        number = int(value)
        if not low <= number <= high:
            # ⚠ The index is into the FLATTENED words, which for an interleaved block is not the
            # gene. Naming the gene and the field is the difference between "element 14,703 is out
            # of range" and "gene 4,901 field 0", only one of which can be looked up.
            where = (
                f"element {index}"
                if fields_per_gene == 1
                else f"gene {index // fields_per_gene} field {index % fields_per_gene}"
            )
            raise GeneBlockError(
                f"{name}: {where} is {number}, outside [{low}, {high}]. "
                "A value that does not fit the column would wrap to a plausible wrong one."
            )
        out.append(number)
    if not _IS_LITTLE_ENDIAN:
        out.byteswap()
    return out.tobytes()


def decode_int32_vector(block: bytes | memoryview, *, name: str = "int32 block", expected_count: int | None = None):
    """A 4 B/gene signed block → an `array("i")`, one element per gene."""
    return _decode_vector("i", block, name=name, stride=INT32_VECTOR_STRIDE, expected_count=expected_count)


def encode_int32_vector(values: Iterable[int], *, name: str = "int32 block") -> bytes:
    """One signed int per gene → a 4 B/gene block."""
    return _encode_vector("i", values, name=name, low=INT32_MIN, high=INT32_MAX)


def positions_of(values: Sequence[int], wanted: int) -> list[int]:
    """Every index at which `values` equals `wanted` — the `flatnonzero` the Sequence tab needs.

    ⭐ This is the one query the gene→locus edge still has to answer, and it is why the assignment
    block is keyed by genome rather than by locus: ~4,900 `int32` scanned in place. A genome at
    ρ > 1 has the same locus at two indices, so this returns a list and never a first hit —
    collapsing it to one would drop the second copy of every duplicated gene, which is exactly the
    case the parity suite counts (`loci_with_copies_above_one`).
    """
    return [index for index, value in enumerate(values) if value == wanted]


# ── gene_spans — 16 B/gene, byte-identical to the `.nseq` gene table ──────────────────────────
def pack_gene_flags(*, is_minus_strand: bool, is_five_prime_partial: bool, gff_phase: int) -> int:
    """The `flags` word: bit 0 `-` strand · bit 1 5'-partial · bits 2-3 GFF phase."""
    if gff_phase not in (0, 1, 2):
        raise GeneBlockError(
            f"gff_phase must be 0, 1 or 2, not {gff_phase!r}. GFF3 defines exactly those three; "
            "a wider value would collide with a flag bit."
        )
    return (
        (FLAG_MINUS_STRAND if is_minus_strand else 0)
        | (FLAG_FIVE_PRIME_PARTIAL if is_five_prime_partial else 0)
        | (gff_phase << FLAG_PHASE_SHIFT)
    )


def encode_gene_spans(
    start_positions: Sequence[int],
    end_positions: Sequence[int],
    contig_indexes: Sequence[int],
    flags: Sequence[int],
) -> bytes:
    """Four parallel columns → `gene_spans`, 16 B/gene, four little-endian `uint32` per gene.

    ⚠ The four columns must be the same length and that is checked: a `zip` over mismatched columns
    stops at the shortest and writes a block that is well-formed, passes its CHECK constraint, and
    is missing the tail of the genome.
    """
    lengths = {
        "start_positions": len(start_positions),
        "end_positions": len(end_positions),
        "contig_indexes": len(contig_indexes),
        "flags": len(flags),
    }
    if len(set(lengths.values())) != 1:
        raise GeneBlockError(f"gene_spans: the four columns have different lengths {lengths}")

    words = array("I")
    for index in range(lengths["start_positions"]):
        row = (start_positions[index], end_positions[index], contig_indexes[index], flags[index])
        for word_index, value in enumerate(row):
            number = int(value)
            if not 0 <= number <= UINT32_MAX:
                raise GeneBlockError(
                    f"gene_spans: gene {index} word {word_index} is {number}, outside "
                    f"[0, {UINT32_MAX}]. The `.nseq` gene table is unsigned; a negative coordinate "
                    "would wrap to a coordinate near 4.29 Gb and render as a plausible gene."
                )
            words.append(number)
        if words[-4] > words[-3]:
            raise GeneBlockError(
                f"gene_spans: gene {index} has start {words[-4]} after end {words[-3]}. "
                "`.nseq` spans are 1-based inclusive and never reversed — strand is a flag bit."
            )
        if words[-1] & ~DEFINED_FLAG_BITS:
            raise GeneBlockError(
                f"gene_spans: gene {index} has flag bits {words[-1]:#x} outside the defined "
                f"{DEFINED_FLAG_BITS:#x}. Use `pack_gene_flags`."
            )
    if not _IS_LITTLE_ENDIAN:
        words.byteswap()
    return words.tobytes()


class GeneSpanTable:
    """`gene_spans` read in place — stride 4 over one `uint32` array, exactly as `app.js` holds it.

    No decode pass and no object per gene: `array("I", block)` is one C-level copy, and each
    accessor is an index into it. `words` is public because the audit census and the round-trip
    against `.nseq` both want the raw table, and hiding it would mean re-deriving it.
    """

    __slots__ = ("gene_count", "words")

    def __init__(self, block: bytes | memoryview, *, gene_count: int | None = None) -> None:
        """Wrap a `gene_spans` block, checking its length against `gene_count` when one is given."""
        self.words = _decode_vector(
            "I", block, name="gene_spans", stride=GENE_SPAN_STRIDE, expected_count=gene_count
        )
        self.gene_count = len(self.words) // 4

    def __len__(self) -> int:
        """How many genes the block holds."""
        return self.gene_count

    def _word(self, flat_index: int, offset: int) -> int:
        if not 0 <= flat_index < self.gene_count:
            raise GeneBlockError(
                f"gene_spans: flat_index {flat_index} is outside this genome's 0..{self.gene_count - 1}. "
                "⚠ A negative index would read from the end of the table and return another gene."
            )
        return self.words[4 * flat_index + offset]

    def start_position(self, flat_index: int) -> int:
        """1-based inclusive start on its contig."""
        return self._word(flat_index, 0)

    def end_position(self, flat_index: int) -> int:
        """Inclusive end, **including the stop codon** — see `.nseq`'s gene-table row."""
        return self._word(flat_index, 1)

    def contig_index(self, flat_index: int) -> int:
        """⚠ A position over the contigs that carry a CDS, NOT the contig number."""
        return self._word(flat_index, 2)

    def flags(self, flat_index: int) -> int:
        """The raw flags word."""
        return self._word(flat_index, 3)

    def length_nt(self, flat_index: int) -> int:
        """`end - start + 1`, derived rather than stored — the column the grain change removes."""
        return self.end_position(flat_index) - self.start_position(flat_index) + 1

    def is_minus_strand(self, flat_index: int) -> bool:
        """Bit 0 of `flags`."""
        return bool(self.flags(flat_index) & FLAG_MINUS_STRAND)

    def strand(self, flat_index: int) -> str:
        """`"+"` or `"-"`, the spelling the `gene.strand` column holds today."""
        return "-" if self.is_minus_strand(flat_index) else "+"

    def is_five_prime_partial(self, flat_index: int) -> bool:
        """Bit 1 of `flags`."""
        return bool(self.flags(flat_index) & FLAG_FIVE_PRIME_PARTIAL)

    def gff_phase(self, flat_index: int) -> int:
        """Bits 2-3 of `flags`."""
        return (self.flags(flat_index) >> FLAG_PHASE_SHIFT) & FLAG_PHASE_MASK

    def to_bytes(self) -> bytes:
        """The block as it would be stored — the other half of the round-trip."""
        if _IS_LITTLE_ENDIAN:
            return self.words.tobytes()
        swapped = array("I", self.words)
        swapped.byteswap()
        return swapped.tobytes()


# ── gene_name_refs — 12 B/gene of references into `gene_text_term` ────────────────────────────
def encode_gene_name_refs(
    symbol_term_ids: Sequence[int],
    product_term_ids: Sequence[int],
    locus_tag_numbers: Sequence[int],
) -> bytes:
    """Three parallel columns → `gene_name_refs`, 12 B/gene, three little-endian `int32`.

    ⚠ **All three words are signed, with `-1` meaning absent.** The design note had the third as
    `uint32`; one row holding two signednesses is a decode trap for no gain — a locus-tag number is
    a per-genome ordinal in the thousands, nowhere near `INT32_MAX`.
    """
    lengths = {
        "symbol_term_ids": len(symbol_term_ids),
        "product_term_ids": len(product_term_ids),
        "locus_tag_numbers": len(locus_tag_numbers),
    }
    if len(set(lengths.values())) != 1:
        raise GeneBlockError(f"gene_name_refs: the three columns have different lengths {lengths}")
    rows: list[int] = []
    for index in range(lengths["symbol_term_ids"]):
        rows.extend((symbol_term_ids[index], product_term_ids[index], locus_tag_numbers[index]))
    return _encode_vector(
        "i", rows, name="gene_name_refs", low=INT32_MIN, high=INT32_MAX, fields_per_gene=3
    )


class GeneNameRefTable:
    """`gene_name_refs` read in place — three `int32` per gene, `-1` for absent."""

    __slots__ = ("gene_count", "values")

    def __init__(self, block: bytes | memoryview, *, gene_count: int | None = None) -> None:
        """Wrap a `gene_name_refs` block, checking its length against `gene_count` when given."""
        self.values = _decode_vector(
            "i", block, name="gene_name_refs", stride=GENE_NAME_REF_STRIDE, expected_count=gene_count
        )
        self.gene_count = len(self.values) // 3

    def __len__(self) -> int:
        """How many genes the block holds."""
        return self.gene_count

    def _value(self, flat_index: int, offset: int) -> int:
        if not 0 <= flat_index < self.gene_count:
            raise GeneBlockError(
                f"gene_name_refs: flat_index {flat_index} is outside 0..{self.gene_count - 1}"
            )
        return self.values[3 * flat_index + offset]

    def symbol_term_id(self, flat_index: int) -> int | None:
        """The `gene_text_term` row for this gene's Bakta symbol, or `None`."""
        found = self._value(flat_index, 0)
        return None if found == ABSENT_TERM_ID else found

    def product_term_id(self, flat_index: int) -> int | None:
        """The `gene_text_term` row for this gene's Bakta product, or `None`."""
        found = self._value(flat_index, 1)
        return None if found == ABSENT_TERM_ID else found

    def locus_tag_number(self, flat_index: int) -> int | None:
        """The numeric half of the locus tag, or `None`; the prefix is one column per genome."""
        found = self._value(flat_index, 2)
        return None if found == ABSENT_LOCUS_TAG_NUMBER else found

    def term_ids_for(self, flat_indexes: Iterable[int]) -> set[int]:
        """The `gene_text_term` ids for JUST these genes — the `WHERE id = ANY(…)` a request needs.

        ⭐ This, not `referenced_term_ids`, is what the Sequence tab wants: a locus holds a handful
        of this genome's genes, so the bounded lookup is over ≤ 2 ids per displayed gene.

        ⛔ Words 0 and 1 only. Word 2 is the locus-tag NUMBER, an ordinal within this genome and not
        an id in any table; folding it in would send integers to `gene_text_term` that match
        unrelated rows and return a product for a gene that has none.
        """
        found: set[int] = set()
        for flat_index in flat_indexes:
            found.add(self._value(flat_index, 0))
            found.add(self._value(flat_index, 1))
        found.discard(ABSENT_TERM_ID)
        return found

    def referenced_term_ids(self) -> set[int]:
        """Every `gene_text_term` id in the WHOLE genome — for ingest and the column audit.

        ⚠ Not a per-request lookup: a ~4,900-gene genome references on the order of thousands of
        distinct products, so passing this to a `WHERE id = ANY(…)` would send thousands of ids to
        answer a question about a handful of genes. Use `term_ids_for` on the serving path.
        """
        found = set(self.values[0::3]) | set(self.values[1::3])
        found.discard(ABSENT_TERM_ID)
        return found

    def to_bytes(self) -> bytes:
        """The block as it would be stored."""
        if _IS_LITTLE_ENDIAN:
            return self.values.tobytes()
        swapped = array("i", self.values)
        swapped.byteswap()
        return swapped.tobytes()
