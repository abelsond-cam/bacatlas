"""The gene layer stored as fixed-stride binary blocks, one row per genome.

Why this package exists, in one paragraph: `gene`, `gene_functional_annotation` and
`gene_locus_membership` are 93 % of the serving database and all three are one row per gene, which
puts the schema at ~347 GB at 72,801 genomes against a design target of 430,000. Postgres charges
~28 B of tuple header per row regardless of columns, so the grain is the cost, not the columns. The
design of record is `nuna/docs/design/storing_the_gene_layer.md`; this package is the encoder and
decoder for the blocks that design stores.

- `packed_gene_blocks` — the codecs. **Stdlib only**: these decode on the serving path, where numpy
  is not installed.
- `nseq_gene_table` — reading nuna's published `.nseq`/`.loci`, which is what makes the claim
  *"`gene_spans` is byte-identical to the `.nseq` gene table"* a measurement rather than an
  assertion.

⛔ Nothing here imports `bac80` or `nuna.tl.cluster`. The ~14× dictionary-encoding precedent cited
in the design doc lives in `src/bac80/marker_hmm/run_extra_copies.py`, which measures the
Bacformer-2 corpora — a different layer, with a hard boundary in `nuna/CLAUDE.md`. Only the idea
travels. The byte-compatibility with `.nseq` is a FORMAT agreement, not a code dependency: this
package reads those files but imports nothing from the package that writes them.
"""

from bacatlas_backend.gene_blocks.nseq_gene_table import (
    NSEQ_MAGIC,
    NseqFile,
    read_loci,
    read_nseq,
)
from bacatlas_backend.gene_blocks.packed_gene_blocks import (
    ABSENT_LOCUS_TAG_NUMBER,
    ABSENT_TERM_ID,
    DEFINED_FLAG_BITS,
    FLAG_FIVE_PRIME_PARTIAL,
    FLAG_MINUS_STRAND,
    FLAG_PHASE_MASK,
    FLAG_PHASE_SHIFT,
    GENE_NAME_REF_STRIDE,
    GENE_SPAN_STRIDE,
    INT32_VECTOR_STRIDE,
    NO_ARRANGEMENT_RANK,
    SPAN_ENCODING_NSEQ_GENE_TABLE,
    UNASSIGNED_LOCUS_ORDINAL,
    GeneBlockError,
    GeneNameRefTable,
    GeneSpanTable,
    block_gene_count,
    byte_length,
    check_block_length,
    decode_int32_vector,
    encode_gene_name_refs,
    encode_gene_spans,
    encode_int32_vector,
    pack_gene_flags,
    positions_of,
)

__all__ = [
    "ABSENT_LOCUS_TAG_NUMBER",
    "ABSENT_TERM_ID",
    "DEFINED_FLAG_BITS",
    "FLAG_FIVE_PRIME_PARTIAL",
    "FLAG_MINUS_STRAND",
    "FLAG_PHASE_MASK",
    "FLAG_PHASE_SHIFT",
    "GENE_NAME_REF_STRIDE",
    "GENE_SPAN_STRIDE",
    "INT32_VECTOR_STRIDE",
    "NO_ARRANGEMENT_RANK",
    "NSEQ_MAGIC",
    "SPAN_ENCODING_NSEQ_GENE_TABLE",
    "UNASSIGNED_LOCUS_ORDINAL",
    "GeneBlockError",
    "GeneNameRefTable",
    "GeneSpanTable",
    "NseqFile",
    "block_gene_count",
    "byte_length",
    "check_block_length",
    "decode_int32_vector",
    "encode_gene_name_refs",
    "encode_gene_spans",
    "encode_int32_vector",
    "pack_gene_flags",
    "positions_of",
    "read_loci",
    "read_nseq",
]
