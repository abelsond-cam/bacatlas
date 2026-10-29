"""The packed gene blocks — the codecs, and the 200-genome differential oracle against `.nseq`.

⭐ **Why a differential oracle and not a round-trip.** `encode(decode(x)) == x` proves the two
halves agree with each other, which is exactly what a symmetric bug gives you for free. The claim
worth testing is that `gene_spans` is byte-identical to nuna's `.nseq` gene table and that both
decoders read the same fields out of it — so the oracle is the published `app.js`'s own `decodeNseq`,
lifted and run in node over the same 200 committed files, and the comparison is field by field on
every gene of every genome. The round-trip is here too, but as the weaker of the two.

⛔ **Coverage is asserted BEFORE any difference is reported.** `4ab35ca`: a comparison that skips
what it cannot compare reports "0 differ" while skipping exactly what changed. Every count this
suite relies on is named in `GENE_BLOCK_COVERAGE` and printed whether it passes or fails, and a
counter that is zero is stated as zero rather than left to look like agreement.

⚠ Skips, with the reason, when node, the nuna checkout (`app.js`) or the committed `.nseq` files are
absent. **No database is involved** — that is the point of this stage.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pytest

from bacatlas_backend.gene_blocks import (
    ABSENT_LOCUS_TAG_NUMBER,
    ABSENT_TERM_ID,
    DEFINED_FLAG_BITS,
    GENE_NAME_REF_STRIDE,
    GENE_SPAN_STRIDE,
    INT32_VECTOR_STRIDE,
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
    read_loci,
    read_nseq,
)
from tests.conftest import PUBLISHED_SITE_CATALOGUE_DIR, require_oracle

SEQUENCE_DIR = PUBLISHED_SITE_CATALOGUE_DIR / "seq"
DUMPER = Path(__file__).with_name("nseq_gene_table_dumper.js")
#: ⚠ Same env var and same default as `test_sequence_parity.py`, so the two suites cannot be
#: pointed at different copies of the page. It is a path into a MUTABLE working tree, not a
#: frozen artifact, which is exactly why it must be overridable.
FROZEN_APP_JS = Path(
    os.environ.get("BACATLAS_FROZEN_APP_JS", Path.home() / "developer/nuna/src/nuna/tl/locus_browser/app.js")
)


# ── the codecs, on their own ──────────────────────────────────────────────────────────────────
def test_the_flags_word_is_the_nseq_bit_layout():
    assert pack_gene_flags(is_minus_strand=False, is_five_prime_partial=False, gff_phase=0) == 0b0000
    assert pack_gene_flags(is_minus_strand=True, is_five_prime_partial=False, gff_phase=0) == 0b0001
    assert pack_gene_flags(is_minus_strand=False, is_five_prime_partial=True, gff_phase=0) == 0b0010
    assert pack_gene_flags(is_minus_strand=False, is_five_prime_partial=False, gff_phase=2) == 0b1000
    assert pack_gene_flags(is_minus_strand=True, is_five_prime_partial=True, gff_phase=1) == 0b0111
    assert DEFINED_FLAG_BITS == 0b1111


def test_a_phase_outside_gff3s_three_values_is_refused():
    """A wider phase would collide with nothing today and with the next flag bit tomorrow."""
    with pytest.raises(GeneBlockError, match="gff_phase"):
        pack_gene_flags(is_minus_strand=False, is_five_prime_partial=False, gff_phase=3)


def test_gene_spans_round_trips_every_field():
    starts, ends, contigs = [1, 500, 4_000_000], [300, 900, 4_000_123], [0, 7, 607]
    flags = [
        pack_gene_flags(is_minus_strand=True, is_five_prime_partial=False, gff_phase=0),
        pack_gene_flags(is_minus_strand=False, is_five_prime_partial=True, gff_phase=2),
        pack_gene_flags(is_minus_strand=True, is_five_prime_partial=True, gff_phase=1),
    ]
    block = encode_gene_spans(starts, ends, contigs, flags)

    assert len(block) == 3 * GENE_SPAN_STRIDE
    table = GeneSpanTable(block, gene_count=3)
    assert [table.start_position(j) for j in range(3)] == starts
    assert [table.end_position(j) for j in range(3)] == ends
    assert [table.contig_index(j) for j in range(3)] == contigs
    assert [table.strand(j) for j in range(3)] == ["-", "+", "-"]
    assert [table.is_five_prime_partial(j) for j in range(3)] == [False, True, True]
    assert [table.gff_phase(j) for j in range(3)] == [0, 2, 1]
    assert [table.length_nt(j) for j in range(3)] == [300, 401, 124]
    assert table.to_bytes() == block


def test_the_block_is_little_endian_whatever_the_cpu_is():
    """⛔ The bytes are a stored format. `start = 1` must be `01 00 00 00`, not the host's order."""
    block = encode_gene_spans([1], [2], [3], [4])
    assert block == bytes([1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0])


def test_columns_of_different_lengths_are_refused():
    """⚠ A `zip` over mismatched columns writes a well-formed block missing the genome's tail."""
    with pytest.raises(GeneBlockError, match="different lengths"):
        encode_gene_spans([1, 2], [3], [0], [0])


def test_a_negative_coordinate_is_refused_rather_than_wrapped():
    with pytest.raises(GeneBlockError, match="outside"):
        encode_gene_spans([-1], [300], [0], [0])


def test_a_reversed_span_is_refused_because_strand_is_a_flag_not_an_order():
    with pytest.raises(GeneBlockError, match="after end"):
        encode_gene_spans([900], [500], [0], [0])


def test_an_undefined_flag_bit_is_refused_on_the_way_in():
    """An unknown bit means the writer knew something this reader does not."""
    with pytest.raises(GeneBlockError, match="outside the defined"):
        encode_gene_spans([1], [300], [0], [1 << 4])


def test_a_block_length_that_is_not_a_whole_number_of_rows_is_refused():
    with pytest.raises(GeneBlockError, match="whole number"):
        block_gene_count(b"\x00" * 17, GENE_SPAN_STRIDE, name="gene_spans")


def test_a_block_that_disagrees_with_gene_count_is_refused():
    """⛔ One row short reads every later gene as the one before it, with no error anywhere."""
    with pytest.raises(GeneBlockError, match="gene_count says"):
        check_block_length(b"\x00" * (2 * GENE_SPAN_STRIDE), GENE_SPAN_STRIDE, 3, name="gene_spans")


def test_an_out_of_range_flat_index_is_refused_rather_than_read_from_the_end():
    table = GeneSpanTable(encode_gene_spans([1], [300], [0], [0]))
    with pytest.raises(GeneBlockError, match="outside"):
        table.start_position(-1)
    with pytest.raises(GeneBlockError, match="outside"):
        table.start_position(1)


def test_gene_name_refs_round_trip_with_absence_distinguished_from_zero():
    """⚠ `0` is a legitimate locus-tag number and term id 0 would be a real row; `-1` is absence."""
    block = encode_gene_name_refs([5, ABSENT_TERM_ID, 0], [7, 7, ABSENT_TERM_ID], [1, 2, ABSENT_LOCUS_TAG_NUMBER])

    assert len(block) == 3 * GENE_NAME_REF_STRIDE
    table = GeneNameRefTable(block, gene_count=3)
    assert [table.symbol_term_id(j) for j in range(3)] == [5, None, 0]
    assert [table.product_term_id(j) for j in range(3)] == [7, 7, None]
    assert [table.locus_tag_number(j) for j in range(3)] == [1, 2, None]
    assert table.to_bytes() == block


def test_referenced_term_ids_excludes_the_locus_tag_number():
    """⛔ Word 2 is an ordinal in this genome, not an id — it must not reach `gene_text_term`."""
    table = GeneNameRefTable(encode_gene_name_refs([5], [7], [999_999]))
    assert table.referenced_term_ids() == {5, 7}






def test_int32_vectors_round_trip_and_keep_the_unassigned_sentinel():
    ordinals = [163, UNASSIGNED_LOCUS_ORDINAL, 889, 163]
    block = encode_int32_vector(ordinals, name="gene_locus_ordinals")

    assert len(block) == len(ordinals) * INT32_VECTOR_STRIDE
    decoded = decode_int32_vector(block, name="gene_locus_ordinals", expected_count=4)
    assert list(decoded) == ordinals


def test_positions_of_returns_every_copy_not_the_first():
    """⭐ A genome at ρ > 1 has the same locus at two indices; a first-hit scan drops the second."""
    ordinals = decode_int32_vector(encode_int32_vector([163, -1, 889, 163, 163]))
    assert positions_of(ordinals, 163) == [0, 3, 4]
    assert positions_of(ordinals, 889) == [2]
    assert positions_of(ordinals, 42) == []


def test_a_buffers_length_is_measured_in_bytes_not_in_items():
    """⛔ `len()` on a `memoryview` counts ITEMS. The miscount goes BOTH ways, which is the danger.

    `memoryview(array("I", block))` has `itemsize` 4, so `len()` is a quarter of the byte length.
    Measured: a 4,900-gene block so wrapped counted as 1,225 genes, and — the quiet direction — a
    19,600-gene block was ACCEPTED as 4,900 and then read as 19,600, which is the "whole genome of
    plausible wrong coordinates with no error anywhere" arriving through the check itself.
    """
    from array import array

    gene_count = 400
    block = encode_gene_spans(
        list(range(1, gene_count + 1)), list(range(2, gene_count + 2)), [0] * gene_count, [0] * gene_count
    )
    assert len(block) == gene_count * GENE_SPAN_STRIDE

    view = memoryview(array("I", block))
    assert len(view) == gene_count * 4, "the view reports WORDS — this is the trap"
    assert byte_length(view) == gene_count * GENE_SPAN_STRIDE
    assert block_gene_count(view, GENE_SPAN_STRIDE, name="gene_spans") == gene_count
    assert check_block_length(view, GENE_SPAN_STRIDE, gene_count, name="gene_spans") == gene_count
    # And the decode agrees with the check, which is the invariant that was broken.
    assert GeneSpanTable(view, gene_count=gene_count).gene_count == gene_count


def test_term_ids_for_is_bounded_by_the_genes_asked_about():
    """⭐ The serving path looks up a handful of genes, not the whole genome's vocabulary."""
    table = GeneNameRefTable(
        encode_gene_name_refs([5, 11, ABSENT_TERM_ID, 5], [7, 12, 13, ABSENT_TERM_ID], [1, 2, 3, 4])
    )
    assert table.term_ids_for([0]) == {5, 7}
    assert table.term_ids_for([0, 3]) == {5, 7}
    assert table.term_ids_for([2]) == {13}
    assert table.term_ids_for([]) == set()
    # The whole-genome view is a different question and stays available for ingest and the audit.
    assert table.referenced_term_ids() == {5, 7, 11, 12, 13}


def test_a_file_that_names_another_sample_is_refused_when_a_name_is_expected():
    """⭐ What makes the header's `sample` key do its job: *"a mis-served file names itself"*.

    ⚠ Unasked, the key is read and compared to nothing, so a file published under the wrong name
    decodes perfectly and serves another genome's coordinates.
    """
    require_oracle(
        SEQUENCE_DIR.is_dir(),
        f"the committed .nseq files are not at {SEQUENCE_DIR}",
        supplies="whether `expect_sample=` actually refuses a file published under another "
                 "genome's name — unasked, such a file decodes perfectly and serves another "
                 "genome's coordinates",
    )
    one = sorted(SEQUENCE_DIR.glob("*.nseq"))[0]

    assert read_nseq(one, expect_sample=one.stem).sample == one.stem
    with pytest.raises(GeneBlockError, match="was asked for"):
        read_nseq(one, expect_sample="SAMEA_SOMETHING_ELSE")


# ── the 200-genome differential oracle ────────────────────────────────────────────────────────
#: ⛔ Every counter this comparison relies on. Printed whether it passes or fails, so a counter
#: that is zero is STATED rather than mistaken for agreement.
GENE_BLOCK_COVERAGE = (
    "genomes",
    "genes",
    "genes_round_tripped_byte_identical",
    "samples_naming_their_own_filename",
    "minus_strand_genes",
    "plus_strand_genes",
    "five_prime_partial_genes",
    "genes_with_phase_above_zero",
    "genes_with_a_flag_bit_outside_the_defined_four",
    "contigs",
    "ambiguous_runs",
    "loci_vectors",
    "loci_entries",
    "loci_unassigned_entries",
)


def _render_genes(table: GeneSpanTable) -> str:
    """The canonical rendering — **identical to `nseq_gene_table_dumper.js`'s `renderGenes`**."""
    return "\n".join(
        f"{table.start_position(j)},{table.end_position(j)},{table.contig_index(j)},"
        f"{table.strand(j)},{1 if table.is_five_prime_partial(j) else 0},{table.gff_phase(j)}"
        for j in range(len(table))
    )


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@pytest.fixture(scope="module")
def page_side() -> dict:
    """What the published `app.js`'s own `decodeNseq` makes of all 200 committed genomes."""
    # ⛔ FAIL, not skip. This is the codec's ONLY independent oracle — our decoder against the
    # published `app.js`'s own `decodeNseq`, on all 200 committed genomes. A missing `node` used to
    # report green, which is the exact shape of defect the comparison exists to catch: an opaque
    # block can be wrong in five ways that all satisfy the length CHECK, and nothing else in the
    # suite would see any of them.
    supplies = (
        "the gene-block codec's only independent oracle — our decoder against the published "
        "app.js's own decodeNseq over all 200 committed genomes"
    )
    node = shutil.which("node")
    require_oracle(bool(node), "node is not on PATH, so the page's own decoder cannot be run",
                   supplies=supplies)  # fmt: skip
    require_oracle(FROZEN_APP_JS.is_file(), f"the nuna checkout's app.js is not at {FROZEN_APP_JS}",
                   supplies=supplies)  # fmt: skip
    require_oracle(SEQUENCE_DIR.is_dir(), f"the committed .nseq files are not at {SEQUENCE_DIR}",
                   supplies=supplies)  # fmt: skip

    genomes = sorted(path.stem for path in SEQUENCE_DIR.glob("*.nseq"))
    require_oracle(bool(genomes), f"no .nseq files in {SEQUENCE_DIR}", supplies=supplies)
    request = {
        "app_js": str(FROZEN_APP_JS),
        "genomes": [
            {
                "sample": sample,
                "nseq": str(SEQUENCE_DIR / f"{sample}.nseq"),
                "loci": str(SEQUENCE_DIR / f"{sample}.loci") if (SEQUENCE_DIR / f"{sample}.loci").is_file() else None,
            }
            for sample in genomes
        ],
    }
    finished = subprocess.run(  # noqa: S603
        [node, str(DUMPER), "-"],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        check=False,
    )
    if finished.returncode != 0:
        pytest.fail(f"the page-side dumper failed: {finished.stderr[-2000:]}")
    return {row["sample"]: row for row in json.loads(finished.stdout)["genomes"]}


@pytest.fixture(scope="module")
def comparison(page_side) -> dict:
    """Decode every committed genome with `gene_blocks` and compare, counting as it goes."""
    coverage: Counter = Counter()
    differences: Counter = Counter()
    examples: list[str] = []

    for sample, page in sorted(page_side.items()):
        found = read_nseq(SEQUENCE_DIR / f"{sample}.nseq")
        coverage["genomes"] += 1
        coverage["contigs"] += len(found.contig_names)
        coverage["ambiguous_runs"] += len(found.ambiguous_runs)

        # ⛔ The filename, which is the ONE thing both decoders cannot agree about by reading the
        # same bytes. Everything else here is two readings of one JSON string, so a file published
        # under the wrong name agrees with itself on every field — `header.sample` against
        # `page["header_sample"]` is a tautology, and only this comparison can fail.
        if found.sample == sample:
            coverage["samples_naming_their_own_filename"] += 1
        else:
            differences["header.sample_is_not_the_filename"] += 1
            if len(examples) < 10:
                examples.append(f"{sample}.nseq: the header names {found.sample!r}")

        # ⛔ Whole tables, not a head and a count. Reversing `contig_lengths` moved every contig
        # base offset and the suite still reported "no differences" — measured 2026-10-26.
        for name, ours, theirs in (
            ("header.sample", found.sample, page["header_sample"]),
            ("header.n_genes", found.gene_count, page["n_genes"]),
            ("header.n_bases", found.base_count, page["n_bases"]),
            ("header.contig_count", len(found.contig_names), page["contig_count"]),
            ("header.contig_names", _digest("\n".join(found.contig_names)), page["contig_names_digest"]),
            (
                "header.contig_lengths",
                _digest("\n".join(str(length) for length in found.contig_lengths)),
                page["contig_lengths_digest"],
            ),
            (
                "header.contig_offsets",
                _digest("\n".join(str(offset) for offset in found.contig_base_offsets())),
                page["contig_offsets_digest"],
            ),
            (
                "header.ambiguous_runs",
                _digest("\n".join(f"{start},{length}" for start, length in found.ambiguous_runs)),
                page["ambiguous_digest"],
            ),
        ):
            if ours != theirs:
                differences[name] += 1
                if len(examples) < 10:
                    examples.append(f"{sample} {name}: ours {ours!r} vs the page's {theirs!r}")

        table = GeneSpanTable(found.gene_table, gene_count=found.gene_count)
        coverage["genes"] += len(table)

        # ⭐ The round-trip: decoded fields re-encoded must be the published bytes, byte for byte.
        # ⛔ The failure is REPORTED, not raised. `encode_gene_spans` refuses an undefined flag bit
        # and a reversed span, so letting it propagate meant a committed file carrying either would
        # crash this fixture — and the assertion downstream that no such flag bit exists could then
        # never fail, because the run never reached it. The sentinel has to be downstream of the
        # failure, not upstream of it.
        try:
            reencoded = encode_gene_spans(
                [table.start_position(j) for j in range(len(table))],
                [table.end_position(j) for j in range(len(table))],
                [table.contig_index(j) for j in range(len(table))],
                [table.flags(j) for j in range(len(table))],
            )
        except GeneBlockError as refusal:
            reencoded = None
            differences["gene_spans_refused_by_the_encoder"] += 1
            if len(examples) < 10:
                examples.append(f"{sample}: the published .nseq will not re-encode — {refusal}")
        if reencoded is None:
            pass
        elif reencoded == found.gene_table:
            coverage["genes_round_tripped_byte_identical"] += len(table)
        else:
            differences["gene_spans_bytes"] += 1
            if len(examples) < 10:
                examples.append(f"{sample}: re-encoded gene table differs from the published .nseq")

        # ⛔ The differential oracle: every field of every gene, against the page's own decoder.
        if _digest(_render_genes(table)) != page["genes_digest"]:
            differences["gene_fields_digest"] += 1
            if len(examples) < 10:
                ours_head = _render_genes(table).split("\n")[:3]
                examples.append(f"{sample}: ours {ours_head} vs the page's {page['genes_head']}")

        counters = page["counters"]
        minus = sum(1 for j in range(len(table)) if table.is_minus_strand(j))
        partial = sum(1 for j in range(len(table)) if table.is_five_prime_partial(j))
        phased = sum(1 for j in range(len(table)) if table.gff_phase(j) > 0)
        outside = sum(1 for j in range(len(table)) if table.flags(j) & ~DEFINED_FLAG_BITS)
        coverage["minus_strand_genes"] += minus
        coverage["plus_strand_genes"] += len(table) - minus
        coverage["five_prime_partial_genes"] += partial
        coverage["genes_with_phase_above_zero"] += phased
        coverage["genes_with_a_flag_bit_outside_the_defined_four"] += outside
        for name, ours, theirs in (
            ("counter.minus_strand", minus, counters["minus_strand"]),
            ("counter.plus_strand", len(table) - minus, counters["plus_strand"]),
            ("counter.five_prime_partial", partial, counters["five_prime_partial"]),
            ("counter.phase_above_zero", phased, sum(counters["phase"][1:])),
            ("counter.flag_bits_outside", outside, counters["flag_bits_outside_the_defined_four"]),
        ):
            if ours != theirs:
                differences[name] += 1
                if len(examples) < 10:
                    examples.append(f"{sample} {name}: ours {ours} vs the page's {theirs}")

        if page.get("loci_digest") is not None:
            ordinals = read_loci(SEQUENCE_DIR / f"{sample}.loci", gene_count=found.gene_count)
            coverage["loci_vectors"] += 1
            coverage["loci_entries"] += len(ordinals)
            unassigned = sum(1 for value in ordinals if value == UNASSIGNED_LOCUS_ORDINAL)
            coverage["loci_unassigned_entries"] += unassigned
            if _digest("\n".join(str(value) for value in ordinals)) != page["loci_digest"]:
                differences["loci_digest"] += 1
            if unassigned != page["loci_unassigned"]:
                differences["loci_unassigned"] += 1
            if encode_int32_vector(ordinals) != (SEQUENCE_DIR / f"{sample}.loci").read_bytes():
                differences["loci_bytes"] += 1

    summary = (
        "gene_blocks vs the page's own decodeNseq over the committed .nseq: "
        + " · ".join(f"{name} {coverage[name]:,}" for name in GENE_BLOCK_COVERAGE)
        + (f" · DIFFERENCES: {dict(differences)} {examples}" if differences else " · no differences")
    )
    print("\n" + summary)
    return {"coverage": coverage, "differences": differences, "summary": summary, "examples": examples}


def _assert_coverage(comparison: dict) -> None:
    """⛔ Called FIRST by every test below: what was compared, before anything about what differed."""
    coverage = comparison["coverage"]
    assert coverage["genomes"] == 200, comparison["summary"]
    assert coverage["genes"] > 900_000, comparison["summary"]
    assert coverage["genes_round_tripped_byte_identical"] == coverage["genes"], comparison["summary"]
    assert coverage["samples_naming_their_own_filename"] == coverage["genomes"], comparison["summary"]
    assert coverage["loci_vectors"] == coverage["genomes"], comparison["summary"]
    assert coverage["loci_entries"] == coverage["genes"], comparison["summary"]
    assert coverage["contigs"] > 10_000, comparison["summary"]
    # Both strands, or the flag bit is untested on the half of the cohort where a flank goes wrong.
    assert coverage["minus_strand_genes"] > 0.4 * coverage["genes"], comparison["summary"]
    assert coverage["plus_strand_genes"] > 0.4 * coverage["genes"], comparison["summary"]


def test_every_committed_genome_is_decoded_and_counted_before_anything_is_reported(comparison):
    _assert_coverage(comparison)


def test_the_header_we_read_is_the_header_the_page_reads(comparison):
    _assert_coverage(comparison)
    named = {name: count for name, count in comparison["differences"].items() if name.startswith("header.")}
    assert not named, comparison["summary"]


def test_every_genes_decoded_fields_match_the_pages_own_decoder(comparison):
    """⭐ The claim under test: `gene_spans` and `.nseq`'s gene table are the same bytes, read the same."""
    _assert_coverage(comparison)
    assert not comparison["differences"]["gene_fields_digest"], comparison["summary"]
    named = {name: count for name, count in comparison["differences"].items() if name.startswith("counter.")}
    assert not named, comparison["summary"]


def test_re_encoding_what_we_decoded_reproduces_the_published_bytes(comparison):
    _assert_coverage(comparison)
    assert not comparison["differences"]["gene_spans_bytes"], comparison["summary"]
    assert not comparison["differences"]["gene_spans_refused_by_the_encoder"], comparison["summary"]


def test_the_loci_vector_is_the_same_layout_as_gene_locus_ordinals(comparison):
    """`.loci` is `n_genes` little-endian `int32` — the column's layout, proven on a published file."""
    _assert_coverage(comparison)
    for name in ("loci_digest", "loci_bytes", "loci_unassigned"):
        assert not comparison["differences"][name], comparison["summary"]


def test_the_thin_coverage_is_named_rather_than_left_to_look_like_agreement(comparison):
    """⚠ Two flag bits are exercised by no committed genome. Named here so stage 4 knows.

    `gff_phase` has one distinct value (0) over all 1.44 M genes in `bacatlas_ingest`, which is the
    design doc's risk 4. If either counter ever becomes non-zero this test is the thing that says
    the cohort changed — and `pack_gene_flags` is unit-tested on both bits regardless.
    """
    _assert_coverage(comparison)
    coverage = comparison["coverage"]
    assert coverage["genes_with_a_flag_bit_outside_the_defined_four"] == 0, comparison["summary"]
    assert coverage["genes_with_phase_above_zero"] == 0, (
        "a committed genome now carries a non-zero GFF phase — the derivation of protein_length_aa "
        f"that assumed phase 0 must be re-checked. {comparison['summary']}"
    )
