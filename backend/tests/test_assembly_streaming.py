"""`ContigSequenceProvider` — the held → cached → ATB chain, and the identity it rests on.

⛔ **The claim this file has to defend.** Streaming is safe only because the ATB assembly and the
``##FASTA`` block of the Bakta GFF we hold are *the same bytes*. BakRep ran Bakta over the ATB
assemblies, so they should be — but "should be" is how a page comes to show plausible wrong bases at
correct-looking coordinates, which nothing downstream can detect. So the identity is MEASURED:
`test_streamed_bases_are_IDENTICAL_to_the_held_ones` fetches real objects and compares them against
the real GFFs, and is the test that would have to fail for this design to be wrong.

⚠ **That test needs the network, so it is gated** behind `BACATLAS_ATB_NETWORK_TESTS=1` and skips
loudly otherwise. Everything else here runs offline with a stub store. Measured 2026-10-23 over all
280 probe genomes outside the suite: **280/280 identical, 0 fetch failures**, mean fetch 0.241 s.
"""

from __future__ import annotations

import gzip
import io
import os
from pathlib import Path

import pytest

from bacatlas_backend.configuration import Configuration
from bacatlas_backend.gff.fasta_reader import read_fasta_block
from bacatlas_backend.gff.gff_cds_parser import parse_genome_annotation
from bacatlas_backend.sequence_source.atb_assembly_store import (
    AssemblyUnavailable,
    AtbAssemblyStore,
    AtbStoreSettings,
)
from bacatlas_backend.sequence_source.contig_sequence_provider import (
    ContigSequenceProvider,
    GenomeSequenceRequest,
    SequenceSourceUnavailable,
)
from bacatlas_backend.sequence_source.provider_factory import build_sequence_provider

HELD_SEQID = "SAMN04393397.contig00001"
HELD_BASES = "ACGTACGTACGTACGTACGTACGT"
STREAMED_BASES = "TTTTCCCCGGGGAAAATTTTCCCC"


class _StubStore:
    """An `AtbAssemblyStore` as far as the provider is concerned, and it counts its calls."""

    def __init__(self, sequences: dict[str, str] | None = None, error: Exception | None = None) -> None:
        self._sequences = sequences
        self._error = error
        self.calls: list[str] = []

    def contig_sequences(self, sample_id: str) -> dict[str, str]:
        self.calls.append(sample_id)
        if self._error is not None:
            raise self._error
        assert self._sequences is not None
        return self._sequences


def _write_held_gff(root: Path, dataset_id: str, sample_id: str, *, seqid: str, bases: str) -> Path:
    """A minimal but REAL Bakta-shaped GFF: one CDS line, then a `##FASTA` block."""
    directory = root / dataset_id / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{sample_id}.bakta.gff3.gz"
    body = (
        "##gff-version 3\n"
        f"{seqid}\tBakta\tCDS\t1\t{len(bases)}\t.\t+\t0\tID=x;locus_tag=TEST_00001\n"
        "##FASTA\n"
        f">{seqid}\n{bases}\n"
    )
    with gzip.open(path, "wt") as handle:
        handle.write(body)
    return path


# ── the order, which is what makes this safe to switch on ──────────────────────────────────────
def test_a_HELD_genome_is_read_from_disk_and_ATB_is_never_asked(tmp_path):
    """⛔⛔ **The property that makes enabling streaming additive rather than a change.** A genome
    whose GFF is on this server must behave exactly as before, so switching streaming on can only
    affect genomes that previously returned `sequence_unavailable`."""
    _write_held_gff(tmp_path, "N043", "SAMN04393397", seqid=HELD_SEQID, bases=HELD_BASES)
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)

    sequences = provider.contig_sequences(
        GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id="N043")
    )

    assert sequences[HELD_SEQID] == HELD_BASES, "the held copy must win"
    assert store.calls == [], "ATB was asked about a genome we hold"


def test_a_genome_we_do_NOT_hold_is_streamed(tmp_path):
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)

    sequences = provider.contig_sequences(
        GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id="N043")
    )

    assert sequences[HELD_SEQID] == STREAMED_BASES
    assert store.calls == ["SAMN04393397"]


def test_a_genome_with_NO_dataset_id_is_streamed_rather_than_refused(tmp_path):
    """⚠ `bakrep_dataset_id` names the directory a held GFF sits in, and its absence used to be
    fatal by definition. For a streamed genome there is no such directory — so absence now means
    "we hold no copy", the ordinary case, and must not read as a fault."""
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)

    sequences = provider.contig_sequences(
        GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None)
    )

    assert sequences[HELD_SEQID] == STREAMED_BASES


# ── the failures, each naming which of the three ways it happened ───────────────────────────────
def test_with_streaming_OFF_a_genome_we_do_not_hold_names_the_setting(tmp_path):
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None)

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(
            GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id="N043")
        )
    assert "not on this server" in str(raised.value)
    assert raised.value.is_retryable is False


def test_with_NOTHING_configured_the_failure_names_BOTH_settings():
    """⚠ "Not configured" and "the file is missing" are indistinguishable unless the message says
    so, and only one of them is ours to fix. This is the case the endpoint used to special-case."""
    provider = ContigSequenceProvider(gff_root=None, atb_store=None)

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(
            GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None)
        )
    detail = str(raised.value)
    assert "BACATLAS_ROOT_GFF" in detail
    assert "BACATLAS_ATB_STREAMING" in detail


def test_the_streams_RETRYABILITY_survives_the_hand_off(tmp_path):
    """⛔ The `Retry-After` header depends on this flag crossing two exception types without being
    flattened — store → provider → service → endpoint."""
    store = _StubStore(error=AssemblyUnavailable("ATB was unreachable", is_retryable=True))
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(
            GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None)
        )
    assert raised.value.is_retryable is True


def test_a_held_gff_with_no_FASTA_block_says_so(tmp_path):
    """⚠ All 280 probe GFFs carry one (measured 2026-09-04); one that does not must fail with a
    sentence rather than a KeyError further down."""
    directory = tmp_path / "N043" / "SAMN04393397"
    directory.mkdir(parents=True)
    with gzip.open(directory / "SAMN04393397.bakta.gff3.gz", "wt") as handle:
        handle.write("##gff-version 3\n")

    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None)
    with pytest.raises(SequenceSourceUnavailable, match="no ##FASTA block"):
        provider.contig_sequences(
            GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id="N043")
        )


# ── the parsed cache ───────────────────────────────────────────────────────────────────────────
def test_a_second_locus_on_the_same_genome_does_not_re_fetch(tmp_path):
    """⚠ A reader walks loci within one genome, so this is the common case, not an optimisation.
    Without it every locus on a streamed genome would cost another 0.24 s."""
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)
    request = GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None)

    for _ in range(5):
        provider.contig_sequences(request)

    assert store.calls == ["SAMN04393397"]


def test_the_parsed_cache_EVICTS_and_keeps_the_hot_entry(tmp_path):
    """⛔ A bacterial genome is ~5.3 MB of Python strings, so the bound is a memory bound. And the
    eviction must be LRU: the first version of this used "append on every touch and dedupe", which
    silently evicted the entry being read on every call."""
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store, cache_size=2)

    hot = GenomeSequenceRequest(sample_id="HOT", bakrep_dataset_id=None)
    provider.contig_sequences(hot)
    for name in ("A", "B", "C"):
        provider.contig_sequences(GenomeSequenceRequest(sample_id=name, bakrep_dataset_id=None))
        provider.contig_sequences(hot)  # keep it hot

    assert store.calls.count("HOT") == 1, f"the hot entry was evicted: {store.calls}"


def test_clear_drops_the_cache(tmp_path):
    store = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=store)
    request = GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None)

    provider.contig_sequences(request)
    provider.clear()
    provider.contig_sequences(request)

    assert store.calls == ["SAMN04393397", "SAMN04393397"]


# ── configuration ──────────────────────────────────────────────────────────────────────────────
def test_streaming_is_OFF_unless_configured():
    """⛔ An existing deployment must behave exactly as it did until someone switches this on."""
    provider = build_sequence_provider(Configuration(database_url="postgresql://x/y"))
    assert provider.can_stream is False


def test_streaming_is_ON_when_configured(tmp_path):
    provider = build_sequence_provider(
        Configuration(
            database_url="postgresql://x/y",
            atb_streaming_enabled=True,
            atb_cache_directory=tmp_path,
        )
    )
    assert provider.can_stream is True


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
def test_the_streaming_flag_accepts_the_same_words_every_other_flag_does(monkeypatch, value):
    """⚠ `BACATLAS_ATB_STREAMING=on` silently meaning *off* is exactly the configuration bug that
    reads as "streaming does not work"."""
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv("BACATLAS_ATB_STREAMING", value)
    assert Configuration.from_environment().atb_streaming_enabled is True


@pytest.mark.parametrize("value", ["", "0", "false", "no", "off", "maybe"])
def test_anything_else_leaves_streaming_off(monkeypatch, value):
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv("BACATLAS_ATB_STREAMING", value)
    assert Configuration.from_environment().atb_streaming_enabled is False


def test_an_unparseable_tuning_value_RAISES_rather_than_falling_back(monkeypatch):
    """⚠ A typo'd knob must not look like one that was never set — otherwise a deployment runs on a
    default while its operator believes it is running on the value they wrote."""
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv("BACATLAS_ATB_MAX_CONCURRENT_FETCHES", "four")
    with pytest.raises(RuntimeError, match="not an integer"):
        Configuration.from_environment()


# ── ⭐ the measured claim, against the real bucket ──────────────────────────────────────────────
NETWORK_TESTS_ENABLED = os.environ.get("BACATLAS_ATB_NETWORK_TESTS", "").strip().lower() in {
    "1",
    "true",
    "yes",
    "on",
}


@pytest.mark.skipif(
    not NETWORK_TESTS_ENABLED,
    reason="BACATLAS_ATB_NETWORK_TESTS is not set — this test fetches real objects from "
    "AllTheBacteria over the network, which a CI run should not do by default",
)
def test_streamed_bases_are_IDENTICAL_to_the_held_ones():
    """⭐⭐ **The whole design rests on this.** For each held GFF, fetch the same genome from ATB and
    compare the contig maps the two produce — not approximately, not by length: equal maps.

    ⛔ Read through the SAME `read_fasta_block` both sides, which is why that rule is one function:
    the defline differs (ATB adds `len=`, `cov=`, `origname=`, `sw=`), the wrapping differs (60 vs
    70 columns) and the case may differ, so two readers would disagree for reasons that have nothing
    to do with the assembly.
    """
    gff_root = os.environ.get("BACATLAS_ROOT_GFF")
    if not gff_root:
        pytest.skip("BACATLAS_ROOT_GFF is not set — there is nothing to compare the stream against")
    held = sorted(Path(gff_root).glob("*/*/*.bakta.gff3.gz"))
    if not held:
        pytest.skip(f"no GFFs under BACATLAS_ROOT_GFF={gff_root}")

    limit = int(os.environ.get("BACATLAS_ATB_NETWORK_TEST_GENOMES", "5"))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=None))

    compared = 0
    for path in held[:limit]:
        sample_id = path.name.removesuffix(".bakta.gff3.gz")
        streamed = store.contig_sequences(sample_id)
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
            else:
                pytest.fail(f"{path.name} carries no ##FASTA block")
            from_gff = read_fasta_block(handle)

        assert set(streamed) == set(from_gff), f"{sample_id}: seqid sets differ"
        assert streamed == from_gff, f"{sample_id}: bases differ"
        compared += 1

    assert compared > 0
    # ⚠ Asserted so a glob that silently matched nothing cannot read as a pass — the repo's own
    # "silence is not success" rule, and the trap `##FASTA`-less files would otherwise spring.
    assert compared == min(limit, len(held))


@pytest.mark.skipif(not NETWORK_TESTS_ENABLED, reason="BACATLAS_ATB_NETWORK_TESTS is not set")
def test_the_gzipped_assembly_round_trips_through_the_real_store(tmp_path):
    """⚠ The real object against the real decompressor, which the synthetic fixtures cannot prove:
    ATB's gzip is produced by a different tool than `gzip.open` writes."""
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=tmp_path))
    sequences = store.contig_sequences("SAMN04393397")

    assert len(sequences) > 1, "a draft short-read assembly has many contigs"
    assert all(set(bases) <= set("ACGTN") for bases in sequences.values())
    assert (tmp_path / "SAMN04393397.fa.gz").exists()
    # The second read comes from the cache and must agree exactly.
    assert AtbAssemblyStore(AtbStoreSettings(cache_directory=tmp_path)).contig_sequences(
        "SAMN04393397"
    ) == sequences


def test_the_held_reader_and_the_fasta_reader_agree_on_a_synthetic_file(tmp_path):
    """⚠ The offline half of the identity claim: one file read both ways must give one answer, so a
    divergence in the two code paths fails here rather than only over the network."""
    path = _write_held_gff(tmp_path, "N043", "SAMN04393397", seqid=HELD_SEQID, bases=HELD_BASES)

    through_gff = parse_genome_annotation(path, want_sequence=True).contig_sequences
    plain = _gzipped_fasta_as_read(f">{HELD_SEQID} len=24 cov=1\n{HELD_BASES}\n")

    assert through_gff == plain


def _gzipped_fasta_as_read(text: str) -> dict[str, str]:
    buffer = io.BytesIO()
    with gzip.open(buffer, "wt") as handle:
        handle.write(text)
    with gzip.open(io.BytesIO(buffer.getvalue()), "rt") as handle:
        return read_fasta_block(handle)
