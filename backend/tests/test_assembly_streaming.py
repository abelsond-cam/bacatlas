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
from bacatlas_backend.models.enumerations import SampleIdentifierKind
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


# ── ⭐ the fourth station: NCBI for RefSeq/GenBank-keyed genomes, and NEVER for BioSample ones ─────
STEM = "GCF_000512165.1_ASM51216v1_genomic"
BARE = "GCF_000512165.1"
#: The namespace split as the live API returns it: `NZ_`-prefixed headers against GenBank-style seqids.
NCBI_HEADERS = {"NZ_CP003999.1": "ACGTACGT", "NZ_CP003996.1": "TTTTCCCC"}


def _ncbi_request(**overrides) -> GenomeSequenceRequest:
    fields = dict(
        sample_id=STEM,
        bakrep_dataset_id=None,
        sample_id_kind=SampleIdentifierKind.ASSEMBLY_STEM,
        bare_accession=BARE,
        wanted_seqids=frozenset({"CP003999.1", "CP003996.1"}),
    )
    fields.update(overrides)
    return GenomeSequenceRequest(**fields)


def test_an_NCBI_keyed_genome_is_fetched_from_NCBI_by_its_BARE_accession_and_ATB_is_never_asked(tmp_path):
    """⛔ The bucket holds no stem under any name (0 of 5,389 found, 2026-10-29), so asking it costs a
    request and a 404 the negative cache would remember. The kind decides, and only NCBI is asked."""
    atb = _StubStore({HELD_SEQID: STREAMED_BASES})
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=atb, ncbi_store=ncbi)

    sequences = provider.contig_sequences(_ncbi_request())

    assert ncbi.calls == [BARE], "NCBI must be asked for the bare accession, not the stem"
    assert atb.calls == []
    assert sequences == {"CP003999.1": "ACGTACGT", "CP003996.1": "TTTTCCCC"}


def test_a_BioSample_genome_is_fetched_from_ATB_and_NCBI_is_never_asked(tmp_path):
    atb = _StubStore({HELD_SEQID: STREAMED_BASES})
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=atb, ncbi_store=ncbi)

    sequences = provider.contig_sequences(
        GenomeSequenceRequest(
            sample_id="SAMN04393397",
            bakrep_dataset_id=None,
            sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
            wanted_seqids=frozenset({HELD_SEQID}),
        )
    )

    assert atb.calls == ["SAMN04393397"]
    assert ncbi.calls == []
    assert sequences[HELD_SEQID] == STREAMED_BASES


def test_a_request_that_does_not_say_its_kind_behaves_as_before_NCBI_existed(tmp_path):
    """⚠ Every caller built before this station leaves the kind unset; they must keep going to ATB."""
    atb = _StubStore({HELD_SEQID: STREAMED_BASES})
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=atb, ncbi_store=ncbi)

    provider.contig_sequences(GenomeSequenceRequest(sample_id="SAMN04393397", bakrep_dataset_id=None))

    assert atb.calls == ["SAMN04393397"]
    assert ncbi.calls == []


def test_an_NCBI_keyed_genome_with_NCBI_OFF_names_the_setting_and_does_NOT_fall_through_to_ATB(tmp_path):
    atb = _StubStore({HELD_SEQID: STREAMED_BASES})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=atb, ncbi_store=None)

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(_ncbi_request())
    assert "BACATLAS_NCBI_STREAMING" in str(raised.value)
    assert raised.value.is_retryable is False
    assert atb.calls == [], "a stem was sent to AllTheBacteria"


def test_an_NCBI_keyed_genome_with_NO_bare_accession_is_refused_by_name(tmp_path):
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=ncbi)

    with pytest.raises(SequenceSourceUnavailable, match="no bare accession") as raised:
        provider.contig_sequences(_ncbi_request(bare_accession=None))
    assert raised.value.is_retryable is False
    assert ncbi.calls == []


def test_a_HELD_copy_wins_for_an_NCBI_keyed_genome_too(tmp_path):
    _write_held_gff(tmp_path, "N043", STEM, seqid="CP003999.1", bases=HELD_BASES)
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=ncbi)

    sequences = provider.contig_sequences(
        _ncbi_request(bakrep_dataset_id="N043", wanted_seqids=frozenset({"CP003999.1"}))
    )

    assert sequences == {"CP003999.1": HELD_BASES}
    assert ncbi.calls == []


# ── the namespace: an NCBI assembly answers as `NZ_…`, the coordinates say `CP…` ───────────────────
def test_an_NCBI_assembly_is_RE_KEYED_to_the_contig_names_the_coordinates_use(tmp_path):
    """⭐ The whole reason `wanted_seqids` is on the request: the service's `sequences.get(row.seqid)`
    must hit for `CP003999.1` although NCBI wrote `NZ_CP003999.1`."""
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=_StubStore(NCBI_HEADERS))

    sequences = provider.contig_sequences(_ncbi_request(wanted_seqids=frozenset({"CP003999.1"})))

    assert sequences == {"CP003999.1": "ACGTACGT"}


def test_an_NCBI_assembly_whose_names_match_EXACTLY_is_returned_as_is(tmp_path):
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=_StubStore(NCBI_HEADERS))

    sequences = provider.contig_sequences(_ncbi_request(wanted_seqids=frozenset({"NZ_CP003996.1"})))

    assert sequences == {"NZ_CP003996.1": "TTTTCCCC"}


def test_an_UNRESOLVABLE_contig_name_is_a_PERMANENT_failure_that_names_it(tmp_path):
    """⛔ Never retryable (the same assembly comes back tomorrow) and never read as a missing file."""
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=_StubStore(NCBI_HEADERS))

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(_ncbi_request(wanted_seqids=frozenset({"CP003999.1", "contig00007"})))
    assert "'contig00007'" in str(raised.value)
    assert STEM in str(raised.value)
    assert raised.value.is_retryable is False


def test_the_parsed_cache_holds_the_assembly_as_NCBI_keyed_it_and_re_keys_PER_CALL(tmp_path):
    """⚠ A reader walks loci within one genome; each locus names different contigs. One fetch, and every
    call gets its own names back."""
    ncbi = _StubStore(NCBI_HEADERS)
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=ncbi)

    first = provider.contig_sequences(_ncbi_request(wanted_seqids=frozenset({"CP003999.1"})))
    second = provider.contig_sequences(_ncbi_request(wanted_seqids=frozenset({"CP003996.1"})))
    everything = provider.contig_sequences(_ncbi_request(wanted_seqids=None))

    assert ncbi.calls == [BARE]
    assert first == {"CP003999.1": "ACGTACGT"}
    assert second == {"CP003996.1": "TTTTCCCC"}
    assert everything == NCBI_HEADERS


def test_a_BioSample_genome_is_NOT_re_keyed_even_when_the_request_names_seqids(tmp_path):
    """⚠ The held GFF and the ATB object write the same first token (280/280 measured), so the service's
    own "contig … is not in the assembly" lookup stays the one that speaks for them — unchanged."""
    atb = _StubStore({HELD_SEQID: STREAMED_BASES, "SAMN04393397.contig00002": "GGGG"})
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=atb, ncbi_store=None)

    sequences = provider.contig_sequences(
        GenomeSequenceRequest(
            sample_id="SAMN04393397",
            bakrep_dataset_id=None,
            sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
            wanted_seqids=frozenset({"nothing.like.this"}),
        )
    )

    assert set(sequences) == {HELD_SEQID, "SAMN04393397.contig00002"}


def test_the_NCBI_failures_RETRYABILITY_survives_the_hand_off(tmp_path):
    ncbi = _StubStore(error=AssemblyUnavailable("NCBI rate-limited this server", is_retryable=True))
    provider = ContigSequenceProvider(gff_root=tmp_path, atb_store=None, ncbi_store=ncbi)

    with pytest.raises(SequenceSourceUnavailable) as raised:
        provider.contig_sequences(_ncbi_request())
    assert raised.value.is_retryable is True


# ── NCBI configuration, with the same discipline as the ATB knobs ─────────────────────────────────
def test_NCBI_streaming_is_OFF_unless_configured_and_the_stations_say_so():
    provider = build_sequence_provider(Configuration(database_url="postgresql://x/y"))
    assert provider.streams_ncbi is False
    assert provider.stations == ()


def test_NCBI_streaming_is_ON_when_configured_and_ATB_stays_as_it_was(tmp_path):
    provider = build_sequence_provider(
        Configuration(
            database_url="postgresql://x/y",
            ncbi_streaming_enabled=True,
            ncbi_cache_directory=tmp_path / "ncbi",
        )
    )
    assert provider.streams_ncbi is True
    assert provider.streams_atb is False
    assert provider.can_stream is True
    assert provider.stations == ("ncbi_datasets",)
    assert (tmp_path / "ncbi").is_dir(), "the cache directory is created at start-up"


def test_both_stations_report_in_chain_order(tmp_path):
    provider = build_sequence_provider(
        Configuration(
            database_url="postgresql://x/y",
            artifact_roots={"gff": tmp_path},
            atb_streaming_enabled=True,
            ncbi_streaming_enabled=True,
        )
    )
    assert provider.stations == ("held_annotation_tree", "allthebacteria", "ncbi_datasets")


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
def test_the_NCBI_flag_accepts_the_same_words_every_other_flag_does(monkeypatch, value):
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv("BACATLAS_NCBI_STREAMING", value)
    assert Configuration.from_environment().ncbi_streaming_enabled is True


def test_the_NCBI_knobs_are_read_from_the_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv("BACATLAS_NCBI_STREAMING", "1")
    monkeypatch.setenv("BACATLAS_NCBI_CACHE_DIR", str(tmp_path / "ncbi"))
    monkeypatch.setenv("BACATLAS_NCBI_BASE_URL", "https://example.test/v2")
    monkeypatch.setenv("BACATLAS_NCBI_TIMEOUT_SECONDS", "12.5")
    monkeypatch.setenv("BACATLAS_NCBI_MAX_CONCURRENT_FETCHES", "2")
    monkeypatch.setenv("BACATLAS_NCBI_CACHE_BUDGET_BYTES", "1024")
    monkeypatch.setenv("BACATLAS_NCBI_REQUESTS_PER_SECOND", "10")
    monkeypatch.setenv("BACATLAS_NCBI_API_KEY", "s3cret")

    configuration = Configuration.from_environment()

    assert configuration.ncbi_cache_directory == tmp_path / "ncbi"
    assert configuration.ncbi_base_url == "https://example.test/v2"
    assert configuration.ncbi_timeout_seconds == 12.5
    assert configuration.ncbi_max_concurrent_fetches == 2
    assert configuration.ncbi_cache_budget_bytes == 1024
    assert configuration.ncbi_requests_per_second == 10.0
    assert configuration.ncbi_api_key == "s3cret"
    provider = build_sequence_provider(configuration)
    assert provider.streams_ncbi is True


@pytest.mark.parametrize(
    ("name", "value", "expected"),
    [
        ("BACATLAS_NCBI_REQUESTS_PER_SECOND", "fast", "not a number"),
        ("BACATLAS_NCBI_TIMEOUT_SECONDS", "soon", "not a number"),
        ("BACATLAS_NCBI_CACHE_BUDGET_BYTES", "big", "not an integer"),
        ("BACATLAS_NCBI_MAX_CONCURRENT_FETCHES", "four", "not an integer"),
    ],
)
def test_an_unparseable_NCBI_knob_RAISES_at_start_up(monkeypatch, name, value, expected):
    """⛔ The same rule as the ATB knobs: a typo must not masquerade as a setting that was never made."""
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql://x/y")
    monkeypatch.setenv(name, value)
    with pytest.raises(RuntimeError, match=expected):
        Configuration.from_environment()


def test_a_non_positive_NCBI_rate_fails_at_start_up_not_on_the_first_fetch():
    with pytest.raises(ValueError, match="requests_per_second"):
        build_sequence_provider(
            Configuration(database_url="postgresql://x/y", ncbi_streaming_enabled=True, ncbi_requests_per_second=0.0)
        )


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
