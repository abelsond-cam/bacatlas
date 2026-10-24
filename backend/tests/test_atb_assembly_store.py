"""`AtbAssemblyStore` — the fetch, the cache, and the failures that must not look alike.

⛔ **`urlopen` is stubbed throughout, and that is the point.** These tests are about OUR behaviour
on each answer ATB can give — a 404, a timeout, a truncated body, two workers racing for the same
genome — not about whether the bucket is up. The one test that really goes to the network is
`test_assembly_streaming.py`'s identity check, which is gated behind an environment flag.

⚠ The gzipped FASTA fixtures here are SYNTHETIC and tiny. Nothing in this file asserts anything
about real sequence; it asserts that the right bytes reach the parser and the wrong ones produce the
right sentence.
"""

from __future__ import annotations

import gzip
import io
import urllib.error
from pathlib import Path

import pytest

from bacatlas_backend.sequence_source import atb_assembly_store as store_module
from bacatlas_backend.sequence_source.atb_assembly_store import (
    AssemblyUnavailable,
    AtbAssemblyStore,
    AtbStoreSettings,
)

#: ⭐ An ATB defline, verbatim in shape: the seqid is the first token and five attributes follow it.
#: Our own Bakta GFF writes the bare seqid, so a reader that kept the whole line would key every
#: streamed contig under a string no `genome_contig` row can match.
ATB_DEFLINE = (
    ">SAMN04393397.contig00001 len=24 cov=30.8 corr=0 "
    "origname=NODE_1_length_296336_cov_30.837355_pilon sw=shovill-spades/1.0.4 date=20190103"
)
BASES = "ACGTACGTACGTACGTACGTACGT"


def _gzipped(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.open(buffer, "wt") as handle:
        handle.write(text)
    return buffer.getvalue()


ONE_CONTIG = _gzipped(f"{ATB_DEFLINE}\n{BASES[:12]}\n{BASES[12:]}\n")


class _Response(io.BytesIO):
    """Just enough of an `http.client.HTTPResponse` for a context-managed `read()`."""

    def __enter__(self):
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def _serving(payload: bytes, *, calls: list[str] | None = None):
    def urlopen(request, timeout=None):  # noqa: ARG001
        if calls is not None:
            calls.append(request.full_url)
        return _Response(payload)

    return urlopen


def _raising(error: Exception, *, calls: list[str] | None = None):
    def urlopen(request, timeout=None):  # noqa: ARG001
        if calls is not None:
            calls.append(request.full_url)
        raise error

    return urlopen


@pytest.fixture()
def cache_directory(tmp_path: Path) -> Path:
    return tmp_path / "seqcache"


# ── the URL, and the defline rule ──────────────────────────────────────────────────────────────
def test_the_object_url_is_the_flat_per_sample_key():
    """⛔ One object per BioSample. The EBI/OSF copies are ~4,000-genome `tar.xz` archives, which
    cannot serve a page — naming the wrong layout here is how that mistake would arrive."""
    store = AtbAssemblyStore(AtbStoreSettings(base_url="https://example.test/bucket/"))
    assert store.object_url("SAMN04393397") == "https://example.test/bucket/SAMN04393397.fa.gz"


def test_the_seqid_is_the_FIRST_TOKEN_and_not_the_whole_ATB_defline(monkeypatch, cache_directory):
    """⛔⛔ The join key. `genome_contig.seqid` holds `SAMN04393397.contig00001`; ATB's defline adds
    `len=`, `cov=`, `origname=` and `sw=`. Keeping them makes every contig lookup miss, and the
    endpoint then says *"this contig is not in the assembly"* — a sentence about the data."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    sequences = store.contig_sequences("SAMN04393397")

    assert list(sequences) == ["SAMN04393397.contig00001"]
    # ⚠ And unwrapped and upper-cased: our GFF wraps at 60, ATB at 70, so a reader that kept the
    # newlines would make two identical assemblies compare unequal.
    assert sequences["SAMN04393397.contig00001"] == BASES


# ── the two failures that must never be the same failure ───────────────────────────────────────
def test_a_404_is_PERMANENT_and_is_remembered_so_it_is_asked_only_once(monkeypatch, cache_directory):
    """⛔ ATB not having a genome is a fact that will not change in a session. Without the negative
    cache every page view of such a genome pays a full round trip to learn the same 404."""
    calls: list[str] = []
    error = urllib.error.HTTPError("https://example.test/x.fa.gz", 404, "Not Found", {}, None)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(error, calls=calls))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    for _ in range(3):
        with pytest.raises(AssemblyUnavailable) as raised:
            store.contig_sequences("SAMN00000000")
        assert raised.value.is_retryable is False
        assert "no assembly for SAMN00000000" in str(raised.value)

    assert len(calls) == 1, "a known-absent genome must not be re-fetched"


def test_an_expired_negative_cache_entry_is_asked_AGAIN(monkeypatch, cache_directory):
    """⚠ Permanent-for-now, not permanent-forever: ATB publishes incremental releases, so a genome
    absent today can exist next month. A process that never re-asks would never see it."""
    calls: list[str] = []
    error = urllib.error.HTTPError("https://example.test/x.fa.gz", 404, "Not Found", {}, None)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(error, calls=calls))
    store = AtbAssemblyStore(
        AtbStoreSettings(cache_directory=cache_directory, absent_ttl_seconds=0.0)
    )

    for _ in range(2):
        with pytest.raises(AssemblyUnavailable):
            store.contig_sequences("SAMN00000000")
    assert len(calls) == 2


@pytest.mark.parametrize(
    ("error", "expected"),
    [
        (urllib.error.URLError("name resolution failed"), "could not be reached"),
        (TimeoutError("timed out"), "could not be reached"),
        (
            urllib.error.HTTPError("https://example.test/x.fa.gz", 503, "Slow Down", {}, None),
            "HTTP 503",
        ),
    ],
)
def test_a_reachability_failure_is_RETRYABLE_and_says_so(monkeypatch, cache_directory, error, expected):
    """⛔ The distinction the `Retry-After` header rests on. A timeout told as "ATB has no such
    genome" sends the reader away from sequence that is in fact available."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(error))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences("SAMN04393397")
    assert raised.value.is_retryable is True
    assert expected in str(raised.value)


def test_a_transient_failure_is_NOT_remembered_as_an_absence(monkeypatch, cache_directory):
    """⛔ Otherwise one timeout poisons the genome for the whole TTL — the page would report a
    genome as missing from ATB because our network blinked once."""
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(TimeoutError("slow")))
    with pytest.raises(AssemblyUnavailable):
        store.contig_sequences("SAMN04393397")

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    assert store.contig_sequences("SAMN04393397")["SAMN04393397.contig00001"] == BASES


# ── the disk cache ─────────────────────────────────────────────────────────────────────────────
def test_a_cached_assembly_is_not_fetched_twice(monkeypatch, cache_directory):
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG, calls=calls))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    store.contig_sequences("SAMN04393397")
    AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory)).contig_sequences("SAMN04393397")

    assert len(calls) == 1, "the second store should have read the first one's cache file"
    assert (cache_directory / "SAMN04393397.fa.gz").exists()


def test_a_TRUNCATED_cache_entry_is_DISCARDED_rather_than_merely_reported(monkeypatch, cache_directory):
    """⛔⛔ The one failure that survives a restart. A half-written cache file fails to decompress
    for every later reader, so a store that only *reports* it breaks that genome permanently."""
    cache_directory.mkdir(parents=True)
    corrupt = cache_directory / "SAMN04393397.fa.gz"
    corrupt.write_bytes(ONE_CONTIG[: len(ONE_CONTIG) // 2])
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences("SAMN04393397")
    assert raised.value.is_retryable is True
    assert not corrupt.exists(), "the corrupt entry must be gone, not just complained about"

    # And the retry the message promises actually works.
    assert store.contig_sequences("SAMN04393397")["SAMN04393397.contig00001"] == BASES


def test_an_assembly_with_no_records_is_a_PERMANENT_failure_not_an_empty_map(monkeypatch, cache_directory):
    """⚠ Returning `{}` would hand back a map in which every contig lookup misses, and the endpoint
    would blame the contig rather than the object."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(_gzipped("")))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences("SAMN04393397")
    assert raised.value.is_retryable is False
    assert "no FASTA records" in str(raised.value)


def test_the_cache_is_held_to_its_BUDGET(monkeypatch, cache_directory):
    """⛔ An unbounded cache on an 80 GB instance disk is a storage bug, not a performance win —
    and the global CLAUDE.md rule is that nothing lands outside a budgeted directory."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    # A budget of one file, so the second fetch must evict the first.
    store = AtbAssemblyStore(
        AtbStoreSettings(cache_directory=cache_directory, cache_budget_bytes=len(ONE_CONTIG))
    )
    store.contig_sequences("SAMN00000001")
    store.contig_sequences("SAMN00000002")

    held = sorted(path.name for path in cache_directory.glob("*.fa.gz"))
    assert len(held) == 1, f"budget not enforced: {held}"


def test_no_part_files_are_left_behind(monkeypatch, cache_directory):
    """⚠ The cache is written to a unique temporary name then renamed, so two workers fetching one
    genome cannot interleave into a half-written file. A leaked `.part` would grow without bound."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))
    store.contig_sequences("SAMN04393397")

    assert list(cache_directory.glob("*.part")) == []


def test_a_sample_id_that_is_not_a_bare_basename_is_NOT_turned_into_a_path(monkeypatch, cache_directory):
    """⛔ `sample_id` reaches the store from the DATABASE, so it is already constrained — but a cache
    path assembled from a value that ever becomes caller-supplied is the classic traversal. The
    guard degrades to "no disk cache for this one", never to a write outside the directory."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory))

    store.contig_sequences("../escape")

    assert not (cache_directory.parent / "escape.fa.gz").exists()
    assert list(cache_directory.glob("**/*.fa.gz")) == []


def test_with_no_cache_directory_the_store_still_works(monkeypatch):
    """⚠ `None` is a legitimate configuration — it is what a test uses — and it must not crash."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=None))
    assert store.contig_sequences("SAMN04393397")["SAMN04393397.contig00001"] == BASES


# ── the worker bound ───────────────────────────────────────────────────────────────────────────
def test_a_request_that_cannot_get_a_fetch_SLOT_is_told_to_retry(monkeypatch, cache_directory):
    """⛔ gunicorn serves this synchronously on 2 vCPU: N concurrent cold opens hold N workers for
    ~0.27 s each. At human concurrency that is invisible; under a crawler it is a denial of service
    against our own API. The honest answer arrives faster than the queue would."""
    monkeypatch.setattr(store_module, "SLOT_WAIT_SECONDS", 0.01)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory, max_concurrent_fetches=1))
    # Take the only slot and do not give it back.
    assert store._slots.acquire(timeout=1)

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences("SAMN04393397")
    assert raised.value.is_retryable is True
    assert "as many assemblies as it will fetch at once" in str(raised.value)


def test_a_failed_fetch_RELEASES_its_slot(monkeypatch, cache_directory):
    """⛔ Otherwise the first N failures permanently exhaust the pool and every later request is
    told the server is busy — with nothing in flight at all."""
    monkeypatch.setattr(store_module, "SLOT_WAIT_SECONDS", 0.01)
    store = AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory, max_concurrent_fetches=1))

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(TimeoutError("slow")))
    with pytest.raises(AssemblyUnavailable):
        store.contig_sequences("SAMN00000001")

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ONE_CONTIG))
    assert store.contig_sequences("SAMN00000002")["SAMN04393397.contig00001"] == BASES


def test_identity_encoding_is_REQUESTED(monkeypatch, cache_directory):
    """⛔ The object is already gzip. Letting the transport add a second `Content-Encoding` is how
    `serving_at_scale.md` §5's measured GitHub Pages trap worked — a range or a body over the
    *compressed* representation, returning plausible wrong bytes with no error."""
    seen: dict[str, str] = {}

    def urlopen(request, timeout=None):  # noqa: ARG001
        seen.update(request.headers)
        return _Response(ONE_CONTIG)

    monkeypatch.setattr(store_module.urllib.request, "urlopen", urlopen)
    AtbAssemblyStore(AtbStoreSettings(cache_directory=cache_directory)).contig_sequences("SAMN04393397")

    assert seen.get("Accept-encoding") == "identity"
