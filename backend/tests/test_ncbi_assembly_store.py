"""`NcbiAssemblyStore` — the fetch, the zip, the cache, the rate limit, and the failures that must not look alike.

⛔ **`urlopen` is stubbed throughout, and that is the point.** These tests are about OUR behaviour on each
answer NCBI Datasets can give — a package with the member, a package without it, the abandoned stream an
unknown accession really returns, a 404, a 429, a timeout — not about whether NCBI is up. The one test
that goes to the network is gated behind `BACATLAS_LIVE_NCBI=1` and skips loudly otherwise.

⚠ The packages here are SYNTHETIC zips built with `zipfile`, shaped like the one the live API returned on
2026-11-13 (README, catalogue, one `ncbi_dataset/data/<accession>/<stem>_genomic.fna`). Nothing asserts
anything about real sequence; it asserts that the right bytes reach the parser and the wrong ones
produce the right sentence.
"""

from __future__ import annotations

import io
import os
import urllib.error
import zipfile
from pathlib import Path

import pytest

from bacatlas_backend.sequence_source import ncbi_assembly_store as store_module
from bacatlas_backend.sequence_source.atb_assembly_store import AssemblyUnavailable
from bacatlas_backend.sequence_source.ncbi_assembly_store import (
    NcbiAssemblyStore,
    NcbiStoreSettings,
    TokenBucket,
)

ACCESSION = "GCF_000512165.1"
MEMBER = f"ncbi_dataset/data/{ACCESSION}/{ACCESSION}_ASM51216v1_genomic.fna"
#: ⭐ Verbatim in shape from the live package: the first token is the RefSeq-prefixed accession and a
#: free-text description follows. Keeping the whole line would key every contig under a string no
#: `genome_contig` row — and no `NZ_` rule — could ever match.
FASTA = (
    ">NZ_CP003999.1 Klebsiella pneumoniae subsp. pneumoniae Kp13 chromosome, complete genome\n"
    "acgtACGTacgt\nACGTACGTACGT\n"
    ">NZ_CP003996.1 Klebsiella pneumoniae subsp. pneumoniae Kp13 plasmid pKP13a, complete sequence\n"
    "TTTTCCCC\n"
)


def _package(members: dict[str, str]) -> bytes:
    """A data package the way NCBI lays one out, with `members` under `ncbi_dataset/data/`."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.md", "# NCBI Datasets\n")
        archive.writestr("ncbi_dataset/data/dataset_catalog.json", '{"apiVersion": "V2"}')
        for name, text in members.items():
            archive.writestr(name, text)
        archive.writestr("md5sum.txt", "")
    return buffer.getvalue()


PACKAGE = _package({MEMBER: FASTA})
#: A valid zip that carries no assembly — one of the two shapes "no such accession" can take.
EMPTY_PACKAGE = _package({})
#: ⛔ The shape the live API ACTUALLY returned for an unknown accession (2026-11-13): HTTP 200, 971 bytes,
#: a local file header for README.md and no central directory — not a zip at all to `zipfile`.
ABANDONED_STREAM = b"PK\x03\x04" + b"\x00" * 26 + b"README.md" + b"# NCBI Datasets\n" + b"\x00" * 900
#: A real download cut off before its central directory: every local header is present, so the member's
#: path IS in the bytes, and that is what tells a truncation from an absence.
TRUNCATED_PACKAGE = PACKAGE[: PACKAGE.index(b"PK\x01\x02")]


class _Response(io.BytesIO):
    """Just enough of an `http.client.HTTPResponse` for a context-managed `read()`."""

    def __enter__(self):
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def _serving(payload: bytes, *, calls: list[str] | None = None, headers: dict[str, str] | None = None):
    def urlopen(request, timeout=None):  # noqa: ARG001
        if calls is not None:
            calls.append(request.full_url)
        if headers is not None:
            headers.update(request.headers)
        return _Response(payload)

    return urlopen


def _raising(error: Exception, *, calls: list[str] | None = None):
    def urlopen(request, timeout=None):  # noqa: ARG001
        if calls is not None:
            calls.append(request.full_url)
        raise error

    return urlopen


def _serving_any_accession(*, calls: list[str] | None = None):
    """A stub NCBI that has every accession: the package is built for whichever one the URL names.

    ⚠ Needed because the member lookup is scoped to the accession's OWN directory, so a fixed package
    for `GCF_000512165.1` is — correctly — "absent" when `GCF_000000002.1` is asked for.
    """

    def urlopen(request, timeout=None):  # noqa: ARG001
        if calls is not None:
            calls.append(request.full_url)
        accession = request.full_url.split("/accession/")[1].split("/")[0]
        return _Response(_package({f"ncbi_dataset/data/{accession}/{accession}_x_genomic.fna": FASTA}))

    return urlopen


def _http_error(code: int, reason: str) -> urllib.error.HTTPError:
    return urllib.error.HTTPError("https://example.test/download", code, reason, {}, None)


class _FakeClock:
    """A clock that moves only when something sleeps on it, so spacing is asserted, not waited for."""

    def __init__(self) -> None:
        self.now = 1000.0
        self.sleeps: list[float] = []

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.sleeps.append(round(seconds, 6))
        self.now += seconds


@pytest.fixture()
def cache_directory(tmp_path: Path) -> Path:
    return tmp_path / "seqcache" / "ncbi"


def _store(cache_directory: Path | None, **overrides) -> NcbiAssemblyStore:
    return NcbiAssemblyStore(NcbiStoreSettings(cache_directory=cache_directory, **overrides))


# ── the URL, the request, and the defline rule ────────────────────────────────────────────────────
def test_the_download_url_is_the_datasets_v2_per_accession_package():
    store = NcbiAssemblyStore(NcbiStoreSettings(base_url="https://example.test/datasets/v2/"))
    assert store.object_url(ACCESSION) == (
        "https://example.test/datasets/v2/genome/accession/GCF_000512165.1/download"
        "?include_annotation_type=GENOME_FASTA"
    )


def test_the_accession_is_url_quoted():
    """⚠ The accession reaches here from the database, but a value with a `/` must not become a path."""
    store = NcbiAssemblyStore(NcbiStoreSettings(base_url="https://example.test"))
    assert "/accession/a%2Fb/download" in store.object_url("a/b")


def test_contigs_are_keyed_by_the_FIRST_TOKEN_of_the_defline_and_upper_cased(monkeypatch, cache_directory):
    """⛔⛔ The join key, before the `NZ_` rule even runs. NCBI's defline is `>NZ_CP003999.1 Klebsiella
    pneumoniae …`; keeping it would make every lookup miss with a sentence about the data."""
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))

    sequences = _store(cache_directory).contig_sequences(ACCESSION)

    assert list(sequences) == ["NZ_CP003999.1", "NZ_CP003996.1"]
    assert sequences["NZ_CP003999.1"] == "ACGTACGTACGTACGTACGTACGT"
    assert sequences["NZ_CP003996.1"] == "TTTTCCCC"


def test_the_request_asks_for_a_zip_with_identity_encoding_and_NO_api_key_by_default(monkeypatch, cache_directory):
    """⛔ The zip is already compressed; a transport that adds a second encoding is the GitHub Pages trap.
    And no `api-key` header is sent unless one is configured — an empty header is a malformed request."""
    seen: dict[str, str] = {}
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE, headers=seen))

    _store(cache_directory).contig_sequences(ACCESSION)

    assert seen.get("Accept") == "application/zip"
    assert seen.get("Accept-encoding") == "identity"
    assert "Api-key" not in seen and "api-key" not in seen


def test_the_api_key_is_sent_as_a_header_when_configured(monkeypatch, cache_directory):
    seen: dict[str, str] = {}
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE, headers=seen))

    _store(cache_directory, api_key="s3cret").contig_sequences(ACCESSION)

    # `urllib` title-cases header names on the way in.
    assert seen.get("Api-key") == "s3cret"


# ── absence, in the three shapes it takes — each permanent, each remembered ───────────────────────
def test_a_package_with_NO_member_for_the_accession_is_an_ABSENCE_and_is_asked_only_once(
    monkeypatch, cache_directory
):
    """⛔ NCBI's "not found" is a 200 — so the negative cache keys on the BODY, not the status."""
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(EMPTY_PACKAGE, calls=calls))
    store = _store(cache_directory)

    for _ in range(3):
        with pytest.raises(AssemblyUnavailable) as raised:
            store.contig_sequences(ACCESSION)
        assert raised.value.is_retryable is False
        assert "no assembly for GCF_000512165.1" in str(raised.value)

    assert len(calls) == 1, "a known-absent accession must not be re-fetched"
    assert list(cache_directory.glob("*.zip")) == [], "an empty package must never be cached as the assembly"


def test_the_ABANDONED_STREAM_an_unknown_accession_really_returns_is_an_absence_not_a_truncation(
    monkeypatch, cache_directory
):
    """⭐ Measured 2026-11-13: `GCF_000000000.1` → 200, 971 bytes, no central directory. It names no
    member, so it is NCBI saying "nothing here" — permanent, remembered, and not a retry loop."""
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(ABANDONED_STREAM, calls=calls))
    store = _store(cache_directory)

    for _ in range(2):
        with pytest.raises(AssemblyUnavailable) as raised:
            store.contig_sequences("GCF_000000000.1")
        assert raised.value.is_retryable is False

    assert len(calls) == 1
    assert list(cache_directory.glob("*")) == []


def test_a_404_is_PERMANENT_and_remembered_too(monkeypatch, cache_directory):
    """⚠ Not what the live API does today, but the obvious thing for it to start doing."""
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(_http_error(404, "Not Found"), calls=calls))
    store = _store(cache_directory)

    for _ in range(3):
        with pytest.raises(AssemblyUnavailable) as raised:
            store.contig_sequences(ACCESSION)
        assert raised.value.is_retryable is False

    assert len(calls) == 1


def test_an_expired_negative_cache_entry_is_asked_AGAIN(monkeypatch, cache_directory):
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(EMPTY_PACKAGE, calls=calls))
    store = _store(cache_directory, absent_ttl_seconds=0.0)

    for _ in range(2):
        with pytest.raises(AssemblyUnavailable):
            store.contig_sequences(ACCESSION)
    assert len(calls) == 2


def test_a_member_filed_under_ANOTHER_accession_is_not_taken(monkeypatch, cache_directory):
    """⛔ Bases under the wrong accession are plausible wrong bases. The package is treated as holding
    nothing for the accession asked — which is the truth about it."""
    other = _package({"ncbi_dataset/data/GCF_999999999.1/GCF_999999999.1_x_genomic.fna": FASTA})
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(other))

    with pytest.raises(AssemblyUnavailable) as raised:
        _store(cache_directory).contig_sequences(ACCESSION)
    assert raised.value.is_retryable is False


def test_two_genomic_members_are_REFUSED_rather_than_guessed(monkeypatch, cache_directory):
    two = _package(
        {
            MEMBER: FASTA,
            f"ncbi_dataset/data/{ACCESSION}/{ACCESSION}_other_genomic.fna": ">x\nAAAA\n",
        }
    )
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(two))

    with pytest.raises(AssemblyUnavailable, match="refusing to guess") as raised:
        _store(cache_directory).contig_sequences(ACCESSION)
    assert raised.value.is_retryable is False


def test_cds_and_rna_FROM_genomic_members_are_not_mistaken_for_the_assembly(monkeypatch, cache_directory):
    package = _package(
        {
            MEMBER: FASTA,
            f"ncbi_dataset/data/{ACCESSION}/cds_from_genomic.fna": ">cds\nATG\n",
            f"ncbi_dataset/data/{ACCESSION}/rna_from_genomic.fna": ">rna\nAUG\n",
        }
    )
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(package))

    assert list(_store(cache_directory).contig_sequences(ACCESSION)) == ["NZ_CP003999.1", "NZ_CP003996.1"]


def test_an_assembly_with_no_records_is_a_PERMANENT_failure_not_an_empty_map(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(_package({MEMBER: ""})))

    with pytest.raises(AssemblyUnavailable) as raised:
        _store(cache_directory).contig_sequences(ACCESSION)
    assert raised.value.is_retryable is False
    assert "no FASTA records" in str(raised.value)


# ── the transient failures, each RETRYABLE and none remembered as an absence ──────────────────────
@pytest.mark.parametrize(
    ("error", "expected"),
    [
        (urllib.error.URLError("name resolution failed"), "could not be reached"),
        (TimeoutError("timed out"), "could not be reached"),
        (_http_error(503, "Service Unavailable"), "HTTP 503"),
        (_http_error(429, "Too Many Requests"), "rate-limited"),
    ],
)
def test_a_reachability_failure_is_RETRYABLE_and_says_so(monkeypatch, cache_directory, error, expected):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(error))

    with pytest.raises(AssemblyUnavailable) as raised:
        _store(cache_directory).contig_sequences(ACCESSION)
    assert raised.value.is_retryable is True
    assert expected in str(raised.value)


def test_a_TRUNCATED_download_is_retryable_and_NOTHING_is_cached(monkeypatch, cache_directory):
    """⛔ The member's path is in the bytes, so this is a cut-off download, not an absence — and the
    half-package must not land in the cache, where it would fail every later reader."""
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(TRUNCATED_PACKAGE, calls=calls))
    store = _store(cache_directory)

    with pytest.raises(AssemblyUnavailable, match="truncated") as raised:
        store.contig_sequences(ACCESSION)
    assert raised.value.is_retryable is True
    assert list(cache_directory.glob("*")) == []

    # Not remembered as absent: the next attempt goes to the network again, and succeeds.
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE, calls=calls))
    assert store.contig_sequences(ACCESSION)["NZ_CP003996.1"] == "TTTTCCCC"
    assert len(calls) == 2


def test_a_transient_failure_is_NOT_remembered_as_an_absence(monkeypatch, cache_directory):
    store = _store(cache_directory)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(TimeoutError("slow")))
    with pytest.raises(AssemblyUnavailable):
        store.contig_sequences(ACCESSION)

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))
    assert store.contig_sequences(ACCESSION)["NZ_CP003996.1"] == "TTTTCCCC"


def test_a_400_is_permanent_but_NOT_remembered_as_an_absence(monkeypatch, cache_directory):
    """⚠ A rejected accession SHAPE says nothing about the assembly's existence."""
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(_http_error(400, "Bad Request"), calls=calls))
    store = _store(cache_directory)

    for _ in range(2):
        with pytest.raises(AssemblyUnavailable) as raised:
            store.contig_sequences(ACCESSION)
        assert raised.value.is_retryable is False
    assert len(calls) == 2


# ── the disk cache ────────────────────────────────────────────────────────────────────────────────
def test_a_cached_package_is_not_fetched_twice_and_is_SEPARATE_from_the_ATB_cache(monkeypatch, cache_directory):
    calls: list[str] = []
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE, calls=calls))

    _store(cache_directory).contig_sequences(ACCESSION)
    again = _store(cache_directory).contig_sequences(ACCESSION)

    assert len(calls) == 1, "the second store should have read the first one's cache file"
    assert again["NZ_CP003999.1"] == "ACGTACGTACGTACGTACGTACGT"
    assert (cache_directory / f"{ACCESSION}.zip").exists()
    # ⛔ Nothing of NCBI's lands beside the ATB objects one level up.
    assert list(cache_directory.parent.glob("*.zip")) == []


def test_a_TRUNCATED_cache_entry_is_DISCARDED_rather_than_merely_reported(monkeypatch, cache_directory):
    cache_directory.mkdir(parents=True)
    corrupt = cache_directory / f"{ACCESSION}.zip"
    corrupt.write_bytes(TRUNCATED_PACKAGE)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))
    store = _store(cache_directory)

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences(ACCESSION)
    assert raised.value.is_retryable is True
    assert not corrupt.exists(), "the corrupt entry must be gone, not just complained about"

    assert store.contig_sequences(ACCESSION)["NZ_CP003996.1"] == "TTTTCCCC"


def test_the_cache_is_held_to_its_BUDGET(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving_any_accession())
    # A budget of one package, so the second fetch must evict the first.
    store = _store(cache_directory, cache_budget_bytes=len(PACKAGE) + 16)

    store.contig_sequences("GCF_000000001.1")
    store.contig_sequences("GCF_000000002.1")

    held = sorted(path.name for path in cache_directory.glob("*.zip"))
    assert len(held) == 1, f"budget not enforced: {held}"


def test_no_part_files_are_left_behind(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))
    _store(cache_directory).contig_sequences(ACCESSION)
    assert list(cache_directory.glob("*.part")) == []


def test_an_accession_that_is_not_a_bare_basename_is_NOT_turned_into_a_path(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(EMPTY_PACKAGE))
    with pytest.raises(AssemblyUnavailable):
        _store(cache_directory).contig_sequences("../escape")

    assert not (cache_directory.parent / "escape.zip").exists()
    assert list(cache_directory.glob("**/*.zip")) == []


def test_with_no_cache_directory_the_store_still_works(monkeypatch):
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))
    assert _store(None).contig_sequences(ACCESSION)["NZ_CP003996.1"] == "TTTTCCCC"


# ── the rate limit ────────────────────────────────────────────────────────────────────────────────
def test_the_token_bucket_SPACES_requests_at_the_rate_after_one_seconds_burst():
    """At 2/s: two tokens are banked, so two requests go at once; the third and fourth each wait 0.5 s."""
    clock = _FakeClock()
    bucket = TokenBucket(2.0, clock=clock, sleep=clock.sleep)

    waits = [bucket.acquire() for _ in range(4)]

    assert waits == [0.0, 0.0, 0.5, 0.5]
    assert clock.sleeps == [0.5, 0.5]


def test_the_token_bucket_REFILLS_while_idle_but_never_beyond_one_second():
    clock = _FakeClock()
    bucket = TokenBucket(3.0, clock=clock, sleep=clock.sleep)
    for _ in range(3):
        bucket.acquire()
    clock.now += 60.0  # a minute idle banks one second's worth, not sixty

    waits = [bucket.acquire() for _ in range(4)]

    assert waits[:3] == [0.0, 0.0, 0.0]
    assert waits[3] == pytest.approx(1 / 3)


def test_the_store_WAITS_on_its_bucket_before_the_second_fetch(monkeypatch, cache_directory):
    """⚠ Through the store, not the bucket alone: the token has to be taken on the fetch path."""
    clock = _FakeClock()
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving_any_accession())
    store = NcbiAssemblyStore(
        NcbiStoreSettings(cache_directory=cache_directory),
        rate_limiter=TokenBucket(1.0, clock=clock, sleep=clock.sleep),
    )

    store.contig_sequences("GCF_000000001.1")
    store.contig_sequences("GCF_000000002.1")
    store.contig_sequences("GCF_000000001.1")  # cached: no token spent

    assert clock.sleeps == [1.0]


def test_a_cache_hit_spends_no_token_and_touches_no_network(monkeypatch, cache_directory):
    calls: list[str] = []
    clock = _FakeClock()
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE, calls=calls))
    store = NcbiAssemblyStore(
        NcbiStoreSettings(cache_directory=cache_directory),
        rate_limiter=TokenBucket(1.0, clock=clock, sleep=clock.sleep),
    )

    for _ in range(5):
        store.contig_sequences(ACCESSION)

    assert calls == [store.object_url(ACCESSION)]
    assert clock.sleeps == []


def test_a_non_positive_rate_is_refused_at_construction():
    """⛔ At start-up, where it belongs — not as a fetch that waits forever."""
    with pytest.raises(ValueError, match="requests_per_second"):
        NcbiAssemblyStore(NcbiStoreSettings(requests_per_second=0))


# ── the worker bound ──────────────────────────────────────────────────────────────────────────────
def test_a_request_that_cannot_get_a_fetch_SLOT_is_told_to_retry(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module, "SLOT_WAIT_SECONDS", 0.01)
    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving(PACKAGE))
    store = _store(cache_directory, max_concurrent_fetches=1)
    assert store._slots.acquire(timeout=1)

    with pytest.raises(AssemblyUnavailable) as raised:
        store.contig_sequences(ACCESSION)
    assert raised.value.is_retryable is True
    assert "as many assemblies from NCBI" in str(raised.value)


def test_a_failed_fetch_RELEASES_its_slot(monkeypatch, cache_directory):
    monkeypatch.setattr(store_module, "SLOT_WAIT_SECONDS", 0.01)
    store = _store(cache_directory, max_concurrent_fetches=1)

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _raising(TimeoutError("slow")))
    with pytest.raises(AssemblyUnavailable):
        store.contig_sequences("GCF_000000001.1")

    monkeypatch.setattr(store_module.urllib.request, "urlopen", _serving_any_accession())
    assert store.contig_sequences("GCF_000000002.1")["NZ_CP003996.1"] == "TTTTCCCC"


# ── ⭐ the one real request, opt-in ───────────────────────────────────────────────────────────────
LIVE_NCBI = os.environ.get("BACATLAS_LIVE_NCBI", "").strip().lower() in {"1", "true", "yes", "on"}


@pytest.mark.skipif(
    not LIVE_NCBI,
    reason="BACATLAS_LIVE_NCBI is not set — this test fetches GCF_000512165.1 from NCBI Datasets over the "
    "network, which a CI run should not do by default",
)
def test_the_live_package_for_GCF_000512165_1_parses_to_uniquely_keyed_contigs(tmp_path):
    """⭐ The real zip against the real reader. Measured 2026-11-13: 1,699,070 bytes, 1.08 s, six contigs
    (one chromosome, five plasmids), every defline `>NZ_CP0039xx.1 …`."""
    store = NcbiAssemblyStore(NcbiStoreSettings(cache_directory=tmp_path))
    sequences = store.contig_sequences(ACCESSION)

    assert len(sequences) >= 1
    assert all(bases and bases == bases.upper() for bases in sequences.values())
    # ⛔ "Every header token is unique" cannot be read off a dict — a collision would have silently
    # overwritten. Count the deflines in the cached member and compare.
    with zipfile.ZipFile(tmp_path / f"{ACCESSION}.zip") as archive:
        member = next(name for name in archive.namelist() if name.endswith("_genomic.fna"))
        deflines = sum(1 for line in archive.open(member) if line.startswith(b">"))
    assert deflines == len(sequences), f"{deflines} deflines but {len(sequences)} distinct first tokens"
    # And the second read is served from the cache and agrees exactly.
    assert NcbiAssemblyStore(NcbiStoreSettings(cache_directory=tmp_path)).contig_sequences(ACCESSION) == sequences
