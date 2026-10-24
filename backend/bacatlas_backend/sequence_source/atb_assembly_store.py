"""One genome's assembly, fetched from AllTheBacteria on demand and cached on disk.

⭐ **Why this exists.** The served sequence tier is ~2.1 MB of gzipped Bakta GFF per genome
(measured, mean over all 280 probe genomes), of which **81 % is the FASTA block**. At the 79,056
Klebsiella genomes that carry a GFF+assembly that is ~162 GB on the instance disk, where AWS charges
$0.10/GB-month — 8× what the same cold bytes cost in object storage. ATB publishes **one object per
BioSample**, so the bytes need not be held at all.

⛔ **The bytes are the SAME bytes, and that is measured rather than assumed.** On 2026-10-23 all 280
probe genomes were fetched from ATB and compared against the ``##FASTA`` block of the Bakta GFF we
hold, normalised identically: **280/280 byte-identical, 0 fetch failures**, and the seqid sets equal
on every genome. BakRep ran Bakta over the ATB assemblies, so this is one upstream object read two
ways — but it is now *checked*, because a different assembly would have produced plausible wrong
bases at correct-looking coordinates, which no page can detect.

⚠ **Latency, measured over those 280 fetches** (UK domestic line → ``eu-west-2``): mean **0.241 s**,
median 0.230 s, p95 0.330 s, max 0.510 s. Parse is a further ~26 ms. So a cold Sequence-tab open is
~0.27 s against ~0.03 s warm. That is a tab a reader clicks, not a hot path.

⛔ **The bucket is open-access and OWNER-pays** — an anonymous GET returns 200 with no
``x-amz-request-charged`` header, and Lightsail inbound transfer is $0.00/GB. Streaming costs us
nothing per request. ⚠ But it is a **third party**, and ATB has already moved once (release 0.2 was
distributed from OSF; the per-sample objects are now on S3). So this is written as a *cache with a
fallback*, never a dependency: a genome whose GFF we hold is never fetched, and a fetch that fails
is a named 503 rather than a wrong answer.

⚠ **No new dependency.** The backend's serving requirements are flask/sqlalchemy/psycopg/alembic;
this uses ``urllib.request`` from the standard library rather than adding ``requests`` or ``boto3``
to a production image for one GET.
"""

from __future__ import annotations

import gzip
import io
import logging
import os
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from bacatlas_backend.gff.fasta_reader import read_fasta_block

logger = logging.getLogger(__name__)

#: The published per-sample object layout: a flat keyspace, one gzipped FASTA per sample accession.
#: ⚠ **Not** the EBI FTP or OSF copies — those are ``tar.xz`` archives of ~4,000 genomes batched by
#: species call, and ATB's own documentation says getting one sample from them means downloading and
#: extracting the whole archive. That shape cannot serve a page.
DEFAULT_ATB_BASE_URL: Final = "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com"

#: ⚠ Generous against the measured p95 of 0.330 s, because the cost of being wrong is asymmetric: a
#: timeout that is too tight turns a slow-but-fine fetch into a 503, while one that is too loose
#: only holds a worker that the semaphore has already bounded.
DEFAULT_TIMEOUT_SECONDS: Final = 10.0

#: ⛔ **The worker bound, and the reason it is not optional.** gunicorn on 2 vCPU serves this
#: synchronously, so N concurrent cold Sequence-tab opens hold N workers for ~0.27 s each. At human
#: concurrency that is invisible; under a crawler it is a denial of service against our own API. A
#: request that cannot get a slot is told to retry — which is true — instead of queueing behind a
#: fetch the reader has already navigated away from.
DEFAULT_MAX_CONCURRENT_FETCHES: Final = 4

#: How long to wait for a fetch slot before saying "busy". Short on purpose: the honest answer
#: arrives faster than the queue would.
SLOT_WAIT_SECONDS: Final = 2.0

#: ⛔ **A negative cache, because absence is a fact about ATB and will not change within a session.**
#: ATB's documentation notes ~10,000 assemblies are unavailable from ENA owing to submission errors;
#: whether the S3 mirror shares those holes is unmeasured for our cohort. Without this, every page
#: view of such a genome pays a full round trip to learn the same 404.
DEFAULT_ABSENT_TTL_SECONDS: Final = 600.0

#: The disk cache budget. ⛔ An unbounded cache on an 80 GB instance disk is a storage bug, not a
#: performance win — 20 GB holds ~11,700 assemblies at the measured 1.703 MB mean, which is far more
#: working set than any reading session has.
DEFAULT_CACHE_BUDGET_BYTES: Final = 20 * 1024**3

#: Evicting one file per fetch would thrash at the boundary; a batch amortises the ``scandir``.
EVICTION_BATCH: Final = 64


class AssemblyUnavailable(LookupError):
    """An assembly could not be obtained — and ``is_retryable`` says whether asking again may help.

    ⛔ The distinction is the whole value of this exception. *"ATB does not have this genome"* is
    permanent and the reader should be told so; *"the fetch timed out"* is transient and a retry is
    the right advice. Collapsing them produces a page that either nags about a genome that will
    never arrive, or gives up on one that was merely slow.
    """

    def __init__(self, message: str, *, is_retryable: bool) -> None:
        super().__init__(message)
        self.is_retryable = is_retryable


@dataclass(frozen=True)
class AtbStoreSettings:
    """Everything the store reads from the outside, resolved once at start-up."""

    base_url: str = DEFAULT_ATB_BASE_URL
    cache_directory: Path | None = None
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS
    max_concurrent_fetches: int = DEFAULT_MAX_CONCURRENT_FETCHES
    absent_ttl_seconds: float = DEFAULT_ABSENT_TTL_SECONDS
    cache_budget_bytes: int = DEFAULT_CACHE_BUDGET_BYTES


class AtbAssemblyStore:
    """Fetch, cache and parse ATB assemblies. Thread-safe; one instance per process.

    ⚠ **The disk cache is keyed on the sample accession and holds the GZIPPED bytes**, not the
    parsed dict: the parse is 26 ms and the gzip is 1.7 MB against ~5.3 MB parsed, so caching the
    compressed form is both smaller on disk and a smaller thing to get wrong. The parsed form is
    cached one layer up, by `ContigSequenceProvider`, where the local-GFF path shares it.
    """

    def __init__(self, settings: AtbStoreSettings | None = None) -> None:
        self._settings = settings or AtbStoreSettings()
        self._slots = threading.BoundedSemaphore(max(1, self._settings.max_concurrent_fetches))
        self._absent: dict[str, float] = {}
        self._absent_lock = threading.Lock()
        if self._settings.cache_directory is not None:
            self._settings.cache_directory.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ public

    def object_url(self, sample_id: str) -> str:
        """The object this accession would be at. Separate so a test can assert it without a GET."""
        return f"{self._settings.base_url.rstrip('/')}/{sample_id}.fa.gz"

    def contig_sequences(self, sample_id: str) -> dict[str, str]:
        """``{seqid: bases}`` for one genome, from the disk cache or from ATB.

        Raises `AssemblyUnavailable` and never returns an empty dict for a genome that exists: an
        assembly with no records is a corrupt object, and saying so beats handing back a map in
        which every contig lookup will miss.
        """
        self._raise_if_known_absent(sample_id)

        cached = self._read_cached_bytes(sample_id)
        payload = cached if cached is not None else self._fetch_and_cache(sample_id)

        try:
            with gzip.open(io.BytesIO(payload), "rt") as handle:
                sequences = read_fasta_block(handle)
        except (OSError, EOFError, gzip.BadGzipFile) as error:
            # ⚠ A truncated cache entry is the one failure that survives a restart, so it is
            # DISCARDED here rather than merely reported — otherwise the genome is broken forever.
            self._discard_cached(sample_id)
            raise AssemblyUnavailable(
                f"{sample_id}'s assembly could not be decompressed ({error}); the cached copy has "
                "been discarded, so retrying will re-fetch it",
                is_retryable=True,
            ) from error

        if not sequences:
            self._discard_cached(sample_id)
            raise AssemblyUnavailable(
                f"{sample_id}'s assembly held no FASTA records, so it carries no bases",
                is_retryable=False,
            )
        return sequences

    # ----------------------------------------------------------------- absence

    def _raise_if_known_absent(self, sample_id: str) -> None:
        with self._absent_lock:
            seen_at = self._absent.get(sample_id)
            if seen_at is None:
                return
            if time.monotonic() - seen_at > self._settings.absent_ttl_seconds:
                del self._absent[sample_id]
                return
        raise AssemblyUnavailable(
            f"AllTheBacteria has no assembly for {sample_id}, so its bases cannot be read",
            is_retryable=False,
        )

    def _remember_absent(self, sample_id: str) -> None:
        with self._absent_lock:
            self._absent[sample_id] = time.monotonic()

    # ------------------------------------------------------------------- fetch

    def _fetch_and_cache(self, sample_id: str) -> bytes:
        url = self.object_url(sample_id)
        if not self._slots.acquire(timeout=SLOT_WAIT_SECONDS):
            raise AssemblyUnavailable(
                "this server is already fetching as many assemblies as it will fetch at once; "
                "the sequence is available, so please retry",
                is_retryable=True,
            )
        started = time.monotonic()
        try:
            request = urllib.request.Request(url, method="GET", headers={"Accept-Encoding": "identity"})
            with urllib.request.urlopen(request, timeout=self._settings.timeout_seconds) as response:
                payload = response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:
                self._remember_absent(sample_id)
                raise AssemblyUnavailable(
                    f"AllTheBacteria has no assembly for {sample_id}, so its bases cannot be read",
                    is_retryable=False,
                ) from error
            raise AssemblyUnavailable(
                f"AllTheBacteria returned HTTP {error.code} for {sample_id}", is_retryable=True
            ) from error
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            # ⛔ Named as a REACHABILITY failure, not as a missing genome. The two look identical on
            # a page that says only "unavailable", and only one of them resolves itself.
            raise AssemblyUnavailable(
                f"AllTheBacteria could not be reached for {sample_id} ({error})", is_retryable=True
            ) from error
        finally:
            self._slots.release()

        logger.info(
            "fetched ATB assembly %s: %d bytes in %.3f s", sample_id, len(payload), time.monotonic() - started
        )
        self._write_cached(sample_id, payload)
        return payload

    # ------------------------------------------------------------- disk cache

    def _cache_path(self, sample_id: str) -> Path | None:
        directory = self._settings.cache_directory
        if directory is None:
            return None
        # ⛔ Basename only. `sample_id` reaches here from the DATABASE, not from the URL, so it is
        # already constrained — but a cache path assembled from a value that ever becomes
        # caller-supplied is the classic traversal, and the guard costs one call.
        safe = os.path.basename(sample_id)
        if not safe or safe != sample_id:
            return None
        return directory / f"{safe}.fa.gz"

    def _read_cached_bytes(self, sample_id: str) -> bytes | None:
        path = self._cache_path(sample_id)
        if path is None:
            return None
        try:
            payload = path.read_bytes()
        except OSError:
            return None
        # Touch so eviction is least-recently-USED rather than least-recently-written.
        try:
            os.utime(path, None)
        except OSError:
            pass
        return payload or None

    def _write_cached(self, sample_id: str, payload: bytes) -> None:
        path = self._cache_path(sample_id)
        if path is None:
            return
        # ⚠ Written to a unique temporary name then renamed, because two workers fetching the same
        # genome at once would otherwise interleave into one half-written file that then fails to
        # decompress for every later reader.
        temporary = path.with_suffix(f".fa.gz.{os.getpid()}.{threading.get_ident()}.part")
        try:
            temporary.write_bytes(payload)
            temporary.replace(path)
        except OSError as error:
            logger.warning("could not cache ATB assembly %s: %s", sample_id, error)
            temporary.unlink(missing_ok=True)
            return
        self._evict_if_over_budget()

    def _discard_cached(self, sample_id: str) -> None:
        path = self._cache_path(sample_id)
        if path is not None:
            path.unlink(missing_ok=True)

    def _evict_if_over_budget(self) -> None:
        directory = self._settings.cache_directory
        if directory is None:
            return
        try:
            entries = [
                (entry.stat().st_mtime, entry.stat().st_size, entry.path)
                for entry in os.scandir(directory)
                if entry.is_file() and entry.name.endswith(".fa.gz")
            ]
        except OSError:
            return
        total = sum(size for _, size, _ in entries)
        if total <= self._settings.cache_budget_bytes:
            return
        entries.sort()  # oldest access first
        for _, size, path in entries[:EVICTION_BATCH]:
            if total <= self._settings.cache_budget_bytes:
                break
            try:
                os.unlink(path)
            except OSError:
                continue
            total -= size
