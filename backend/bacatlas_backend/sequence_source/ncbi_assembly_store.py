"""One RefSeq/GenBank assembly, fetched from NCBI Datasets on demand and cached on disk.

⭐ **Why this exists.** 4,806 of the 71,435-genome Klebsiella cohort are RefSeq/GenBank long-read
assemblies (`genome.sample_id_kind = assembly_stem`, `genome.bare_accession = GCF_000512165.1`), and their
bases are not on the server and never will be — David, 2026-10-30: *"If the anchoring option to an NCBI
genome is chosen, then the assembly will be streamed so that the sequence tab can be seen. That is all."*
So this is the fourth station of `ContigSequenceProvider`'s chain — held GFF → cached → ATB → NCBI — and it
is **serving-only**: the ingest never reads it.

⚠ **Measured against the live API, 2026-11-13, one request** (UK domestic line):
``GET /datasets/v2/genome/accession/GCF_000512165.1/download?include_annotation_type=GENOME_FASTA``
answered **200 ``application/zip``, 1,699,070 bytes in 1.08 s**. The zip held ``README.md``, ``md5sum.txt``,
``ncbi_dataset/data/dataset_catalog.json``, ``ncbi_dataset/data/assembly_data_report.jsonl`` and exactly ONE
``ncbi_dataset/data/GCF_000512165.1/GCF_000512165.1_ASM51216v1_genomic.fna`` (5,812,292 bytes), whose
deflines read ``>NZ_CP003999.1 Klebsiella pneumoniae subsp. pneumoniae Kp13 chromosome, complete genome``.
So the first token is the **RefSeq-prefixed** accession — the ``NZ_`` namespace split the record describes
is exactly what this store hands back, and `contig_namespace.resolve_contig_names` is what reconciles it.

⛔ **An unknown accession is NOT a 404.** Measured in the same session: ``GCF_000000000.1`` answered **200
``application/zip``** with a **971-byte body that is not a valid zip at all** — a streamed archive abandoned
after the README's local header, with no central directory. So "NCBI has no such assembly" is read off the
BODY: a 200 carrying no ``ncbi_dataset/data/<accession>/*_genomic.fna`` member, valid zip or not, is the
absence, and it is that case the negative cache remembers. A 404 is handled the same way as a second,
defensive signal for the day the API changes. ⚠ A body that is not a zip but DOES name the member is a
truncated download, which is transient — the two are told apart by looking for the member's path in the
bytes, because a zip writes each member's name in cleartext ahead of its data.

⛔ **A third party, like ATB; a cache with a fallback, never a dependency.** A genome whose bases are held
or cached is never fetched, and a fetch that fails is a named 503 rather than a wrong answer. NCBI meters
anonymous clients at **3 requests/s (10 with an API key)** and answers an overrun with 429, so this store
spaces its own requests with a token bucket, per process, rather than discovering the limit from the
error. ⚠ Per PROCESS: gunicorn's N workers each hold a bucket, so the instance's ceiling is N × the rate.
With 4 workers at the default 3/s that is 12/s on paper — but a worker only spends a token on a COLD open
of an NCBI-keyed genome, which at human concurrency is rare, and the semaphore bounds it further.

⚠ **No new dependency.** ``urllib.request`` and ``zipfile`` from the standard library, as the ATB store
uses ``urllib.request`` and ``gzip``.
"""

from __future__ import annotations

import io
import logging
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from bacatlas_backend.gff.fasta_reader import read_fasta_block
from bacatlas_backend.sequence_source.atb_assembly_store import AssemblyUnavailable

logger = logging.getLogger(__name__)

#: NCBI Datasets v2 REST. The per-accession download returns a zip "data package"; the member we want
#: is selected by `include_annotation_type=GENOME_FASTA`, which keeps the package to the genomic FASTA
#: plus its catalogue — ~1.7 MB for a 5.8 MB Klebsiella genome, against the full annotation package.
DEFAULT_NCBI_BASE_URL: Final = "https://api.ncbi.nlm.nih.gov/datasets/v2"

#: ⚠ Three times the ATB timeout, because the one measured fetch took 1.08 s against ATB's p95 of
#: 0.330 s — NCBI assembles the package on request where S3 serves a stored object. Too tight turns a
#: slow-but-fine fetch into a 503; too loose only holds a worker the semaphore has already bounded.
DEFAULT_TIMEOUT_SECONDS: Final = 30.0

#: The worker bound, for the same reason the ATB store has one: gunicorn serves this synchronously, and
#: N concurrent cold opens hold N workers for the whole fetch. See `atb_assembly_store`.
DEFAULT_MAX_CONCURRENT_FETCHES: Final = 4

#: How long to wait for a fetch slot before saying "busy". Short on purpose.
SLOT_WAIT_SECONDS: Final = 2.0

#: ⛔ The negative cache. "NCBI has no assembly under this accession" is a fact that will not change
#: within a session, and without this every page view of such a genome pays a round trip (and a rate
#: token) to learn it again.
DEFAULT_ABSENT_TTL_SECONDS: Final = 600.0

#: ⛔ A SEPARATE budget from ATB's, in a separate directory, so the two caches never evict each other.
#: 10 GiB holds the whole 4,806-genome slice at the measured 1.7 MB mean with room to spare — which is
#: a bound, not a target: the working set of any reading session is a handful of genomes.
DEFAULT_CACHE_BUDGET_BYTES: Final = 10 * 1024**3

#: NCBI's published anonymous limit. 10/s with an API key — set `requests_per_second` to match when
#: one is configured; the store does not raise the rate on its own because the key may be wrong.
DEFAULT_REQUESTS_PER_SECOND: Final = 3.0

#: Evicting one file per fetch would thrash at the boundary; a batch amortises the ``scandir``.
EVICTION_BATCH: Final = 64

#: Where the data package keeps its files. The accession's own directory sits directly under it.
_DATA_PREFIX: Final = "ncbi_dataset/data/"

#: The cache file suffix — the fetched object is kept as-is, exactly as the ATB store keeps its gzip.
_CACHE_SUFFIX: Final = ".zip"


class TokenBucket:
    """A token bucket: ``rate_per_second`` tokens refill continuously, up to ``capacity``; one per request.

    ⚠ ``clock`` and ``sleep`` are injectable so a test can prove the spacing without waiting for it. A
    caller that finds the bucket empty RESERVES its token (the count goes negative) and sleeps for the
    deficit outside the lock, so concurrent callers queue at the right intervals instead of all waking
    at once and racing for the same refill.
    """

    def __init__(
        self,
        rate_per_second: float,
        *,
        capacity: float | None = None,
        clock: Callable[[], float] = time.monotonic,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        if rate_per_second <= 0:
            raise ValueError(f"requests_per_second must be positive, not {rate_per_second!r}")
        self._rate = float(rate_per_second)
        #: One second's allowance, so a burst after idleness never exceeds what the sustained rate
        #: would have permitted in that second anyway.
        self._capacity = float(capacity) if capacity is not None else max(1.0, self._rate)
        self._tokens = self._capacity
        self._clock = clock
        self._sleep = sleep
        self._updated = clock()
        self._lock = threading.Lock()

    def acquire(self) -> float:
        """Take one token, sleeping until one is due. Returns the seconds waited."""
        with self._lock:
            now = self._clock()
            self._tokens = min(self._capacity, self._tokens + (now - self._updated) * self._rate)
            self._updated = now
            self._tokens -= 1.0
            if self._tokens >= 0.0:
                return 0.0
            wait = -self._tokens / self._rate
        self._sleep(wait)
        return wait


@dataclass(frozen=True)
class NcbiStoreSettings:
    """Everything the store reads from the outside, resolved once at start-up."""

    base_url: str = DEFAULT_NCBI_BASE_URL
    cache_directory: Path | None = None
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS
    max_concurrent_fetches: int = DEFAULT_MAX_CONCURRENT_FETCHES
    absent_ttl_seconds: float = DEFAULT_ABSENT_TTL_SECONDS
    cache_budget_bytes: int = DEFAULT_CACHE_BUDGET_BYTES
    requests_per_second: float = DEFAULT_REQUESTS_PER_SECOND
    #: Sent as the ``api-key`` header when set. Raises NCBI's limit to 10/s — raise the knob above to
    #: use it; the store will not assume the key is valid.
    api_key: str | None = None


class NcbiAssemblyStore:
    """Fetch, cache and parse NCBI assemblies by accession. Thread-safe; one instance per process.

    ⚠ **The disk cache is keyed on the bare accession and holds the ZIP as fetched**, not the parsed
    dict — the same discipline as the ATB store's gzip: the compressed object is smaller on disk and a
    smaller thing to get wrong, and the parse (one member, one FASTA read) is tens of milliseconds. The
    parsed form is cached one layer up, by `ContigSequenceProvider`, where every source shares it.

    ⛔ Returns contigs keyed by the FIRST TOKEN of each defline, which for a RefSeq assembly is the
    ``NZ_``-prefixed accession. This store does not know which names the caller's coordinates use; the
    provider does, and resolves them with `contig_namespace.resolve_contig_names`.
    """

    def __init__(self, settings: NcbiStoreSettings | None = None, *, rate_limiter: TokenBucket | None = None) -> None:
        self._settings = settings or NcbiStoreSettings()
        self._slots = threading.BoundedSemaphore(max(1, self._settings.max_concurrent_fetches))
        self._bucket = rate_limiter or TokenBucket(self._settings.requests_per_second)
        self._absent: dict[str, float] = {}
        self._absent_lock = threading.Lock()
        if self._settings.cache_directory is not None:
            self._settings.cache_directory.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ public

    def object_url(self, accession: str) -> str:
        """The download this accession would be at. Separate so a test can assert it without a GET."""
        return (
            f"{self._settings.base_url.rstrip('/')}/genome/accession/"
            f"{urllib.parse.quote(accession, safe='')}/download?include_annotation_type=GENOME_FASTA"
        )

    def contig_sequences(self, accession: str) -> dict[str, str]:
        """``{header token: bases}`` for one assembly, from the disk cache or from NCBI.

        Raises `AssemblyUnavailable` and never returns an empty dict for an assembly that exists: a
        package with no records is a corrupt object, and saying so beats handing back a map in which
        every contig lookup will miss.
        """
        self._raise_if_known_absent(accession)

        cached = self._read_cached_bytes(accession)
        if cached is not None:
            sequences = self._parse_cached(accession, cached)
        else:
            payload = self._fetch(accession)
            sequences = self._parse_fetched(accession, payload)
            # ⛔ Written only once the member has been found and read: a package that carries no
            # assembly is NCBI's "not found", and caching it would turn a negative answer into a
            # permanent one that survives the negative cache's TTL.
            self._write_cached(accession, payload)

        if not sequences:
            self._discard_cached(accession)
            raise AssemblyUnavailable(
                f"{accession}'s assembly held no FASTA records, so it carries no bases", is_retryable=False
            )
        return sequences

    # ----------------------------------------------------------------- absence

    def _raise_if_known_absent(self, accession: str) -> None:
        with self._absent_lock:
            seen_at = self._absent.get(accession)
            if seen_at is None:
                return
            if time.monotonic() - seen_at > self._settings.absent_ttl_seconds:
                del self._absent[accession]
                return
        raise self._absent_error(accession)

    def _remember_absent(self, accession: str) -> None:
        with self._absent_lock:
            self._absent[accession] = time.monotonic()

    @staticmethod
    def _absent_error(accession: str) -> AssemblyUnavailable:
        return AssemblyUnavailable(
            f"NCBI Datasets has no assembly for {accession}, so its bases cannot be read", is_retryable=False
        )

    # ------------------------------------------------------------------- fetch

    def _fetch(self, accession: str) -> bytes:
        url = self.object_url(accession)
        if not self._slots.acquire(timeout=SLOT_WAIT_SECONDS):
            raise AssemblyUnavailable(
                "this server is already fetching as many assemblies from NCBI as it will fetch at once; "
                "the sequence is available, so please retry",
                is_retryable=True,
            )
        started = time.monotonic()
        try:
            # ⚠ The token is taken INSIDE the slot, so at most `max_concurrent_fetches` requests can ever
            # be waiting on the bucket — the longest wait is slots ÷ rate, ~1.3 s at the defaults — and a
            # burst of cold opens is told "busy" by the semaphore rather than queued without bound.
            waited = self._bucket.acquire()
            headers = {"Accept": "application/zip", "Accept-Encoding": "identity"}
            if self._settings.api_key:
                headers["api-key"] = self._settings.api_key
            request = urllib.request.Request(url, method="GET", headers=headers)
            with urllib.request.urlopen(request, timeout=self._settings.timeout_seconds) as response:
                payload = response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:
                self._remember_absent(accession)
                raise self._absent_error(accession) from error
            if error.code == 429:
                raise AssemblyUnavailable(
                    f"NCBI rate-limited this server while fetching {accession} (HTTP 429); please retry",
                    is_retryable=True,
                ) from error
            if error.code == 400:
                # ⚠ Malformed, not missing: NCBI rejected the accession's SHAPE, which no retry fixes —
                # but it is not negative-cached, because it says nothing about the assembly's existence.
                raise AssemblyUnavailable(
                    f"NCBI Datasets rejected {accession!r} as an accession (HTTP 400)", is_retryable=False
                ) from error
            raise AssemblyUnavailable(
                f"NCBI Datasets returned HTTP {error.code} for {accession}", is_retryable=True
            ) from error
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            # ⛔ Named as a REACHABILITY failure, not as a missing genome.
            raise AssemblyUnavailable(
                f"NCBI Datasets could not be reached for {accession} ({error})", is_retryable=True
            ) from error
        finally:
            self._slots.release()

        logger.info(
            "fetched NCBI assembly %s: %d bytes in %.3f s (rate-limit wait %.3f s)",
            accession,
            len(payload),
            time.monotonic() - started,
            waited,
        )
        return payload

    # ------------------------------------------------------------------- parse

    def _parse_fetched(self, accession: str, payload: bytes) -> dict[str, str]:
        """A fresh download → contigs, or the named failure. Absence is REMEMBERED here, not above."""
        try:
            sequences = self._genomic_fasta(accession, payload)
        except zipfile.BadZipFile as error:
            if self._names_the_member(accession, payload):
                raise AssemblyUnavailable(
                    f"the download of {accession} from NCBI was truncated ({error}); nothing was cached, "
                    "so retrying will fetch it again",
                    is_retryable=True,
                ) from error
            # ⛔ The measured shape of "no such accession": a 200 whose body is an abandoned stream with
            # no member for the accession anywhere in it.
            self._remember_absent(accession)
            raise self._absent_error(accession) from error
        if sequences is None:
            self._remember_absent(accession)
            raise self._absent_error(accession)
        return sequences

    def _parse_cached(self, accession: str, payload: bytes) -> dict[str, str]:
        """A cached package → contigs. ⚠ A cache entry that will not read is DISCARDED, not reported."""
        try:
            sequences = self._genomic_fasta(accession, payload)
        except zipfile.BadZipFile as error:
            # The one failure that survives a restart: a half-written entry would otherwise fail for
            # every later reader, so it is removed and the retry the message promises re-fetches.
            self._discard_cached(accession)
            raise AssemblyUnavailable(
                f"{accession}'s cached package could not be read ({error}); the cached copy has been "
                "discarded, so retrying will re-fetch it",
                is_retryable=True,
            ) from error
        if sequences is None:
            self._discard_cached(accession)
            raise AssemblyUnavailable(
                f"{accession}'s cached package carried no assembly; the cached copy has been discarded, "
                "so retrying will re-fetch it",
                is_retryable=True,
            )
        return sequences

    def _genomic_fasta(self, accession: str, payload: bytes) -> dict[str, str] | None:
        """The package's one genomic FASTA for ``accession``, parsed — or ``None`` if it holds none.

        Raises ``zipfile.BadZipFile`` for a body that is not a readable zip (a bad CRC on the member
        counts: ``zipfile`` raises it on read, which is how a truncated member is caught).
        """
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            member = self._member_for(accession, archive.namelist())
            if member is None:
                return None
            with archive.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle:
                # ⛔ The same reader as the held GFF and ATB: first token of the defline, upper-cased.
                return read_fasta_block(handle)

    @staticmethod
    def _member_for(accession: str, names: list[str]) -> str | None:
        """The ``*_genomic.fna`` under THIS accession's directory — exactly one, or none, or a refusal.

        ⛔ Scoped to ``ncbi_dataset/data/<accession>/``, never to any ``_genomic.fna`` in the package:
        bases filed under another accession are plausible wrong bases, which is the failure nothing
        downstream can see. ``cds_from_genomic.fna`` / ``rna_from_genomic.fna`` are excluded by suffix.
        """
        prefix = f"{_DATA_PREFIX}{accession}/"
        candidates = [
            name
            for name in names
            if name.startswith(prefix) and name.endswith("_genomic.fna") and not name.endswith("from_genomic.fna")
        ]
        if not candidates:
            return None
        if len(candidates) > 1:
            raise AssemblyUnavailable(
                f"NCBI's package for {accession} holds {len(candidates)} genomic FASTA members "
                f"({', '.join(sorted(candidates))}); refusing to guess which is the assembly",
                is_retryable=False,
            )
        return candidates[0]

    @staticmethod
    def _names_the_member(accession: str, payload: bytes) -> bool:
        """Whether the bytes carry a path under this accession's data directory — zip names are cleartext."""
        return f"{_DATA_PREFIX}{accession}/".encode() in payload

    # ------------------------------------------------------------- disk cache

    def _cache_path(self, accession: str) -> Path | None:
        directory = self._settings.cache_directory
        if directory is None:
            return None
        # ⛔ Basename only — see the ATB store. The accession reaches here from the database, but a
        # cache path assembled from a value that ever becomes caller-supplied is the classic traversal.
        safe = os.path.basename(accession)
        if not safe or safe != accession:
            return None
        return directory / f"{safe}{_CACHE_SUFFIX}"

    def _read_cached_bytes(self, accession: str) -> bytes | None:
        path = self._cache_path(accession)
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

    def _write_cached(self, accession: str, payload: bytes) -> None:
        path = self._cache_path(accession)
        if path is None:
            return
        # ⚠ Unique temporary name then rename, so two workers fetching one genome cannot interleave.
        temporary = path.with_suffix(f"{_CACHE_SUFFIX}.{os.getpid()}.{threading.get_ident()}.part")
        try:
            temporary.write_bytes(payload)
            temporary.replace(path)
        except OSError as error:
            logger.warning("could not cache NCBI assembly %s: %s", accession, error)
            temporary.unlink(missing_ok=True)
            return
        self._evict_if_over_budget()

    def _discard_cached(self, accession: str) -> None:
        path = self._cache_path(accession)
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
                if entry.is_file() and entry.name.endswith(_CACHE_SUFFIX)
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
