"""`Configuration` → a wired `ContigSequenceProvider`, and the key it is stored under.

⛔ **Here rather than in `application_factory` because both the factory and the species blueprint
need it**, and putting it in the factory made the blueprint import the module that imports the
blueprint — a circular import that fails at start-up. This module imports only configuration and
the sequence sources, so nothing can import its way back round.
"""

from __future__ import annotations

import logging

from bacatlas_backend.configuration import Configuration
from bacatlas_backend.sequence_source.atb_assembly_store import AtbAssemblyStore, AtbStoreSettings
from bacatlas_backend.sequence_source.contig_sequence_provider import ContigSequenceProvider
from bacatlas_backend.sequence_source.ncbi_assembly_store import NcbiAssemblyStore, NcbiStoreSettings

logger = logging.getLogger(__name__)

#: ⭐ Where the Sequence tab's resolver lives on `application.extensions`. ⛔ **One per application,
#: not one per request** — it holds the parsed-genome cache that used to be a module-level
#: `lru_cache`, so building it per request would silently turn every tab open into a fresh parse
#: (and, with streaming on, a fresh 0.24 s fetch) while every test still passed.
SEQUENCE_PROVIDER_KEY = "bacatlas_sequence_provider"


def build_sequence_provider(configuration: Configuration) -> ContigSequenceProvider:
    """The held → cached → ATB | NCBI chain, wired from configuration.

    ⚠ **A provider is always built, even with no GFF root and no streaming.** It then raises a named
    failure on every request, which the endpoint already turns into a 503 naming the missing
    setting. Returning `None` here instead would push a second "is sequence configured at all?"
    branch into the endpoint, and that branch is how "not configured" came to be indistinguishable
    from "the file is missing" in the first place.
    """
    store: AtbAssemblyStore | None = None
    if configuration.atb_streaming_enabled:
        defaults = AtbStoreSettings()
        store = AtbAssemblyStore(
            AtbStoreSettings(
                base_url=configuration.atb_base_url or defaults.base_url,
                cache_directory=configuration.atb_cache_directory,
                timeout_seconds=configuration.atb_timeout_seconds or defaults.timeout_seconds,
                max_concurrent_fetches=(
                    configuration.atb_max_concurrent_fetches or defaults.max_concurrent_fetches
                ),
                cache_budget_bytes=(configuration.atb_cache_budget_bytes or defaults.cache_budget_bytes),
            )
        )
        logger.info(
            "assembly streaming ENABLED (cache: %s)",
            configuration.atb_cache_directory or "disabled — every miss re-fetches",
        )
    ncbi_store: NcbiAssemblyStore | None = None
    if configuration.ncbi_streaming_enabled:
        ncbi_defaults = NcbiStoreSettings()
        # ⚠ `requests_per_second` ≤ 0 raises here (ValueError from the token bucket) — at start-up,
        # which is where a knob that would otherwise stall every fetch forever should fail.
        ncbi_store = NcbiAssemblyStore(
            NcbiStoreSettings(
                base_url=configuration.ncbi_base_url or ncbi_defaults.base_url,
                cache_directory=configuration.ncbi_cache_directory,
                timeout_seconds=configuration.ncbi_timeout_seconds or ncbi_defaults.timeout_seconds,
                max_concurrent_fetches=(
                    configuration.ncbi_max_concurrent_fetches or ncbi_defaults.max_concurrent_fetches
                ),
                cache_budget_bytes=(configuration.ncbi_cache_budget_bytes or ncbi_defaults.cache_budget_bytes),
                requests_per_second=(
                    configuration.ncbi_requests_per_second
                    if configuration.ncbi_requests_per_second is not None
                    else ncbi_defaults.requests_per_second
                ),
                api_key=configuration.ncbi_api_key,
            )
        )
        logger.info(
            "NCBI assembly streaming ENABLED (cache: %s; %s/s; api key %s)",
            configuration.ncbi_cache_directory or "disabled — every miss re-fetches",
            configuration.ncbi_requests_per_second or ncbi_defaults.requests_per_second,
            "set" if configuration.ncbi_api_key else "not set",
        )
    return ContigSequenceProvider(
        gff_root=configuration.artifact_roots.get("gff"),
        atb_store=store,
        ncbi_store=ncbi_store,
    )
