"""Environment-driven settings — the only place a path or a URL is read from the outside.

Two deployment profiles, and the difference between them is entirely *where the files are*:
the Mac has the cluster artifacts under ``~/developer/nuna/data``, the med school server has
whatever its sysadmin gave us. Nothing else in the package reads ``os.environ``, so a wrong
path is one exception at start-up rather than a 404 on one endpoint three weeks later.

⛔ ``artifact_roots`` is a MAP, not a single path. Every ``artifact_pointer`` row stores a root
KEY plus a relative path, so moving the store is a config change and not a database migration.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

#: Root keys an artifact_pointer row may name. Adding one is a schema decision, not a config
#: convenience, so the set is closed here and validated at start-up.
ARTIFACT_ROOT_KEYS = ("gff", "assemblies", "embeddings", "analysis")


#: ⛔ The same truthy words `BACATLAS_SQL_ECHO` already accepted, in one place so a new flag cannot
#: quietly recognise a different set — `BACATLAS_ATB_STREAMING=on` silently meaning *off* is exactly
#: the kind of configuration bug that reads as "streaming does not work".
TRUTHY_FLAG_VALUES = frozenset({"1", "true", "yes", "on"})


def _flag(name: str) -> bool:
    """A boolean environment flag. Absent is False."""
    return os.environ.get(name, "").strip().lower() in TRUTHY_FLAG_VALUES


def _optional_int(name: str) -> int | None:
    """An integer override, or ``None``.

    ⚠ Raises on a value that is set but unparseable, rather than falling back to the default — a
    typo'd tuning knob must not look like one that was never set.
    """
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return None
    try:
        return int(raw)
    except ValueError as error:
        raise RuntimeError(f"{name}={raw!r} is not an integer") from error


def _optional_float(name: str) -> float | None:
    """As `_optional_int`, for a duration."""
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return None
    try:
        return float(raw)
    except ValueError as error:
        raise RuntimeError(f"{name}={raw!r} is not a number") from error


def _require(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(
            f"{name} is not set. The API refuses to start without it rather than failing on the "
            f"first request that needs it — see bacatlas_backend/configuration.py."
        )
    return value


@dataclass(frozen=True)
class Configuration:
    """Everything the application needs from its environment, resolved once."""

    database_url: str
    artifact_roots: dict[str, Path] = field(default_factory=dict)
    deployment_profile: str = "development"
    sql_echo: bool = False

    #: ⭐ **Fetch a genome's assembly from AllTheBacteria when we do not hold its GFF.** Off by
    #: default, so an existing deployment behaves exactly as it did until it is switched on. ⚠ It is
    #: additive even when on: the held GFF always wins, so this can only change the answer for a
    #: genome that previously returned `sequence_unavailable`.
    atb_streaming_enabled: bool = False
    #: Where fetched assemblies are cached. ⛔ **Never `$HOME` and never the repo** — on the instance
    #: this is a directory under `/srv/bacatlas/`, on the budget below. `None` disables the disk
    #: cache entirely, which is correct for a test and wrong for a server.
    atb_cache_directory: Path | None = None
    #: Overridable only so a test can point at a local stub; there is no second real value.
    atb_base_url: str | None = None
    atb_timeout_seconds: float | None = None
    atb_max_concurrent_fetches: int | None = None
    atb_cache_budget_bytes: int | None = None

    @classmethod
    def from_environment(cls) -> Configuration:
        """Build from ``BACATLAS_*`` variables. Raises on anything missing or unreadable."""
        roots: dict[str, Path] = {}
        for key in ARTIFACT_ROOT_KEYS:
            raw = os.environ.get(f"BACATLAS_ROOT_{key.upper()}")
            if raw:
                roots[key] = Path(raw).expanduser()
        cache_directory = os.environ.get("BACATLAS_ATB_CACHE_DIR")
        return cls(
            database_url=_require("BACATLAS_DATABASE_URL"),
            artifact_roots=roots,
            deployment_profile=os.environ.get("BACATLAS_PROFILE", "development"),
            sql_echo=_flag("BACATLAS_SQL_ECHO"),
            atb_streaming_enabled=_flag("BACATLAS_ATB_STREAMING"),
            atb_cache_directory=Path(cache_directory).expanduser() if cache_directory else None,
            atb_base_url=os.environ.get("BACATLAS_ATB_BASE_URL") or None,
            atb_timeout_seconds=_optional_float("BACATLAS_ATB_TIMEOUT_SECONDS"),
            atb_max_concurrent_fetches=_optional_int("BACATLAS_ATB_MAX_CONCURRENT_FETCHES"),
            atb_cache_budget_bytes=_optional_int("BACATLAS_ATB_CACHE_BUDGET_BYTES"),
        )

    def resolve_artifact(self, root_key: str, relative_path: str) -> Path:
        """``(root_key, relative_path)`` → an absolute path, or a named failure.

        ⚠ The error names the root key AND the configured roots, because the failure this
        actually produces in practice is a profile that was never given the root at all —
        which is indistinguishable from a missing file unless the message says so.
        """
        if root_key not in self.artifact_roots:
            raise KeyError(
                f"artifact root {root_key!r} is not configured for profile "
                f"{self.deployment_profile!r} (have: {sorted(self.artifact_roots) or 'none'}). "
                f"Set BACATLAS_ROOT_{root_key.upper()}."
            )
        return self.artifact_roots[root_key] / relative_path
