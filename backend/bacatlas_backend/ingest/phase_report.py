"""One line per phase of a long ingest stage: its seconds, and the process's memory when it ended.

⛔ **A stage that runs for hours in silence can only be sized after it dies.** The pangenome stage of
a 71,435-genome cohort builds every frame of the catalogue in one process on one node; when the
kernel kills it for memory there is no traceback, and without a trail the only facts left are the
exit code and the wall clock. A line per phase turns that into a located measurement — which phase,
how long, how far the peak had climbed — and for a run that succeeds it is the per-phase record the
next cohort is sized from.

⚠ **To stderr, flushed, and never part of a report.** The loaders' stdout is their report and is
compared between runs; a timing line there would make two identical loads differ. stderr is where a
scheduler's log wants it, and `flush=True` is what makes the last line before a kill the true one.

⭐ `peak` is the kernel's high-water mark for this process (`ru_maxrss`) — the number a scheduler
checks its memory limit against, monotone, never falling. `rss` is what is resident now, read from
`/proc/self/statm` where that exists; it falls when a frame is released, which is how a release is
seen to have worked.
"""

from __future__ import annotations

import os
import resource
import sys
import time

_last: float | None = None


def _peak_bytes() -> int:
    """`ru_maxrss` in bytes: Linux reports kilobytes, macOS bytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak) if sys.platform == "darwin" else int(peak) * 1024


def _resident_bytes() -> int | None:
    """Resident set now, or `None` where the platform has no `/proc`."""
    try:
        with open("/proc/self/statm") as handle:
            return int(handle.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError, IndexError):
        return None


def phase(name: str) -> None:
    """Say that ``name`` has just finished: seconds since the previous mark, memory now and at peak."""
    global _last  # noqa: PLW0603 — one clock per process is the point
    now = time.monotonic()
    elapsed = 0.0 if _last is None else now - _last
    _last = now
    resident = _resident_bytes()
    now_text = "" if resident is None else f"  rss={resident / 2**30:7.2f}G"
    print(
        f"  [phase] {name:<34} {elapsed:9.1f}s{now_text}  peak={_peak_bytes() / 2**30:7.2f}G",
        file=sys.stderr,
        flush=True,
    )
