"""The phase trail goes to stderr, flushed, one line per phase — and never into a loader's report."""

from __future__ import annotations

import re

from bacatlas_backend.ingest import phase_report


def test_a_phase_line_names_the_phase_its_seconds_and_the_peak_and_goes_to_stderr_only(capsys):
    phase_report.phase("first")
    phase_report.phase("frames: gap table")
    captured = capsys.readouterr()
    assert captured.out == "", "stdout is the loaders' report and is compared between runs"
    lines = captured.err.splitlines()
    assert len(lines) == 2
    pattern = r"  \[phase\] frames: gap table\s+\d+\.\ds(  rss=\s*\d+\.\d\dG)?  peak=\s*\d+\.\d\dG"
    assert re.fullmatch(pattern, lines[1]), lines[1]


def test_the_peak_is_in_bytes_on_both_platforms(monkeypatch):
    class Usage:
        ru_maxrss = 3 * 2**20

    monkeypatch.setattr(phase_report.resource, "getrusage", lambda _who: Usage)
    monkeypatch.setattr(phase_report.sys, "platform", "linux")
    assert phase_report._peak_bytes() == 3 * 2**30, "Linux reports kilobytes"
    monkeypatch.setattr(phase_report.sys, "platform", "darwin")
    assert phase_report._peak_bytes() == 3 * 2**20, "macOS reports bytes"
