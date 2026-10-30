"""The frame-peak instrument's own failure modes — measured nothing must never print as measured zero.

⛔ Its B/row figure comes from `RssAnon`, which `/proc/self/status` supplies on Linux and nothing supplies
on macOS. So on a laptop this instrument CAN run and CANNOT measure, and the one unacceptable outcome is a
plausible number. These tests pin that, and the roster slice's coverage refusal, both of which are pure and
testable anywhere; the measurement itself is only meaningful on CSD3.
"""

from __future__ import annotations

import re

import pandas as pd
import pytest

from bacatlas_backend.instruments.frame_peak_census import Census, Mark, _roster_filter, report


def _census(**anon):
    c = Census(categorical=False, n_genomes=2, n_genes=10, n_loci=4)
    for label, value in anon.items():
        c.marks.append(Mark(label, hwm_gb=1.0, anon_gb=value, rows=100))
    return c


def _source(samples, rows):
    frame = pd.DataFrame({
        "sample_id": [s for s, _ in rows],
        "flat_index": [f for _, f in rows],
        "node": ["n0"] * len(rows),
    })  # fmt: skip
    return {
        "samples": list(samples),
        "assign": frame,
        "products": frame.drop(columns="node").assign(product="p"),
        "coords": frame.drop(columns="node").assign(start=1, end=2, contig_idx=0, strand="+", rank_on_contig=0),
        "noncoding": None,
        "dbxrefs": None,
        "pfam": None,
        "n_genes": len(frame),
    }


def test_the_roster_slice_filters_every_per_gene_frame_not_some():
    """⛔ A coords of 1,000 genomes merged against an assign of 71,435 still yields a 1,000-genome window
    frame — so a partial filter would report a peak that silently carries the big frame too."""
    source = _source(["A", "B", "C"], [("A", 0), ("A", 1), ("B", 0), ("C", 0)])
    got = _roster_filter(source, _roster(["A", "B"]))
    assert got["samples"] == ["A", "B"]
    for key in ("assign", "products", "coords"):
        assert set(got[key]["sample_id"]) == {"A", "B"}, f"{key} was not sliced"
    assert got["n_genes"] == 3


def test_the_roster_slice_refuses_a_genome_this_catalogue_does_not_have():
    source = _source(["A", "B"], [("A", 0), ("B", 0)])
    with pytest.raises(SystemExit, match="would measure a different cohort"):
        _roster_filter(source, _roster(["A", "ZZ"]))


def test_no_roster_leaves_the_source_untouched():
    source = _source(["A", "B"], [("A", 0), ("B", 0)])
    assert _roster_filter(source, None) is source


def test_a_missing_rssanon_reports_not_measurable_rather_than_a_number():
    """⛔ The macOS path. A B/row printed from a None delta is the exact shape of a measured-zero lie."""
    c = _census(sources_loaded=None, annotation_frame=None, assign_coords=None, oriented_windows=None)
    text = report(c)
    assert not re.search(r"[\d.]+ B/row", text), "a B/row NUMBER from a None delta is the measured-zero lie"
    assert text.count("RssAnon unavailable") == 2, "both frames must say so, not just the first"
    assert "VmHWM over the whole run" in text, "the peak is still reported — it comes from a different field"


def test_bytes_per_row_is_computed_from_the_heap_delta_over_that_frames_own_rows():
    c = Census(categorical=True, n_genomes=2, n_genes=10, n_loci=4)
    c.marks += [
        Mark("sources_loaded", hwm_gb=1.0, anon_gb=1.0, rows=10),
        Mark("annotation_frame", hwm_gb=2.0, anon_gb=2.0, rows=10),
        Mark("assign_coords", hwm_gb=2.0, anon_gb=2.0, rows=10),
        # +1 GiB of heap over 10,000,000 rows = 107.4 B/row
        Mark("oriented_windows", hwm_gb=3.0, anon_gb=3.0, rows=10_000_000),
        Mark("windows_freed", hwm_gb=3.0, anon_gb=2.0),
    ]
    text = report(c)
    assert "107.4 B/row" in text
    assert "returned 1.000 GB" in text, "freeing the frame must be reported, or a leak reads as a cost"
    assert "VmHWM over the whole run: 3.000 GB" in text


def _roster(names):
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "roster.txt"
    path.write_text("\n".join(names) + "\n")
    return path
