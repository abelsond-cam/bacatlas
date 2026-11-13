"""The window phase, reduced by genome chunk, against the unchunked build it replaced — exactly.

⛔ **The oracle is the OLD code, kept here and nowhere else.** `_unchunked_window_frames` and
`_unchunked_gene_memberships` below are `build_catalogue_frames`' window block and `_gene_memberships` as they
stood at bacatlas 73e59b0: one ~10-rows-per-gene window frame, the strand merge, two `gene_slot_vectors` passes
and a merge of every gene's vector back onto the arrangement frame. At 71,435 genomes that phase measured
0.8-1.0 TB (preflight lens "frames"), which is why it left `src`; it stays here because it is the definition the
chunked path must reproduce — values, dtypes and row order — at ANY chunk size.

Every comparison asserts its own coverage first: same columns, same row count, non-empty. A comparison that
skips what changed reports "0 differ" having checked nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bacatlas_backend.ingest import catalogue_frames
from bacatlas_backend.ingest.catalogue_frames import OFFSETS, _genome_chunks, _window_frames
from tests.conftest import artifacts_for, requires_artifacts


# ── the oracle: the unchunked window block, verbatim in substance ─────────────────────────────────
def _unchunked_gene_memberships(assign, windows, arrangement_frame, slot_of):
    from nuna.tl.locus_browser.export_payload import gene_slot_vectors

    base = assign[["sample_id", "flat_index", "node"]].copy()
    base["sample_id"] = base["sample_id"].astype(str)
    if not len(windows) or not len(arrangement_frame):
        return base.assign(arrangement_rank=None)[["node", "sample_id", "flat_index", "arrangement_rank"]]
    columns = [f"s{index}" for index in range(len(slot_of))]
    genes = gene_slot_vectors(windows, slot_of)
    ranks = arrangement_frame[["node", *columns, "rank"]].rename(columns={"rank": "arrangement_rank"})
    windowed = genes.merge(ranks, on=["node", *columns], how="left")[
        ["node", "sample_id", "flat_index", "arrangement_rank"]
    ]
    windowed["sample_id"] = windowed["sample_id"].astype(str)
    return base.merge(
        windowed[["sample_id", "flat_index", "arrangement_rank"]],
        on=["sample_id", "flat_index"],
        how="left",
    )[["node", "sample_id", "flat_index", "arrangement_rank"]]


def _unchunked_window_frames(assign, assign_coords, node_position, sample_position):
    from nuna.eval.metrics.syntology import NEAR_SHELLS, near_synteny_agreement
    from nuna.tl.locus_browser.export_payload import TOP_NEIGHBOURS, arrangements, neighbour_counts, oriented_windows

    windows = oriented_windows(assign_coords, node_position, radius=NEAR_SHELLS)
    occupants = neighbour_counts(windows, top=TOP_NEIGHBOURS)
    slot_of = {offset: index for index, offset in enumerate(OFFSETS)}
    arrangement_frame = arrangements(windows, slot_of, top=None, sample_pos=sample_position)
    synteny = near_synteny_agreement(
        windows.rename(columns={"neigh_pos": "neigh_id"}), radius=NEAR_SHELLS
    ).set_index("node")["syntenic_A5"]
    memberships = _unchunked_gene_memberships(assign, windows, arrangement_frame, slot_of)
    return occupants, arrangement_frame, synteny, memberships


# ── comparison ──────────────────────────────────────────────────────────────────────────────────
def _assert_frames_identical(got, want, what):
    assert isinstance(want, pd.DataFrame) and len(want), f"{what}: the oracle is empty — nothing compared"
    assert list(got.columns) == list(want.columns), f"{what}: columns {list(got.columns)} vs {list(want.columns)}"
    assert len(got) == len(want), f"{what}: {len(got)} rows vs {len(want)}"
    arrays = [c for c in want.columns if want[c].dtype == object and isinstance(want[c].iloc[0], np.ndarray)]
    pd.testing.assert_frame_equal(
        got.drop(columns=arrays), want.drop(columns=arrays), check_exact=True, check_dtype=True
    )
    for column in arrays:  # `assert_frame_equal` compares array VALUES but not their dtype
        for mine, theirs in zip(got[column], want[column], strict=True):
            assert type(mine) is type(theirs) and mine.dtype == theirs.dtype and np.array_equal(mine, theirs)


def _assert_window_frames_identical(got, want):
    names = ("offset_occupants", "arrangements", "syntenic_A5", "gene_memberships")
    assert len(got) == len(want) == len(names)
    for name, mine, theirs in zip(names, got, want, strict=True):
        if isinstance(theirs, pd.Series):
            assert len(theirs), f"{name}: the oracle is empty"
            pd.testing.assert_series_equal(mine, theirs, check_exact=True, check_dtype=True)
        else:
            _assert_frames_identical(mine, theirs, name)


# ── fixtures ────────────────────────────────────────────────────────────────────────────────────
def _cohort(seed: int):
    """An assignment and its coordinates holding the cases the chunked path could get wrong.

    Genomes whose rows interleave in input order; 1-gene contigs (no window), contigs shorter than the radius and
    long ones; tied starts; both strands; a conserved backbone contig, sometimes twice in one genome (ρ > 1); node
    labels whose string order is not their numeric order; assigned genes with NO coordinates (a NULL rank); a
    genome with no coordinates at all; a non-Range index on the assignment.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for genome in range(19):
        sample = f"SAMEA{(genome * 7919) % 100_003:07d}"
        flat = 0
        for copy in range(1 + int(rng.random() < 0.4)):
            for i in range(12):
                rows.append((sample, flat, 100 + copy, 1_000 * i, "+", str(20 + i)))
                flat += 1
        for contig in range(int(rng.integers(1, 5))):
            size = int(rng.choice([1, 1, 2, 3, 9, 25]))
            starts = np.sort(rng.integers(0, 30, size)) * 100
            for i in range(size):
                node = (i + contig * 3) % 12 if rng.random() < 0.8 else int(rng.integers(0, 14))
                rows.append((sample, flat, contig, int(starts[i]), "+" if rng.random() < 0.6 else "-", str(node)))
                flat += 1
    coords = pd.DataFrame(rows, columns=["sample_id", "flat_index", "contig_idx", "start", "strand", "node"])
    coords["end"] = coords["start"] + 90
    assign = coords[["sample_id", "flat_index", "node"]].copy()
    # genes assigned but never given coordinates, and one genome with none at all
    extra = pd.DataFrame({"sample_id": ["SAMEA9999999"] * 3 + [coords["sample_id"].iloc[0]],
                          "flat_index": [0, 1, 2, 10_000], "node": ["3", "4", "3", "5"]})  # fmt: skip
    assign = pd.concat([assign, extra], ignore_index=True).sample(frac=1.0, random_state=seed)
    assign.index = assign.index * 2 + 7
    coords = coords.drop(columns="node").sample(frac=1.0, random_state=seed + 1)
    assign_coords = coords.merge(assign, on=["sample_id", "flat_index"], how="inner")
    labels = sorted(assign["node"].unique(), key=lambda label: (len(label), label))
    node_position = {label: index for index, label in enumerate(reversed(labels))}
    samples = sorted(assign["sample_id"].unique())
    return assign, assign_coords, node_position, {sample: index for index, sample in enumerate(samples)}


# ── the tests ───────────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("chunk_genes", ["1", "40", "250", "100000000"])
def test_the_chunked_window_frames_are_the_unchunked_ones_exactly(monkeypatch, seed, chunk_genes):
    """`"1"` makes every genome its own chunk; `"40"` puts a genome bigger than the bound in a chunk alone."""
    monkeypatch.setenv("BACATLAS_WINDOW_CHUNK_GENES", chunk_genes)
    assign, assign_coords, node_position, sample_position = _cohort(seed)
    want = _unchunked_window_frames(assign, assign_coords, node_position, sample_position)
    got = _window_frames(assign, assign_coords, node_position, sample_position)
    _assert_window_frames_identical(got, want)


def test_the_fixture_reaches_every_case_it_claims():
    assign, assign_coords, node_position, sample_position = _cohort(0)
    *_, memberships = _unchunked_window_frames(assign, assign_coords, node_position, sample_position)
    ranks = memberships["arrangement_rank"]
    assert ranks.isna().any() and ranks.notna().any(), "no gene without a window, or none with one"
    assert assign_coords.duplicated(["sample_id", "contig_idx", "start"]).any(), "no tied start"
    sizes = np.diff(np.searchsorted(np.sort(assign_coords["sample_id"].map(sample_position)), range(20)))
    assert sizes.max() > 40, "no genome bigger than the 40-gene chunk bound"
    assert len(_genome_chunks(sizes, 250)) > 3, "the 250-gene bound does not split the cohort"


def test_every_rank_is_placed_gives_int64_as_the_join_did():
    """The join's dtype rule, kept: all placed → int64; any unplaced → float64 with NaN."""
    assign, assign_coords, node_position, sample_position = _cohort(1)
    windowed = assign_coords.groupby(["sample_id", "contig_idx"])["flat_index"].transform("size") > 1
    kept = assign_coords[windowed]
    assign = assign.merge(kept[["sample_id", "flat_index"]], on=["sample_id", "flat_index"])
    want = _unchunked_window_frames(assign, kept, node_position, sample_position)
    got = _window_frames(assign, kept, node_position, sample_position)
    assert want[3]["arrangement_rank"].dtype == np.int64, "the fixture did not reach the all-placed case"
    _assert_window_frames_identical(got, want)


def test_a_cohort_with_no_window_at_all_keeps_the_all_none_rank():
    assign = pd.DataFrame({"sample_id": ["A", "B"], "flat_index": [0, 0], "node": ["1", "2"]})
    coords = pd.DataFrame({"sample_id": ["A", "B"], "flat_index": [0, 0], "contig_idx": [0, 0],
                           "start": [1, 1], "end": [9, 9], "strand": ["+", "-"]})  # fmt: skip
    assign_coords = coords.merge(assign, on=["sample_id", "flat_index"])
    node_position, sample_position = {"1": 0, "2": 1}, {"A": 0, "B": 1}
    want = _unchunked_window_frames(assign, assign_coords, node_position, sample_position)
    got = _window_frames(assign, assign_coords, node_position, sample_position)
    for mine, theirs in zip(got[:2], want[:2], strict=True):
        assert list(mine.columns) == list(theirs.columns) and not len(mine)
        pd.testing.assert_frame_equal(mine, theirs, check_exact=True)
    pd.testing.assert_series_equal(got[2], want[2], check_exact=True)
    assert want[3]["arrangement_rank"].isna().all() and want[3]["arrangement_rank"].dtype == object
    pd.testing.assert_frame_equal(got[3], want[3], check_exact=True)


@pytest.mark.parametrize("bound", [1, 7, 40, 10_000])
def test_a_chunk_is_whole_genomes_in_order_and_never_splits_one(bound):
    sizes = np.array([0, 5, 40, 3, 3, 3, 120, 0, 0, 7, 1])
    chunks = _genome_chunks(sizes, bound)
    assert chunks[0][0] == 0 and chunks[-1][1] == len(sizes)
    assert all(stop == first for (_, stop), (first, _) in zip(chunks, chunks[1:], strict=False)), "gap or overlap"
    for first, stop in chunks:
        assert stop > first
        held = int(sizes[first:stop].sum())
        assert held <= bound or stop - first == 1 or held == int(sizes[first:stop].max()), (first, stop, held)


def test_the_chunk_size_comes_from_the_environment_and_refuses_nonsense(monkeypatch):
    monkeypatch.delenv("BACATLAS_WINDOW_CHUNK_GENES", raising=False)
    assert catalogue_frames._window_chunk_genes() == catalogue_frames.WINDOW_CHUNK_GENES
    monkeypatch.setenv("BACATLAS_WINDOW_CHUNK_GENES", "0")
    with pytest.raises(catalogue_frames.CatalogueFrameError):
        catalogue_frames._window_chunk_genes()


# ── on the real probe ───────────────────────────────────────────────────────────────────────────
@requires_artifacts
@pytest.mark.parametrize("species_key", ["kp", "ecoli"])
def test_the_probe_window_frames_are_identical_with_a_genome_per_chunk(monkeypatch, species_key):
    """Both published probe catalogues, every genome its own chunk — the most boundaries the path can have."""
    from nuna.eval.metrics.oversubscription import node_membership_profile
    from nuna.eval.prevalence import categorise
    from nuna.tl.locus_browser.export_payload import node_order

    source = catalogue_frames._load_source_frames(artifacts_for(species_key))
    assign, samples = source["assign"], source["samples"]
    profile = node_membership_profile(assign, len(samples))
    profile["band"] = categorise(profile["genomes_present"].to_numpy(), len(samples))
    node_position = {node: index for index, node in enumerate(node_order(profile)["node"].tolist())}
    sample_position = {sample: index for index, sample in enumerate(samples)}
    assign_coords = source["coords"].merge(assign, on=["sample_id", "flat_index"], how="inner")

    want = _unchunked_window_frames(assign, assign_coords, node_position, sample_position)
    monkeypatch.setenv("BACATLAS_WINDOW_CHUNK_GENES", "1")
    got = _window_frames(assign, assign_coords, node_position, sample_position)
    _assert_window_frames_identical(got, want)
