"""A similarity artifact may DECLARE two of the card's three views unmeasured — and only by declaring it.

nuna's lean producer (`build_bacatlas_artifacts`, the one that runs at 71,435 genomes) cannot compute the
weakest-link and own-neighbours views: both are gene-kNN facts, and the gene kNN is what does not scale.
Its `.meta` says `omitted=weak_intra,weak_inter,nn_own_frac` and its CSV leaves the columns out. The
ingest demanded all five and refused the artifact — found only by ingesting a lean build of the probe.

The rule these tests hold: a view is absent from a catalogue **because the artifact said so**, whole, and
from every locus. An empty column nobody declared is still the disagreement it always was.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas
import pytest

from bacatlas_backend.ingest.artifact_locator import REPRESENTATIONS, CatalogueArtifacts
from bacatlas_backend.ingest.catalogue_frames import (
    CatalogueFrameError,
    _cluster_similarity,
    _omitted_similarity_views,
    build_catalogue_frames,
)
from bacatlas_backend.ingest.derive_locus_ranking import similarity_index
from tests.conftest import NUNA_DATA_ROOT, artifacts_for, requires_artifacts

LEAN = "weak_intra,weak_inter,nn_own_frac"
BOTH = ("weak_margin", "own_fraction")

# five loci: two singletons (no within-locus pair), three measurable
WITHIN = [0.9, 0.8, None, 0.7, None]
NEAR = [0.4, 0.5, 0.3, 0.2, 0.6]
WEAK_OWN = [0.6, 0.5, None, 0.4, None]
WEAK_OTHER = [0.2, 0.3, None, 0.1, None]
OWN = [1.0, 0.9, None, 0.8, None]
NONE = [None] * 5


# ── the index: three views, or the ones the artifact measured ───────────────────────────────────
def test_a_full_artifact_ranks_all_three_views_over_one_measurable_set():
    index = similarity_index(WITHIN, NEAR, WEAK_OWN, WEAK_OTHER, OWN)
    assert index.measurable_count == 3
    for percentiles in (index.separation_percentile, index.weak_margin_percentile, index.own_fraction_percentile):
        assert [value is None for value in percentiles] == [False, False, True, False, True]


def test_two_views_declared_unmeasured_are_NULL_on_every_locus_and_separation_is_untouched():
    full = similarity_index(WITHIN, NEAR, WEAK_OWN, WEAK_OTHER, OWN)
    lean = similarity_index(WITHIN, NEAR, NONE, NONE, NONE, omitted_views=BOTH)
    assert lean.measurable_count == 3
    assert lean.separation_percentile == full.separation_percentile
    assert lean.weak_margin_percentile == NONE and lean.own_fraction_percentile == NONE


def test_an_EMPTY_view_nobody_declared_is_still_refused():
    with pytest.raises(ValueError, match="weak margin 0, own fraction 0"):
        similarity_index(WITHIN, NEAR, NONE, NONE, NONE)


def test_a_view_declared_unmeasured_that_carries_a_value_is_refused():
    with pytest.raises(ValueError, match=r"weak margin 1.*declared omitted: weak_margin, own_fraction"):
        similarity_index(WITHIN, NEAR, [0.6, *NONE[1:]], [0.2, *NONE[1:]], NONE, omitted_views=BOTH)


@pytest.mark.parametrize(
    "omitted,weak_own,weak_other,own",
    [(("own_fraction",), WEAK_OWN, WEAK_OTHER, NONE), (("weak_margin",), NONE, NONE, OWN)],
)
def test_one_view_may_be_unmeasured_while_the_other_is_ranked(omitted, weak_own, weak_other, own):
    index = similarity_index(WITHIN, NEAR, weak_own, weak_other, own, omitted_views=omitted)
    assert index.measurable_count == 3
    measured = index.weak_margin_percentile if "own_fraction" in omitted else index.own_fraction_percentile
    absent = index.own_fraction_percentile if "own_fraction" in omitted else index.weak_margin_percentile
    assert sum(value is not None for value in measured) == 3 and absent == NONE


# ── the declaration: which stems may go, and only whole views ───────────────────────────────────
def test_no_declaration_omits_nothing():
    assert _omitted_similarity_views({}, "x.meta") == ()
    assert _omitted_similarity_views({"omitted": ""}, "x.meta") == ()


def test_the_lean_producers_declaration_is_both_views():
    assert _omitted_similarity_views({"omitted": LEAN}, "x.meta") == BOTH
    assert _omitted_similarity_views({"omitted": "nn_own_frac"}, "x.meta") == ("own_fraction",)


@pytest.mark.parametrize("declared", ["within", "near,nn_own_frac", "weak_intra,weak_inter,cross_p75"])
def test_the_headline_pair_and_unknown_stems_may_not_be_omitted(declared):
    with pytest.raises(CatalogueFrameError, match="may not be omitted"):
        _omitted_similarity_views({"omitted": declared}, "x.meta")


def test_half_of_the_weakest_link_view_is_refused():
    with pytest.raises(CatalogueFrameError, match="omits weak_intra but not weak_inter"):
        _omitted_similarity_views({"omitted": "weak_intra"}, "x.meta")


# ── the reader: a lean-shaped artifact on disk ──────────────────────────────────────────────────
def _write_similarity(root: Path, *, omitted: str | None, weak_columns: bool, weak_value=None) -> CatalogueArtifacts:
    artifacts = CatalogueArtifacts(data_root=root, set_key="kp", model_label="kp_lean", run_id="r", species_key="kp")
    for representation in REPRESENTATIONS:
        csv_path = artifacts.cluster_similarity(representation)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        frame = pandas.DataFrame({"node": ["0", "1", "2"], "within_raw": [0.9, 0.8, None], "near_raw": [0.4, 0.5, 0.3]})
        if weak_columns:
            frame["weak_intra_raw"] = [weak_value, None, None]
            frame["weak_inter_raw"] = [weak_value, None, None]
            frame["nn_own_frac_raw"] = [weak_value, None, None]
        frame.to_csv(csv_path, index=False)
        meta = [f"rep={representation}", "n_loci=3", "n_measurable=2", "forms=raw", "floor_raw=0.06", "k=3"]
        if omitted is not None:
            meta.append(f"omitted={omitted}")
        artifacts.cluster_similarity_metadata(representation).write_text("\n".join(meta) + "\n")
        pandas.DataFrame(
            {"node": ["0", "1"], "form": ["raw", "raw"], "rank": [1, 1], "neighbour": ["1", "0"], "cross": [0.4, 0.5]}
        ).to_csv(artifacts.cluster_nearest(representation), index=False)
    return artifacts


POSITIONS = {"0": 0, "1": 1, "2": 2}


def test_a_lean_artifact_loads_with_the_two_views_NULL_and_says_which(tmp_path):
    blocks = _cluster_similarity(_write_similarity(tmp_path, omitted=LEAN, weak_columns=False), POSITIONS, 3, pandas)
    for representation in REPRESENTATIONS:
        block = blocks[representation]
        assert block["omitted_views"] == BOTH
        for column in ("weak_intra_raw", "weak_inter_raw", "nn_own_frac_raw"):
            assert block["loci"][column].isna().all()
        assert block["loci"]["within_raw"].notna().sum() == 2


def test_the_same_columns_missing_with_NO_declaration_are_refused_as_before(tmp_path):
    with pytest.raises(CatalogueFrameError, match="has no weak_intra_raw, weak_inter_raw, nn_own_frac_raw"):
        _cluster_similarity(_write_similarity(tmp_path, omitted=None, weak_columns=False), POSITIONS, 3, pandas)


def test_a_declared_omission_contradicted_by_values_in_the_column_is_refused(tmp_path):
    with pytest.raises(CatalogueFrameError, match="declares .* omitted, but .* carries values in"):
        _cluster_similarity(
            _write_similarity(tmp_path, omitted=LEAN, weak_columns=True, weak_value=0.5), POSITIONS, 3, pandas
        )


def test_declared_omitted_columns_that_are_present_but_EMPTY_are_the_same_omission(tmp_path):
    blocks = _cluster_similarity(_write_similarity(tmp_path, omitted=LEAN, weak_columns=True), POSITIONS, 3, pandas)
    assert all(block["omitted_views"] == BOTH for block in blocks.values())


# ── the whole catalogue: the omission changes the similarity block and NOTHING else ─────────────
@requires_artifacts
def test_a_lean_shaped_copy_of_the_published_kp_similarity_changes_no_other_frame(tmp_path):
    """The real Kp probe artifacts with the similarity CSVs rewritten the way the lean producer writes them."""
    real = artifacts_for("kp")
    proc = tmp_path / "proc"
    (proc / "analysis").mkdir(parents=True)
    for entry in (NUNA_DATA_ROOT / "proc").iterdir():
        if entry.name != "analysis":
            os.symlink(entry, proc / entry.name)
    for entry in (NUNA_DATA_ROOT / "proc" / "analysis").iterdir():
        if entry.name != "locus_browser":
            os.symlink(entry, proc / "analysis" / entry.name)
    (proc / "analysis" / "locus_browser").mkdir()
    lean = CatalogueArtifacts(
        data_root=tmp_path, set_key=real.set_key, model_label=real.model_label, run_id=real.run_id,
        species_key=real.species_key,
    )  # fmt: skip
    for representation in REPRESENTATIONS:
        frame = pandas.read_csv(real.cluster_similarity(representation), dtype={"node": str})
        frame.drop(columns=["weak_intra_raw", "weak_inter_raw", "nn_own_frac_raw"]).to_csv(
            lean.cluster_similarity(representation), index=False
        )
        lean.cluster_similarity_metadata(representation).write_text(
            real.cluster_similarity_metadata(representation).read_text().rstrip("\n") + f"\nomitted={LEAN}\n"
        )
        os.symlink(real.cluster_nearest(representation), lean.cluster_nearest(representation))
        if real.cluster_similarity_audit(representation).exists():
            os.symlink(real.cluster_similarity_audit(representation), lean.cluster_similarity_audit(representation))

    full_frames, lean_frames = build_catalogue_frames(real), build_catalogue_frames(lean)
    for representation in REPRESENTATIONS:
        assert lean_frames.similarity[representation]["omitted_views"] == BOTH
        assert full_frames.similarity[representation]["omitted_views"] == ()
        sentence = lean_frames.omitted[f"similarity_views[{representation}]"]
        assert "weakest-link view" in sentence and "own-neighbours view" in sentence
        full_index = similarity_index(
            *(list(full_frames.similarity[representation]["loci"][f"{stem}_raw"]) for stem in
              ("within", "near", "weak_intra", "weak_inter", "nn_own_frac"))
        )  # fmt: skip
        lean_index = similarity_index(
            *(list(lean_frames.similarity[representation]["loci"][f"{stem}_raw"]) for stem in
              ("within", "near", "weak_intra", "weak_inter", "nn_own_frac")),
            omitted_views=BOTH,
        )  # fmt: skip
        assert lean_index.measurable_count == full_index.measurable_count > 0
        assert lean_index.separation_percentile == full_index.separation_percentile
    assert set(lean_frames.omitted) - set(full_frames.omitted) == {
        f"similarity_views[{representation}]" for representation in REPRESENTATIONS
    }
    for name in ("loci", "annotation_entries", "arrangements", "offset_occupants", "gene_memberships", "gaps"):
        pandas.testing.assert_frame_equal(getattr(full_frames, name), getattr(lean_frames, name), check_exact=True)


# ── the estimator: the card's sentence says "median over every gene pair", and a lean build is neither ──
def test_a_full_audit_artifact_needs_no_note_because_the_cards_sentence_is_true_of_it():
    from bacatlas_backend.ingest.catalogue_frames import _similarity_estimator_note

    assert _similarity_estimator_note("esm", {"rank_by": "sim", "n_measurable": "11272"}) is None
    assert _similarity_estimator_note("esm", {}) is None


def test_a_lean_artifact_says_MEAN_and_over_how_many_of_the_genes():
    from bacatlas_backend.ingest.catalogue_frames import _similarity_estimator_note

    note = _similarity_estimator_note(
        "esm",
        {"rank_by": "mean_cross_pair_cosine", "n_genes_scored": "97094905", "n_genes_unembedded": "276886110"},
    )
    assert note == (
        "a median over every gene pair on the esm similarity card — its within-locus and nearest-locus numbers "
        "are MEANS on this build, taken over the 97,094,905 of 373,981,015 genes that carry a esm embedding"
    )
    full = _similarity_estimator_note(
        "bacformer", {"rank_by": "mean_cross_pair_cosine", "n_genes_scored": "532851", "n_genes_unembedded": "0"}
    )
    assert full.endswith("are MEANS on this build, taken over all 532,851 genes")


@requires_artifacts
def test_the_published_probe_catalogue_carries_no_estimator_note(tmp_path):
    frames = build_catalogue_frames(artifacts_for("kp"))
    assert not [key for key in frames.omitted if key.startswith("similarity_")], frames.omitted
