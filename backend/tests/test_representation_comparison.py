"""ESM against Bacformer on the same pairs — the two guards that keep the comparison honest.

⛔ Both failures this pins are silent. A column labelled `ESM` fed a Bacformer map inverts the whole
comparison and every number stays in range; and comparing each map's AUC over its *own* shortlist
rewards whichever shortlist carries more easy negatives, which is a statement about candidate
generation rather than about the embeddings.

Runs without a database: both guards fire before any query.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


def _load(name: str):
    path = Path(__file__).resolve().parent.parent / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cmp = _load("compare_representations_cross_species")

COLUMNS = [
    "kp_node",
    "ecoli_node",
    "n_kp_genes",
    "n_ecoli_genes",
    "n_pairs",
    "cross",
    "cross_p75",
    "cross_max",
    "rank_cross",
    "rank_cross_p75",
    "rank_cross_max",
]


def _map(tmp_path: Path, name: str, rep: str, rows) -> Path:
    path = tmp_path / f"{name}.tsv"
    path.write_text("\t".join(COLUMNS) + "\n" + "".join("\t".join(str(value) for value in row) + "\n" for row in rows))
    path.with_suffix(".json").write_text(
        json.dumps({"donor": "ecoli", "recipient": "kp", "rep": rep, "k": 256, "cap": 30, "recall": []})
    )
    return path


def test_a_column_labelled_with_the_wrong_representation_is_refused(tmp_path, monkeypatch):
    """⛔ `ESM=<a bacformer map>` inverts the comparison, and every number stays in range.

    The map records what it was built from, so the label is checked against it rather than trusted.
    """
    bacformer = _map(tmp_path, "b", "bacformer", [["1", "2", 7, 9, 63, 0.7, 0.7, 0.7, 1, 1, 1]])
    esm = _map(tmp_path, "e", "esm", [["1", "2", 7, 9, 63, 0.99, 0.99, 0.99, 1, 1, 1]])
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql+psycopg://unused/unused")
    monkeypatch.setattr("sys.argv", ["x", f"ESM={bacformer}", f"BACFORMER={esm}"])
    with pytest.raises(SystemExit, match="was built from BACFORMER but is labelled 'ESM'"):
        cmp.main()


def test_maps_sharing_no_pair_are_refused_rather_than_compared_on_nothing(tmp_path, monkeypatch):
    """⛔ An empty intersection would print an AUC of NaN per axis and read as 'no difference'."""
    a = _map(tmp_path, "a", "esm", [["1", "2", 7, 9, 63, 0.99, 0.99, 0.99, 1, 1, 1]])
    b = _map(tmp_path, "b", "bacformer", [["9", "8", 7, 9, 63, 0.70, 0.70, 0.70, 1, 1, 1]])
    monkeypatch.setenv("BACATLAS_DATABASE_URL", "postgresql+psycopg://unused/unused")
    monkeypatch.setattr("sys.argv", ["x", f"ESM={a}", f"BACFORMER={b}"])
    with pytest.raises(SystemExit, match="share no node pair"):
        cmp.main()


def test_one_map_is_not_a_comparison(tmp_path, monkeypatch):
    """The script exists to compare, so a single map is an error rather than a degenerate run."""
    a = _map(tmp_path, "a", "esm", [["1", "2", 7, 9, 63, 0.99, 0.99, 0.99, 1, 1, 1]])
    monkeypatch.setattr("sys.argv", ["x", f"ESM={a}"])
    with pytest.raises(SystemExit, match="at least two maps"):
        cmp.main()


def test_the_scores_come_through_the_SHARED_reader_so_labels_stay_strings(tmp_path):
    """⛔ Node labels look numeric; the shared `read_map` is what keeps them text on both sides."""
    path = _map(tmp_path, "a", "esm", [["1098", "02811", 7, 9, 63, 0.991, 0.993, 0.998, 1, 1, 1]])
    scores, provenance = cmp.read_scores(str(path), donor="ecoli", recipient="kp")
    assert list(scores) == [("1098", "02811")]
    assert scores[("1098", "02811")] == pytest.approx(0.991)
    assert provenance["rep"] == "esm"
