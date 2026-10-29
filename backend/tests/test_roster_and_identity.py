"""The ingest pointed at a cohort that does not follow the probe's conventions.

⛔ **Every case here is a SILENT SUCCESS that was reachable before these changes.** None of them
raised, none of them logged, and each produced a database that looked like a smaller correct load:

| what | what it did | what it looked like |
|---|---|---|
| no nested GFF tree | roster globbed to `[]` | `0 loaded, 0 refused` — a cohort with no genomes |
| `_meta` with no `species` column | a bare `KeyError` out of pyarrow | the whole run aborted |
| one unreadable parquet | `ArrowInvalid` is a `ValueError`, uncaught | the whole run aborted |
| `{set}_close_seq.parquet` not the universe's name | `FileNotFoundError` | — (the one that was loud) |
| `_strand.parquet` short | that genome's whole gene set forced to `+` | every arrow forward, A5 → 0.5 |

⚠ The three that are loud are here too, because the fix for each was an *override*, and an
override that resolves to the wrong file is the next silent success.
"""

from __future__ import annotations

from pathlib import Path

import pandas
import pytest

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_command_line import (
    EmptyRoster,
    annotation_index,
    build_parser,
    discover_genomes,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import GenomeIngestError, ingest_one_genome
from bacatlas_backend.ingest.published_catalogues import PUBLISHED_CATALOGUES
from tests.conftest import NUNA_DATA_ROOT, artifacts_for


def _universe(path: Path, samples: list[str], genes_each: int = 3) -> Path:
    """A close_seq parquet: one row per gene, which is all the roster needs of it."""
    pandas.DataFrame(
        {
            "sample_id": [s for s in samples for _ in range(genes_each)],
            "flat_index": [i for _ in samples for i in range(genes_each)],
        }
    ).to_parquet(path)
    return path


def _meta(tmp_path: Path, sample_id: str, *, species: str | None = None) -> Path:
    """A `{BS}_meta.parquet`. ⭐ `species=None` is the kp71k schema: ten columns, no `species`."""
    directory = tmp_path / "proc" / "embeddings" / "meta"
    directory.mkdir(parents=True, exist_ok=True)
    frame = pandas.DataFrame(
        {
            "flat_index": [0, 1], "contig_idx": [0, 0], "start": [1, 100], "end": [60, 160],
            "protein_sequence": ["MA", "MB"], "gene_name": ["a", "b"],
        }
    )
    if species is not None:
        frame["species"] = species
    path = directory / f"{sample_id}_meta.parquet"
    frame.to_parquet(path)
    return path


def _artifacts(tmp_path: Path, **overrides) -> CatalogueArtifacts:
    fields = dict(
        data_root=tmp_path, set_key="kp", model_label="kp_probe", run_id="kp_probe_run"
    )
    fields.update(overrides)
    return CatalogueArtifacts(**fields)


# ── R1 · the roster is the gene universe, and an empty one is refused ─────────────────────────
def test_the_roster_comes_from_the_universe_and_NOT_from_the_gff_tree(tmp_path):
    """⛔ The defect: with no `<datasetID>/<BS>/<BS>.bakta.gff3.gz` tree the glob matched nothing,
    so the genome layer loaded zero genomes and reported success."""
    artifacts = _artifacts(
        tmp_path, gene_universe_path=_universe(tmp_path / "u.parquet", ["SAMN1", "SAMN2"])
    )
    assert not artifacts.annotation_root.exists(), "the point of the case is that there is no tree"

    found = discover_genomes(artifacts)

    assert [sample_id for sample_id, _ in found] == ["SAMN1", "SAMN2"]  # the universe, in its order
    assert [dataset_id for _, dataset_id in found] == [None, None], (
        "a genome the tree does not hold must be present-and-unlocated, not absent"
    )


def test_an_empty_roster_is_REFUSED_rather_than_loaded(tmp_path):
    """⭐ Loading nothing succeeds. That is the whole bug, and it needs its own refusal."""
    artifacts = _artifacts(
        tmp_path, gene_universe_path=_universe(tmp_path / "u.parquet", [], genes_each=0)
    )
    with pytest.raises(EmptyRoster, match="names a genome"):
        discover_genomes(artifacts)


def test_a_genome_in_the_tree_and_in_no_universe_is_STILL_in_the_roster(tmp_path):
    """⚠ The projection depends on exactly these. `--stage projection` places genomes the model
    never clustered and `projected_gene_placement` is keyed onto `gene (genome_id, flat_index)`,
    so a roster narrowed to the catalogue's own universe would silently stop them loading.
    Measured on the probe store: 280 GFFs, 200 genomes across both universes, and 10 of the other
    80 are the placed genomes."""
    directory = tmp_path / "raw" / "gff" / "d1" / "SAMN_projectable"
    directory.mkdir(parents=True)
    (directory / "SAMN_projectable.bakta.gff3.gz").write_bytes(b"")
    artifacts = _artifacts(
        tmp_path, gene_universe_path=_universe(tmp_path / "u.parquet", ["SAMN1"])
    )

    found = dict(discover_genomes(artifacts))

    assert found == {"SAMN1": None, "SAMN_projectable": "d1"}


def test_the_union_is_empty_only_when_BOTH_sources_are(tmp_path):
    """⛔ The refusal is on the union, because either source alone may legitimately be empty: a
    cohort served from streamed assemblies has no tree, and a store mirrored before its catalogue
    has no universe."""
    directory = tmp_path / "raw" / "gff" / "d1" / "SAMN_only"
    directory.mkdir(parents=True)
    (directory / "SAMN_only.bakta.gff3.gz").write_bytes(b"")
    empty = _universe(tmp_path / "u.parquet", [], genes_each=0)
    assert discover_genomes(_artifacts(tmp_path, gene_universe_path=empty)) == [
        ("SAMN_only", "d1")
    ]


def test_the_roster_order_is_the_universes_own_sorted_order(tmp_path):
    """⚠ `collection_genome_ordinal` is a position in exactly this list, and `arr.gid` indexes into
    it — so the genome layer and the roster must derive it by ONE rule, which is why both call
    `read_genome_vocabulary` rather than sorting independently."""
    from bacatlas_backend.ingest.ingest_genome_collection_roster import read_genome_vocabulary

    artifacts = _artifacts(
        tmp_path,
        gene_universe_path=_universe(tmp_path / "u.parquet", ["SAMN9", "SAMN1", "SAMN5"]),
    )
    roster, n_genes = read_genome_vocabulary(artifacts)
    assert [sample_id for sample_id, _ in discover_genomes(artifacts)] == roster
    assert (roster, n_genes) == (["SAMN1", "SAMN5", "SAMN9"], 9)


def test_a_genome_under_two_dataset_directories_is_refused(tmp_path):
    """One genome must resolve to one GFF, or the sequence endpoint serves whichever won the glob."""
    for dataset_id in ("d1", "d2"):
        directory = tmp_path / "raw" / "gff" / dataset_id / "SAMN1"
        directory.mkdir(parents=True)
        (directory / "SAMN1.bakta.gff3.gz").write_bytes(b"")
    with pytest.raises(RuntimeError, match="more than one dataset directory"):
        annotation_index(_artifacts(tmp_path))


def test_a_roster_genome_with_no_gff_is_refused_BY_NAME(tmp_path):
    """⚠ Missing and absent must not look the same: it is refused, not quietly left out."""
    _meta(tmp_path, "SAMN404")
    artifacts = _artifacts(tmp_path)
    with pytest.raises(GenomeIngestError, match="SAMN404.*under no dataset directory"):
        ingest_one_genome(
            None, artifacts, sample_id="SAMN404", bakrep_dataset_id=None,
            pathogen_species_ids={"kp": 1}, species_key="kp",
        )


# ── R3 · the species is named, because 300 of 300 kp71k `_meta` files have no such column ─────
def test_a_named_species_key_means_the_species_column_is_never_read(tmp_path):
    """⛔ The parquet is read with an explicit column list, so requesting an absent `species`
    raises inside pyarrow before a single gene is seen. Naming the species must not merely
    override the value — it must stop the column being asked for."""
    _meta(tmp_path, "SAMN1")
    artifacts = _artifacts(tmp_path)

    # Without the key, the absent column is what fails — and it fails in pyarrow, not in a check.
    with pytest.raises(Exception, match="species"):
        ingest_one_genome(
            None, artifacts, sample_id="SAMN1", bakrep_dataset_id="d1",
            pathogen_species_ids={"kp": 1},
        )
    # With it, the read succeeds and the run reaches the GFF — the next real check.
    with pytest.raises(GenomeIngestError, match="no GFF at"):
        ingest_one_genome(
            None, artifacts, sample_id="SAMN1", bakrep_dataset_id="d1",
            pathogen_species_ids={"kp": 1}, species_key="kp",
        )


def test_a_species_key_with_no_row_is_refused_naming_what_is_known(tmp_path):
    _meta(tmp_path, "SAMN1")
    artifacts = _artifacts(tmp_path)
    with pytest.raises(GenomeIngestError, match="has no row"):
        ingest_one_genome(
            None, artifacts, sample_id="SAMN1", bakrep_dataset_id="d1",
            pathogen_species_ids={"kp": 1}, species_key="not_a_species",
        )


# ── R4/R5 · the two overrides, and that neither changes the derived path ──────────────────────
def test_the_derived_paths_are_unchanged_when_no_override_is_given(tmp_path):
    """⛔ An override must be an override. A second naming rule would move every probe catalogue."""
    artifacts = _artifacts(tmp_path)
    assert artifacts.gene_universe == tmp_path / "proc" / "mmseqs" / "kp_close_seq.parquet"
    assert artifacts.assignment == (
        tmp_path / "proc" / "analysis" / "step4_global_merge" / "assignments" / "kp_probe_run.tsv"
    )


def test_each_override_names_its_file_outright(tmp_path):
    """⭐ The 71k universe is `kleb_close_seq_kp71k_id0.98.parquet` and its assignment is under
    `comparison/nuna/<model>/<run>/step4/` — neither spelling is reachable by any rule."""
    universe = tmp_path / "kleb_close_seq_kp71k_id0.98.parquet"
    assignment = tmp_path / "comparison" / "nuna" / "hypergeo" / "kp71k_s0" / "step4" / "a.parquet"
    artifacts = _artifacts(tmp_path, gene_universe_path=universe, assignment_path=assignment)
    assert artifacts.gene_universe == universe
    assert artifacts.assignment == assignment


def test_a_tilde_in_an_override_is_expanded_like_every_other_path(tmp_path):
    artifacts = _artifacts(tmp_path, gene_universe_path=Path("~/u.parquet"))
    assert "~" not in str(artifacts.gene_universe)


def test_the_parquet_sidecar_resolves_ONLY_when_the_tsv_is_absent(tmp_path):
    """⛔ Preferring it where both exist would change which file every published catalogue was
    built from, silently, as a side effect of a sidecar appearing beside it."""
    directory = tmp_path / "proc" / "analysis" / "step4_global_merge" / "assignments"
    directory.mkdir(parents=True)
    artifacts = _artifacts(tmp_path)

    (directory / "kp_probe_run.parquet").write_bytes(b"")
    assert artifacts.assignment.suffix == ".parquet", "a parquet-only run must be addressable"

    (directory / "kp_probe_run.tsv").write_text("")
    assert artifacts.assignment.suffix == ".tsv", "the TSV is the primary artifact when it exists"


# ── the CLI surface the four overrides are reached through ────────────────────────────────────
@pytest.mark.parametrize(
    "flag,attribute",
    [
        ("--gene-universe", "gene_universe"),
        ("--assignment", "assignment"),
        ("--species-key", "species_key"),
    ],
)
def test_every_override_is_reachable_from_the_command_line(flag, attribute):
    """⚠ `--species-key` already existed and fed `CatalogueArtifacts` only — never the genome
    loader — which is exactly why the kp71k schema was not survivable before this."""
    args = build_parser().parse_args([flag, "value"])
    assert getattr(args, attribute.replace("-", "_")) is not None


@pytest.mark.parametrize(
    "flag,attribute",
    [
        ("--no-functional-annotations", "no_functional_annotations"),
        ("--no-noncoding-features", "no_noncoding_features"),
        ("--allow-missing-strand", "allow_missing_strand"),
    ],
)
def test_every_omission_is_a_flag_that_must_be_TYPED(flag, attribute):
    """⭐ Each of these leaves a table out. Default off, so an absence is always a decision."""
    assert getattr(build_parser().parse_args([]), attribute) is False
    assert getattr(build_parser().parse_args([flag]), attribute) is True


# ── R2 · one unreadable parquet is a refusal, not an aborted run ──────────────────────────────
def test_arrow_invalid_IS_a_value_error(tmp_path):
    """⛔ The fact the handler rests on, pinned because it is pyarrow's and not ours.

    `ArrowInvalid` subclasses `ValueError`, which is why it slipped past a handler catching
    `(GenomeIngestError, GeneAlignmentError, KeyError)` — and why adding `ValueError` fixes it.
    """
    import pyarrow

    corrupt = tmp_path / "corrupt.parquet"
    corrupt.write_bytes(b"this is not a parquet file")
    with pytest.raises(ValueError):
        pandas.read_parquet(corrupt)
    assert issubclass(pyarrow.lib.ArrowInvalid, ValueError)


@pytest.mark.usefixtures("session")
def test_one_unreadable_parquet_refuses_THAT_genome_and_the_run_continues(tmp_path, session):
    """⛔ Before this, an `ArrowInvalid` escaped the per-genome handler and killed the whole load —
    the opposite of what a per-genome unit of work is for. The run must end with a report."""
    from bacatlas_backend.ingest.ingest_command_line import load_genome_layer

    _universe(tmp_path / "u.parquet", ["SAMN_corrupt", "SAMN_missing"])
    directory = tmp_path / "proc" / "embeddings" / "meta"
    directory.mkdir(parents=True)
    (directory / "SAMN_corrupt_meta.parquet").write_bytes(b"not a parquet")

    report = load_genome_layer(
        session,
        _artifacts(tmp_path, gene_universe_path=tmp_path / "u.parquet"),
        species_key="kp",
        progress_every=0,
    )

    assert (report.genomes_loaded, report.genomes_refused) == (0, 2)
    refused = dict(report.refusals)
    assert set(refused) == {"SAMN_corrupt", "SAMN_missing"}, "every refusal is named, not counted"
    assert "parquet" in refused["SAMN_corrupt"].lower()
    assert report.genome_ids == []


# ── S1/S3 · strand, where partial coverage is worse than none ─────────────────────────────────
def _strand_store(tmp_path: Path, samples: list[str]) -> Path:
    directory = tmp_path / "embeddings" / "meta"
    directory.mkdir(parents=True, exist_ok=True)
    for sample_id in samples:
        pandas.DataFrame({"flat_index": [0], "strand": ["+"]}).to_parquet(
            directory / f"{sample_id}_strand.parquet"
        )
    return tmp_path


def test_a_short_strand_store_is_refused_and_the_missing_genomes_are_NAMED(tmp_path):
    """⛔ `load_gene_coords` fills an absent sidecar with `+` for that genome's whole gene set and
    only PRINTS a counter — every arrow forward, every window unreversed, A5 pulled to 0.5."""
    from bacatlas_backend.ingest.catalogue_frames import MissingStrand, _gate_strand

    processed = _strand_store(tmp_path, ["SAMN1"])
    with pytest.raises(MissingStrand) as refusal:
        _gate_strand(processed, ["SAMN1", "SAMN2"], {}, allow_missing=False)
    assert "SAMN2" in str(refusal.value), "a count is not a miss list"
    assert "1 of 2" in str(refusal.value).replace(",", "")


def test_a_complete_strand_store_passes_silently(tmp_path):
    from bacatlas_backend.ingest.catalogue_frames import _gate_strand

    omitted: dict = {}
    _gate_strand(_strand_store(tmp_path, ["SAMN1", "SAMN2"]), ["SAMN1", "SAMN2"], omitted,
                 allow_missing=False)  # fmt: skip
    assert omitted == {}


def test_proceeding_past_a_short_store_RECORDS_the_shortfall(tmp_path):
    """⭐ `omitted` is where an absence is stated, so it is never read as a measured zero."""
    from bacatlas_backend.ingest.catalogue_frames import _gate_strand

    omitted: dict = {}
    _gate_strand(_strand_store(tmp_path, ["SAMN1"]), ["SAMN1", "SAMN2"], omitted, allow_missing=True)
    assert "SAMN2" in omitted["strand"] and "PROCEEDING" in omitted["strand"]


def test_sidecars_that_all_read_plus_are_refused_TOO(tmp_path):
    """⚠ The second silent-failure path: this passes `strand_coverage` by construction and then
    produces exactly the uniform map an absent sidecar does."""
    from bacatlas_backend.ingest.catalogue_frames import MissingStrand, _gate_strand_balance

    degenerate = pandas.DataFrame({"strand": ["+"] * 100})
    with pytest.raises(MissingStrand, match="degenerate"):
        _gate_strand_balance(degenerate, {}, allow_missing=False)

    balanced = pandas.DataFrame({"strand": ["+"] * 50 + ["-"] * 50})
    omitted: dict = {}
    _gate_strand_balance(balanced, omitted, allow_missing=False)
    assert omitted == {}


# ── the gate: nothing about the two published catalogues moved ────────────────────────────────
@pytest.mark.parametrize("entry", list(PUBLISHED_CATALOGUES), ids=lambda e: e.species_key)
def test_the_published_rosters_are_PAIR_FOR_PAIR_what_the_glob_produced(entry):
    """⛔ The regression guard for this whole change, against the real store.

    Both published catalogues were built from a roster the old glob derived. The union must
    reproduce it exactly — same genomes, same dataset ids — or every `collection_genome_ordinal`
    moves, and `arr.gid` indexes into that list, so every arrangement would name a different
    genome while each name stayed a real one.

    ⚠ The reference is the OLD implementation written out here, not a recorded count: a count
    agreeing proves only that the two lists are the same length.
    """
    if not (NUNA_DATA_ROOT / "raw" / "gff").is_dir():
        pytest.skip(f"the cluster artifact mirror is not present at {NUNA_DATA_ROOT}")
    artifacts = artifacts_for(entry.species_key)

    before = [
        (path.parent.name, path.parent.parent.name)
        for path in sorted(artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz"))
    ]
    after = discover_genomes(artifacts)

    assert len(before) == 280, "the probe store is 280 GFFs — the premise of this test moved"
    assert sorted(after) == sorted(before)
    assert not [sample_id for sample_id, dataset_id in after if dataset_id is None], (
        "every probe genome is located; an unlocated one would be refused and the load short"
    )
