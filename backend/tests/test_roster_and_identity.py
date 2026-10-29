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
        ("--allow-missing-strand", "allow_missing_strand"),
    ],
)
def test_every_omission_is_a_flag_that_must_be_TYPED(flag, attribute):
    """⭐ Each of these leaves a table out. Default off, so an absence is always a decision."""
    assert getattr(build_parser().parse_args([]), attribute) is False
    assert getattr(build_parser().parse_args([flag]), attribute) is True


def test_the_noncoding_table_is_OPT_IN_because_nothing_reads_it():
    """⛔ The polarity is inverted from the other two, and the inversion is the point.

    `--no-functional-annotations` omits a table the page *would* read if it were cheap, so leaving it
    out is a decision that has to be typed. `genome_noncoding_feature` has no reader at all — traced
    2026-10-29 across both repos — so WRITING it is the decision, and at 71,435 genomes it is the
    largest table in the database. A test rather than a comment, because the previous default was
    justified in the help text by a renderer that does not exist.
    """
    assert build_parser().parse_args([]).with_noncoding_features is False
    assert build_parser().parse_args(["--with-noncoding-features"]).with_noncoding_features is True


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


# ── T1 · the coordinate convention, which `need_length=False` would otherwise have retired ────
def _convention_store(tmp_path: Path, samples: list[str], *, inclusive: bool = True) -> Path:
    """A `_meta` store whose coordinates do, or do not, follow the GFF's 1-based-inclusive rule.

    ⚠ The wrong case is the REAL one, not a scribble: 0-based half-open is what every other
    genomics convention uses, so `end - start` rather than `end - start + 1` is the mistake an
    importer actually makes. It is one base out on every gene and nothing raises.
    """
    directory = tmp_path / "embeddings" / "meta"
    directory.mkdir(parents=True, exist_ok=True)
    aa = [100, 250, 60, 410]
    for sample_id in samples:
        start = [1, 2_000, 5_000, 9_000]
        span = [3 * (n + 1) for n in aa]
        pairs = list(zip(start, span, strict=True))
        end = [s + w - 1 for s, w in pairs] if inclusive else [s + w for s, w in pairs]
        pandas.DataFrame({
            "sample_id": [sample_id] * len(aa),
            "flat_index": list(range(len(aa))),
            "contig_idx": [0] * len(aa),
            "start": start,
            "end": end,
            "protein_sequence": ["M" * n for n in aa],
        }).to_parquet(directory / f"{sample_id}_meta.parquet")  # fmt: skip
    return tmp_path


def test_a_1_based_inclusive_cohort_passes_the_convention_check_silently(tmp_path):
    from bacatlas_backend.ingest.catalogue_frames import _gate_length_convention

    omitted: dict = {}
    _gate_length_convention(_convention_store(tmp_path, ["SAMN1", "SAMN2"]), ["SAMN1", "SAMN2"], omitted)
    assert omitted == {}


def test_a_half_open_cohort_is_RECORDED_because_every_block_width_would_be_wrong(tmp_path, capsys):
    """⛔ The check `need_length=False` would have silently retired. One base out on every gene:
    no exception, no empty table, just a track drawn slightly wrong forever."""
    from bacatlas_backend.ingest.catalogue_frames import _gate_length_convention

    omitted: dict = {}
    processed = _convention_store(tmp_path, ["SAMN1", "SAMN2"], inclusive=False)
    _gate_length_convention(processed, ["SAMN1", "SAMN2"], omitted)

    assert "length_convention" in omitted, "a collapsed convention must not pass unrecorded"
    assert "0.0%" in omitted["length_convention"], omitted["length_convention"]
    assert "1-based" in omitted["length_convention"]
    assert "⚠" in capsys.readouterr().out


def test_a_sample_that_yields_no_protein_length_says_UNCHECKED_not_nothing(tmp_path):
    """⚠ `not measured` and `measured fine` are the same silence unless one of them speaks."""
    from bacatlas_backend.ingest.catalogue_frames import _gate_length_convention

    (tmp_path / "embeddings" / "meta").mkdir(parents=True)
    omitted: dict = {}
    _gate_length_convention(tmp_path, ["SAMN1", "SAMN2"], omitted)
    assert "UNCHECKED" in omitted["length_convention"]


def test_the_convention_check_READS_A_BOUNDED_SEEDED_SAMPLE_not_the_cohort(tmp_path, monkeypatch):
    """⛔ The gate exists inside the change that stopped reading `protein_sequence` for the cohort.
    If it reads every genome it has given the whole saving back, and nothing else would notice —
    the numbers it prints would be right either way. So the genomes it asks for are the assertion.
    """
    import nuna.eval.reference_annotation as reference_annotation

    from bacatlas_backend.ingest.catalogue_frames import (
        LENGTH_CONVENTION_SAMPLE,
        _gate_length_convention,
    )

    samples = [f"SAMN{i:04d}" for i in range(LENGTH_CONVENTION_SAMPLE * 4)]
    processed = _convention_store(tmp_path, samples)

    asked: list[list[str]] = []
    real = reference_annotation.load_gene_coords

    def recording(proc, drawn, **kwargs):
        asked.append(list(drawn))
        return real(proc, drawn, **kwargs)

    monkeypatch.setattr(reference_annotation, "load_gene_coords", recording)
    _gate_length_convention(processed, samples, {})
    _gate_length_convention(processed, samples, {})

    assert len(asked) == 2
    assert len(asked[0]) == LENGTH_CONVENTION_SAMPLE, f"read {len(asked[0])} of {len(samples)}"
    assert asked[0] == asked[1], "seeded — two runs of one cohort must check the same genomes"
    assert set(asked[0]) < set(samples)


def test_a_cohort_smaller_than_the_sample_is_read_WHOLE_and_still_checked(tmp_path):
    """⚠ `random.sample` raises when k exceeds the population — the probe cohorts are 100 and 280."""
    from bacatlas_backend.ingest.catalogue_frames import _gate_length_convention

    omitted: dict = {}
    _gate_length_convention(_convention_store(tmp_path, ["SAMN1"]), ["SAMN1"], omitted)
    assert omitted == {}


def test_the_cohort_read_ASKS_FOR_NO_PROTEIN_SEQUENCE_which_is_the_whole_input_saving(
    tmp_path, monkeypatch
):
    """⛔ T1, gated where it lives: `_load_source_frames` must call `load_gene_coords` LEAN.

    `protein_sequence` is **92.59 %** of the `_meta` store's column bytes (measured over all 280
    probe genomes, 458.7 MiB of column chunks, 1,436,421 genes) and `length_aa` — the only thing
    derived from it — is read **nowhere** in this repository. At 71,435 genomes the keyword is the
    difference between ~124.5 GiB and ~6.9 GiB of input.

    ⚠ **It is gated here because losing it is INVISIBLE.** The lean call returns the identical frame
    minus one unread column, so every other test, every row count and the whole database are
    byte-for-byte unchanged either way. Nothing but the clock would notice a revert, and the clock
    is only read at a scale where the run is already hours long.

    ⭐ This is also the only test that executes `_load_source_frames` at all — the orchestration
    above seven tested loaders, which until now was reached only through a full database ingest.
    """
    import nuna.eval.assignments as assignments
    import nuna.eval.reference_annotation as reference_annotation

    from bacatlas_backend.ingest import catalogue_frames

    samples = ["SAMN1", "SAMN2"]
    universe = pandas.DataFrame({
        "sample_id": samples * 2, "flat_index": [0, 0, 1, 1], "node": ["n0", "n0", "n1", "n1"],
    })  # fmt: skip
    coords = pandas.DataFrame({
        "sample_id": samples * 2, "flat_index": [0, 0, 1, 1], "contig_idx": [0] * 4,
        "start": [1, 1, 100, 100], "end": [60, 60, 160, 160],
        "strand": ["+", "-", "+", "-"], "rank_on_contig": [0, 0, 1, 1],
    })  # fmt: skip

    seen: dict = {}

    def fake_load_gene_coords(proc, requested, **kwargs):
        seen["kwargs"] = kwargs
        return coords

    monkeypatch.setattr(assignments, "gene_universe", lambda *a, **k: universe)
    monkeypatch.setattr(assignments, "load_assignments", lambda *a, **k: universe)
    monkeypatch.setattr(reference_annotation, "load_gene_coords", fake_load_gene_coords)
    for empty in ("load_bakta_products", "load_bakta_dbxrefs", "load_bakta_noncoding"):
        monkeypatch.setattr(reference_annotation, empty, lambda *a, **k: pandas.DataFrame())
    monkeypatch.setattr(catalogue_frames, "_gate_strand", lambda *a, **k: None)
    monkeypatch.setattr(catalogue_frames, "_gate_length_convention", lambda *a, **k: None)

    source = catalogue_frames._load_source_frames(_artifacts(tmp_path))

    assert seen["kwargs"].get("need_length") is False, (
        f"the cohort read must be lean; it asked {seen['kwargs']!r}"
    )
    assert "length_aa" not in source["coords"].columns
    assert source["samples"] == samples and source["n_genes"] == 4


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


# ── T3 · the categorical window knob, and the merge that would NULL a column rather than raise ──
# ⭐ nuna measured the knob on its 100-genome probe: 132.8 → 16.4 B/row, 614.4 → 75.7 MB, **8.12×**,
# with `assert_frame_equal` exact. That is why `build_catalogue_frames` can reach it. ⚠ It is default
# OFF — it buys nothing at 280 genomes and costs ~1.5× the build — but the 71k ingest will want it, and
# the thing standing between here and there is this file's own `_gene_memberships`.
#
# ⛔ **Its hazard does not raise.** `_gene_memberships` factorizes a MultiIndex over the window frame and
# then LEFT-joins the result back onto the assignment on `(sample_id, flat_index)`. Under
# `categorical=True` that join has a **categorical `flat_index` on one side and a plain `int64` on the
# other**. pandas 2.3.3 matches them by value — but had it not, every row would still be returned, each
# with a NULL `arrangement_rank`, which is a *documented legitimate value* (a gene alone on its contig
# has no window). The table would be the right length, the column would be all nulls, and the only thing
# that could tell the difference is a test that asserts the ranks are actually there.


def _window_inputs():
    """A cohort big enough that most genes have a full ±5 window and several share an arrangement."""
    import pandas as pd

    rows = [
        {"sample_id": f"SAMEA{genome:07d}", "flat_index": index, "contig_idx": 0,
         "start": 1_000 * index, "strand": "+" if index % 3 else "-", "node": f"n{index % 7}"}
        for genome in (1, 2, 3) for index in range(14)
    ]  # fmt: skip
    assign_coords = pd.DataFrame(rows)
    node_position = {n: i for i, n in enumerate(sorted(assign_coords["node"].unique()))}
    sample_position = {s: i for i, s in enumerate(sorted(assign_coords["sample_id"].unique()))}
    return assign_coords, node_position, sample_position


def _memberships(categorical: bool):
    from nuna.tl.locus_browser.export_payload import (
        OFFSETS,
        arrangements,
        oriented_windows,
    )

    from bacatlas_backend.ingest import catalogue_frames

    assign_coords, node_position, sample_position = _window_inputs()
    windows = oriented_windows(assign_coords, node_position, radius=5, categorical=categorical)
    slot_of = {offset: index for index, offset in enumerate(OFFSETS)}
    arrangement_frame = arrangements(windows, slot_of, top=None, sample_pos=sample_position)
    assign = assign_coords[["sample_id", "flat_index", "node"]]
    got = catalogue_frames._gene_memberships(assign, windows, arrangement_frame, slot_of)
    return got, windows


def test_categorical_windows_give_the_IDENTICAL_gene_memberships():
    """⛔ Same rows, same order, same ranks — the only claim a memory knob may make of this table."""
    import pandas as pd

    plain, plain_w = _memberships(categorical=False)
    cat, cat_w = _memberships(categorical=True)
    assert str(cat_w["sample_id"].dtype) == "category", "the fixture never reached the categorical path"
    assert str(plain_w["sample_id"].dtype) != "category"

    back = cat.copy()
    for column in ("node", "sample_id", "flat_index"):
        back[column] = back[column].astype(plain[column].dtype)
    pd.testing.assert_frame_equal(plain, back)


def test_the_arrangement_RANKS_ARE_ACTUALLY_THERE_not_merely_a_column_of_the_right_length():
    """⛔ The assertion the identity test cannot make for itself.

    If the categorical/int64 join silently matched nothing, BOTH paths' ranks would be null together
    only if the plain path were broken too — but a reader comparing two frames would see them agree and
    call it proved. So the rank column is asserted non-null on its own terms, in both paths: a NULL rank
    is a *legitimate* value here (a gene with no window), which is exactly what makes an all-NULL column
    indistinguishable from a correct one without this line.
    """
    for categorical in (False, True):
        got, _ = _memberships(categorical=categorical)
        placed = got["arrangement_rank"].notna().sum()
        assert placed > 0, f"categorical={categorical}: every arrangement_rank is NULL — the join missed"
        assert placed == len(got), (
            f"categorical={categorical}: {len(got) - placed} of {len(got)} genes lost their rank; "
            "every gene in this fixture has a full window, so none should be unplaced"
        )


def test_build_catalogue_frames_ACTUALLY_FORWARDS_the_knob_to_the_window_build():
    """⛔ A switch nothing passes on is a switch that silently does nothing — the T1 lesson again.

    ⚠ Asserted on the SOURCE, and the reason is worth stating rather than hiding: `build_catalogue_frames`
    is called by **no test in this repository** — only by `ingest_locus_catalogue`, i.e. only inside a full
    database ingest — so there is no harness in which to observe the call. Until there is, the source is
    the only thing that can tell `categorical=categorical_windows` from a parameter that was added to the
    signature and then dropped on the floor, because with the default OFF the two behave identically.
    """
    import inspect

    from bacatlas_backend.ingest.catalogue_frames import build_catalogue_frames

    source = inspect.getsource(build_catalogue_frames)
    assert "categorical_windows: bool = False" in source or "categorical_windows" in str(
        inspect.signature(build_catalogue_frames)
    ), "the parameter is gone"
    assert "categorical=categorical_windows" in source, (
        "`oriented_windows` is called without the knob, so `categorical_windows=True` would do nothing"
    )


def test_gene_memberships_is_EXACTLY_ONE_ROW_PER_ASSIGNED_GENE_in_both_dtypes():
    """⛔ The function's own documented contract — *"exactly one row per gene, always"* — asserted.

    ⚠ **Added because a mutation survived both tests above, and the reason it survived is the lesson.**
    Replace the `MultiIndex.factorize` with an enumeration of the index LEVELS — the categorical trap, and
    a change someone might make believing it equivalent — and the intermediate gains a row per *possible*
    (node, genome, gene) triple. Those extra rows then give `windowed` several rows per
    `(sample_id, flat_index)`, and the LEFT join multiplies the output. Neither test above sees it: the
    identity test compares the two dtypes against **each other**, and they break together; the rank test
    scales with the row count. Only the contract catches it, and only stated against the ASSIGNMENT.
    """
    assign_coords, _, _ = _window_inputs()
    n_genes = len(assign_coords)
    for categorical in (False, True):
        got, _ = _memberships(categorical=categorical)
        assert len(got) == n_genes, (
            f"categorical={categorical}: {len(got)} rows for {n_genes} assigned genes"
        )
        assert not got.duplicated(subset=["sample_id", "flat_index"]).any(), (
            f"categorical={categorical}: a gene appears more than once"
        )
