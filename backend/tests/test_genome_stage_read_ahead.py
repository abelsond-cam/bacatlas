"""The genome stage's reads moved onto a thread pool — and nothing the database receives may change.

⛔ **The standard is bit-exact, and it is proved against the loop as it was.** `_reference_load` below is
the genome-stage loop from bacatlas `73e59b0`, kept HERE as the oracle: serial, a cold `is_file()` per
genome, every file read inside `ingest_one_genome`. The new loop reads each genome's `_meta` columns
and `.nseq` gene table on a pool, takes them in ROSTER order, and must write the same `genome`,
`genome_assembly`, `genome_contig` and `genome_genes` rows, in the same physical order, with the same
ids, and refuse the same genomes with the same messages — at every pool size, including none.

The synthetic tree is built to hold the shapes a reordering or a mis-keyed read would get wrong: more
genomes than the read-ahead window, a `.nseq` naming another sample, a gene-count mismatch, a corrupt
`_meta`, a `_meta` missing a column, a named genome with no table, a genome with neither file, an
assembly stem with NCBI seqids, and a genome only the `.nseq` tree holds.

⚠ Two differences from `73e59b0` are by design, and both are confined to assembly stems — every probe
genome is a BioSample: `genome.bare_accession` is now written, and `genome_contig.contig_name` is the
seqid where it carries no `{sample}.` prefix. Both are asserted directly, not hidden in the comparison.
"""

from __future__ import annotations

import json
import os
import pathlib
import struct
import threading
from collections import Counter

import numpy as np
import pandas
import pytest
from sqlalchemy import func, select, text

from bacatlas_backend.ingest import artifact_locator
from bacatlas_backend.ingest import ingest_genome_gene_table as gene_table
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.genome_gene_alignment import GeneAlignmentError
from bacatlas_backend.ingest.ingest_command_line import (
    NSEQ_GATE_REFUSED_REASON,
    GenomeMetaMissing,
    discover_genomes,
    load_genome_layer,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import (
    GenomeIngestError,
    PrefetchedReads,
    _contig_name,
    bare_accession_for,
    ingest_one_genome,
)
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.models.genome import Genome

GENES, CONTIGS = 40, 6
STEM = "GCF_000512165.1_ASM51216v1_genomic"
NAMED = "SAMN09999990"
NEITHER = "SAMN09999991"
TREE_ONLY = "SAMN09999992"
WRONG_SAMPLE = "SAMN00000007"
COUNT_MISMATCH = "SAMN00000013"
CORRUPT_META = "SAMN00000021"
SHORT_META = "SAMN00000034"
TABLES = ("genome", "genome_assembly", "genome_contig", "genome_genes")


def _nseq_bytes(sample: str, cname: list[str], starts, ends, contigs) -> bytes:
    clen = [1000] * len(cname)
    header = json.dumps(
        {"sample": sample, "cname": cname, "clen": clen, "n_genes": len(starts),
         "n_bases": sum(clen), "ambiguous": []},
        separators=(",", ":"),
    ).encode()
    table = np.stack([starts, ends, contigs, np.zeros(len(starts), int)], axis=1).astype("<u4").tobytes()
    return b"NUNASEQ\x01" + struct.pack("<I", len(header)) + header + table + b"\0" * (-(-sum(clen) // 4))


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    """A `.nseq`-only data root: 61 BioSamples created in REVERSE roster order, plus the edge cases."""
    root = tmp_path_factory.mktemp("read_ahead")
    meta_dir = root / "proc" / "embeddings" / "meta"
    seq = root / "seq"
    meta_dir.mkdir(parents=True)
    seq.mkdir()
    rng = np.random.default_rng(7)
    biosamples = [f"SAMN{i:08d}" for i in range(61)]
    universe = [*reversed(biosamples), STEM, NAMED, NEITHER]
    for sample in [*universe, TREE_ONLY]:
        contigs = np.sort(rng.integers(0, CONTIGS, GENES))
        starts = np.arange(GENES) * 900 + 101 + int(rng.integers(0, 50))
        ends = starts + 899
        meta = pandas.DataFrame({"flat_index": np.arange(GENES), "contig_idx": contigs, "start": starts, "end": ends,
                                 "protein_sequence": ["M"] * GENES, "gene_name": ["g"] * GENES})
        if sample == SHORT_META:
            meta = meta.drop(columns=["end"])
        if sample == CORRUPT_META:
            (meta_dir / f"{sample}_meta.parquet").write_bytes(b"not a parquet")
        else:
            meta.sample(frac=1, random_state=1).to_parquet(meta_dir / f"{sample}_meta.parquet")  # unsorted rows
        if sample in (NAMED, NEITHER):
            continue
        cname = ([f"NZ_CP{j:06d}.1" for j in range(CONTIGS)] if sample == STEM
                 else [f"{sample}.contig{j:05d}" for j in range(CONTIGS)])
        if sample == COUNT_MISMATCH:  # one gene short of its `_meta`
            body = _nseq_bytes(sample, cname, starts[:-1], ends[:-1], contigs[:-1])
        elif sample == WRONG_SAMPLE:  # a file whose header names ANOTHER genome of the roster
            body = _nseq_bytes("SAMN00000008", cname, starts, ends, contigs)
        else:
            body = _nseq_bytes(sample, cname, starts, ends, contigs)
        (seq / f"{sample}.nseq").write_bytes(body)
    pandas.DataFrame({"sample_id": [s for s in universe for _ in range(2)]}).to_parquet(root / "universe.parquet")
    return root


def _artifacts(root) -> CatalogueArtifacts:
    return CatalogueArtifacts(
        data_root=root, set_key="kp71k", model_label="kp71k_hypergeo", run_id="r", species_key="kp",
        gene_universe_path=root / "universe.parquet", nseq_root=root / "seq",
    )


def _reset(session) -> None:
    session.execute(text(
        "TRUNCATE genome_genes, genome_contig, genome_assembly, genome, pathogen_species RESTART IDENTITY CASCADE"
    ))
    session.commit()


def _rows(session) -> dict:
    """Every row of the four tables, in PHYSICAL order — every column, values and types."""
    out = {}
    for table in TABLES:
        columns = session.execute(text(
            "select a.attname, format_type(a.atttypid, a.atttypmod) from pg_attribute a "
            f"where a.attrelid = 'public.{table}'::regclass and a.attnum > 0 and not a.attisdropped order by a.attnum"
        )).all()
        rows = session.execute(text(f"select * from {table} order by ctid")).all()
        out[table] = (columns, [tuple((type(v).__name__, v) for v in row) for row in rows])
    return out


def _reference_load(session, artifacts, *, species_key, named) -> list[tuple[str, str]]:
    """⛔ THE ORACLE — the genome-stage loop of bacatlas 73e59b0, reduced to what writes rows.

    Serial; `ingest_one_genome` called WITHOUT `prefetched`, so it stats the `_meta` and reads both
    files itself; the original except tuple. Returns the refusals as the report lists them.
    """
    species_ids = ingest_pathogen_species(session)
    session.commit()
    present = artifacts.nseq_present()
    refusals = []
    for sample_id, dataset_id in discover_genomes(artifacts):
        try:
            ingest_one_genome(
                session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
                pathogen_species_ids=species_ids, species_key=species_key,
                nseq_path=artifacts.per_genome_nseq(sample_id) if sample_id in present else None,
                gene_table_omitted_reason=NSEQ_GATE_REFUSED_REASON if sample_id in named else None,
            )
            session.commit()
        except (GenomeIngestError, GeneAlignmentError, KeyError, ValueError) as error:
            session.rollback()
            refusals.append((sample_id, str(error).replace("\n", " ")[:300]))
    return refusals


@pytest.fixture()
def reference(session, tree):
    _reset(session)
    refusals = _reference_load(session, _artifacts(tree), species_key="kp", named=frozenset({NAMED}))
    rows = _rows(session)
    _reset(session)
    yield refusals, rows
    _reset(session)


@pytest.mark.parametrize("workers", [0, 1, 2, 8])
def test_the_read_ahead_writes_exactly_the_rows_the_serial_loop_wrote(session, tree, reference, workers):
    """⭐ Same rows, same physical order, same ids, same refusals — at every pool size."""
    expected_refusals, expected_rows = reference
    report = load_genome_layer(
        session, _artifacts(tree), species_key="kp", genomes_without_gene_table=frozenset({NAMED}),
        progress_every=0, read_workers=workers,
    )
    session.commit()
    rows = _rows(session)

    # ⛔ Coverage first: the comparison must have something in it, column for column.
    assert set(rows) == set(expected_rows) == set(TABLES)
    for table in TABLES:
        assert rows[table][0] == expected_rows[table][0], f"{table}: columns or types differ"
        assert len(rows[table][1]) == len(expected_rows[table][1]) > 0, table
    assert len(rows["genome"][1]) == 60, "61 BioSamples + stem + named + neither + tree-only − 5 refused"
    assert report.refusals == expected_refusals
    assert {sample for sample, _ in report.refusals} == {
        NEITHER, WRONG_SAMPLE, COUNT_MISMATCH, CORRUPT_META, SHORT_META,
    }
    for table in TABLES:
        assert rows[table][1] == expected_rows[table][1], f"{table}: a row differs from the serial loop's"


def test_the_two_designed_differences_land_only_on_the_assembly_stem(session, tree):
    """`bare_accession` written for the stem and NULL for every BioSample; its NCBI seqids kept whole."""
    _reset(session)
    load_genome_layer(session, _artifacts(tree), species_key="kp", genomes_without_gene_table=frozenset({NAMED}),
                      progress_every=0)
    session.commit()
    accessions = dict(session.execute(select(Genome.sample_id, Genome.bare_accession)).all())
    assert accessions.pop(STEM) == "GCF_000512165.1"
    assert set(accessions.values()) == {None}
    names = session.execute(text(
        "select g.sample_id, c.seqid, c.contig_name from genome_contig c join genome g using (genome_id)"
    )).all()
    assert {(seqid, name) for sample, seqid, name in names if sample == STEM} == {
        (f"NZ_CP{j:06d}.1", f"NZ_CP{j:06d}.1") for j in range(CONTIGS)
    }
    assert all(name == seqid.split(".", 1)[1] for sample, seqid, name in names if sample != STEM)
    _reset(session)


def test_no_meta_is_stat_ed_per_genome_and_no_nseq_path_is_stat_ed_at_all(session, tree, monkeypatch):
    """⛔ The regression guard for the audit's ~3 cold path stats per genome: the `_meta` presence is
    one scan, and the `.nseq` length comes from `fstat` on the handle already open."""
    _reset(session)
    calls = {"is_file": 0, "nseq_stat": 0, "scans": 0}
    real_is_file, real_stat, real_scandir = pathlib.Path.is_file, os.stat, artifact_locator.os.scandir

    def is_file(self):
        calls["is_file"] += self.name.endswith("_meta.parquet")
        return real_is_file(self)

    def stat(path, *args, **kwargs):
        calls["nseq_stat"] += str(path).endswith(".nseq")
        return real_stat(path, *args, **kwargs)

    def scandir(path="."):
        calls["scans"] += str(path).endswith(os.path.join("embeddings", "meta"))
        return real_scandir(path)

    monkeypatch.setattr(pathlib.Path, "is_file", is_file)
    monkeypatch.setattr(os, "stat", stat)
    monkeypatch.setattr(artifact_locator.os, "scandir", scandir)
    report = load_genome_layer(session, _artifacts(tree), species_key="kp",
                               genomes_without_gene_table=frozenset({NAMED}), progress_every=0, read_workers=4)
    session.commit()
    assert report.genomes_loaded == 60
    assert calls == {"is_file": 0, "nseq_stat": 0, "scans": 1}
    _reset(session)


def witness_reads(monkeypatch) -> dict:
    """Record every `_meta` / `.nseq` read: what the consumer took, what it found prefetched, and on which thread.

    ⚠ **The rows cannot see a prefetch that misses.** `PrefetchedReads.take` falls through to the loader
    when the key it is asked for was not prefetched — by design, so a drift costs a read and never a wrong
    frame — and the genome loads exactly as before, only with its reads back on the serial SQL thread:
    the ~2.6–4.4 h the pool exists to remove, invisible to every row comparison. This is the witness for it.
    ``taken`` counts the consumer's requests by kind, ``missed`` lists each that fell through, ``loaded``
    counts the loaders' runs by kind, and ``loaded_on_sql_thread`` those run on the calling thread.
    """
    seen: dict = {"taken": Counter(), "missed": [], "loaded": Counter(), "loaded_on_sql_thread": []}
    sql_thread = threading.current_thread()
    real_take = PrefetchedReads.take

    def take(self, key, loader, *args):
        seen["taken"][key[0]] += 1
        if key not in self._results:
            seen["missed"].append(key)
        return real_take(self, key, loader, *args)

    monkeypatch.setattr(PrefetchedReads, "take", take)
    for kind, name in (("meta", "_load_meta"), ("nseq", "_load_nseq")):
        real = getattr(gene_table, name)

        def load(*args, _kind=kind, _real=real):
            seen["loaded"][_kind] += 1
            if threading.current_thread() is sql_thread:
                seen["loaded_on_sql_thread"].append((_kind, args))
            return _real(*args)

        monkeypatch.setattr(gene_table, name, load)
    return seen


@pytest.mark.parametrize("workers", [0, 1, 4])
def test_every_read_is_TAKEN_from_the_prefetch_and_none_runs_on_the_sql_thread(session, tree, monkeypatch, workers):
    """⛔ The speed-up's own gate: no request falls through to the loader, and with a pool none reads serially.

    Counts, not "some": 64 `_meta` requests (every roster genome but NEITHER, refused before a read; the
    corrupt and short `_meta` included — their stored errors are taken, not re-read) and 61 `.nseq`
    requests (those 64 less the two bad `_meta` and NAMED, which reads only `flat_index`). Each file is
    loaded exactly once — 63 `.nseq`, because the pool also reads the two whose genome the bad `_meta` then
    refuses, which is a read on a worker and nothing on the SQL thread — so a miss shows twice: in
    ``missed`` and as a load beyond these counts.
    """
    _reset(session)
    seen = witness_reads(monkeypatch)
    report = load_genome_layer(session, _artifacts(tree), species_key="kp",
                               genomes_without_gene_table=frozenset({NAMED}), progress_every=0, read_workers=workers)
    session.commit()
    assert report.genomes_loaded == 60
    assert seen["taken"] == {"meta": 64, "nseq": 61}
    assert seen["missed"] == []
    assert seen["loaded"] == {"meta": 64, "nseq": 63}
    if workers > 1:
        assert seen["loaded_on_sql_thread"] == []
    else:  # inline by design: the same reads, on the caller's thread
        assert len(seen["loaded_on_sql_thread"]) == 64 + 63
    _reset(session)


def test_a_roster_genome_with_no_meta_refuses_the_RUN_by_name_before_anything_is_written(session, tree, tmp_path):
    """⛔ In seconds, not hours in — and a NAMED genome is not exempt: its gene count is read from it."""
    _reset(session)
    meta_dir = tmp_path / "proc" / "embeddings" / "meta"
    meta_dir.mkdir(parents=True)
    samples = ["SAMN1", "SAMN2", NAMED]
    for sample in samples[:1]:
        pandas.DataFrame({"flat_index": [0]}).to_parquet(meta_dir / f"{sample}_meta.parquet")
    (meta_dir / "SAMN2_meta.parquet").mkdir()  # a DIRECTORY of that name is not a file — as `is_file()` said
    pandas.DataFrame({"sample_id": samples}).to_parquet(tmp_path / "universe.parquet")
    (tmp_path / "seq").mkdir()
    with pytest.raises(GenomeMetaMissing, match=r"2 of 3 roster genome\(s\) have no `_meta.parquet`.*"
                                                 rf"\['{NAMED}', 'SAMN2'\]"):
        load_genome_layer(session, _artifacts(tmp_path), species_key="kp",
                          genomes_without_gene_table=frozenset({NAMED}), progress_every=0)
    session.rollback()
    assert session.execute(select(func.count()).select_from(Genome)).scalar_one() == 0
    _reset(session)


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads a mode-000 file")
def test_an_unreadable_meta_is_ONE_named_refusal_and_the_run_continues(session, tree, tmp_path):
    """⚠ An `OSError` on a `_meta` (EIO, ESTALE, EACCES on Lustre) used to escape the handler and end the stage."""
    _reset(session)
    import shutil

    root = tmp_path / "copy"
    shutil.copytree(tree, root)
    locked = root / "proc" / "embeddings" / "meta" / "SAMN00000003_meta.parquet"
    locked.chmod(0)
    try:
        report = load_genome_layer(session, _artifacts(root), species_key="kp",
                                   genomes_without_gene_table=frozenset({NAMED}), progress_every=0)
    finally:
        locked.chmod(0o644)
    session.commit()
    refused = dict(report.refusals)
    assert "Permission denied" in refused["SAMN00000003"] and str(locked) in refused["SAMN00000003"]
    assert report.genomes_loaded == 59
    _reset(session)


# ── the two rules, unit by unit, against the definitions they replace or reuse ─────────────────
def _old_contig_name(seqid: str) -> str:
    """⛔ ORACLE — `_contig_rows`' rule at 73e59b0: everything after the first dot."""
    return seqid.split(".", 1)[1] if "." in seqid else seqid


@pytest.mark.parametrize(
    "sample,seqid",
    [
        ("SAMEA103923484", "SAMEA103923484.contig00324"),
        ("SAMN12345678", "SAMN12345678.contig00001"),
        ("SAMD00012345", "SAMD00012345.contig.with.dots"),
        ("SAMN12345678", "SAMN12345678."),
        ("SAMN12345678", "contig_without_a_dot"),
    ],
)
def test_the_contig_name_is_unchanged_for_every_bakrep_seqid(sample, seqid):
    assert _contig_name(seqid, sample) == _old_contig_name(seqid)


@pytest.mark.parametrize(
    "sample,seqid,name",
    [
        (STEM, "NZ_CP003999.1", "NZ_CP003999.1"),
        (STEM, "CP004000.1", "CP004000.1"),
        (STEM, f"{STEM}.contig1", "contig1"),
        ("SAMN12345678", "SAMN87654321.contig00001", "SAMN87654321.contig00001"),
    ],
)
def test_the_contig_name_is_never_a_fragment_of_an_accession(sample, seqid, name):
    """The old rule gave `1`, `1`, `1_ASM51216v1_genomic.contig1` and `contig00001` for these."""
    assert _contig_name(seqid, sample) == name


@pytest.mark.parametrize(
    "sample,expected",
    [
        (STEM, "GCF_000512165.1"),
        ("GCA_003010635.1_INF156_genomic", "GCA_003010635.1"),
        ("GCF_000016305.1_KSB2_2B_genomic", "GCF_000016305.1"),
        ("gcf_000016305_ASM1_genomic", "gcf_000016305"),
        ("SAMN12345678", None),
        ("SAMEA103923484", None),
        ("SAMD00012345", None),
        ("GCF_000512165.1", None),  # not a stem: no tail, nothing to strip — NULL, and the gate names it
        ("ERR123456", None),
    ],
)
def test_the_bare_accession_is_nunas_rule_and_NULL_where_that_rule_does_not_apply(sample, expected):
    from nuna.tl.stability.build_genome_set import bare_accession

    assert bare_accession_for(sample) == expected
    if expected is not None:
        assert bare_accession(sample) == expected, "one rule, nuna's"
