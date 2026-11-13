"""A genome holds ONE gene table — and a catalogue that reads it by position pins it.

`genome_genes` is keyed by the genome alone; a catalogue's `genome_locus_assignment` block says
"gene *i* of that table is in locus L". Two stores that annotate one assembly differently therefore
cannot both be loaded: the second genome stage would replace the table under the first catalogue's
positions. The guard refuses exactly that — *different bytes under a catalogue that reads them* — and
must refuse nothing else: not an identical reload, not a change to a genome no catalogue reads.
"""

from __future__ import annotations

import json
import struct

import numpy as np
import pandas
import pytest
from sqlalchemy import select

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_genome_gene_table import GenomeIngestError, ingest_one_genome
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment

SAMPLE = "SAMN03892119"
GENES = 12
REASON = "the .nseq acceptance gate refused this genome"


def _nseq_bytes(sample: str, starts, ends) -> bytes:
    """A minimal `.nseq`: one contig, `len(starts)` genes."""
    header = json.dumps(
        {"sample": sample, "cname": [f"{sample}.contig00001"], "clen": [50_000], "n_genes": len(starts),
         "n_bases": 50_000, "ambiguous": []},
        separators=(",", ":"),
    ).encode()  # fmt: skip
    zeros = np.zeros(len(starts), int)
    table = np.stack([starts, ends, zeros, zeros], axis=1).astype("<u4").tobytes()
    return b"NUNASEQ\x01" + struct.pack("<I", len(header)) + header + table + b"\0" * (-(-50_000 // 4))


def _store(root, *, shift: int = 0, genes: int = GENES):
    """Write one genome's `_meta` and `.nseq` — `shift` moves every gene, `genes` changes how many."""
    meta_dir, seq = root / "proc" / "embeddings" / "meta", root / "seq"
    meta_dir.mkdir(parents=True, exist_ok=True)
    seq.mkdir(exist_ok=True)
    starts = np.arange(genes) * 900 + 101 + shift
    ends = starts + 899
    pandas.DataFrame({
        "flat_index": np.arange(genes), "contig_idx": np.zeros(genes, int), "start": starts, "end": ends,
        "protein_sequence": ["M"] * genes, "gene_name": ["g"] * genes,
    }).to_parquet(meta_dir / f"{SAMPLE}_meta.parquet")  # fmt: skip
    (seq / f"{SAMPLE}.nseq").write_bytes(_nseq_bytes(SAMPLE, starts, ends))
    return CatalogueArtifacts(
        data_root=root, set_key="kp71k", model_label="kp71k_hypergeo", run_id="r", species_key="kp",
        nseq_root=seq,
    )  # fmt: skip


def _load(session, artifacts, *, named: bool = False):
    """One genome through the `.nseq`-only route, inside a SAVEPOINT — the unit the real loop rolls back."""
    species_ids = ingest_pathogen_species(session)
    with session.begin_nested():
        return ingest_one_genome(
            session, artifacts, sample_id=SAMPLE, bakrep_dataset_id=None, pathogen_species_ids=species_ids,
            species_key="kp", nseq_path=None if named else artifacts.per_genome_nseq(SAMPLE),
            gene_table_omitted_reason=REASON if named else None,
        )  # fmt: skip


def _block(session):
    session.expire_all()
    return session.execute(
        select(GenomeGenes.gene_count, GenomeGenes.gene_spans)
        .join(Genome, Genome.genome_id == GenomeGenes.genome_id)
        .where(Genome.sample_id == SAMPLE)
    ).one_or_none()


def _catalogue_reads(session, seeded, genome_id: int, gene_count: int) -> None:
    """A loaded catalogue that addresses this genome's table by position."""
    session.add(GenomeLocusAssignment(
        pangenome_id=seeded["pangenome"].pangenome_id, genome_id=genome_id, gene_count=gene_count,
        placed_gene_count=0, distinct_locus_count=0,
        gene_locus_ordinals=b"\xff" * (4 * gene_count), gene_arrangement_ranks=b"\xff" * (4 * gene_count),
    ))  # fmt: skip
    session.flush()


def test_a_different_table_under_a_catalogue_that_reads_it_is_refused_by_name(session, seeded, tmp_path):
    first = _load(session, _store(tmp_path))
    before = _block(session)
    _catalogue_reads(session, seeded, first.genome_id, GENES)

    with pytest.raises(GenomeIngestError) as refusal:
        _load(session, _store(tmp_path, shift=3))  # the same genes, three bases along: another store's table
    message = str(refusal.value)
    assert SAMPLE in message and "probe-nuna4" in message and "by position" in message
    assert tuple(_block(session)) == tuple(before), "the refused load must leave the stored table untouched"


def test_a_table_with_another_gene_count_is_refused_the_same_way(session, seeded, tmp_path):
    first = _load(session, _store(tmp_path))
    _catalogue_reads(session, seeded, first.genome_id, GENES)
    with pytest.raises(GenomeIngestError, match=f"{GENES + 3:,} genes"):
        _load(session, _store(tmp_path, genes=GENES + 3))
    assert _block(session).gene_count == GENES


def test_the_SAME_table_reloads_under_a_catalogue_because_a_resume_must(session, seeded, tmp_path):
    artifacts = _store(tmp_path)
    first = _load(session, artifacts)
    before = _block(session)
    _catalogue_reads(session, seeded, first.genome_id, GENES)
    again = _load(session, artifacts)
    assert again.genome_id == first.genome_id and again.gene_spans_written == 1
    assert tuple(_block(session)) == tuple(before)


def test_a_genome_NO_catalogue_reads_is_free_to_change(session, tmp_path):
    _load(session, _store(tmp_path))
    before = _block(session)
    _load(session, _store(tmp_path, shift=3))
    after = _block(session)
    assert after.gene_count == before.gene_count and bytes(after.gene_spans) != bytes(before.gene_spans)


def test_naming_the_genome_as_holding_no_table_cannot_delete_one_a_catalogue_reads(session, seeded, tmp_path):
    artifacts = _store(tmp_path)
    first = _load(session, artifacts)
    before = _block(session)
    _catalogue_reads(session, seeded, first.genome_id, GENES)
    with pytest.raises(GenomeIngestError, match="named as holding none"):
        _load(session, artifacts, named=True)
    assert tuple(_block(session)) == tuple(before)
    contigs = session.execute(select(GenomeContig).where(GenomeContig.genome_id == first.genome_id)).all()
    assert len(contigs) == 1, "the contig rows go with the block, so they must survive the refusal too"


def test_naming_a_genome_no_catalogue_reads_still_removes_its_table(session, tmp_path):
    artifacts = _store(tmp_path)
    _load(session, artifacts)
    report = _load(session, artifacts, named=True)
    assert report.gene_table_omitted_reason == REASON
    assert _block(session) is None
