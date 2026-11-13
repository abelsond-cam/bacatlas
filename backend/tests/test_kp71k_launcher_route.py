"""The 71,435-genome launcher's route through `main()` — driven with its EXACT flag set.

⛔ **The probe could not see this class of bug, by construction.** Both probe cohorts have
`set_key == species_key` (`ecoli`/`ecoli`, `kp`/`kp`), so a key that silently defaulted from one to the
other resolved to the right cohort anyway. The kp71k launcher (`nuna setup/csd3/ingest_bacatlas_kp71k.sbatch`)
passes `--set-key kp71k --species-key kp` and NO `--collection-key` / `--collection-label`; `main()` built
the artifacts without the species key, so the roster asked `published_cohorts` for `('kp71k', 'kp71k')`,
which nobody declares, and `--stage pangenome` died with `UndeclaredCohort` hours after the genome stage
had committed — and the launcher's EXIT trap then deleted the database.

The second half is the NCBI station's one input: an assembly stem's `bare_accession`, now gated.
"""

from __future__ import annotations

import pandas
import pytest
from sqlalchemy import delete

from bacatlas_backend.ingest import ingest_command_line
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.ingest.publish_pangenome import verify_pangenome_is_servable
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from tests.conftest import PROBE_URL

SET, SPECIES, LABEL = "kp71k", "kp", "kp71k_hypergeo"
RUN_ID = (
    "kp71k_s0_bacformer_clever_exploded_preclusterbfglobal0.5rhoPAIRMAX-p_hypergeo-rhofree0.05"
    "_k100_g71435_res0.6_seed0"
)
SAMPLES = ["GCF_000512165.1_ASM51216v1_genomic", "SAMEA103923484", "SAMN12345678"]


class ReachedTheRun(Exception):
    """Raised in place of `ingest_pangenome_run`: the roster resolved, and this is how far the test drives."""


def _launcher_tree(tmp_path):
    """A mirror `main()` accepts with the launcher's overrides: universe, assignment, .nseq dir, the list."""
    universe = tmp_path / "kleb_close_seq_kp71k_id0.98.parquet"
    pandas.DataFrame({"sample_id": [s for s in SAMPLES for _ in range(3)]}).to_parquet(universe)
    assignment = tmp_path / "comparison" / "nuna" / "hypergeo" / "kp71k_s0" / "step4" / "assignments"
    assignment = assignment / f"{RUN_ID}.parquet"
    nseq_dir = tmp_path / "seq_kp71k"
    nseq_dir.mkdir()
    no_gene_table = tmp_path / "genomes_without_gene_table.txt"
    no_gene_table.write_text(f"{SAMPLES[2]}\n")
    mirror = tmp_path / "bacatlas_mirror" / "kp71k"
    located = CatalogueArtifacts(
        data_root=mirror, set_key=SET, model_label=LABEL, run_id=RUN_ID, species_key=SPECIES,
        gene_universe_path=universe, assignment_path=assignment, nseq_root=nseq_dir,
    )
    for _role, path, _what in located.required():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch(exist_ok=True)
    return [  # ⛔ the launcher's COMMON, flag for flag, in its order
        "--data-root", str(mirror), "--set-key", SET, "--species-key", SPECIES, "--model-label", LABEL,
        "--run-id", RUN_ID, "--gene-universe", str(universe), "--assignment", str(assignment),
        "--nseq-dir", str(nseq_dir), "--genomes-without-gene-table", str(no_gene_table),
    ]


def test_the_pangenome_stage_with_the_launchers_flags_resolves_cohort_71k_and_key_kp_71k(
    session, tmp_path, monkeypatch
):
    """⛔ Fails on 73e59b0 with `UndeclaredCohort: no cohort declared for species 'kp71k' + dataset 'kp71k'`."""
    common = _launcher_tree(tmp_path)
    species_id = ingest_pathogen_species(session)[SPECIES]
    session.add_all([
        Genome(sample_id=sample, pathogen_species_id=species_id, strand_is_observed=True,
               sample_id_kind=SampleIdentifierKind.ASSEMBLY_STEM if sample.startswith("GCF_")
               else SampleIdentifierKind.BIOSAMPLE)
        for sample in SAMPLES
    ])
    session.commit()

    seen: dict = {}

    def run_recorder(session, artifacts, **kwargs):
        seen.update(kwargs, artifacts_species_key=artifacts.species_key)
        raise ReachedTheRun

    monkeypatch.setattr(ingest_command_line, "ingest_pangenome_run", run_recorder)
    monkeypatch.setenv("BACATLAS_DATABASE_URL", PROBE_URL)
    try:
        with pytest.raises(ReachedTheRun):
            ingest_command_line.main(["--stage", "pangenome", "--offer", *common])
    finally:
        session.rollback()
        session.execute(delete(Genome).where(Genome.sample_id.in_(SAMPLES)))
        session.commit()

    assert seen["artifacts_species_key"] == SPECIES
    assert seen["catalogue_key"] == "kp-71k"
    assert seen["genome_count"] == len(SAMPLES)
    assert seen["gene_count"] == 3 * len(SAMPLES)


# ── the NCBI station's one input, gated ───────────────────────────────────────────────────────
STEM_CHECK = "every assembly-stem genome has a bare accession for the NCBI station"


def _stem_check(results: list[str]) -> list[str]:
    return [line for line in results if line.startswith(STEM_CHECK)]


def test_publish_names_an_assembly_stem_without_a_bare_accession_and_passes_once_it_has_one(session, seeded):
    pangenome, species = seeded["pangenome"], seeded["species"]
    passed, failed = verify_pangenome_is_servable(session, pangenome)
    assert _stem_check(passed) == [f"{STEM_CHECK} — 0 assembly stems in the collection"], "BioSamples only: a pass"
    assert not _stem_check(failed)

    stems = [
        Genome(sample_id=sample, pathogen_species_id=species.pathogen_species_id, strand_is_observed=True,
               sample_id_kind=SampleIdentifierKind.ASSEMBLY_STEM, bare_accession=accession)
        for sample, accession in (
            ("GCF_000512165.1_ASM51216v1_genomic", None),
            ("GCA_003010635.1_INF156_genomic", "GCA_003010635.1"),
        )
    ]
    session.add_all(stems)
    session.flush()
    session.add_all([
        GenomeCollectionMembership(
            genome_collection_id=seeded["collection"].genome_collection_id, genome_id=genome.genome_id,
            collection_genome_ordinal=900 + index, requested_sample_id=genome.sample_id, matched_by="sample_id",
        )
        for index, genome in enumerate(stems)
    ])
    session.flush()

    passed, failed = verify_pangenome_is_servable(session, pangenome)
    assert not _stem_check(passed)
    (line,) = _stem_check(failed)
    assert "2 assembly stems in the collection" in line
    assert "1 of 2 assembly-stem genomes have no bare_accession" in line
    assert "['GCF_000512165.1_ASM51216v1_genomic']" in line

    stems[0].bare_accession = "GCF_000512165.1"
    session.flush()
    passed, failed = verify_pangenome_is_servable(session, pangenome)
    assert _stem_check(passed) == [f"{STEM_CHECK} — 2 assembly stems in the collection"]
    assert not _stem_check(failed)
