"""`python -m bacatlas_backend.ingest` — the offline loader.

⛔ **Writes are OURS ALONE, and never on a request path.** This is the only thing in the package
that inserts, and it runs on a machine that has the cluster artifacts. The serving install has
neither them nor `nuna`.

⭐ **Every stage reports what it did, in counts.** A loader that prints "done" has told you nothing
you can check; one that prints `280 genomes · 1,436,421 genes · 0 refused` can be reconciled against
the store it read. A refusal is named with the genome and the check that caught it, and the run
continues — so one bad genome costs one line, not the whole load.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

import bacatlas_backend.models  # noqa: F401  (registers every table on the metadata)
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.genome_gene_alignment import GeneAlignmentError
from bacatlas_backend.ingest.ingest_genome_collection_roster import (
    genome_id_by_ordinal,
    ingest_genome_collection,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import GenomeIngestError, ingest_one_genome
from bacatlas_backend.ingest.ingest_locus_catalogue import ingest_locus_catalogue
from bacatlas_backend.ingest.ingest_nuna_model_registry import (
    ingest_model_registry,
    model_key_for_audit_label,
)
from bacatlas_backend.ingest.ingest_pangenome_run import catalogue_key_for, ingest_pangenome_run
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.ingest.ingest_projected_genomes import ProjectionRefused, load_projection
from bacatlas_backend.ingest.ingest_reference_vocabularies import load_reference_vocabularies
from bacatlas_backend.ingest.publish_pangenome import PublishRefused, publish_pangenome
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation, GenomeNoncodingFeature
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.nuna_model import NunaModel
from bacatlas_backend.models.pangenome import Pangenome


@dataclass
class GenomeLayerReport:
    """The whole genome-layer load, in numbers that can be reconciled against the store."""

    genomes_loaded: int = 0
    genomes_refused: int = 0
    genes: int = 0
    contigs: int = 0
    functional_annotations: int = 0
    noncoding_features: int = 0
    seconds: float = 0.0
    refusals: list[tuple[str, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def render(self) -> str:
        """The report as the lines the CLI prints — counts first, then every refusal by name."""
        lines = [
            f"genome layer: {self.genomes_loaded:,} loaded, {self.genomes_refused:,} refused, "
            f"{self.seconds:.0f}s",
            f"  genes                  {self.genes:,}",
            f"  contigs                {self.contigs:,}",
            f"  functional annotations {self.functional_annotations:,}",
            f"  non-coding features    {self.noncoding_features:,}",
        ]
        for sample_id, reason in self.refusals:
            lines.append(f"  REFUSED {sample_id}: {reason}")
        for note in self.notes:
            lines.append(f"  NOTE {note}")
        return "\n".join(lines)


def discover_genomes(artifacts: CatalogueArtifacts) -> list[tuple[str, str]]:
    """`(sample_id, bakrep_dataset_id)` for every GFF in the store.

    ⚠ The dataset id is a PATH component and appears in no parquet — captured here or lost. The GFF
    sits one level deeper than `<BS>/`, which is why the glob has three segments and not two.
    """
    pairs = [
        (path.parent.name, path.parent.parent.name)
        for path in sorted(artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz"))
    ]
    duplicates = {name for name, _ in pairs if sum(1 for other, _ in pairs if other == name) > 1}
    if duplicates:
        raise RuntimeError(
            f"{len(duplicates)} sample id(s) appear under more than one dataset directory: "
            f"{sorted(duplicates)[:5]}. One genome must resolve to one GFF, or the sequence "
            "endpoint would serve whichever the glob happened to reach first."
        )
    return pairs


def load_genome_layer(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    limit: int | None = None,
    only: set[str] | None = None,
    progress_every: int = 25,
) -> GenomeLayerReport:
    """Load every genome in the store, one atomic transaction each."""
    report = GenomeLayerReport()
    species_ids = ingest_pathogen_species(session)
    session.commit()

    genomes = discover_genomes(artifacts)
    if only:
        genomes = [pair for pair in genomes if pair[0] in only]
    if limit:
        genomes = genomes[:limit]

    started = time.time()
    for index, (sample_id, dataset_id) in enumerate(genomes, start=1):
        try:
            one = ingest_one_genome(
                session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
                pathogen_species_ids=species_ids,
            )
            session.commit()
        except (GenomeIngestError, GeneAlignmentError, KeyError) as error:
            # ⛔ Roll back THIS genome and carry on. A genome is the unit of work, so one refusal
            # must not discard the 200 already loaded — and it must be named, not counted.
            session.rollback()
            report.genomes_refused += 1
            report.refusals.append((sample_id, str(error).replace("\n", " ")[:300]))
            continue
        report.genomes_loaded += 1
        report.genes += one.genes_written
        report.contigs += one.contigs_written
        report.functional_annotations += one.functional_annotations_written
        report.noncoding_features += one.noncoding_features_written
        report.notes.extend(one.notes)
        if progress_every and index % progress_every == 0:
            print(
                f"  {index}/{len(genomes)}  {report.genes:,} genes  {time.time() - started:.0f}s",
                file=sys.stderr, flush=True,
            )
    report.seconds = time.time() - started
    return report


def reconcile(session: Session, report: GenomeLayerReport) -> list[str]:
    """⛔ What the loader WROTE against what the table HOLDS — the two are not the same claim.

    A loader reporting its own writes confirms itself. This asks the database, and a difference means
    something else wrote (or deleted) rows underneath the run.
    """
    held = {
        "genes": session.execute(select(func.count()).select_from(Gene)).scalar_one(),
        "contigs": session.execute(select(func.count()).select_from(GenomeContig)).scalar_one(),
        "functional_annotations": session.execute(
            select(func.count()).select_from(GeneFunctionalAnnotation)
        ).scalar_one(),
        "noncoding_features": session.execute(
            select(func.count()).select_from(GenomeNoncodingFeature)
        ).scalar_one(),
        "genomes": session.execute(select(func.count()).select_from(Genome)).scalar_one(),
    }
    written = {
        "genes": report.genes, "contigs": report.contigs,
        "functional_annotations": report.functional_annotations,
        "noncoding_features": report.noncoding_features, "genomes": report.genomes_loaded,
    }
    return [
        f"{key}: wrote {written[key]:,} but the table holds {held[key]:,}"
        for key in written
        if written[key] != held[key]
    ]


def load_pangenome_layer(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    species_key: str,
    catalogue_key: str | None = None,
    collection_key: str | None = None,
) -> list[str]:
    """Model registry → roster → run → catalogue, in the one order the foreign keys allow.

    ⚠ **The genome layer must already be loaded**, and this refuses rather than loading a partial
    roster: every arrangement's membership is a position in the genome vocabulary, so a hole in it
    would address the wrong genome and nothing downstream could detect that.
    """
    from nuna.tl.cluster.nuna_pipeline import MODELS

    species_ids = ingest_pathogen_species(session)
    session.flush()
    species_id = species_ids[species_key]

    model_key = model_key_for_audit_label(artifacts.model_label, artifacts.set_key, MODELS)
    from bacatlas_backend.ingest.ingest_pangenome_run import nuna_git_sha

    ingest_model_registry(session, models=MODELS, registry_git_sha=nuna_git_sha())
    session.flush()
    model = session.execute(
        select(NunaModel).where(NunaModel.model_key == model_key)
    ).scalar_one()

    collection_id, roster = ingest_genome_collection(
        session, artifacts, pathogen_species_id=species_id, collection_key=collection_key
    )
    session.flush()

    pangenome_id, run_report = ingest_pangenome_run(
        session,
        artifacts,
        pathogen_species_id=species_id,
        genome_collection_id=collection_id,
        nuna_model_id=model.nuna_model_id,
        registry_steps=list(model.steps),
        genome_count=roster.genome_count,
        gene_count=roster.genes_in_universe,
        # ⚠ Zero until the loci exist; `ingest_locus_catalogue` sets the real count in the same
        # transaction, so a committed row never claims a locus count it does not hold.
        locus_count=0,
        # ⭐ How a reader addresses this catalogue. Defaults to `{species}-{model}`, which is also
        # what nuna's exporter writes as `--dset`, so one key spans both surfaces.
        catalogue_key=catalogue_key_for(species_key, model_key, catalogue_key),
    )
    catalogue_report = ingest_locus_catalogue(
        session,
        artifacts,
        pangenome_id=pangenome_id,
        pathogen_species_id=species_id,
        genome_id_by_ordinal=genome_id_by_ordinal(session, collection_id),
    )
    return [
        f"catalogue: {catalogue_key_for(species_key, model_key, catalogue_key)}",
        f"model: {model_key} ({model.step_count} steps, {model.exclusivity_form.value})",
        f"collection: {roster.collection_key} — {roster.genome_count} genomes, "
        f"{roster.genes_in_universe:,} genes in the universe",
        run_report.render(),
        catalogue_report.render(),
    ]


def reconcile_pangenome(session: Session, artifacts: CatalogueArtifacts) -> list[str]:
    """⛔ Ask the DATABASE what it holds, against what the checked-in triple says it should.

    The published counts are a fact of record (`published_catalogues.py`), so this compares the load
    against them rather than against itself.
    """
    from bacatlas_backend.ingest.published_catalogues import PUBLISHED_CATALOGUES

    entry = next((e for e in PUBLISHED_CATALOGUES if e.run_id == artifacts.run_id), None)
    if entry is None:
        return []
    # ⚠ `run_id` alone is NOT unique — `uq_pangenome__run_id_ingest_generation` is the key, and a
    # re-ingest bumps the generation rather than mutating rows. Unqualified, this raised
    # `MultipleResultsFound` the first time anyone re-ingested. Newest generation is the live one.
    pangenome = session.execute(
        select(Pangenome)
        .where(Pangenome.run_id == artifacts.run_id)
        .order_by(Pangenome.ingest_generation.desc())
        .limit(1)
    ).scalar_one()
    held_loci = session.execute(
        select(func.count()).select_from(Locus).where(Locus.pangenome_id == pangenome.pangenome_id)
    ).scalar_one()
    differences = []
    for name, held, expected in (
        ("genome_count", pangenome.genome_count, entry.genome_count),
        ("gene_count", pangenome.gene_count, entry.gene_count),
        ("locus_count", held_loci, entry.locus_count),
    ):
        if held != expected:
            differences.append(
                f"{name}: the database holds {held:,} but the published catalogue is {expected:,}"
            )
    return differences


def build_parser() -> argparse.ArgumentParser:
    """The CLI surface. `--limit` and `--only` are smoke-run conveniences, not modes to load in."""
    parser = argparse.ArgumentParser(prog="python -m bacatlas_backend.ingest", description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("~/developer/nuna/data"),
                        help="the local mirror of the cluster artifacts")
    parser.add_argument("--set-key", default="ecoli", help="the dataset token in filenames")
    # ⚠ NOT `required=True`: `--stage reference` loads the public vocabularies, which belong to no
    # model and no run, and demanding a run id there would mean naming an arbitrary catalogue to
    # load a table that has nothing to do with it. `main` requires them for every other stage.
    parser.add_argument("--model-label", default=None, help="the audit's --model-label")
    parser.add_argument("--run-id", default=None, help="the assignment stem")
    parser.add_argument("--database-url", default=os.environ.get("BACATLAS_DATABASE_URL"))
    parser.add_argument("--stage", choices=("reference", "genomes", "pangenome", "projection", "all"),
                        default="all",
                        help="`reference` is the public vocabularies, global to every species and "
                             "model; `genomes` is model-INDEPENDENT and a new pangenome must never "
                             "rewrite it; `pangenome` needs it already loaded; `projection` places "
                             "genomes the model never clustered and is DELIBERATELY not in `all` — "
                             "it is an addition to a published catalogue, not part of building one")
    parser.add_argument("--projection-dir", type=Path, default=None,
                        help="where nuna wrote *_placed_windows.tsv (default: the locator's own path)")
    parser.add_argument("--publish", action="store_true",
                        help="after loading, VERIFY the catalogue and point its species at it. The "
                             "pointer is what the service reads, so this is a separate, explicit "
                             "step: a load that is not published changes nothing a reader can see")
    parser.add_argument("--species-key", default=None,
                        help="the browser key (`ecoli` | `kp`); defaults to the artifacts' own")
    parser.add_argument("--offer", action="store_true",
                        help="after loading, SHOW this catalogue in the picker without making it the "
                             "default. ⭐ The three decisions are separate: load it, offer it, serve "
                             "it by default. This is the middle one — it is what lets nuna4 stay the "
                             "default and the rollback while nuna5 is selectable beside it")  # fmt: skip
    parser.add_argument("--collection-key", default=None,
                        help="the genome roster this catalogue is built on; defaults to the "
                             "artifacts' set key. ⛔ Two catalogues of one species SHARE a collection, "
                             "and the membership ordinal is what arr.gid indexes into — so a model "
                             "over a genuinely different roster needs its own key here, or the "
                             "roster guard will refuse the load")  # fmt: skip
    parser.add_argument("--catalogue-key", default=None,
                        help="how a reader addresses THIS catalogue (`ecoli-nuna5`). Defaults to "
                             "`{species}-{model}`, which is also what nuna's exporter writes as "
                             "`--dset`. ⛔ Hyphen, never underscore: the keys must stay prefix-free "
                             "or one catalogue's payloads match another's key")  # fmt: skip
    parser.add_argument("--limit", type=int, help="stop after N genomes (a smoke run, not a mode)")
    parser.add_argument("--only", nargs="*", help="load only these sample ids")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run a load and return a shell exit code: 0 clean, 1 refusals or unreconciled, 2 misconfigured."""
    args = build_parser().parse_args(argv)
    if not args.database_url:
        print("BACATLAS_DATABASE_URL is not set and --database-url was not given", file=sys.stderr)
        return 2

    # ⛔ The reference stage reads no artifacts at all — the vocabularies are global — so it must
    # not be made to name a catalogue, and must not verify one.
    if args.stage == "reference":
        engine = create_engine(args.database_url, future=True)
        with Session(engine) as session:
            print(load_reference_vocabularies(session).render())
            session.commit()
        return 0

    if not args.model_label or not args.run_id:
        print("--model-label and --run-id are required for this stage", file=sys.stderr)
        return 2

    artifacts = CatalogueArtifacts(
        data_root=args.data_root, set_key=args.set_key,
        model_label=args.model_label, run_id=args.run_id,
    )
    optional = artifacts.verify()
    print(f"artifacts: all {len(artifacts.required())} required present; optional {optional}")

    engine = create_engine(args.database_url, future=True)
    report, differences = None, []
    with Session(engine) as session:
        # ⚠ First in `all`: the vocabularies are global and idempotent, so re-running costs one
        # upsert and an absent reference degrades chips rather than failing — but the count is
        # printed either way, so a silently missing table is visible here rather than inferred from
        # blank chips three screens into the page.
        if args.stage == "all":
            print(load_reference_vocabularies(session).render())
            session.commit()
        if args.stage in ("genomes", "all"):
            report = load_genome_layer(
                session, artifacts, limit=args.limit, only=set(args.only) if args.only else None
            )
            session.commit()
            print(report.render())
            differences += reconcile(session, report)
        # ⛔ NOT part of `all`. A projection is an addition to a catalogue that already exists, and
        # folding it into the build would make every load depend on an optional feature's artifacts.
        if args.stage == "projection":
            # ⚠ Newest generation, for the reason recorded in `reconcile_pangenome`: `run_id` is
            # unique only WITH `ingest_generation`, so this threw once a second generation existed.
            pangenome_id = session.execute(
                select(Pangenome.pangenome_id)
                .where(Pangenome.run_id == artifacts.run_id)
                .order_by(Pangenome.ingest_generation.desc())
                .limit(1)
            ).scalar_one_or_none()
            if pangenome_id is None:
                print(f"no pangenome with run_id {artifacts.run_id} — load it first", file=sys.stderr)
                return 2
            try:
                print(load_projection(
                    session,
                    pangenome_id=int(pangenome_id),
                    projection_root=args.projection_dir or artifacts.projection_root,
                ).render())
                session.commit()
            except ProjectionRefused as error:
                session.rollback()
                print(f"PROJECTION REFUSED: {error}", file=sys.stderr)
                return 1
            return 0
        if args.stage in ("pangenome", "all"):
            for line in load_pangenome_layer(
                session,
                artifacts,
                species_key=args.species_key or artifacts.species_key,
                catalogue_key=args.catalogue_key,
                collection_key=args.collection_key,
            ):
                print(line)
            session.commit()
            differences += reconcile_pangenome(session, artifacts)
            if args.offer:
                from bacatlas_backend.ingest.publish_pangenome import offer_catalogue

                # ⚠ Read the key back off the row that was just written rather than recomputing it.
                # Recomputing would be a second implementation of the same rule, and the two could
                # disagree about exactly the string a reader addresses this catalogue by.
                just_loaded = session.execute(
                    select(Pangenome.catalogue_key)
                    .where(Pangenome.run_id == artifacts.run_id)
                    .order_by(Pangenome.ingest_generation.desc())
                    .limit(1)
                ).scalar_one()
                print(offer_catalogue(session, catalogue_key=just_loaded))
                session.commit()
            if args.publish:
                try:
                    print(publish_pangenome(session, run_id=artifacts.run_id).render())
                    session.commit()
                except PublishRefused as error:
                    session.rollback()
                    print(f"PUBLISH REFUSED: {error}", file=sys.stderr)
                    differences.append("the catalogue was loaded but not published")

    # ⚠ Reconciliation is only meaningful on a FULL load: a --limit or --only run writes a subset
    # into a table that may already hold more, so a difference there is expected, not a fault.
    if args.limit or args.only:
        print("reconciliation skipped: a partial load cannot be reconciled against the whole table")
    elif differences:
        print("RECONCILIATION FAILED:", file=sys.stderr)
        for difference in differences:
            print(f"  {difference}", file=sys.stderr)
        return 1
    else:
        print("reconciled: every count the loader reported is what the tables hold")
    return 1 if (report and report.refusals) else 0


if __name__ == "__main__":
    raise SystemExit(main())
