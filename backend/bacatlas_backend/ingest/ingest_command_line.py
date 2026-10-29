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

from sqlalchemy import ARRAY, Integer, any_, create_engine, func, literal, select
from sqlalchemy.orm import Session

import bacatlas_backend.models  # noqa: F401  (registers every table on the metadata)
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.genome_gene_alignment import GeneAlignmentError
from bacatlas_backend.ingest.ingest_genome_collection_roster import (
    genome_id_by_ordinal,
    ingest_genome_collection,
    read_genome_vocabulary,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import GenomeIngestError, ingest_one_genome
from bacatlas_backend.ingest.ingest_locus_catalogue import ingest_locus_catalogue
from bacatlas_backend.ingest.ingest_pangenome_run import catalogue_key_for, ingest_pangenome_run
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.ingest.ingest_projected_genomes import ProjectionRefused, load_projection
from bacatlas_backend.ingest.ingest_reference_vocabularies import load_reference_vocabularies
from bacatlas_backend.ingest.publish_pangenome import PublishRefused, publish_pangenome
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation, GenomeNoncodingFeature
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.locus import Locus
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
    #: ⭐ The genomes THIS run wrote — `reconcile`'s denominator. Without it the only available
    #: comparison is against the whole table, which equals the run's own writes only when the
    #: table started empty.
    genome_ids: list[int] = field(default_factory=list)

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


class EmptyRoster(RuntimeError):
    """A roster with no genomes in it — refused, because loading nothing succeeds."""


def annotation_index(artifacts: CatalogueArtifacts) -> dict[str, str]:
    """`sample_id → bakrep_dataset_id`, from ONE walk of the annotation tree.

    ⚠ The dataset id is a PATH component and appears in no parquet — captured here or lost. The GFF
    sits one level deeper than `<BS>/`, which is why the glob has three segments and not two.

    ⛔ **Locating a genome, which is not the same question as which genomes exist.** An absent tree
    returns an empty index and every genome is then refused BY NAME, where it used to produce an
    empty roster and a successful load of nothing.
    """
    index: dict[str, str] = {}
    duplicates: set[str] = set()
    for path in sorted(artifacts.annotation_root.glob("*/*/*.bakta.gff3.gz")):
        sample_id, dataset_id = path.parent.name, path.parent.parent.name
        if index.setdefault(sample_id, dataset_id) != dataset_id:
            duplicates.add(sample_id)
    if duplicates:
        raise RuntimeError(
            f"{len(duplicates)} sample id(s) appear under more than one dataset directory: "
            f"{sorted(duplicates)[:5]}. One genome must resolve to one GFF, or the sequence "
            "endpoint would serve whichever the glob happened to reach first."
        )
    return index


def discover_genomes(artifacts: CatalogueArtifacts) -> list[tuple[str, str | None]]:
    """`(sample_id, bakrep_dataset_id | None)` for every genome this data root should hold.

    ⛔⛔ **THE GENE UNIVERSE IS NOW A ROSTER SOURCE, WHERE THE GFF TREE USED TO BE THE ONLY ONE.**
    This function derived the entire roster by globbing `annotation_root/*/*/*.bakta.gff3.gz`.
    Where that nested tree does not exist — a cohort whose Bakta GFFs are flat
    `<BioSample>.bakta.gff3.gz`, or one served from streamed assemblies — the glob matched
    nothing, the genome layer loaded **zero genomes**, printed `0 loaded, 0 refused`, and the run
    reported success. Nothing downstream could tell that apart from a cohort that genuinely had no
    genomes, because the two produce the same output.

    ⭐ **The roster is the UNION of the two, and neither alone is right:**

    * the **gene universe** is what the catalogue *requires*. Every one of its genomes must have a
      `genome` row or `ingest_genome_collection` refuses the whole load, so a roster that does not
      cover it cannot build a catalogue at all. It is also read by `read_genome_vocabulary`, the
      same function the collection builds its ordinals from — one reader, so the genome layer and
      the membership cannot disagree about which genomes a cohort has.
    * the **annotation tree** is what this data root *holds*, and the genome layer is explicitly
      model-independent. ⚠ A genome in the tree and in no universe is exactly a **projectable**
      genome — `--stage projection` places genomes the model never clustered, and
      `projected_gene_placement` is keyed onto `gene (genome_id, flat_index)`, so it needs rows
      the catalogue's own roster will never name. Measured on the probe store: 280 GFFs, 200
      genomes across the two universes, and **10 of the other 80 are the placed genomes**.

    ⛔ An empty UNION is refused by name. That is the silent success, and it is the only case
    where loading nothing and finding nothing look the same.

    ⚠ A universe genome the tree does not hold gets `None` and is **refused by name** by
    `ingest_one_genome` rather than dropped. Missing and absent must not look the same.
    """
    roster, _ = read_genome_vocabulary(artifacts)
    index = annotation_index(artifacts)
    in_universe = set(roster)
    roster += sorted(sample_id for sample_id in index if sample_id not in in_universe)
    if not roster:
        raise EmptyRoster(
            f"neither the gene universe at {artifacts.gene_universe} nor the annotation tree at "
            f"{artifacts.annotation_root} names a genome. A roster of zero loads zero genomes and "
            "reports success, which is indistinguishable from a cohort that has none — so it is "
            "refused here instead."
        )
    return [(sample_id, index.get(sample_id)) for sample_id in roster]


def load_genome_layer(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    limit: int | None = None,
    only: set[str] | None = None,
    progress_every: int = 25,
    species_key: str | None = None,
    write_functional_annotations: bool = True,
    write_noncoding_features: bool = False,
) -> GenomeLayerReport:
    """Load every genome in the roster, one atomic transaction each.

    ``species_key`` names the species for every genome instead of reading it per parquet — see
    `ingest_one_genome`, which is where the reason lives. The two ``write_*`` switches exist so a
    table can be left out of a load deliberately and *said*, rather than arriving empty.
    """
    report = GenomeLayerReport()
    species_ids = ingest_pathogen_species(session)
    session.commit()
    if species_key is not None and species_key not in species_ids:
        raise KeyError(
            f"--species-key {species_key!r} has no `pathogen_species` row (known: "
            f"{sorted(species_ids)}). Named rather than defaulted: a load under the wrong species "
            "is a load nothing downstream can detect."
        )

    genomes = discover_genomes(artifacts)
    # ⭐ The frame, stated. Two sources make one roster, and a reader cannot otherwise tell a
    # 280-genome load of a 100-genome cohort from a 280-genome cohort.
    located = sum(1 for _, dataset_id in genomes if dataset_id is not None)
    in_universe, _ = read_genome_vocabulary(artifacts)
    report.notes.append(
        f"roster: {len(genomes):,} genomes — {len(in_universe):,} in the gene universe, "
        f"{len(genomes) - len(in_universe):,} in the annotation tree only; "
        f"{located:,} located under a dataset directory"
    )
    if only:
        genomes = [pair for pair in genomes if pair[0] in only]
    if limit:
        genomes = genomes[:limit]

    started = time.time()
    for index, (sample_id, dataset_id) in enumerate(genomes, start=1):
        try:
            one = ingest_one_genome(
                session, artifacts, sample_id=sample_id, bakrep_dataset_id=dataset_id,
                pathogen_species_ids=species_ids, species_key=species_key,
                write_functional_annotations=write_functional_annotations,
                write_noncoding_features=write_noncoding_features,
            )
            session.commit()
        # ⚠ `ValueError` is in the tuple because `pyarrow.lib.ArrowInvalid` IS one: a single
        # unreadable parquet used to escape this handler and abort the whole run, which is the
        # opposite of what a per-genome unit of work is for.
        except (GenomeIngestError, GeneAlignmentError, KeyError, ValueError) as error:
            # ⛔ Roll back THIS genome and carry on. A genome is the unit of work, so one refusal
            # must not discard the 200 already loaded — and it must be named, not counted.
            session.rollback()
            report.genomes_refused += 1
            report.refusals.append((sample_id, str(error).replace("\n", " ")[:300]))
            continue
        report.genomes_loaded += 1
        report.genome_ids.append(one.genome_id)
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

    ⭐ **Scoped to the genomes THIS run loaded**, which is the only frame in which the claim is
    true. Counted over the whole table it agrees with the run's own writes exactly when the table
    started empty — so the two probe species in one database, or a second cohort beside the first,
    made every reconciliation report a difference that was not one. ⚠ `any_` over an array rather
    than `IN (…)`: at 71,435 genomes an `IN` list exceeds Postgres' 65,535 bind parameters, which
    fails at the one scale this exists to serve.
    """
    scope = literal(report.genome_ids, ARRAY(Integer))
    held = {
        "genes": session.execute(
            select(func.count()).select_from(Gene).where(Gene.genome_id == any_(scope))
        ).scalar_one(),
        "contigs": session.execute(
            select(func.count()).select_from(GenomeContig)
            .where(GenomeContig.genome_id == any_(scope))
        ).scalar_one(),
        "functional_annotations": session.execute(
            select(func.count()).select_from(GeneFunctionalAnnotation)
            .where(GeneFunctionalAnnotation.genome_id == any_(scope))
        ).scalar_one(),
        "noncoding_features": session.execute(
            select(func.count()).select_from(GenomeNoncodingFeature)
            .where(GenomeNoncodingFeature.genome_id == any_(scope))
        ).scalar_one(),
        "genomes": session.execute(
            select(func.count()).select_from(Genome).where(Genome.genome_id == any_(scope))
        ).scalar_one(),
    }
    written = {
        "genes": report.genes, "contigs": report.contigs,
        "functional_annotations": report.functional_annotations,
        "noncoding_features": report.noncoding_features, "genomes": report.genomes_loaded,
    }
    return [
        f"{key}: wrote {written[key]:,} but the table holds {held[key]:,} "
        f"for the {len(report.genome_ids):,} genome(s) this run loaded"
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
    collection_label: str | None = None,
    allow_missing_strand: bool = False,
) -> list[str]:
    """Model registry → roster → run → catalogue, in the one order the foreign keys allow.

    ⚠ **The genome layer must already be loaded**, and this refuses rather than loading a partial
    roster: every arrangement's membership is a position in the genome vocabulary, so a hole in it
    would address the wrong genome and nothing downstream could detect that.
    """
    species_ids = ingest_pathogen_species(session)
    session.flush()
    species_id = species_ids[species_key]

    collection_id, roster = ingest_genome_collection(
        session,
        artifacts,
        pathogen_species_id=species_id,
        collection_key=collection_key,
        collection_label=collection_label,
    )
    session.flush()

    pangenome_id, run_report = ingest_pangenome_run(
        session,
        artifacts,
        pathogen_species_id=species_id,
        genome_collection_id=collection_id,
        genome_count=roster.genome_count,
        gene_count=roster.genes_in_universe,
        # ⚠ Zero until the loci exist; `ingest_locus_catalogue` sets the real count in the same
        # transaction, so a committed row never claims a locus count it does not hold.
        locus_count=0,
        # ⭐ How a reader addresses this catalogue: `{species}-{cohort}`, from the roster this
        # load actually wrote. ⛔ The second half used to be the MODEL key. Taking it from
        # `roster.collection_key` is what makes the key and the picker's option name the same
        # thing, and `published_cohorts` is what stops that key defaulting to the dataset token
        # (which produced `ecoli-ecoli`: legal under all three guards and nobody's cohort).
        catalogue_key=catalogue_key_for(species_key, roster.collection_key, catalogue_key),
    )
    catalogue_report = ingest_locus_catalogue(
        session,
        artifacts,
        pangenome_id=pangenome_id,
        pathogen_species_id=species_id,
        genome_id_by_ordinal=genome_id_by_ordinal(session, collection_id),
        allow_missing_strand=allow_missing_strand,
    )
    return [
        f"catalogue: {catalogue_key_for(species_key, roster.collection_key, catalogue_key)}",
        f"cohort: {roster.collection_key} — {roster.collection_label!r}, "
        f"{roster.genome_count} genomes, "
        f"{roster.genes_in_universe:,} genes in the universe",
        run_report.render(),
        catalogue_report.render(),
    ]


def reconcile_pangenome(session: Session, artifacts: CatalogueArtifacts) -> list[str]:
    """⛔ Ask the DATABASE what it holds, against what the checked-in triple says it should.

    The published counts are a fact of record (`published_catalogues.py`), so this compares the load
    against them rather than against itself.

    ⚠ **An unknown run_id is NOT a clean reconciliation**, and saying so is the whole of
    `pangenome_is_reconcilable`. Only two runs have a checked-in triple, so every cohort loaded from
    here on is unknown to this function — and an empty difference list made `main` print *"reconciled:
    every count the loader reported is what the tables hold"* having compared nothing at all. A
    reader cannot tell that from a load that was checked and agreed (`silence-is-not-success`).
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


def pangenome_is_reconcilable(artifacts: CatalogueArtifacts) -> bool:
    """Is there a checked-in triple to reconcile THIS run against?

    ⛔ The denominator for `reconcile_pangenome`'s zero. `PUBLISHED_CATALOGUES` holds the two static
    publications and nothing else, so a load of any other run — every cohort from the browser's
    cohort axis onward — has no record to be compared with. *Nothing differed* and *nothing was
    compared* must not print the same line.
    """
    from bacatlas_backend.ingest.published_catalogues import PUBLISHED_CATALOGUES

    return any(entry.run_id == artifacts.run_id for entry in PUBLISHED_CATALOGUES)


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
                        help="the browser key (`ecoli` | `kp`); defaults to the artifacts' own. "
                             "\u2b50 It now also names the species of every GENOME, instead of "
                             "reading `meta.species` per parquet \u2014 300 of 300 sampled kp71k "
                             "`_meta` files have no such column, and a constant is correct for a "
                             "roster built as `kpsc_final \u2227 species == K. pneumoniae`")  # fmt: skip
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
    parser.add_argument("--collection-label", default=None,
                        help="what the PICKER SHOWS for this cohort (`E. coli probe`, `Kp`). "
                             "\u26d4 Do NOT put the genome count in it — the page renders "
                             "`{label} — {count} genomes` from this and the roster's own count, so a "
                             "label that spelled its own number could contradict the roster it "
                             "names. Defaults to the collection key, and `--offer` then refuses: an "
                             "unlabelled cohort may be loaded but not shown")  # fmt: skip
    parser.add_argument("--catalogue-key", default=None,
                        help="how a reader addresses THIS catalogue (`ecoli-nuna5`). Defaults to "
                             "`{species}-{model}`, which is also what nuna's exporter writes as "
                             "`--dset`. ⛔ Hyphen, never underscore: the keys must stay prefix-free "
                             "or one catalogue's payloads match another's key")  # fmt: skip
    parser.add_argument("--gene-universe", type=Path, default=None,
                        help="the close_seq parquet, where the derived "
                             "`proc/mmseqs/{set}_close_seq.parquet` does not name it. \u26d4 Pass "
                             "this rather than symlinking: the set key also pins the model label, "
                             "the Pfam path, the payload path and the species key")  # fmt: skip
    parser.add_argument("--assignment", type=Path, default=None,
                        help="the step-4 assignment, where the derived "
                             "`analysis/step4_global_merge/assignments/{run_id}.tsv` does not name "
                             "it \u2014 a big run writes under `comparison/nuna/<model>/<run>/step4/`. "
                             "Either the .tsv or its .parquet sidecar")  # fmt: skip
    parser.add_argument("--no-functional-annotations", action="store_true",
                        help="\u2b50 do not write `gene_functional_annotation`. The per-node rollups "
                             "in `locus_annotation_entry` are what the page reads; the per-gene "
                             "table is 372 M rows at 71,435 genomes and nothing renders it. Left "
                             "OUT deliberately and reported, never arriving empty")  # fmt: skip
    parser.add_argument("--with-noncoding-features", action="store_true",
                        help="\u2b50 write `genome_noncoding_feature`, which is OFF by default. "
                             "\u26d4 The earlier `--no-noncoding-features` said suppressing it was "
                             "'a decision to take and state, not a default' because 'it renders a "
                             "track'. That was wrong, and the polarity followed the error: traced "
                             "2026-10-29, the table has NO reader \u2014 no service, no route, no "
                             "frontend field \u2014 only the ingest that writes it and the reconcile "
                             "that counts it. The track's intergenic slot is drawn from "
                             "`intergenic_gap`, whose lengths come from gene adjacencies and whose "
                             "features come from the `_noncoding` PARQUET, never from this table. "
                             "It is also the largest single table at 71,435 genomes (~5.8 GB), so "
                             "writing it is what now has to be asked for")  # fmt: skip
    parser.add_argument("--allow-missing-strand", action="store_true",
                        help="proceed when genomes have no `{BS}_strand.parquet`, recording the "
                             "shortfall in `omitted`. \u26d4 Refused by default: `load_gene_coords` "
                             "fills an absent sidecar with `+` for that genome's WHOLE gene set, "
                             "which un-reverses its windows, points every track arrow forward and "
                             "pulls syntenic_A5 toward 0.5 with no error anywhere")  # fmt: skip
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
        gene_universe_path=args.gene_universe, assignment_path=args.assignment,
    )
    optional = artifacts.verify()
    print(f"artifacts: all {len(artifacts.required())} required present; optional {optional}")

    engine = create_engine(args.database_url, future=True)
    report, differences, reconcilable = None, [], True
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
                session, artifacts, limit=args.limit, only=set(args.only) if args.only else None,
                species_key=args.species_key,
                write_functional_annotations=not args.no_functional_annotations,
                write_noncoding_features=args.with_noncoding_features,
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
                collection_label=args.collection_label,
                allow_missing_strand=args.allow_missing_strand,
            ):
                print(line)
            session.commit()
            reconcilable = pangenome_is_reconcilable(artifacts)
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
                # ⛔ CAUGHT, exactly as `--publish` below catches it. `offer_catalogue` refuses a
                # cohort still wearing its key, and that refusal is the DESIGNED way an operator
                # learns the label is missing — so it has to arrive as a line they can read, not as
                # a traceback that kills the process twelve lines after `session.commit()` already
                # wrote the load. Uncaught, the load succeeded, the catalogue is loaded and
                # addressable, and the exit status says the whole run failed.
                try:
                    print(offer_catalogue(session, catalogue_key=just_loaded))
                    session.commit()
                except PublishRefused as error:
                    session.rollback()
                    print(f"OFFER REFUSED: {error}", file=sys.stderr)
                    differences.append("the catalogue was loaded but not offered in the picker")
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
    elif reconcilable:
        print("reconciled: every count the loader reported is what the tables hold")
    else:
        # ⛔ NOT "reconciled". This run has no checked-in triple, so nothing was compared — and the
        # line that used to print here claimed every count had been checked.
        print(
            "NOT reconciled: this run_id has no row in `published_catalogues.py`, so the load was "
            "compared against nothing. The loader's own counts are above; they are unchecked."
        )
    return 1 if (report and report.refusals) else 0


if __name__ == "__main__":
    raise SystemExit(main())
