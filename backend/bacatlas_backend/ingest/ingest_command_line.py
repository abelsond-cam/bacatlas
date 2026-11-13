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
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy import ARRAY, Integer, any_, create_engine, func, literal, select
from sqlalchemy.orm import Session

import bacatlas_backend.models  # noqa: F401  (registers every table on the metadata)
from bacatlas_backend.ingest.artifact_locator import DEFAULT_READ_WORKERS, CatalogueArtifacts, read_ahead_in_order
from bacatlas_backend.ingest.genome_gene_alignment import GeneAlignmentError
from bacatlas_backend.ingest.ingest_genome_collection_roster import (
    genome_id_by_ordinal,
    ingest_genome_collection,
    read_genome_vocabulary,
)
from bacatlas_backend.ingest.ingest_genome_gene_table import (
    GENE_SPANS_SOURCES,
    GenomeIngestError,
    ingest_one_genome,
    read_genome_inputs,
)
from bacatlas_backend.ingest.ingest_locus_catalogue import ingest_locus_catalogue
from bacatlas_backend.ingest.ingest_pangenome_run import catalogue_key_for, ingest_pangenome_run
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.ingest.ingest_projected_genomes import ProjectionRefused, load_projection
from bacatlas_backend.ingest.ingest_reference_vocabularies import load_reference_vocabularies
from bacatlas_backend.ingest.publish_pangenome import PublishRefused, publish_pangenome
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome

#: ⭐ The reason stored on `genome.gene_table_omitted_reason` for every genome named by
#: `--genomes-without-gene-table`. One sentence, because the Sequence tab shows it to a reader and the
#: catalogue's `omitted_sections["gene_table"]` quotes it.
NSEQ_GATE_REFUSED_REASON = (
    "the .nseq acceptance gate refused it: its stored proteins do not reproduce from the assembly"
)


class GeneTableListRefused(ValueError):
    """`--genomes-without-gene-table` names something the data root contradicts — the list is wrong."""


class GenomeMetaMissing(ValueError):
    """Roster genomes with no `_meta` parquet — the RUN is refused, by name, before a genome is written."""


def refuse_unless_every_genome_has_meta(
    genomes: list[tuple[str, str | None]], meta_present: frozenset[str], per_genome_root: Path
) -> None:
    """⛔ Every roster genome — a named one included — needs its `_meta` parquet, so ask ONCE, up front.

    The parquet is what every route checks its coordinates against and what a named genome's gene
    count is read from. It used to be asked per genome, as a cold `is_file()` inside the serial loop,
    and a gap surfaced as one refusal hours in. Now it is one scan of the directory and a refusal of
    the whole run in seconds, naming the genomes — a cohort missing a `_meta` is a mirror that is
    wrong, not a genome that is.
    """
    missing = [sample_id for sample_id, _ in genomes if sample_id not in meta_present]
    if missing:
        raise GenomeMetaMissing(
            f"{len(missing):,} of {len(genomes):,} roster genome(s) have no `_meta.parquet` under "
            f"{per_genome_root}: {missing[:10]}. Every route reads it, so nothing is loaded."
        )


def read_genomes_without_gene_table(path: Path) -> frozenset[str]:
    """The named allowlist: one `sample_id` per line, blank lines ignored, nothing else interpreted.

    ⛔ A name repeated is refused rather than collapsed: a list written twice over is a list nobody
    checked, and the count it implies is the count the report states.
    """
    names = [line.strip() for line in Path(path).expanduser().read_text().splitlines() if line.strip()]
    repeated = sorted(name for name, seen in Counter(names).items() if seen > 1)
    if repeated:
        raise GeneTableListRefused(f"{path}: names {repeated[:5]} more than once")
    return frozenset(names)


def refuse_unless_gene_table_list_matches_roster(
    named: frozenset[str],
    genomes: list[tuple[str, str | None]],
    present: frozenset[str],
) -> None:
    """⛔ Before a genome is loaded: every named genome is in the roster and holds NEITHER file.

    The two refusals are the two ways the list can be stale — a name the roster does not contain
    (a typo, or another cohort's list) and a name the data root holds a GFF or `.nseq` for (the list
    claims the gate refused a genome that has a table). Either refuses the RUN, by name, because the
    list is one claim about the cohort and a list wrong in one place is not trusted in the others.
    """
    roster = {sample_id: dataset_id for sample_id, dataset_id in genomes}
    not_in_roster = sorted(named - set(roster))
    holds_a_table = sorted(
        sample_id for sample_id in named
        if sample_id in roster and (roster[sample_id] is not None or sample_id in present)
    )
    problems = []
    if not_in_roster:
        problems.append(f"{len(not_in_roster)} named genome(s) are not in the roster: {not_in_roster[:10]}")
    if holds_a_table:
        problems.append(
            f"{len(holds_a_table)} named genome(s) have a GFF or a .nseq in this data root, so the list "
            f"is stale — it claims the gate refused a genome that has a table: {holds_a_table[:10]}"
        )
    if problems:
        raise GeneTableListRefused("--genomes-without-gene-table: " + "; ".join(problems))


@dataclass
class GenomeLayerReport:
    """The whole genome-layer load, in numbers that can be reconciled against the store."""

    genomes_loaded: int = 0
    genomes_refused: int = 0
    #: Genes across the `genome_genes` blocks written — the sum of their `gene_count`.
    genes: int = 0
    contigs: int = 0
    #: ⭐ `genome_genes` rows written, and by which ROUTE. The three route counts plus
    #: `without_gene_table` add up to `genomes_loaded`, and a reader needs them apart: a load where
    #: every genome reads `encoded` is a load that checked no block against nuna's file, which one
    #: total could not say.
    gene_spans: int = 0
    gene_spans_routes: dict[str, int] = field(default_factory=lambda: dict.fromkeys(GENE_SPANS_SOURCES, 0))
    #: ⛔ Genomes loaded WITHOUT a gene table because they were NAMED (`--genomes-without-gene-table`)
    #: — loaded, not refused, and listed by name so the count can be checked against the list.
    without_gene_table: list[str] = field(default_factory=list)
    #: Their genes — the `_meta` row counts stored as `coding_gene_count` — kept OUT of `genes`, which
    #: is the blocks' total and is what `verify_pangenome_is_servable` adds these to.
    genes_without_gene_table: int = 0
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
            f"  gene_spans blocks      {self.gene_spans:,}  ("
            + " · ".join(f"{route} {self.gene_spans_routes.get(route, 0):,}" for route in GENE_SPANS_SOURCES)
            + ")",
            f"  no gene table (named)  {len(self.without_gene_table):,}"
            + (
                f"  ({self.genes_without_gene_table:,} genes: {', '.join(sorted(self.without_gene_table))})"
                if self.without_gene_table
                else ""
            ),
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

    * the **`.nseq` tree** (`artifacts.nseq_root`, when one is named) is the third thing this data
      root can hold, and the ONLY thing it holds for the 71,435-genome cohort, which has no GFFs for
      the gene layer. A genome there is LOCATED — `ingest_one_genome` loads it from the file — so
      it joins the roster on the same footing as a genome in the annotation tree. ⚠ Read from ONE
      directory listing (`nseq_present`), never a per-genome `exists()`.

    ⛔ An empty UNION is refused by name. That is the silent success, and it is the only case
    where loading nothing and finding nothing look the same.

    ⚠ A universe genome neither tree holds gets `None` here and is **refused by name** by
    `ingest_one_genome` rather than dropped. Missing and absent must not look the same. The
    dataset id is the GFF's; whether a `.nseq` exists is a separate fact the loader resolves from
    the same listing, so one pair still means one thing.
    """
    return _discover_genomes(artifacts)[0]


def _discover_genomes(artifacts: CatalogueArtifacts) -> tuple[list[tuple[str, str | None]], int]:
    """`discover_genomes`, and how many of its genomes are the universe's — from ONE read of the universe.

    ⚠ The genome stage used to read the universe's 374 M-row `sample_id` column a second time only to
    count its genomes (~10 s and ~13 GB transient at 71,435 genomes). The count is `len(roster)` as it
    stands before the trees' extras are appended, so it is taken there.
    """
    roster, _ = read_genome_vocabulary(artifacts)
    universe_count = len(roster)
    index = annotation_index(artifacts)
    present = artifacts.nseq_present()
    in_universe = set(roster)
    roster += sorted((set(index) | present) - in_universe)
    if not roster:
        raise EmptyRoster(
            f"neither the gene universe at {artifacts.gene_universe}, nor the annotation tree at "
            f"{artifacts.annotation_root}, nor the .nseq tree at {artifacts.nseq_root} names a "
            "genome. A roster of zero loads zero genomes and reports success, which is "
            "indistinguishable from a cohort that has none — so it is refused here instead."
        )
    return [(sample_id, index.get(sample_id)) for sample_id in roster], universe_count


def load_genome_layer(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    limit: int | None = None,
    only: set[str] | None = None,
    progress_every: int = 25,
    species_key: str | None = None,
    genomes_without_gene_table: frozenset[str] | None = None,
    read_workers: int = DEFAULT_READ_WORKERS,
) -> GenomeLayerReport:
    """Load every genome in the roster, one atomic transaction each.

    ``species_key`` names the species for every genome instead of reading it per parquet — see
    `ingest_one_genome`, which is where the reason lives.

    ``genomes_without_gene_table`` NAMES the genomes that hold neither a GFF nor a `.nseq` because
    nuna's `.nseq` acceptance gate refused them. Each is loaded as a `genome` row with no block; the
    list is checked against the whole roster first and raises `GeneTableListRefused` — before any
    genome is written — if it names a genome the roster lacks or one this data root holds a table for.

    ⭐ ``read_workers`` threads read each genome's files (the `_meta` columns, the `.nseq` gene table)
    ahead of the loop, which stays serial and takes them in ROSTER ORDER — so the rows written, their
    order and every id are what a serial read gives; only the Lustre latency moves off the SQL thread.
    ``read_workers <= 1`` reads inline. Every roster genome's `_meta` is established first by ONE scan
    (`GenomeMetaMissing` refuses the run, by name, before any genome is written).
    """
    named = genomes_without_gene_table or frozenset()
    report = GenomeLayerReport()
    species_ids = ingest_pathogen_species(session)
    session.commit()
    if species_key is not None and species_key not in species_ids:
        raise KeyError(
            f"--species-key {species_key!r} has no `pathogen_species` row (known: "
            f"{sorted(species_ids)}). Named rather than defaulted: a load under the wrong species "
            "is a load nothing downstream can detect."
        )

    genomes, universe_count = _discover_genomes(artifacts)
    # ⭐ ONE listing of the .nseq tree for the whole run — never a stat per genome.
    present = artifacts.nseq_present()
    # ⭐ The frame, stated. Three sources make one roster, and a reader cannot otherwise tell a
    # 280-genome load of a 100-genome cohort from a 280-genome cohort — nor a load that copied
    # every block from nuna's files from one that encoded every block itself.
    with_gff = sum(1 for _, dataset_id in genomes if dataset_id is not None)
    with_nseq = sum(1 for sample_id, _ in genomes if sample_id in present)
    with_both = sum(1 for sample_id, dataset_id in genomes if dataset_id is not None and sample_id in present)
    located = with_gff + with_nseq - with_both
    report.notes.append(
        f"roster: {len(genomes):,} genomes — {universe_count:,} in the gene universe, "
        f"{len(genomes) - universe_count:,} held by this data root only; "
        f"located {located:,}: {with_gff:,} with a GFF, {with_nseq:,} with a .nseq, "
        f"{with_both:,} with both, {len(genomes) - located:,} with neither — "
        f"{len(named):,} of them named as holding no gene table, "
        f"{len(genomes) - located - len(named):,} refused by name"
    )
    # ⛔ Against the WHOLE roster, before `--only`/`--limit` narrow it and before anything is written.
    refuse_unless_gene_table_list_matches_roster(named, genomes, present)
    refuse_unless_every_genome_has_meta(genomes, artifacts.per_genome_present("meta"), artifacts.per_genome_root)
    if only:
        genomes = [pair for pair in genomes if pair[0] in only]
    if limit:
        genomes = genomes[:limit]

    def inputs_of(pair: tuple[str, str | None]) -> dict:
        sample_id, dataset_id = pair
        return dict(
            sample_id=sample_id, bakrep_dataset_id=dataset_id, species_key=species_key,
            nseq_path=artifacts.per_genome_nseq(sample_id) if sample_id in present else None,
            gene_table_omitted_reason=NSEQ_GATE_REFUSED_REASON if sample_id in named else None,
        )

    def read(pair: tuple[str, str | None]):
        return read_genome_inputs(artifacts, **inputs_of(pair))

    started = time.time()
    in_order = read_ahead_in_order(genomes, read, read_workers)
    for index, ((sample_id, dataset_id), prefetched) in enumerate(in_order, start=1):
        try:
            one = ingest_one_genome(
                session, artifacts, pathogen_species_ids=species_ids, prefetched=prefetched,
                **inputs_of((sample_id, dataset_id)),
            )
            session.commit()
        # ⚠ `ValueError` is in the tuple because `pyarrow.lib.ArrowInvalid` IS one: a single
        # unreadable parquet used to escape this handler and abort the whole run, which is the
        # opposite of what a per-genome unit of work is for. ⚠ `OSError` for the same reason: pandas
        # opens the `_meta` itself, so an EIO, ESTALE or EACCES on one of 71,435 files on Lustre was a
        # traceback that ended the stage. A database error is NOT caught: it is not a fact about one
        # genome, and the session it leaves is not one to carry on in.
        except (GenomeIngestError, GeneAlignmentError, KeyError, ValueError, OSError) as error:
            # ⛔ Roll back THIS genome and carry on. A genome is the unit of work, so one refusal
            # must not discard the 200 already loaded — and it must be named, not counted.
            session.rollback()
            report.genomes_refused += 1
            report.refusals.append((sample_id, str(error).replace("\n", " ")[:300]))
            continue
        report.genomes_loaded += 1
        report.genome_ids.append(one.genome_id)
        if one.gene_table_omitted_reason is None:
            report.genes += one.gene_count
        else:
            report.without_gene_table.append(sample_id)
            report.genes_without_gene_table += one.gene_count
        report.contigs += one.contigs_written
        report.gene_spans += one.gene_spans_written
        if one.gene_spans_source is not None:
            report.gene_spans_routes[one.gene_spans_source] = report.gene_spans_routes.get(one.gene_spans_source, 0) + 1
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
        "contigs": session.execute(
            select(func.count()).select_from(GenomeContig)
            .where(GenomeContig.genome_id == any_(scope))
        ).scalar_one(),
        "genomes": session.execute(
            select(func.count()).select_from(Genome).where(Genome.genome_id == any_(scope))
        ).scalar_one(),
        "gene_spans": session.execute(
            select(func.count()).select_from(GenomeGenes).where(GenomeGenes.genome_id == any_(scope))
        ).scalar_one(),
    }
    written = {
        "contigs": report.contigs, "genomes": report.genomes_loaded, "gene_spans": report.gene_spans,
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
    genomes_without_gene_table: frozenset[str] | None = None,
    read_workers: int = DEFAULT_READ_WORKERS,
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
        genomes_without_gene_table=genomes_without_gene_table,
        read_workers=read_workers,
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
    parser.add_argument("--nseq-dir", type=Path, default=None,
                        help="where nuna published `{BS}.nseq` / `{BS}.loci`. \u2b50 A genome with a "
                             ".nseq has its gene_spans block COPIED from the file (and, where it also "
                             "has a GFF, checked byte for byte against the GFF's encoding); a genome "
                             "with a .nseq and NO GFF is loaded from the header alone. "
                             "\u26d4 Listed ONCE, never stat-ed per genome. Unset, nothing is read "
                             "from one and the block is encoded from the GFF as before")  # fmt: skip
    parser.add_argument("--allow-missing-strand", action="store_true",
                        help="proceed when genomes have no `{BS}_strand.parquet`, recording the "
                             "shortfall in `omitted`. \u26d4 Refused by default: `load_gene_coords` "
                             "fills an absent sidecar with `+` for that genome's WHOLE gene set, "
                             "which un-reverses its windows, points every track arrow forward and "
                             "pulls syntenic_A5 toward 0.5 with no error anywhere")  # fmt: skip
    parser.add_argument("--genomes-without-gene-table", type=Path, default=None,
                        help="a text file, one sample_id per line: the modelled genomes nuna's .nseq "
                             "acceptance gate refused (their stored proteins do not reproduce from the "
                             "assembly). \u2b50 `genomes`: each is loaded as a genome row with no gene "
                             "table instead of being refused. `pangenome`: the collection's genomes "
                             "without a gene table must be EXACTLY this set. \u26d4 A named genome the "
                             "roster lacks, or one with a GFF or .nseq, refuses the run. Unset, a "
                             "genome with neither file is refused by name as before")  # fmt: skip
    parser.add_argument("--read-workers", type=int, default=DEFAULT_READ_WORKERS,
                        help="threads reading each genome's `_meta` and `.nseq` (genomes) and its `.loci` "
                             "(pangenome) ahead of the serial load, in roster order. \u2b50 The rows "
                             "written are the same at any value; it only moves Lustre's per-file "
                             "latency off the SQL thread. 0 reads inline")  # fmt: skip
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

    # ⛔ `species_key` MUST reach the artifacts: the roster resolves its cohort from
    # `(artifacts.species_key, artifacts.set_key)`, and left out it defaults to the set key — so the
    # kp71k launcher's `--set-key kp71k --species-key kp` looked up `('kp71k', 'kp71k')`, which no
    # cohort declares, and `--stage pangenome` died with `UndeclaredCohort` after the whole genome
    # stage had run. The probe never saw it: there the two keys are the same string.
    artifacts = CatalogueArtifacts(
        data_root=args.data_root, set_key=args.set_key,
        model_label=args.model_label, run_id=args.run_id,
        species_key=args.species_key,
        gene_universe_path=args.gene_universe, assignment_path=args.assignment,
        nseq_root=args.nseq_dir,
    )
    optional = artifacts.verify()
    print(f"artifacts: all {len(artifacts.required())} required present; optional {optional}")

    try:
        without_gene_table = (
            read_genomes_without_gene_table(args.genomes_without_gene_table)
            if args.genomes_without_gene_table is not None
            else None
        )
    except (OSError, GeneTableListRefused) as error:
        print(f"--genomes-without-gene-table: {error}", file=sys.stderr)
        return 2

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
            try:
                report = load_genome_layer(
                    session, artifacts, limit=args.limit, only=set(args.only) if args.only else None,
                    species_key=args.species_key, genomes_without_gene_table=without_gene_table,
                    read_workers=args.read_workers,
                )
            except (GeneTableListRefused, GenomeMetaMissing) as error:
                # ⛔ Raised before any genome is written: the list is a misconfiguration, not a refusal
                # of one genome, so nothing is loaded and the exit code says so.
                session.rollback()
                print(f"GENOME LAYER REFUSED: {error}", file=sys.stderr)
                return 2
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
                genomes_without_gene_table=without_gene_table,
                read_workers=args.read_workers,
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
