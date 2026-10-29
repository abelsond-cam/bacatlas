"""One genome at a time: its contigs, its genes, and their functional labels.

⛔ **This is the model-INDEPENDENT half.** Coordinates are a property of the genome and outlive
every run; a locus is a property of a model. A re-ingest of a new pangenome must not rewrite these
rows — 412 M of them at 80,000 genomes.

⭐ **The unit of work is one genome, and it is atomic.** A load that fails on genome 200 of 280 must
be resumable without dropping the database, and a re-run after a fix must not double a row. Each
genome is delete-then-COPY inside the caller's transaction, so a crash between the two leaves the
previous rows intact and a re-run replaces exactly what it replaced before.

⚠ **Nothing is written until the alignment gate passes.** `genome_gene_alignment.check_against_meta`
verifies the numbering against the parquet on four independent things; a genome that fails is
refused whole rather than written partially.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy.orm import Session

from bacatlas_backend.gff.gff_cds_parser import parse_genome_annotation
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.genome_gene_alignment import (
    GenomeAlignment,
    check_against_meta,
    flatten_to_extractor_order,
)
from bacatlas_backend.ingest.ingest_pathogen_species import species_key_for_parquet_value
from bacatlas_backend.ingest.staging_table_loader import replace_rows_for
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.gene import Gene, GeneFunctionalAnnotation, GenomeNoncodingFeature
from bacatlas_backend.models.genome import Genome, GenomeAssembly, GenomeContig

#: Columns written per gene, in COPY order.
GENE_COLUMNS = (
    "genome_id", "flat_index", "pathogen_species_id", "contig_index",
    "start_position", "end_position", "length_nt", "strand", "strand_is_observed",
    "gff_phase", "is_five_prime_partial", "protein_length_aa",
    "bakta_gene_symbol", "bakta_product", "locus_tag",
)

FUNCTIONAL_COLUMNS = (
    "genome_id", "flat_index", "uniref50_accession", "kegg_orthology_id",
    "cog_id", "cog_categories", "ec_numbers", "gene_ontology_terms",
)

CONTIG_COLUMNS = ("genome_id", "seqid", "contig_index", "contig_name", "length_bases")

NONCODING_COLUMNS = (
    "genome_id", "contig_index", "noncoding_feature_index", "start_position", "end_position",
    "strand", "feature_type", "feature_gene", "feature_product", "overlaps_coding_sequence",
)

#: ⛔ The Bakta locus-tag shape — `nuna.eval.assignments.strip_locus_tags`' regex, and the reason a
#: tag must never reach `bakta_gene_symbol`: tags are unique per genome, so a family holding only
#: them would report one distinct symbol and read as unanimous.
LOCUS_TAG_SHAPE = re.compile(r"^[A-Z0-9]{4,}_\d+$")

#: ⛔ `ec` and `go` are comma-joined SORTED-UNIQUE SETS, not scalars. Measured over 40 genomes:
#: 12.3 % of non-null `ec` values carry more than one term, up to 4 terms and 39 characters.
SET_VALUED_DBXREF_COLUMNS = ("ec", "go")


class GenomeIngestError(RuntimeError):
    """A genome could not be loaded, named with the reason and the check that caught it."""


@dataclass
class GenomeIngestReport:
    """What one genome's load actually did — counts, not a boolean."""

    sample_id: str
    genome_id: int
    genes_written: int = 0
    contigs_written: int = 0
    functional_annotations_written: int = 0
    noncoding_features_written: int = 0
    genes_deleted: int = 0
    #: ⚠ Coverage per label, so *not annotated* can be told from *not measured*. A column absent
    #: from the parquet arrives as `None` here and is ANNOUNCED, never counted as zero coverage.
    label_coverage: dict[str, int | None] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _sample_identifier_kind(sample_id: str) -> SampleIdentifierKind:
    """⚠ BioSample accessions have ONE form; RefSeq/GenBank assembly stems have two."""
    return (
        SampleIdentifierKind.BIOSAMPLE
        if sample_id.startswith(("SAMN", "SAMEA", "SAMD"))
        else SampleIdentifierKind.ASSEMBLY_STEM
    )


def split_set_value(value) -> list[str] | None:
    """A comma-joined sorted-unique set → a list. ⛔ NULL stays NULL; `[]` is never manufactured.

    *Not annotated* and *annotated with nothing* are different claims, and the parquets contain no
    empty strings — so an empty list here could only have been invented by this function.
    """
    if value is None or value != value:  # NaN compares unequal to itself
        return None
    text = str(value).strip()
    return text.split(",") if text else None


def _split_cog_categories(value) -> list[str] | None:
    """`CP` → `['C', 'P']`. ⛔ A concatenated SET of single letters, not one category.

    Measured over 40 probe genomes: 14,812 of 129,865 non-null values (11.4 %) carry more than one,
    up to 4. ⚠ NULL stays NULL — *not annotated* is not *annotated with nothing*.
    """
    if value is None or value != value:
        return None
    text = str(value).strip()
    return list(text) if text else None


def scalar_or_none(value):
    """A pandas cell → a Python scalar or None, without importing pandas into this signature.

    ⭐ **Public, with `split_set_value` above it, because these two DEFINE the stored shape of a
    functional annotation value** — what `gene_functional_annotation` holds for one gene, from the
    one parquet column it holds it from. `ingest/catalogue_frames.py` reads the same parquet through
    nuna's annotation frame and must reach the same values, or a rollup computed there and a rate
    measured over the gene rows differ by a parse rather than by a fact. One definition, two
    callers.
    """
    if value is None or value != value:
        return None
    return value


def ingest_one_genome(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    sample_id: str,
    bakrep_dataset_id: str | None,
    pathogen_species_ids: dict[str, int],
    species_key: str | None = None,
    write_functional_annotations: bool = True,
    write_noncoding_features: bool = False,
) -> GenomeIngestReport:
    """Load one genome's contigs, genes and labels. Raises before writing anything if it cannot.

    ``pathogen_species_ids`` maps `species_key` → id, from `ingest_pathogen_species`.

    ⭐ **``species_key`` names the species instead of reading it, and then ``species`` is not
    requested from the parquet at all.** The probe store's `_meta` files carry that column; **300
    of 300 sampled kp71k files do not** — they share a 10-column schema without it — so a load of
    that cohort failed on a `KeyError` from pyarrow before it read a single gene. ⚠ A constant is
    *correct* there and not a convenience: `selected_genomes_kp_final` is
    `kpsc_final_list ∧ species == "Klebsiella pneumoniae"` exactly, and its report JSON records the
    6,618 non-*pneumoniae* genomes excluded when the roster was built. Absent, the parquet is read
    as before, so the probe path is untouched.

    ⚠ ``bakrep_dataset_id`` may be ``None`` — a roster genome the annotation tree does not hold.
    It is refused **by name** rather than dropped, because the roster is the gene universe now and
    a genome missing from it is a different fact than a genome missing a file.
    """
    import pandas  # the `ingest` extra, deliberately not a serving dependency

    meta_path = artifacts.per_genome(sample_id, "meta")
    if not meta_path.is_file():
        raise GenomeIngestError(f"{sample_id}: no meta parquet at {meta_path}")

    # `gene_name` is read for the symbol FALLBACK only — see `_read_products`. It is never the
    # primary source, because it carries the locus tag where there is no symbol.
    columns = ["flat_index", "contig_idx", "start", "end", "protein_sequence", "gene_name"]
    if species_key is None:
        columns.insert(5, "species")
    meta = pandas.read_parquet(meta_path, columns=columns).sort_values("flat_index")
    if meta.empty:
        raise GenomeIngestError(f"{sample_id}: the meta parquet is empty")

    if species_key is None:
        species_key = species_key_for_parquet_value(str(meta["species"].iloc[0]))
    if species_key not in pathogen_species_ids:
        raise GenomeIngestError(
            f"{sample_id}: species {species_key!r} has no row. Run ingest_pathogen_species first."
        )
    species_id = pathogen_species_ids[species_key]

    if bakrep_dataset_id is None:
        raise GenomeIngestError(
            f"{sample_id}: in the roster but under no dataset directory in "
            f"{artifacts.annotation_root} — no `<datasetID>/<BS>/<BS>.bakta.gff3.gz`"
        )
    annotation_path = artifacts.genome_annotation(sample_id, bakrep_dataset_id)
    if not annotation_path.is_file():
        raise GenomeIngestError(f"{sample_id}: no GFF at {annotation_path}")

    alignment = flatten_to_extractor_order(parse_genome_annotation(annotation_path))
    check_against_meta(
        alignment, sample_id, meta["flat_index"], meta["contig_idx"],
        meta["start"], meta["end"], meta["protein_sequence"],
    )

    genome_id = _upsert_genome(
        session, sample_id=sample_id, species_id=species_id, alignment=alignment,
        bakrep_dataset_id=bakrep_dataset_id, annotation_path=annotation_path,
    )
    report = GenomeIngestReport(sample_id=sample_id, genome_id=genome_id)

    _, report.contigs_written = replace_rows_for(
        session, GenomeContig.__table__, CONTIG_COLUMNS,
        _contig_rows(genome_id, alignment),
        where=GenomeContig.__table__.c.genome_id == genome_id,
    )

    products = _read_products(pandas, artifacts, sample_id, meta)
    report.genes_deleted, report.genes_written = replace_rows_for(
        session, Gene.__table__, GENE_COLUMNS,
        _gene_rows(genome_id, species_id, alignment, meta, products),
        where=Gene.__table__.c.genome_id == genome_id,
    )

    if write_functional_annotations:
        frame, coverage, note = _read_dbxrefs(pandas, artifacts, sample_id, len(alignment.genes))
        report.label_coverage = coverage
        if note:
            report.notes.append(note)
        _, report.functional_annotations_written = replace_rows_for(
            session, GeneFunctionalAnnotation.__table__, FUNCTIONAL_COLUMNS,
            _functional_rows(genome_id, frame),
            where=GeneFunctionalAnnotation.__table__.c.genome_id == genome_id,
        )

    if write_noncoding_features:
        _, report.noncoding_features_written = replace_rows_for(
            session, GenomeNoncodingFeature.__table__, NONCODING_COLUMNS,
            _noncoding_rows(pandas, artifacts, sample_id, genome_id),
            where=GenomeNoncodingFeature.__table__.c.genome_id == genome_id,
        )

    return report


def _upsert_genome(
    session: Session,
    *,
    sample_id: str,
    species_id: int,
    alignment: GenomeAlignment,
    bakrep_dataset_id: str,
    annotation_path: Path,
) -> int:
    """One row per genome, ever. Re-ingest refreshes the measured counts and nothing else."""
    from sqlalchemy import select

    genome = session.execute(select(Genome).where(Genome.sample_id == sample_id)).scalar_one_or_none()
    if genome is None:
        genome = Genome(sample_id=sample_id, pathogen_species_id=species_id)
        session.add(genome)

    genome.sample_id_kind = _sample_identifier_kind(sample_id)
    genome.bakrep_dataset_id = bakrep_dataset_id
    genome.contig_count_with_coding_sequence = len(alignment.contig_names_in_index_order)
    genome.total_base_count = alignment.total_base_count
    genome.coding_gene_count = len(alignment.genes)
    # ⛔ True because the strand was READ from GFF column 7, which every CDS line carries — not
    # because a strand parquet happens to exist. See `Genome.strand_is_observed`.
    genome.strand_is_observed = True
    session.flush()

    assembly = next(
        (row for row in session.query(GenomeAssembly).filter_by(genome_id=genome.genome_id)), None
    )
    if assembly is None:
        assembly = GenomeAssembly(genome_id=genome.genome_id)
        session.add(assembly)
    # ⚠ `assembly_accession` and `assembly_file_path` stay NULL here, and that is a fact about what
    # is on this machine, not a gap in the load: the probe assembly manifest is not mirrored
    # locally, and every one of its 281 rows reads `source = atb` with a blank `ncbi_accession`
    # anyway. The annotation path IS known, and it is the one the sequence endpoint opens.
    assembly.annotation_file_path = str(annotation_path)
    session.flush()
    return genome.genome_id


def _contig_rows(genome_id: int, alignment: GenomeAlignment):
    """⛔ Only contigs carrying a kept CDS — that is what `contig_index` enumerates."""
    for contig_index, seqid in enumerate(alignment.contig_names_in_index_order):
        yield (
            genome_id,
            seqid,
            contig_index,
            # The contig NAME, as distinct from the seqid: BakRep writes `{BioSample}.contig00001`,
            # and the part after the dot is what a reader recognises. Where there is no dot the two
            # are the same string, which is the honest answer rather than an empty column.
            seqid.split(".", 1)[1] if "." in seqid else seqid,
            alignment.contig_lengths_by_seqid[seqid],
        )


def _read_products(pandas, artifacts: CatalogueArtifacts, sample_id: str, meta):
    """`flat_index → (symbol, product)`, from `product.gene` with a TAG-FILTERED `gene_name` fallback.

    ⛔ **`product.gene` first, and never `meta.gene_name` unfiltered.** `gene_name` falls back to the
    genome-private Bakta locus tag where there is no symbol — 809 of 4,876 rows (16.6 %) on
    SAMEA103923484 — and a per-genome tag in a cross-genome symbol column makes a family holding
    only those report one distinct symbol and read as unanimous.

    ⭐ **But `product.gene` alone loses real symbols, and the fallback is what recovers them.**
    `bakta_products.join_products` drops every `(contig_idx, start, end)` key that is not unique in
    the genome (`drop_duplicates(keep=False)`), so an ambiguous coordinate yields NaN in the product
    parquet while `meta.gene_name` still carries the symbol. Measured over 60 probe genomes: **25
    real symbols** present in `gene_name` and absent from `product.gene` — `perC`, `yfdN`, `dicB`,
    `traW`, `traC` — against 249,797 that both carry. Taking `product.gene` as an unqualified rule
    silently discards them.

    ⚠ The fallback is admitted ONLY where the value is not locus-tag-shaped, so it recovers symbols
    without ever admitting a tag. Where both sources have a real symbol they agree on every row.
    """
    path = artifacts.per_genome(sample_id, "product")
    fallback = {
        int(row.flat_index): value
        for row in meta.itertuples(index=False)
        if (value := scalar_or_none(getattr(row, "gene_name", None))) and not LOCUS_TAG_SHAPE.fullmatch(str(value))
    }
    if not path.is_file():
        return {index: (symbol, None) for index, symbol in fallback.items()}
    frame = pandas.read_parquet(path, columns=["flat_index", "gene", "product"])
    out = {}
    for row in frame.itertuples(index=False):
        index = int(row.flat_index)
        symbol = scalar_or_none(row.gene)
        out[index] = (symbol if symbol is not None else fallback.get(index), scalar_or_none(row.product))
    return out


def _gene_rows(genome_id: int, species_id: int, alignment: GenomeAlignment, meta, products: dict):
    protein_lengths = [len(str(value)) for value in meta["protein_sequence"]]
    for gene in alignment.genes:
        symbol, product = products.get(gene.flat_index, (None, None))
        yield (
            genome_id, gene.flat_index, species_id, gene.contig_index,
            gene.start_position, gene.end_position, gene.length_nt,
            gene.strand, True, gene.phase, gene.is_five_prime_partial,
            # ⚠ From the parquet's own protein string, not from `(end-start+1)/3 - 1`. That
            # arithmetic is a coordinate-convention CHECK and holds here only because `phase` is 0
            # on every CDS in this cohort; a non-zero phase breaks it and not the measurement.
            protein_lengths[gene.flat_index],
            symbol, product, gene.locus_tag,
        )


def _read_dbxrefs(pandas, artifacts: CatalogueArtifacts, sample_id: str, gene_count: int):
    """The six label columns, with per-label coverage and an ANNOUNCEMENT for any absent column.

    ⛔ A column missing from the file arrives as `None` coverage, **not 0**. *Not measured* and
    *zero coverage* are different findings, and files written before a column existed are exactly
    the case that makes them look alike.
    """
    wanted = ["flat_index", "uniref50", "kegg_ko", "cog_id", "cog_cat", "ec", "go"]
    path = artifacts.per_genome(sample_id, "dbxref")
    if not path.is_file():
        return None, dict.fromkeys(wanted[1:]), f"{sample_id}: no dbxref parquet — no label was measured"

    frame = pandas.read_parquet(path)
    absent = [name for name in wanted[1:] if name not in frame.columns]
    frame = frame.reindex(columns=wanted)
    coverage = {
        name: (None if name in absent else int(frame[name].notna().sum())) for name in wanted[1:]
    }
    note = None
    if absent:
        note = (
            f"{sample_id}: {absent} absent from the dbxref parquet — recorded as NOT MEASURED, "
            f"which is not the same as zero coverage over {gene_count:,} genes"
        )
    return frame, coverage, note


def _functional_rows(genome_id: int, frame):
    if frame is None:
        return
    for row in frame.itertuples(index=False):
        yield (
            genome_id, int(row.flat_index),
            scalar_or_none(row.uniref50), scalar_or_none(row.kegg_ko),
            scalar_or_none(row.cog_id), _split_cog_categories(row.cog_cat),
            split_set_value(row.ec), split_set_value(row.go),
        )


def _noncoding_rows(pandas, artifacts: CatalogueArtifacts, sample_id: str, genome_id: int):
    """tRNA, rRNA, ncRNA, CRISPR — keyed `(genome, contig, nc_index)` and never by `flat_index`.

    ⛔ Admitting a non-coding row into the gene numbering would renumber every gene and shift every
    ±5 neighbourhood, which is why the parquet deliberately carries no `flat_index` at all.
    """
    path = artifacts.per_genome(sample_id, "noncoding")
    if not path.is_file():
        return
    frame = pandas.read_parquet(path)
    for row in frame.itertuples(index=False):
        yield (
            genome_id, int(row.contig_idx), int(row.nc_index),
            int(row.start), int(row.end), scalar_or_none(row.strand),
            str(row.ftype), scalar_or_none(row.gene), scalar_or_none(row.product),
            bool(row.overlaps_cds) if scalar_or_none(row.overlaps_cds) is not None else None,
        )
