"""One genome at a time: its contigs, its gene-span block, and — where a GFF is held — its genes and labels.

⛔ **This is the model-INDEPENDENT half.** Coordinates are a property of the genome and outlive
every run; a locus is a property of a model. A re-ingest of a new pangenome must not rewrite these
rows — 412 M of them at 80,000 genomes.

⭐ **The unit of work is one genome, and it is atomic.** A load that fails on genome 200 of 280 must
be resumable without dropping the database, and a re-run after a fix must not double a row. Each
genome is delete-then-COPY inside the caller's transaction, so a crash between the two leaves the
previous rows intact and a re-run replaces exactly what it replaced before.

⚠ **Nothing is written until the gate passes.** On the GFF routes `genome_gene_alignment.check_against_meta`
verifies the numbering against the parquet on four independent things; on the `.nseq` routes the
file's gene table is checked against the GFF's encoding, or against the parquet, before a row is
written. A genome that fails is refused whole rather than written partially.

⭐ **Three routes, decided per genome from the two facts the roster already resolved — a GFF under a
dataset directory, and a `.nseq` in the named tree (`CatalogueArtifacts.nseq_root`):**

| GFF | `.nseq` | written | `gene_spans_source` |
|---|---|---|---|
| yes | yes | today's load; the GFF's encoding must equal the file's block byte for byte | `nseq+encoded_checked` |
| yes | no | today's load, the block encoded from the GFF | `encoded` |
| no | yes | the genome layer from the `.nseq` header, the block copied verbatim, NO per-gene rows | `nseq` |
| no | no | nothing — refused by name | — |

The third route is the 71,435-genome cohort, which holds no GFFs for the gene layer. ⛔ On it `gene`,
`gene_functional_annotation` and `genome_noncoding_feature` are not written: those tables are being
dropped, and the per-gene annotations reach the page as locus rollups built at the pangenome stage
from the parquets. The `_meta` parquet is still required — it is what the `.nseq` is checked
against — and the species must be named, because that cohort's `_meta` files carry no `species` column.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from sqlalchemy import delete
from sqlalchemy.orm import Session

from bacatlas_backend.gene_blocks import (
    GENE_SPAN_STRIDE,
    SPAN_ENCODING_NSEQ_GENE_TABLE,
    GeneBlockError,
    GeneSpanTable,
    NseqFile,
    encode_gene_spans,
    pack_gene_flags,
    read_nseq_gene_table,
)
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
from bacatlas_backend.models.genome_genes import GenomeGenes

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

#: One row per genome: the whole gene table as a block. ⚠ `gene_spans` is `bytes`, never `str` —
#: `copy_rows` refuses a `str` in a bytea column because bytea's input function reinterprets it.
GENOME_GENES_COLUMNS = ("genome_id", "gene_count", "span_encoding", "gene_spans")

#: Where a genome's `gene_spans` block came from — the three routes above, named so a layer report
#: can say how many genomes took each rather than one count a GFF-only load and a `.nseq`-only load
#: would both satisfy.
GeneSpansSource = Literal["nseq+encoded_checked", "encoded", "nseq"]
GENE_SPANS_SOURCES: tuple[GeneSpansSource, ...] = ("nseq+encoded_checked", "encoded", "nseq")

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
    #: ⭐ 0 or 1, and WHERE the block came from. A report that said only "1 block" could not tell a
    #: load that checked every block against nuna's file from one that read no file at all.
    gene_spans_written: int = 0
    gene_spans_source: GeneSpansSource | None = None
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


def gene_spans_from_alignment(alignment: GenomeAlignment) -> bytes:
    """The `gene_spans` block a GFF implies — `encode_gene_spans` over the alignment, in `flat_index` order.

    ⭐ Public because it is one half of the differential: the 200-genome test encodes every probe GFF
    through this and compares the bytes to nuna's `.nseq`, and the dual route refuses a genome where
    the two disagree. The flags word is `pack_gene_flags` and nothing else — the same bit layout nuna
    writes, from the same three GFF facts.
    """
    genes = alignment.genes
    return encode_gene_spans(
        [gene.start_position for gene in genes],
        [gene.end_position for gene in genes],
        [gene.contig_index for gene in genes],
        [
            pack_gene_flags(
                is_minus_strand=gene.strand == "-",
                is_five_prime_partial=gene.is_five_prime_partial,
                gff_phase=gene.phase,
            )
            for gene in genes
        ],
    )


def ingest_one_genome(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    sample_id: str,
    bakrep_dataset_id: str | None,
    pathogen_species_ids: dict[str, int],
    species_key: str | None = None,
    nseq_path: Path | None = None,
    write_functional_annotations: bool = True,
    write_noncoding_features: bool = False,
) -> GenomeIngestReport:
    """Load one genome: contigs, the gene-span block and, with a GFF, its genes and labels — or raise before writing.

    ``pathogen_species_ids`` maps `species_key` → id, from `ingest_pathogen_species`.

    ⭐ **``species_key`` names the species instead of reading it, and then ``species`` is not
    requested from the parquet at all.** The probe store's `_meta` files carry that column; **300
    of 300 sampled kp71k files do not** — they share a 10-column schema without it — so a load of
    that cohort failed on a `KeyError` from pyarrow before it read a single gene. ⚠ A constant is
    *correct* there and not a convenience: `selected_genomes_kp_final` is
    `kpsc_final_list ∧ species == "Klebsiella pneumoniae"` exactly, and its report JSON records the
    6,618 non-*pneumoniae* genomes excluded when the roster was built. Absent, the parquet is read
    as before, so the probe path is untouched — except on the `.nseq`-only route, which REQUIRES it.

    ⚠ ``bakrep_dataset_id`` may be ``None`` — a roster genome the annotation tree does not hold —
    and ``nseq_path`` may be ``None`` — a genome the `.nseq` tree does not hold. The caller resolved
    both from one walk each (`discover_genomes`, `nseq_present`); this never stats a tree to find
    out. ⛔ Both ``None`` is refused **by name** rather than dropped, because the roster is the gene
    universe now and a genome missing from it is a different fact than a genome missing a file.
    """
    import pandas  # the `ingest` extra, deliberately not a serving dependency

    meta_path = artifacts.per_genome(sample_id, "meta")
    if not meta_path.is_file():
        raise GenomeIngestError(f"{sample_id}: no meta parquet at {meta_path}")

    # ⭐ The route, decided before a byte is read — so the `.nseq`-only route never asks the parquet
    # for the protein strings it has no use for.
    has_gff = bakrep_dataset_id is not None
    if not has_gff and nseq_path is None:
        raise GenomeIngestError(
            f"{sample_id}: in the roster but under no dataset directory in "
            f"{artifacts.annotation_root} — no `<datasetID>/<BS>/<BS>.bakta.gff3.gz` — and in no "
            "`.nseq` tree either (`--nseq-dir`). Nothing holds this genome's coordinates."
        )
    if not has_gff and species_key is None:
        raise GenomeIngestError(
            f"{sample_id}: the `.nseq`-only route needs `species_key` named. The cohort it exists for "
            "carries no `species` column in its `_meta` parquets, so there is nothing to read one from."
        )

    if has_gff:
        # `gene_name` is read for the symbol FALLBACK only — see `_read_products`. It is never the
        # primary source, because it carries the locus tag where there is no symbol.
        columns = ["flat_index", "contig_idx", "start", "end", "protein_sequence", "gene_name"]
        if species_key is None:
            columns.insert(5, "species")
    else:
        columns = ["flat_index", "contig_idx", "start", "end"]
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

    source: GeneSpansSource
    if has_gff:
        annotation_path = artifacts.genome_annotation(sample_id, bakrep_dataset_id)
        if not annotation_path.is_file():
            raise GenomeIngestError(f"{sample_id}: no GFF at {annotation_path}")

        alignment = flatten_to_extractor_order(parse_genome_annotation(annotation_path))
        check_against_meta(
            alignment, sample_id, meta["flat_index"], meta["contig_idx"],
            meta["start"], meta["end"], meta["protein_sequence"],
        )
        encoded = gene_spans_from_alignment(alignment)
        contig_names = alignment.contig_names_in_index_order
        contig_lengths = tuple(alignment.contig_lengths_by_seqid[seqid] for seqid in contig_names)
        gene_count = len(alignment.genes)
        if nseq_path is not None:
            # ⭐ The V1 oracle, at ingest: nuna's encoder and this one must agree byte for byte, or
            # the genome is refused. The block written is the FILE's — nuna's encoder is the only
            # encoder — and the two are identical by this check, so which is copied is not a choice.
            found = _read_nseq_or_refuse(nseq_path, sample_id)
            _refuse_unless_nseq_matches_gff(sample_id, found, contig_names, contig_lengths, encoded)
            gene_spans, source = found.gene_table, "nseq+encoded_checked"
        else:
            gene_spans, source = encoded, "encoded"
        genome_id = _upsert_genome(
            session, sample_id=sample_id, species_id=species_id,
            bakrep_dataset_id=bakrep_dataset_id, annotation_path=annotation_path,
            contig_count=len(contig_names), total_base_count=alignment.total_base_count,
            coding_gene_count=gene_count,
        )
    else:
        found = _read_nseq_or_refuse(nseq_path, sample_id)
        _refuse_unless_nseq_matches_meta(sample_id, found, meta)
        gene_spans, source = found.gene_table, "nseq"
        contig_names, contig_lengths, gene_count = found.contig_names, found.contig_lengths, found.gene_count
        # ⚠ `total_base_count` is `n_bases` here: the sum over the contigs that CARRY A CDS, because
        # those are the only contigs a `.nseq` holds. The GFF route stores the whole `##FASTA`
        # extent, which is larger on every genome with a CDS-free contig (SAMEA103923484: 348
        # contigs, 270 with a CDS). Nothing reads the column today; the frame is named here so that
        # the day something does, a genome loaded by this route is not read as a smaller assembly.
        genome_id = _upsert_genome(
            session, sample_id=sample_id, species_id=species_id,
            bakrep_dataset_id=None, annotation_path=None,
            contig_count=len(contig_names), total_base_count=found.base_count,
            coding_gene_count=gene_count,
        )
    report = GenomeIngestReport(sample_id=sample_id, genome_id=genome_id)

    _, report.contigs_written = replace_rows_for(
        session, GenomeContig.__table__, CONTIG_COLUMNS,
        _contig_rows(genome_id, contig_names, contig_lengths),
        where=GenomeContig.__table__.c.genome_id == genome_id,
    )

    if has_gff:
        products = _read_products(pandas, artifacts, sample_id, meta)
        report.genes_deleted, report.genes_written = replace_rows_for(
            session, Gene.__table__, GENE_COLUMNS,
            _gene_rows(genome_id, species_id, alignment, meta, products),
            where=Gene.__table__.c.genome_id == genome_id,
        )

        if write_functional_annotations:
            frame, coverage, note = _read_dbxrefs(pandas, artifacts, sample_id, gene_count)
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
    else:
        # ⛔ Replace-by-genome holds on this route too. A genome loaded from its GFF once and from
        # its `.nseq` now must not keep per-gene rows nothing wrote this time — three deletes that
        # usually match nothing, so a re-run replaces exactly what it replaced before.
        for table in (GeneFunctionalAnnotation.__table__, GenomeNoncodingFeature.__table__):
            session.execute(delete(table).where(table.c.genome_id == genome_id))
        report.genes_deleted = session.execute(
            delete(Gene.__table__).where(Gene.__table__.c.genome_id == genome_id)
        ).rowcount

    _, report.gene_spans_written = replace_rows_for(
        session, GenomeGenes.__table__, GENOME_GENES_COLUMNS,
        [(genome_id, gene_count, SPAN_ENCODING_NSEQ_GENE_TABLE, gene_spans)],
        where=GenomeGenes.__table__.c.genome_id == genome_id,
    )
    report.gene_spans_source = source
    return report


def _read_nseq_or_refuse(path: Path, sample_id: str) -> NseqFile:
    """The `.nseq`'s header and gene table, with any failure named for the GENOME and not just the file.

    `expect_sample=sample_id` is what makes the header's own name do its job: a file published under
    the wrong name decodes perfectly and would serve another genome's coordinates.
    """
    try:
        return read_nseq_gene_table(path, expect_sample=sample_id)
    except (GeneBlockError, OSError, KeyError, ValueError) as error:
        raise GenomeIngestError(f"{sample_id}: the .nseq at {path} could not be read — {error}") from error


def _span_row(table: GeneSpanTable, flat_index: int) -> str:
    if flat_index >= len(table):
        return "no such row"
    return (
        f"start {table.start_position(flat_index)} end {table.end_position(flat_index)} "
        f"contig {table.contig_index(flat_index)} flags {table.flags(flat_index):#x}"
    )


def _refuse_unless_nseq_matches_gff(
    sample_id: str,
    found: NseqFile,
    contig_names: Sequence[str],
    contig_lengths: Sequence[int],
    encoded: bytes,
) -> None:
    """⛔ Byte-identical, or refused — naming WHAT differed, never a warning.

    The `.nseq` is nuna's encoding of this genome and `encoded` is ours from the same GFF. They are
    supposed to be the same bytes; where they are not, one encoder is wrong about this genome, and
    serving either would put coordinates under the wrong gene with nothing on any page looking odd.
    """
    differences: list[str] = []
    encoded_count = len(encoded) // GENE_SPAN_STRIDE
    if found.gene_count != encoded_count:
        differences.append(f"gene count: the .nseq holds {found.gene_count:,}, the GFF aligned {encoded_count:,}")
    if tuple(found.contig_names) != tuple(contig_names):
        first = _first_difference(found.contig_names, contig_names)
        differences.append(
            f"contig names: {len(found.contig_names)} in the .nseq vs {len(contig_names)} from the GFF, "
            f"first difference at contig_index {first}: "
            f"{_at(found.contig_names, first)!r} vs {_at(contig_names, first)!r}"
        )
    if tuple(found.contig_lengths) != tuple(contig_lengths):
        first = _first_difference(found.contig_lengths, contig_lengths)
        differences.append(
            f"contig lengths: first difference at contig_index {first}: "
            f"{_at(found.contig_lengths, first)!r} vs {_at(contig_lengths, first)!r}"
        )
    if found.gene_table != encoded:
        theirs = GeneSpanTable(found.gene_table)
        ours = GeneSpanTable(encoded)
        first = next(
            (
                j for j in range(min(len(theirs), len(ours)))
                if found.gene_table[j * GENE_SPAN_STRIDE:(j + 1) * GENE_SPAN_STRIDE]
                != encoded[j * GENE_SPAN_STRIDE:(j + 1) * GENE_SPAN_STRIDE]
            ),
            min(len(theirs), len(ours)),
        )
        differences.append(
            f"gene table: first difference at flat_index {first} — "
            f".nseq {_span_row(theirs, first)} vs GFF {_span_row(ours, first)}"
        )
    if differences:
        raise GenomeIngestError(
            f"{sample_id}: the .nseq at {found.path} does not match the GFF — " + "; ".join(differences)
            + ". ⛔ Refused whole: nuna's encoder and this one disagree about this genome, and whichever "
            "is wrong, serving either would put coordinates under the wrong gene."
        )


def _first_difference(ours: Sequence, theirs: Sequence) -> int:
    for index, (a, b) in enumerate(zip(ours, theirs, strict=False)):
        if a != b:
            return index
    return min(len(ours), len(theirs))


def _at(values: Sequence, index: int):
    return values[index] if index < len(values) else None


def _refuse_unless_nseq_matches_meta(sample_id: str, found: NseqFile, meta) -> None:
    """⛔ The `.nseq`-only route's gate: the file's gene table against the parquet every artifact is keyed on.

    No GFF means no `check_against_meta`, so this is the one independent witness the route has:
    `start`, `end` and `contig_idx` for every gene, vectorised, and `flat_index` as `0..n-1` so row
    *j* of the table IS gene *j*. A stale or mis-built `.nseq` would put every gene's coordinates
    under the parquet's numbering with nothing downstream able to see it.
    """
    import numpy  # the `ingest` extra, function-local like every numpy import in this package

    gene_count = len(meta)
    if found.gene_count != gene_count:
        raise GenomeIngestError(
            f"{sample_id}: the .nseq at {found.path} holds {found.gene_count:,} genes but the meta "
            f"parquet has {gene_count:,} rows. ⛔ Refused whole: row j of the gene table must be "
            "flat_index j, and a count that differs means it is not."
        )
    flat_index = meta["flat_index"].to_numpy()
    if not numpy.array_equal(flat_index, numpy.arange(gene_count)):
        raise GenomeIngestError(
            f"{sample_id}: the meta parquet's flat_index is not 0..{gene_count - 1} — the gene table "
            "is positional and cannot be aligned to it"
        )
    table = numpy.frombuffer(found.gene_table, dtype="<u4").reshape(gene_count, 4)
    for name, word in (("start", 0), ("end", 1), ("contig_idx", 2)):
        parquet = meta[name].to_numpy().astype("int64")
        in_file = table[:, word].astype("int64")
        bad = numpy.flatnonzero(parquet != in_file)
        if bad.size:
            first = int(bad[0])
            raise GenomeIngestError(
                f"{sample_id}: the .nseq at {found.path} disagrees with the meta parquet on `{name}` "
                f"for {bad.size:,} of {gene_count:,} genes — first at flat_index {first}: parquet "
                f"{int(parquet[first])}, .nseq {int(in_file[first])}. ⛔ Refused whole."
            )
    if gene_count and int(table[:, 2].max()) >= len(found.contig_names):
        raise GenomeIngestError(
            f"{sample_id}: the .nseq gene table names contig index {int(table[:, 2].max())} but its "
            f"header lists {len(found.contig_names)} contigs"
        )


def _upsert_genome(
    session: Session,
    *,
    sample_id: str,
    species_id: int,
    bakrep_dataset_id: str | None,
    annotation_path: Path | None,
    contig_count: int,
    total_base_count: int,
    coding_gene_count: int,
) -> int:
    """One row per genome, ever. Re-ingest refreshes the measured counts and nothing else.

    ⚠ Takes the measured values, not the thing they were measured from, because two routes measure
    them from different things — a GFF alignment, or a `.nseq` header — and one row has to come out.
    """
    from sqlalchemy import select

    genome = session.execute(select(Genome).where(Genome.sample_id == sample_id)).scalar_one_or_none()
    if genome is None:
        genome = Genome(sample_id=sample_id, pathogen_species_id=species_id)
        session.add(genome)

    genome.sample_id_kind = _sample_identifier_kind(sample_id)
    genome.bakrep_dataset_id = bakrep_dataset_id
    genome.contig_count_with_coding_sequence = contig_count
    genome.total_base_count = total_base_count
    genome.coding_gene_count = coding_gene_count
    # ⛔ True because the strand was READ from GFF column 7, which every CDS line carries — not
    # because a strand parquet happens to exist. On the `.nseq` route it is the flags word, which
    # nuna set from the same GFF column when it built the file. See `Genome.strand_is_observed`.
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
    # anyway. The annotation path IS known where a GFF was read, and is NULL where none was.
    assembly.annotation_file_path = None if annotation_path is None else str(annotation_path)
    session.flush()
    return genome.genome_id


def _contig_rows(genome_id: int, seqids: Sequence[str], lengths: Sequence[int]):
    """⛔ Only contigs carrying a kept CDS — that is what `contig_index` enumerates, on both routes."""
    for contig_index, (seqid, length) in enumerate(zip(seqids, lengths, strict=True)):
        yield (
            genome_id,
            seqid,
            contig_index,
            # The contig NAME, as distinct from the seqid: BakRep writes `{BioSample}.contig00001`,
            # and the part after the dot is what a reader recognises. Where there is no dot the two
            # are the same string, which is the honest answer rather than an empty column.
            seqid.split(".", 1)[1] if "." in seqid else seqid,
            length,
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
