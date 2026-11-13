"""One genome at a time: its `genome` / assembly / contig rows and its `genome_genes` block.

⛔ **This is the model-INDEPENDENT half.** Coordinates are a property of the genome and outlive
every run; a locus is a property of a model. A re-ingest of a new pangenome must not rewrite these
rows.

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
| yes | yes | the file's block; the GFF's encoding must equal it byte for byte | `nseq+encoded_checked` |
| yes | no | the block encoded from the GFF | `encoded` |
| no | yes | the genome layer from the `.nseq` header, the block copied verbatim | `nseq` |
| no | no | nothing — refused by name | — |

The third route is the 71,435-genome cohort, which holds no GFFs for the gene layer. ⛔ **No route
writes per-gene rows**: the per-gene tables were dropped (D8, migration `013b739e04ef`), and the
per-gene annotations reach the page as locus rollups built at the pangenome stage from the parquets.
The `_meta` parquet is still required — it is what the GFF and the `.nseq` are checked against — and
on the third route the species must be named, because that cohort's `_meta` files carry no
`species` column.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

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
from bacatlas_backend.models.genome import Genome, GenomeAssembly, GenomeContig
from bacatlas_backend.models.genome_genes import GenomeGenes

CONTIG_COLUMNS = ("genome_id", "seqid", "contig_index", "contig_name", "length_bases")

#: One row per genome: the whole gene table as a block. ⚠ `gene_spans` is `bytes`, never `str` —
#: `copy_rows` refuses a `str` in a bytea column because bytea's input function reinterprets it.
GENOME_GENES_COLUMNS = ("genome_id", "gene_count", "span_encoding", "gene_spans")

#: Where a genome's `gene_spans` block came from — the three routes above, named so a layer report
#: can say how many genomes took each rather than one count a GFF-only load and a `.nseq`-only load
#: would both satisfy.
GeneSpansSource = Literal["nseq+encoded_checked", "encoded", "nseq"]
GENE_SPANS_SOURCES: tuple[GeneSpansSource, ...] = ("nseq+encoded_checked", "encoded", "nseq")


class GenomeIngestError(RuntimeError):
    """A genome could not be loaded, named with the reason and the check that caught it."""


@dataclass
class GenomeIngestReport:
    """What one genome's load actually did — counts, not a boolean."""

    sample_id: str
    genome_id: int
    #: Genes in the `genome_genes` block written — its `gene_count`, on every route.
    gene_count: int = 0
    contigs_written: int = 0
    #: ⭐ 0 or 1, and WHERE the block came from. A report that said only "1 block" could not tell a
    #: load that checked every block against nuna's file from one that read no file at all.
    gene_spans_written: int = 0
    gene_spans_source: GeneSpansSource | None = None
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

    ⭐ **Public, with `split_set_value` above it, because these two DEFINE the parsed shape of a
    functional annotation value** read from the dbxref parquet. `ingest/catalogue_frames.py` builds
    every per-locus rollup through them, so one definition fixes what a value is.
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
) -> GenomeIngestReport:
    """Load one genome: its `genome` row, assembly, contigs and gene-span block — or raise before writing.

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
        # `protein_sequence` is one of `check_against_meta`'s four independent witnesses.
        columns = ["flat_index", "contig_idx", "start", "end", "protein_sequence"]
        if species_key is None:
            columns.append("species")
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
    report = GenomeIngestReport(sample_id=sample_id, genome_id=genome_id, gene_count=gene_count)

    _, report.contigs_written = replace_rows_for(
        session, GenomeContig.__table__, CONTIG_COLUMNS,
        _contig_rows(genome_id, contig_names, contig_lengths),
        where=GenomeContig.__table__.c.genome_id == genome_id,
    )

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

