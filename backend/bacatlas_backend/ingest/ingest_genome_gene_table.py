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
| no | no, and NAMED as having no gene table | the `genome` row only — no block, no contigs | — |
| no | no | nothing — refused by name | — |

⛔ **The fourth row is a named allowlist, never a default.** nuna's `.nseq` acceptance gate refuses a
genome whose stored proteins do not reproduce from its assembly — 17 of the 71,435-genome cohort.
They are still MODELLED genomes: the clustering, the assignment and every per-locus number include
them, so the catalogue needs their `genome` row and their gene count. The caller names them
(`--genomes-without-gene-table`); this then reads only `flat_index` from the `_meta` parquet, requires
it to be exactly `0..n-1`, stores `coding_gene_count = n` and the reason in
`genome.gene_table_omitted_reason`, and writes no block. A genome with neither file that is NOT named
is still refused by name, and a named genome that HAS a GFF or a `.nseq` is refused — the name claims
the gate refused a genome that has a table, so the list is stale.

The third route is the 71,435-genome cohort, which holds no GFFs for the gene layer. ⛔ **No route
writes per-gene rows**: the per-gene tables were dropped (D8, migration `013b739e04ef`), and the
per-gene annotations reach the page as locus rollups built at the pangenome stage from the parquets.
The `_meta` parquet is still required — it is what the GFF and the `.nseq` are checked against — and
on the third route the species must be named, because that cohort's `_meta` files carry no
`species` column.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from sqlalchemy import func, select
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
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
from bacatlas_backend.models.pangenome import Pangenome

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


class PrefetchedReads:
    """One genome's input files, read OFF the serial SQL path — each a value, or the error its read raised.

    ⭐ **A cache keyed on exactly what was read, never a second reader.** The key is the path AND the
    columns (or the sample a `.nseq` must name), and a consumer asking for anything else reads it
    itself, so a prefetch that guessed wrong costs a read and can never hand back a different frame.
    ⛔ A failed read is stored, not raised: it surfaces when `ingest_one_genome` reaches the point
    where it would have read that file, so a genome is refused by the same check, with the same
    message, in the same place as before — and a genome refused EARLIER for another reason never
    sees it.
    """

    def __init__(self) -> None:
        self._results: dict = {}

    def run(self, key, loader, *args) -> None:
        """Read ``loader(*args)`` now, and keep its value — or the error it raised — under ``key``."""
        try:
            self._results[key] = (True, loader(*args))
        except Exception as error:  # noqa: BLE001 — re-raised, unchanged, where the consumer reads it
            self._results[key] = (False, error)

    def take(self, key, loader, *args):
        """The value kept under ``key`` (raising its error, if the read failed), or a fresh ``loader(*args)``."""
        if key in self._results:
            succeeded, value = self._results.pop(key)
            if succeeded:
                return value
            raise value
        return loader(*args)


def _meta_columns(*, has_gff: bool, species_key: str | None) -> tuple[str, ...]:
    """The `_meta` columns a route reads — ONE definition, for the prefetch and the route alike."""
    if has_gff:
        # `protein_sequence` is one of `check_against_meta`'s four independent witnesses.
        columns = ("flat_index", "contig_idx", "start", "end", "protein_sequence")
        return columns + ("species",) if species_key is None else columns
    return ("flat_index", "contig_idx", "start", "end")


#: The named route reads only this — `0..n-1` is checked and `n` stored.
_NAMED_ROUTE_META_COLUMNS = ("flat_index",)


def _load_meta(meta_path: Path, columns: tuple[str, ...]):
    import pandas  # the `ingest` extra, deliberately not a serving dependency

    return pandas.read_parquet(meta_path, columns=list(columns))


def _load_nseq(path: Path, sample_id: str) -> NseqFile:
    return read_nseq_gene_table(path, expect_sample=sample_id)


def _read_meta(meta_path: Path, columns: tuple[str, ...], prefetched: PrefetchedReads | None):
    if prefetched is None:
        return _load_meta(meta_path, columns)
    return prefetched.take(("meta", meta_path, columns), _load_meta, meta_path, columns)


def read_genome_inputs(
    artifacts: CatalogueArtifacts,
    *,
    sample_id: str,
    bakrep_dataset_id: str | None,
    species_key: str | None = None,
    nseq_path: Path | None = None,
    gene_table_omitted_reason: str | None = None,
) -> PrefetchedReads:
    """Read what `ingest_one_genome` would read for this genome — and nothing it would refuse before reading.

    Safe on a worker thread: no session, no shared state, every failure captured per file. The route
    is decided by the same four facts `ingest_one_genome` decides it by, and the columns come from
    `_meta_columns`, so the two cannot drift apart; where they did, the consumer would read for itself.
    """
    reads = PrefetchedReads()
    meta_path = artifacts.per_genome(sample_id, "meta")
    has_gff = bakrep_dataset_id is not None
    if gene_table_omitted_reason is not None:
        if not has_gff and nseq_path is None:
            columns = _NAMED_ROUTE_META_COLUMNS
            reads.run(("meta", meta_path, columns), _load_meta, meta_path, columns)
        return reads
    if not has_gff and (nseq_path is None or species_key is None):
        return reads
    columns = _meta_columns(has_gff=has_gff, species_key=species_key)
    reads.run(("meta", meta_path, columns), _load_meta, meta_path, columns)
    if nseq_path is not None:
        reads.run(("nseq", nseq_path, sample_id), _load_nseq, nseq_path, sample_id)
    return reads


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
    #: ⭐ Non-None exactly when this genome was loaded WITHOUT a gene table, named by the caller —
    #: the reason stored on `genome.gene_table_omitted_reason`. Then `gene_count` is the `_meta`
    #: parquet's row count and nothing was written to `genome_genes`.
    gene_table_omitted_reason: str | None = None
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
    gene_table_omitted_reason: str | None = None,
    prefetched: PrefetchedReads | None = None,
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

    ⭐ ``gene_table_omitted_reason`` non-None is the caller NAMING this genome as one that holds no
    gene table, and why (module table, fourth row). It is honoured only where both files are absent;
    a named genome with a GFF or a `.nseq` is refused, because the name is then stale.

    ⭐ ``prefetched`` is this genome's files already read (`read_genome_inputs`, on a thread pool),
    and then the caller has ALSO established from one directory scan that the `_meta` parquet exists
    — so the per-genome `is_file()` is not asked again. A file that vanished in between fails its read
    with an `OSError`, which the caller refuses by name. Absent, this reads for itself, as it always has.
    """
    meta_path = artifacts.per_genome(sample_id, "meta")
    if prefetched is None and not meta_path.is_file():
        raise GenomeIngestError(f"{sample_id}: no meta parquet at {meta_path}")

    # ⭐ The route, decided before a byte is read — so the `.nseq`-only route never asks the parquet
    # for the protein strings it has no use for.
    has_gff = bakrep_dataset_id is not None
    if gene_table_omitted_reason is not None:
        if has_gff or nseq_path is not None:
            held = " and ".join(
                name for name, present in (("a GFF", has_gff), ("a .nseq", nseq_path is not None)) if present
            )
            raise GenomeIngestError(
                f"{sample_id}: named as holding no gene table ({gene_table_omitted_reason}), but this "
                f"data root holds {held} for it. ⛔ The list is stale — it claims the gate refused a "
                "genome that has a table — so it is refused rather than either half being believed."
            )
        return _ingest_named_without_gene_table(
            session, meta_path, sample_id=sample_id, species_key=species_key,
            pathogen_species_ids=pathogen_species_ids, reason=gene_table_omitted_reason,
            prefetched=prefetched,
        )
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

    columns = _meta_columns(has_gff=has_gff, species_key=species_key)
    meta = _read_meta(meta_path, columns, prefetched).sort_values("flat_index")
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
            found = _read_nseq_or_refuse(nseq_path, sample_id, prefetched)
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
        found = _read_nseq_or_refuse(nseq_path, sample_id, prefetched)
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
        _contig_rows(genome_id, contig_names, contig_lengths, sample_id=sample_id),
        where=GenomeContig.__table__.c.genome_id == genome_id,
    )

    _refuse_to_replace_a_gene_table_a_catalogue_reads(
        session, sample_id=sample_id, genome_id=genome_id, gene_count=gene_count, gene_spans=gene_spans
    )
    _, report.gene_spans_written = replace_rows_for(
        session, GenomeGenes.__table__, GENOME_GENES_COLUMNS,
        [(genome_id, gene_count, SPAN_ENCODING_NSEQ_GENE_TABLE, gene_spans)],
        where=GenomeGenes.__table__.c.genome_id == genome_id,
    )
    report.gene_spans_source = source
    return report


def _refuse_to_replace_a_gene_table_a_catalogue_reads(
    session: Session, *, sample_id: str, genome_id: int, gene_count: int | None, gene_spans: bytes | None
) -> None:
    """Refuse to change — or delete — a gene table that a loaded catalogue addresses by position.

    ⛔ **ONE GENOME, ONE GENE TABLE, and the schema cannot hold two.** `genome` is unique on
    `sample_id` and `genome_genes` is keyed by the genome alone. A catalogue's
    `genome_locus_assignment` block is positional: element *i* is the locus of gene *i* **of the
    table it was built against**. Replace that table with another store's and every catalogue built
    on the old one names the wrong gene at every position after the first difference — with every
    length CHECK satisfied and nothing on the page looking wrong.

    ⚠ **Measured, not hypothetical.** 63 genomes sit in both the 100-genome Kp probe and the
    71,435-genome cohort, and the two stores' gene tables differ in **all 63** (SAMN03892119: 5,341
    genes against 5,344 — one store keeps the internal-stop and selenoprotein genes). Loading the
    second cohort's genome stage over the first would have overwritten 63 tables in silence.

    ⭐ **Keyed to the cause, so it never blocks a legitimate reload.** The same bytes are an
    idempotent re-run (a resumed load) and pass. A genome no catalogue reads is free to change. Only
    the combination that breaks something is refused: *different bytes* under *a catalogue that
    reads them by position*. The digest is taken in the database, so a resume does not pull every
    stored block through the socket to compare it.

    ``gene_spans=None`` is the named no-gene-table route, which DELETES the block: the same harm,
    refused on the same terms.
    """
    held = session.execute(
        select(GenomeGenes.gene_count, func.sha256(GenomeGenes.gene_spans)).where(GenomeGenes.genome_id == genome_id)
    ).one_or_none()
    if held is None:
        return
    held_count, held_digest = held
    same_bytes = gene_spans is not None and bytes(held_digest) == hashlib.sha256(gene_spans).digest()
    if same_bytes and held_count == gene_count:
        return
    readers = (
        session.execute(
            select(Pangenome.catalogue_key)
            .join(GenomeLocusAssignment, GenomeLocusAssignment.pangenome_id == Pangenome.pangenome_id)
            .where(GenomeLocusAssignment.genome_id == genome_id)
            .order_by(Pangenome.catalogue_key)
        )
        .scalars()
        .all()
    )
    if not readers:
        return
    if gene_spans is None:
        offered = "no gene table (it is named as holding none)"
    else:
        offered = f"a different one ({gene_count:,} genes)"
    raise GenomeIngestError(
        f"{sample_id}: this genome already holds a gene table ({held_count:,} genes) and this load offers "
        f"{offered}. ⛔ Refused: {len(readers)} loaded catalogue(s) address the stored table by position "
        f"({', '.join(readers[:5])}{', …' if len(readers) > 5 else ''}), and a genome can hold only one. "
        "Two cohorts whose stores annotate one assembly differently cannot share a database; to "
        "change the table, remove those catalogues first and load them again afterwards."
    )


def _ingest_named_without_gene_table(
    session: Session,
    meta_path: Path,
    *,
    sample_id: str,
    species_key: str | None,
    pathogen_species_ids: dict[str, int],
    reason: str,
    prefetched: PrefetchedReads | None = None,
) -> GenomeIngestReport:
    """The fourth route: a `genome` row, the gene count, the reason — and NO block.

    ⛔ The gene count is the `_meta` parquet's, and `flat_index` must be exactly `0..n-1`: that count
    is what `genome_locus_assignment.gene_count` is written at and what every membership's
    `flat_index` is range-checked against, so a parquet with a hole would size the block wrong.
    ⚠ The species must be NAMED, exactly as on the `.nseq`-only route — the cohort this exists for
    has no `species` column to read one from.

    ⚠ Any `genome_genes` / `genome_contig` row left by an earlier load is DELETED: a stale block would
    otherwise sit beside the mark, and the publish gate refuses a genome carrying both.
    """
    import numpy  # the `ingest` extra, function-local like every numpy import in this package

    if species_key is None:
        raise GenomeIngestError(
            f"{sample_id}: a genome named as holding no gene table needs `species_key` named. The "
            "cohort it exists for carries no `species` column in its `_meta` parquets."
        )
    if species_key not in pathogen_species_ids:
        raise GenomeIngestError(
            f"{sample_id}: species {species_key!r} has no row. Run ingest_pathogen_species first."
        )
    flat_index = _read_meta(meta_path, _NAMED_ROUTE_META_COLUMNS, prefetched)["flat_index"].to_numpy()
    gene_count = len(flat_index)
    if gene_count == 0:
        raise GenomeIngestError(f"{sample_id}: the meta parquet is empty")
    if not numpy.array_equal(numpy.sort(flat_index), numpy.arange(gene_count)):
        raise GenomeIngestError(
            f"{sample_id}: the meta parquet's flat_index is not exactly 0..{gene_count - 1} "
            f"({gene_count:,} rows), so its row count is not this genome's gene count. ⛔ Refused: "
            "the assignment block is sized from it."
        )
    genome_id = _upsert_genome(
        session, sample_id=sample_id, species_id=pathogen_species_ids[species_key],
        bakrep_dataset_id=None, annotation_path=None,
        contig_count=None, total_base_count=None, coding_gene_count=gene_count,
        # ⚠ False: nothing was read from which a strand could be observed. No reader reaches it —
        # the strand is served only from a gene table, which this genome does not have.
        strand_is_observed=False,
        gene_table_omitted_reason=reason,
    )
    _refuse_to_replace_a_gene_table_a_catalogue_reads(
        session, sample_id=sample_id, genome_id=genome_id, gene_count=None, gene_spans=None
    )
    for table in (GenomeGenes.__table__, GenomeContig.__table__):
        session.execute(table.delete().where(table.c.genome_id == genome_id))
    session.flush()
    return GenomeIngestReport(
        sample_id=sample_id, genome_id=genome_id, gene_count=gene_count, gene_table_omitted_reason=reason
    )


def _read_nseq_or_refuse(path: Path, sample_id: str, prefetched: PrefetchedReads | None = None) -> NseqFile:
    """The `.nseq`'s header and gene table, with any failure named for the GENOME and not just the file.

    `expect_sample=sample_id` is what makes the header's own name do its job: a file published under
    the wrong name decodes perfectly and would serve another genome's coordinates. A read done ahead
    on a thread pool (``prefetched``) fails into the same handler as one done here.
    """
    try:
        if prefetched is None:
            return _load_nseq(path, sample_id)
        return prefetched.take(("nseq", path, sample_id), _load_nseq, path, sample_id)
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
    contig_count: int | None,
    total_base_count: int | None,
    coding_gene_count: int,
    strand_is_observed: bool = True,
    gene_table_omitted_reason: str | None = None,
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
    # ⛔ Written on EVERY route, because the NCBI sequence station asks for nothing else: an assembly
    # stem's bases are fetched by its bare accession, and a NULL there is refused at serving time as
    # "an ingest fault". NULL for a BioSample id, which has one form.
    genome.bare_accession = bare_accession_for(sample_id)
    genome.bakrep_dataset_id = bakrep_dataset_id
    genome.contig_count_with_coding_sequence = contig_count
    genome.total_base_count = total_base_count
    genome.coding_gene_count = coding_gene_count
    # ⛔ True because the strand was READ from GFF column 7, which every CDS line carries — not
    # because a strand parquet happens to exist. On the `.nseq` route it is the flags word, which
    # nuna set from the same GFF column when it built the file. See `Genome.strand_is_observed`.
    genome.strand_is_observed = strand_is_observed
    # ⛔ Written on EVERY route: NULL wherever a block is about to be written, so a genome that was
    # once named as holding no gene table and has since gained one loses the mark.
    genome.gene_table_omitted_reason = gene_table_omitted_reason
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


def bare_accession_for(sample_id: str) -> str | None:
    """`GCF_000512165.1_ASM51216v1_genomic` → `GCF_000512165.1`; NULL for a BioSample id.

    ⭐ nuna's `build_genome_set.bare_accession` and nothing else — the rule that anchors on the
    ACCESSION, not on the `_ASM…` tail (`_INF156_genomic` and `_KSB2_2B_genomic` occur). That function
    hands back anything it does not recognise unchanged; here an unrecognised assembly stem is NULL
    instead, because a stem is not an accession and asking NCBI for one is a request that cannot be
    answered — `verify_pangenome_is_servable` then names the genome rather than the page failing on it.
    """
    if _sample_identifier_kind(sample_id) is SampleIdentifierKind.BIOSAMPLE:
        return None
    from nuna.tl.stability.build_genome_set import bare_accession

    accession = bare_accession(sample_id)
    return None if accession == sample_id else accession


def _contig_name(seqid: str, sample_id: str) -> str:
    """The contig NAME, as distinct from the seqid — the BakRep prefix `{sample_id}.` removed, and only that.

    BakRep writes `{BioSample}.contig00001`, and the part after the prefix is what a reader recognises.
    ⛔ It used to be "everything after the first dot". That is the same string for every BakRep seqid
    — measured: the seqid on every one of the 1,568,072 feature lines of the 280 probe GFFs begins
    with its own `{BS}.` — and the VERSION for an NCBI accession: `NZ_CP003999.1` was stored, and
    served, as contig `1`. A seqid without the prefix is its own name, which is the honest answer
    rather than a fragment of it.
    """
    prefix = f"{sample_id}."
    return seqid[len(prefix):] if seqid.startswith(prefix) else seqid


def _contig_rows(genome_id: int, seqids: Sequence[str], lengths: Sequence[int], *, sample_id: str):
    """⛔ Only contigs carrying a kept CDS — that is what `contig_index` enumerates, on both routes."""
    for contig_index, (seqid, length) in enumerate(zip(seqids, lengths, strict=True)):
        yield (genome_id, seqid, contig_index, _contig_name(seqid, sample_id), length)

