"""Genome, its assemblies and its contigs — **model-independent** facts that outlive every run.

⛔ This is the half of the gene layer that a re-ingest must never rewrite. Coordinates are a
property of the genome; a locus is a property of the model. The split is the same one that keeps
`.nseq` and `.loci` as two files, one level up.
"""

from __future__ import annotations

from sqlalchemy import BigInteger, Boolean, ForeignKey, Index, Integer, SmallInteger, String, Text, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import SampleIdentifierKind


class Genome(Base):
    """One assembled genome, ever. Not scoped to a collection or a model."""

    __tablename__ = "genome"
    __table_args__ = (
        UniqueConstraint("sample_id"),
        Index("ix_genome__bare_accession", "bare_accession"),
    )

    genome_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    pathogen_species_id: Mapped[int] = mapped_column(
        ForeignKey("pathogen_species.pathogen_species_id"), nullable=False
    )

    #: ⭐ **The store key.** Every parquet, `.npy`, GFF and assignment file is named by this.
    #: ⚠ The parquets' own `species` column is NOT the species key: it reads `ecoli | kpneumoniae`
    #: (the only two values over all 280 genomes), while `pathogen_species.species_key` is
    #: `ecoli | kp` (`published.tsv`). A direct string join drops every Klebsiella genome silently.
    sample_id: Mapped[str] = mapped_column(String(128), nullable=False)
    sample_id_kind: Mapped[SampleIdentifierKind] = mapped_column(nullable=False)

    #: The second key form, and ONLY for RefSeq/GenBank stems. ⚠ `GCF_000512165.1` against
    #: `GCF_000512165.1_ASM51216v1_genomic`. The tail is **not always** `_ASM…_genomic` — NCBI
    #: writes the submitter's strain name there, so `_INF156_genomic` and `_KSB2_2B_genomic` (two
    #: underscores) occur. Populated by `build_genome_set.bare_accession`, which anchors on the
    #: ACCESSION and not on the tail; an `_ASM`-only rule silently lost 27 genomes.
    #: NULL for BioSample-keyed short-read assemblies, which have one form.
    bare_accession: Mapped[str | None] = mapped_column(String(128), nullable=True)

    #: The BakRep dataset directory (`gff/<datasetID>/<BS>/`). No parquet carries it, so it is
    #: captured from the path at ingest or it is lost.
    bakrep_dataset_id: Mapped[str | None] = mapped_column(String(64), nullable=True)

    sublineage: Mapped[str | None] = mapped_column(String(64), nullable=True)
    lincode: Mapped[str | None] = mapped_column(String(64), nullable=True)
    phylogroup: Mapped[str | None] = mapped_column(String(64), nullable=True)

    #: ⚠ NULL means "not a Klebsiella genome / flag not applicable", NOT false.
    is_kpsc_final_list: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    is_complete_assembly: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    is_hybrid_assembly: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    #: Contigs that carry at least one KEPT CDS — the set `contig_index` enumerates, and strictly
    #: smaller than the assembly's contig count. Measured on SAMEA103923484: **270 of 348**.
    contig_count_with_coding_sequence: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: ⚠ **Every contig in the assembly**, including the 78 that carry no CDS — this is a property
    #: of the genome, not of the annotation. So it is NOT the sum of the contigs `contig_index`
    #: enumerates, and the two differ on every genome measured (5,025,405 vs 4,880,150 bases on
    #: SAMEA103923484). Taken from the GFF's `##FASTA` block, which carries them all.
    total_base_count: Mapped[int | None] = mapped_column(BigInteger, nullable=True)
    #: The CDS the extractor KEPT — the meta parquet's row count, not the GFF's CDS line count.
    #: The GFF has 6–18 more (mean 11.2 over 25 genomes): pseudo, missing-contig, empty-translation
    #: and internal-stop features are dropped BEFORE `flat_index` is numbered.
    coding_gene_count: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: ⛔ **Whether this genome's strands were READ FROM AN ANNOTATION, or defaulted to `+`.**
    #: ⚠ NOT "does a `{BS}_strand.parquet` exist" — that was the original wording and it is wrong
    #: now that strand is taken from GFF column 7, which is present on 100 % of CDS lines. Under
    #: the old reading the entire 13,573-genome Klebsiella cohort would load as `False` (it has no
    #: strand parquets) while its strand came from its GFFs and was fully observed. The flag exists
    #: for the forced-`+` path — `reference_annotation.load_gene_coords(ignore_strand=True)`, and
    #: `extract_strand.strand_for_genome`, which is permitted to default up to 20 genes or 1 % of a
    #: genome — so it must record what the INGEST did, not which file happened to be on disk.
    #: ⚠ It is honest only at GENOME granularity: the parquet's `strand` column has no nulls and
    #: `n_defaulted` lives only in an in-memory `df.attrs`, so there is no per-gene observation flag
    #: to be had from that source at all.
    strand_is_observed: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    #: ⛔ **Why this genome has NO `genome_genes` row — or NULL, meaning it has one.** A modelled
    #: genome can be clustered, assigned and counted in every per-locus number and still hold no
    #: verified gene table: nuna's `.nseq` acceptance gate refuses a genome whose stored proteins do
    #: not reproduce from its assembly (17 of the 71,435-genome *K. pneumoniae* cohort). Such a genome
    #: is NAMED at ingest (`--genomes-without-gene-table`), never inferred from an absent file, and the
    #: reason is stored here so the publish gate and the Sequence tab can tell *named as carrying none*
    #: from *the genome layer was never loaded for it* — the second is still a defect.
    #: ⚠ Every route that writes a `genome_genes` row sets this back to NULL, so a genome that later
    #: gains a table loses the mark. A row with both is refused by `verify_pangenome_is_servable`.
    gene_table_omitted_reason: Mapped[str | None] = mapped_column(Text, nullable=True)

    def __repr__(self) -> str:
        return f"<Genome {self.sample_id}>"


class GenomeAssembly(Base):
    """Where a genome's assembly came from — accession, source and the file we hold."""

    __tablename__ = "genome_assembly"
    __table_args__ = (UniqueConstraint("genome_id", "assembly_accession"),)

    genome_assembly_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), nullable=False)

    #: ⚠ Every probe assembly is an **AllTheBacteria** assembly, not an NCBI one: all 281 manifest
    #: rows read `source = atb` and `ncbi_accession` is blank on all of them. So this is nullable
    #: and its absence is a fact about the cohort, not a gap in the load.
    assembly_accession: Mapped[str | None] = mapped_column(String(64), nullable=True)
    assembly_source: Mapped[str | None] = mapped_column(String(32), nullable=True)
    assembly_file_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    annotation_file_path: Mapped[str | None] = mapped_column(String(512), nullable=True)


class GenomeContig(Base):
    """`contig_index` ↔ seqid ↔ contig name and length — the map, stored because it is measured.

    ⛔ **`contig_index` is not the contig number.** It is a position over contigs that have a CDS.
    Index 269 on SAMEA103923484 is `contig00324`. `contig_index + 1` is right for the low indices
    and silently wrong in the tail.

    ⛔ **Two GFF seqids can map to one `contig_index`, and it is not a broken join.** A short-read
    assembler emits byte-identical duplicate contigs and the extractor keeps only one, so the
    coordinate map assigns both seqids to it. That is why the map is stored as measured at ingest
    rather than re-derived per request — and why the unique constraint is on (genome, seqid) and
    NOT on (genome, contig_index).
    """

    __tablename__ = "genome_contig"
    __table_args__ = (
        UniqueConstraint("genome_id", "seqid"),
        Index("ix_genome_contig__genome_id__contig_index", "genome_id", "contig_index"),
    )

    genome_contig_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), nullable=False)

    #: The GFF's own column-1 identifier, e.g. `SAMEA103923484.contig00001`.
    seqid: Mapped[str] = mapped_column(String(128), nullable=False)
    contig_index: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    contig_name: Mapped[str] = mapped_column(String(128), nullable=False)
    length_bases: Mapped[int] = mapped_column(Integer, nullable=False)
