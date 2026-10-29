"""The gene layer, split into what belongs to the GENOME and what belongs to the MODEL.

⛔ **`gene` is model-independent and `gene_locus_membership` is not.** Coordinates are a property of
the genome and outlive every run; a locus is a property of the model. Folding them into one table
would mean rewriting 412 M coordinate rows every time a model changes. This is the same split that
keeps a genome's DNA and its locus index as two artifacts, one level up.

⚠ Everything the Sequence tab shows is a column here EXCEPT the base letters themselves, which are
parsed on demand from the original gzipped Bakta GFF — see `bacatlas_backend.gff`.
"""

from __future__ import annotations

from sqlalchemy import (
    BigInteger,
    Boolean,
    CheckConstraint,
    ForeignKey,
    Index,
    Integer,
    SmallInteger,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base


class Gene(Base):
    """One coding sequence in one genome. Keyed `(genome_id, flat_index)`.

    ⛔ **`flat_index` is a POSITIONAL row offset**, shared by the embedding `.npy`, every
    `{BS}_*.parquet` and the assignment files. Nothing here may renumber it: a DB that did would
    break all three at once. It is a running counter over CDS lines and **cannot be inferred from
    coordinates**.
    """

    __tablename__ = "gene"
    __table_args__ = (
        # ⭐ BRIN, not btree, on genome_id: ingest loads one genome at a time so the table is
        # physically clustered by genome, and a BRIN over 412 M rows is a few hundred kB against
        # ~13 GB for a btree.
        Index("ix_gene__genome_id", "genome_id", postgresql_using="brin"),
        Index("ix_gene__genome_id__contig_index__start_position", "genome_id", "contig_index", "start_position"),
        CheckConstraint("start_position >= 1 AND end_position >= start_position", name="span_is_forward"),
        CheckConstraint("strand IN ('+', '-')", name="strand_is_plus_or_minus"),
    )

    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)
    flat_index: Mapped[int] = mapped_column(Integer, primary_key=True)
    pathogen_species_id: Mapped[int] = mapped_column(SmallInteger, nullable=False)

    #: ⛔ Order of first appearance over contigs that HAVE A CDS — not the contig number.
    contig_index: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    #: 1-based inclusive, the GFF's own numbers.
    start_position: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Inclusive, and **includes the stop codon** — so protein length is `(end-start+1)/3 - 1`.
    end_position: Mapped[int] = mapped_column(Integer, nullable=False)
    #: `end - start + 1`. Stored because the track draws a NEIGHBOUR's block at its own gene length.
    length_nt: Mapped[int] = mapped_column(Integer, nullable=False)

    strand: Mapped[str] = mapped_column(String(1), nullable=False)
    #: ⚠ Whether a strand parquet existed for this genome, or every gene was forced to `+`. A forced
    #: `+` must never read as an observation.
    strand_is_observed: Mapped[bool] = mapped_column(Boolean, nullable=False, default=True)

    #: The GFF phase column. Non-zero means the CDS does not begin at a start codon, and it is what
    #: `translate_cds` drops FIRST — before trimming to a codon boundary.
    gff_phase: Mapped[int] = mapped_column(SmallInteger, nullable=False, default=0)
    #: ⛔ An alternative initiator is promoted to `M` **only when the CDS is not 5'-partial**.
    #: Measured on `nhaA`: GTG is stored as `M`, and the plain table gives `V`.
    is_five_prime_partial: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    # ⛔ **`gc_percent` was DROPPED 2026-10-26** (David: *"I don't think we are using this anywhere
    # in BacAtlas right now. We don't need to hold it in the database ingest at all."*). It was
    # justified here by *"stored so `WHERE gc_percent > 60` never has to open a file"* — and that
    # query was never written. Measured across the repository: the only page that prints GC computes
    # it from the sliced bases (`GeneSequenceCard.vue`), the T6 parity suite compares that computed
    # value on both sides, no `locus` rollup holds GC, and nothing in `queries/`, `api/` or
    # `instruments/` referred to it. GC is a property of the bases, so it is derived where the bases
    # are. ⚠ This table now has no measurement column and therefore no `nan_guards`.
    protein_length_aa: Mapped[int | None] = mapped_column(Integer, nullable=True)

    bakta_gene_symbol: Mapped[str | None] = mapped_column(String(128), nullable=True)
    bakta_product: Mapped[str | None] = mapped_column(Text, nullable=True)
    #: The genome-private tag Bakta emits where it has no symbol (`AAOCBP_22210`).
    #: ⚠ Unique per genome, so it must never be treated as a gene NAME — a family holding only these
    #: would report one distinct symbol and read as unanimous.
    locus_tag: Mapped[str | None] = mapped_column(String(64), nullable=True)


class GeneFunctionalAnnotation(Base):
    """The six Bakta `Dbxref` label columns, one row per gene.

    ⚠ Ingesting this at all is a **collection-level policy**: at 80,000 genomes it is ~406 M rows
    for information whose only consumer is a per-locus rollup that ingest already computes.
    Measured coverage (ecoli/kp, job 34023134): UniRef50 98.5/99.4 %, COG 64.1/70.1 %,
    GO 71.9/**41.4** %, EC 37.3/27.9 %, KEGG 8.3/8.2 %.
    ⛔ A column absent from the parquet arrives NULL and is ANNOUNCED — files written before a
    column existed would otherwise report 0 % coverage for it, and *not measured* is not *zero*.
    """

    __tablename__ = "gene_functional_annotation"
    __table_args__ = (
        Index("ix_gene_functional_annotation__uniref50_accession", "uniref50_accession"),
        Index("ix_gene_functional_annotation__cog_id", "cog_id"),
    )

    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)
    flat_index: Mapped[int] = mapped_column(Integer, primary_key=True)

    uniref50_accession: Mapped[str | None] = mapped_column(String(64), nullable=True)
    kegg_orthology_id: Mapped[str | None] = mapped_column(String(32), nullable=True)
    #: Single-valued: `bakta_dbxrefs` takes the first `^COG:(COG\d+)$` hit and stops.
    cog_id: Mapped[str | None] = mapped_column(String(32), nullable=True)
    #: ⛔ **A SET of single-letter COG functional categories, not one category.** Bakta writes them
    #: concatenated (`CP`, `DZ`, `EFG`), and measured over 40 probe genomes **14,812 of 129,865
    #: non-null values (11.4 %) carry more than one letter**, up to 4. Stored as a scalar it fits —
    #: which is why the error would have been silent — but `WHERE cog_category = 'C'` then misses
    #: every gene whose categories are `CO`, `CP` or `CR`, and a per-category census under-counts by
    #: 11.4 %. This is the same set-in-a-string trap `ec` and `go` are modelled around, and it is
    #: modelled the same way here.
    #: ⚠ `locus.modal_cog_category` stays a STRING deliberately: it is the modal *value* the page
    #: prints via `strings.cog[cog_cat]`, and `CP` is what that value is.
    cog_categories: Mapped[list[str] | None] = mapped_column(ARRAY(String), nullable=True)
    #: ⛔ **A SET, not a scalar — and `String(32)` was too narrow to hold it.** Bakta's `ec` column
    #: is built exactly like `go`: `",".join(sorted(ec)) or None` (`bakta_dbxrefs.py:103`). Measured
    #: over 40 probe genomes: **8,660 of 70,477 non-null values (12.3 %) carry more than one term,
    #: up to 4 terms and 39 characters** — so the original scalar column would have raised on 12 %
    #: of the rows it was given, or silently truncated them under a driver that did not.
    #: Modelled as an array for the same reason `gene_ontology_terms` is: the reverse query ("which
    #: genes carry EC 2.7.7.7") is real, and a comma-joined string answers it with `LIKE '%…%'`.
    #: ⚠ NULL = not annotated; `{}` is never written.
    ec_numbers: Mapped[list[str] | None] = mapped_column(ARRAY(String), nullable=True)
    #: The comma-joined sorted-unique `go` column, exploded on `,` and nothing else. The `GO:`
    #: prefix is retained deliberately — it is part of the identifier, unlike the others' prefixes.
    #: Up to 26 terms / 285 characters measured. NULL = not annotated; `{}` is never written,
    #: because an empty list and an absent one are different claims.
    gene_ontology_terms: Mapped[list[str] | None] = mapped_column(ARRAY(String), nullable=True)


class GeneLocusMembership(Base):
    """Which locus each gene landed in, for one pangenome. **Exactly one row per gene.**

    The contract invariant (`every gene appears exactly once`) is enforced by the primary key here
    rather than by a `.duplicated()` pass at load time.
    """

    __tablename__ = "gene_locus_membership"
    __table_args__ = (
        Index("ix_gene_locus_membership__locus_id", "locus_id"),
        # ⛔⛔ **AN UNINDEXED FOREIGN KEY ON THIS TABLE IS NOT A TUNING QUESTION**, and the history is
        # kept here because the lesson outlives the column that taught it. A `locus_arrangement_id`
        # FK once sat on this table with no index, and deleting the pangenome's 70,519 arrangements
        # made Postgres run one sequential scan of 486,794 membership rows PER ARRANGEMENT to satisfy
        # `ON DELETE SET NULL`: the first load took 49 s into empty tables, the second sat in
        # `DELETE FROM locus_arrangement` for over 8 minutes and had to be killed. At 80,000 genomes
        # this table is ~412 M rows and 49 M arrangements.
        # ⭐ **That column and its FK are now GONE** (2026-10-26), so the stall is structurally
        # impossible rather than index-mitigated — and the index that mitigated it, which carried
        # 1,040,201 scans of pure FK maintenance and 45 MB, went with it. Any new FK added here needs
        # its own index in the same commit.
        # ⚠ `genome_id` FIRST, deliberately. The FK to `genome` is `ON DELETE CASCADE`, and a
        # composite index can only serve it from the LEADING column — `(pangenome_id, genome_id)`
        # could not, so deleting one genome scanned the whole table. Nothing is lost by the swap:
        # the primary key is `(pangenome_id, genome_id, flat_index)`, so `pangenome_id`-leading
        # lookups already have an index and the old one was redundant with it.
        Index("ix_gene_locus_membership__genome_id__pangenome_id", "genome_id", "pangenome_id"),
        # ⭐ The anchor query: which genes does THIS genome have at THIS locus, and in which
        # neighbourhood. Two index lookups instead of a linear scan over a per-genome int32 array.
        Index("ix_gene_locus_membership__locus_id__genome_id", "locus_id", "genome_id"),
    )

    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), primary_key=True
    )
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True)
    flat_index: Mapped[int] = mapped_column(Integer, primary_key=True)

    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)

    # ⛔ **`locus_arrangement_id` and `is_representative` were DROPPED 2026-10-26**, and neither is
    # coming back here. Measured before removing them: neither column was SELECTed by any query, any
    # serialiser or any test — `locus_arrangement_id` was written by `_load_gene_memberships` and
    # read by nothing, and `is_representative` existed only as this declaration. The docstring that
    # justified the first claimed it "keeps the `member_genome_ids` array honest", but nothing
    # enforced that: no test compared the two, so the honesty was an intention rather than a check.
    #
    # ⚠ The fact the first column carried — *"this genome has a gene here with no recorded
    # neighbourhood"*, which is a different sentence from *"this genome has no gene here"* — is
    # still available: `genes_without_a_window` is counted at ingest from the frame's own
    # `arrangement_rank`, not from this column, and `gene_arrangement_ranks` in
    # `genome_locus_assignment` carries it per gene once that table lands.
    # Design of record: `nuna/docs/design/storing_the_gene_layer.md`.


class GenomeNoncodingFeature(Base):
    """tRNA, rRNA, ncRNA, oriC, CRISPR — keyed `(genome, contig, index)`.

    ⛔⛔ **NOTHING READS THIS TABLE, AND IT IS NO LONGER WRITTEN BY DEFAULT** (David, 2026-10-29;
    traced the same day across both repos). There is no service, no route and no frontend field for
    any column here — only the ingest that writes it and the reconcile that counts it. The page's
    intergenic slot is drawn from `intergenic_gap`, whose length statistics come from GENE
    adjacencies and whose named features come from the `_noncoding` **parquet** via
    `catalogue_frames.load_bakta_noncoding`, never from this table. At 71,435 genomes it would be
    the single largest table in the database (~5.8 GB of ~19.7 GB), so `--with-noncoding-features`
    is now opt-in. ⚠ The table is kept in the schema only until D8's drop migration lands; do not
    build on it, and if the IGR work needs a per-genome feature store, design one for that work
    rather than reviving this.

    ⛔ **Deliberately carries no `flat_index`.** Admitting a non-coding row into that numbering would
    renumber every gene and shift every ±5 neighbourhood. It meets the gene table by coordinate
    containment only.

    **Cohort census, read from `data/proc/embeddings/meta/*_noncoding.parquet` (280 files,
    2026-09-04): 81,683 features, median 275/genome.** ncRNA 43.1 %, tRNA 26.2 %,
    regulatory_region 18.6 %, crispr-repeat 3.6 %, crispr-spacer 3.6 %, rRNA 1.9 %, oriC 1.6 %,
    oriT 0.7 %, CRISPR 0.3 %, tmRNA 0.3 %; **23.0 % overlap a CDS**.

    ⚠ Job 34060534's census reported **81,819** features and a median of 276 for the same 280
    genomes. Every proportion above matches it to a tenth of a point, so the shape is the same set —
    but the totals differ by 136, and the parquets are what the ingest actually loaded (verified:
    the load wrote 81,683 and the table holds 81,683). Both figures are kept because only one names
    the artifact it was read from, and that is the one of record.
    """

    __tablename__ = "genome_noncoding_feature"
    __table_args__ = (
        Index("ix_genome_noncoding_feature__genome_id__contig_index", "genome_id", "contig_index"),
        UniqueConstraint("genome_id", "contig_index", "noncoding_feature_index"),
    )

    genome_noncoding_feature_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.genome_id", ondelete="CASCADE"), nullable=False)
    contig_index: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    noncoding_feature_index: Mapped[int] = mapped_column(Integer, nullable=False)

    start_position: Mapped[int] = mapped_column(Integer, nullable=False)
    end_position: Mapped[int] = mapped_column(Integer, nullable=False)
    strand: Mapped[str | None] = mapped_column(String(1), nullable=True)
    feature_type: Mapped[str] = mapped_column(String(64), nullable=False)
    feature_gene: Mapped[str | None] = mapped_column(String(128), nullable=True)
    feature_product: Mapped[str | None] = mapped_column(Text, nullable=True)
    overlaps_coding_sequence: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
