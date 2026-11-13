"""How many loci each genome carries in one pangenome — the anchor dropdown's number, stored.

⭐ **Precomputed at ingest, never aggregated per request.** The published page derived this on first
opening of the anchor box (`app.js::genomeCounts`) by one pass over the ~490k-entry membership index it
already held. A server holds no such index in memory; the equivalent query was an aggregation over
the per-gene membership rows, **~412 M per pangenome at the 80,000-genome design target**, on every
keystroke in a picker. So it is written once, by the pangenome layer, in the same transaction as
the rows it counts.
"""

from __future__ import annotations

from sqlalchemy import BigInteger, CheckConstraint, ForeignKey, Index, Integer
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base


class PangenomeGenomeLocusCount(Base):
    """One genome of a pangenome's collection, and the two different numbers of loci it is in.

    ⛔ **Two counts, named apart, because they are different facts and only one is the page's.**

    * `locus_count` — distinct loci where the genome has **at least one gene** (`genome_locus_assignment`).
    * `arrangement_locus_count` — distinct loci where the genome sits in **some arrangement's**
      `member_genome_ids`. This is exactly `app.js::genomeCounts`, which is what the dropdown printed.

    They differ wherever a genome's gene reached no ±5 window — `coords` is an inner join in the export,
    so a gene alone on its contig is counted present at its locus and appears in no arrangement. That is
    the same gap that makes 6.26 % of ecoli loci "incomplete" (`locus.arrangement_member_genome_count`),
    seen from the genome side instead of the locus side, and it touches **every** probe genome.

    ⛔ **Both are `COUNT(DISTINCT locus)`, never `COUNT(*)`.** A genome at ρ > 1 has two genes at one
    locus and can sit in two arrangements there; counting rows instead of loci returns a plausible,
    slightly-too-large number — per-genome gene totals, not locus totals — that nothing on the page
    could contradict.

    ⚠ Every member of the pangenome's collection has a row, including one that contributed no gene at
    all (both counts 0). The dropdown lists every genome it can anchor; a missing row would drop a genome
    from the picker rather than say it is in nothing. `publish_pangenome` refuses a catalogue without one
    row per collection genome.
    """

    __tablename__ = "pangenome_genome_locus_count"
    __table_args__ = (
        # ⚠ The FK to `genome` is `ON DELETE CASCADE`, and the primary key leads with `pangenome_id`,
        # so without this a genome retraction scans the whole table — the omission that once turned a
        # routine re-ingest into an 8-minute stall on `gene_locus_membership`.
        Index("ix_pangenome_genome_locus_count__genome_id", "genome_id"),
        # ⛔ A genome cannot sit in an arrangement at a locus where it has no gene, so the arranged loci
        # are a subset of the present ones. A violation means the two counts were read from rows that
        # describe different catalogues, and the write fails rather than serving both.
        # ⚠ The name is short on purpose: with the `ck_<table>__` prefix it must fit Postgres' 63-byte
        # identifier limit, or the server truncates it and the models drift from the database.
        CheckConstraint(
            "arrangement_locus_count >= 0 AND arrangement_locus_count <= locus_count",
            name="arranged_within_present",
        ),
    )

    #: ⚠ BigInteger, because `pangenome.pangenome_id` is.
    pangenome_id: Mapped[int] = mapped_column(
        BigInteger, ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), primary_key=True
    )
    genome_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("genome.genome_id", ondelete="CASCADE"), primary_key=True
    )

    #: Distinct loci where this genome has ≥ 1 gene. `0` = the genome contributed no gene to this
    #: pangenome at all — a measured zero, never "not counted".
    locus_count: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Distinct loci where this genome is in some arrangement's `member_genome_ids` — the published
    #: page's `genomeCounts()`, which the anchor dropdown prints as "N loci".
    arrangement_locus_count: Mapped[int] = mapped_column(Integer, nullable=False)
