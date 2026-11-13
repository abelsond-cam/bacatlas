"""Each genome's two locus counts for one pangenome — written once, read on every picker keystroke.

⭐ **Counted from the rows just written, in the same transaction, not from the frames.** Both counts
are read back out of `genome_locus_assignment.distinct_locus_count` and
`locus_arrangement.member_genome_ids` — the rows the API serves — so the number the picker prints
cannot describe a different catalogue from the one the locus view draws, and a membership a loader
dropped cannot be counted anyway.

⚠ **The `locus_count` half is O(genomes) now, not O(genes).** `distinct_locus_count` is computed
when the genome's assignment block is built (`ingest_locus_catalogue._assignment_block`, from the
same memberships frame the old `count(DISTINCT locus_id)` aggregated), so this reads one integer per
genome rather than grouping ~412 M membership rows. The `arrangement_locus_count` half is unchanged
and still unnests `locus_arrangement.member_genome_ids` — minutes, at ingest — where the page's
`genomeCounts()` was the same pass on first opening of the box.
"""

from __future__ import annotations

from sqlalchemy import text
from sqlalchemy.orm import Session

#: ⭐ The definition of record. `pangenome_id` is the only parameter.
#:
#: ⛔ DISTINCT loci on BOTH sides, never a gene or row count: a genome at ρ > 1 has two genes at one
#: locus and can occupy two of its arrangements, so a row count is the genome's gene total — a
#: plausible number a little too large, and one nothing on the page could contradict.
#: `distinct_locus_count` is exactly `count(DISTINCT locus)` over the genome's placed genes, and
#: `test_genome_listing_endpoint` proves it equal to the old membership aggregate genome for genome.
#:
#: ⚠ Driven from the COLLECTION, with outer joins: a genome that contributed no gene still gets a row
#: with two zeros, because the picker lists every genome it can anchor and a missing row would drop
#: one silently rather than say it is in nothing.
#:
#: ⚠ The migration that created the table carries a frozen copy of the ORIGINAL (membership-based)
#: SELECT for its backfill — a migration must not import code that keeps changing, and that copy is
#: deliberately left as it was; `test_genome_listing_endpoint` recomputes this against the stored
#: rows so the definition and the table cannot drift apart unseen.
GENOME_LOCUS_COUNT_SELECT = """
    SELECT member.genome_id,
           coalesce(assigned.distinct_locus_count, 0) AS locus_count,
           coalesce(arranged.arrangement_locus_count, 0) AS arrangement_locus_count
      FROM pangenome p
      JOIN genome_collection_membership member
        ON member.genome_collection_id = p.genome_collection_id
      LEFT JOIN genome_locus_assignment assigned
        ON assigned.pangenome_id = :pangenome_id
       AND assigned.genome_id = member.genome_id
      LEFT JOIN (SELECT arranged_genome_id AS genome_id,
                        count(DISTINCT a.locus_id) AS arrangement_locus_count
                   FROM locus_arrangement a
                  CROSS JOIN LATERAL unnest(a.member_genome_ids) AS arranged_genome_id
                  WHERE a.pangenome_id = :pangenome_id
                  GROUP BY arranged_genome_id) arranged
        ON arranged.genome_id = member.genome_id
     WHERE p.pangenome_id = :pangenome_id
"""


def write_genome_locus_counts(session: Session, pangenome_id: int) -> int:
    """One row per collection genome for this pangenome, in ONE statement. Returns rows written.

    ⚠ Call it AFTER the arrangements and the genome assignments are loaded — it counts what is there.
    The caller deletes the pangenome's previous rows first (`_delete_pangenome_layer`), so a re-ingest
    replaces them wholesale like every other table in the layer.
    """
    result = session.execute(
        text(
            "INSERT INTO pangenome_genome_locus_count "
            "(pangenome_id, genome_id, locus_count, arrangement_locus_count) "
            "SELECT :pangenome_id, counted.genome_id, counted.locus_count, "
            "counted.arrangement_locus_count "
            f"FROM ({GENOME_LOCUS_COUNT_SELECT}) AS counted"
        ),
        {"pangenome_id": pangenome_id},
    )
    return result.rowcount
