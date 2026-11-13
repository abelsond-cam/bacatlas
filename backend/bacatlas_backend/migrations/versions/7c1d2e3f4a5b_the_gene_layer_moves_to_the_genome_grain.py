"""the gene layer moves to the genome grain

Two new tables and one new column, beside the three per-gene tables they replace — which stay for
one more migration so the two can be written together and compared row for row before anything is
dropped. Design of record: `nuna/docs/design/storing_the_gene_layer.md`.

* `genome_genes` — one row per genome: `gene_count`, `span_encoding`, `gene_spans` (16 B/gene,
  byte-identical to nuna's `.nseq` gene table).
* `genome_locus_assignment` — one row per (pangenome, genome): `gene_locus_ordinals` (4 B/gene,
  ≡ nuna's `.loci`) and `gene_arrangement_ranks`, with the three counts the loader checks.
* `locus.present_genome_ids int[] NOT NULL` + a GIN index — the one per-locus array kept on
  purpose, because *"which genomes carry locus L"* has no other answer once
  `gene_locus_membership` is gone.

⚠ `present_genome_ids` is added NOT NULL to a table that already holds rows (four catalogues on the
build database). It is created with a server default of the empty array, which every existing row
takes, and the default is then dropped so the models and the migration agree (the drift check
compares server defaults). An empty array on a pre-existing locus is honest: nothing has computed
its genomes yet, and the re-ingest that writes the new tables fills it. ⛔ Nothing reads the column
until the serving rewrite lands, so an empty array cannot be served as "no genome carries this".

Revision ID: 7c1d2e3f4a5b
Revises: 2ed24dc9ee53
Create Date: 2026-11-13 12:30:00.000000
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = "7c1d2e3f4a5b"
down_revision: str | None = "2ed24dc9ee53"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "genome_genes",
        sa.Column("genome_id", sa.Integer(), nullable=False),
        sa.Column("gene_count", sa.Integer(), nullable=False),
        sa.Column("span_encoding", sa.SmallInteger(), nullable=False),
        sa.Column("gene_spans", sa.LargeBinary(), nullable=False),
        sa.CheckConstraint("gene_count >= 0", name="gene_count_is_not_negative"),
        sa.CheckConstraint("octet_length(gene_spans) = gene_count * 16", name="gene_spans_is_16_bytes_per_gene"),
        sa.ForeignKeyConstraint(["genome_id"], ["genome.genome_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("genome_id"),
    )
    op.create_table(
        "genome_locus_assignment",
        sa.Column("pangenome_id", sa.BigInteger(), nullable=False),
        sa.Column("genome_id", sa.Integer(), nullable=False),
        sa.Column("gene_count", sa.Integer(), nullable=False),
        sa.Column("placed_gene_count", sa.Integer(), nullable=False),
        sa.Column("distinct_locus_count", sa.Integer(), nullable=False),
        sa.Column("gene_locus_ordinals", sa.LargeBinary(), nullable=False),
        sa.Column("gene_arrangement_ranks", sa.LargeBinary(), nullable=False),
        sa.CheckConstraint(
            "distinct_locus_count >= 0 AND distinct_locus_count <= placed_gene_count",
            name="distinct_loci_within_placed_genes",
        ),
        sa.CheckConstraint(
            "octet_length(gene_arrangement_ranks) = gene_count * 4", name="gene_arrangement_ranks_is_4_bytes_per_gene"
        ),
        sa.CheckConstraint(
            "octet_length(gene_locus_ordinals) = gene_count * 4", name="gene_locus_ordinals_is_4_bytes_per_gene"
        ),
        sa.CheckConstraint(
            "placed_gene_count >= 0 AND placed_gene_count <= gene_count", name="placed_genes_within_gene_count"
        ),
        sa.ForeignKeyConstraint(["genome_id"], ["genome.genome_id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["pangenome_id"], ["pangenome.pangenome_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("pangenome_id", "genome_id"),
    )
    op.add_column(
        "locus",
        sa.Column(
            "present_genome_ids",
            postgresql.ARRAY(sa.Integer()),
            nullable=False,
            server_default=sa.text("'{}'::integer[]"),
        ),
    )
    op.alter_column("locus", "present_genome_ids", server_default=None)
    op.create_index(
        "ix_locus__present_genome_ids", "locus", ["present_genome_ids"], unique=False, postgresql_using="gin"
    )


def downgrade() -> None:
    op.drop_index("ix_locus__present_genome_ids", table_name="locus", postgresql_using="gin")
    op.drop_column("locus", "present_genome_ids")
    op.drop_table("genome_locus_assignment")
    op.drop_table("genome_genes")
