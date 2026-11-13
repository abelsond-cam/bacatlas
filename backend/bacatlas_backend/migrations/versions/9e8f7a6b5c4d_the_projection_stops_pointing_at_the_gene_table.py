"""the projection stops pointing at the gene table

`projected_gene_placement` loses its composite foreign key `(genome_id, flat_index) →
gene (genome_id, flat_index)`. The gene layer is one row per genome now (`genome_genes`, migration
`7c1d2e3f4a5b`), so once `gene` is dropped there is no row for the constraint to reference — and it
was the ONLY inbound foreign key onto the per-gene table (`storing_the_gene_layer.md` §6 item 3).

The guarantee moves into the loader: `ingest_projected_genomes._check_the_genes_are_this_genomes`
reads `genome_genes.gene_count` and refuses any placed `flat_index` outside `[0, gene_count)`. A
check in a loader is weaker than a constraint in a schema, and that is a chosen cost, not an
oversight.

⚠ Nothing else about the table changes: the primary key, the FK onto `projected_genome` (which is
what cascades a re-ingest), both indexes and every CHECK stay as `b41f7c9e2a08` created them.
`downgrade` recreates the constraint under its original name so the round-trip is exact — which
requires every placement to reference an existing `gene` row at that moment, as it did before.

Revision ID: 9e8f7a6b5c4d
Revises: 7c1d2e3f4a5b
Create Date: 2026-11-13 18:00:00.000000
"""

from __future__ import annotations

from alembic import op

revision: str = "9e8f7a6b5c4d"
down_revision: str | None = "7c1d2e3f4a5b"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.drop_constraint("fk_projected_gene_placement__gene", "projected_gene_placement", type_="foreignkey")


def downgrade() -> None:
    op.create_foreign_key(
        "fk_projected_gene_placement__gene",
        "projected_gene_placement",
        "gene",
        ["genome_id", "flat_index"],
        ["genome_id", "flat_index"],
        ondelete="CASCADE",
    )
