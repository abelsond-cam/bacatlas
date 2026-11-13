"""a modelled genome can be named as holding no gene table

`genome.gene_table_omitted_reason` — nullable text. NULL means the genome has a `genome_genes` row;
non-NULL is WHY it has none (nuna's `.nseq` acceptance gate refused 17 genomes of the 71,435-genome
*K. pneumoniae* cohort because their stored proteins do not reproduce from the assembly). Such a
genome is still modelled — clustered, assigned and counted — and has no Sequence tab.

⛔ Read by hand: one nullable column, no default, no index, no constraint, so `downgrade` is exactly
the `drop_column` and recreates nothing else. No rows exist to backfill: NULL is the truth for every
genome loaded before this column existed, because every one of them was loaded with a gene table.

Revision ID: a0d837efeaad
Revises: 013b739e04ef
Create Date: 2026-11-13 20:26:16.555490
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision: str = "a0d837efeaad"
down_revision: str | None = "013b739e04ef"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("genome", sa.Column("gene_table_omitted_reason", sa.Text(), nullable=True))


def downgrade() -> None:
    op.drop_column("genome", "gene_table_omitted_reason")
