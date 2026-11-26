"""the locus card states the audit's numbers, not its verdicts

Revision ID: 9601b9c2d400
Revises: 493f46d0cc95
Create Date: 2026-11-26 19:08:36.542675

⭐ **David's 2026-11-26 ruling.** At 71,435 genomes the audit's per-locus verdicts were decided by stray
genes — kp-71k locus 4008 (Wzi, PF14052 on 69,825 of 69,884 genes) was shown as "grouped on context alone"
and "conflicting architecture". The redesigned audit states numbers instead, and `locus` gains seven
columns to carry them: the dominant Pfam domain SET and its gene count, the count of annotated genes on any
other set, the audit's own `pfam_shared` (that set on more than 95 % of ALL genes), the share of sampled
members MMseqs groups at the resolved identity, and the 5th percentile and median of member → medoid ESM
cosine.

⚠ **All seven are NULL on every catalogue already loaded, and that is the truth for each of them**: those
audits predate 2026-11-26 and measured none of it. No row is written here. The ingest refuses an audit that
lacks the new columns rather than filling them.

⛔ **ONE hand adjustment to what autogenerate wrote, and it is a failure autogenerate cannot see**: it
omits every CHECK on an existing table (Alembic compares none), so the four below are added by hand and
`test_migration_round_trip` asserts each exists after `upgrade head`. Their downgrade needs nothing beyond
the columns: Postgres drops a CHECK with the columns it references, and no enum type is added.
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision: str = '9601b9c2d400'
down_revision: str | None = '493f46d0cc95'
branch_labels = None
depends_on = None

#: `models/served_locus_annotation.COSINE_BOUND`, written out: a migration must not import a model
#: module, whose value may move after this revision is applied.
COSINE_BOUND = 1.001


def upgrade() -> None:
    op.add_column('locus', sa.Column('pfam_dominant_domains', sa.Text(), nullable=True))
    op.add_column('locus', sa.Column('pfam_dominant_domain_gene_count', sa.Integer(), nullable=True))
    op.add_column('locus', sa.Column('pfam_other_architecture_gene_count', sa.Integer(), nullable=True))
    op.add_column('locus', sa.Column('pfam_architecture_shared', sa.Boolean(), nullable=True))
    op.add_column('locus', sa.Column('mmseqs_group_share', sa.Float(), nullable=True))
    op.add_column('locus', sa.Column('esm_member_similarity_p05', sa.Float(), nullable=True))
    op.add_column('locus', sa.Column('esm_member_similarity_median', sa.Float(), nullable=True))
    # ⛔ The hand adjustment — see the module docstring.
    op.create_check_constraint(
        op.f('ck_locus__mmseqs_group_share_is_a_share'),
        'locus',
        "mmseqs_group_share IS NULL OR mmseqs_group_share BETWEEN 0 AND 1",
    )
    op.create_check_constraint(
        op.f('ck_locus__esm_member_similarities_are_cosines'),
        'locus',
        f"(esm_member_similarity_p05 IS NULL OR esm_member_similarity_p05 BETWEEN -{COSINE_BOUND} AND "
        f"{COSINE_BOUND}) AND (esm_member_similarity_median IS NULL OR esm_member_similarity_median "
        f"BETWEEN -{COSINE_BOUND} AND {COSINE_BOUND})",
    )
    op.create_check_constraint(
        op.f('ck_locus__pfam_architecture_counts_are_counts'),
        'locus',
        "(pfam_dominant_domain_gene_count IS NULL OR pfam_dominant_domain_gene_count >= 0) AND "
        "(pfam_other_architecture_gene_count IS NULL OR pfam_other_architecture_gene_count >= 0)",
    )
    op.create_check_constraint(
        op.f('ck_locus__pfam_architectures_within_annotated'),
        'locus',
        "pfam_dominant_domain_gene_count IS NULL OR pfam_other_architecture_gene_count IS NULL "
        "OR pfam_annotated_member_count IS NULL OR pfam_dominant_domain_gene_count + "
        "pfam_other_architecture_gene_count <= pfam_annotated_member_count",
    )


def downgrade() -> None:
    # The four CHECKs go with their columns (see the module docstring).
    op.drop_column('locus', 'esm_member_similarity_median')
    op.drop_column('locus', 'esm_member_similarity_p05')
    op.drop_column('locus', 'mmseqs_group_share')
    op.drop_column('locus', 'pfam_architecture_shared')
    op.drop_column('locus', 'pfam_other_architecture_gene_count')
    op.drop_column('locus', 'pfam_dominant_domain_gene_count')
    op.drop_column('locus', 'pfam_dominant_domains')
