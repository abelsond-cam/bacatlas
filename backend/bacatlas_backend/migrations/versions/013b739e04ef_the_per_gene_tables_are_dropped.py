"""the per-gene tables are dropped

D8. `gene`, `gene_functional_annotation`, `gene_locus_membership` and `genome_noncoding_feature` go —
372 M rows each at 71,435 genomes. Everything the page reads now comes from the genome grain
(`7c1d2e3f4a5b`): `genome_genes` (packed spans), `genome_locus_assignment` (packed locus ordinals and
arrangement ranks) and `locus.present_genome_ids`. The one inbound foreign key, from
`projected_gene_placement` onto `gene`, was dropped first (`9e8f7a6b5c4d`).

⛔ `downgrade` recreates each table EXACTLY as the earlier migrations left it — primary key, foreign
keys, CHECKs, the unique constraint and every index, including `ix_gene__genome_id`'s BRIN — and
recreates no rows. Checked by hand against a `pg_dump --schema-only` taken at `9e8f7a6b5c4d`.

Revision ID: 013b739e04ef
Revises: 9e8f7a6b5c4d
Create Date: 2026-11-13 12:36:05.665295
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = "013b739e04ef"
down_revision: str | None = "9e8f7a6b5c4d"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.drop_index(op.f('ix_genome_noncoding_feature__genome_id__contig_index'), table_name='genome_noncoding_feature')
    op.drop_table('genome_noncoding_feature')
    op.drop_index(op.f('ix_gene_locus_membership__genome_id__pangenome_id'), table_name='gene_locus_membership')
    op.drop_index(op.f('ix_gene_locus_membership__locus_id'), table_name='gene_locus_membership')
    op.drop_index(op.f('ix_gene_locus_membership__locus_id__genome_id'), table_name='gene_locus_membership')
    op.drop_table('gene_locus_membership')
    op.drop_index(op.f('ix_gene__genome_id'), table_name='gene', postgresql_using='brin')
    op.drop_index(op.f('ix_gene__genome_id__contig_index__start_position'), table_name='gene')
    op.drop_table('gene')
    op.drop_index(op.f('ix_gene_functional_annotation__cog_id'), table_name='gene_functional_annotation')
    op.drop_index(op.f('ix_gene_functional_annotation__uniref50_accession'), table_name='gene_functional_annotation')
    op.drop_table('gene_functional_annotation')


def downgrade() -> None:
    op.create_table('gene_functional_annotation',
    sa.Column('genome_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('flat_index', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('uniref50_accession', sa.VARCHAR(length=64), autoincrement=False, nullable=True),
    sa.Column('kegg_orthology_id', sa.VARCHAR(length=32), autoincrement=False, nullable=True),
    sa.Column('cog_id', sa.VARCHAR(length=32), autoincrement=False, nullable=True),
    sa.Column('cog_categories', postgresql.ARRAY(sa.VARCHAR()), autoincrement=False, nullable=True),
    sa.Column('ec_numbers', postgresql.ARRAY(sa.VARCHAR()), autoincrement=False, nullable=True),
    sa.Column('gene_ontology_terms', postgresql.ARRAY(sa.VARCHAR()), autoincrement=False, nullable=True),
    sa.ForeignKeyConstraint(['genome_id'], ['genome.genome_id'], name=op.f('fk_gene_functional_annotation__genome_id__genome'), ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('genome_id', 'flat_index', name=op.f('pk_gene_functional_annotation'))
    )
    op.create_index(op.f('ix_gene_functional_annotation__uniref50_accession'), 'gene_functional_annotation', ['uniref50_accession'], unique=False)
    op.create_index(op.f('ix_gene_functional_annotation__cog_id'), 'gene_functional_annotation', ['cog_id'], unique=False)
    op.create_table('gene',
    sa.Column('genome_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('flat_index', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('pathogen_species_id', sa.SMALLINT(), autoincrement=False, nullable=False),
    sa.Column('contig_index', sa.SMALLINT(), autoincrement=False, nullable=False),
    sa.Column('start_position', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('end_position', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('length_nt', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('strand', sa.VARCHAR(length=1), autoincrement=False, nullable=False),
    sa.Column('strand_is_observed', sa.BOOLEAN(), autoincrement=False, nullable=False),
    sa.Column('gff_phase', sa.SMALLINT(), autoincrement=False, nullable=False),
    sa.Column('is_five_prime_partial', sa.BOOLEAN(), autoincrement=False, nullable=False),
    sa.Column('protein_length_aa', sa.INTEGER(), autoincrement=False, nullable=True),
    sa.Column('bakta_gene_symbol', sa.VARCHAR(length=128), autoincrement=False, nullable=True),
    sa.Column('bakta_product', sa.TEXT(), autoincrement=False, nullable=True),
    sa.Column('locus_tag', sa.VARCHAR(length=64), autoincrement=False, nullable=True),
    # ⚠ The initial migration's own spelling, so the stored definition is byte-identical.
    sa.CheckConstraint("strand IN ('+', '-')", name=op.f('ck_gene__strand_is_plus_or_minus')),
    sa.CheckConstraint('start_position >= 1 AND end_position >= start_position', name=op.f('ck_gene__span_is_forward')),
    sa.ForeignKeyConstraint(['genome_id'], ['genome.genome_id'], name=op.f('fk_gene__genome_id__genome'), ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('genome_id', 'flat_index', name=op.f('pk_gene'))
    )
    op.create_index(op.f('ix_gene__genome_id__contig_index__start_position'), 'gene', ['genome_id', 'contig_index', 'start_position'], unique=False)
    op.create_index(op.f('ix_gene__genome_id'), 'gene', ['genome_id'], unique=False, postgresql_using='brin')
    op.create_table('gene_locus_membership',
    sa.Column('pangenome_id', sa.BIGINT(), autoincrement=False, nullable=False),
    sa.Column('genome_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('flat_index', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('locus_id', sa.BIGINT(), autoincrement=False, nullable=False),
    sa.ForeignKeyConstraint(['genome_id'], ['genome.genome_id'], name=op.f('fk_gene_locus_membership__genome_id__genome'), ondelete='CASCADE'),
    sa.ForeignKeyConstraint(['locus_id'], ['locus.locus_id'], name=op.f('fk_gene_locus_membership__locus_id__locus'), ondelete='CASCADE'),
    sa.ForeignKeyConstraint(['pangenome_id'], ['pangenome.pangenome_id'], name=op.f('fk_gene_locus_membership__pangenome_id__pangenome'), ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('pangenome_id', 'genome_id', 'flat_index', name=op.f('pk_gene_locus_membership'))
    )
    op.create_index(op.f('ix_gene_locus_membership__locus_id__genome_id'), 'gene_locus_membership', ['locus_id', 'genome_id'], unique=False)
    op.create_index(op.f('ix_gene_locus_membership__locus_id'), 'gene_locus_membership', ['locus_id'], unique=False)
    op.create_index(op.f('ix_gene_locus_membership__genome_id__pangenome_id'), 'gene_locus_membership', ['genome_id', 'pangenome_id'], unique=False)
    op.create_table('genome_noncoding_feature',
    sa.Column('genome_noncoding_feature_id', sa.BIGINT(), autoincrement=True, nullable=False),
    sa.Column('genome_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('contig_index', sa.SMALLINT(), autoincrement=False, nullable=False),
    sa.Column('noncoding_feature_index', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('start_position', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('end_position', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('strand', sa.VARCHAR(length=1), autoincrement=False, nullable=True),
    sa.Column('feature_type', sa.VARCHAR(length=64), autoincrement=False, nullable=False),
    sa.Column('feature_gene', sa.VARCHAR(length=128), autoincrement=False, nullable=True),
    sa.Column('feature_product', sa.TEXT(), autoincrement=False, nullable=True),
    sa.Column('overlaps_coding_sequence', sa.BOOLEAN(), autoincrement=False, nullable=True),
    sa.ForeignKeyConstraint(['genome_id'], ['genome.genome_id'], name=op.f('fk_genome_noncoding_feature__genome_id__genome'), ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('genome_noncoding_feature_id', name=op.f('pk_genome_noncoding_feature')),
    sa.UniqueConstraint('genome_id', 'contig_index', 'noncoding_feature_index', name=op.f('uq_genome_noncoding_feature__genome_id_contig_index_non_c8d5'))
    )
    op.create_index(op.f('ix_genome_noncoding_feature__genome_id__contig_index'), 'genome_noncoding_feature', ['genome_id', 'contig_index'], unique=False)
