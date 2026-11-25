"""annotation-transfer predictions are served from their own table

Revision ID: 493f46d0cc95
Revises: a0d837efeaad
Create Date: 2026-11-14 22:35:31.700054

⭐ `served_locus_annotation` holds what annotation-transfer PREDICTS for a locus — one row per (locus,
claim), the class set, its calibrated accuracy, its donor and the two medoid cosines. On `pangenome`,
three columns say where the catalogue stands: `served_annotation_state` (not loaded · none produced ·
loaded), `annotation_model_version` and `served_annotation_sha256`, the last two set exactly when loaded.
⛔ Never `locus_annotation_entry`, which holds what a locus's own genes were observed to carry
(`nuna/docs/model_evaluation/annotation_transfer_gains.md` §8). No rows are written here: every catalogue
already loaded becomes NOT LOADED through the column default, which is the truth for each of them.

⛔ The donor key is DEFERRABLE INITIALLY DEFERRED. Postgres fires referential checks per deleted ROW, so a
plain one failed a catalogue reload whenever a donor locus went before its recipient (review of 9cfa481).

⚠ This revision was revised on its branch before it was applied anywhere outside a probe database; a
database at the earlier draft must `alembic downgrade a0d837efeaad` with that draft before upgrading.

⛔ **FOUR hand adjustments to what autogenerate wrote. Each is a failure it cannot see.**

1. **It omits every CHECK on `pangenome`** — `annotation_model_version_is_tagged`,
   `served_state_names_its_version`, `served_state_names_its_file`. Alembic does not compare CHECK
   constraints on an existing table, so `test_the_migration_still_matches_the_models` cannot see the gap
   either; `test_migration_round_trip` asserts all three exist after `upgrade head`. Their downgrade needs
   nothing: Postgres drops a CHECK with the columns it references.
2. **The downgrade did not drop `servedannotationkind`** (`postgres-schema-traps` trap 4), and
3. **nor `servedannotationstate`**: dropping the last column using an enum leaves the type, so down → up
   would die on `CREATE TYPE ... already exists`. Both are dropped by hand.
4. **The state column's server default is written out**, because existing rows need a value before NOT
   NULL can hold, and autogenerate does not compare server defaults in this repo (`compare_server_default`
   is off), so it would not notice the default missing either.
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = '493f46d0cc95'
down_revision: str | None = 'a0d837efeaad'
branch_labels = None
depends_on = None

SERVED_ANNOTATION_KIND = postgresql.ENUM(
    'COG_FUNCTIONAL_CATEGORY', 'GENE_ONTOLOGY_SLIM', name='servedannotationkind', create_type=False
)
SERVED_ANNOTATION_STATE = postgresql.ENUM(
    'NOT_LOADED', 'NONE_PRODUCED', 'LOADED', name='servedannotationstate', create_type=False
)


def upgrade() -> None:
    SERVED_ANNOTATION_KIND.create(op.get_bind(), checkfirst=False)
    SERVED_ANNOTATION_STATE.create(op.get_bind(), checkfirst=False)
    op.create_table('served_locus_annotation',
    sa.Column('served_locus_annotation_id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('locus_id', sa.BigInteger(), nullable=False),
    sa.Column('annotation_kind', SERVED_ANNOTATION_KIND, nullable=False),
    sa.Column('gene_ontology_namespace', sa.SmallInteger(), nullable=True),
    sa.Column('term_values', postgresql.ARRAY(sa.Text()), nullable=False),
    sa.Column('estimated_accuracy', sa.Float(), nullable=False),
    sa.Column('donor_locus_id', sa.BigInteger(), nullable=False),
    sa.Column('esm_medoid_similarity', sa.Float(), nullable=False),
    sa.Column('bacformer_medoid_similarity', sa.Float(), nullable=False),
    sa.Column('annotation_model_version', sa.String(length=32), nullable=False),
    sa.CheckConstraint("(annotation_kind = 'GENE_ONTOLOGY_SLIM') = (gene_ontology_namespace IS NOT NULL)", name=op.f('ck_served_locus_annotation__namespace_matches_kind')),
    sa.CheckConstraint("annotation_model_version ~ '^at-v[0-9]+\\.[0-9]+$'", name=op.f('ck_served_locus_annotation__version_is_tagged')),
    sa.CheckConstraint('cardinality(term_values) >= 1 AND array_position(term_values, NULL) IS NULL', name=op.f('ck_served_locus_annotation__terms_are_stated')),
    sa.CheckConstraint('esm_medoid_similarity BETWEEN -1.001 AND 1.001 AND bacformer_medoid_similarity BETWEEN -1.001 AND 1.001', name=op.f('ck_served_locus_annotation__similarities_are_cosines')),
    sa.CheckConstraint('estimated_accuracy >= 0 AND estimated_accuracy <= 1', name=op.f('ck_served_locus_annotation__accuracy_is_a_probability')),
    sa.CheckConstraint('gene_ontology_namespace IS NULL OR gene_ontology_namespace BETWEEN 0 AND 2', name=op.f('ck_served_locus_annotation__namespace_is_known')),
    sa.ForeignKeyConstraint(['donor_locus_id'], ['locus.locus_id'], name=op.f('fk_served_locus_annotation__donor_locus_id__locus'), initially='DEFERRED', deferrable=True),
    sa.ForeignKeyConstraint(['locus_id'], ['locus.locus_id'], name=op.f('fk_served_locus_annotation__locus_id__locus'), ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('served_locus_annotation_id', name=op.f('pk_served_locus_annotation'))
    )
    op.create_index('ix_served_locus_annotation__donor_locus_id', 'served_locus_annotation', ['donor_locus_id'], unique=False)
    op.create_index('uq_served_locus_annotation__locus_kind_namespace', 'served_locus_annotation', ['locus_id', 'annotation_kind', sa.literal_column('COALESCE(gene_ontology_namespace, -1)')], unique=True)
    # ⛔ adjustment 4 — every existing catalogue is NOT LOADED, which is true of each.
    op.add_column('pangenome', sa.Column(
        'served_annotation_state', SERVED_ANNOTATION_STATE, server_default='NOT_LOADED', nullable=False
    ))
    op.add_column('pangenome', sa.Column('annotation_model_version', sa.String(length=32), nullable=True))
    op.add_column('pangenome', sa.Column('served_annotation_sha256', sa.String(length=64), nullable=True))
    # ⛔ adjustment 1 — see the module docstring.
    op.create_check_constraint(
        op.f('ck_pangenome__annotation_model_version_is_tagged'),
        'pangenome',
        "annotation_model_version IS NULL OR annotation_model_version ~ '^at-v[0-9]+\\.[0-9]+$'",
    )
    op.create_check_constraint(
        op.f('ck_pangenome__served_state_names_its_version'),
        'pangenome',
        "(served_annotation_state = 'LOADED') = (annotation_model_version IS NOT NULL)",
    )
    op.create_check_constraint(
        op.f('ck_pangenome__served_state_names_its_file'),
        'pangenome',
        "(served_annotation_state = 'LOADED') = (served_annotation_sha256 IS NOT NULL)",
    )


def downgrade() -> None:
    # The CHECKs go with their columns (adjustment 1).
    op.drop_column('pangenome', 'served_annotation_sha256')
    op.drop_column('pangenome', 'annotation_model_version')
    op.drop_column('pangenome', 'served_annotation_state')
    op.drop_index('uq_served_locus_annotation__locus_kind_namespace', table_name='served_locus_annotation')
    op.drop_index('ix_served_locus_annotation__donor_locus_id', table_name='served_locus_annotation')
    op.drop_table('served_locus_annotation')
    # ⛔ adjustments 2 and 3 — dropping a table or a column does NOT drop its enum type.
    SERVED_ANNOTATION_KIND.drop(op.get_bind(), checkfirst=False)
    SERVED_ANNOTATION_STATE.drop(op.get_bind(), checkfirst=False)
