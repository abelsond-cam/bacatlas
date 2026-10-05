"""the slim is a selection not a hierarchy

⭐⭐ The table that fixes the GO comparison. `goslim_metagenomics` holds `membrane`, `plasma
membrane` and `outer membrane` as three separate classes, and all three namespace roots besides — so
folding a term onto it put that term and its own parent in different buckets, and the set algebra
that followed called them `disjoint`. Klebsiella `gumC` read "classes differ" for 16 genes saying
plasma membrane and one saying membrane. 11 of the 43 `disjoint` verdicts across both catalogues
involved a namespace root, which is an ancestor of everything.

`gene_ontology_slim_ancestor` carries the slim's own hierarchy, TRANSITIVELY and with no root. 111
classes, 108 of them under another. `gene_ontology_term` carries the raw term → (namespace, slim
classes) map the within-node propagation side needs, because it reads
`gene_functional_annotation.gene_ontology_terms` — raw accessions, all three namespaces pooled.

⛔ No row names `GO:0003674`, `GO:0008150` or `GO:0005575`. Under full ancestor closure every pair
shares its namespace root and `disjoint` becomes unreachable — `vendor_go`'s own standing objection,
and it is right. Restricting the closure to the slim and dropping the roots is what keeps two
genuinely different branches disjoint.

⚠ Global, not per pangenome, and no catalogue is re-ingested: `ingest --stage reference` is the
whole data migration, and a future change to the measure needs no re-ingest either.

Revision ID: b7d4e21a9c56
Revises: f3a1c08e5d72
Create Date: 2026-10-06 00:10:00.000000
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


revision: str = 'b7d4e21a9c56'
down_revision: str | None = 'f3a1c08e5d72'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        'gene_ontology_term',
        sa.Column('go_id', sa.String(length=16), nullable=False),
        sa.Column('name', sa.Text(), nullable=False, server_default=''),
        # ⛔ An INDEX, not a namespace string: `enumerations.GENE_ONTOLOGY_NAMESPACE_NAMES` already
        # owns the index↔name contract and this joins straight to
        # `locus_annotation_entry.gene_ontology_namespace`. A second spelling is how two parts of one
        # response come to disagree about which namespace a term is in.
        sa.Column('namespace_index', sa.SmallInteger(), nullable=False),
        # ⚠ EVERY nearest slim class, not one — dropping the others would make two genes that share a
        # dropped ancestor read as disjoint. May be empty: such a term self-maps on read.
        sa.Column('slim_go_ids', postgresql.ARRAY(sa.String(length=16)), nullable=False),
        sa.PrimaryKeyConstraint('go_id', name=op.f('pk_gene_ontology_term')),
    )
    op.create_index(
        op.f('ix_gene_ontology_term__namespace_index'),
        'gene_ontology_term', ['namespace_index'], unique=False,
    )
    op.create_table(
        'gene_ontology_slim_ancestor',
        sa.Column('go_id', sa.String(length=16), nullable=False),
        sa.Column('ancestor_go_id', sa.String(length=16), nullable=False),
        # ⚠ The FKs are the point of an edge table over an array column: "an ancestor that is not a
        # known class" becomes impossible rather than merely unlikely.
        sa.ForeignKeyConstraint(
            ['go_id'], ['gene_ontology_term.go_id'],
            name=op.f('fk_gene_ontology_slim_ancestor__go_id__gene_ontology_term'),
        ),
        sa.ForeignKeyConstraint(
            ['ancestor_go_id'], ['gene_ontology_term.go_id'],
            name=op.f('fk_gene_ontology_slim_ancestor__ancestor_go_id__gene_ontology_term'),
        ),
        sa.PrimaryKeyConstraint(
            'go_id', 'ancestor_go_id', name=op.f('pk_gene_ontology_slim_ancestor')
        ),
    )
    # The reverse direction — which classes are under `membrane`? — is a real query.
    op.create_index(
        op.f('ix_gene_ontology_slim_ancestor__ancestor'),
        'gene_ontology_slim_ancestor', ['ancestor_go_id'], unique=False,
    )


def downgrade() -> None:
    op.drop_index(
        op.f('ix_gene_ontology_slim_ancestor__ancestor'), table_name='gene_ontology_slim_ancestor'
    )
    op.drop_table('gene_ontology_slim_ancestor')
    op.drop_index(op.f('ix_gene_ontology_term__namespace_index'), table_name='gene_ontology_term')
    op.drop_table('gene_ontology_term')
