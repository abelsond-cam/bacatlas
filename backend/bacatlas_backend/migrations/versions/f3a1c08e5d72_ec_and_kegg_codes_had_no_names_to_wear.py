"""ec and kegg codes had no names to wear

`term_name` is populated for 100 % of the catalogue's COG (20,823) and GO (59,879) entries and
**0 %** of its EC (9,398) and KEGG (5,309) ones, because the names were never vendored. So the
Function tab showed `2.7.10.-` under a seven-word gloss — *"transferase · 4th level not stated"* —
and `K01991` as a bare link.

⭐ The EC gloss was a workaround for a licence restriction that does not exist. ExPASy ENZYME is
**CC BY 4.0**, read from `enzuser.txt`, not CC BY-ND. `2.7.10.-` is *Protein-tyrosine kinases*.

⛔ KEGG's licence is NOT resolved and this migration does not resolve it: KEGG requires an academic
service provider licence from anyone offering services. The table exists on David's decision of
2026-10-05 for an internal prototype, with the question carried on the page. Nothing ships to
BacAtlas.org until the subscription is confirmed.

⚠ Global, not per pangenome — an EC code is the same reaction in every species and model. Both
tables are populated by `ingest --stage reference`, which is idempotent, so this needs no catalogue
re-ingest.

Revision ID: f3a1c08e5d72
Revises: dd9e1b5d60b2
Create Date: 2026-10-05 23:50:00.000000
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision: str = 'f3a1c08e5d72'
down_revision: str | None = 'dd9e1b5d60b2'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        'enzyme_class',
        # ⛔ Dash-padded exactly as the annotation carries it (`2.7.10.-`), never normalised to
        # `2.7.10`: 126 of the 1,410 codes in the two catalogues state fewer than four fields, and
        # the dashes are how the depth is stated. Width 32 leaves room for the longest sub-subclass.
        sa.Column('ec_code', sa.String(length=32), nullable=False),
        # ⚠ NOT NULL with an empty default, as `pfam_family` does: "" means the reference had no
        # value, and a NULL here would invite exactly the null-vs-empty confusion the rest of this
        # schema spends its effort avoiding.
        sa.Column('name', sa.Text(), nullable=False, server_default=''),
        sa.PrimaryKeyConstraint('ec_code', name=op.f('pk_enzyme_class')),
    )
    op.create_table(
        'kegg_orthology',
        sa.Column('ko_id', sa.String(length=16), nullable=False),
        # ⚠ Plural and optional: `K01991` is `wza, gfcE`, and 1,186 of 28,516 KOs have no symbol.
        sa.Column('symbol', sa.String(length=256), nullable=False, server_default=''),
        sa.Column('definition', sa.Text(), nullable=False, server_default=''),
        sa.PrimaryKeyConstraint('ko_id', name=op.f('pk_kegg_orthology')),
    )


def downgrade() -> None:
    op.drop_table('kegg_orthology')
    op.drop_table('enzyme_class')
