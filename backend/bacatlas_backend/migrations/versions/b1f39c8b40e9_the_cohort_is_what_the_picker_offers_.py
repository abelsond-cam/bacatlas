"""the cohort is what the picker offers and the model is one label

Revision ID: b1f39c8b40e9
Revises: c99e925f7341
Create Date: 2026-10-27 20:34:40.366266

Two columns, both NOT NULL, both on tables that already hold rows. ⛔ **Autogenerate's form of this
is a migration that cannot run** — `ADD COLUMN … NOT NULL` with no default is an error on a
non-empty table, and `bacatlas_ingest` holds 2 collections and 4 pangenomes. So each column is added
nullable, backfilled from data **already in this database**, and only then made NOT NULL. On an
empty database every UPDATE is a no-op and the SET NOT NULL still succeeds, which is what
`test_migration_round_trip` exercises.

⭐ **Neither backfill invents a value.**

- `pangenome.model_label` comes from the `pangenome_evaluation` row with `metric_name = 'label'` —
  nuna's own `export_payload --model-label`, and already the source `meta.model_label` is served
  from (`payload_serialiser.py:819-822`). Measured here before writing: all 4 catalogues have one.
  The COALESCE chain behind it (`nuna_model.label`, then `run_id`) exists so the SET NOT NULL cannot
  fail on a database where the evaluation row is missing — never because a guess is acceptable.
- `genome_collection.collection_label` is backfilled **mechanically from `collection_key`**, which
  is deliberately not a display string. ⚠ The real labels ("E. coli probe", the Kp cohorts' names)
  are user-facing text and are supplied by the ingest's `--collection-label`; a schema migration is
  the wrong place to choose words a reader will see. Until a catalogue is re-ingested, the picker
  would show a key — which is the visible, correctable state, not a silent one.
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision: str = "b1f39c8b40e9"
down_revision: str | None = "c99e925f7341"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("genome_collection", sa.Column("collection_label", sa.String(length=128), nullable=True))
    op.execute("UPDATE genome_collection SET collection_label = collection_key WHERE collection_label IS NULL")
    op.alter_column("genome_collection", "collection_label", nullable=False)

    op.add_column("pangenome", sa.Column("model_label", sa.String(length=256), nullable=True))
    op.execute(
        """
        UPDATE pangenome AS p
           SET model_label = COALESCE(
                 (SELECT e.detail
                    FROM pangenome_evaluation AS e
                   WHERE e.pangenome_id = p.pangenome_id
                     AND e.metric_name = 'label'
                   LIMIT 1),
                 (SELECT m.label FROM nuna_model AS m WHERE m.nuna_model_id = p.nuna_model_id),
                 p.run_id)
         WHERE p.model_label IS NULL
        """
    )
    op.alter_column("pangenome", "model_label", nullable=False)


def downgrade() -> None:
    op.drop_column("pangenome", "model_label")
    op.drop_column("genome_collection", "collection_label")
