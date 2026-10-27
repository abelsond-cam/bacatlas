"""the two probe cohorts get the names a reader sees

⭐ **`b1f39c8b40e9` deliberately backfilled `collection_label` to the collection KEY**, because a key
in a picker reads as unfinished, where a blank or a guessed species name would read as finished and
wrong. The placeholder has done its job: David named the cohorts on 2026-10-27, so this replaces it.

    ecoli / ecoli  ->  'E. coli AMR'
    kp    / kp     ->  'Kp AMR'

⛔ **THE LABEL DOES NOT CARRY THE COUNT, and the words David gave did.** He wrote "E.coli AMR 100,
Kp AMR 100, Kp 71,500". The picker renders `"{label} — {count:,} genomes"` from this column and
`genome_collection.genome_count`, which is *passed, never inferred* — so storing "Kp 71,500" would
render "Kp 71,500 — 71,435 genomes": a label contradicting the roster it names, on the first day,
with nothing anywhere able to say which half was right. The words are stored here, the number comes
from the roster, and the rendered option reads as he asked for it.

⛔ **A DATA migration, and load-bearing rather than cosmetic.** `verify_pangenome_is_servable` and
`offer_catalogue` both refuse a cohort whose label is still its own key, so until this runs the live
`ecoli-nuna4` and `kp-nuna4` FAIL the publish gate. Measured on `bacatlas_ingest`: 14 checks, 13
passed, 1 failed before; 14 passed after. A guard that the production data violates is not a guard.

⚠ **It replaces only the PLACEHOLDER** (`WHERE collection_label = collection_key`). A label somebody
has already chosen is never overwritten, so re-running this — or running it after a better name has
landed — changes nothing.

⚠ **And the label is not yet STICKY.** `ingest_genome_collection` assigns
`collection_label = collection_label or key` unconditionally, so the next `--stage pangenome` run
that omits `--collection-label` puts the key back; `offer_catalogue`'s guard fires only on the way
in, so a catalogue already published would not re-check and the machine token would reach a reader
with nothing refusing it. Making the checked-in cohort declaration the ingest's own default is what
closes that, and it belongs beside the declaration rather than here.

Revision ID: 166560ca71e2
Revises: b1f39c8b40e9
Create Date: 2026-10-27 23:03:20.594924
"""

from __future__ import annotations

from alembic import op

revision: str = '166560ca71e2'
down_revision: str | None = 'b1f39c8b40e9'
branch_labels = None
depends_on = None


#: `(species_key, collection_key, label)` — David, 2026-10-27.
#: ⚠ Keyed on the SPECIES KEY, never on `genome_collection_id`: the ids are whatever this particular
#: database happened to assign, and a migration that names them is correct on one machine.
COHORT_LABELS = (
    ("ecoli", "ecoli", "E. coli AMR"),
    ("kp", "kp", "Kp AMR"),
)

_RELABEL = """
    UPDATE genome_collection AS gc
       SET collection_label = :new_label
      FROM pathogen_species AS ps
     WHERE ps.pathogen_species_id = gc.pathogen_species_id
       AND ps.species_key = :species_key
       AND gc.collection_key = :collection_key
       AND gc.collection_label = :only_if_label_is
"""


def _relabel(species_key: str, collection_key: str, new_label: str, only_if_label_is: str) -> None:
    import sqlalchemy as sa

    op.get_bind().execute(
        sa.text(_RELABEL),
        {
            "new_label": new_label,
            "species_key": species_key,
            "collection_key": collection_key,
            "only_if_label_is": only_if_label_is,
        },
    )


def upgrade() -> None:
    for species_key, collection_key, label in COHORT_LABELS:
        # ⚠ `only_if_label_is` is the KEY — the placeholder — so a chosen name survives.
        _relabel(species_key, collection_key, label, collection_key)


def downgrade() -> None:
    """Back to the placeholder, which is a true inverse here because the placeholder IS the key.

    ⚠ Only for a row still wearing the label this migration wrote. A name chosen since is left
    alone, for the same reason `upgrade` will not overwrite one.
    """
    for species_key, collection_key, label in COHORT_LABELS:
        _relabel(species_key, collection_key, collection_key, label)
