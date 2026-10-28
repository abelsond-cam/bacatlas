"""the cohort keys the declaration names

⛔ **`published_cohorts.py` is load-bearing now, and the database did not agree with it.** The two
collections were keyed `ecoli` and `kp` — SPECIES tokens, from the days when `collection_key`
defaulted to the dataset token. The declaration names them `probe100`, scoped to their species.

**Left as it was, the next `--stage pangenome` run would have CREATED A SECOND COLLECTION.**
`ingest_genome_collection` is idempotent on `(pathogen_species_id, collection_key)`, so resolving
`probe100` from the declaration against a row keyed `kp` finds nothing and inserts. The old row
stays behind — still carrying whichever catalogues were not re-ingested, still offered, still
labelled — and `ingest_pangenome_run` REPOINTS `pangenome.genome_collection_id` on the one that
was, so a cohort-by-cohort rollout ends with one species showing two options over two collections
that are the same hundred genomes. Making the declaration authoritative and leaving the keys stale
would have produced exactly the duplicate `2b1faa871117` had just finished removing.

⚠ **A plain UPDATE, and that is safe here for a reason worth stating.** `collection_key` is a
lookup handle, not a position: nothing stores it, no foreign key spells it, and the ordinal
contract lives in `genome_collection_membership`, which is untouched. The roster guard that refuses
a changed membership order is not involved — no membership row moves. ⛔ What is NOT safe is the
reverse shape, a key kept while the ROSTER changes, which is what that guard exists for.

⚠ `pangenome.catalogue_key` does NOT move. It stays `ecoli-nuna4` / `kp-nuna4` until 2.5d
re-ingests the cohorts, so every shared `?catalogue=` link keeps working and the six scripts that
hardcode those keys keep running. The wire's `cohort.key` does change, from `kp` to `probe100` —
a data change with no shape change, which is exactly the case the list routes' content-digest ETag
was chosen over a shape-only tag to cover.

Revision ID: 2ed24dc9ee53
Revises: 49f61160a9f5
Create Date: 2026-10-27 23:58:00.000000
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision: str = '2ed24dc9ee53'
down_revision: str | None = '49f61160a9f5'
branch_labels = None
depends_on = None

#: `(species_key, from, to)` — the declaration's `collection_key` for each live cohort.
REKEY = (
    ("ecoli", "ecoli", "probe100"),
    ("kp", "kp", "probe100"),
)

_REKEY = sa.text(
    """
    UPDATE genome_collection AS gc
       SET collection_key = :to_key
      FROM pathogen_species AS ps
     WHERE ps.pathogen_species_id = gc.pathogen_species_id
       AND ps.species_key = :species_key
       AND gc.collection_key = :from_key
    """
)


def _apply(pairs) -> None:
    bind = op.get_bind()
    for species_key, from_key, to_key in pairs:
        bind.execute(
            _REKEY, {"species_key": species_key, "from_key": from_key, "to_key": to_key}
        )


def upgrade() -> None:
    _apply(REKEY)


def downgrade() -> None:
    """Back to the species tokens — a true inverse, since only the handle moved."""
    _apply((species, to_key, from_key) for species, from_key, to_key in REKEY)
