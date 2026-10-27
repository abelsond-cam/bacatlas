"""two options with the same name are not a picker

⛔ **THE COHORT AXIS CANNOT REPRESENT TWO CATALOGUES OVER ONE COHORT, and the live database holds
two pairs of them.** Measured on `bacatlas_ingest` the moment the picker started naming cohorts:

    pangenome 1  ecoli-nuna4  is_published=t  cohort "E. coli AMR"  (default)
    pangenome 4  ecoli-nuna5  is_published=t  cohort "E. coli AMR"
    pangenome 2  kp-nuna4     is_published=t  cohort "Kp AMR"       (default)
    pangenome 5  kp-nuna5     is_published=t  cohort "Kp AMR"

The picker would offer each species two options reading `E. coli AMR — 100 genomes`, distinguished
only by ` (default)` on one of them. That is not a degraded picker, it is one that lies: the thing
a reader is choosing between stops being named anywhere on the control.

⭐ **So the nuna5 pair is HIDDEN — `is_published = false` — and nothing else about it changes.**
`is_published` means *offer this in the picker*; `pathogen_species.default_pangenome_id` means
*serve this for a bare species request*. Neither nuna5 row is a default (checked before this was
written), so hiding them strands nothing: both stay loaded, both stay addressable by
`?catalogue=ecoli-nuna5`, and every locus of theirs still resolves. Only the menu changes.

⚠ **This is editorial, and it is in a migration because the schema change forces it.** "Load it ·
offer it · serve it by default are three decisions" is this repo's rule, and the second normally
belongs to `offer_catalogue`, not to a migration. But the release that makes the picker show
cohorts is the release that makes these four rows ambiguous, so a database at this revision with
all four offered is internally inconsistent. The downgrade restores them exactly.

⚑ It is also what the models of record say: every four- and five-step chain — `nuna4`, `nuna5` and
their variants — is SUPERSEDED (`nuna/CLAUDE.md`, 2026-10-24), and the three cohorts at launch are
`ecoli-probe100` on nuna4, `kp-probe100` on hypergeo and `kp-71k` on hypergeo. No nuna5 among them.

⚠ Keyed on `catalogue_key`, so it is a no-op on a database that does not hold these two rows — a
fresh one, or one already re-ingested onto the cohort keys.

Revision ID: 2b1faa871117
Revises: 166560ca71e2
Create Date: 2026-10-27 23:25:29.723831
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision: str = '2b1faa871117'
down_revision: str | None = '166560ca71e2'
branch_labels = None
depends_on = None

#: The catalogues whose cohort is already named by another OFFERED catalogue of the same species.
HIDDEN = ("ecoli-nuna5", "kp-nuna5")

_SET_OFFERED = sa.text(
    "UPDATE pangenome SET is_published = :offered WHERE catalogue_key = ANY(:keys)"
)


def upgrade() -> None:
    op.get_bind().execute(_SET_OFFERED, {"offered": False, "keys": list(HIDDEN)})


def downgrade() -> None:
    """Back into the picker — a true inverse, because hiding wrote one boolean and nothing else.

    ⚠ It restores the ambiguity along with the rows, which is correct: a downgrade returns the
    database to the state the previous revision describes, and at that revision the picker named
    models, where four distinct options is exactly right.
    """
    op.get_bind().execute(_SET_OFFERED, {"offered": True, "keys": list(HIDDEN)})
