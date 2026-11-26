"""The migration must build the models exactly, and must be reversible AND re-appliable.

⛔ Three properties, and each was broken on the first attempt — none of which reading the generated
file would have shown:

1. Autogenerate **omitted the circular foreign key** `pathogen_species.default_pangenome_id →
   pangenome`. It can only be created by ALTER after both tables exist, and Alembic dropped that
   edge on the floor. The column is the publish-flip pointer, so without the constraint nothing
   stops it naming a pangenome that does not exist.
2. `downgrade` dropped the tables but **not the eleven ENUM types**, so `downgrade` → `upgrade`
   died on `CREATE TYPE ... already exists`. A migration that cannot be re-applied is not
   reversible.
3. Nothing checks that the migration still matches the models unless something asks — so
   `test_no_drift` does.
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest
from sqlalchemy import create_engine, text

BACKEND_ROOT = Path(__file__).resolve().parents[1]
#: ⚠ The RUNNING interpreter's alembic, never a path into one particular venv. The old
#: `BACKEND_ROOT.parent / ".venv/bin/alembic"` silently SKIPPED this whole module from a git worktree
#: or a CI runner — neither has that directory — and "skipped" is not "the migration is fine".
#: Running `-m alembic` from `BACKEND_ROOT` also means `env.py` imports THIS checkout's models.
ALEMBIC = [sys.executable, "-m", "alembic"]
PROBE_URL = os.environ.get(
    "BACATLAS_MIGRATION_DATABASE_URL",
    f"postgresql+psycopg://{os.environ.get('USER', 'postgres')}@localhost:5432/bacatlas_migrate_probe",
)

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("alembic") is None, reason="alembic is not installed in this interpreter"
)


def _alembic(*args):
    return subprocess.run(
        [*ALEMBIC, *args],
        cwd=BACKEND_ROOT,
        env={**os.environ, "BACATLAS_DATABASE_URL": PROBE_URL},
        capture_output=True,
        text=True,
    )


@pytest.fixture(scope="module")
def probe_engine():
    engine = create_engine(PROBE_URL, future=True)
    try:
        with engine.connect() as connection:
            connection.execute(text("SELECT 1"))
    except Exception as error:  # noqa: BLE001
        pytest.skip(f"migration probe database unavailable: {error}")
    _alembic("downgrade", "base")
    return engine


def _count(engine, sql):
    with engine.connect() as connection:
        return connection.execute(text(sql)).scalar_one()


TABLES = "SELECT count(*) FROM information_schema.tables WHERE table_schema='public'"
ENUM_TYPES = (
    "SELECT count(*) FROM pg_type t JOIN pg_namespace n ON n.oid=t.typnamespace "
    "WHERE n.nspname='public' AND t.typtype='e'"
)


def test_upgrade_then_downgrade_then_upgrade_again(probe_engine):
    assert _alembic("upgrade", "head").returncode == 0
    after_first = _count(probe_engine, TABLES)
    assert after_first > 20, f"only {after_first} tables after upgrade"

    assert _alembic("downgrade", "base").returncode == 0
    # alembic_version survives a downgrade to base; everything else must not.
    assert _count(probe_engine, TABLES) == 1
    # ⛔ The defect this test exists for: eleven enum types survived here on the first attempt.
    assert _count(probe_engine, ENUM_TYPES) == 0, "downgrade left orphaned enum types"

    second = _alembic("upgrade", "head")
    assert second.returncode == 0, f"re-upgrade failed:\n{second.stderr[-1500:]}"
    assert _count(probe_engine, TABLES) == after_first


def test_the_circular_publish_pointer_foreign_key_exists(probe_engine):
    # ⛔ Autogenerate omitted this. It is the publish-flip pointer.
    _alembic("upgrade", "head")
    with probe_engine.connect() as connection:
        definition = connection.execute(
            text(
                "SELECT pg_get_constraintdef(c.oid) FROM pg_constraint c "
                "JOIN pg_class t ON t.oid = c.conrelid "
                "WHERE t.relname = 'pathogen_species' AND c.contype = 'f'"
            )
        ).scalars().all()
    assert any("pangenome" in d and "default_pangenome_id" in d for d in definition), (
        f"the publish-flip pointer has no referential integrity: {definition}"
    )
    # DEFERRABLE is required: the flip and the partition attach share one transaction.
    assert any("DEFERRABLE" in d for d in definition)


def test_the_served_annotation_constraints_autogenerate_cannot_see_exist_after_upgrade(probe_engine):
    """⛔ Hand-added in `493f46d0cc95`, and invisible to the drift check below: Alembic compares no CHECK on
    an existing table. Production databases are built from migrations, so this is the only place a missing
    or mistyped one would show (review of 9cfa481)."""
    _alembic("upgrade", "head")
    with probe_engine.connect() as connection:
        checks = dict(
            connection.execute(
                text(
                    "SELECT c.conname, pg_get_constraintdef(c.oid) FROM pg_constraint c "
                    "JOIN pg_class t ON t.oid = c.conrelid WHERE t.relname = 'pangenome' AND c.contype = 'c'"
                )
            ).all()
        )
        donor = connection.execute(
            text(
                "SELECT condeferrable, condeferred FROM pg_constraint "
                "WHERE conname = 'fk_served_locus_annotation__donor_locus_id__locus'"
            )
        ).one()
    assert "'^at-v[0-9]+\\.[0-9]+$'" in checks["ck_pangenome__annotation_model_version_is_tagged"]
    assert "annotation_model_version IS NOT NULL" in checks["ck_pangenome__served_state_names_its_version"]
    assert "served_annotation_sha256 IS NOT NULL" in checks["ck_pangenome__served_state_names_its_file"]
    # ⛔ Deferred to commit: a reload deleting donor and recipient together must not fail per row.
    assert tuple(donor) == (True, True)
    with probe_engine.begin() as connection:
        default = connection.execute(
            text(
                "SELECT column_default FROM information_schema.columns "
                "WHERE table_name = 'pangenome' AND column_name = 'served_annotation_state'"
            )
        ).scalar_one()
    assert default is not None and "NOT_LOADED" in default


def test_the_audit_fact_constraints_autogenerate_cannot_see_exist_after_upgrade(probe_engine):
    """⛔ Hand-added in `9601b9c2d400` (David's 2026-11-26 audit facts), and invisible to the drift check
    below for the same reason as the served-annotation ones: Alembic compares no CHECK on an existing
    table. A database built from migrations without them would store a share of 1.4 or a negative count.
    """
    _alembic("upgrade", "head")
    with probe_engine.connect() as connection:
        checks = dict(
            connection.execute(
                text(
                    "SELECT c.conname, pg_get_constraintdef(c.oid) FROM pg_constraint c "
                    "JOIN pg_class t ON t.oid = c.conrelid WHERE t.relname = 'locus' AND c.contype = 'c'"
                )
            ).all()
        )
    share = checks["ck_locus__mmseqs_group_share_is_a_share"]
    assert "mmseqs_group_share >= (0)::double precision" in share
    assert "mmseqs_group_share <= (1)::double precision" in share
    cosines = checks["ck_locus__esm_member_similarities_are_cosines"]
    for column in ("esm_member_similarity_p05", "esm_member_similarity_median"):
        assert f"{column} <= (1.001)::double precision" in cosines
        assert f"{column} >= ('-1.001'::numeric)::double precision" in cosines
    counts = checks["ck_locus__pfam_architecture_counts_are_counts"]
    assert "pfam_dominant_domain_gene_count >= 0" in counts
    assert "pfam_other_architecture_gene_count >= 0" in counts
    assert (
        "(pfam_dominant_domain_gene_count + pfam_other_architecture_gene_count) <= pfam_annotated_member_count"
        in checks["ck_locus__pfam_architectures_within_annotated"]
    )


def test_the_migration_still_matches_the_models(probe_engine):
    """Autogenerate against the migrated database must find NOTHING to do."""
    _alembic("upgrade", "head")
    result = _alembic("revision", "--autogenerate", "-m", "drift check")
    generated = None
    try:
        assert result.returncode == 0, result.stderr[-1000:]
        for line in result.stdout.splitlines():
            if line.startswith("Generating "):
                generated = Path(line.split("Generating ", 1)[1].split(" ...")[0].strip())
        assert generated and generated.exists(), f"no migration generated:\n{result.stdout}"
        body = generated.read_text()
        upgrade_body = body.split("def upgrade() -> None:")[1].split("def downgrade")[0]
        operations = [
            line.strip()
            for line in upgrade_body.splitlines()
            if line.strip().startswith("op.") or line.strip().startswith("sa.")
        ]
        assert not operations, (
            "the migration has drifted from the models — autogenerate wants:\n  "
            + "\n  ".join(operations)
        )
    finally:
        if generated and generated.exists():
            generated.unlink()


def test_the_ec_rename_preserves_the_values_it_retypes(probe_engine):
    """⭐ The claim the second migration actually makes, exercised with real rows.

    Autogenerate proposed drop-and-add for `ec_number` → `ec_numbers`. With no deployment yet the
    two are indistinguishable — which is exactly why the difference has to be pinned now rather
    than discovered at 406 M rows, where a drop-and-add is a table rewrite and a data loss.

    ⚠ Three values, three different fates: a single term, a multi-term set (the 12.3 % case that
    overflowed the old `String(32)`), and NULL — which must stay NULL and never become `{}`.
    """
    assert _alembic("upgrade", "3bb9b090284b").returncode == 0
    assert _alembic("downgrade", "167f47ef6e68").returncode == 0

    with probe_engine.connect() as connection:
        connection.execute(
            text(
                "INSERT INTO pathogen_species (species_key, scientific_name)"
                " VALUES ('ecrt', 'Escherichia coli') RETURNING pathogen_species_id"
            )
        )
        species_id = connection.execute(
            text("SELECT pathogen_species_id FROM pathogen_species WHERE species_key = 'ecrt'")
        ).scalar_one()
        genome_id = connection.execute(
            text(
                "INSERT INTO genome (pathogen_species_id, sample_id, sample_id_kind, strand_is_observed)"
                " VALUES (:s, 'SAMEA_EC_RENAME', 'BIOSAMPLE', true) RETURNING genome_id"
            ),
            {"s": species_id},
        ).scalar_one()
        connection.execute(
            text(
                "INSERT INTO gene_functional_annotation (genome_id, flat_index, ec_number) VALUES"
                " (:g, 0, '2.7.7.7'), (:g, 1, '1.1.1.1,2.7.7.7,3.4.21.-'), (:g, 2, NULL)"
            ),
            {"g": genome_id},
        )
        connection.commit()

    assert _alembic("upgrade", "3bb9b090284b").returncode == 0

    with probe_engine.connect() as connection:
        rows = dict(
            connection.execute(
                text("SELECT flat_index, ec_numbers FROM gene_functional_annotation ORDER BY flat_index")
            ).all()
        )
    assert rows[0] == ["2.7.7.7"]
    assert rows[1] == ["1.1.1.1", "2.7.7.7", "3.4.21.-"], "the multi-term set did not survive the cast"
    assert rows[2] is None, "NULL became an empty array — 'not annotated' is not 'annotated with nothing'"

    with probe_engine.connect() as connection:
        connection.execute(text("DELETE FROM gene_functional_annotation WHERE genome_id = :g"), {"g": genome_id})
        connection.execute(text("DELETE FROM genome WHERE genome_id = :g"), {"g": genome_id})
        connection.execute(text("DELETE FROM pathogen_species WHERE species_key = 'ecrt'"))
        connection.commit()
