"""Fixtures for the backend suite.

The minimal chain species → collection → model → pangenome → locus, built once per module. It
exists because the constraints worth testing are on the REAL tables: a hand-written temp table with
a lookalike CHECK proves only that Postgres can evaluate a CHECK, which was never in doubt.
"""

import json
import os
from pathlib import Path

import pytest
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.orm import Session

import bacatlas_backend.models  # noqa: F401  (registers every table on the metadata)
from bacatlas_backend.database import Base
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.published_catalogues import catalogue as published_catalogue
from bacatlas_backend.models.enumerations import (
    PrevalenceBand,
    RosterLineage,
    SampleIdentifierKind,
)
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollection
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies

PROBE_URL = os.environ.get(
    "BACATLAS_PROBE_DATABASE_URL",
    f"postgresql+psycopg://{os.environ.get('USER', 'postgres')}@localhost:5432/bacatlas_schema_probe",
)

#: The ten signed offsets, `0` deliberately absent.
TEN_ZEROS = [0] * 10


@pytest.fixture(scope="module")
def engine():
    engine = create_engine(PROBE_URL, future=True)
    try:
        with engine.connect() as connection:
            connection.execute(text("SELECT 1"))
    except Exception as error:  # noqa: BLE001
        pytest.skip(f"probe database unavailable at {PROBE_URL}: {error}")
    # ⛔⛔ **`drop_all` CANNOT SHED A TABLE THE MODELS NO LONGER DECLARE, and that is a trap with
    # a long fuse.** It drops exactly what `Base.metadata` holds, so the moment a table is DELETED
    # from the models its row in the probe database becomes an orphan no run can ever remove.
    # Measured when `nuna_model_step` went (2026-10-27): the orphan still depended on the
    # `embeddingrepresentation` enum, so the next `create_all` died on `cannot drop type
    # embeddingrepresentation because other objects depend on it` — which reads as "the probe
    # database is broken", not "a table you removed is still there". Orphaned enum TYPES are the
    # same shape: `drop_all` removes the ones the metadata still knows about and leaves the rest,
    # and `create_all` then cannot recreate one that is already there.
    #
    # ⚠ **Surgical, not `DROP SCHEMA public CASCADE`.** That was the first fix and it was worse:
    # this fixture is MODULE-scoped and several modules share one probe database, so a statement
    # that removes everything makes the fixture destructive across modules. It produced eight
    # `relation "pathogen_species" does not exist` errors in one full run and none in the next —
    # and a flaky probe is a worse failure than the one being fixed, because it teaches people to
    # re-run. This touches only tables and types the models no longer declare, so for a module
    # with no orphans it does exactly what `drop_all` always did.
    with engine.begin() as connection:
        # ⛔ `pg_trgm` BEFORE `create_all`: `ix_locus__search_text_trigram` uses `gin_trgm_ops`,
        # and without the extension the whole schema build fails on that one index — which reads
        # as "the probe database is broken" rather than "an extension is missing".
        connection.execute(text("CREATE EXTENSION IF NOT EXISTS pg_trgm"))

    declared_tables = set(Base.metadata.tables)
    orphan_tables = [
        name
        for name in inspect(engine).get_table_names()
        if name not in declared_tables and name != "alembic_version"
    ]
    with engine.begin() as connection:
        for name in orphan_tables:
            connection.execute(text(f'DROP TABLE IF EXISTS "{name}" CASCADE'))

    Base.metadata.drop_all(engine)

    # ⚠ AFTER `drop_all`, so every type the models still declare has gone with its columns and
    # whatever enum types remain are by definition orphans of a deleted model.
    with engine.begin() as connection:
        remaining = connection.execute(
            text("SELECT typname FROM pg_type WHERE typtype = 'e'")
        ).scalars().all()
        for name in remaining:
            connection.execute(text(f'DROP TYPE IF EXISTS "{name}" CASCADE'))

    Base.metadata.create_all(engine)
    return engine


@pytest.fixture()
def session(engine):
    with Session(engine) as session:
        yield session
        session.rollback()


@pytest.fixture(scope="module")
def _seed_ids(engine):
    """Commit the shared seed ONCE per module, and hand back plain ids.

    ⚠ Ids, not ORM objects: an object created in a closed session is detached, and every attribute
    access then raises somewhere unrelated to the test that is actually failing.
    """
    with Session(engine) as session:
        ids = _build_seed(session)
        session.commit()
    return ids


@pytest.fixture()
def seeded(session, _seed_ids):
    """The shared seed, re-attached to this test's session."""
    return {
        "species": session.get(PathogenSpecies, _seed_ids["species"]),
        "genome": session.get(Genome, _seed_ids["genome"]),
        "collection": session.get(GenomeCollection, _seed_ids["collection"]),
        "pangenome": session.get(Pangenome, _seed_ids["pangenome"]),
    }


def _build_seed(session):
    """One species, one genome, one collection, one pangenome.

    ⚠ There is no model row any more: `nuna_model` and `nuna_model_step` were dropped with the
    model axis (2026-10-27). A pangenome is one RUN over one cohort, and the model it used is the
    free-text `model_label` below.
    """
    species = PathogenSpecies(species_key="probe", scientific_name="Escherichia coli")
    session.add(species)
    session.flush()

    genome = Genome(
        pathogen_species_id=species.pathogen_species_id,
        sample_id="SAMEA0000001",
        sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
        strand_is_observed=True,
    )
    collection = GenomeCollection(
        pathogen_species_id=species.pathogen_species_id,
        # ⚠ No underscore: the collection key is now half of `catalogue_key`, which must be
        # prefix-free, and `ingest_genome_collection` runs `catalogue_key_for`'s guards on it.
        collection_key="probe100",
        # ⚠ Deliberately NOT equal to the key: `offer_catalogue` refuses a cohort still wearing
        # its key, and a fixture that tripped that guard would make every picker test test the guard.
        collection_label="E. coli probe",
        roster_lineage=RosterLineage.PROBE_BIOSAMPLE,
        genome_count=100,
    )
    session.add_all([genome, collection])
    session.flush()

    pangenome = Pangenome(
        run_id="probe_run_res0.1_seed0",
        # The key a reader addresses this catalogue by — `{species}-{model}`, hyphenated so the keys
        # stay prefix-free. The seed species is `probe`.
        catalogue_key="probe-nuna4",
        model_label="nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
        ingest_generation=1,
        pathogen_species_id=species.pathogen_species_id,
        genome_collection_id=collection.genome_collection_id,
        genome_count=100,
        gene_count=489146,
        locus_count=1,
    )
    session.add(pangenome)
    session.flush()
    return {
        "species": species.pathogen_species_id,
        "genome": genome.genome_id,
        "collection": collection.genome_collection_id,
        "pangenome": pangenome.pangenome_id,
    }


def make_locus(seeded, *, ordinal=0, label="0", **overrides):
    """A minimally-valid Locus, so a test can vary exactly the column it is about."""
    fields = dict(
        pangenome_id=seeded["pangenome"].pangenome_id,
        pathogen_species_id=seeded["species"].pathogen_species_id,
        node_label=label,
        catalogue_ordinal=ordinal,
        member_gene_count=99,
        member_genome_count=97,
        prevalence_band=PrevalenceBand.CORE,
        display_name="wzi",
        display_name_source="bakta_symbol",
        named_member_count=63,
        context_observed_member_counts=TEN_ZEROS,
        total_arrangement_count=4,
    )
    fields.update(overrides)
    return Locus(**fields)


# ── the local artifact mirror ──────────────────────────────────────────────────────────────────
#: Where the cluster artifacts are mirrored on this machine. ⚠ Overridable, because a machine that
#: does not have them must SKIP rather than fail — the mirror is ~150 MB of pulled files and is not
#: a checked-in fixture.
NUNA_DATA_ROOT = Path(os.environ.get("BACATLAS_NUNA_DATA_ROOT", "~/developer/nuna/data")).expanduser()

#: ⭐ The SITE catalogue, not the export. `data/{species}.json` in this repo is what `render_site`
#: shipped, so it alone carries `meta.landing`, `meta.examples` and the vendored `pfam_names` —
#: `render` MUTATES the payload and the site build is the only caller that serialises afterwards.
#: The export under `data/browser/` has none of the three, which is why it cannot be this oracle.
PUBLISHED_SITE_CATALOGUE_DIR = Path(__file__).resolve().parents[2] / "data"

requires_artifacts = pytest.mark.skipif(
    not (NUNA_DATA_ROOT / "proc" / "embeddings" / "meta").is_dir(),
    reason=f"the cluster artifact mirror is not present at {NUNA_DATA_ROOT}",
)


def artifacts_for(species_key: str) -> "CatalogueArtifacts":
    """The published catalogue's artifacts, addressed from the checked-in triple — never guessed."""
    entry = published_catalogue(species_key)
    return CatalogueArtifacts(
        data_root=NUNA_DATA_ROOT,
        set_key=entry.set_key,
        model_label=entry.model_label,
        run_id=entry.run_id,
        species_key=entry.species_key,
    )


@pytest.fixture(scope="session")
def ecoli_artifacts():
    if not (NUNA_DATA_ROOT / "proc" / "embeddings" / "meta").is_dir():
        pytest.skip(f"the cluster artifact mirror is not present at {NUNA_DATA_ROOT}")
    return artifacts_for("ecoli")


@pytest.fixture(scope="session")
def published_ecoli_site_catalogue():
    """The shipped `ecoli.json`, decoded — the acceptance oracle for everything derived."""
    path = PUBLISHED_SITE_CATALOGUE_DIR / "ecoli.json"
    if not path.exists():
        pytest.skip(f"the published site catalogue is not present at {path}")
    return json.loads(path.read_text())


@pytest.fixture(scope="session")
def exported_ecoli_payload():
    """⭐ The EXPORT — what `build_payload` emitted, and the only correct oracle for a rebuild.

    ⛔ Not `published_ecoli_site_catalogue` above. That is this same payload *after* `render_page`
    mutated it: it carries three extra blocks (`cog_names`, `go_names`, `pfam_names`) and two extra
    `meta` keys (`landing`, `examples`) which `build_payload` never produced. Diffing a rebuild
    against it reports five differences that are all correct behaviour, and no real signal.
    """
    if not (NUNA_DATA_ROOT / "browser").is_dir():
        pytest.skip(f"the browser export mirror is not present at {NUNA_DATA_ROOT / 'browser'}")
    return json.loads(artifacts_for("ecoli").published_payload.read_text())


@pytest.fixture(scope="session")
def exported_kp_payload():
    """The kp export — the second species, so a rule cannot be fitted to one catalogue."""
    if not (NUNA_DATA_ROOT / "browser").is_dir():
        pytest.skip(f"the browser export mirror is not present at {NUNA_DATA_ROOT / 'browser'}")
    return json.loads(artifacts_for("kp").published_payload.read_text())


# ── ⛔ silent gates ────────────────────────────────────────────────────────────────────────────
#
# **A check that skips has not passed, and a suite that prints `36 skipped` has not said which 36.**
# This matters more the further the real work moves off this machine: a laptop reduced to a conduit
# for a dump still runs the suite, and that suite can report green having checked the parity
# oracles, the route costs and the column-meaning audit *not at all*. Two mechanisms here:
#
# * `require_oracle` turns an absent oracle into a FAILURE, not a skip. An oracle is the only
#   independent witness a piece of code has; losing it silently is the whole failure mode.
# * `pytest_terminal_summary` prints a CENSUS of every skip, grouped by reason, on every run. The
#   21 files gating on `BACATLAS_DATABASE_URL` and the 9 on the artifact mirror then announce
#   themselves without 30 edits, and a number that grows is visible the day it grows.

#: ⭐ The explicit, typed opt-out. ⛔ Default OFF: a machine that genuinely cannot run an oracle
#: says so out loud, once, rather than every suite everywhere quietly measuring less.
ALLOW_UNORACLED = os.environ.get("BACATLAS_ALLOW_UNORACLED", "").strip().lower() in {"1", "true", "yes"}


def require_oracle(condition: bool, what: str, *, supplies: str) -> None:
    """Fail — not skip — when an independent oracle is unavailable, naming what it would have checked.

    ``what`` is the thing that is missing, in the imperative ("node is not on PATH"); ``supplies``
    is what the suite therefore did not check. Both appear in the message, because *not looked at*
    and *no difference* are indistinguishable in a test report unless you make them different.
    """
    if condition:
        return
    message = (
        f"ORACLE UNAVAILABLE: {what}. This test is the independent witness for {supplies}, so the "
        "suite cannot report on it at all. Install what is missing, or set "
        "BACATLAS_ALLOW_UNORACLED=1 to accept a suite that does not check it."
    )
    if ALLOW_UNORACLED:
        pytest.skip(message)
    pytest.fail(message, pytrace=False)


def pytest_configure(config):  # noqa: ARG001
    """⭐ Default `BACATLAS_ROOT_GFF` to the mirror, because the census showed what it was costing.

    ⛔ **34 of this suite's 36 skips were the whole GFF-backed sequence path** — 22 API tests and 12
    sequence-endpoint tests — on a machine that had the GFFs all along, simply because one variable
    was unset. Nobody chose that; it is what `36 skipped` hides. Setting it here runs them, and
    they pass.

    ⚠ Only when the directory exists, and never over an explicit value: a machine without the
    mirror still skips, and the census then NAMES the path it did not check rather than leaving a
    reader to infer it from a number.
    """
    if os.environ.get("BACATLAS_ROOT_GFF"):
        return
    mirrored = NUNA_DATA_ROOT / "raw" / "gff"
    if mirrored.is_dir():
        os.environ["BACATLAS_ROOT_GFF"] = str(mirrored)
        print(f"BACATLAS_ROOT_GFF defaulted to the artifact mirror at {mirrored}")


def pytest_terminal_summary(terminalreporter, exitstatus, config):  # noqa: ARG001
    """⭐ Every skip, grouped by reason and counted — printed on every run, pass or fail.

    `-rs` would do this on request, and nobody requests it. The point is that the census is
    unconditional: a reader sees *what was not checked* beside what was.
    """
    skipped = terminalreporter.stats.get("skipped", [])
    if not skipped:
        return
    reasons: dict[str, int] = {}
    for report in skipped:
        # longrepr for a skip is (path, lineno, "Skipped: <reason>")
        raw = report.longrepr[2] if isinstance(report.longrepr, tuple) else str(report.longrepr)
        reasons[raw.removeprefix("Skipped: ").strip()] = reasons.get(raw.removeprefix("Skipped: ").strip(), 0) + 1
    terminalreporter.write_sep("-", f"what this run did NOT check ({len(skipped)} skipped)")
    for reason, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
        terminalreporter.write_line(f"  {count:4d}  {reason}")
