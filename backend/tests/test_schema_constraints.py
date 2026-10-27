"""The constraints that encode a trap must fire — and must not fire on valid data.

⛔ A guard that never fires is worse than no guard: it reads as protection. Every constraint here is
tested in BOTH directions on the REAL tables, because a hand-written temp table with a lookalike
CHECK proves only that Postgres can evaluate a CHECK, which was never in doubt.
"""

import pytest
from sqlalchemy import text
from sqlalchemy.exc import IntegrityError

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import ExclusivityForm, ExclusivityFormSource
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant
from bacatlas_backend.models.pangenome import Pangenome

from .conftest import make_locus

OUR_TABLES = tuple(Base.metadata.tables)


def test_the_whole_schema_builds_in_postgres(engine):
    with engine.connect() as connection:
        present = {
            row[0]
            for row in connection.execute(
                text("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'")
            )
        }
    missing = set(OUR_TABLES) - present
    assert not missing, f"a model did not reach the database: {sorted(missing)}"


def test_postgres_really_does_say_nan_equals_nan(engine):
    # This is WHY the guard is `<> 'NaN'` and not the obvious IEEE-754 `col = col`. If this ever
    # stops being true the guard's rationale changes, and someone should notice here first.
    with engine.connect() as connection:
        assert connection.execute(text("SELECT 'NaN'::float8 = 'NaN'::float8")).scalar_one() is True
        assert connection.execute(text("SELECT 'NaN'::float8 <> 'NaN'::float8")).scalar_one() is False


@pytest.mark.parametrize("value", [0.0, 0.9436, 1e-9, 1.0])
def test_a_measurement_column_ACCEPTS_ordinary_values_INCLUDING_ZERO(session, seeded, value):
    # The half that matters most. 0.0 is a MEASUREMENT and must survive: a guard that rejected it
    # would be indistinguishable from one that worked, until a real measured zero arrived.
    session.add(make_locus(seeded, ordinal=100, label="m0", syntenic_a5=value))
    session.flush()
    assert session.query(type(session.query(Pangenome).first())).count() >= 1


def test_a_measurement_column_REJECTS_nan(session, seeded):
    session.add(make_locus(seeded, ordinal=101, label="m1", syntenic_a5=float("nan")))
    with pytest.raises(IntegrityError, match="syntenic_a5_is_not_nan"):
        session.flush()


def test_a_measurement_column_ACCEPTS_null_because_null_means_not_measured(session, seeded):
    # NULL is the whole reason these columns are nullable: `u50_impurity` is NULL below 5 labelled
    # members, and the seqid columns are NULL under --skip-seqid-to-medoid.
    session.add(make_locus(seeded, ordinal=102, label="m2", syntenic_a5=None, uniref50_impurity=None))
    session.flush()


def test_the_exclusivity_constraint_REJECTS_a_token_that_disagrees_with_its_form(session, seeded):
    # ⛔ `-excl` is a PREFIX of `-exclLOGP` and four readers got this wrong silently in one day.
    # After ingest no query looks at a run_id string again, because this row cannot exist.
    session.add(
        Pangenome(
            run_id="wrong_token_run",
            catalogue_key="probe-wrong-token",
            model_label="nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
            pathogen_species_id=seeded["species"].pathogen_species_id,
            genome_collection_id=seeded["collection"].genome_collection_id,
            exclusivity_form=ExclusivityForm.EXCLUSION,       # standard …
            exclusivity_form_source=ExclusivityFormSource.RUN_ID_TOKEN,
            run_id_exclusivity_token="-excl",                  # … but the DAMPED token
            genome_count=100,
            gene_count=1,
            locus_count=1,
        )
    )
    with pytest.raises(IntegrityError, match="exclusivity_token_agrees"):
        session.flush()


@pytest.mark.parametrize(
    ("token", "form"),
    [("-exclLOGP", ExclusivityForm.EXCLUSION), ("-excl", ExclusivityForm.DAMPED_EXCLUSION)],
)
def test_the_exclusivity_constraint_ACCEPTS_each_token_with_its_own_form(session, seeded, token, form):
    session.add(
        Pangenome(
            run_id=f"ok_{token}",
            catalogue_key=f"probe-ok{token}",
            model_label="nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
            pathogen_species_id=seeded["species"].pathogen_species_id,
            genome_collection_id=seeded["collection"].genome_collection_id,
            exclusivity_form=form,
            exclusivity_form_source=ExclusivityFormSource.RUN_ID_TOKEN,
            run_id_exclusivity_token=token,
            genome_count=100,
            gene_count=1,
            locus_count=1,
        )
    )
    session.flush()


def test_an_offset_of_zero_is_rejected(session, seeded):
    # `0` is the focal locus; a row at 0 would be a locus claiming to be its own neighbour.
    locus = make_locus(seeded, ordinal=200, label="z0")
    session.add(locus)
    session.flush()
    session.add(
        LocusOffsetOccupant(
            locus_id=locus.locus_id,
            pangenome_id=seeded["pangenome"].pangenome_id,
            signed_offset=0,
            rank_within_offset=0,
            neighbour_locus_id=locus.locus_id,
            member_gene_count=1,
            same_strand_member_count=1,
        )
    )
    with pytest.raises(IntegrityError, match="signed_offset_excludes_zero"):
        session.flush()


def test_our_constraint_names_all_come_from_the_convention(engine):
    # A regression guard for a real bug: an explicit short `name=` on UniqueConstraint BYPASSES the
    # naming convention, and two tables both got a constraint literally called `species_key`, which
    # Postgres rejected on the second CREATE TABLE. Scoped to OUR tables — pg_constraint is full of
    # system catalog entries that follow no convention of ours.
    with engine.connect() as connection:
        names = [
            row[0]
            for row in connection.execute(
                text(
                    "SELECT conname FROM pg_constraint c "
                    "JOIN pg_class t ON t.oid = c.conrelid "
                    "JOIN pg_namespace n ON n.oid = t.relnamespace "
                    "WHERE n.nspname = 'public' AND c.contype IN ('u','c','p','f')"
                )
            )
        ]
    assert names, "no constraints found — the scoping is wrong, not the schema"
    stray = [n for n in names if n.split("_")[0] not in {"uq", "ck", "pk", "fk"}]
    assert not stray, f"constraints bypassing the naming convention: {sorted(stray)}"


def test_there_is_NO_unique_constraint_on_locus_and_genome(engine):
    # ⛔ A genome at ρ > 1 occupies TWO arrangements at one locus. A uniqueness constraint on
    # (locus, genome) would assert something false about the biology. This test exists so nobody
    # adds one for tidiness.
    # ⚠ Compare COLUMN NAMES, not a substring of the definition: `"genome_id" in definition` is
    # true of `member_genome_ids` and this test passed for the wrong reason on the first attempt.
    with engine.connect() as connection:
        rows = connection.execute(
            text(
                "SELECT c.conname, array_agg(a.attname ORDER BY a.attname) FROM pg_constraint c "
                "JOIN pg_class t ON t.oid = c.conrelid "
                "JOIN unnest(c.conkey) AS k(attnum) ON true "
                "JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = k.attnum "
                "WHERE t.relname = 'locus_arrangement' AND c.contype = 'u' GROUP BY c.conname"
            )
        ).all()
    for name, columns in rows:
        assert not {"locus_id", "genome_id"} <= set(columns), (
            f"{name} forbids a genome occupying two arrangements at one locus: {columns}"
        )


def test_rho_ceiling_can_hold_step_twos_literal_1e9(session, seeded):
    """⛔ `nuna_pipeline.py:410` passes the string `"1e9"`, and CLAUDE.md's table says so too.

    At `Numeric(12, 4)` this raised `numeric field overflow`, so ingesting the step 2 of EVERY
    multi-step model would have failed. The value is stored, not merely accepted, because a silent
    truncation to `99999999.9999` would still describe a ρ rail — just not this one.
    """
    from bacatlas_backend.models.enumerations import RhoRule
    from bacatlas_backend.models.nuna_model import NunaModelStep

    step = NunaModelStep(
        nuna_model_id=seeded["model"].nuna_model_id,
        step_ordinal=2, step_name="step2", gamma=0.98,
        rho_rule=RhoRule.OFF, rho_ceiling=1e9, node_mode="mmseq", uses_exclusivity=False,
    )
    session.add(step)
    session.flush()
    session.expire(step)
    assert float(step.rho_ceiling) == 1e9
