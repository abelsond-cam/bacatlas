"""The constraints that encode a trap must fire — and must not fire on valid data.

⛔ A guard that never fires is worse than no guard: it reads as protection. Every constraint here is
tested in BOTH directions on the REAL tables, because a hand-written temp table with a lookalike
CHECK proves only that Postgres can evaluate a CHECK, which was never in doubt.
"""

import pytest
from sqlalchemy import text
from sqlalchemy.exc import IntegrityError

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import RhoRule
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


# ⛔ **THE TWO `exclusivity_token_agrees` TESTS ARE DELETED WITH THE CONSTRAINT (2026-10-27)**, and
# the loss is named rather than quietly absorbed. They pinned that `-excl`, a PREFIX of
# `-exclLOGP`, could not be stored against the wrong decoded form — the trap four readers fell into
# in one day. Nothing in this service decodes a run_id any anymore, so there is no derived form to be
# inconsistent with: the id is stored verbatim and the token is inside it, where nuna's own suite
# owns the vocabulary. A test over a vocabulary this repo no longer reads is theatre.
#
# ⚠ The structurally identical trap that IS still live is in the CATALOGUE KEYS — `ecoli` globbing
# `locus_browser_ecoli_nuna5_*.json` — and `test_catalogue_addressing.py` is where it is pinned.

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


def test_pangenome_step_rho_ceiling_can_hold_step_twos_literal_1e9(session, seeded):
    """⛔ `nuna_pipeline.py` passes the string `"1e9"`, and every multi-step model's step 2 carries it.

    The test this replaced pinned the same value on `nuna_model_step`, whose `rho_ceiling` is
    `Numeric(18, 4)` — widened from `Numeric(12, 4)` after `numeric field overflow` would have
    failed the ingest of EVERY multi-step model. That table is gone with the registry, but the
    value is not: it arrives from the run's own manifest and lands in `pangenome_step`, which is
    kept. ⭐ That column is a `Float`, so the overflow hazard does not transfer — asserted rather
    than assumed, because "a different type, therefore fine" is the kind of reasoning this suite
    exists to replace with a measurement.
    """
    from bacatlas_backend.models.pangenome import PangenomeStep

    step = PangenomeStep(
        pangenome_id=seeded["pangenome"].pangenome_id,
        step_ordinal=2, step_name="step2", gamma=0.98,
        rho_rule=RhoRule.OFF, rho_ceiling=1e9, node_mode="mmseq", uses_exclusivity=False,
    )
    session.add(step)
    session.flush()
    session.expire(step)
    assert float(step.rho_ceiling) == 1e9


# ── the 2026-11-26 audit facts ─────────────────────────────────────────────────────────────────
# ⛔ Each CHECK is exercised in BOTH directions: a boundary value it must accept, and the value it
# exists to refuse. The three float columns carry no `nan_guards` — the range CHECK refuses a NaN
# because Postgres orders NaN above every number — so the NaN case is asserted here, per column.
@pytest.mark.parametrize("value", [0.0, 0.95, 1.0])
def test_the_mmseqs_group_share_ACCEPTS_a_share_including_both_ends(session, seeded, value):
    session.add(make_locus(seeded, ordinal=300, label="g0", mmseqs_group_share=value))
    session.flush()


@pytest.mark.parametrize("value", [1.0001, -0.0001, float("nan")])
def test_the_mmseqs_group_share_REJECTS_anything_that_is_not_a_share(session, seeded, value):
    session.add(make_locus(seeded, ordinal=301, label="g1", mmseqs_group_share=value))
    with pytest.raises(IntegrityError, match="mmseqs_group_share_is_a_share"):
        session.flush()


@pytest.mark.parametrize("value", [1.000001, -0.2, 0.0])
@pytest.mark.parametrize("column", ["esm_member_similarity_p05", "esm_member_similarity_median"])
def test_an_esm_member_similarity_ACCEPTS_a_cosine_with_float32_slack(session, seeded, column, value):
    # ⚠ 1.000001 is a CORRECT float32 cosine of two L2-normalised vectors (`COSINE_BOUND`'s note).
    session.add(make_locus(seeded, ordinal=302, label="e0", **{column: value}))
    session.flush()


@pytest.mark.parametrize("value", [1.01, -1.01, float("nan")])
@pytest.mark.parametrize("column", ["esm_member_similarity_p05", "esm_member_similarity_median"])
def test_an_esm_member_similarity_REJECTS_what_no_cosine_can_be(session, seeded, column, value):
    session.add(make_locus(seeded, ordinal=303, label="e1", **{column: value}))
    with pytest.raises(IntegrityError, match="esm_member_similarities_are_cosines"):
        session.flush()


@pytest.mark.parametrize("column", ["pfam_dominant_domain_gene_count", "pfam_other_architecture_gene_count"])
def test_a_pfam_architecture_count_REJECTS_a_negative(session, seeded, column):
    session.add(make_locus(seeded, ordinal=304, label="p0", **{column: -1}))
    with pytest.raises(IntegrityError, match="pfam_architecture_counts_are_counts"):
        session.flush()


def test_dominant_plus_other_may_EQUAL_the_annotated_genes(session, seeded):
    session.add(
        make_locus(
            seeded, ordinal=305, label="p1", pfam_annotated_member_count=90,
            pfam_dominant_domain_gene_count=82, pfam_other_architecture_gene_count=8,
        )
    )  # fmt: skip
    session.flush()


def test_dominant_plus_other_may_NOT_exceed_the_annotated_genes(session, seeded):
    # ⛔ An unannotated gene is "no Pfam", never "another architecture": counting it as the second
    # would push the sum past the annotated genes, which is exactly what this refuses.
    session.add(
        make_locus(
            seeded, ordinal=306, label="p2", pfam_annotated_member_count=90,
            pfam_dominant_domain_gene_count=82, pfam_other_architecture_gene_count=9,
        )
    )  # fmt: skip
    with pytest.raises(IntegrityError, match="pfam_architectures_within_annotated"):
        session.flush()
