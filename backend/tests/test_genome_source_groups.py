"""The anchor control's two lists — the derived source group, its census, and the filter.

⭐ **Why there are two lists.** The Klebsiella cohort holds two populations with different
properties: ~4,000 complete long-read NCBI genomes, short enough to browse, and ~75,000
AllTheBacteria short-read assemblies, which need a type-to-narrow search. One combobox over both
serves neither, so the picker offers a group and then a list within it.

⛔ **The group is DERIVED from `sample_id_kind` + `is_complete_assembly`, never stored.** No
migration, no re-ingest — both columns are already written by `build_genome_set`. The cost is that
the derivation has exactly one home (`SOURCE_GROUP_EXPRESSION`) and every test below goes through
it rather than restating the rule.

⚠ **The loaded cohort cannot exercise two of the three groups.** All 280 probe genomes are
BioSample-keyed ATB assemblies with `is_complete_assembly` NULL, so `ncbi_complete` and
`ncbi_draft` are empty there — measured, not assumed, by
`test_every_probe_genome_is_an_ATB_assembly`. The derivation tests therefore build genomes on the
schema probe, which is the only place an NCBI-keyed row exists at all today.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.pangenome_genome_locus_count import PangenomeGenomeLocusCount
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.services.genome_listing_service import (
    GenomeSourceGroup,
    count_genome_source_groups,
    list_genomes,
)


@pytest.fixture(scope="module")
def application():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the endpoint is tested on a loaded database")
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        published = session.execute(
            select(func.count())
            .select_from(PathogenSpecies)
            .where(PathogenSpecies.default_pangenome_id.is_not(None))
        ).scalar_one()
        if published < 2:
            pytest.skip(f"only {published} species are published in {url}")
        if not session.execute(select(func.count()).select_from(PangenomeGenomeLocusCount)).scalar_one():
            pytest.skip(f"{url} holds no per-genome locus counts — run `alembic upgrade head` there")
    return create_application(Configuration(database_url=url))


@pytest.fixture(scope="module")
def client(application):
    return application.test_client()


def _listing(client, species_key="ecoli", **query):
    response = client.get(f"/api/v1/species/{species_key}/genomes", query_string=query)
    assert response.status_code == 200, response.get_data(as_text=True)[:300]
    return response.get_json()


def _join_catalogue(session, pangenome, genome, *, ordinal: int) -> None:
    """Make `genome` a LISTABLE member of `pangenome`'s collection.

    ⛔ **Both rows, because the listing needs both.** Membership alone puts a genome in the
    collection; the `pangenome_genome_locus_count` row is what the listing joins for its two
    per-genome counts. A test that added only the first would build the exact state the publish gate
    forbids — and that state is what proved the census was counting the wrong population.
    """
    session.add(
        GenomeCollectionMembership(
            genome_collection_id=pangenome.genome_collection_id,
            genome_id=genome.genome_id,
            collection_genome_ordinal=ordinal,
        )
    )
    session.add(
        PangenomeGenomeLocusCount(
            pangenome_id=pangenome.pangenome_id,
            genome_id=genome.genome_id,
            locus_count=1,
            arrangement_locus_count=1,
        )
    )
    session.flush()


def _shell(client, species_key="ecoli"):
    response = client.get(f"/api/v1/species/{species_key}")
    assert response.status_code == 200
    return response.get_json()


# ── the derivation, on genomes the loaded cohort does not contain ──────────────────────────────
@pytest.mark.parametrize(
    ("sample_id", "kind", "is_complete", "expected"),
    [
        # A BioSample-keyed assembly is an ATB one — `NEXT_anchor_to_genome.md`: all 281 probe
        # manifest rows read `source = atb` and `ncbi_accession` is blank on every one.
        ("SAMEA9100001", SampleIdentifierKind.BIOSAMPLE, None, GenomeSourceGroup.ATB),
        # ⛔ And a BioSample row stays ATB even if the completeness flag is set, because the KEY
        # form is what says where the assembly came from. Reversing this precedence would file an
        # ATB assembly under "complete NCBI genomes", which is a claim about provenance.
        ("SAMEA9100002", SampleIdentifierKind.BIOSAMPLE, True, GenomeSourceGroup.ATB),
        ("GCF_910000001.1", SampleIdentifierKind.ASSEMBLY_STEM, True, GenomeSourceGroup.NCBI_COMPLETE),
        # ⚠ A `GCA_` draft long-read assembly: NCBI-keyed and NOT complete. `metadata_v2` holds
        # 5,519 LRA rows against 4,017 `is_complete`, so ~1,500 genomes land here — and they must
        # appear in a third group rather than vanish from a control that lists every genome.
        ("GCA_910000002.1", SampleIdentifierKind.ASSEMBLY_STEM, False, GenomeSourceGroup.NCBI_DRAFT),
        # ⛔ NULL is "not applicable / not recorded", NEVER false — but for GROUPING it cannot be
        # treated as complete either, so it falls to draft. Stated as a test because the nullable
        # boolean is exactly where a `== True` and an `is not False` silently disagree.
        ("GCF_910000003.1", SampleIdentifierKind.ASSEMBLY_STEM, None, GenomeSourceGroup.NCBI_DRAFT),
    ],
)
def test_the_source_group_is_DERIVED_from_the_key_form_and_the_completeness_flag(
    session, seeded, sample_id, kind, is_complete, expected
):
    pangenome = seeded["pangenome"]
    genome = Genome(
        pathogen_species_id=seeded["species"].pathogen_species_id,
        sample_id=sample_id,
        sample_id_kind=kind,
        is_complete_assembly=is_complete,
        strand_is_observed=True,
    )
    session.add(genome)
    session.flush()
    _join_catalogue(session, pangenome, genome, ordinal=0)

    counts = count_genome_source_groups(session, pangenome)

    assert counts[expected] == 1, f"{sample_id} was not grouped as {expected.value}: {counts}"
    # ⚠ Exactly one, because the shared seed's own genome is a collection member of a DIFFERENT
    # collection — asserted so a seed that later joined this one cannot make this test vacuous.
    assert sum(counts.values()) == 1


def test_the_census_reports_a_ZERO_for_a_group_with_no_members(session, seeded):
    """⛔ Zeros are PRESENT, not absent. An absent key cannot be told from one the server forgot to
    compute — and "this catalogue holds no complete NCBI genomes" is the answer today, so it is the
    sentence the control has to be able to say."""
    counts = count_genome_source_groups(session, seeded["pangenome"])

    assert set(counts) == set(GenomeSourceGroup)
    assert all(count == 0 for count in counts.values())


def test_the_row_and_the_FILTER_cannot_disagree(session, seeded):
    """⛔ One SQL expression serves the reported group, the filter and the census. A Python-side
    second copy of the rule is exactly how a filter and its tally drift apart."""
    pangenome = seeded["pangenome"]
    for ordinal, (sample_id, kind, complete) in enumerate(
        [
            ("SAMEA9200001", SampleIdentifierKind.BIOSAMPLE, None),
            ("GCF_920000001.1", SampleIdentifierKind.ASSEMBLY_STEM, True),
            ("GCA_920000002.1", SampleIdentifierKind.ASSEMBLY_STEM, False),
        ]
    ):
        genome = Genome(
            pathogen_species_id=seeded["species"].pathogen_species_id,
            sample_id=sample_id,
            sample_id_kind=kind,
            is_complete_assembly=complete,
            strand_is_observed=True,
        )
        session.add(genome)
        session.flush()
        _join_catalogue(session, pangenome, genome, ordinal=ordinal)

    census = count_genome_source_groups(session, pangenome)
    for group in GenomeSourceGroup:
        listing = list_genomes(session, pangenome, query="", limit=1000, source_group=group)
        # ⚠ `matched_genome_count` is the count BEFORE the limit, so it is the comparable number.
        assert listing.matched_genome_count == census[group], (
            f"{group.value}: the filter matched {listing.matched_genome_count} but the census "
            f"says {census[group]}"
        )
        assert all(row.source_group is group for row in listing.genomes)


# ── on the loaded catalogues ───────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("species_key", ["ecoli", "kp"])
def test_every_probe_genome_is_an_ATB_assembly(client, species_key):
    """⚠ **Measured, not assumed.** The probe cohort is entirely BioSample-keyed ATB assemblies, so
    the complete-genome list is empty on every catalogue loaded today. That is a fact about the
    cohort and the reason the NCBI group renders as a stated absence rather than an empty box — and
    it is the thing that will change the day the 4,017 complete genomes are ingested."""
    groups = _shell(client, species_key)["pangenome"]["genome_source_groups"]

    assert groups["ncbi_complete"] == 0
    assert groups["ncbi_draft"] == 0
    assert groups["atb"] > 0
    # ⛔ And the census adds up to the catalogue, so no genome is missing from every group.
    assert sum(groups.values()) == _listing(client, species_key)["genome_count"]


@pytest.mark.parametrize("species_key", ["ecoli", "kp"])
def test_every_listed_row_carries_its_group(client, species_key):
    rows = _listing(client, species_key, limit=25)["genomes"]

    assert rows
    assert {row["source_group"] for row in rows} == {"atb"}


def test_the_filter_narrows_the_list(client):
    everything = _listing(client, limit=1000)
    atb = _listing(client, limit=1000, source_group="atb")
    empty = _listing(client, limit=1000, source_group="ncbi_complete")

    assert atb["matched_genome_count"] == everything["matched_genome_count"]
    assert atb["source_group"] == "atb"
    assert empty["matched_genome_count"] == 0
    assert empty["genomes"] == []
    # ⛔ Still a 200 naming the collection size — an empty group is an ANSWER, not a 404.
    assert empty["genome_count"] == everything["genome_count"]


def test_an_unfiltered_listing_echoes_a_NULL_group(client):
    """⚠ Echoed so the client can discard a late answer to a superseded filter, exactly as `query`
    already is — the picker switches group and keystrokes race."""
    assert _listing(client, limit=1)["source_group"] is None


def test_the_filter_and_the_query_COMPOSE(client):
    """⚠ Both narrow, and neither replaces the other. A filter that silently dropped `q` would list
    every ATB genome under a query the reader had typed."""
    row = _listing(client, limit=1, source_group="atb")["genomes"][0]
    narrowed = _listing(client, limit=1000, source_group="atb", q=row["sample_id"])

    assert [hit["sample_id"] for hit in narrowed["genomes"]] == [row["sample_id"]]


def test_the_group_census_does_NOT_move_with_the_query(client):
    """⛔⛔ The distinction the control depends on: *"nothing matches what you typed"* against
    *"this catalogue holds none of these at all"*. The census is a property of the catalogue, so it
    rides on the shell and is independent of `q` by construction — this pins that it stays so."""
    groups = _shell(client)["pangenome"]["genome_source_groups"]

    # However specific the query, the census is unchanged — it is not computed from the listing.
    assert _shell(client)["pangenome"]["genome_source_groups"] == groups
    assert _listing(client, q="ZZZZZZZZ-no-such-genome")["matched_genome_count"] == 0


@pytest.mark.parametrize("value", ["nonsense", "NCBI_COMPLETE", "complete", "ncbi", ""])
def test_an_unknown_group_is_a_NAMED_400_and_never_a_silently_WIDER_list(client, value):
    """⛔ A typo'd filter that quietly widened the list would show complete-genome rows under an
    "AllTheBacteria" heading — a claim about provenance, made by a fallback.

    ⚠ `""` is the exception and is deliberately "no filter": an absent query parameter and an empty
    one arrive identically from a form, and refusing the empty string would break the unfiltered
    list the moment a client sent `?source_group=`.
    """
    response = client.get("/api/v1/species/ecoli/genomes", query_string={"source_group": value})

    if value == "":
        assert response.status_code == 200
        assert response.get_json()["source_group"] is None
    else:
        assert response.status_code == 400
        assert response.get_json()["error"] == "unknown_source_group"
        # The message names what IS allowed, so the fix does not need the source.
        assert "atb" in response.get_json()["detail"]


def test_the_shell_census_is_cacheable_with_the_catalogue(client):
    """⚠ It rides on the shell precisely because the shell is fetched once and is cacheable, where a
    census on the keystroke route would be re-run on every keypress.

    ⛔ Asserts the HEADERS, not a 304 round-trip: `_immutable` tags and sets `max-age` but does not
    call `make_conditional`, so revalidation is the cache's job and not this app's. Asserting a 304
    here would be asserting behaviour nothing implements — it failed for that reason when first
    written, which is the honest version of this test.
    """
    shell = client.get("/api/v1/species/ecoli")

    assert shell.headers["Cache-Control"] == "public, max-age=86400"
    assert shell.headers["ETag"]
    # The census is part of the tagged representation, so it cannot drift from the catalogue.
    assert "genome_source_groups" in shell.get_json()["pangenome"]
