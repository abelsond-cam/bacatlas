"""`GET /species/{key}/genomes` — the anchor control's list, and the two counts it carries.

⭐ It replaces three things the published page held resident — `meta.genomes`, `anchorSearch` and
`GENOME_N` — and the property worth pinning is that the COUNTS are read, never aggregated per
request: at the 80,000-genome design target the aggregate is over ~412 M membership rows, on every
keystroke. The route's statement count is pinned in `test_route_cost.py`; what the counts MEAN is
pinned here, and against the published page in `test_anchoring_parity.py` (T4).

⚠ The endpoint tests run against `BACATLAS_DATABASE_URL` with both catalogues loaded, published and
migrated to the table; the ingest tests at the bottom run on the schema probe, where a hand-built
catalogue can put a genome at ρ > 1 on purpose.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, func, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.ingest.ingest_genome_locus_counts import (
    GENOME_LOCUS_COUNT_SELECT,
    write_genome_locus_counts,
)
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.gene import Gene, GeneLocusMembership
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pangenome_genome_locus_count import PangenomeGenomeLocusCount
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from tests.conftest import make_locus


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
        counted = session.execute(
            select(func.count()).select_from(PangenomeGenomeLocusCount)
        ).scalar_one()
        if not counted:
            pytest.skip(f"{url} holds no per-genome locus counts — run `alembic upgrade head` there")
    return create_application(Configuration(database_url=url))


@pytest.fixture(scope="module")
def client(application):
    return application.test_client()


def _listing(client, species_key="ecoli", **query):
    response = client.get(f"/api/v1/species/{species_key}/genomes", query_string=query)
    assert response.status_code == 200, response.get_data(as_text=True)[:300]
    return response.get_json()


# ── the contract ───────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("species_key", ["ecoli", "kp"])
def test_an_empty_query_lists_EVERY_collection_genome_in_the_order_meta_genomes_published(
    client, species_key
):
    body = _listing(client, species_key, limit=1000)
    assert body["species_key"] == species_key
    assert body["query"] == ""
    assert body["genome_count"] == body["matched_genome_count"] == len(body["genomes"]) == 100
    assert body["truncated"] is False
    # ⛔ The ordinal is the contract `arr.gid` indexed; the list is in it, never in accession order.
    assert [g["collection_genome_ordinal"] for g in body["genomes"]] == list(range(100))
    # ⛔ An EXACT key set, not a subset: a row that quietly grew a field is a contract change, and
    # this assertion is what makes someone come here and say so. `source_group` was added
    # 2026-10-24 for the anchor control's two lists — `test_genome_source_groups.py` owns what it
    # MEANS; this line only records that the row carries it.
    assert set(body["genomes"][0]) == {
        "sample_id", "collection_genome_ordinal", "locus_count", "arrangement_locus_count",
        "source_group",
    }


def test_the_default_limit_is_100_and_a_cut_list_SAYS_it_was_cut(client):
    assert len(_listing(client)["genomes"]) == 100
    cut = _listing(client, limit=3)
    assert len(cut["genomes"]) == 3
    # ⛔ Matches BEFORE the limit. Without it three rows read as "three genomes match".
    assert cut["matched_genome_count"] == 100
    assert cut["truncated"] is True
    assert [g["collection_genome_ordinal"] for g in cut["genomes"]] == [0, 1, 2]


@pytest.mark.parametrize(("asked", "served"), [(0, 1), (-5, 1), (5000, 100), ("many", 100)])
def test_the_limit_is_HELD_to_between_one_and_a_thousand(client, asked, served):
    """A nonsense limit is clamped or defaulted, never an error and never an empty 200."""
    assert len(_listing(client, limit=asked)["genomes"]) == served


def test_the_query_is_a_CASE_INSENSITIVE_substring_and_is_echoed_trimmed(client):
    everything = _listing(client, limit=1000)["genomes"]
    probe = everything[37]["sample_id"]
    fragment = probe[3:9]
    expected = [g["sample_id"] for g in everything if fragment.upper() in g["sample_id"].upper()]
    for query in (fragment, fragment.lower(), f"  {fragment.upper()}  "):
        body = _listing(client, q=query, limit=1000)
        assert body["query"] == query.strip()
        assert [g["sample_id"] for g in body["genomes"]] == expected
        assert body["matched_genome_count"] == len(expected)
    assert 1 <= len(expected) < 100, f"{fragment!r} matched {len(expected)} — the filter is untested"


@pytest.mark.parametrize("query", ["%", "_", "SAMEA_", "%2204%"])
def test_LIKE_metacharacters_are_LITERAL_as_indexOf_reads_them(client, query):
    """⛔ `anchorSearch` is `indexOf`, where `%` and `_` are characters. No accession contains one."""
    body = _listing(client, q=query)
    assert body["genomes"] == [] and body["matched_genome_count"] == 0
    assert body["truncated"] is False


def test_a_query_matching_nothing_is_an_EMPTY_200_that_still_names_the_collection_size(client):
    body = _listing(client, q="NOT-A-GENOME")
    assert body["genomes"] == []
    assert body["matched_genome_count"] == 0
    assert body["genome_count"] == 100


@pytest.mark.parametrize("species_key", ["tuberculosis", "ECOLI"])
def test_an_unknown_species_is_a_NAMED_404_and_never_an_empty_200(client, species_key):
    response = client.get(f"/api/v1/species/{species_key}/genomes")
    assert response.status_code == 404
    assert response.get_json()["error"] == "not_found"
    assert species_key in response.get_json()["detail"]


def test_the_listing_is_cacheable_and_tagged_by_its_pangenome_whatever_the_query(client):
    """⚠ Two queries, two bodies, ONE tag — and that is right, not an oversight.

    An entity tag distinguishes representations of one resource, and the query string is part of the
    resource's URI (RFC 9110 §8.8.3): each query is cached and revalidated under its own URI, and for
    any one URI the body changes only with the pangenome.
    """
    first = client.get("/api/v1/species/ecoli/genomes", query_string={"q": "SAMEA22"})
    second = client.get("/api/v1/species/ecoli/genomes", query_string={"q": "SAMEA23"})
    assert first.headers["Cache-Control"] == second.headers["Cache-Control"] == "public, max-age=86400"
    assert first.get_json() != second.get_json()
    assert first.headers["ETag"] == second.headers["ETag"]
    species = client.get("/api/v1/species/ecoli")
    assert first.headers["ETag"] == species.headers["ETag"]


# ── the two counts ─────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("species_key", ["ecoli", "kp"])
def test_the_two_counts_are_DIFFERENT_facts_and_differ_for_every_probe_genome(client, species_key):
    """⛔ Loci with a gene, and loci with a recorded neighbourhood — the second is the page's number.

    A gene alone on its contig is counted present and reaches no window, so the arranged loci are a
    strict subset of the present ones. Measured: every one of the 100 genomes in both species has
    such a gene somewhere, so collapsing the two fields would be wrong for the whole picker.
    """
    genomes = _listing(client, species_key, limit=1000)["genomes"]
    assert len(genomes) == 100
    assert all(g["arrangement_locus_count"] <= g["locus_count"] for g in genomes)
    differing = [g for g in genomes if g["arrangement_locus_count"] < g["locus_count"]]
    assert len(differing) == 100


def test_the_per_genome_counts_are_COUNT_DISTINCT_and_add_up_to_the_locus_side(application):
    """⛔ `COUNT(*)` where `COUNT(DISTINCT locus)` is meant gives each genome's GENE total instead.

    The two sides of one bipartite count must agree: summed over genomes, "loci this genome is in"
    is exactly "genomes this locus is in" summed over loci — `locus.member_genome_count`, written by a
    different code path from different frames. A row count would exceed it by every ρ > 1 gene.
    """
    engine = create_engine(application.config["BACATLAS"].database_url, future=True)
    with Session(engine) as session:
        for pangenome_id in session.execute(
            select(PathogenSpecies.default_pangenome_id)
        ).scalars():
            per_genome = session.execute(
                select(func.sum(PangenomeGenomeLocusCount.locus_count)).where(
                    PangenomeGenomeLocusCount.pangenome_id == pangenome_id
                )
            ).scalar_one()
            per_locus = session.execute(
                select(func.sum(Locus.member_genome_count)).where(Locus.pangenome_id == pangenome_id)
            ).scalar_one()
            gene_rows = session.execute(
                select(func.count())
                .select_from(GeneLocusMembership)
                .where(GeneLocusMembership.pangenome_id == pangenome_id)
            ).scalar_one()
            assert per_genome == per_locus, (pangenome_id, per_genome, per_locus)
            # And the naive count really is different here, so the assertion above has teeth.
            assert gene_rows > per_locus, (pangenome_id, gene_rows, per_locus)


def test_the_stored_counts_are_the_ingest_definition_recomputed(application):
    """⭐ The migration's backfill and the ingest statement are two copies; this holds them together.

    The rows in the loaded database came from the migration's frozen backfill; a fresh ingest writes
    them with `GENOME_LOCUS_COUNT_SELECT`. Recomputing the latter over the same rows must reproduce
    the former exactly, genome for genome, on both species.
    """
    engine = create_engine(application.config["BACATLAS"].database_url, future=True)
    with Session(engine) as session:
        examined = 0
        for pangenome_id in session.execute(
            select(PathogenSpecies.default_pangenome_id)
        ).scalars():
            stored = {
                genome_id: (present, arranged)
                for genome_id, present, arranged in session.execute(
                    select(
                        PangenomeGenomeLocusCount.genome_id,
                        PangenomeGenomeLocusCount.locus_count,
                        PangenomeGenomeLocusCount.arrangement_locus_count,
                    ).where(PangenomeGenomeLocusCount.pangenome_id == pangenome_id)
                )
            }
            recomputed = {
                genome_id: (present, arranged)
                for genome_id, present, arranged in session.execute(
                    text(GENOME_LOCUS_COUNT_SELECT), {"pangenome_id": pangenome_id}
                )
            }
            assert stored == recomputed, pangenome_id
            examined += len(stored)
    assert examined == 200, f"examined {examined} genomes"


@pytest.mark.parametrize("pangenome_id", [1, 2])
def test_the_publish_gate_NAMES_the_genome_count_check_and_it_passes(application, pangenome_id):
    from bacatlas_backend.ingest.publish_pangenome import verify_pangenome_is_servable

    engine = application.extensions["bacatlas_database"].engine
    with Session(engine) as session:
        passed, failed = verify_pangenome_is_servable(session, session.get(Pangenome, pangenome_id))
    assert failed == []
    assert "every collection genome has a locus count row" in passed


# ── the ingest statement, on a catalogue built to break it ─────────────────────────────────────
@pytest.fixture()
def small_catalogue(session, seeded):
    """Three genomes, two loci, and every shape the count can get wrong.

    * genome A: TWO genes at locus 0 sitting in TWO arrangements there (ρ > 1), plus a gene at
      locus 1 in an arrangement.
    * genome B: a gene at locus 0 in an arrangement, and a gene at locus 1 with NO window.
    * genome C: in the collection, and no gene anywhere.
    """
    species_id = seeded["species"].pathogen_species_id
    pangenome = seeded["pangenome"]
    genome_a = seeded["genome"]
    genome_b = Genome(
        pathogen_species_id=species_id, sample_id="SAMEA0000002",
        sample_id_kind=SampleIdentifierKind.BIOSAMPLE, strand_is_observed=True,
    )
    genome_c = Genome(
        pathogen_species_id=species_id, sample_id="SAMEA0000003",
        sample_id_kind=SampleIdentifierKind.BIOSAMPLE, strand_is_observed=True,
    )
    session.add_all([genome_b, genome_c])
    session.flush()
    for ordinal, genome in enumerate((genome_a, genome_b, genome_c)):
        session.add(
            GenomeCollectionMembership(
                genome_collection_id=pangenome.genome_collection_id,
                genome_id=genome.genome_id,
                collection_genome_ordinal=ordinal,
            )
        )
    loci = [make_locus(seeded, ordinal=index, label=str(index)) for index in (0, 1)]
    session.add_all(loci)
    session.flush()
    a, b = genome_a.genome_id, genome_b.genome_id
    arrangements = [
        LocusArrangement(
            locus_id=loci[0].locus_id, pangenome_id=pangenome.pangenome_id, rank_within_locus=0,
            member_gene_count=2, member_genome_count=2, neighbour_slot_codes=[-1] * 10,
            member_genome_ids=[a, b],
        ),
        LocusArrangement(
            locus_id=loci[0].locus_id, pangenome_id=pangenome.pangenome_id, rank_within_locus=1,
            member_gene_count=1, member_genome_count=1, neighbour_slot_codes=[-1] * 10,
            member_genome_ids=[a],
        ),
        LocusArrangement(
            locus_id=loci[1].locus_id, pangenome_id=pangenome.pangenome_id, rank_within_locus=0,
            member_gene_count=1, member_genome_count=1, neighbour_slot_codes=[-1] * 10,
            member_genome_ids=[a],
        ),
    ]
    session.add_all(arrangements)
    session.flush()
    # ⚠ No arrangement is paired with a gene here. `arrangement_locus_count` is counted from
    # `locus_arrangement.member_genome_ids` above, never from a per-gene link — which is why
    # dropping `gene_locus_membership.locus_arrangement_id` left this fixture's subject intact.
    genes = [
        # (genome, flat_index, locus)
        (a, 0, loci[0]),
        (a, 1, loci[0]),
        (a, 2, loci[1]),
        (b, 0, loci[0]),
        (b, 1, loci[1]),
    ]
    for genome_id, flat_index, _locus in genes:
        session.add(
            Gene(
                genome_id=genome_id, flat_index=flat_index, pathogen_species_id=species_id,
                contig_index=0, start_position=1 + 100 * flat_index,
                end_position=90 + 100 * flat_index, length_nt=90, strand="+",
            )
        )
    session.flush()
    for genome_id, flat_index, locus in genes:
        session.add(
            GeneLocusMembership(
                pangenome_id=pangenome.pangenome_id, genome_id=genome_id, flat_index=flat_index,
                locus_id=locus.locus_id,
            )
        )
    session.flush()
    return pangenome, {"a": a, "b": b, "c": genome_c.genome_id}


def test_the_ingest_counts_DISTINCT_loci_so_rho_above_one_is_counted_once(session, small_catalogue):
    pangenome, genome = small_catalogue
    written = write_genome_locus_counts(session, pangenome.pangenome_id)
    rows = {
        row.genome_id: (row.locus_count, row.arrangement_locus_count)
        for row in session.execute(
            select(PangenomeGenomeLocusCount).where(
                PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id
            )
        ).scalars()
    }
    assert written == 3
    # ⛔ A: three genes and three arrangement memberships, over TWO loci. `count(*)` would say 3 / 3.
    assert rows[genome["a"]] == (2, 2)
    # B: present at both loci, arranged at one — the gene with no window is a present locus only.
    assert rows[genome["b"]] == (2, 1)
    # ⚠ C: in the collection and in nothing. A row of zeros, never a missing row.
    assert rows[genome["c"]] == (0, 0)


def test_the_publish_gate_REFUSES_a_collection_genome_with_no_count_row(session, small_catalogue):
    from bacatlas_backend.ingest.publish_pangenome import verify_pangenome_is_servable

    pangenome, genome = small_catalogue
    write_genome_locus_counts(session, pangenome.pangenome_id)
    passed, _ = verify_pangenome_is_servable(session, pangenome)
    assert "every collection genome has a locus count row" in passed

    session.execute(
        PangenomeGenomeLocusCount.__table__.delete().where(
            PangenomeGenomeLocusCount.genome_id == genome["c"]
        )
    )
    passed, failed = verify_pangenome_is_servable(session, pangenome)
    assert "every collection genome has a locus count row" not in passed
    assert any(
        failure.startswith("every collection genome has a locus count row: 1 collection genomes")
        for failure in failed
    ), failed
