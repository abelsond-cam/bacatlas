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
from bacatlas_backend.gene_blocks import encode_int32_vector
from bacatlas_backend.ingest.ingest_genome_locus_counts import (
    GENOME_LOCUS_COUNT_SELECT,
    write_genome_locus_counts,
)
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
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
                select(func.sum(GenomeLocusAssignment.placed_gene_count)).where(
                    GenomeLocusAssignment.pangenome_id == pangenome_id
                )
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
    a, b = genome_a.genome_id, genome_b.genome_id
    # ⚠ The counts are the TRUE ones for the genes below, so the publish gate's cross-table sums
    # (placed genes vs member genes, present genomes vs member_genome_count) hold on this fixture.
    loci = [
        make_locus(
            seeded, ordinal=index, label=str(index),
            member_gene_count=member_genes, member_genome_count=2, present_genome_ids=[a, b],
        )
        for index, member_genes in ((0, 3), (1, 2))
    ]
    session.add_all(loci)
    session.flush()
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
    # `locus_arrangement.member_genome_ids` above, never from a per-gene link.
    # ⭐ The genes, at the genome grain — `locus_count` is read from `distinct_locus_count`, and the
    # publish gate checks `gene_count` against `genome_genes`.
    # A: three genes over loci 0, 0, 1 → placed 3, distinct 2. B: 0, 1 → 2, 2. C: nothing.
    blocks = {
        a: ([0, 0, 1], 3, 2),
        b: ([0, 1], 2, 2),
        genome_c.genome_id: ([], 0, 0),
    }
    for genome_id, (ordinals, placed, distinct) in blocks.items():
        session.add(
            GenomeGenes(genome_id=genome_id, gene_count=len(ordinals), gene_spans=b"\x00" * (16 * len(ordinals)))
        )
        session.add(
            GenomeLocusAssignment(
                pangenome_id=pangenome.pangenome_id, genome_id=genome_id, gene_count=len(ordinals),
                placed_gene_count=placed, distinct_locus_count=distinct,
                gene_locus_ordinals=encode_int32_vector(ordinals),
                gene_arrangement_ranks=encode_int32_vector([-1] * len(ordinals)),
            )
        )
    session.flush()
    return pangenome, {"a": a, "b": b, "c": genome_c.genome_id}


def test_the_assignment_based_locus_count_counts_distinct_loci_genome_for_genome(session, small_catalogue):
    """⭐ The picker's count is DISTINCT loci per genome — ρ > 1 counted once, a genome with no gene a
    zero and not a gap. Coverage is asserted first: three genomes, before a value is compared."""
    pangenome, genome = small_catalogue
    from_assignments = {
        row.genome_id: row.locus_count
        for row in session.execute(text(GENOME_LOCUS_COUNT_SELECT), {"pangenome_id": pangenome.pangenome_id})
    }
    assert set(from_assignments) == set(genome.values())
    assert from_assignments == {genome["a"]: 2, genome["b"]: 2, genome["c"]: 0}


def test_a_genome_whose_assignment_row_is_missing_counts_ZERO_loci_not_a_missing_row(session, small_catalogue):
    """The join is OUTER: an inner join would write two rows and the picker would lose C silently.

    ⚠ C, not A: `ck_pangenome_genome_locus_count__arranged_within_present` refuses a genome that is
    arranged somewhere but assigned nowhere, which is the schema saying that state is not data.
    """
    pangenome, genome = small_catalogue
    session.execute(
        GenomeLocusAssignment.__table__.delete().where(GenomeLocusAssignment.genome_id == genome["c"])
    )
    written = write_genome_locus_counts(session, pangenome.pangenome_id)
    assert written == 3
    stored = session.execute(
        select(PangenomeGenomeLocusCount.locus_count, PangenomeGenomeLocusCount.arrangement_locus_count).where(
            PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id,
            PangenomeGenomeLocusCount.genome_id == genome["c"],
        )
    ).one()
    assert tuple(stored) == (0, 0)


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


# ── the publish gate on the genome-grain layer ────────────────────────────────────────────────
ASSIGNMENT_CHECKS = (
    "every collection genome has a locus assignment row",
    "every assignment's gene_count equals genome_genes.gene_count, or genome.coding_gene_count where the "
    "genome is named as holding none",
    "the placed genes add up to the loci's member genes",
    "every locus lists exactly member_genome_count present genomes",
)


def _gate(session, pangenome):
    from bacatlas_backend.ingest.publish_pangenome import verify_pangenome_is_servable

    return verify_pangenome_is_servable(session, pangenome)


def test_the_publish_gate_NAMES_all_four_assignment_checks_and_they_pass_on_a_consistent_fixture(
    session, small_catalogue
):
    pangenome, _genome = small_catalogue
    passed, failed = _gate(session, pangenome)
    for name in ASSIGNMENT_CHECKS:
        assert name in passed, (name, failed)


def test_the_gate_REFUSES_a_collection_genome_with_no_assignment_row(session, small_catalogue):
    pangenome, genome = small_catalogue
    session.execute(
        GenomeLocusAssignment.__table__.delete().where(GenomeLocusAssignment.genome_id == genome["c"])
    )
    _passed, failed = _gate(session, pangenome)
    assert any(
        f.startswith("every collection genome has a locus assignment row: 1 collection genomes have no")
        for f in failed
    ), failed


def test_the_gate_REFUSES_a_gene_count_that_disagrees_with_genome_genes(session, small_catalogue):
    pangenome, genome = small_catalogue
    # The genome stage re-annotated A with one more gene; the assignment block is now too short.
    session.execute(
        GenomeGenes.__table__.update().where(GenomeGenes.genome_id == genome["a"]).values(
            gene_count=4, gene_spans=b"\x00" * 64
        )
    )
    _passed, failed = _gate(session, pangenome)
    assert any(
        f.startswith(f"{ASSIGNMENT_CHECKS[1]}: 1 assignment rows disagree")
        for f in failed
    ), failed


def test_the_gate_REFUSES_when_the_placed_genes_do_not_add_up_and_STATES_both_numbers(session, small_catalogue):
    pangenome, _genome = small_catalogue
    session.execute(Locus.__table__.update().where(Locus.catalogue_ordinal == 0).values(member_gene_count=4))
    _passed, failed = _gate(session, pangenome)
    assert any(
        "sum(genome_locus_assignment.placed_gene_count) = 5 against sum(locus.member_gene_count) = 6" in f
        for f in failed
    ), failed


def test_the_gate_REFUSES_a_locus_whose_present_genomes_are_not_member_genome_count(session, small_catalogue):
    pangenome, genome = small_catalogue
    session.execute(
        Locus.__table__.update().where(Locus.catalogue_ordinal == 1).values(present_genome_ids=[genome["a"]])
    )
    _passed, failed = _gate(session, pangenome)
    assert any(
        f.startswith("every locus lists exactly member_genome_count present genomes: 1 loci have")
        for f in failed
    ), failed


# ── the gene layer in the publish gate: the blocks themselves ─────────────────────────────────
#: The three checks by name, as `verify_pangenome_is_servable` spells them. ⚠ `small_catalogue` seeds
#: the blocks these read — a: 3 genes, b: 2, c: 0 (a row, not a missing row) — so the total is 5, and
#: each test below tells the pangenome row the number it wants compared against.
#: ⚠ The first check names its own counts when it passes — all three genomes with a block, none named.
BLOCK_OR_NAMED = "every collection genome has a gene_spans block or is named as holding none"
BLOCK_TOTAL = "the gene_spans blocks, with the named genomes' coding_gene_count, hold the pangenome's gene_count"
GENE_LAYER_CHECKS = (
    f"{BLOCK_OR_NAMED} — 3 with a block, 0 named as carrying none, of 3",
    "every gene_spans block is 16 bytes per gene",
    BLOCK_TOTAL,
)


def _gene_layer_lines(passed, failed):
    return (
        [line for line in passed if line.startswith(GENE_LAYER_CHECKS) or "NOT compared" in line],
        [line for line in failed if line.startswith(GENE_LAYER_CHECKS)],
    )


def test_the_publish_gate_PASSES_the_gene_layer_when_the_blocks_add_up_to_the_pangenomes_gene_count(
    session, small_catalogue
):
    """The positive case, so a check that can only ever fail is told from one that is failing unnoticed."""
    pangenome, _genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()

    passed_lines, failed_lines = _gene_layer_lines(*_gate(session, pangenome))
    assert sorted(passed_lines) == sorted(GENE_LAYER_CHECKS), (passed_lines, failed_lines)
    assert failed_lines == []


def test_the_publish_gate_REFUSES_a_collection_genome_with_no_gene_spans_block(session, small_catalogue):
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()
    session.execute(GenomeGenes.__table__.delete().where(GenomeGenes.genome_id == genome["c"]))

    passed, failed = _gate(session, pangenome)
    assert not any(line.startswith(BLOCK_OR_NAMED) for line in passed)
    assert any(
        failure.startswith(
            f"{BLOCK_OR_NAMED} — 2 with a block, 0 named as carrying none, of 3: 1 of 3 collection genomes "
            "have no genome_genes row and no gene_table_omitted_reason"
        )
        and "['SAMEA0000003']" in failure
        for failure in failed
    ), failed
    # The other two still hold: the blocks that remain are well-formed, and C held no genes.
    assert "every gene_spans block is 16 bytes per gene" in passed
    assert BLOCK_TOTAL in passed


def test_the_publish_gate_REFUSES_a_block_total_that_is_not_the_pangenomes_gene_count_and_SAYS_both(
    session, small_catalogue
):
    pangenome, _genome = small_catalogue
    pangenome.gene_count = 6
    session.flush()

    passed, failed = _gate(session, pangenome)
    assert BLOCK_TOTAL not in passed
    assert (
        f"{BLOCK_TOTAL}: the blocks hold 5 genes and the 0 genome(s) named as carrying none 0 — 5 — and "
        "the pangenome row claims 6"
    ) in failed, failed


def _give_genome_b_a_gene_in_no_locus(session, pangenome, genome) -> None:
    """B's table gains a third gene the catalogue never clustered: `-1` in its block, placed still 2."""
    session.execute(
        GenomeGenes.__table__.update().where(GenomeGenes.genome_id == genome["b"]).values(
            gene_count=3, gene_spans=b"\x00" * 48
        )
    )
    session.execute(
        GenomeLocusAssignment.__table__.update()
        .where(
            GenomeLocusAssignment.pangenome_id == pangenome.pangenome_id,
            GenomeLocusAssignment.genome_id == genome["b"],
        )
        .values(
            gene_count=3, placed_gene_count=2,
            gene_locus_ordinals=encode_int32_vector([0, 1, -1]),
            gene_arrangement_ranks=encode_int32_vector([-1, -1, -1]),
        )
    )
    session.flush()


def test_the_publish_gate_sets_aside_a_gene_the_tables_hold_in_NO_locus_and_says_how_many(
    session, small_catalogue
):
    """⛔ The 71,435-genome cohort: 1,482 genes sit in gene tables and in no locus. The universe is still 5."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    _give_genome_b_a_gene_in_no_locus(session, pangenome, genome)

    passed, failed = _gate(session, pangenome)
    assert f"{BLOCK_TOTAL} — 1 gene(s) in 1 genome(s) are in no locus and are set aside" in passed, (passed, failed)
    assert not any(line.startswith(BLOCK_TOTAL) for line in failed)


def test_the_publish_gate_still_REFUSES_a_total_that_only_matches_because_of_a_gene_in_no_locus(
    session, small_catalogue
):
    """The tables now hold 6 genes. A gate that only asked for 'at least' — or for the raw total — would
    pass a universe of 6; it is 5, and the sixth gene is in no locus."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 6
    _give_genome_b_a_gene_in_no_locus(session, pangenome, genome)

    passed, failed = _gate(session, pangenome)
    assert not any(line.startswith(BLOCK_TOTAL) for line in passed)
    assert (
        f"{BLOCK_TOTAL} — 1 gene(s) in 1 genome(s) are in no locus and are set aside: the blocks hold 6 genes "
        "and the 0 genome(s) named as carrying none 0 — 6 — and the pangenome row claims 6"
    ) in failed, failed


def test_the_publish_gate_SAYS_when_the_block_total_was_not_compared(session, small_catalogue):
    """⚠ A zero `gene_count` is nothing to compare against — stated in the report, never a silent pass."""
    pangenome, _genome = small_catalogue
    pangenome.gene_count = 0
    session.flush()

    passed, failed = _gate(session, pangenome)
    assert any(
        line == "the gene_spans block total (5 genes) was NOT compared: pangenome.gene_count is 0" for line in passed
    ), passed
    assert BLOCK_TOTAL not in passed
    assert not any(line.startswith(BLOCK_TOTAL) for line in failed)


def test_the_publish_gate_REFUSES_a_block_that_is_not_16_bytes_per_gene(session, small_catalogue):
    """The table's own CHECK forbids such a row, so it is dropped for THIS transaction only — DDL is
    transactional in Postgres and the session fixture rolls back — to prove the gate re-asks it."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 6
    session.flush()
    # ⚠ Looked up by what it CHECKS, not by name: the metadata's naming convention prefixes it.
    (constraint_name,) = session.execute(text(
        "SELECT conname FROM pg_constraint WHERE conrelid = 'genome_genes'::regclass AND contype = 'c' "
        "AND pg_get_constraintdef(oid) LIKE '%octet_length%'"
    )).one()
    session.execute(text(f'ALTER TABLE genome_genes DROP CONSTRAINT "{constraint_name}"'))
    # C now claims one gene and carries fifteen bytes.
    session.execute(
        GenomeGenes.__table__.update().where(GenomeGenes.genome_id == genome["c"]).values(
            gene_count=1, gene_spans=b"\x00" * 15
        )
    )

    passed, failed = _gate(session, pangenome)
    assert "every gene_spans block is 16 bytes per gene" not in passed
    assert any(
        failure.startswith("every gene_spans block is 16 bytes per gene: 1 block(s) whose octet_length")
        for failure in failed
    ), failed
    # Every genome HAS a row and the counts add to 6 — only the length is wrong, and only that is said.
    assert GENE_LAYER_CHECKS[0] in passed
    assert BLOCK_TOTAL in passed


# ── the publish gate on a genome NAMED as holding no gene table ───────────────────────────────
NAMED_REASON = "the .nseq acceptance gate refused it: its stored proteins do not reproduce from the assembly"


def _name_b_without_gene_table(session, genome, *, coding_gene_count=2):
    """B as the genome stage leaves a named genome: no block, its `_meta` row count, the reason."""
    session.execute(GenomeGenes.__table__.delete().where(GenomeGenes.genome_id == genome["b"]))
    session.execute(
        Genome.__table__.update().where(Genome.genome_id == genome["b"]).values(
            gene_table_omitted_reason=NAMED_REASON, coding_gene_count=coding_gene_count
        )
    )


def test_the_publish_gate_PASSES_a_genome_NAMED_as_holding_no_gene_table_and_says_how_many(
    session, small_catalogue
):
    """⭐ The positive case: the blocks hold 3, B is named with 2, the pangenome claims 5 — every check passes,
    and the first one states its own counts."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()
    _name_b_without_gene_table(session, genome)

    passed, failed = _gate(session, pangenome)
    # ⚠ The fixture has no similarity layer, so other checks fail here; only the gene layer's are asked.
    assert not any(line.startswith((BLOCK_OR_NAMED, BLOCK_TOTAL, ASSIGNMENT_CHECKS[1])) for line in failed), failed
    assert f"{BLOCK_OR_NAMED} — 2 with a block, 1 named as carrying none, of 3" in passed
    assert BLOCK_TOTAL in passed
    assert ASSIGNMENT_CHECKS[1] in passed


def test_the_publish_gate_REFUSES_a_genome_with_NEITHER_a_block_nor_a_reason_by_name(session, small_catalogue):
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()
    session.execute(GenomeGenes.__table__.delete().where(GenomeGenes.genome_id == genome["b"]))

    passed, failed = _gate(session, pangenome)
    assert not any(line.startswith(BLOCK_OR_NAMED) for line in passed)
    assert any(
        line.startswith(BLOCK_OR_NAMED) and "no gene_table_omitted_reason" in line and "['SAMEA0000002']" in line
        for line in failed
    ), failed
    # ⚠ And the assignment row now has nothing to equal: counted as a disagreement, never skipped.
    assert any(line.startswith(f"{ASSIGNMENT_CHECKS[1]}: 1 assignment rows disagree") for line in failed), failed


def test_the_publish_gate_REFUSES_a_genome_with_BOTH_a_block_and_a_reason_by_name(session, small_catalogue):
    """⛔ A stale mark: the genome stage clears it whenever it writes a block, so both is a defect."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()
    # ⚠ With its true count: the fixture's genomes carry none, and a NULL would hide a total that
    # counted a stale-marked genome's genes twice — once in its block, once as "named".
    session.execute(
        Genome.__table__.update().where(Genome.genome_id == genome["a"]).values(
            gene_table_omitted_reason=NAMED_REASON, coding_gene_count=3
        )
    )

    passed, failed = _gate(session, pangenome)
    assert not any(line.startswith(BLOCK_OR_NAMED) for line in passed)
    assert any(
        line.startswith(BLOCK_OR_NAMED) and "the mark is stale" in line and "['SAMEA0000001']" in line
        for line in failed
    ), failed
    # Its genes are counted ONCE, in its block — the stale mark refuses, it does not double the total.
    assert BLOCK_TOTAL in passed, failed


def test_the_publish_gate_REFUSES_a_named_genome_whose_coding_gene_count_is_not_its_assignment_s(
    session, small_catalogue
):
    """B's assignment block holds 2 genes; a stored count of 3 disagrees with it and with the total."""
    pangenome, genome = small_catalogue
    pangenome.gene_count = 5
    session.flush()
    _name_b_without_gene_table(session, genome, coding_gene_count=3)

    passed, failed = _gate(session, pangenome)
    assert any(line.startswith(f"{ASSIGNMENT_CHECKS[1]}: 1 assignment rows disagree") for line in failed), failed
    assert any(
        line.startswith(f"{BLOCK_TOTAL}: the blocks hold 3 genes and the 1 genome(s) named as carrying none 3")
        for line in failed
    ), failed
    assert f"{BLOCK_OR_NAMED} — 2 with a block, 1 named as carrying none, of 3" in passed
