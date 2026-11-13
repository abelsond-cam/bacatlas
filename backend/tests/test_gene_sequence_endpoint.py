"""The Sequence tab's endpoint — the genome-grain blocks decoded for one locus, the bases sliced from the GFF.

⛔ **The defect this file exists for is a flank taken from the wrong end.** On a minus-strand gene —
about half of them — the upstream flank sits at HIGHER contig coordinates and is reverse-complemented
with the gene. Slicing `start - flank` unconditionally returns the *downstream* flank instead, and it
renders as a perfectly plausible 100 bases of DNA. Nothing on the page could contradict it, so the
checks below go back to the contig and read the bases themselves.

**Two sections, two databases.**

- *The probe database* (`BACATLAS_PROBE_DATABASE_URL`, the conftest fixtures): a hand-built genome
  with two copies at one locus, a minus-strand gene, a contig whose name is not its index + 1, an
  unplaced gene, a second genome the catalogue never placed, and three deliberately malformed
  genomes — served through a fake provider holding the bases in memory. Runs wherever the schema
  probe runs, which is every machine.
- *The loaded database* (`BACATLAS_DATABASE_URL` + `BACATLAS_ROOT_GFF`): the real catalogues and the
  real GFFs. Skips with the reason when either is absent. ⚠ Its candidate genes are chosen by
  DECODING THE BLOCKS — the same tables the service reads — never from a per-gene table, so the
  selection and the thing selected cannot drift apart.

⚠ Four fields left the row and the response on 2026-11-13 — `protein_length_aa`, `gene_symbol`,
`product`, `locus_tag` — and the assertions that read them went with them: nothing renders them
(David, 2026-10-26 / 2026-10-29). The protein length is checked against the SPAN instead.
"""

from __future__ import annotations

import os
import random
from collections import Counter

import pytest
from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.api.blueprint_species import RESPONSE_SHAPE_VERSION
from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.gene_blocks import GeneBlockError, encode_gene_spans, encode_int32_vector, pack_gene_flags
from bacatlas_backend.gff.gene_sequence_reader import reverse_complement
from bacatlas_backend.instruments.sql_cost_oracle import SqlCostOracle
from bacatlas_backend.models.genome import Genome, GenomeContig
from bacatlas_backend.models.genome_collection import GenomeCollection, GenomeCollectionMembership
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.sequence_source.contig_sequence_provider import (
    ContigSequenceProvider,
    GenomeSequenceRequest,
)
from bacatlas_backend.sequence_source.provider_factory import SEQUENCE_PROVIDER_KEY
from bacatlas_backend.services.gene_sequence_service import SequenceUnavailable, load_gene_sequences
from tests.conftest import PROBE_URL, make_locus
from tests.sequence_parity_comparison import GeneTableRow, decode_gene_table_rows, genome_gene_rows

# ══════════════════════════════════════════════════════════════════════════════════════════════
# The probe database — a genome built by hand, so every case is one the fixture put there on purpose
# ══════════════════════════════════════════════════════════════════════════════════════════════

HOST = "SAMEA0000001"  # the conftest seed genome: two copies at LOCUS, one at OTHER, one unplaced
OUTSIDE = "SAMEA0000002"  # in the collection and in `genome_genes`, never placed by this catalogue
UNOBSERVED = "SAMEA0000003"  # placed, with `strand_is_observed` False
MISMATCHED = "SAMEA0000004"  # `genome_genes` says 2 genes, the assignment says 3 — both CHECKs pass
UNLOADED = "SAMEA0000005"  # placed, but no `genome_genes` row at all
ENCODED_OTHERWISE = "SAMEA0000006"  # `span_encoding` 2, which no decoder knows
NAMED = "SAMEA0000007"  # placed, no `genome_genes` row, and NAMED as holding none — with the reason why
NAMED_REASON = "the .nseq acceptance gate refused it: its stored proteins do not reproduce from the assembly"

LOCUS = "2811"  # catalogue_ordinal 5 — HOST has genes 0 (+) and 1 (−) here
OTHER = "3126"  # catalogue_ordinal 7 — HOST has gene 2 here, on the renamed contig
LACKED = "9999"  # catalogue_ordinal 9 — in the catalogue, HOST has no gene there

PLUS = pack_gene_flags(is_minus_strand=False, is_five_prime_partial=False, gff_phase=0)
MINUS = pack_gene_flags(is_minus_strand=True, is_five_prime_partial=False, gff_phase=0)


def _bases(length: int, seed: int) -> str:
    """Pseudo-random bases, seeded — NOT a repeat or a palindrome, so a flank from the wrong end differs."""
    return "".join(random.Random(seed).choices("ACGT", k=length))


CONTIG_1 = _bases(1000, 1)  # HOST contig_index 0 → `contig00001`, which IS index + 1
CONTIG_7 = _bases(600, 7)  # HOST contig_index 1 → `contig00007`, which is NOT index + 1

#: HOST's four genes, 1-based inclusive, 90 bases each: `(start, end, contig_index, flags, ordinal)`.
HOST_GENES = (
    (101, 190, 0, PLUS, 5),
    (301, 390, 0, MINUS, 5),
    (51, 140, 1, PLUS, 7),
    (201, 290, 1, MINUS, -1),
)


class _HeldBases:
    """`ContigSequenceProvider`'s one method over bases held in memory — and it records what it was asked.

    `requests` is the witness that a "no gene here" answer never opened an assembly at all.
    """

    def __init__(self, genomes: dict[str, dict[str, str]]) -> None:
        self.genomes = genomes
        self.requests: list[GenomeSequenceRequest] = []

    def contig_sequences(self, request: GenomeSequenceRequest) -> dict[str, str]:
        self.requests.append(request)
        return self.genomes[request.sample_id]


def _gene_layer(genome_id: int, genes, *, span_encoding: int = 1) -> GenomeGenes:
    return GenomeGenes(
        genome_id=genome_id,
        gene_count=len(genes),
        span_encoding=span_encoding,
        gene_spans=encode_gene_spans(
            [g[0] for g in genes], [g[1] for g in genes], [g[2] for g in genes], [g[3] for g in genes]
        ),
    )


def _assignment(pangenome_id: int, genome_id: int, ordinals: list[int]) -> GenomeLocusAssignment:
    placed = [o for o in ordinals if o != -1]
    return GenomeLocusAssignment(
        pangenome_id=pangenome_id,
        genome_id=genome_id,
        gene_count=len(ordinals),
        placed_gene_count=len(placed),
        distinct_locus_count=len(set(placed)),
        gene_locus_ordinals=encode_int32_vector(ordinals, name="gene_locus_ordinals"),
        gene_arrangement_ranks=encode_int32_vector([-1] * len(ordinals), name="gene_arrangement_ranks"),
    )


def _contig(genome_id: int, sample_id: str, index: int, name: str, length: int) -> GenomeContig:
    return GenomeContig(
        genome_id=genome_id, seqid=f"{sample_id}.{name}", contig_index=index, contig_name=name, length_bases=length
    )


@pytest.fixture(scope="module")
def gene_layer(engine, _seed_ids):
    """The probe catalogue's gene layer, COMMITTED so the Flask client's own session can see it.

    ⚠ Committed rather than flushed, like `_seed_ids`: the HTTP tests go through the application's
    engine, which cannot see an uncommitted row of this module's session.
    """
    with Session(engine) as session:
        species = session.get(PathogenSpecies, _seed_ids["species"])
        collection = session.get(GenomeCollection, _seed_ids["collection"])
        pangenome = session.get(Pangenome, _seed_ids["pangenome"])
        host = session.get(Genome, _seed_ids["genome"])
        seeded = {"species": species, "collection": collection, "pangenome": pangenome, "genome": host}
        pangenome_id = pangenome.pangenome_id
        # So `/api/v1/species/probe/…` resolves: the species serves this pangenome.
        species.default_pangenome_id = pangenome_id

        session.add_all(
            [
                make_locus(seeded, ordinal=5, label=LOCUS),
                make_locus(seeded, ordinal=7, label=OTHER),
                make_locus(seeded, ordinal=9, label=LACKED, present_genome_ids=[], member_genome_count=0),
            ]
        )
        session.add_all(
            [
                _contig(host.genome_id, HOST, 0, "contig00001", len(CONTIG_1)),
                _contig(host.genome_id, HOST, 1, "contig00007", len(CONTIG_7)),
                _gene_layer(host.genome_id, HOST_GENES),
                _assignment(pangenome_id, host.genome_id, [g[4] for g in HOST_GENES]),
            ]
        )

        others = {}
        for _ordinal, (sample_id, strand_is_observed) in enumerate(
            [(OUTSIDE, True), (UNOBSERVED, False), (MISMATCHED, True), (UNLOADED, True), (ENCODED_OTHERWISE, True)],
            start=1,
        ):
            genome = Genome(
                pathogen_species_id=species.pathogen_species_id,
                sample_id=sample_id,
                sample_id_kind=host.sample_id_kind,
                strand_is_observed=strand_is_observed,
            )
            session.add(genome)
            session.flush()
            others[sample_id] = genome.genome_id
            session.add(_contig(genome.genome_id, sample_id, 0, "contig00001", 400))
        memberships = [(host.genome_id, 0)] + [(genome_id, n) for n, genome_id in enumerate(others.values(), 1)]
        session.add_all(
            [
                GenomeCollectionMembership(
                    genome_collection_id=collection.genome_collection_id,
                    genome_id=genome_id,
                    collection_genome_ordinal=ordinal,
                )
                for genome_id, ordinal in memberships
            ]
        )
        one_gene = ((11, 100, 0, PLUS, 7),)
        session.add_all(
            [
                # OUTSIDE: a gene table, no assignment row — this catalogue never placed it.
                _gene_layer(others[OUTSIDE], one_gene),
                # UNOBSERVED: placed at OTHER, strands defaulted rather than read.
                _gene_layer(others[UNOBSERVED], one_gene),
                _assignment(pangenome_id, others[UNOBSERVED], [7]),
                # MISMATCHED: two genes in one table, three in the other. Each block passes its own CHECK.
                _gene_layer(others[MISMATCHED], ((11, 100, 0, PLUS, 5), (111, 200, 0, PLUS, 5))),
                _assignment(pangenome_id, others[MISMATCHED], [5, 5, 5]),
                # UNLOADED: placed, and no gene table at all.
                _assignment(pangenome_id, others[UNLOADED], [5]),
                # ENCODED_OTHERWISE: a layout this decoder does not know.
                _gene_layer(others[ENCODED_OTHERWISE], ((11, 100, 0, PLUS, 5),), span_encoding=2),
                _assignment(pangenome_id, others[ENCODED_OTHERWISE], [5]),
            ]
        )
        # NAMED: what the genome stage leaves for a genome nuna's `.nseq` gate refused — a `genome` row
        # with the reason and its `_meta` gene count, no block, no contigs — placed at LOCUS only.
        named = Genome(
            pathogen_species_id=species.pathogen_species_id,
            sample_id=NAMED,
            sample_id_kind=host.sample_id_kind,
            strand_is_observed=False,
            coding_gene_count=2,
            gene_table_omitted_reason=NAMED_REASON,
        )
        session.add(named)
        session.flush()
        others[NAMED] = named.genome_id
        session.add_all(
            [
                GenomeCollectionMembership(
                    genome_collection_id=collection.genome_collection_id,
                    genome_id=named.genome_id,
                    collection_genome_ordinal=len(memberships),
                ),
                _assignment(pangenome_id, named.genome_id, [5, -1]),
            ]
        )
        session.commit()
        return {"pangenome_id": pangenome_id, "host_id": host.genome_id, **others}


@pytest.fixture()
def held():
    return _HeldBases(
        {
            HOST: {f"{HOST}.contig00001": CONTIG_1, f"{HOST}.contig00007": CONTIG_7},
            **{
                sample: {f"{sample}.contig00001": _bases(400, 3)}
                for sample in (OUTSIDE, UNOBSERVED, MISMATCHED, UNLOADED, ENCODED_OTHERWISE)
            },
        }
    )


@pytest.fixture(scope="module")
def probe_application(gene_layer):
    application = create_application(Configuration(database_url=PROBE_URL))
    # ⛔ As in production: an unhandled error is a 500 here, not an exception in the test runner.
    application.config["PROPAGATE_EXCEPTIONS"] = False
    return application


@pytest.fixture()
def probe_client(probe_application, held):
    probe_application.extensions[SEQUENCE_PROVIDER_KEY] = held
    return probe_application.test_client()


def _load(session, gene_layer, held, *, sample_id: str, label: str):
    return load_gene_sequences(
        session, pangenome_id=gene_layer["pangenome_id"], node_label=label, sample_id=sample_id, sequence_provider=held
    )


def _url(sample_id: str, label: str) -> str:
    return f"/api/v1/species/probe/genomes/{sample_id}/loci/{label}/sequence"


def test_two_copies_at_one_locus_return_copy_count_2_in_flat_index_order(session, gene_layer, held, probe_client):
    """⛔ `positions_of`, not a first hit: ρ > 1 is one genome at one locus twice, and both are served."""
    rows = _load(session, gene_layer, held, sample_id=HOST, label=LOCUS)
    assert [row.flat_index for row in rows] == [0, 1]
    assert [row.copy_ordinal for row in rows] == [1, 2]
    assert {row.copy_count for row in rows} == {2}
    assert [row.strand for row in rows] == ["+", "-"]
    assert [row.length_nt for row in rows] == [90, 90]
    assert rows[0].sequence.coding_sequence == CONTIG_1[100:190]
    assert rows[1].sequence.coding_sequence == reverse_complement(CONTIG_1[300:390])
    # ⚠ One assembly, however many copies.
    assert [request.sample_id for request in held.requests] == [HOST]

    genes = probe_client.get(_url(HOST, LOCUS)).get_json()["genes"]
    assert [(gene["copy_ordinal"], gene["copy_count"], gene["flat_index"]) for gene in genes] == [(1, 2, 0), (2, 2, 1)]


def test_a_MINUS_strand_gene_takes_its_flanks_from_the_HIGHER_end_and_a_PLUS_gene_from_the_LOWER(
    session, gene_layer, held
):
    """⛔⛔ The bug of record, against bases the fixture wrote — so the expected flank is not a guess."""
    plus, minus = _load(session, gene_layer, held, sample_id=HOST, label=LOCUS)
    assert minus.strand == "-"
    assert minus.sequence.upstream_flank_sequence == reverse_complement(CONTIG_1[390:490])
    assert minus.sequence.downstream_flank_sequence == reverse_complement(CONTIG_1[200:300])
    assert tuple(minus.sequence.upstream_flank_span) == (391, 490)
    assert tuple(minus.sequence.downstream_flank_span) == (201, 300)
    # ⚠ And they are NOT each other: on a palindromic region this test would pass either way.
    assert minus.sequence.upstream_flank_sequence != minus.sequence.downstream_flank_sequence

    assert plus.strand == "+"
    assert plus.sequence.upstream_flank_sequence == CONTIG_1[0:100]
    assert plus.sequence.downstream_flank_sequence == CONTIG_1[190:290]
    assert tuple(plus.sequence.upstream_flank_span) == (1, 100)


def test_a_genome_the_catalogue_never_placed_is_an_EMPTY_LIST_and_not_an_error(session, gene_layer, held, probe_client):
    """⛔ No `genome_locus_assignment` row → `[]`, an ANSWER. And no assembly is opened to give it."""
    assert _load(session, gene_layer, held, sample_id=OUTSIDE, label=LOCUS) == []
    assert held.requests == []
    response = probe_client.get(_url(OUTSIDE, LOCUS))
    assert response.status_code == 200
    assert response.get_json()["genes"] == []


def test_a_locus_the_catalogue_lacks_or_the_genome_lacks_is_an_EMPTY_LIST(session, gene_layer, held, probe_client):
    """Two absences, both `[]` from the service; the endpoint NAMES the first as a 404 before asking."""
    assert _load(session, gene_layer, held, sample_id=HOST, label="not-a-locus") == []
    assert _load(session, gene_layer, held, sample_id=HOST, label=LACKED) == []
    assert held.requests == []

    lacked = probe_client.get(_url(HOST, LACKED))
    assert lacked.status_code == 200 and lacked.get_json()["genes"] == []
    missing = probe_client.get(_url(HOST, "not-a-locus"))
    assert missing.status_code == 404 and "not-a-locus" in missing.get_json()["detail"]


def test_the_contig_name_is_the_seqid_suffix_after_the_first_dot_and_not_index_plus_one(session, gene_layer, held):
    """⛔ `contig_index` 1 is `contig00007` here. `contig{index + 1}` would print `contig00002`."""
    [row] = _load(session, gene_layer, held, sample_id=HOST, label=OTHER)
    assert row.flat_index == 2
    assert row.seqid == f"{HOST}.contig00007"
    assert row.contig_name == "contig00007" == row.seqid.split(".", 1)[1]
    assert row.contig_name != "contig00002"
    assert row.sequence.coding_sequence == CONTIG_7[50:140]


def test_strand_is_observed_is_the_GENOME_S_flag(session, gene_layer, held):
    """⚠ A forced `+` must never read as an observation — the flag travels from `genome`, per row."""
    [unobserved] = _load(session, gene_layer, held, sample_id=UNOBSERVED, label=OTHER)
    assert unobserved.strand_is_observed is False
    assert all(row.strand_is_observed for row in _load(session, gene_layer, held, sample_id=HOST, label=LOCUS))


def test_the_response_has_NONE_of_the_four_retired_fields_and_the_shape_version_is_5(gene_layer, probe_client):
    """⛔ A removal needs the bump as much as an addition: a warm client would otherwise replay a
    day-old body carrying fields the server no longer has, and nothing anywhere would error."""
    response = probe_client.get(_url(HOST, LOCUS))
    assert response.status_code == 200
    genes = response.get_json()["genes"]
    assert len(genes) == 2
    retired = {"gene_symbol", "product", "locus_tag", "protein_length_aa"}
    for gene in genes:
        assert retired.isdisjoint(gene), sorted(retired & set(gene))
    assert set(genes[0]) == {
        "copy_ordinal",
        "copy_count",
        "flat_index",
        "contig_name",
        "seqid",
        "start_position",
        "end_position",
        "strand",
        "strand_is_observed",
        "gff_phase",
        "is_five_prime_partial",
        "length_nt",
        "coding_sequence",
        "protein_sequence",
        "upstream_flank_sequence",
        "downstream_flank_sequence",
        "upstream_flank_is_truncated_by_contig_end",
        "downstream_flank_is_truncated_by_contig_end",
        "upstream_flank_span",
        "downstream_flank_span",
        "gc_percent",
    }
    assert RESPONSE_SHAPE_VERSION == 5
    assert response.headers["ETag"] == f'"{gene_layer["pangenome_id"]}-v5"'


@pytest.mark.parametrize(
    "sample_id,reason",
    [
        (MISMATCHED, "genome_locus_assignment says 3 genes and genome_genes says 2"),
        (UNLOADED, "no genome_genes row"),
        (ENCODED_OTHERWISE, "span_encoding is 2"),
    ],
)
def test_a_malformed_gene_layer_is_REFUSED_by_name_and_NEVER_served_as_no_gene_here(
    session, gene_layer, held, probe_client, sample_id, reason
):
    """⛔ Three defects the CHECK constraints cannot see. Each is a `GeneBlockError` naming the genome
    and the fault, and at the route a 500 — never a 200 with `genes: []`, which a reader would take
    as "this genome has nothing here"."""
    with pytest.raises(GeneBlockError, match=reason):
        _load(session, gene_layer, held, sample_id=sample_id, label=LOCUS)
    assert held.requests == []
    response = probe_client.get(_url(sample_id, LOCUS))
    assert response.status_code == 500
    assert response.get_json() == {"error": "server_error", "detail": "the server could not answer this request"}


def test_a_NAMED_genome_with_a_gene_here_is_a_NAMED_non_retryable_503_and_never_a_500(
    session, gene_layer, held, probe_client
):
    """⭐ A modelled genome with no gene table: the gene exists (its assignment says so) and its bases
    cannot be shown, for the stored reason. The existing "sequence unavailable" path, without
    `Retry-After` — no retry will produce a table — and without opening an assembly."""
    with pytest.raises(SequenceUnavailable) as raised:
        _load(session, gene_layer, held, sample_id=NAMED, label=LOCUS)
    assert raised.value.is_retryable is False
    detail = (
        f"{NAMED} has a gene at this locus, but BacAtlas holds no gene table for this genome — "
        f"{NAMED_REASON} — so its coordinates and bases cannot be shown"
    )
    assert str(raised.value) == detail

    response = probe_client.get(_url(NAMED, LOCUS))
    assert response.status_code == 503
    assert "Retry-After" not in response.headers
    assert response.get_json() == {"error": "sequence_unavailable", "detail": detail}
    assert held.requests == []


def test_a_NAMED_genome_with_NO_gene_here_is_still_the_ANSWER_no_gene_here(session, gene_layer, held, probe_client):
    """Its assignment block is real, so *"this genome has no gene at this locus"* is still true and still said."""
    assert _load(session, gene_layer, held, sample_id=NAMED, label=OTHER) == []
    response = probe_client.get(_url(NAMED, OTHER))
    assert response.status_code == 200
    assert response.get_json()["genes"] == []
    assert held.requests == []


def test_the_sequence_view_is_ONE_statement_whatever_it_answers(engine, gene_layer, held):
    """⚠ Two copies, no gene here, and a locus the catalogue lacks: one statement each.

    The contig map rides on the same statement as the blocks, so the count does not grow with the
    copies, the contigs, or the catalogue.
    """
    for sample_id, label, copies in [(HOST, LOCUS, 2), (HOST, LACKED, 0), (OUTSIDE, LOCUS, 0), (HOST, "none", 0)]:
        with Session(engine) as session, SqlCostOracle(engine) as report:
            rows = _load(session, gene_layer, held, sample_id=sample_id, label=label)
        assert len(rows) == copies
        report.assert_at_most(1, what=f"the sequence view for {sample_id} at {label}")
    # ⚠ The named state reads its reason off the same statement — no second one to learn why.
    with Session(engine) as session, SqlCostOracle(engine) as report, pytest.raises(SequenceUnavailable):
        _load(session, gene_layer, held, sample_id=NAMED, label=LOCUS)
    report.assert_at_most(1, what=f"the sequence view for the named genome {NAMED}")


def test_the_sequence_ROUTE_is_the_view_plus_the_three_named_lookups(probe_application, probe_client):
    """The species, the genome-in-collection and the locus — each a named 404 when it misses — then
    the view's one statement. ⚠ Held to the same number in `test_route_cost.SEQUENCE_TAB`."""
    for sample_id, status in ((HOST, 200), (NAMED, 503)):
        with SqlCostOracle(probe_application.extensions["bacatlas_database"].engine) as report:
            response = probe_client.get(_url(sample_id, LOCUS))
        assert response.status_code == status
        report.assert_at_most(4, what=f"the Sequence tab, at the route, for {sample_id}")


def test_an_unreadable_assembly_RAISES_rather_than_returning_no_genes(session, gene_layer, tmp_path):
    """⛔ The distinction the whole endpoint turns on: *we could not read it* is not *there is none*.

    ⚠ The probe genome has no `bakrep_dataset_id`, so the held-only chain's named reason is "no
    recorded annotation directory" rather than "not on this server" — the message differs, the
    verdict must not: a raise, never `[]`.
    """
    unreadable = ContigSequenceProvider(gff_root=tmp_path, atb_store=None)
    with pytest.raises(SequenceUnavailable, match="cannot be located"):
        _load(session, gene_layer, unreadable, sample_id=HOST, label=LOCUS)


# ── the parity suite's database side, which decodes the same blocks ───────────────────────────
def test_decode_gene_table_rows_reads_the_blocks_the_way_the_page_reads_its_nseq():
    """⭐ No database: the pure helper, on blocks the test encoded, with every column pinned."""
    spans = encode_gene_spans([101, 301], [190, 390], [0, 1], [PLUS, MINUS])
    contigs = {0: ("S.contig00001", "contig00001", 1000), 1: ("S.contig00007", "contig00007", 600)}
    rows = decode_gene_table_rows(
        gene_count=2,
        gene_spans=spans,
        gene_locus_ordinals=encode_int32_vector([5, -1]),
        contigs_by_index=contigs,
        label_by_ordinal={5: "2811"},
    )
    assert rows == [
        GeneTableRow(0, 101, 190, 0, "+", 0, False, "S.contig00001", "contig00001", 1000, "2811", 5),
        GeneTableRow(1, 301, 390, 1, "-", 0, False, "S.contig00007", "contig00007", 600, None, None),
    ]
    # No assignment row at all → every gene unplaced, and nothing is looked up.
    unplaced = decode_gene_table_rows(
        gene_count=2, gene_spans=spans, gene_locus_ordinals=None, contigs_by_index=contigs, label_by_ordinal={}
    )
    assert [row.node_label for row in unplaced] == [None, None]
    # ⛔ A contig the map lacks and an ordinal the catalogue lacks are refused, never a row with a None.
    with pytest.raises(LookupError, match="contig_index 1"):
        decode_gene_table_rows(
            gene_count=2,
            gene_spans=spans,
            gene_locus_ordinals=None,
            contigs_by_index={0: contigs[0]},
            label_by_ordinal={},
        )
    with pytest.raises(LookupError, match="catalogue_ordinal 5"):
        decode_gene_table_rows(
            gene_count=2,
            gene_spans=spans,
            gene_locus_ordinals=encode_int32_vector([5, -1]),
            contigs_by_index=contigs,
            label_by_ordinal={},
        )


def test_genome_gene_rows_decodes_the_probe_genome_exactly_as_the_fixture_wrote_it(session, gene_layer):
    """The database-backed helper, against the rows this module put there — labels resolved, unplaced kept."""
    rows = genome_gene_rows(session, genome_id=gene_layer["host_id"], pangenome_id=gene_layer["pangenome_id"])
    assert [(r.flat_index, r.start_position, r.end_position, r.contig_index, r.strand) for r in rows] == [
        (0, 101, 190, 0, "+"),
        (1, 301, 390, 0, "-"),
        (2, 51, 140, 1, "+"),
        (3, 201, 290, 1, "-"),
    ]
    assert [r.node_label for r in rows] == [LOCUS, LOCUS, OTHER, None]
    assert [r.catalogue_ordinal for r in rows] == [5, 5, 7, None]
    assert [r.contig_name for r in rows] == ["contig00001", "contig00001", "contig00007", "contig00007"]
    assert [r.length_bases for r in rows] == [1000, 1000, 600, 600]
    # A genome the catalogue never placed still has its gene table, all of it unassigned.
    [outside] = genome_gene_rows(session, genome_id=gene_layer[OUTSIDE], pangenome_id=gene_layer["pangenome_id"])
    assert outside.node_label is None and outside.start_position == 11
    with pytest.raises(LookupError, match="no genome_genes row"):
        genome_gene_rows(session, genome_id=gene_layer[UNLOADED], pangenome_id=gene_layer["pangenome_id"])


def test_an_empty_flank_has_a_NULL_span_rather_than_a_zero_length_one():
    """⚠ `[1, 0]` is not a span; it is arithmetic leaking. A flank that does not exist says so."""
    from bacatlas_backend.gff.gene_sequence_reader import read_gene_sequence

    view = read_gene_sequence("ATGAAATAG", start_position=1, end_position=9, strand="+", flank_length=100)
    assert view.upstream_flank_span is None
    assert view.upstream_flank_sequence == ""
    assert view.upstream_flank_is_truncated_by_contig_end is True


# ══════════════════════════════════════════════════════════════════════════════════════════════
# The loaded database — the real catalogues and the real GFFs; skips, with the reason, without them
# ══════════════════════════════════════════════════════════════════════════════════════════════

PANGENOME = 1


@pytest.fixture(scope="module")
def loaded_application():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set")
    root = os.environ.get("BACATLAS_ROOT_GFF")
    if not root:
        pytest.skip("BACATLAS_ROOT_GFF is not set — the sequence endpoint reads the original GFFs")
    return create_application(Configuration.from_environment())


@pytest.fixture(scope="module")
def loaded_client(loaded_application):
    return loaded_application.test_client()


@pytest.fixture(scope="module")
def loaded_engine(loaded_application):
    return loaded_application.extensions["bacatlas_database"].engine


@pytest.fixture(scope="module")
def gff_root(loaded_application):
    return loaded_application.config["BACATLAS"].artifact_roots["gff"]


@pytest.fixture()
def provider(gff_root):
    """The HELD-only chain — no ATB store, so these tests cannot reach the network.

    ⛔ Function-scoped on purpose: a module-scoped provider would carry its parsed-genome cache
    between tests, and the "one statement and one file" assertion below measures a cold read.
    """
    return ContigSequenceProvider(gff_root=gff_root, atb_store=None)


def _catalogue_genomes(session) -> list[tuple[str, int]]:
    """Every genome the published pangenome placed — an assignment row is what membership IS now."""
    return [
        (str(sample_id), int(genome_id))
        for sample_id, genome_id in session.execute(
            select(Genome.sample_id, Genome.genome_id)
            .join(GenomeLocusAssignment, GenomeLocusAssignment.genome_id == Genome.genome_id)
            .where(GenomeLocusAssignment.pangenome_id == PANGENOME)
            .order_by(Genome.genome_id)
        )
    ]


def _placed_rows(session, genome_id: int) -> list[GeneTableRow]:
    return [
        row
        for row in genome_gene_rows(session, genome_id=genome_id, pangenome_id=PANGENOME)
        if row.node_label is not None
    ]


def _first_gene(session, predicate) -> tuple[str, str, GeneTableRow]:
    """The first placed gene of the catalogue satisfying `predicate`, decoding one genome at a time."""
    for sample_id, genome_id in _catalogue_genomes(session):
        for row in _placed_rows(session, genome_id):
            if predicate(row):
                return sample_id, row.node_label, row
    raise AssertionError("no placed gene in the catalogue satisfies the selection")


def _a_gene(session, *, strand: str):
    """One (sample_id, node_label) whose gene sits on `strand`, away from either contig end."""
    sample_id, label, _ = _first_gene(
        session,
        lambda row: (
            row.strand == strand and row.start_position > 500 and row.end_position - row.start_position + 1 > 900
        ),
    )
    return sample_id, label


def _dataset_id(engine, sample_id: str) -> str:
    with Session(engine) as session:
        return session.execute(select(Genome.bakrep_dataset_id).where(Genome.sample_id == sample_id)).scalar_one()


# ── the flank orientation, which is the whole point ────────────────────────────────────────────
def test_a_MINUS_strand_gene_takes_its_upstream_flank_from_the_HIGHER_coordinates(loaded_engine, provider):
    """⛔⛔ **The bug of record, checked against the contig itself.**

    `load_meta_flanks` orients on strand and `app.js:4376-4386` states the convention. A reader
    following a guide designed against it would otherwise be pointed at the opposite end of the gene.
    """
    with Session(loaded_engine) as session:
        sample_id, label = _a_gene(session, strand="-")
        rows = load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    assert rows, f"{sample_id} has no gene at {label}"
    gene = rows[0]
    assert gene.strand == "-"

    contig = provider.contig_sequences(
        GenomeSequenceRequest(sample_id=sample_id, bakrep_dataset_id=_dataset_id(loaded_engine, sample_id))
    )[gene.seqid]
    higher = contig[gene.end_position : gene.end_position + 100]
    lower = contig[max(0, gene.start_position - 101) : gene.start_position - 1]

    assert gene.sequence.upstream_flank_sequence == reverse_complement(higher)
    assert gene.sequence.downstream_flank_sequence == reverse_complement(lower)
    # ⚠ And they are NOT each other: on a palindromic region this test would pass either way.
    assert gene.sequence.upstream_flank_sequence != gene.sequence.downstream_flank_sequence


def test_a_PLUS_strand_gene_takes_its_upstream_flank_from_the_LOWER_coordinates(loaded_engine, provider):
    with Session(loaded_engine) as session:
        sample_id, label = _a_gene(session, strand="+")
        rows = load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    gene = rows[0]
    contig = provider.contig_sequences(
        GenomeSequenceRequest(sample_id=sample_id, bakrep_dataset_id=_dataset_id(loaded_engine, sample_id))
    )[gene.seqid]
    assert gene.sequence.upstream_flank_sequence == contig[gene.start_position - 101 : gene.start_position - 1]
    assert gene.sequence.downstream_flank_sequence == contig[gene.end_position : gene.end_position + 100]


# ── the sliced values against the span ─────────────────────────────────────────────────────────
def test_the_sliced_sequence_LENGTHS_agree_with_the_span(loaded_engine, provider):
    """⭐ Two computations of one thing: the span from the block, the bases sliced from the file.

    ⚠ **This test compared `protein_length_aa` against the sliced protein until 2026-11-13.** That
    column went with the per-gene table — nothing rendered it (David, 2026-10-26 / 2026-10-29) — and
    the comparison it carried survives without it: the span INCLUDES the stop codon, so a phase-0 CDS
    is `length_nt / 3 - 1` residues, which is the derivation the design records (`storing_the_gene_layer.md`
    §6 item 4) and the one a non-zero phase would break. `stored_gc_percent` went the same way on
    2026-10-26; T6 compares GC from the page's own decoded bases instead.
    """
    checked = derived = 0
    with Session(loaded_engine) as session:
        for sample_id, genome_id in _catalogue_genomes(session)[:3]:
            for row in _placed_rows(session, genome_id)[:20]:
                for gene in load_gene_sequences(
                    session,
                    pangenome_id=PANGENOME,
                    node_label=row.node_label,
                    sample_id=sample_id,
                    sequence_provider=provider,
                ):
                    assert len(gene.sequence.coding_sequence) == gene.length_nt, (sample_id, gene.flat_index)
                    checked += 1
                    if gene.gff_phase == 0 and gene.length_nt % 3 == 0:
                        assert len(gene.sequence.protein_sequence) == gene.length_nt // 3 - 1, (
                            sample_id,
                            gene.flat_index,
                        )
                        derived += 1
    # ⛔ Coverage before the verdict: a selection that matched nothing would report 0 green assertions.
    assert checked >= 40, f"only {checked} genes were examined"
    assert derived >= 30, f"only {derived} phase-0 proteins were derived"


# ── the three different answers ────────────────────────────────────────────────────────────────
def test_a_genome_with_NO_gene_here_is_an_ANSWER_and_not_a_failure(loaded_client, loaded_engine):
    """⛔ `genes: []` with a 200. A 404 would read as "that locus does not exist"."""
    with Session(loaded_engine) as session:
        label, present = session.execute(
            select(Locus.node_label, Locus.present_genome_ids)
            .where(Locus.pangenome_id == PANGENOME, Locus.member_genome_count == 1)
            .limit(1)
        ).one()
        carriers = set(present)
        absent = next(sample_id for sample_id, genome_id in _catalogue_genomes(session) if genome_id not in carriers)
    response = loaded_client.get(f"/api/v1/species/ecoli/genomes/{absent}/loci/{label}/sequence")
    assert response.status_code == 200
    assert response.get_json()["genes"] == []


def test_an_unknown_genome_and_an_unknown_locus_are_NAMED_404s(loaded_client):
    missing_genome = loaded_client.get("/api/v1/species/ecoli/genomes/SAMEA000000/loci/2811/sequence")
    assert missing_genome.status_code == 404
    assert "SAMEA000000" in missing_genome.get_json()["detail"]

    missing_locus = loaded_client.get("/api/v1/species/ecoli/genomes/SAMEA103923484/loci/not-a-locus/sequence")
    assert missing_locus.status_code == 404
    assert "not-a-locus" in missing_locus.get_json()["detail"]


def test_an_unreadable_annotation_file_RAISES_rather_than_returning_no_genes(loaded_engine, tmp_path):
    """⛔ The distinction the whole endpoint turns on: *we could not read it* is not *there is none*."""
    with Session(loaded_engine) as session:
        sample_id, label = _a_gene(session, strand="+")
        with pytest.raises(SequenceUnavailable, match="not on this server"):
            load_gene_sequences(
                session,
                pangenome_id=PANGENOME,
                node_label=label,
                sample_id=sample_id,
                sequence_provider=ContigSequenceProvider(gff_root=tmp_path, atb_store=None),
            )


# ── the things that are easy to get plausibly wrong ────────────────────────────────────────────
def test_the_contig_is_named_by_its_NAME_and_not_by_its_index(loaded_client, loaded_engine):
    """⛔ `contig_index` enumerates contigs that HAVE a CDS, so `index + 1` is a different contig.

    ⚠ **And it coincides most of the time, which is exactly why this is dangerous.** Measured on the
    loaded catalogues: the derived name is wrong on **7,696 of 26,878 contigs (28.6 %)** and right on
    the rest — so a page printing `contig{index + 1}` looks correct on seven contigs in ten and names
    a real contig the gene is not on for the other three. The test therefore selects from the
    DIVERGING set rather than taking whatever comes first, which is how this assertion passed
    vacuously the first time it was written.
    """
    with Session(loaded_engine) as session:
        sample_id, label, row = _first_gene(session, lambda r: r.contig_name != f"contig{r.contig_index + 1:05d}")
    payload = loaded_client.get(f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence").get_json()
    gene = next(gene for gene in payload["genes"] if gene["flat_index"] == row.flat_index)
    assert gene["contig_name"] == row.contig_name
    assert gene["contig_name"] != f"contig{row.contig_index + 1:05d}"
    assert gene["seqid"].endswith(gene["contig_name"])


def test_rho_above_one_returns_EVERY_copy_with_its_ordinal(loaded_client, loaded_engine):
    """⛔ `COUNT(DISTINCT)` territory: one genome really can hold two genes at one locus."""
    found = None
    with Session(loaded_engine) as session:
        for sample_id, genome_id in _catalogue_genomes(session):
            copies = Counter(row.node_label for row in _placed_rows(session, genome_id))
            label, count = max(copies.items(), key=lambda item: item[1])
            if count > 1:
                found = (sample_id, label, count)
                break
    if found is None:
        pytest.skip("no genome holds two genes at one locus in this catalogue")
    sample_id, label, copies = found
    genes = loaded_client.get(f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence").get_json()["genes"]
    assert len(genes) == copies
    assert [gene["copy_ordinal"] for gene in genes] == list(range(1, copies + 1))
    assert {gene["copy_count"] for gene in genes} == {copies}
    # ⚠ And they are different genes, not one row repeated — in flat_index order.
    flat_indexes = [gene["flat_index"] for gene in genes]
    assert len(set(flat_indexes)) == copies and flat_indexes == sorted(flat_indexes)


def test_a_flank_that_RUNS_OFF_the_contig_says_so(loaded_client, loaded_engine):
    """⚠ Draft assemblies: mean contig 14.4 kb, ~14 genes, so this is common rather than exotic.

    A short flank that does not say it is short reads as a complete one.
    """
    with Session(loaded_engine) as session:
        sample_id, label, row = _first_gene(session, lambda r: r.start_position < 50)
    genes = loaded_client.get(f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence").get_json()["genes"]
    gene = next(gene for gene in genes if gene["flat_index"] == row.flat_index)
    truncated = (
        gene["upstream_flank_is_truncated_by_contig_end"]
        if gene["strand"] == "+"
        else gene["downstream_flank_is_truncated_by_contig_end"]
    )
    assert truncated is True
    assert len(gene["upstream_flank_sequence"]) < 100 or len(gene["downstream_flank_sequence"]) < 100


def test_the_sequence_view_is_ONE_statement_and_ONE_file(loaded_engine, provider):
    """⚠ However many copies the genome holds. The file is opened once and cached per process."""
    with Session(loaded_engine) as session:
        sample_id, label = _a_gene(session, strand="+")
    with Session(loaded_engine) as session, SqlCostOracle(loaded_engine) as report:
        load_gene_sequences(
            session, pangenome_id=PANGENOME, node_label=label, sample_id=sample_id, sequence_provider=provider
        )
    report.assert_at_most(1, what="one gene sequence view")


def test_each_flank_NAMES_the_contig_coordinates_it_came_from(loaded_client, loaded_engine):
    """⭐ So the page can say *"contig 1,388–1,487, reverse-complemented"* without re-deriving it.

    ⛔ A client working out which end a flank came from would have to re-implement the strand rule,
    which is the exact thing the reader exists to hold in one place — and getting it wrong there
    mislabels a correct sequence, which is worse than showing the wrong one.
    """
    with Session(loaded_engine) as session:
        sample_id, label = _a_gene(session, strand="-")
    genes = loaded_client.get(f"/api/v1/species/ecoli/genomes/{sample_id}/loci/{label}/sequence").get_json()["genes"]
    gene = next(gene for gene in genes if gene["strand"] == "-")

    upstream = gene["upstream_flank_span"]
    downstream = gene["downstream_flank_span"]
    # ⛔ On a MINUS gene the upstream flank is at the HIGHER coordinates.
    assert upstream[0] > gene["end_position"]
    assert downstream[1] < gene["start_position"]
    # The spans describe exactly the strings that were returned.
    assert upstream[1] - upstream[0] + 1 == len(gene["upstream_flank_sequence"])
    assert downstream[1] - downstream[0] + 1 == len(gene["downstream_flank_sequence"])
