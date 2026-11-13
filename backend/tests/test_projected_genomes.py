"""Genomes placed on a pangenome they were never part of.

Every test here pins something that would otherwise produce a page that looks entirely ordinary:
a modelled genome placed on its own model; a locus label from another run resolving to a different
locus; a projected genome's matched arrangement missing from the response because it sat past the
display cap; an agreement count served without the denominator that bounds it; and the three
absences — no vector, no window, no matching arrangement — collapsed into one.

⛔ And the one that guards everything else: **loading a projection changes no count the page already
quotes.** The genome list, an unanchored locus and the shell must be byte-identical before and after.
"""

from __future__ import annotations

import json

import pytest
from sqlalchemy import text

from bacatlas_backend.gene_blocks import encode_gene_spans, pack_gene_flags
from bacatlas_backend.ingest.ingest_projected_genomes import ProjectionRefused, load_projection
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.projected_genome import ProjectedGenePlacement, ProjectedGenome
from bacatlas_backend.services.locus_detail_service import load_locus_detail
from bacatlas_backend.services.projected_genome_service import (
    list_projected_genomes,
    resolve_projected_genome,
)
from tests.conftest import TEN_ZEROS, make_locus

OUTSIDER = "SAMN09999999"
RULE = "the nearest modelled gene decides the locus; the first n check it"
WINDOW = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]


@pytest.fixture()
def catalogue(session, seeded):
    """One locus with two arrangements, and an outsider genome with three genes.

    The second arrangement is rank 9 — **past the display cap of 8** — because that is the case the
    projection has to handle through its own clause: a projected genome is in no `member_genome_ids`,
    so the anchor's array test can never find it.
    """
    locus = make_locus(seeded, ordinal=0, label="2811")
    other = make_locus(seeded, ordinal=1, label="3126")
    session.add_all([locus, other])
    session.flush()

    common = LocusArrangement(
        locus_id=locus.locus_id,
        pangenome_id=seeded["pangenome"].pangenome_id,
        rank_within_locus=0,
        member_gene_count=90,
        member_genome_count=90,
        neighbour_slot_codes=TEN_ZEROS,
        member_genome_ids=[seeded["genome"].genome_id],
    )
    rare = LocusArrangement(
        locus_id=locus.locus_id,
        pangenome_id=seeded["pangenome"].pangenome_id,
        rank_within_locus=9,
        member_gene_count=1,
        member_genome_count=1,
        neighbour_slot_codes=WINDOW,
        member_genome_ids=[seeded["genome"].genome_id],
    )
    session.add_all([common, rare])

    outsider = Genome(
        pathogen_species_id=seeded["species"].pathogen_species_id,
        sample_id=OUTSIDER,
        sample_id_kind=seeded["genome"].sample_id_kind,
        strand_is_observed=True,
    )
    session.add(outsider)
    session.flush()
    # ⭐ The outsider's three genes as ONE `genome_genes` row — the gene layer is a block per genome,
    # and `gene_count` is what the projection loader checks every placed `flat_index` against.
    plus = pack_gene_flags(is_minus_strand=False, is_five_prime_partial=False, gff_phase=0)
    session.add(
        GenomeGenes(
            genome_id=outsider.genome_id,
            gene_count=3,
            gene_spans=encode_gene_spans(
                [100 * j + 1 for j in range(3)], [100 * j + 90 for j in range(3)], [0, 0, 0], [plus] * 3
            ),
        )
    )
    session.flush()
    return {"locus": locus, "other": other, "common": common, "rare": rare, "outsider": outsider}


def write_placements(tmp_path, rows: list[dict], *, sample_id: str = OUTSIDER, assignment_sha: str | None = None):
    """A `*_placed_windows.tsv` + `placement_summary.json` pair, exactly as nuna writes them."""
    columns = [
        "sample_id", "flat_index", "locus_label", "nearest_cosine", "agreeing_neighbours",
        "locus_modelled_genes", "available_neighbours", "placed_summed_cosine", "runner_up_label",
        "runner_up_summed_cosine", "is_contested", "copies_at_this_locus", "copy_ordinal",
        *[f"s{index}" for index in range(10)],
        "matched_arrangement_rank", "matched_arrangement_genomes", "observed_slots", "has_a_window",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        full = {
            "sample_id": sample_id, "nearest_cosine": 0.99, "agreeing_neighbours": 10,
            "locus_modelled_genes": 99, "available_neighbours": 10, "placed_summed_cosine": 9.9,
            "runner_up_label": "", "runner_up_summed_cosine": "", "is_contested": False,
            "copies_at_this_locus": 1, "copy_ordinal": 1, "matched_arrangement_rank": -1,
            "matched_arrangement_genomes": 0, "observed_slots": 10, "has_a_window": True,
            **{f"s{index}": 0 for index in range(10)},
            **row,
        }
        lines.append("\t".join("" if full[column] == "" else str(full[column]) for column in columns))
    (tmp_path / f"{sample_id}_placed_windows.tsv").write_text("\n".join(lines) + "\n")
    (tmp_path / "placement_summary.json").write_text(
        json.dumps(
            {
                "rule": RULE,
                "checking_neighbours": 10,
                "neighbours_searched": 100,
                "assignment": {"sha256": assignment_sha} if assignment_sha else {},
                "nuna_git_sha": "abc123",
            }
        )
    )
    return tmp_path


def test_a_projection_loads_and_is_counted_per_genome(session, seeded, catalogue, tmp_path):
    """The per-genome row carries the rule, the evidence and the distribution — not a single score."""
    write_placements(
        tmp_path,
        [
            {"flat_index": 0, "locus_label": "2811", "nearest_cosine": 0.99},
            {"flat_index": 1, "locus_label": "2811", "nearest_cosine": 0.95, "copies_at_this_locus": 2,
             "copy_ordinal": 2},
            {"flat_index": 2, "locus_label": "3126", "nearest_cosine": 0.42, "is_contested": True,
             "runner_up_label": "2811", "runner_up_summed_cosine": 5.5, "agreeing_neighbours": 1},
        ],
    )
    report = load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()

    assert len(report.genomes) == 1
    written = report.genomes[0]
    assert written["placed"] == 3
    assert written["loci"] == 2
    assert written["contested"] == 1

    row = session.get(ProjectedGenome, (seeded["pangenome"].pangenome_id, catalogue["outsider"].genome_id))
    assert row.rule_label == RULE
    assert row.neighbours_searched == 100 and row.neighbours_reported == 10
    assert row.nearest_cosine_minimum == pytest.approx(0.42)
    assert row.nearest_cosine_median == pytest.approx(0.95)


def test_a_modelled_genome_cannot_be_projected_onto_its_own_model(session, seeded, catalogue, tmp_path):
    """⛔ Its nearest neighbour would be itself at cosine 1.0 — a perfect score meaning nothing."""
    session.add(
        GenomeCollectionMembership(
            genome_collection_id=seeded["collection"].genome_collection_id,
            genome_id=catalogue["outsider"].genome_id,
            collection_genome_ordinal=99,
        )
    )
    session.flush()
    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "2811"}])
    with pytest.raises(ProjectionRefused, match="modelled genomes"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


def test_a_locus_label_from_another_run_is_refused_rather_than_resolved(session, seeded, catalogue, tmp_path):
    """⛔ Labels are model-private: a label from another model names a DIFFERENT locus, not nothing."""
    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "99999"}])
    with pytest.raises(ProjectionRefused, match="not in this pangenome"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


def test_a_flat_index_this_genome_does_not_have_is_refused(session, seeded, catalogue, tmp_path):
    """`flat_index` is positional and shared with the embedding store — a gap means another gene set."""
    write_placements(tmp_path, [{"flat_index": 77, "locus_label": "2811"}])
    with pytest.raises(ProjectionRefused, match="flat_index this genome does not have"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


@pytest.mark.parametrize("flat_index", [3, -1])
def test_the_range_check_is_HALF_OPEN_at_gene_count_and_closed_at_zero(
    session, seeded, catalogue, tmp_path, flat_index
):
    """⛔ The boundary the dropped foreign key used to hold. The genome has genes 0, 1, 2 — so `3`
    (`== gene_count`) and `-1` are both outside, and a `<=` or a missing lower bound would let one in.
    Each is tested alone, because a file holding both would be refused for the first and say nothing
    about the second."""
    write_placements(tmp_path, [{"flat_index": flat_index, "locus_label": "2811"}])
    with pytest.raises(ProjectionRefused, match="flat_index this genome does not have"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


def test_a_genome_with_NO_genome_genes_row_is_refused_by_NAME_not_as_an_empty_gene_set(
    session, seeded, catalogue, tmp_path
):
    """⛔ `count(*)` over an absent genome returned 0 and then called every placement unknown — a true
    verdict with a misleading reason. The unloaded genome layer is named as such."""
    session.query(GenomeGenes).filter(GenomeGenes.genome_id == catalogue["outsider"].genome_id).delete()
    session.flush()
    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "2811"}])
    with pytest.raises(ProjectionRefused, match="no genome_genes row"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


def test_a_projection_built_against_a_different_assignment_is_refused(session, seeded, catalogue, tmp_path):
    """⛔ The failure nothing on the page could show: every label would name a different locus."""
    seeded["pangenome"].assignment_sha256 = "a" * 64
    session.flush()
    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "2811"}], assignment_sha="b" * 64)
    with pytest.raises(ProjectionRefused, match="computed against assignment"):
        load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)


def test_the_window_is_matched_by_its_WHOLE_vector_not_by_a_rank(session, seeded, catalogue, tmp_path):
    """⚠ Nine slots agreeing and one differing is a different neighbourhood.

    The file claims rank 9 for both genes; only the one whose ten slots actually equal that
    arrangement's is matched, because the database checks the vector rather than the claim.
    """
    nearly = dict(zip([f"s{index}" for index in range(10)], [*WINDOW[:9], 99], strict=True))
    write_placements(
        tmp_path,
        [
            {"flat_index": 0, "locus_label": "2811", "matched_arrangement_rank": 9,
             **dict(zip([f"s{index}" for index in range(10)], WINDOW, strict=True))},
            {"flat_index": 1, "locus_label": "2811", "matched_arrangement_rank": 9, **nearly},
        ],
    )
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()
    placements = {
        placement.flat_index: placement
        for placement in session.query(ProjectedGenePlacement).order_by(ProjectedGenePlacement.flat_index)
    }
    assert placements[0].matched_locus_arrangement_id == catalogue["rare"].locus_arrangement_id
    assert placements[1].matched_locus_arrangement_id is None


def test_a_gene_alone_on_its_contig_has_no_window_which_is_not_a_window_that_matched_nothing(
    session, seeded, catalogue, tmp_path
):
    """⛔ Three outcomes; collapsing any two reports a neighbourhood comparison for a gene with none."""
    write_placements(
        tmp_path,
        [
            {"flat_index": 0, "locus_label": "2811", "has_a_window": False, "observed_slots": 0},
            {"flat_index": 1, "locus_label": "2811",
             **dict(zip([f"s{index}" for index in range(10)], WINDOW, strict=True))},
        ],
    )
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()
    placements = {
        placement.flat_index: placement for placement in session.query(ProjectedGenePlacement)
    }
    assert placements[0].neighbour_slot_codes is None
    assert placements[1].neighbour_slot_codes == WINDOW
    row = session.get(ProjectedGenome, (seeded["pangenome"].pangenome_id, catalogue["outsider"].genome_id))
    assert row.gene_without_window_count == 1
    #: the third gene of this genome has no placement row at all — a different absence again
    assert row.gene_count == 3 and row.placed_gene_count == 2 and row.gene_without_vector_count == 1


def test_the_matched_arrangement_is_offered_even_though_it_is_past_the_display_cap(
    session, seeded, catalogue, tmp_path
):
    """⭐ The projected genome's own version of `arrShown`'s rule.

    Its window matches rank 9 and the cap is 8, so without the extra clause the reader is told their
    genome's neighbourhood is #10 and the response does not contain #10.
    """
    write_placements(
        tmp_path,
        [{"flat_index": 0, "locus_label": "2811",
          **dict(zip([f"s{index}" for index in range(10)], WINDOW, strict=True))}],
    )
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()

    detail = load_locus_detail(
        session,
        pangenome_id=seeded["pangenome"].pangenome_id,
        node_label="2811",
        projected_genome_id=catalogue["outsider"].genome_id,
    )
    assert 9 in [arrangement.rank_within_locus for arrangement in detail.arrangements]
    assert detail.anchor_kind == "projected"
    assert len(detail.projected_placements) == 1
    # ⛔ And it is NOT reported as a member: `arrangement_ranks` is the membership relation.
    assert detail.anchor_arrangement_ranks == []


def test_the_agreement_is_served_with_the_denominator_that_bounds_it(session, seeded, catalogue, tmp_path):
    """⛔ A locus with m modelled genes can supply at most min(n, m) checkers.

    Serving `agreeing_neighbours` alone reports a bounded numerator as a score, and at singleton loci
    that number is 1 however good the placement is.
    """
    from bacatlas_backend.serialisers.locus_serialiser import serialise_projected_placement

    write_placements(
        tmp_path,
        [{"flat_index": 0, "locus_label": "2811", "agreeing_neighbours": 1, "available_neighbours": 1,
          "locus_modelled_genes": 1, "is_contested": True}],
    )
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()
    detail = load_locus_detail(
        session,
        pangenome_id=seeded["pangenome"].pangenome_id,
        node_label="2811",
        projected_genome_id=catalogue["outsider"].genome_id,
    )
    serialised = serialise_projected_placement(detail.projected_placements[0], detail)
    assert serialised["agreeing_neighbours"] == 1
    assert serialised["available_neighbours"] == 1
    assert serialised["is_contested"] is True


def test_a_projection_changes_no_count_the_page_already_quotes(session, seeded, catalogue, tmp_path):
    """⛔ The guard for everything else: the modelled catalogue must be untouched by a projection."""
    before = _catalogue_counts(session, seeded)
    write_placements(
        tmp_path,
        [{"flat_index": 0, "locus_label": "2811"}, {"flat_index": 1, "locus_label": "3126"}],
    )
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()
    assert _catalogue_counts(session, seeded) == before


def _catalogue_counts(session, seeded) -> dict:
    """Every count a projection must not move, read from the tables the page reads."""
    pangenome_id = seeded["pangenome"].pangenome_id
    return {
        "loci": session.execute(
            text("SELECT count(*) AS n, coalesce(sum(member_gene_count),0) AS genes, "
                 "coalesce(sum(member_genome_count),0) AS genomes FROM locus WHERE pangenome_id = :p"),
            {"p": pangenome_id},
        ).mappings().one(),
        "arrangements": session.execute(
            text("SELECT count(*) AS n, coalesce(sum(member_genome_count),0) AS genomes "
                 "FROM locus_arrangement WHERE pangenome_id = :p"),
            {"p": pangenome_id},
        ).mappings().one(),
        "collection": session.execute(
            text("SELECT count(*) AS n FROM genome_collection_membership WHERE genome_collection_id = :c"),
            {"c": seeded["collection"].genome_collection_id},
        ).mappings().one(),
        # ⚠ The gene→locus edge at the genome grain: rows and placed genes of `genome_locus_assignment`,
        # not `gene_locus_membership`, which the gene-layer redesign is retiring.
        "assignments": session.execute(
            text(
                "SELECT count(*) AS n, coalesce(sum(placed_gene_count), 0) AS placed "
                "FROM genome_locus_assignment WHERE pangenome_id = :p"
            ),
            {"p": pangenome_id},
        ).mappings().one(),
    }


def test_the_listing_names_the_rule_and_only_projected_genomes(session, seeded, catalogue, tmp_path):
    """A genome that exists but was never projected must not appear, and must not resolve."""
    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "2811"}])
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()

    rows = list_projected_genomes(session, seeded["pangenome"].pangenome_id)
    assert [row.sample_id for row in rows] == [OUTSIDER]
    assert rows[0].rule_label == RULE
    assert resolve_projected_genome(session, seeded["pangenome"].pangenome_id, OUTSIDER) is not None
    assert resolve_projected_genome(session, seeded["pangenome"].pangenome_id, "SAMEA0000001") is None


def test_re_ingesting_the_catalogue_takes_the_projection_with_it_and_SAYS_so(
    session, seeded, catalogue, tmp_path, capsys
):
    """⛔ `projected_gene_placement.locus_id` cascades from `locus`, so a silent re-ingest would leave
    a `projected_genome` row whose counts describe an empty table — a summary reading "3 genes placed
    on 2 loci" over nothing. Both go, and the loader says which."""
    from bacatlas_backend.ingest.ingest_locus_catalogue import _delete_pangenome_layer

    write_placements(tmp_path, [{"flat_index": 0, "locus_label": "2811"}])
    load_projection(session, pangenome_id=seeded["pangenome"].pangenome_id, projection_root=tmp_path)
    session.flush()

    _delete_pangenome_layer(session, seeded["pangenome"].pangenome_id)
    session.flush()
    assert session.query(ProjectedGenome).count() == 0
    assert session.query(ProjectedGenePlacement).count() == 0
    assert "re-run `--stage projection`" in capsys.readouterr().out
