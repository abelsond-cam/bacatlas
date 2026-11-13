"""`genome_locus_assignment` and `locus.present_genome_ids` — the genome-grain gene layer, at ingest.

⭐ Both are written by `ingest_locus_catalogue` from the ONE memberships frame, so the property worth
pinning is that the blocks decode to exactly the memberships — every `(genome, flat_index)` to its locus's catalogue
ordinal, `-1` everywhere else — and that every refusal names the genome, the index or the node it is
about. Design of record: `nuna/docs/design/storing_the_gene_layer.md` §3B.

⚠ The loaders are driven directly with hand-built frames on the schema probe, because a full
`CatalogueFrames` needs the cluster artifact mirror and these are the cases a real catalogue does
not contain on purpose: a flat_index past the end, a repeated one, a `.loci` that disagrees.

⛔ Mutation-minded: the two loci are numbered 1 and 0 — the REVERSE of their row order in the frame —
so a loader that read the ordinal as the row position would write the wrong locus into every block
and the `.loci` comparison, the decode assertions and the `present_genome_ids` check would all see it.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from sqlalchemy import select

from bacatlas_backend.gene_blocks import (
    NO_ARRANGEMENT_RANK,
    UNASSIGNED_LOCUS_ORDINAL,
    decode_int32_vector,
    encode_int32_vector,
    positions_of,
)
from bacatlas_backend.ingest.catalogue_frames import CatalogueFrameError
from bacatlas_backend.ingest.ingest_locus_catalogue import (
    LOCUS_COLUMNS,
    CatalogueLoadReport,
    _delete_pangenome_layer,
    _gene_counts_by_genome,
    _genome_ids_by_sample,
    _load_gene_memberships,
    _locus_rows,
    _omitted_sections,
    _present_genome_ids_by_node,
    _refuse_unless_named_set_matches,
)
from bacatlas_backend.ingest.staging_table_loader import copy_rows
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
from bacatlas_backend.models.locus import Locus
from tests.conftest import TEN_ZEROS, make_locus

SAMPLE_A, SAMPLE_B, SAMPLE_C = "SAMEA0000001", "SAMEA0000002", "SAMEA0000003"
GENE_COUNTS = {SAMPLE_A: 5, SAMPLE_B: 4, SAMPLE_C: 3}

#: ⛔ Node `n0` is ordinal 1 and `n1` is ordinal 0 — reversed on purpose (module note).
LOCI_FRAME = pd.DataFrame({"node": ["n0", "n1"], "catalogue_ordinal": [1, 0]})
ORDINAL_OF = {"n0": 1, "n1": 0}


def _memberships(rank_absent=math.nan) -> pd.DataFrame:
    """Five genes over two loci, every shape a block can get wrong.

    * A: TWO genes at `n0` (ρ > 1) with ranks 0 and 1, one gene at `n1` with NO window; genes 2 and
      4 of its five are unassigned.
    * B: one gene at each locus; listed BEFORE A at `n1`, so a `present_genome_ids` that kept frame
      order instead of sorting would read `[B, A]`.
    * C: in the collection and in nothing at all.

    ⛔ `n0` appears FIRST in the frame and is ordinal 1. A loader numbering nodes by first appearance
    (what `pd.factorize` hands back) would give `n0` → 0 — measured: with `n1` first, that mutant
    coincided with the reversed ordinals and survived. The frame order and the ordinal order must
    differ from each other AND from the loci frame's row order for all three mutants to die.
    """
    return pd.DataFrame(
        {
            "node": ["n0", "n1", "n0", "n1", "n0"],
            "sample_id": [SAMPLE_A, SAMPLE_B, SAMPLE_A, SAMPLE_A, SAMPLE_B],
            "flat_index": [0, 0, 1, 3, 2],
            "arrangement_rank": [0.0, 2.0, 1.0, rank_absent, 0.0],
        }
    )


#: What the blocks must decode to, per sample: `(ordinals, ranks, placed, distinct)`.
EXPECTED = {
    SAMPLE_A: ([1, 1, -1, 0, -1], [0, 1, -1, -1, -1], 3, 2),
    SAMPLE_B: ([0, -1, 1, -1], [2, -1, 0, -1], 2, 2),
    SAMPLE_C: ([-1, -1, -1], [-1, -1, -1], 0, 0),
}


def _frames(memberships=None, loci=LOCI_FRAME):
    return SimpleNamespace(
        gene_memberships=_memberships() if memberships is None else memberships,
        loci=loci,
        samples=[SAMPLE_A, SAMPLE_B, SAMPLE_C],
    )


@pytest.fixture()
def cohort(session, seeded):
    """Three genomes with `genome_genes` rows, and the two loci — ids keyed by sample."""
    species_id = seeded["species"].pathogen_species_id
    genomes = {SAMPLE_A: seeded["genome"]}
    for sample in (SAMPLE_B, SAMPLE_C):
        genomes[sample] = Genome(
            pathogen_species_id=species_id,
            sample_id=sample,
            sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
            strand_is_observed=True,
        )
    session.add_all([genomes[SAMPLE_B], genomes[SAMPLE_C]])
    session.flush()
    for sample, genome in genomes.items():
        session.add(
            GenomeGenes(
                genome_id=genome.genome_id,
                gene_count=GENE_COUNTS[sample],
                gene_spans=b"\x00" * (16 * GENE_COUNTS[sample]),
            )
        )
    loci = {
        label: make_locus(
            seeded, ordinal=ordinal, label=label, member_gene_count=3 if label == "n0" else 2,
            member_genome_count=2, present_genome_ids=[genomes[SAMPLE_A].genome_id, genomes[SAMPLE_B].genome_id],
        )
        for label, ordinal in ORDINAL_OF.items()
    }
    session.add_all(loci.values())
    session.flush()
    return {
        "pangenome": seeded["pangenome"],
        "species_id": species_id,
        "genome_id_by_sample": {sample: genome.genome_id for sample, genome in genomes.items()},
        "locus_id_by_label": {label: locus.locus_id for label, locus in loci.items()},
    }


def _load(session, cohort, frames=None, *, loci_root=None):
    return _load_gene_memberships(
        session,
        frames or _frames(),
        cohort["pangenome"].pangenome_id,
        cohort["genome_id_by_sample"],
        loci_root=loci_root,
    )


def _stored(session, cohort) -> dict:
    """`{sample: row}` of what the table holds for this pangenome."""
    sample_of = {genome_id: sample for sample, genome_id in cohort["genome_id_by_sample"].items()}
    return {
        sample_of[row.genome_id]: row
        for row in session.execute(
            select(GenomeLocusAssignment).where(
                GenomeLocusAssignment.pangenome_id == cohort["pangenome"].pangenome_id
            )
        ).scalars()
    }


# ── the blocks decode to exactly the memberships ──────────────────────────────────────────────
def test_every_genome_gets_a_row_whose_blocks_decode_to_exactly_its_memberships(session, cohort):
    written = _load(session, cohort)
    assert written.memberships == 5, "every membership row of the frame is encoded"
    assert written.assignments == 3
    assert written.without_window == 1

    stored = _stored(session, cohort)
    assert set(stored) == set(EXPECTED), "every collection genome has a row, C included"
    for sample, (ordinals, ranks, placed, distinct) in EXPECTED.items():
        row = stored[sample]
        assert row.gene_count == GENE_COUNTS[sample]
        assert list(decode_int32_vector(row.gene_locus_ordinals, expected_count=row.gene_count)) == ordinals, sample
        assert list(decode_int32_vector(row.gene_arrangement_ranks, expected_count=row.gene_count)) == ranks, sample
        assert (row.placed_gene_count, row.distinct_locus_count) == (placed, distinct), sample


def test_a_genome_at_rho_above_one_carries_the_same_ordinal_at_BOTH_indexes(session, cohort):
    _load(session, cohort)
    ordinals = decode_int32_vector(_stored(session, cohort)[SAMPLE_A].gene_locus_ordinals)
    # ⛔ Both copies, never a first hit — `positions_of` is what the Sequence tab asks.
    assert positions_of(ordinals, ORDINAL_OF["n0"]) == [0, 1]


@pytest.mark.parametrize("absent", [math.nan, None], ids=["nan", "none"])
def test_a_gene_with_no_window_carries_NO_ARRANGEMENT_RANK_and_is_counted(session, cohort, absent):
    """`arrangement_rank` reaches the loader as float NaN or as object None depending on the frame's
    dtype path (`catalogue_frames._gene_memberships`); both must become `-1` and be counted once."""
    written = _load(session, cohort, _frames(_memberships(rank_absent=absent)))
    assert written.without_window == 1
    ranks = decode_int32_vector(_stored(session, cohort)[SAMPLE_A].gene_arrangement_ranks)
    assert ranks[3] == NO_ARRANGEMENT_RANK
    assert list(ranks) == EXPECTED[SAMPLE_A][1]


def test_a_genome_with_no_membership_at_all_gets_an_all_unassigned_row(session, cohort):
    _load(session, cohort)
    row = _stored(session, cohort)[SAMPLE_C]
    assert row.gene_count == 3
    assert list(decode_int32_vector(row.gene_locus_ordinals)) == [UNASSIGNED_LOCUS_ORDINAL] * 3
    assert (row.placed_gene_count, row.distinct_locus_count) == (0, 0)


def test_the_encoding_is_the_codec_of_record_byte_for_byte(session, cohort):
    _load(session, cohort)
    row = _stored(session, cohort)[SAMPLE_A]
    assert row.gene_locus_ordinals == encode_int32_vector(EXPECTED[SAMPLE_A][0])
    assert row.gene_arrangement_ranks == encode_int32_vector(EXPECTED[SAMPLE_A][1])


# ── refusals, each naming what it is about ────────────────────────────────────────────────────
@pytest.mark.parametrize("bad_index", [4, -1], ids=["past_the_end", "negative"])
def test_a_flat_index_outside_the_genome_refuses_by_sample_and_index(session, cohort, bad_index):
    memberships = _memberships()
    memberships.loc[4, "flat_index"] = bad_index  # B's gene at n0; B has 4 genes
    expected = rf"{SAMPLE_B}: membership flat_index {bad_index} is outside \[0, 4\)"
    with pytest.raises(CatalogueFrameError, match=expected):
        _load(session, cohort, _frames(memberships))


def test_a_repeated_flat_index_refuses_rather_than_overwriting(session, cohort):
    memberships = _memberships()
    memberships.loc[4, "flat_index"] = 0  # B now has two rows at flat_index 0
    with pytest.raises(CatalogueFrameError, match=rf"{SAMPLE_B}: flat_index 0 has 2 membership rows"):
        _load(session, cohort, _frames(memberships))


def test_a_collection_genome_with_no_genome_genes_row_refuses_by_name(session, cohort):
    session.execute(
        GenomeGenes.__table__.delete().where(
            GenomeGenes.genome_id == cohort["genome_id_by_sample"][SAMPLE_C]
        )
    )
    with pytest.raises(
        KeyError,
        match=rf"1 of 3 collection genomes have no `genome_genes` row and are not named as holding none "
        rf"\(\['{SAMPLE_C}'\]\)",
    ):
        _load(session, cohort)


def test_a_membership_node_the_locus_frame_does_not_hold_refuses_by_label(session, cohort):
    memberships = _memberships()
    memberships.loc[0, "node"] = "n9"
    with pytest.raises(CatalogueFrameError, match=r"1 membership node\(s\) are not in the locus frame \(\['n9'\]\)"):
        _load(session, cohort, _frames(memberships))


def test_a_fractional_rank_is_refused_rather_than_truncated(session, cohort):
    memberships = _memberships()
    memberships.loc[1, "arrangement_rank"] = 1.5
    expected = r"arrangement_rank holds values that are not non-negative whole numbers \(\[1\.5\]\)"
    with pytest.raises(CatalogueFrameError, match=expected):
        _load(session, cohort, _frames(memberships))


# ── the .loci oracle ──────────────────────────────────────────────────────────────────────────
def _write_loci(root, sample, ordinals) -> None:
    (root / f"{sample}.loci").write_bytes(encode_int32_vector(ordinals))


def test_the_loci_oracle_counts_checked_and_absent_separately(session, cohort, tmp_path):
    _write_loci(tmp_path, SAMPLE_A, EXPECTED[SAMPLE_A][0])
    _write_loci(tmp_path, SAMPLE_B, EXPECTED[SAMPLE_B][0])
    written = _load(session, cohort, loci_root=tmp_path)
    assert (written.loci_checked, written.loci_absent) == (2, 1)
    assert written.loci_root == str(tmp_path)
    assert written.assignments == 3


def test_a_loci_file_that_disagrees_refuses_naming_the_sample_and_the_first_index(session, cohort, tmp_path):
    _write_loci(tmp_path, SAMPLE_A, EXPECTED[SAMPLE_A][0])
    disagreeing = list(EXPECTED[SAMPLE_B][0])
    disagreeing[2] = 0  # ours says ordinal 1 at flat_index 2
    _write_loci(tmp_path, SAMPLE_B, disagreeing)
    with pytest.raises(
        CatalogueFrameError,
        match=rf"{SAMPLE_B}: gene_locus_ordinals disagrees with the published {SAMPLE_B}\.loci at 1 of 4 genes, "
        r"first at flat_index 2 \(ours 1, the file's 0\)",
    ):
        _load(session, cohort, loci_root=tmp_path)


def test_a_loci_file_of_the_wrong_length_refuses_naming_the_sample(session, cohort, tmp_path):
    _write_loci(tmp_path, SAMPLE_A, EXPECTED[SAMPLE_A][0][:-1])
    with pytest.raises(CatalogueFrameError, match=rf"{SAMPLE_A}: the published {SAMPLE_A}\.loci does not hold 5 genes"):
        _load(session, cohort, loci_root=tmp_path)


def test_without_a_root_the_oracle_is_reported_as_NOT_looked_at(session, cohort):
    written = _load(session, cohort)
    assert written.loci_root is None
    assert (written.loci_checked, written.loci_absent) == (0, 0)
    report = CatalogueLoadReport(pangenome_id=1, genome_locus_assignments=3)
    assert ".loci oracle NOT looked at" in report.render()
    report = CatalogueLoadReport(
        pangenome_id=1, genome_locus_assignments=3, loci_oracle_checked=2, loci_oracle_absent=1, loci_oracle_root="/x"
    )
    assert "2 checked byte-for-byte against a .loci, 1 with no .loci under /x" in report.render()


# ── present_genome_ids ────────────────────────────────────────────────────────────────────────
def test_present_genome_ids_are_sorted_distinct_and_one_per_member_genome(cohort):
    ids = cohort["genome_id_by_sample"]
    a, b = ids[SAMPLE_A], ids[SAMPLE_B]
    assert a < b, "the fixture's point: B is listed before A in the frame, so sorting is observable"
    present = _present_genome_ids_by_node(_frames(), ids)
    # n0 has A twice (ρ > 1) and B once → two distinct genomes; n1 has B then A in the frame → sorted.
    assert present == {"n0": [a, b], "n1": [a, b]}


def _loci_frame_with_every_column(member_genome_counts: dict) -> pd.DataFrame:
    """A loci frame with every column `_locus_rows` reads, so the real row builder can be driven."""
    rows = []
    for label, ordinal in ORDINAL_OF.items():
        row = {
            "node": label,
            "catalogue_ordinal": ordinal,
            "member_gene_count": 3 if label == "n0" else 2,
            "member_genome_count": member_genome_counts[label],
            "prevalence_band": "core",
            "bakta_gene_symbol": None,
            "named_member_count": 0,
            "medoid_sample_id": SAMPLE_A,
            "medoid_flat_index": 0,
            "modal_cog_category": None,
            "context_observed_member_counts": TEN_ZEROS,
            "total_arrangement_count": 1,
            "arrangement_member_gene_count": 1,
            "arrangement_member_genome_count": 1,
        }
        for name in (
            "uniref50_family_count", "uniref50_major_family_count", "uniref50_labelled_member_count",
            "pfam_annotated_member_count", "pfam_architecture_count", "pfam_concordance_class",
            "syntenic_a5", "collapse_tier", "collapse_bucket", "resolved_threshold", "seqid_coverage",
            "uniref50_impurity", "uniref50_coverage", "cog_annotated_member_count", "cog_distinct_id_count",
            "ec_annotated_member_count", "kegg_annotated_member_count", "go_annotated_member_count_any",
            "propagation_verdict_cog", "propagation_verdict_ec", "propagation_verdict_kegg",
            "propagation_verdict_go", "median_gene_length_nt", "gene_length_interquartile_range_nt",
        ):
            row[name] = None
        for namespace in ("molecular_function", "biological_process", "cellular_component"):
            row[f"go_annotated_member_count_{namespace}"] = None
            row[f"go_claiming_member_count_{namespace}"] = None
            row[f"go_verdict_{namespace}"] = None
        rows.append(row)
    return pd.DataFrame(rows)


def _derived(n: int) -> dict:
    return {
        "display_name": [f"name{i}" for i in range(n)],
        "display_name_source": ["node_label"] * n,
        "display_name_source_accession": [None] * n,
        "best_product": [None] * n,
        "search_text": [""] * n,
        "interest_score": [0.0] * n,
    }


def test_the_locus_rows_carry_present_genome_ids_and_load(session, seeded):
    """Driven through the REAL row builder and `copy_rows`, so the column is proven at the table."""
    species_id = seeded["species"].pathogen_species_id
    genome_b = Genome(
        pathogen_species_id=species_id, sample_id=SAMPLE_B,
        sample_id_kind=SampleIdentifierKind.BIOSAMPLE, strand_is_observed=True,
    )
    session.add(genome_b)
    session.flush()
    ids = {SAMPLE_A: seeded["genome"].genome_id, SAMPLE_B: genome_b.genome_id, SAMPLE_C: genome_b.genome_id + 10_000}
    frames = _frames(loci=_loci_frame_with_every_column({"n0": 2, "n1": 2}))
    rows = list(_locus_rows(frames, _derived(2), seeded["pangenome"].pangenome_id, species_id, ids))
    position = LOCUS_COLUMNS.index("present_genome_ids")
    assert [row[position] for row in rows] == [[ids[SAMPLE_A], ids[SAMPLE_B]]] * 2

    copy_rows(session, Locus.__table__, LOCUS_COLUMNS, rows)
    stored = dict(
        session.execute(
            select(Locus.node_label, Locus.present_genome_ids).where(
                Locus.pangenome_id == seeded["pangenome"].pangenome_id
            )
        ).all()
    )
    assert stored == {"n0": [ids[SAMPLE_A], ids[SAMPLE_B]], "n1": [ids[SAMPLE_A], ids[SAMPLE_B]]}


def test_a_member_genome_count_the_memberships_do_not_support_refuses_by_node(seeded):
    ids = {SAMPLE_A: 1, SAMPLE_B: 2, SAMPLE_C: 3}
    frames = _frames(loci=_loci_frame_with_every_column({"n0": 2, "n1": 3}))
    with pytest.raises(
        CatalogueFrameError,
        match=r"locus 'n1': the memberships name 2 distinct genome\(s\) but member_genome_count says 3",
    ):
        list(_locus_rows(frames, _derived(2), seeded["pangenome"].pangenome_id, 1, ids))


def test_present_genome_ids_cannot_be_built_without_memberships(seeded):
    frames = _frames(memberships=None, loci=_loci_frame_with_every_column({"n0": 2, "n1": 2}))
    frames.gene_memberships = None
    with pytest.raises(CatalogueFrameError, match="the frames carry no gene memberships"):
        list(_locus_rows(frames, _derived(2), seeded["pangenome"].pangenome_id, 1, {}))


# ── the layer delete ──────────────────────────────────────────────────────────────────────────
def test_deleting_the_pangenome_layer_takes_the_assignment_rows_with_it(session, cohort):
    _load(session, cohort)
    assert len(_stored(session, cohort)) == 3
    _delete_pangenome_layer(session, cohort["pangenome"].pangenome_id)
    assert _stored(session, cohort) == {}


def test_the_vectors_are_numpy_int32_before_encoding():
    """The block builder's dtype is what keeps a 374 M-row frame off the Python heap."""
    from bacatlas_backend.ingest.ingest_locus_catalogue import _assignment_block

    block = _assignment_block("s", 4, np.array([0, 2]), np.array([7, 7]), np.array([1, -1], dtype=np.int32))
    assert block["ordinals"].dtype == np.int32 and block["ranks"].dtype == np.int32
    assert block["ordinals"].tolist() == [7, -1, 7, -1]
    assert (block["placed_gene_count"], block["distinct_locus_count"]) == (2, 1)


# ── a genome NAMED as holding no gene table ────────────────────────────────────────────────────
NAMED_REASON = "the .nseq acceptance gate refused it: its stored proteins do not reproduce from the assembly"


def _name_without_gene_table(session, cohort, sample, *, coding_gene_count):
    """What the genome stage leaves for a named genome: no block, the `_meta` row count, the reason."""
    genome_id = cohort["genome_id_by_sample"][sample]
    session.execute(GenomeGenes.__table__.delete().where(GenomeGenes.genome_id == genome_id))
    session.execute(
        Genome.__table__.update().where(Genome.genome_id == genome_id).values(
            gene_table_omitted_reason=NAMED_REASON, coding_gene_count=coding_gene_count
        )
    )


def test_a_NAMED_genome_is_assigned_like_any_other_at_its_coding_gene_count(session, cohort):
    """⭐ B carries two memberships; with its block gone and the reason recorded, its row is written at
    `genome.coding_gene_count` and decodes to exactly the memberships, as it did with a block."""
    _name_without_gene_table(session, cohort, SAMPLE_B, coding_gene_count=GENE_COUNTS[SAMPLE_B])
    written = _load(session, cohort)
    assert written.assignments == 3
    row = _stored(session, cohort)[SAMPLE_B]
    assert row.gene_count == GENE_COUNTS[SAMPLE_B]
    assert list(decode_int32_vector(row.gene_locus_ordinals, expected_count=row.gene_count)) == EXPECTED[SAMPLE_B][0]


def test_an_UNNAMED_genome_with_no_block_is_refused_EVEN_WITH_a_coding_gene_count(session, cohort):
    """⛔ The cause is the reason, not the count: the cohort fixture's genomes carry no count, so without
    one set here a loader that used ANY genome's `coding_gene_count` would pass the refusal test above."""
    genome_id = cohort["genome_id_by_sample"][SAMPLE_C]
    session.execute(GenomeGenes.__table__.delete().where(GenomeGenes.genome_id == genome_id))
    session.execute(Genome.__table__.update().where(Genome.genome_id == genome_id).values(coding_gene_count=3))
    with pytest.raises(KeyError, match=rf"1 of 3 collection genomes have no `genome_genes` row .*{SAMPLE_C}"):
        _load(session, cohort)


def test_a_NAMED_genome_with_no_coding_gene_count_is_still_refused(session, cohort):
    """A reason with no count is not a count: the block would have no length."""
    _name_without_gene_table(session, cohort, SAMPLE_C, coding_gene_count=None)
    with pytest.raises(KeyError, match=rf"1 of 3 collection genomes have no `genome_genes` row .*{SAMPLE_C}"):
        _load(session, cohort)


def test_the_gene_counts_name_exactly_the_genomes_without_a_block_and_their_reason(session, cohort):
    _name_without_gene_table(session, cohort, SAMPLE_C, coding_gene_count=3)
    samples = [SAMPLE_A, SAMPLE_B, SAMPLE_C]
    counts = _gene_counts_by_genome(session, [cohort["genome_id_by_sample"][s] for s in samples], samples)
    assert counts.omitted_reasons == {SAMPLE_C: NAMED_REASON}
    assert counts.by_genome == {cohort["genome_id_by_sample"][s]: GENE_COUNTS[s] for s in samples}


@pytest.mark.parametrize(
    "listed,expected",
    [
        (frozenset(), rf"1 collection genome\(s\) hold no gene table but are not on the list: \['{SAMPLE_C}'\]"),
        (
            frozenset({SAMPLE_C, SAMPLE_B}),
            rf"1 listed genome\(s\) are not collection genomes without a gene table .*\['{SAMPLE_B}'\]",
        ),
    ],
    ids=["actual_not_listed", "listed_not_actual"],
)
def test_the_named_set_must_equal_the_genomes_without_a_block_both_ways(listed, expected):
    with pytest.raises(KeyError, match=expected):
        _refuse_unless_named_set_matches(listed, {SAMPLE_C: NAMED_REASON})
    _refuse_unless_named_set_matches(frozenset({SAMPLE_C}), {SAMPLE_C: NAMED_REASON})


def test_omitted_sections_STATES_the_count_the_reason_and_the_genomes(session, cohort):
    """⭐ The footer prints it as "Not shown on this build: …", so it must read as a sentence on its own."""
    _name_without_gene_table(session, cohort, SAMPLE_C, coding_gene_count=3)
    samples = [SAMPLE_A, SAMPLE_B, SAMPLE_C]
    counts = _gene_counts_by_genome(session, [cohort["genome_id_by_sample"][s] for s in samples], samples)
    sections = _omitted_sections({"strand": "kept"}, counts, roster_size=3)
    assert sections["strand"] == "kept"
    assert sections["gene_table"] == (
        f"the gene table (Sequence tab) of 1 of 3 genomes — {NAMED_REASON} ({SAMPLE_C}). They are "
        "clustered, assigned and counted like every other genome"
    )
    nothing_named = _gene_counts_by_genome(session, [cohort["genome_id_by_sample"][SAMPLE_A]], [SAMPLE_A])
    assert _omitted_sections({}, nothing_named, roster_size=1) is None


def test_the_genome_lookups_survive_MORE_genomes_than_postgres_has_bind_parameters(session, cohort):
    """⛔ 71,435 genomes. `IN (…)` sends one bind parameter per genome and Postgres carries at most
    65,535, so the pangenome stage could not resolve its own roster. Measured: psycopg refuses it."""
    many = [f"SAMEA_ABSENT_{index}" for index in range(70_000)]
    assert _genome_ids_by_sample(session, [SAMPLE_A, *many]) == {SAMPLE_A: cohort["genome_id_by_sample"][SAMPLE_A]}
    with pytest.raises(KeyError, match="70000 of 70000 collection genomes have no `genome_genes` row"):
        _gene_counts_by_genome(session, list(range(-70_000, 0)), many)


# ── the list reaches the pangenome stage, and is checked there BEFORE a locus row is written ──
@pytest.mark.parametrize(
    "listed,refusal",
    [(frozenset({SAMPLE_C}), None), (frozenset(), "not on the list"), (frozenset({SAMPLE_C, SAMPLE_A}), "are not")],
    ids=["matches", "named_genome_unlisted", "listed_genome_has_a_block"],
)
def test_ingest_locus_catalogue_CHECKS_the_list_before_writing_a_locus(session, cohort, monkeypatch, listed, refusal):
    """⛔ The real entry function, its heavy derivations stubbed: the named-set check runs right after the
    gene counts, and a matching list gets past it to the first locus write — which the stub stops."""
    from bacatlas_backend.ingest import ingest_locus_catalogue as module

    class ReachedTheLocusRows(Exception):
        pass

    def stop(*_args, **_kwargs):
        raise ReachedTheLocusRows

    monkeypatch.setattr(module, "gene_ontology_fold", lambda session: None)
    monkeypatch.setattr(module, "_policy", dict)
    monkeypatch.setattr(module, "load_pfam_name_table", dict)
    monkeypatch.setattr(module, "_derive_locus_columns", lambda *args: {})
    monkeypatch.setattr(module, "_locus_rows", stop)
    _name_without_gene_table(session, cohort, SAMPLE_C, coding_gene_count=3)
    frames = SimpleNamespace(**vars(_frames()), coverage={}, omitted={})
    expected = pytest.raises(KeyError, match=refusal) if refusal else pytest.raises(ReachedTheLocusRows)
    with expected:
        module.ingest_locus_catalogue(
            session, artifacts=None, pangenome_id=cohort["pangenome"].pangenome_id,
            pathogen_species_id=cohort["species_id"], genome_id_by_ordinal=[], frames=frames,
            genomes_without_gene_table=listed,
        )


def test_the_list_is_THREADED_from_the_command_line_to_the_catalogue_load(monkeypatch, tmp_path):
    """⛔ A flag parsed and never threaded is a silent success: the stage would accept any genome the
    genome stage marked, list or no list. Driven through `main`, every stage function stubbed."""
    from bacatlas_backend.ingest import ingest_command_line as cli

    seen: dict = {}
    listed = tmp_path / "without.txt"
    listed.write_text(f"{SAMPLE_C}\n")

    class Artifacts:
        run_id, species_key = "r", "probe"

        def __init__(self, **_kwargs):
            pass

        def verify(self):
            return {}

        def required(self):
            return []

    def record(session, artifacts, **kwargs):
        seen.update(kwargs)
        return SimpleNamespace(render=lambda: "catalogue")

    monkeypatch.setattr(cli, "CatalogueArtifacts", Artifacts)
    monkeypatch.setattr(cli, "ingest_pathogen_species", lambda session: {"probe": 1})
    monkeypatch.setattr(cli, "ingest_genome_collection", lambda *a, **k: (1, SimpleNamespace(
        genome_count=3, genes_in_universe=12, collection_key="c", collection_label="C")))
    monkeypatch.setattr(cli, "ingest_pangenome_run", lambda *a, **k: (1, SimpleNamespace(render=lambda: "run")))
    monkeypatch.setattr(cli, "genome_id_by_ordinal", lambda *a: [])
    monkeypatch.setattr(cli, "ingest_locus_catalogue", record)
    monkeypatch.setattr(cli, "reconcile_pangenome", lambda *a: [])
    from tests.conftest import PROBE_URL

    code = cli.main([
        "--database-url", PROBE_URL, "--model-label", "m", "--run-id", "r", "--stage", "pangenome",
        "--species-key", "probe", "--genomes-without-gene-table", str(listed),
    ])
    assert code == 0
    assert seen["genomes_without_gene_table"] == frozenset({SAMPLE_C})


# ── the .loci oracle, read ahead on a pool — against the serial oracle it replaced ─────────────
def _old_refuse_unless_loci_agree(path, sample, gene_count, ordinals):
    """⛔ ORACLE — `_refuse_unless_loci_agree` at bacatlas 73e59b0, verbatim but for the docstring."""
    from bacatlas_backend.gene_blocks import GeneBlockError, read_loci

    try:
        published = read_loci(path, gene_count=gene_count)
    except GeneBlockError as error:
        raise CatalogueFrameError(
            f"{sample}: the published {path.name} does not hold {gene_count:,} genes "
            f"(`genome_genes.gene_count`): {error}"
        ) from error
    expected = np.frombuffer(published, dtype=np.int32)
    differing = np.flatnonzero(expected != ordinals)
    if len(differing):
        first = int(differing[0])
        raise CatalogueFrameError(
            f"{sample}: gene_locus_ordinals disagrees with the published {path.name} at "
            f"{len(differing):,} of {gene_count:,} genes, first at flat_index {first} "
            f"(ours {int(ordinals[first])}, the file's {int(expected[first])}). The block and the "
            "file are two encodings of one assignment, so this catalogue was not built from the "
            "frame the file was."
        )


def _old_assignment_rows(session, frames, pangenome_id, genome_id_by_sample, *, loci_root):
    """⛔ ORACLE — `_load_gene_memberships`' block loop at 73e59b0: a serial `is_file()` and a read per genome.

    Returns `(rows, (checked, absent))` instead of COPYing, so the two can be compared row for row.
    """
    from bacatlas_backend.ingest.ingest_locus_catalogue import (
        _arrangement_ranks_of,
        _assignment_block,
        _catalogue_ordinals_of,
        _genome_ids_of_samples,
    )

    memberships = frames.gene_memberships
    genome_ids_in_collection = [genome_id_by_sample[sample] for sample in frames.samples]
    gene_count_by_genome = _gene_counts_by_genome(session, genome_ids_in_collection, frames.samples).by_genome
    ordinals = _catalogue_ordinals_of(memberships["node"], frames.loci)
    resolved = _genome_ids_of_samples(memberships["sample_id"], genome_id_by_sample)
    sample_codes, sample_of_code = resolved["codes"], resolved["samples"]
    flat_indexes = np.asarray(memberships["flat_index"].to_numpy(), dtype=np.int64)
    ranks, _ = _arrangement_ranks_of(memberships["arrangement_rank"])
    order = np.argsort(sample_codes, kind="stable")
    sorted_codes = sample_codes[order]
    starts = np.searchsorted(sorted_codes, np.arange(len(sample_of_code)), side="left")
    ends = np.searchsorted(sorted_codes, np.arange(len(sample_of_code)), side="right")
    slice_of_sample = {
        sample_of_code[code]: (int(starts[code]), int(ends[code])) for code in range(len(sample_of_code))
    }
    rows, checked, absent = [], 0, 0
    for sample, genome_id in zip(frames.samples, genome_ids_in_collection, strict=True):
        gene_count = gene_count_by_genome[genome_id]
        start, end = slice_of_sample.get(sample, (0, 0))
        picked = order[start:end]
        block = _assignment_block(sample, gene_count, flat_indexes[picked], ordinals[picked], ranks[picked])
        oracle = loci_root / f"{sample}.loci"
        if oracle.is_file():
            _old_refuse_unless_loci_agree(oracle, sample, gene_count, block["ordinals"])
            checked += 1
        else:
            absent += 1
        rows.append((
            pangenome_id, genome_id, gene_count, block["placed_gene_count"], block["distinct_locus_count"],
            encode_int32_vector(block["ordinals"], name=f"{sample} gene_locus_ordinals"),
            encode_int32_vector(block["ranks"], name=f"{sample} gene_arrangement_ranks"),
        ))
    return rows, (checked, absent)


def _outcome(call):
    """`("ok", value)` or `("raised", type, message)` — so a refusal is compared as exactly as a row."""
    try:
        return ("ok", call())
    except Exception as error:  # noqa: BLE001
        return ("raised", type(error).__name__, str(error))


def _shape_loci_tree(root, case: str) -> None:
    """Every way a `.loci` can be present, absent, or wrong — each answered by `is_file()` + a read before."""
    _write_loci(root, SAMPLE_A, EXPECTED[SAMPLE_A][0])
    _write_loci(root, SAMPLE_B, EXPECTED[SAMPLE_B][0])
    if case == "c_is_a_directory":
        (root / f"{SAMPLE_C}.loci").mkdir()
    elif case == "c_is_a_dangling_link":
        (root / f"{SAMPLE_C}.loci").symlink_to(root / "nowhere.loci")
    elif case == "c_is_a_live_link":
        target = root / "elsewhere"
        target.mkdir()
        _write_loci(target, SAMPLE_C, EXPECTED[SAMPLE_C][0])
        (root / f"{SAMPLE_C}.loci").symlink_to(target / f"{SAMPLE_C}.loci")
    elif case == "b_disagrees":
        _write_loci(root, SAMPLE_B, [0, 0, 0, 0])
    elif case == "a_too_short":
        _write_loci(root, SAMPLE_A, EXPECTED[SAMPLE_A][0][:-1])
    elif case == "b_unreadable":
        (root / f"{SAMPLE_B}.loci").chmod(0)
    elif case == "no_directory":
        for path in root.iterdir():
            path.unlink()
        root.rmdir()


@pytest.mark.parametrize(
    "case",
    ["all_agree", "c_is_a_directory", "c_is_a_dangling_link", "c_is_a_live_link", "b_disagrees", "a_too_short",
     "b_unreadable", "no_directory"],
)
@pytest.mark.parametrize("workers", [0, 1, 3])
def test_the_loci_oracle_read_ahead_is_the_serial_oracle_row_for_row(session, cohort, tmp_path, monkeypatch,
                                                                      case, workers):
    """⭐ Same rows in the same order, same checked/absent counts, and the same refusal with the same words."""
    import os

    from bacatlas_backend.ingest import ingest_locus_catalogue

    if case == "b_unreadable" and os.geteuid() == 0:
        pytest.skip("root reads a mode-000 file")
    root = tmp_path / "loci"
    root.mkdir()
    _shape_loci_tree(root, case)
    frames = _frames()
    args = (session, frames, cohort["pangenome"].pangenome_id, cohort["genome_id_by_sample"])
    expected = _outcome(lambda: _old_assignment_rows(*args, loci_root=root))

    collected: list = []
    monkeypatch.setattr(
        ingest_locus_catalogue, "copy_rows", lambda _s, _t, _c, rows: len([collected.append(r) for r in rows])
    )

    def new():
        written = _load_gene_memberships(*args, loci_root=root, read_workers=workers)
        return collected, (written.loci_checked, written.loci_absent)

    try:
        actual = _outcome(new)
    finally:
        if case == "b_unreadable":
            (root / f"{SAMPLE_B}.loci").chmod(0o644)
    # ⛔ Coverage first: an outcome of each kind, and three rows of seven columns where it loaded.
    assert expected[0] in ("ok", "raised")
    if expected[0] == "ok":
        assert len(expected[1][0]) == 3 and all(len(row) == 7 for row in expected[1][0])
        assert sum(expected[1][1]) == 3
    assert actual == expected
    if case in ("b_disagrees", "a_too_short", "b_unreadable"):
        assert expected[0] == "raised", "the edge case must actually refuse, or it tests nothing"


def test_the_arrangement_count_is_the_rows_the_old_read_back_found(session, cohort):
    """⛔ `_load_arrangements` returned a `{(node, rank): id}` dict read back for `len()`; it now returns the
    COPY's count. The oracle is that read-back, run here against the rows the call just wrote."""
    from bacatlas_backend.ingest.ingest_locus_catalogue import _load_arrangements
    from bacatlas_backend.models.locus_arrangement import LocusArrangement

    genome_ids = [cohort["genome_id_by_sample"][sample] for sample in (SAMPLE_A, SAMPLE_B, SAMPLE_C)]
    arrangements = pd.DataFrame({
        "node": ["n0", "n0", "n1"], "rank": [0, 1, 0], "genes": [2, 1, 2], "genomes": [1, 1, 2],
        "flip": [False, True, False], **{f"s{i}": [i, 0, 1] for i in range(10)},
        "gset": [[0], [0], [0, 1]],
    })
    pangenome_id = cohort["pangenome"].pangenome_id
    count = _load_arrangements(
        session, SimpleNamespace(arrangements=arrangements), pangenome_id, cohort["locus_id_by_label"], genome_ids
    )
    read_back = {
        (label, rank): arrangement_id
        for label, rank, arrangement_id in session.execute(
            select(Locus.node_label, LocusArrangement.rank_within_locus, LocusArrangement.locus_arrangement_id)
            .join(Locus, Locus.locus_id == LocusArrangement.locus_id)
            .where(LocusArrangement.pangenome_id == pangenome_id)
        ).all()
    }
    assert count == len(read_back) == 3
    assert _load_arrangements(
        session, SimpleNamespace(arrangements=arrangements.iloc[:0]), pangenome_id, {}, genome_ids
    ) == 0
