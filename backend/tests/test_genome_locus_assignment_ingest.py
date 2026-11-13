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
    _load_gene_memberships,
    _locus_rows,
    _present_genome_ids_by_node,
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
    with pytest.raises(KeyError, match=rf"1 of 3 collection genomes have no `genome_genes` row \(\['{SAMPLE_C}'\]\)"):
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
