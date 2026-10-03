"""Validate the decoder against the shipped catalogue BEFORE anything uses it as an oracle.

⛔ An oracle is only worth what its own correctness is worth. Every assertion here is a property the
payload's *own* contract states, checked over the whole 17,531-locus catalogue — so a decode error
fails here rather than surfacing later as a database that "disagrees with the page".

⛔ Every test reports the coverage it actually examined. This repo's own rule: *a diff loop that
`continue`s on a missing column reports "0 differ" while silently skipping exactly the columns that
changed*. A suite that examined 300 of 17,531 loci must say 300, not "pass".
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.known_parity_exceptions import RETIRED_TIER
from tests.payload_oracle import CONTIG_END, OFFSETS, _b64_int16, _b64_int32, load_catalogue

CATALOGUE_DIR = Path(__file__).resolve().parents[2] / "data"

pytestmark = pytest.mark.skipif(
    not (CATALOGUE_DIR / "ecoli.json").is_file(),
    reason=f"published catalogue not present at {CATALOGUE_DIR}",
)


@pytest.fixture(scope="module")
def ecoli():
    return load_catalogue(CATALOGUE_DIR / "ecoli.json")


def test_the_catalogue_is_the_one_the_live_page_serves(ecoli):
    """Pin what this oracle actually is, so a re-export cannot silently change the baseline.

    ⭐ **Re-baselined 14 → 16 on 2026-10-03, and this test fired exactly as designed.** Both pages
    were re-exported at schema 16 and republished; the schema-14 files are archived at
    `~/developer/nuna/data/browser/retired_schema14/*.schema14.json`. ⚠ The pin was the ONLY thing
    that noticed — and when it moved it also exposed that the rebuild's own
    `INGESTED_PAYLOAD_SCHEMA` had read 14 since the 2026-09-25 re-ingest made the rows schema 16.
    Two stale numbers agreeing is why neither was visible until one of them moved.
    """
    assert ecoli.schema_version == 16
    assert ecoli.meta["model_label"] == "ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL"
    assert ecoli.meta["n_genomes"] == 100
    assert ecoli.n_loci == 17_531
    assert len(ecoli.meta["genomes"]) == ecoli.meta["n_genomes"]
    assert tuple(ecoli.meta["offsets"]) == OFFSETS


def test_every_run_length_array_sums_to_the_length_of_what_it_indexes(ecoli):
    """The CSR invariant. If this is wrong, every list on every locus is shifted."""
    checked = 0
    for key, block in ecoli.raw["lists"].items():
        assert sum(block["n"]) == len(block["idx"]), key
        assert len(block["n"]) == ecoli.n_loci, key
        for column in set(block) - {"n"}:
            assert len(block[column]) == len(block["idx"]), f"{key}.{column}"
        checked += 1

    arr = ecoli.raw["arr"]
    assert sum(arr["n"]) == len(arr["cnt"]) == len(arr["gen"]) == len(arr["flip"])
    assert len(arr["vec"]) == len(arr["cnt"]) * len(OFFSETS)
    # ⭐ `gid` rides on `gen`, not on `n` — a second prefix sum over a DIFFERENT array.
    assert sum(arr["gen"]) == len(arr["gid"])

    ctx = ecoli.raw["ctx"]
    assert len(ctx["n"]) == len(ctx["obs"]) == ecoli.n_loci * len(OFFSETS)
    assert sum(ctx["n"]) == len(ctx["nid"]) == len(ctx["cnt"]) == len(ctx["same"])
    assert checked == 8, f"only {checked} list blocks checked"


def test_the_union_of_a_locus_genome_sets_is_its_genome_count(ecoli):
    """⛔ The ρ>1 rule, over every locus.

    *"A genome at ρ > 1 can occupy two arrangements at one locus and then appears in both runs. So
    `sum(gen)` over a locus may exceed `nodes.genomes[i]`, while the union of its runs equals it."*
    A consumer using `COUNT(*)` where `COUNT(DISTINCT)` was meant gets a plausible, slightly-too-large
    number — which is exactly what this asserts is not the same thing.
    """
    unions_matched = 0
    sums_exceeded = 0
    for i in range(ecoli.n_loci):
        arrangements = ecoli.arrangements(i)
        if not arrangements:
            continue
        union = set()
        total = 0
        for arrangement in arrangements:
            union.update(arrangement.genome_ordinals)
            total += arrangement.genome_count
        expected = ecoli.scalar("genomes", i)
        # ⚠ Only loci whose members ALL have a recorded neighbourhood can match exactly; a member
        # with none contributes to `nodes.genomes` and to no arrangement. So the union is a subset,
        # and the strict claim is that it never EXCEEDS the count.
        assert len(union) <= expected, f"locus {i}: union {len(union)} > genomes {expected}"
        if len(union) == expected:
            unions_matched += 1
        if total > len(union):
            sums_exceeded += 1

    assert unions_matched > 0.9 * ecoli.n_loci, f"only {unions_matched:,}/{ecoli.n_loci:,} matched exactly"
    assert sums_exceeded > 0, (
        "no locus had sum(gen) > |union| — either ρ>1 double membership does not occur in this "
        "catalogue (in which case this oracle cannot test for it) or the decoder is collapsing it"
    )


def test_a_contig_end_survives_as_a_value_and_is_never_folded_away(ecoli):
    """⛔ `-1` in `arr.vec` is *the contig ends here* — a real observation, not missing data."""
    contig_end_slots = 0
    total_slots = 0
    for i in range(ecoli.n_loci):
        for arrangement in ecoli.arrangements(i):
            total_slots += len(arrangement.slot_codes)
            contig_end_slots += sum(1 for code in arrangement.slot_codes if code == CONTIG_END)

    assert contig_end_slots > 0, "no contig-end slot decoded — the sentinel is being lost"
    assert total_slots == len(ecoli.raw["arr"]["vec"])


def test_a_flipped_row_reads_slot_nine_minus_j_and_complements_the_strand_bit(ecoli):
    """⛔ `obsSlot()`. Reversing without complementing — or complementing a `-1` — renders perfectly."""
    flipped = [
        (i, arrangement)
        for i in range(ecoli.n_loci)
        for arrangement in ecoli.arrangements(i)
        if arrangement.is_flipped
    ]
    assert flipped, "no flipped arrangement in the catalogue — the schema says ~0.7 % of them are"

    for _, arrangement in flipped[:500]:
        displayed = [arrangement.displayed_slot(j) for j in range(len(OFFSETS))]
        recorded = list(arrangement.slot_codes)
        assert displayed != recorded or all(code == CONTIG_END for code in recorded)
        for j, code in enumerate(displayed):
            source = recorded[len(recorded) - 1 - j]
            if source == CONTIG_END:
                assert code == CONTIG_END, "a contig end was strand-complemented into locus 0"
            else:
                assert code >> 1 == source >> 1, "the locus changed, not just its strand"
                assert code & 1 != source & 1, "the strand bit was not complemented"

    assert len(flipped) / len(ecoli.raw["arr"]["cnt"]) < 0.05, "flip rate far above the stated 0.7 %"


def test_observed_counts_never_exceed_the_locus_size_and_expose_truncation(ecoli):
    """`obs < size` is the *only* signal that a member ran off a contig edge."""
    truncated_loci = 0
    for i in range(ecoli.n_loci):
        size = ecoli.scalar("size", i)
        observed = ecoli.observed_member_counts(i)
        assert len(observed) == len(OFFSETS)
        assert max(observed) <= size, f"locus {i}: obs {max(observed)} > size {size}"
        if min(observed) < size:
            truncated_loci += 1

    assert truncated_loci > 0, "no locus shows contig-edge truncation — obs is being read as size"


def test_every_occupant_names_a_real_locus_and_both_directions_agree(ecoli):
    """⭐ The reverse index must be the forward one read the other way, not a second computation."""
    forward: dict[tuple[int, int, int], tuple[int, int]] = {}
    reverse_count: dict[int, int] = {}
    for i in range(ecoli.n_loci):
        for occupant in ecoli.offset_occupants(i):
            assert 0 <= occupant.neighbour_locus_index < ecoli.n_loci
            assert occupant.same_strand_count <= occupant.gene_count
            forward[(i, occupant.signed_offset, occupant.rank)] = (
                occupant.neighbour_locus_index,
                occupant.gene_count,
            )
            reverse_count[occupant.neighbour_locus_index] = reverse_count.get(occupant.neighbour_locus_index, 0) + 1

    assert len(forward) == len(ecoli.raw["ctx"]["nid"])
    assert sum(reverse_count.values()) == len(forward)


def test_the_joint_and_the_marginal_are_consistent_in_the_direction_that_holds(ecoli):
    """⛔ One direction only, and the converse is DELIBERATELY not asserted.

    Every occupant in the marginal must appear in at least one arrangement slot at that offset. The
    converse — that the marginal's modes co-occur in some genome — is exactly the misreading the two
    views exist to prevent, and a test demanding it would encode the error.
    """
    checked = 0
    for i in range(0, ecoli.n_loci, 7):  # a seventh of the catalogue; the loop is O(arrangements)
        by_slot: list[set[int]] = [set() for _ in OFFSETS]
        for arrangement in ecoli.arrangements(i):
            for slot, code in enumerate(arrangement.slot_codes):
                if code != CONTIG_END:
                    by_slot[slot].add(code >> 1)
        for occupant in ecoli.offset_occupants(i):
            slot = OFFSETS.index(occupant.signed_offset)
            assert occupant.neighbour_locus_index in by_slot[slot], (
                f"locus {i} offset {occupant.signed_offset}: marginal names locus "
                f"{occupant.neighbour_locus_index}, which sits in no arrangement at that slot"
            )
            checked += 1

    assert checked > 30_000, f"only {checked:,} occupant rows checked"


def test_an_ABSENT_nearest_locus_drops_its_slot_and_is_not_a_similarity_of_zero(ecoli):
    """⛔ The rule the retired `map_reps` tests could only check SYNTHETICALLY, now on live data.

    Those three tests asserted `dropped_seen == 0` and `absent == {bacformer: 0, esm: 0}` — they
    recorded a non-occurrence, and said so honestly. The rule moved to `sim.near_i` when the medoid
    geometry was withdrawn on 2026-09-24, and there it **happens**: thousands of loci have fewer than
    `k` neighbours in the catalogue.

    ⚠ **And the hazard is sharper than it was.** `near_v` at a dropped slot is `0`, while Bacformer's
    real cross-locus medians go NEGATIVE — so a consumer that reads the value without consulting
    `near_i` reads an absent neighbour as *more similar* than a real one. That is the `gaps.vi` trap
    again, in the one block schema 16 added and nothing had tested.
    """
    absent = {}
    for representation in ("bacformer", "esm"):
        block = ecoli.sim_block(representation)
        k, scale = block["k"], block["cos_scale"]
        assert k == 5
        gone = 0
        for index in range(ecoli.n_loci):
            neighbours = ecoli.sim_neighbours(representation, index)
            slots = [slot for slot, _, _ in neighbours]
            assert slots == sorted(slots), "slots must stay in rank order"
            assert all(0 <= slot < k for slot in slots)
            assert len({target for _, target, _ in neighbours}) == len(neighbours), "a repeated neighbour"
            assert all(0 <= target < ecoli.n_loci for _, target, _ in neighbours), "an out-of-catalogue index"
            assert all(-1.0001 <= cosine <= 1.0001 for _, _, cosine in neighbours)
            gone += k - len(neighbours)
        absent[representation] = gone
        # ⛔ The rule is EXERCISED, not merely stated: if this ever reaches 0 the test has gone back
        # to proving nothing and should be re-pointed at whatever block carries the sentinel then.
        assert gone > 0, f"{representation} has no absent slot — this test is synthetic again"

        # ⚠ Every absent slot's VALUE is 0, and 0 is not below the live range.
        values = [value / scale for value in _b64_int16(block["near_v"])]
        indices = _b64_int32(block["near_i"])
        assert {value for value, target in zip(values, indices, strict=True) if target == -1} == {0.0}
        live = [value for value, target in zip(values, indices, strict=True) if target != -1]
        assert min(live) < 0.0 if representation == "bacformer" else min(live) > 0.0

    # Recorded so a change in the catalogue is visible rather than absorbed.
    assert absent == {"bacformer": 1_622, "esm": 3_525}, absent


def test_the_sim_block_measures_the_same_catalogue_the_tier_column_does(ecoli):
    """⭐ A three-way cross-check the medoid tests had no way to make.

    `within` is `null` for exactly the loci with no within-locus pair — the singletons — and the tier
    is absent for exactly the loci the audit could not place. Those are the same population, and the
    block DECLARES its own count. All three must agree or the payload is describing two catalogues.
    """
    for representation in ("bacformer", "esm"):
        block = ecoli.sim_block(representation)
        measurable = [value is not None for value in block["within"]]
        assert len(measurable) == ecoli.n_loci
        assert sum(measurable) == block["n_measurable"] == 12_104
        tiered = [value >= 0 for value in ecoli.nodes["tier"]]
        assert measurable == tiered, "the sim block and the tier column disagree about which loci exist"
    # ⚠ The two representations must agree with each other too — one catalogue, measured twice.
    assert [v is None for v in ecoli.sim_block("esm")["within"]] == [
        v is None for v in ecoli.sim_block("bacformer")["within"]
    ]


def test_absence_from_the_sparse_variance_index_is_a_MEASURED_ZERO(ecoli):
    """⛔ The rule that inverts, and whose reversal is silent.

    `app.js::gapVarRaw` has two clauses and they mean opposite things: the whole block absent is
    *not measured* (`null`), but a **gap** absent from `vi` scored exactly zero — every genome
    agrees. 86.1 % of this catalogue's gaps are in that second case, so reading absence as `None`
    would blank the majority of the block and turn *"identical in every genome"* into
    *"we never looked"*.
    """
    gaps = ecoli.intergenic_gaps
    assert gaps, "no gaps block decoded"
    indexed = len(ecoli.raw["gaps"]["vi"])
    zero = [gap for gap in gaps if gap.variance_score == 0.0]
    not_measured = [gap for gap in gaps if gap.variance_score is None]

    assert not_measured == [], (
        f"{len(not_measured):,} gaps decoded as not-measured, but this payload CARRIES the variance "
        "block — every gap in it has a score, and most of them are zero"
    )
    # The absent majority, plus the handful that are indexed and happen to score zero.
    assert len(zero) >= len(gaps) - indexed, f"{len(zero):,} zeros against {len(gaps) - indexed:,} absent"
    assert 0.80 < (len(gaps) - indexed) / len(gaps) < 0.92, "the absent share moved far from the stated 86 %"

    # ⚠ `vd` is the score × 1000. A decoder comparing it raw against a 0..1 threshold reads every
    # varying gap as saturated, and the track goes uniformly dark.
    assert max(ecoli.raw["gaps"]["vd"]) == 2000
    assert max(gap.variance_score for gap in gaps) == 2.0

    # ⚠ The modal length is NOT recovered for an absent gap — one absence, two different fates.
    absent_with_mode = [g for g in gaps if g.variance_score == 0.0 and g.modal_length_nt is not None]
    assert len(absent_with_mode) <= indexed



def test_the_flanking_pair_is_sorted_by_NODE_LABEL_and_not_by_payload_index(ecoli):
    """⛔ The finding this oracle was written to be capable of making.

    `gene_adjacencies` sorts the pair with ``left <= right`` on the **node ids**, which are TEXT.
    `intergenic_block` then maps both through `node_pos` — `node_order`'s band-and-size ordering —
    so the shipped `a`/`b` are label-sorted and their payload indices are in no particular order.
    On this catalogue **32.3 % of gaps have `a > b` by index**.

    That matters because `app.js::gapBetween` keys its lookup by INDEX
    (``(x < y ? x : y) + "|" + (x < y ? y : x)``) against a map built from the shipped order — so
    those gaps cannot be found. The database therefore canonicalises on `node_label`, matching the
    science, and the API resolves the pair by label before it queries.
    """
    labels = ecoli.nodes["label"]
    by_label = sum(
        1 for gap in ecoli.intergenic_gaps
        if str(labels[gap.low_locus_index]) <= str(labels[gap.high_locus_index])
    )
    by_index = sum(
        1 for gap in ecoli.intergenic_gaps if gap.low_locus_index <= gap.high_locus_index
    )
    total = len(ecoli.intergenic_gaps)

    assert by_label == total, f"only {by_label:,}/{total:,} pairs are label-sorted"
    assert by_index < total, (
        "every pair is ALSO index-sorted, so the two orderings agree on this catalogue and the "
        "distinction below is untestable here"
    )
    unfindable = total - by_index
    assert unfindable / total > 0.25, f"{unfindable:,}/{total:,} index-inverted — expected ~32 %"


def test_an_overlap_is_a_negative_gap_and_is_not_clamped(ecoli):
    """*"Clamping this at zero made abuts exactly and overlaps by 190 bases the same number."*"""
    overlaps = [gap for gap in ecoli.intergenic_gaps if gap.median_signed_length_nt < 0]
    share = len(overlaps) / len(ecoli.intergenic_gaps)
    assert 0.05 < share < 0.40, f"{share:.1%} of adjacencies overlap — the schema measured 18.8 %"


def test_only_min_equals_max_certifies_that_every_genome_agrees(ecoli):
    """⛔ `q1 == q3` means the MIDDLE HALF agrees. `[1,1,1,1,100]` has an IQR of zero."""
    ranged = [g for g in ecoli.intergenic_gaps if g.minimum is not None and g.quartile1 is not None]
    assert ranged, "no gap carries both a range and quartiles"
    misleading = [g for g in ranged if g.quartile1 == g.quartile3 and g.minimum != g.maximum]
    assert misleading, (
        "no gap has a degenerate IQR over a real spread — the distinction the schema records as a "
        "reporting bug does not occur here, so this catalogue cannot test it"
    )


def test_an_absent_string_index_is_none_and_never_the_zeroth_string(ecoli):
    """⛔ `-1` is absent. Resolving it as `strings[pool][-1]` returns the LAST string, plausibly."""
    assert ecoli.string("sym", -1) is None
    unnamed = [i for i in range(ecoli.n_loci) if ecoli.nodes["name"][i] < 0]
    assert unnamed, "every locus is named — the absent path is untested"
    assert ecoli.string("sym", ecoli.nodes["name"][unnamed[0]]) is None


def test_the_tier_pool_round_trips_the_audit_artifact_it_was_exported_from(ecoli):
    """⛔ Parity is an invariant now, not an allowlist — and the PREMISE changed, not just the list.

    This used to say: *"the local audit artifacts were regenerated after the pages were published,
    so they are not the files the pages were built from"*, and allowed the 2 loci that drift caused.
    That is no longer true. The pages were **re-exported on 2026-10-03 from these very artifacts**,
    so the right assertion is that the interned tier pool round-trips the CSV it came from, on every
    locus — which is a stronger claim than the allowlist it replaces, and one that would fail if a
    future re-export ever interned the pool against a different audit.

    ⚠ Narrowing the allowlist to an empty set would have left `assert set() == frozenset()`, green
    forever. The denominator is asserted first for the same reason: a comparison that examined
    nothing also finds nothing differing.
    """
    import csv
    import os
    from pathlib import Path

    data_root = Path(os.environ.get("BACATLAS_NUNA_DATA_ROOT", "~/developer/nuna/data")).expanduser()
    label = "ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL"
    waterfall = data_root / "proc" / "analysis" / "accessory_audit" / f"{label}_homology_waterfall.csv"
    if not waterfall.is_file():
        pytest.skip(f"the audit artifacts are not mirrored at {waterfall}")

    with waterfall.open() as handle:
        current = {row["node"]: row["bucket_esm"] for row in csv.DictReader(handle)}

    tiers = ecoli.strings["tier"]
    compared, differing = 0, {}
    for i in range(ecoli.n_loci):
        index = ecoli.nodes["tier"][i]
        if index < 0:
            continue  # -1 = no waterfall row at all, which is every singleton
        label_text = ecoli.label(i)
        if label_text not in current or not current[label_text]:
            continue
        compared += 1
        if tiers[index] != current[label_text]:
            differing[label_text] = (tiers[index], current[label_text])

    assert compared == 12_104, f"only {compared:,} loci carried a tier on both sides"
    assert differing == {}, f"the page's tier disagrees with the audit on {len(differing)} loci: {differing}"
    # ⭐ And the retired name is in neither the pool nor the artifact — the retirement as an
    # invariant rather than as a tolerated difference.
    assert RETIRED_TIER not in tiers
    assert RETIRED_TIER not in set(current.values())
