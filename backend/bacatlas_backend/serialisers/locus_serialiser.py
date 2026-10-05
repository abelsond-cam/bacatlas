"""ORM rows → plain dicts. **No ORM object escapes a serialiser.**

⛔ **A slot is `null` plus an `absence_reason`, never a bare `-1`.** `occ(code)` is deleted: that
packed form is exactly where the "−1 means five different things" trap lives, and a serialiser that
mapped every `-1` to `null` would be right twice and destroy three meanings.

⛔ **Every remainder is its own named field, one per sentence.** `observed_not_listed` is a
*display* cut; `members_without_an_observation` is *missing data* (a contig end);
`arrangements_not_listed` is what the arrangement cap left out; `members_in_arrangements_not_listed`
is the members inside those; `members_without_a_neighbourhood` is members with no recorded window at
all. A rewrite that computes "the rest" once collapses them and the page starts making a claim it
cannot support — so each is emitted, and each is emitted even when it is zero.

⚠ **There are five here and four in the plan, and the fifth is the API's own doing.** The published
payload is uncapped, so *members past the cap* did not exist as a category and `size − Σ(listed)`
really did mean "no coordinates for the gene, so no window". The moment the API capped the
arrangement list that stopped being true, and splitting the two is what keeps the fourth honest.

⚠ **`null` is *not measured* and `0.0` is *measured zero*, everywhere below.** A serialiser that
emitted `0` for an absent float would turn "we did not look" into "we looked and found nothing".
"""

from __future__ import annotations

from bacatlas_backend.models.enumerations import gene_ontology_namespace_name
from bacatlas_backend.services.locus_detail_service import (
    SIGNED_OFFSETS,
    LocusDetail,
    membership_is_complete,
)
from bacatlas_backend.services.reference_name_service import ReferenceNames

#: What a slot's absence means, when it has one. ⛔ A contig end is an OBSERVATION — the member
#: genuinely has no gene there — and is not the same as a neighbour that fell outside the catalogue.
ABSENCE_CONTIG_END = "contig_end"


def serialise_annotation_entry(entry, names: ReferenceNames | None = None) -> dict:
    """One (vocabulary, term) row. `name` is `null` where nothing can name this term.

    ⛔ **`term_name` first, the reference second — never the other way round.** Bakta's own name
    is what this catalogue's genes were annotated with; a vendored reference is a later lookup of
    the same code. Where both exist, preferring the reference would silently replace the recorded
    annotation with a different release's wording and make the page disagree with the GFFs.

    In practice they never collide: `term_name` is populated for 100 % of COG and GO entries and
    0 % of EC and KEGG ones, which is exactly why `names` exists. The ordering is the guarantee that
    it stays a gap-filler if that ever changes.
    """
    return {
        "rank": entry.rank_within_locus,
        "term": entry.term_value,
        "name": entry.term_name
        or (names.name_for(entry.annotation_kind, entry.term_value) if names else None),
        "gene_count": entry.member_gene_count,
        # ⛔ The NAME, never the stored 0/1/2. See `gene_ontology_namespace_name`: the coverage
        # block of this very response keys its namespaces by name, and a client grouping entries by
        # an integer renders three GO cards with empty term lists under coverage lines that promise
        # otherwise. Found by building the tab, invisible to every test before it existed.
        **(
            {"gene_ontology_namespace": gene_ontology_namespace_name(entry.gene_ontology_namespace)}
            if entry.gene_ontology_namespace is not None
            else {}
        ),
    }


def serialise_uniref_family(family) -> dict:
    """One UniRef50 family's row of the cross-tab.

    ⚠ `distinct_symbol_count` is a COUNT and not a list: the page shows the modal plus "+N" and
    cannot name the others, so the response must not appear to offer them.

    ⚠ `named_gene_count` is that column's coverage, and `null` is *not measured* while `0` is
    *measured and none* — the same distinction `pfam_annotated_gene_count` carries.
    """
    return {
        "rank": family.rank_within_locus,
        "uniref50_accession": family.uniref50_accession,
        "gene_count": family.member_gene_count,
        "modal_product": family.modal_bakta_product,
        "modal_architecture": family.modal_pfam_architecture,
        "pfam_annotated_gene_count": family.pfam_annotated_member_count,
        "modal_symbol": family.modal_bakta_gene_symbol,
        "distinct_symbol_count": family.distinct_real_symbol_count,
        # ⛔ The symbol column's DENOMINATOR, and it travels with the modal symbol always. Without
        # it a family of 9 genes naming 2 reads as nine genes agreeing — kp locus 3992.
        "named_gene_count": family.named_member_count,
    }


def serialise_arrangement(arrangement, neighbours) -> dict:
    """One whole ±5 neighbourhood — the JOINT view, with each slot resolved or explicitly absent."""
    slots = []
    for position, offset in enumerate(SIGNED_OFFSETS):
        code = arrangement.neighbour_slot_codes[position]
        if code < 0:
            slots.append(
                {
                    "signed_offset": offset,
                    "locus": None,
                    # ⛔ Named, not implied by the null. A member truncated at a contig edge has a
                    # genuinely different OBSERVED neighbourhood; folding it into the untruncated
                    # group would claim an observation that was never made.
                    "absence_reason": ABSENCE_CONTIG_END,
                    "same_strand": None,
                }
            )
            continue
        row = neighbours.by_catalogue_ordinal.get(code // 2)
        slots.append(
            {
                "signed_offset": offset,
                "locus": row.node_label if row else None,
                "absence_reason": None if row else "outside_catalogue",
                "same_strand": bool(code % 2),
            }
        )
    return {
        "rank": arrangement.rank_within_locus,
        "gene_count": arrangement.member_gene_count,
        "genome_count": arrangement.member_genome_count,
        # ⚠ INTRINSIC to the arrangement, not display mirroring. The two are different slot spaces,
        # and conflating them puts a count from one position under a gene from another.
        "is_recorded_reverse_complement": arrangement.is_recorded_reverse_complement,
        "slots": slots,
    }


def serialise_offset_occupants(detail: LocusDetail) -> list[dict]:
    """The MARGINAL view — one entry per offset, each with its own honest denominator."""
    observed = list(detail.locus.context_observed_member_counts)
    out = []
    for position, offset in enumerate(SIGNED_OFFSETS):
        occupants = detail.offset_occupants.get(offset, [])
        listed = sum(occupant.member_gene_count for occupant in occupants)
        observed_here = observed[position] if position < len(observed) else 0
        out.append(
            {
                "signed_offset": offset,
                "observed_member_count": observed_here,
                # ⛔ Two different remainders, both emitted, both named.
                "observed_not_listed": max(0, observed_here - listed),
                "members_without_an_observation": max(
                    0, detail.locus.member_gene_count - observed_here
                ),
                "occupants": [
                    {
                        "rank": occupant.rank_within_offset,
                        "locus": (
                            detail.neighbour_display_rows.by_locus_id[
                                occupant.neighbour_locus_id
                            ].node_label
                            if occupant.neighbour_locus_id in detail.neighbour_display_rows.by_locus_id
                            else None
                        ),
                        "gene_count": occupant.member_gene_count,
                        "same_strand_gene_count": occupant.same_strand_member_count,
                    }
                    for occupant in occupants
                ],
            }
        )
    return out


def serialise_intergenic_gap(gap, labels: dict) -> dict:
    """One gap, keyed by the two flanking LABELS — never by an index.

    ⛔ The published page keys its gap lookup by payload INDEX against a map built from
    **label**-sorted pairs, and 7,379 of 22,838 *E. coli* gaps therefore cannot be found on it: the
    track draws no block, which reads as *"these genes are adjacent with nothing between them"*. The
    API answers with the labels, so a client sorts two strings and the miss cannot recur.
    """
    return {
        "flanking_loci": [labels.get(gap.flanking_locus_id_a), labels.get(gap.flanking_locus_id_b)],
        "observed_genome_count": gap.observed_genome_count,
        # ⛔ SIGNED. Negative means the two genes OVERLAP — 18.8 % of adjacent pairs do, at a median
        # of −4 bases. Clamping at zero made *abuts exactly* and *overlaps by 190* the same number.
        "median_signed_length_nt": gap.median_signed_length_nt,
        "quartile1_signed_length_nt": gap.quartile1_signed_length_nt,
        "quartile3_signed_length_nt": gap.quartile3_signed_length_nt,
        "minimum_signed_length_nt": gap.minimum_signed_length_nt,
        "maximum_signed_length_nt": gap.maximum_signed_length_nt,
        # ⛔ `null` = not measured; `0.0` = every genome agrees. White on the track has to mean
        # "identical in every genome", never "small".
        "length_variance_score": gap.length_variance_score,
        "modal_length_nt": gap.modal_length_nt,
        "distinct_named_feature_count": gap.distinct_named_feature_count,
        # ⚠ `q1 == q3` means the MIDDLE HALF agrees; only `mn == mx` certifies the strong claim.
        "every_genome_agrees": (
            gap.minimum_signed_length_nt is not None
            and gap.minimum_signed_length_nt == gap.maximum_signed_length_nt
        ),
    }



def serialise_projected_placement(placement, detail) -> dict:
    """One gene of a projected genome at this locus — the evidence, with its own denominator.

    ⛔ **`agreeing_neighbours` is meaningless without `available_neighbours`.** A locus with *m*
    modelled genes can supply at most min(n, m) of the n checkers, so the raw count tracks the
    locus's SIZE as much as the evidence: measured over the first ten placed genomes, `is_contested`
    runs at 0.39 % where the locus has ten genes to offer and 42.6 % where it has fewer, reaching
    96.3 % at singletons. Both numbers are serialised together for that reason, and a client that
    shows one without the other is reporting a bounded numerator as a score.

    ⚠ **Three absences, kept apart.** `matched_arrangement_rank` is null with
    `has_a_neighbourhood` true → the window was built and matches none of the published
    arrangements, which is an observation. `has_a_neighbourhood` false → the gene is alone on its
    contig and has no window at all. A gene with no Bacformer vector has no placement row anywhere.
    """
    matched = None
    for arrangement in detail.arrangements:
        if arrangement.locus_arrangement_id == placement.matched_locus_arrangement_id:
            matched = arrangement
            break
    return {
        "flat_index": placement.flat_index,
        "copy_ordinal": placement.copy_ordinal,
        "copies_at_locus": placement.copies_at_locus,
        "nearest_cosine": placement.nearest_cosine,
        "agreeing_neighbours": placement.agreeing_neighbour_count,
        "available_neighbours": placement.available_neighbour_count,
        "placed_summed_cosine": placement.placed_summed_cosine,
        "is_contested": placement.is_contested,
        "runner_up_summed_cosine": placement.runner_up_summed_cosine,
        "has_a_neighbourhood": placement.neighbour_slot_codes is not None,
        "matched_arrangement_rank": matched.rank_within_locus if matched is not None else None,
        "matched_arrangement_genome_count": matched.member_genome_count if matched is not None else None,
    }


def serialise_locus_detail(detail: LocusDetail, names: ReferenceNames | None = None) -> dict:
    """The whole locus response — one round trip, and the popover is then offline.

    `names` fills in the EC and KEGG names Bakta never emitted. Optional, because an absent
    reference must degrade a chip to a bare accession rather than fail the request.

    ⚠ `cosine_scale_factor` went with the map on 2026-09-24. The 6×6 cosines were stored as scaled
    int16 and had to be divided back; the similarities are plain floats, so nothing is scaled and
    nothing needs a factor travelling beside it to be read.
    """
    locus = detail.locus
    neighbours = detail.neighbour_display_rows
    # ⛔ The card's nearest loci are stored as real `locus_id`s but SERVED as catalogue ordinals, so
    # the client resolves them through `neighbour_display_rows` like every other neighbour address —
    # one index on the page rather than two over the same small integers. A neighbour the fan-out did
    # not resolve is DROPPED rather than served with a null address: the list is already ragged, so a
    # missing entry is a shape the client handles, and a null id is one it would have to guess at.
    ordinal_of = {row.locus_id: row.catalogue_ordinal for row in neighbours.all_rows()}
    # ⛔ EVERY resolved row, both key spaces: a gap between two arrangement occupants that are not
    # marginal modes is flanked by loci reached only by catalogue ordinal, and keyed on `by_locus_id`
    # alone its labels come back `null` — so the gap is dropped client-side and two genes read as
    # adjacent with nothing between them.
    labels = {row.locus_id: row.node_label for row in neighbours.all_rows()}
    labels.setdefault(locus.locus_id, locus.node_label)

    listed_members = sum(
        arrangement.member_gene_count for arrangement in detail.arrangements
    )
    return {
        "locus": {
            "label": locus.node_label,
            "catalogue_ordinal": locus.catalogue_ordinal,
            "display_name": locus.display_name,
            "display_name_source": locus.display_name_source,
            "display_name_source_accession": locus.display_name_source_accession,
            "best_product": locus.best_product,
            "bakta_gene_symbol": locus.bakta_gene_symbol,
            "gene_count": locus.member_gene_count,
            "genome_count": locus.member_genome_count,
            "named_gene_count": locus.named_member_count,
            "prevalence_band": locus.prevalence_band.value,
            "median_gene_length_nt": locus.median_gene_length_nt,
            "gene_length_interquartile_range_nt": locus.gene_length_interquartile_range_nt,
            "uniref50": {
                "family_count": locus.uniref50_family_count,
                # ⚠ What the card CLAIMS with, while the list below shows every family.
                "major_family_count": locus.uniref50_major_family_count,
                "labelled_gene_count": locus.uniref50_labelled_member_count,
                "impurity": locus.uniref50_impurity,
                "coverage": locus.uniref50_coverage,
            },
            "pfam": {
                # ⛔ `gene_count − annotated_gene_count` is MISSING COVERAGE, not a competing
                # architecture. Without this denominator the page cannot tell the two apart.
                "annotated_gene_count": locus.pfam_annotated_member_count,
                "architecture_count": locus.pfam_architecture_count,
                # ⛔ The audit's own verdict, read and never re-derived.
                "concordance_class": locus.pfam_concordance_class,
            },
            "evidence": {
                "syntenic_a5": locus.syntenic_a5,
                "collapse_tier": locus.collapse_tier,
                "collapse_bucket": locus.collapse_bucket,
                # ⚠ Over at most 50 members, longest-first (`FAMILY_MEMBER_CAP`) — 4,138 of 12,104
                # ecoli loci sit at that cap, so this is "the identity at which the 50 longest
                # members first group", a different and weaker claim than the whole locus.
                "resolved_threshold": locus.resolved_threshold,
                "resolved_threshold_is_capped_at_50_members": True,
            },
            # ⛔ `geometry` was REPLACED, not renamed. Every field it carried measured a locus by
            # reducing it to its MEDOID — `within_medoid_distance` was its members' distance to that
            # one gene, `cosine_matrix` the six medoids' pairwise geometry, `map_position` where that
            # gene sat on a UMAP of medoids. These are measured over the whole SET, or anchored on a
            # POINT, and a client that read one as the other would be plausible and wrong.
            #
            # ⚠ `separation` and `weak_margin` are NOT here although the card prints both: each is
            # the difference of two fields that ARE, and the client subtracts. A separately-rounded
            # difference drifting from the two rounded numbers printed beside it is the exact class
            # of quiet disagreement this card was rebuilt to end.
            "similarity": {
                representation: (
                    None
                    if row is None
                    else {
                        "within_similarity": row.within_similarity,
                        "nearest_similarity": row.nearest_similarity,
                        "weak_own_similarity": row.weak_own_similarity,
                        "weak_other_similarity": row.weak_other_similarity,
                        "own_neighbour_fraction": row.own_neighbour_fraction,
                        # ⭐ One midrank per VIEW, because the three views are not three pictures of
                        # one thing — only 18–34 % of flagged loci are flagged by all three.
                        # ⚠ Each is NULL for a singleton and must read *not measurable*, never 0.000.
                        "separation_percentile": row.separation_percentile,
                        "weak_margin_percentile": row.weak_margin_percentile,
                        # ⛔ At an own fraction of exactly 1.0 this is a midrank inside a tie block
                        # covering 88.5 % of the catalogue — "p53" reading as *better than half the
                        # catalogue* when it means *tied with nearly all of it*. Served because below
                        # 1.0 it means what it looks like; the CARD is what must decline to print it.
                        "own_fraction_percentile": row.own_fraction_percentile,
                        # ⭐ The five nearest OTHER loci in THIS representation — a different five
                        # from the other's (their separations correlate at rho ~0.47). Ragged: a
                        # locus whose shortlist held fewer simply has fewer entries, never padding.
                        # ⚠ `catalogue_ordinal`, so the client resolves them through
                        # `neighbour_display_rows` like every other neighbour address.
                        "nearest_loci": [
                            {
                                "rank": nearest.rank,
                                "catalogue_ordinal": ordinal_of.get(nearest.neighbour_locus_id),
                                "cross_similarity": nearest.cross_similarity,
                            }
                            for nearest in detail.nearest_loci.get(representation, ())
                            if nearest.neighbour_locus_id in ordinal_of
                        ],
                    }
                )
                for representation, row in (
                    ("esm", detail.similarity.get("esm")),
                    ("bacformer", detail.similarity.get("bacformer")),
                )
            },
            "interest_score": locus.interest_score,
        },
        "annotations": {
            kind: [serialise_annotation_entry(entry, names) for entry in entries]
            for kind, entries in detail.card_annotations.items()
        },
        "uniref50_families": [serialise_uniref_family(family) for family in detail.uniref_families],
        "arrangements": {
            "listed": [
                serialise_arrangement(arrangement, neighbours) for arrangement in detail.arrangements
            ],
            # ⛔ Never moved by any display cap, and never conflated with the number listed.
            "total": locus.total_arrangement_count,
            "arrangements_not_listed": max(
                0, locus.total_arrangement_count - detail.arrangements_listed
            ),
            # ⛔⛔ TWO remainders, and the API needs both where the published page needed only one.
            # That payload is UNCAPPED, so `size − Σ(listed)` there really did mean "no coordinates
            # for the gene, so no window". Capping the list at 8 makes the same subtraction sweep in
            # every member sitting past the cap — 15,912 E. coli genes over 2,340 loci and 10,437 kp
            # genes over 1,548 — each of which the page would then have told a reader has no
            # coordinates. `arrangement_member_gene_count` is the denominator that separates them.
            "members_in_arrangements_not_listed": max(
                0, locus.arrangement_member_gene_count - listed_members
            ),
            "members_without_a_neighbourhood": max(
                0, locus.member_gene_count - locus.arrangement_member_gene_count
            ),
            # ⛔ Whether EVERY genome present reaches an arrangement — a GENOME question, and the
            # only thing that settles which of two sentences the anchor line may say. Empty anchor
            # ranks mean *"has no gene at this locus"* only where membership is complete; where it
            # is not, the true sentence is *"has no recorded neighbourhood at this locus"*, and at
            # 6.26 % of ecoli loci the first one is false. ⚠ Not derivable from the gene counts
            # above: a genome at ρ > 1 can lose one gene's window and keep its arrangement.
            "membership_is_complete": membership_is_complete(locus),
        },
        # ⛔ A LIST, and it is not decoration: without it a client cannot draw the anchored
        # arrangement by default, and cannot offer the button that the display cap's whole rule
        # exists for — "otherwise the reader is told in words that their genome sits in #37 and has
        # no button to go back to it". Empty means the anchored genome has no gene here, which is a
        # different statement from there being no anchor; `is_anchored` separates them.
        "anchor": {
            "is_anchored": detail.is_anchored,
            "arrangement_ranks": list(detail.anchor_arrangement_ranks),
            # ⛔ `catalogue` | `projected` | null. A projected genome MATCHES an arrangement; a
            # catalogue genome is COUNTED in one. The page must not say "your genome is in
            # arrangement #3" about a genome that is in no arrangement at all, and the kind is what
            # stops it. `arrangement_ranks` stays empty for a projected genome — its relation is in
            # `projected_copies` instead, so a client that ignores the kind cannot accidentally
            # render a projected genome as a member.
            "kind": detail.anchor_kind,
            "projected_copies": [
                serialise_projected_placement(placement, detail) for placement in detail.projected_placements
            ],
        },
        "offsets": serialise_offset_occupants(detail),
        "intergenic_gaps": [
            serialise_intergenic_gap(gap, labels) for gap in detail.intergenic_gaps
        ],
        # ⭐ Every Pfam family the response MENTIONS, resolved here rather than by shipping the
        # 833 kB vendored table to the browser. Keyed VERSION-STRIPPED — a chip holding `PF00126.29`
        # must cut at the dot before looking itself up, or it silently misses and falls back to a
        # bare accession, which reads as "this family has no name" rather than as a failed join.
        # ⚠ An empty string is what the vendored table HAD, not a missing key: the reader pads short
        # rows. A caller must test membership, never truthiness of `short_name`.
        "pfam_reference": {
            accession: {
                "short_name": family.short_name,
                "description": family.description,
                # ⭐ Preferred over the Pfam entry for the link: it is the integrated record, and the
                # page a reader following a domain actually wants. "" where there is none.
                "interpro_accession": family.interpro_accession,
                "interpro_name": family.interpro_name,
                # ⚠ "" is CLANLESS and is not an identity — only ~46 % of families are in a clan, so
                # treating "" as a shared clan would make any two clanless families look like the
                # same superfamily.
                "clan_accession": family.clan_accession,
                "clan_name": family.clan_name,
            }
            for accession, family in sorted(detail.pfam_families.items())
        },
        # ⭐ The fan-out, answered in this response rather than in 15–303 more.
        "neighbour_display_rows": [
            {
                "label": row.node_label,
                # ⛔⛔ The row's OTHER address, and the map needs exactly this one. A slot code
                # carries `catalogue_ordinal * 2 + strand` and `nearest_locus_ordinals` is in the
                # same space, while `locus_offset_occupant.neighbour_locus_id` is a surrogate id.
                # Both are small integers over the same range, so looking one up in the other's
                # index draws one locus where another belongs — with a page that still looks right.
                # The client must resolve the map's nearest loci through THIS field.
                "catalogue_ordinal": row.catalogue_ordinal,
                # ⭐ For the MAP legend, not the track — the locus number is not what tells a reader
                # whether a neighbour belongs there; the product is.
                "best_product": row.best_product,
                # ⛔ `within_medoid_distance` and `map_position` went on 2026-09-24 with the map
                # they served — the first was its ring radius, the second its sprite position. A
                # neighbour row carries what NAMES a locus; the similarity relating it to the focal
                # one is on `similarity[rep].nearest_loci`, because it is a property of the PAIR.
                "display_name": row.display_name,
                "display_name_source": row.display_name_source,
                "genome_count": row.member_genome_count,
                "median_gene_length_nt": row.median_gene_length_nt,
                "prevalence_band": row.prevalence_band,
            }
            for row in sorted(neighbours.all_rows(), key=lambda row: row.node_label)
        ],
    }


def serialise_gene_sequence(row) -> dict:
    """One gene copy for the Sequence tab — where it is, what it reads, and what is missing.

    ⭐ **The stored values travel beside the sliced ones.** `gc_percent` and `protein_length_aa` are
    columns, written at ingest from the same coordinates; sending both lets a reader (and a test)
    see them agree rather than taking one on trust. They are separate fields for that reason and
    must not be collapsed into one.

    ⛔ **Truncation is reported per flank, not implied by a short string.** These are draft
    assemblies — mean contig 14.4 kb, about 14 genes — so a gene within 100 bp of a contig end is
    common rather than exotic, and a flank that is short without saying so reads as a complete one.
    """
    sequence = row.sequence
    return {
        # ⭐ `n of m` — rho > 1 puts one genome in a locus twice, and that is something to show
        # rather than something to average away.
        "copy_ordinal": row.copy_ordinal,
        "copy_count": row.copy_count,
        "flat_index": row.flat_index,
        # ⛔ The contig's NAME, never `contig_index + 1`: the index enumerates contigs that HAVE a
        # CDS, so index 269 is `contig00324` on SAMEA103923484.
        "contig_name": row.contig_name,
        "seqid": row.seqid,
        "start_position": row.start_position,
        "end_position": row.end_position,
        "strand": row.strand,
        # ⚠ Whether a strand parquet existed at all. A forced `+` must never read as an observation.
        "strand_is_observed": row.strand_is_observed,
        "gff_phase": row.gff_phase,
        "is_five_prime_partial": row.is_five_prime_partial,
        # Inclusive of the stop codon, so the protein is `length_nt / 3 - 1` residues.
        "length_nt": row.length_nt,
        "protein_length_aa": row.stored_protein_length_aa,
        "gene_symbol": row.bakta_gene_symbol,
        "product": row.bakta_product,
        # ⚠ Genome-private (`AAOCBP_22210`) and never a gene NAME.
        "locus_tag": row.locus_tag,
        "coding_sequence": sequence.coding_sequence,
        "protein_sequence": sequence.protein_sequence,
        # ⛔ In the GENE's reading direction, not the contig's. On a minus-strand gene the upstream
        # flank sits at HIGHER contig coordinates and is reverse-complemented with it. Slicing
        # `start - flank` unconditionally returns the DOWNSTREAM flank for about half of all genes
        # and looks entirely plausible on screen.
        "upstream_flank_sequence": sequence.upstream_flank_sequence,
        "downstream_flank_sequence": sequence.downstream_flank_sequence,
        "upstream_flank_is_truncated_by_contig_end": sequence.upstream_flank_is_truncated_by_contig_end,
        "downstream_flank_is_truncated_by_contig_end": sequence.downstream_flank_is_truncated_by_contig_end,
        # ⭐ Where each flank CAME FROM, in contig coordinates — so the page can say "contig
        # 1,388–1,487, reverse-complemented" without re-implementing the strand rule to work out
        # which end it was. `null` where the flank is empty.
        "upstream_flank_span": list(sequence.upstream_flank_span) if sequence.upstream_flank_span else None,
        "downstream_flank_span": (
            list(sequence.downstream_flank_span) if sequence.downstream_flank_span else None
        ),
        "gc_percent": sequence.gc_percent,
        "stored_gc_percent": row.stored_gc_percent,
    }
