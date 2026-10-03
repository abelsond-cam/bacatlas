/**
 * The `/api/v1` contract, as TypeScript. **Mirrors `serialisers/` exactly** — a change on either
 * side is breaking, which is the analogue of *"a change here is breaking for `app.js`"*.
 *
 * ⚠ **`null` means NOT MEASURED and `0` means MEASURED ZERO, throughout.** Every nullable number
 * below is nullable for that reason and not for tidiness, so `v-if` on any of them is wrong: a
 * measured zero is falsy, and for gap variance it is the *majority* case — 80–90 % of gaps vary not
 * at all. White on the track has to mean "identical in every genome", never "small". Test
 * `!== null`.
 */

import type { NeighbourSlot } from "@/lib/slotSpaces";

export type { NeighbourSlot };

/** Which of the two embeddings a geometric claim is made in. */
export type Representation = "esm" | "bacformer";

/**
 * How common a locus is across the collection — **five** bands, `prevalence.categorise`'s own.
 *
 * ⛔ `rare` was missing from this union while 9,877 of 33,201 loci — **30 %** — carry it. TypeScript
 * cannot catch that: the value arrives as a string and only a lookup keyed on the band, or an
 * exhaustive switch, ever notices, and then as `undefined` at render time rather than as an error.
 * `apiContract.test.ts` now asserts the union covers every band the real responses carry.
 *
 * ⚠ The order below is the order the bands are *drawn* in, commonest first. `PREVALENCE_BANDS` in
 * `lib/prevalence` is the runtime copy; the two are pinned to each other by a test.
 */
export type PrevalenceBand = "core" | "soft_core" | "shell" | "cloud" | "rare";

/** The three GO namespaces, spelled as the API spells them. */
export type GeneOntologyNamespace =
  | "molecular_function"
  | "biological_process"
  | "cellular_component";

/**
 * ⛔ `no_coverage` is a VALUE, not an absence. Fewer than two annotated members is neither
 * agreement nor disagreement, and counting it as either invents a finding.
 */
/**
 * ⛔⛔ **The Pfam ladder, not a yes/no** — and this type said `"agree" | "disagree" | "no_coverage"`
 * until the function tab was built against it.
 *
 * The server stores and sends `GeneOntologyAgreementVerdict`, six values, because GO agreement means
 * the same thing Pfam agreement does and one vocabulary serves both. Nothing caught the error: no
 * component read `go_verdicts` until now, and when one did every chip would have rendered blank —
 * `VERDICTS["single"]` is `undefined` in a table keyed on `"agree"` — with the page showing a locus
 * whose members agree as a locus with no verdict at all.
 *
 * ⚠ `no_coverage` is a VALUE, never a null: fewer than two annotated members is neither agreement
 * nor disagreement and must not be counted as either.
 */
export type GoVerdict =
  | "no_coverage"
  | "single"
  | "same_domains"
  | "nested"
  | "overlapping"
  | "disjoint";

/**
 * ⛔ **The annotation vocabulary, pinned — `backend/bacatlas_backend/models/enumerations.py`'s
 * `AnnotationKind`, all seven members.**
 *
 * It exists because `annotations` was `Record<string, …>`, and a bare string index means a key that
 * matches nothing type-checks, renders, and shows an empty card. `FunctionTab` asked for
 * `enzyme_commission`; the backend has always emitted `ec_number`; the string appeared exactly once
 * in the whole front end and in no test. **4,544 loci across the two live catalogues carry an EC
 * number and not one of them displayed it** — 3,348 drew no EC/KEGG card at all and 1,196 drew it
 * with the KEGG row alone (measured 2026-10-03). Six of the seven keys were right, and the compiler
 * could not tell them apart from the seventh.
 *
 * ⚠ So the fix is the TYPE, not the one string: `entriesOf` now takes an `AnnotationKind`, and the
 * next typo is a build failure rather than a blank card. `[[build-an-instrument-not-another-fix]]`.
 */
export const ANNOTATION_KINDS = [
  "gene_symbol",
  "protein_product",
  "pfam_architecture",
  "cog_orthogroup",
  "gene_ontology_slim",
  "ec_number",
  "kegg_orthology",
] as const;

export type AnnotationKind = (typeof ANNOTATION_KINDS)[number];

/** What the response carries per vocabulary. `Partial`: an absent kind is absent, not empty. */
export type AnnotationsByKind = Readonly<
  Partial<Record<AnnotationKind, readonly AnnotationEntry[]>>
>;

export interface AnnotationEntry {
  readonly rank: number;
  readonly term: string;
  /** `null` where the vocabulary has no name to give — KEGG ids are present and never named. */
  readonly name: string | null;
  readonly gene_count: number;
  readonly gene_ontology_namespace?: GeneOntologyNamespace;
}

export interface UnirefFamily {
  readonly rank: number;
  readonly uniref50_accession: string;
  readonly gene_count: number;
  readonly modal_product: string | null;
  readonly modal_architecture: string | null;
  /**
   * ⚠ **Nullable in the schema**, though populated for all 38,672 crosstab rows in the two loaded
   * catalogues. Typed as the column allows rather than as this data happens to be: a contract that
   * promises non-null where the schema permits null is a latent crash one ingest away, and the card
   * must test `!== null` — a measured **zero** here is the loudest case, a family where no gene is
   * annotated at all.
   */
  readonly pfam_annotated_gene_count: number | null;
  /** ⚠ `null` for 20,600 of 38,672 rows — over half the families carry no modal symbol. */
  readonly modal_symbol: string | null;
  /** ⚠ A COUNT, not a list. The card shows the modal plus "+N" and cannot name the others. */
  readonly distinct_symbol_count: number | null;
  /**
   * How many of THIS family's genes carry any gene name — the modal symbol's DENOMINATOR.
   *
   * ⛔ Without it a family of 9 genes of which **2** are named `rfbX` reports
   * `distinct_symbol_count = 1`, draws no "+N", and reads as nine genes agreeing on a name (kp
   * locus 3992). ⚠ Nullable like `pfam_annotated_gene_count` and for the same reason: `null` is
   * *not measured*, `0` is *measured and none named*, and the card must test `!== null`.
   */
  readonly named_gene_count: number | null;
}

export interface Arrangement {
  readonly rank: number;
  readonly gene_count: number;
  readonly genome_count: number;
  /**
   * ⚠ INTRINSIC to the arrangement — what the badge and the "recorded at" footnote describe. It is
   * NOT display mirroring, and composing the two is `displayMirrorApplied`'s job.
   */
  readonly is_recorded_reverse_complement: boolean;
  /** Ten, always, in RECORDED order. */
  readonly slots: readonly NeighbourSlot[];
}

export interface OffsetOccupant {
  readonly rank: number;
  readonly locus: string | null;
  readonly gene_count: number;
  readonly same_strand_gene_count: number;
}

export interface OffsetMarginal {
  readonly signed_offset: number;
  /** The honest denominator: members with a gene at this position at all. */
  readonly observed_member_count: number;
  /** ⛔ A DISPLAY cut — the occupants past the exported top-N. */
  readonly observed_not_listed: number;
  /** ⛔ MISSING DATA — members with no gene here, typically a contig end. A different sentence. */
  readonly members_without_an_observation: number;
  readonly occupants: readonly OffsetOccupant[];
}

export interface IntergenicGap {
  /** ⛔ Keyed by LABELS, never by index — the published page's index lookup missed 7,379 of 22,838. */
  readonly flanking_loci: readonly [string | null, string | null];
  readonly observed_genome_count: number;
  /** ⛔ SIGNED. Negative means the two genes OVERLAP; 18.8 % of adjacencies do. */
  readonly median_signed_length_nt: number | null;
  readonly quartile1_signed_length_nt: number | null;
  readonly quartile3_signed_length_nt: number | null;
  readonly minimum_signed_length_nt: number | null;
  readonly maximum_signed_length_nt: number | null;
  /** ⛔ `null` = not measured, `0.0` = every genome agrees. */
  readonly length_variance_score: number | null;
  readonly modal_length_nt: number | null;
  readonly distinct_named_feature_count: number;
  /** ⚠ Only `mn == mx` certifies this. `q1 == q3` is the weaker middle-half claim. */
  readonly every_genome_agrees: boolean;
}

/**
 * One of the five nearest OTHER loci in one representation.
 *
 * ⚠ **The list is RAGGED** — a locus whose shortlist held fewer than five simply has fewer rows, so
 * a rank is never a position in an array. Rows rather than the retired `nearest_locus_ordinals`
 * array because each carries its own similarity, and an array would have needed a sentinel: a
 * sentinel read as an index names a real locus that looks entirely plausible.
 */
export interface NearestLocus {
  /** ⚠ 1-based, as the artifact writes it. Rank 1 is BY CONSTRUCTION `nearest_similarity`. */
  readonly rank: number;
  /**
   * ⛔⛔ A **catalogue ordinal**, resolved through `neighbour_display_rows` — the same key space an
   * arrangement slot code carries, and *not* the surrogate `neighbour_locus_id` the offset
   * occupants use. Both are small integers over the same range, so resolving through the other
   * index names one locus where another belongs, on a page that still looks entirely right
   * (`2b99bb4`). The server DROPS a neighbour the fan-out did not resolve rather than sending a
   * null address, so every ordinal here has a row.
   */
  readonly catalogue_ordinal: number;
  readonly cross_similarity: number;
}

/**
 * One locus's **set-to-set** similarity in one representation — three pairs and a share.
 *
 * ⛔ **This REPLACED `LocusGeometry`; it is a different measurement, not a renaming.** The medoid
 * geometry reduced a locus to ONE member and then measured that single point: `within` was its
 * members' distance to that gene, `nearest` that gene's distance to another locus's medoid. Every
 * number here is a median over whole gene SETS, or is anchored on the one member least attached to
 * its locus. A client reading one as the other would be plausible and wrong.
 *
 * ⚠ **`separation` and `margin` are NOT here although the card prints both.** Each is the
 * difference of two fields that are, and the client subtracts — a separately-rounded difference
 * drifting from the two rounded numbers printed beside it is the exact class of quiet disagreement
 * this card was rebuilt to end.
 */
export interface LocusSimilarity {
  /**
   * ⛔ `null` on a SINGLETON, which has no pair inside its locus — 5,427 of *E. coli*'s 17,531
   * loci. A `0.0` here would claim its members are unrelated to each other.
   */
  readonly within_similarity: number | null;
  /**
   * The highest such median against another locus. ⚠ **Present for a singleton**, where the three
   * fields below are not: that question is well posed for one gene. So this is the one field a card
   * must NOT test to decide whether a locus is measurable.
   */
  readonly nearest_similarity: number | null;
  /**
   * The member least attached to its own locus: its nearest neighbour INSIDE the locus, then that
   * same gene's nearest gene anywhere else. ⭐ The point is invariant where a median moves with the
   * set — and the clustering joins points, not medians.
   *
   * ⚠ `weak_own_similarity` implies the own fraction by ARITHMETIC: if the weakest member's nearest
   * gene is a stranger then the fraction cannot be 1. The two agreeing is never evidence.
   */
  readonly weak_own_similarity: number | null;
  readonly weak_other_similarity: number | null;
  /**
   * The share of members whose nearest gene in the whole species is another member of this locus.
   *
   * ⛔ **A SHARE, not a cosine.** It is never drawn against the random gene-pair floor, and never
   * formatted at 0 dp: `0.9999` rounds to "100%", and below 1 is the entire point of the number.
   */
  readonly own_neighbour_fraction: number | null;
  /**
   * ⭐ One precomputed MIDRANK per VIEW, over MEASURABLE loci only — never over the catalogue. The
   * three are not three ranks of one thing: only 18–34 % of flagged loci are flagged by all three.
   * ⚠ `null` reads "not measurable" and must never render as `0.000`.
   */
  readonly separation_percentile: number | null;
  readonly weak_margin_percentile: number | null;
  /**
   * ⛔ At an own fraction of exactly 1.0 this is a midrank inside a tie block covering **88.5 %** of
   * the catalogue, so it reads "p53" — *better than half the catalogue* — when it means *tied with
   * nearly all of it*. Served because below 1.0 the whole tie block is above it and it then means
   * what it looks like; **the card is what must decline to print it at 1.0.**
   */
  readonly own_fraction_percentile: number | null;
  readonly nearest_loci: readonly NearestLocus[];
}

/**
 * One Pfam-A family, resolved server-side from the vendored public reference.
 *
 * ⚠ **Every field is a string and `""` means the table HAD no value** — never that the key is
 * missing. The vendored reader pads short rows rather than skipping them, so a caller tests
 * membership of the map, not the truthiness of `short_name`.
 */
export interface PfamFamilyReference {
  /** `Sigma70_r2` — what the chip says when it has one. */
  readonly short_name: string;
  readonly description: string;
  /**
   * ⭐ Preferred over the Pfam entry for the link: it is the integrated record, and the page a
   * reader following a domain actually wants. `""` where there is no integrated entry.
   */
  readonly interpro_accession: string;
  readonly interpro_name: string;
  /**
   * ⚠ `""` is CLANLESS and is **not** an identity. Only ~46 % of families are in a clan, so
   * treating `""` as a shared clan would make any two clanless families look like the same
   * superfamily — the common case, not the corner.
   */
  readonly clan_accession: string;
  readonly clan_name: string;
}

export interface NeighbourDisplayRow {
  readonly label: string;
  /**
   * ⛔⛔ **The row's OTHER address, and `similarity[rep].nearest_loci` needs exactly this one.** An
   * arrangement slot code carries `catalogue_ordinal * 2 + strand` and a nearest-locus row is in the
   * same space — while the marginal occupants are addressed by `label`. Both are small integers
   * over the same range, so resolving the nearest loci through anything else names one locus where
   * another belongs, on a page that still looks entirely right. That exact merge already cost this
   * project once (`2b99bb4`).
   */
  readonly catalogue_ordinal: number;
  readonly display_name: string;
  /**
   * ⭐ For the NEAREST-LOCI list, not the track. *"The locus NUMBER is not what tells you whether a
   * neighbour belongs here — the product is."* The track has no room for it and does not ask.
   */
  readonly best_product: string | null;
  /**
   * ⛔ `map_position` and `within_medoid_distance` went on 2026-09-24 with the map they served — the
   * first was a sprite position, the second a ring radius. A neighbour row now carries only what
   * NAMES a locus; the similarity relating it to the focal one lives on `similarity[rep].nearest_loci`,
   * because it is a property of the PAIR rather than of this row.
   */
  readonly display_name_source: string;
  readonly genome_count: number;
  readonly median_gene_length_nt: number | null;
  readonly prevalence_band: PrevalenceBand;
}

export interface Locus {
  readonly label: string;
  readonly catalogue_ordinal: number;
  readonly display_name: string;
  readonly display_name_source: string;
  readonly display_name_source_accession: string | null;
  readonly best_product: string | null;
  readonly bakta_gene_symbol: string | null;
  readonly gene_count: number;
  readonly genome_count: number;
  readonly named_gene_count: number;
  readonly prevalence_band: PrevalenceBand;
  readonly median_gene_length_nt: number | null;
  readonly gene_length_interquartile_range_nt: number | null;
  readonly uniref50: {
    readonly family_count: number;
    /** ⚠ What the card CLAIMS with, while the list shows every family. */
    readonly major_family_count: number;
    readonly labelled_gene_count: number;
    readonly impurity: number | null;
    readonly coverage: number | null;
  };
  readonly pfam: {
    /** ⛔ `gene_count − this` is MISSING COVERAGE, not a competing architecture. */
    readonly annotated_gene_count: number;
    readonly architecture_count: number;
    /** ⛔ The audit's own verdict, READ and never re-derived — the page once disagreed on 2 of 22,624. */
    readonly concordance_class: string | null;
  };
  readonly evidence: {
    readonly syntenic_a5: number | null;
    readonly collapse_tier: string | null;
    readonly collapse_bucket: string | null;
    /** ⚠ "the identity at which the 50 LONGEST members first group" — weaker than the whole locus. */
    readonly resolved_threshold: number | null;
    readonly resolved_threshold_is_capped_at_50_members: boolean;
  };
  /**
   * ⛔ **`null` for a representation with no row at all**, which is not the same absence as a row
   * whose `within_similarity` is null — the second is a singleton, measured and found to have no
   * pair. Both must read as sentences rather than as numbers.
   */
  readonly similarity: Readonly<Record<Representation, LocusSimilarity | null>>;
  readonly interest_score: number | null;
}

export interface LocusDetailResponse {
  readonly locus: Locus;
  readonly annotations: AnnotationsByKind;
  readonly uniref50_families: readonly UnirefFamily[];
  readonly arrangements: {
    readonly listed: readonly Arrangement[];
    /** ⛔ Never moved by any display cap, and never conflated with the number listed. */
    readonly total: number;
    /** ⛔ Arrangements the cap left out — a COUNT of arrangements. */
    readonly arrangements_not_listed: number;
    /**
     * ⛔ Members sitting INSIDE those uncapped arrangements. A different remainder from the one
     * below, and the API's own creation: the published payload is uncapped, so this category did
     * not exist there.
     */
    readonly members_in_arrangements_not_listed: number;
    /**
     * ⛔ Members with no recorded neighbourhood **at all** — no coordinates for the gene, so no
     * window. ⚠ Folding this together with the field above tells a reader that 15,912 *E. coli*
     * genes have no coordinates when they simply sit past the display cap.
     */
    readonly members_without_a_neighbourhood: number;
    /**
     * ⛔ Whether EVERY genome present reaches an arrangement — a **genome** question, and the only
     * thing that settles which of two sentences the anchor line may say. With no anchored ranks,
     * `true` means *"has no gene at this locus"* and `false` means *"has no recorded neighbourhood
     * at this locus"*; they are different claims and only one is ever true. Measured: **6.26 % of
     * ecoli loci and 3.69 % of kp loci** are incomplete (worst case 64 genomes of 100).
     *
     * ⚠ Not derivable from the two gene remainders above. A genome at ρ > 1 has two genes here, and
     * one of them losing its window leaves the genome fully present in an arrangement.
     */
    readonly membership_is_complete: boolean;
  };
  /**
   * ⛔ Which arrangement RANKS the anchored genome carries — a **list**, because rho > 1 puts one
   * genome in two arrangements at one locus and there is no uniqueness constraint on
   * (locus, genome) anywhere.
   *
   * ⚠ An empty list carries two different facts and `is_anchored` is what separates them: *"your
   * genome has no gene at this locus"* against *"you have not anchored one"*. They are different
   * sentences on the page.
   */
  readonly anchor: {
    readonly is_anchored: boolean;
    readonly arrangement_ranks: readonly number[];
    /**
     * ⛔ `catalogue` | `projected` | `null`. A catalogue genome is COUNTED in the arrangements it
     * occupies; a projected genome MATCHES one without being in it, so `arrangement_ranks` stays
     * empty for it and its relation is in `projected_copies`. A client that ignores the kind cannot
     * accidentally render a projected genome as a member — but one that reads `arrangement_ranks`
     * alone would say "no gene at this locus" about a gene that is plainly placed there.
     */
    readonly kind: "catalogue" | "projected" | null;
    readonly projected_copies: readonly ProjectedCopy[];
  };
  readonly offsets: readonly OffsetMarginal[];
  readonly intergenic_gaps: readonly IntergenicGap[];
  /**
   * ⭐ Every Pfam family this response MENTIONS, keyed by **version-stripped** accession. The same
   * move as `neighbour_display_rows`: a bounded block in one round trip, instead of the page
   * carrying an 833 kB vendored reference to render a chip.
   *
   * ⛔ A chip holding `PF00126.29` must cut at the dot before looking itself up, or it silently
   * misses and falls back to a bare accession — which reads as *"this family has no name"* rather
   * than as a failed join.
   */
  readonly pfam_reference: Readonly<Record<string, PfamFamilyReference>>;
  /** ⭐ The 15–303-locus fan-out, answered in THIS response rather than in that many more. */
  readonly neighbour_display_rows: readonly NeighbourDisplayRow[];
  readonly resolved_neighbour_count: number;
}

/**
 * One gene of a PROJECTED genome at this locus — placed after the model was built.
 *
 * ⛔ **`agreeing_neighbours` is meaningless without `available_neighbours`.** A locus with *m*
 * modelled genes can supply at most min(n, m) of the n checkers, so the raw count tracks the
 * locus's SIZE as much as the evidence: measured over the first ten placed genomes, `is_contested`
 * runs at 0.39 % where the locus has ten genes to offer and 42.6 % where it has fewer, reaching
 * 96.3 % at singleton loci. Render them together or not at all.
 */
export interface ProjectedCopy {
  readonly flat_index: number;
  readonly copy_ordinal: number;
  readonly copies_at_locus: number;
  readonly nearest_cosine: number | null;
  readonly agreeing_neighbours: number;
  readonly available_neighbours: number;
  readonly placed_summed_cosine: number | null;
  readonly is_contested: boolean;
  readonly runner_up_summed_cosine: number | null;
  /**
   * ⚠ `false` means the gene is ALONE ON ITS CONTIG and has no ±5 window at all — not that its
   * window matched nothing, which is `matched_arrangement_rank === null` with this `true`. Three
   * outcomes, and collapsing two of them reports a neighbourhood comparison for a gene with no
   * neighbours.
   */
  readonly has_a_neighbourhood: boolean;
  readonly matched_arrangement_rank: number | null;
  readonly matched_arrangement_genome_count: number | null;
}

/** A genome placed on this catalogue after the model was built — never a member of it. */
export interface ProjectedGenomeEntry {
  readonly sample_id: string;
  readonly rule: string;
  readonly neighbours_searched: number;
  readonly neighbours_reported: number;
  readonly gene_count: number;
  readonly placed_gene_count: number;
  /** ⛔ Four counts that are not the same number — see the backend serialiser. */
  readonly genes_without_a_vector: number;
  readonly genes_without_a_neighbourhood: number;
  readonly contested_gene_count: number;
  readonly distinct_locus_count: number;
  readonly multi_copy_locus_count: number;
  readonly window_matched_gene_count: number;
  /** ⚠ The distribution, not a mean: there is no novelty threshold, so the spread is the signal. */
  readonly nearest_cosine: {
    readonly median: number | null;
    readonly fifth_percentile: number | null;
    readonly minimum: number | null;
  };
  readonly nuna_git_sha: string | null;
}

export interface ProjectedGenomesResponse {
  readonly species_key: string;
  readonly genomes: readonly ProjectedGenomeEntry[];
  /** The sentence the page must repeat wherever these genomes appear. */
  readonly caveat: string;
}

export interface ArrangementPageResponse {
  readonly arrangements: readonly Arrangement[];
  readonly offset: number;
  readonly total: number;
}

/**
 * ⛔ **The vocabularies the inference ladder covers, pinned to the backend's `SUPPORTED_KINDS`.**
 * A string here would let a typo compile — the defect that hid 4,544 loci's EC numbers behind the
 * key `enzyme_commission` for as long as `annotations` was a `Record<string, …>`.
 */
export const INFERENCE_KINDS = ["cog_orthogroup", "ec_number"] as const;
export type InferenceKind = (typeof INFERENCE_KINDS)[number];

/** One (similarity tier, depth) cell of the measured ladder. */
export interface CalibrationCell {
  readonly tier: string;
  readonly level: number;
  readonly level_label: string;
  readonly pairs: number;
  /** ⛔ `null` below the backend's pair floor — a cell with too few pairs has a count, not a rate. */
  readonly agreement: number | null;
  readonly interval_low: number | null;
  readonly interval_high: number | null;
  readonly chance: number;
  /** Observed ÷ chance. ⛔ NOT an odds ratio — see `annotation_transfer.Cell.lift`. */
  readonly lift: number | null;
}

export interface CalibrationLadder {
  readonly annotation_kind: InferenceKind;
  readonly representation: Representation;
  readonly annotated_locus_count: number;
  readonly locus_count: number;
  readonly min_pairs: number;
  readonly cells: readonly CalibrationCell[];
}

/** What this node's OWN call rests on. */
export interface OwnSupport {
  readonly term: string;
  readonly name: string | null;
  /** How many member genes carry the modal term — the vote, not the node size. */
  readonly gene_count: number;
  readonly annotated_gene_count: number;
  readonly member_gene_count: number;
  /**
   * ⛔ `false` where exactly ONE gene carries a call. Such a node is unanimous *by construction*,
   * which is the absence of evidence rather than evidence — 518 of 4,993 *E. coli* COG nodes, `gumC`
   * among them. A card that rendered it like any other would claim a check that never happened.
   */
  readonly checkable: boolean;
}

/** One step of the outward walk, including the ranks that carried nothing. */
export interface WalkStep {
  readonly rank: number;
  readonly cosine: number | null;
  readonly tier: string;
  readonly node_label: string;
  readonly catalogue_ordinal: number;
  readonly display_name: string | null;
  readonly carries_annotation: boolean;
}

export interface TransferCandidate {
  readonly rank: number;
  readonly cosine: number | null;
  readonly tier: string;
  readonly donor: OwnSupport & {
    readonly node_label: string;
    readonly catalogue_ordinal: number;
    readonly display_name: string | null;
  };
  /** `null` where the tier is too remote to call, or the donor states nothing this deep. */
  readonly level: number | null;
  readonly value: readonly string[] | null;
  readonly calibration: CalibrationCell | null;
}

export interface VocabularyInference {
  readonly annotation_kind: InferenceKind;
  readonly own: OwnSupport | null;
  /** ⛔ Empty where the node has its own call — a suggestion would compete with a real annotation. */
  readonly walk: readonly WalkStep[];
  readonly candidate: TransferCandidate | null;
}

export interface FunctionResponse {
  readonly annotations: AnnotationsByKind;
  readonly inference: {
    readonly representation: Representation;
    readonly vocabularies: readonly VocabularyInference[];
  };
  readonly calibration: Readonly<Record<InferenceKind, CalibrationLadder>>;
  readonly coverage: {
    readonly gene_count: number;
    readonly cog_annotated_gene_count: number;
    /** ⛔ A COUNT, never a relation — two COG ids in one locus is ordinary above family level. */
    readonly cog_distinct_id_count: number;
    readonly modal_cog_categories: readonly string[] | null;
    readonly ec_annotated_gene_count: number;
    readonly kegg_annotated_gene_count: number;
    readonly go_annotated_gene_count: Readonly<Record<GeneOntologyNamespace, number>>;
  };
  readonly go_verdicts: Readonly<Record<GeneOntologyNamespace, GoVerdict | null>>;
}

export interface ModelStep {
  readonly ordinal: number;
  readonly name: string;
  readonly representation: Representation | null;
  readonly gamma: number | null;
  readonly rho_rule: string | null;
  readonly rho_ceiling: number | null;
  readonly uses_exclusivity: boolean;
  readonly stage_name: string | null;
  readonly stage_detail: string | null;
}

/**
 * What makes a cosine on the similarity card mean anything: the **random GENE-pair floor** for one
 * representation, and its spread.
 *
 * ⛔ **This REPLACED `map_projections`, it did not rename it.** That carried a UMAP over every
 * locus's MEDOID, its extent and a whole-catalogue sprite, and a null sampled over random pairs of
 * **medoids** — all of it describing a construction that reduced a locus to one member. On the
 * published *E. coli* catalogue the medoid null and this gene-pair floor sit at **0.0651** and
 * **0.0587**: close enough to look interchangeable, and not interchangeable.
 *
 * ⚠ Without a floor a cosine has no scale at all — ESM's is ~0.742 and Bacformer's ~0.059, so the
 * same 0.41 reads oppositely in the two.
 */
export interface SimilarityBaseline {
  readonly representation: Representation;
  /**
   * `raw` or `centred`. ⛔ Raw is what is published; a centred number captioned as a raw one is
   * indistinguishable from the real thing — ESM's floor moves 0.7417 → ~0.005.
   */
  readonly form: string;
  readonly floor_median: number | null;
  /**
   * ⭐ The floor's spread, so the strip can draw a box rather than a bare tick: a median alone
   * cannot say whether a locus's 0.41 sits far outside random or inside its shoulder. `null` where
   * the run's audit JSON was not beside its CSV — then the strip falls back to the median alone.
   */
  readonly floor_p25: number | null;
  readonly floor_p75: number | null;
  readonly floor_p99: number | null;
  /** ⭐ The other half of "p12 of 12,104 loci" — the denominator every midrank was taken over. */
  readonly measurable_locus_count: number | null;
  /** The kNN width the shortlist and the point measures were computed at. */
  readonly neighbour_knn_k: number | null;
}

/**
 * ⛔ **A catalogue key, distinguishable by the compiler from a species key.**
 *
 * Both are strings, both name something in the URL, and they were interchangeable right up to the
 * moment a species held two catalogues. Passing `ecoli` where `ecoli-nuna5` belongs does not throw:
 * it silently serves the species' DEFAULT clustering under a page that says it is showing another
 * one, and it renders perfectly. The type system could not see it, because nothing about the shape
 * of the string changed — which is why this is a brand rather than a comment.
 *
 * ⚠ Mint it only with {@link asCatalogueKey}, and only from something that really is one: a
 * `catalogue_key` the server sent, an entry in `GET /catalogues`, or a URL parameter about to be
 * handed straight back to the server for resolution. Never from `species.key`.
 */
export type CatalogueKey = string & { readonly __catalogue: unique symbol };

/** Mint a {@link CatalogueKey}. The single place the brand is applied, so it is greppable. */
export function asCatalogueKey(key: string): CatalogueKey {
  return key as CatalogueKey;
}

/**
 * One row of `GET /api/v1/catalogues` — everything the model picker needs.
 *
 * ⚠ A species appears ONCE PER CATALOGUE it holds, so E. coli on nuna4 and nuna5 is two rows.
 * `is_default` marks the one a bare `?species=` resolves to; `key` is what an address should carry
 * to pin this clustering.
 *
 * ⛔ Only catalogues the server OFFERS are listed. One that is loaded but staged is absent from
 * here and still reachable by key, so never gate a load on membership of this list.
 */
export interface CatalogueEntry {
  readonly key: CatalogueKey;
  readonly species: { readonly key: string; readonly scientific_name: string };
  /** ⚠ Nullable: a pangenome need not be tied to a registry model. */
  readonly model: {
    readonly key: string;
    readonly label: string | null;
    readonly step_count: number | null;
    readonly exclusivity_form: string;
  } | null;
  readonly is_default: boolean;
  readonly genome_count: number;
  readonly gene_count: number;
  readonly locus_count: number;
  readonly run_id: string;
}

export interface CataloguesResponse {
  readonly catalogues: readonly CatalogueEntry[];
}

export interface SpeciesCatalogueResponse {
  readonly species: { readonly key: string; readonly scientific_name: string };
  readonly pangenome: {
    /**
     * ⭐ How this catalogue is addressed — `ecoli-nuna4`, `ecoli-nuna5`.
     *
     * A page that followed a bare `?species=` needs this to name what it was actually given, so the
     * address can be pinned to the clustering being read. ⛔ It cannot be rebuilt from
     * `species.key` + `model.key`: that pair is not unique, and a catalogue may carry a deliberate
     * key such as `ecoli-sensitive`.
     */
    readonly catalogue_key: CatalogueKey;
    readonly run_id: string;
    readonly genome_count: number;
    readonly gene_count: number;
    readonly locus_count: number;
    readonly built_at: string | null;
    /**
     * When this stack INGESTED the catalogue. ⛔ Not the same claim as `built_at`, which is when the
     * model was exported — and which is NULL on every catalogue loaded so far, because nothing in
     * `ingest` writes it. This one is populated, so it is what the footer can actually print.
     */
    readonly ingested_at: string | null;
    readonly git_sha: string | null;
    readonly exclusivity_form: string;
    /** ⚠ An absent section is NAMED, so an omission is never mistaken for a measured zero. */
    readonly omitted_sections: Readonly<Record<string, string>>;
  };
  readonly model: {
    readonly key: string;
    readonly label: string | null;
    readonly exclusivity_form: string;
    readonly knn_k: number | null;
    readonly step_count: number | null;
  } | null;
  readonly steps: readonly ModelStep[];
  /** The footer's `<dl>`, one `[label, value]` pair per row, printed verbatim (`render_page._provenance`). */
  readonly provenance_rows: readonly (readonly [string, string])[];
  /** Loci per band. ⛔ Every band is present — a band with no loci is a measured `0`, never absent. */
  readonly prevalence_census: Readonly<Record<PrevalenceBand, number>>;
  /**
   * ⭐ The SAME catalogue partitioned by GENE. Divided by the genome count it is the census's
   * "per genome" line — a different partition, not a restatement: a third of the loci are singletons
   * while a typical genome carries about fifty singleton genes.
   */
  readonly prevalence_gene_census: Readonly<Record<PrevalenceBand, number>>;
  /** Read verbatim from the audit summary at ingest. Keys absent from an older summary are absent here. */
  readonly audit_headline: Readonly<Record<string, number | string | null>>;
  readonly similarity_baselines: readonly SimilarityBaseline[];
  readonly landing_locus: string | null;
  readonly example_loci: readonly string[];
  /** The example chips as drawn: name, and the UniRef50 families the chip quotes. */
  readonly example_locus_rows: readonly ExampleLocusRow[];
}

export interface ExampleLocusRow {
  readonly label: string;
  readonly display_name: string;
  readonly uniref50_family_count: number | null;
}

export interface SpeciesListResponse {
  readonly species: readonly {
    readonly key: string;
    readonly scientific_name: string;
    readonly ncbi_taxonomy_id: number | null;
    /** ⚠ A species with nothing published is LISTED with `false`, never filtered out. */
    readonly published: boolean;
    readonly genome_count: number | null;
    readonly gene_count: number | null;
    readonly locus_count: number | null;
    /**
     * The species' DEFAULT catalogue — what a bare `?species=` resolves to, and the key to
     * address once it has been followed. `null` when the species serves nothing.
     *
     * ⚠ Replaced a `model_label` that the server hardcoded to `null`. The full menu of
     * catalogues is `GET /catalogues`, not this.
     */
    readonly catalogue_key: CatalogueKey | null;
  }[];
}

export interface SearchHit {
  readonly label: string;
  readonly display_name: string;
  /** What the page prints under the name — the number never tells a reader whether a hit is theirs. */
  readonly best_product: string | null;
  readonly gene_count: number;
  readonly genome_count: number;
  readonly prevalence_band: PrevalenceBand;
  readonly rank_band: string | null;
}

export interface SearchResponse {
  readonly query: string;
  /**
   * ⚠ Which mode answered. A 1–2 character query cannot use the trigram index and searches less of
   * the haystack; a reader is entitled to know that rather than to conclude the catalogue is empty.
   */
  readonly mode: string;
  readonly truncated: boolean;
  readonly hits: readonly SearchHit[];
}

/**
 * One gene copy's bases — the Sequence tab's payload.
 *
 * ⭐ **~5 kB out, where the published page pulled 1.3 MB in.** It had no choice: GitHub Pages
 * applies `Range` to the *compressed* stream, so byte offsets return plausible wrong bytes with a
 * 206 and no error, and the browser had to fetch whole `.nseq` files and decode DNA itself. A server
 * slices the file Bakta wrote, so that custom format is retired rather than ported.
 */
export interface GeneSequence {
  /** ⭐ ρ > 1 puts one genome in a locus twice, so a copy is `n of m` and the page says so. */
  readonly copy_ordinal: number;
  readonly copy_count: number;
  readonly flat_index: number;
  /**
   * ⛔ The contig's own NAME, never `contig_index + 1` — the index enumerates contigs that HAVE a
   * CDS. Measured: the derived name is wrong on 7,696 of 26,878 contigs (28.6 %) and right on the
   * rest, so a page deriving it looks correct seven times in ten and names a real contig the gene is
   * not on for the other three.
   */
  readonly contig_name: string;
  readonly seqid: string;
  /** 1-based inclusive, and the end INCLUDES the stop codon. */
  readonly start_position: number;
  readonly end_position: number;
  readonly strand: "+" | "-";
  /** ⚠ Whether a strand parquet existed at all. A forced `+` must never read as an observation. */
  readonly strand_is_observed: boolean;
  readonly gff_phase: number;
  readonly is_five_prime_partial: boolean;
  readonly length_nt: number;
  readonly protein_length_aa: number | null;
  readonly gene_symbol: string | null;
  readonly product: string | null;
  /** ⚠ Genome-private (`AAOCBP_22210`) and never a gene NAME. */
  readonly locus_tag: string | null;
  readonly coding_sequence: string;
  readonly protein_sequence: string;
  /**
   * ⛔ In the GENE's reading direction, not the contig's. On a minus-strand gene the upstream flank
   * sits at HIGHER contig coordinates and is reverse-complemented with it.
   */
  readonly upstream_flank_sequence: string;
  readonly downstream_flank_sequence: string;
  readonly upstream_flank_is_truncated_by_contig_end: boolean;
  readonly downstream_flank_is_truncated_by_contig_end: boolean;
  /** Contig coordinates the flank came from, 1-based inclusive. `null` where the flank is empty. */
  readonly upstream_flank_span: readonly [number, number] | null;
  readonly downstream_flank_span: readonly [number, number] | null;
  readonly gc_percent: number;
  /** What ingest wrote from the same coordinates — sent so the two can be seen to agree. */
  readonly stored_gc_percent: number | null;
}

export interface GeneSequenceResponse {
  readonly genome: { readonly sample_id: string };
  readonly locus: { readonly label: string };
  /**
   * ⛔ An EMPTY list is the answer *"this genome has no gene at this locus"*, not a failure. A
   * request that failed is a `Failure`, and the two must never render alike.
   */
  readonly genes: readonly GeneSequence[];
}

/**
 * One genome the anchor control can offer.
 *
 * ⚠ **Two counts, because they are two facts.** `locus_count` is the loci where this genome has a
 * gene; `arrangement_locus_count` is the loci where it appears in some recorded neighbourhood — the
 * published page's `genomeCounts()`, which could only see the latter. They differ wherever a gene
 * got no window. Both are per LOCUS: a genome at ρ > 1 is counted once, or the number would exceed
 * the catalogue and mean nothing.
 */
export interface GenomeRow {
  readonly sample_id: string;
  readonly collection_genome_ordinal: number;
  readonly locus_count: number;
  readonly arrangement_locus_count: number;
}

/** `GET /species/{key}/genomes?q=` — replaces the resident `meta.genomes` + `anchorSearch`. */
export interface GenomeListResponse {
  readonly species_key: string;
  readonly query: string;
  readonly genome_count: number;
  /** How many match `q` before the limit — so a truncated list can say how much it left out. */
  readonly matched_genome_count: number;
  readonly truncated: boolean;
  readonly genomes: readonly GenomeRow[];
}

/** One row of the footer's residual lists — the evidence beside the name, to judge before clicking. */
export interface ResidualLocusRow {
  readonly label: string;
  readonly display_name: string;
  readonly prevalence_band: PrevalenceBand;
  readonly gene_count: number;
  readonly uniref50_family_count: number | null;
  /** `null` where Pfam could not judge the locus — not "no architectures". */
  readonly pfam_architecture_count: number | null;
  readonly syntenic_a5: number | null;
  /**
   * ⛔ **SIMILARITIES, and nothing is left to convert.** These were medoid DISTANCES and the footer
   * turned each into `1 − d`; they are now the same set-to-set numbers the locus card shows — a
   * median over every within-locus gene pair, and the highest such median against another locus. A
   * surviving `1 − d` would silently print the complement of a real similarity.
   */
  readonly esm_within_similarity: number | null;
  readonly esm_nearest_similarity: number | null;
}

/**
 * `GET /species/{key}/audit/residual-loci` — ⛔ TWO lists, never merged: grouped on context alone is
 * a statement about the evidence, a Pfam conflict is evidence against the merge, and neither is a
 * verdict. A locus can be in both, and then appears in both.
 */
export interface AuditResidualsResponse {
  readonly grouped_on_context_alone: readonly ResidualLocusRow[];
  readonly pfam_conflicts: readonly ResidualLocusRow[];
}
