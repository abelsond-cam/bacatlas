/**
 * The three views of a locus's set-to-set similarity, and the banding each one is read through.
 *
 * ⭐ **Three views, and they REPLACE each other rather than stacking** (David, 2026-09-24: *"a click
 * box option for median (default to show) and then to look at the weak scores and the own fraction
 * as secondary options … when you click on other options, they replace the median view"*). They are
 * not three pictures of one thing, which is why all three are stored: flagging by each — worst
 * decile of separation, negative weak margin, own fraction below 1 — only **34 %** (ecoli/Bacformer)
 * and **18 %** (ESM) of flagged loci are flagged by all three, and of the median's worst decile only
 * 50 % / 22 % have a negative margin. Rank ρ 0.79–0.93: they agree overall and part at the tail,
 * which is exactly where flagging happens.
 *
 * ⛔ **The difference is subtracted HERE and is deliberately not served.** `separation` and `margin`
 * are each the difference of two fields the response carries, and a separately-rounded difference
 * drifting from the two rounded numbers printed beside it is the exact class of quiet disagreement
 * this card was rebuilt to end.
 */

import type { LocusSimilarity, PairEstimator, SimilarityBaseline } from "@/api/types";

/** ⚠ The id is what the store keeps, so it must never be a position in the array below. */
export type SimilarityViewId = "median" | "weak" | "own";

/** Every row key a view can show. ⛔ `own_neighbour_fraction` is a SHARE; the rest are cosines. */
export type SimilarityRowKey =
  | "within_similarity"
  | "nearest_similarity"
  | "weak_own_similarity"
  | "weak_other_similarity"
  | "own_neighbour_fraction";

export type SimilarityPercentileKey =
  | "separation_percentile"
  | "weak_margin_percentile"
  | "own_fraction_percentile";

export interface SimilarityRow {
  readonly label: string;
  readonly key: SimilarityRowKey;
}

export interface SimilarityView {
  readonly id: SimilarityViewId;
  /** What the strip's button says. */
  readonly tab: string;
  readonly rows: readonly SimilarityRow[];
  /**
   * What the difference row is CALLED, or `null` where this view has no difference.
   *
   * ⭐ It doubles as the test for "is this view a pair of cosines?" — which is what decides whether
   * the random gene-pair floor may be drawn and whether the difference row is a value or a rank.
   */
  readonly difference: string | null;
  readonly percentileKey: SimilarityPercentileKey;
  /** `[clean, marginal, bad]`, one parallel construction so a reader compares one thing across views. */
  readonly bands: readonly [string, string, string];
  readonly note: string;
  /**
   * What the tab and the note say on a build whose numbers are MEANS (`SimilarityBaseline.estimator`).
   * Absent where the view reads the same either way. ⛔ Read through `viewTab` / `viewNote`, never directly.
   */
  readonly mean?: { readonly tab: string; readonly note: string };
}

/**
 * ⚠ **The pair view is FIRST and is the default** — a median or a mean as the build says (`viewTab`). The other two are secondary readings of the
 * same locus, and opening on either would make the card answer a question nobody asked.
 *
 * ⚠ The band wording qualifies the RELATION between this locus and its nearest rival — never the
 * quality of the locus, which the rows above show is excellent in almost every one of these cases.
 * The two difference views cut at zero and at p5: the zero crossing already sits near p2.5 on the
 * median view, so `Poor` catches the bottom ~2 % alone. The own view has no difference and is cut on
 * the value itself — anything below 1.0 means some member's nearest gene in the whole species is not
 * a member of this locus, which is the thing worth seeing.
 */
export const SIMILARITY_VIEWS: readonly SimilarityView[] = [
  {
    id: "median",
    tab: "median pair",
    rows: [
      { label: "within cluster", key: "within_similarity" },
      { label: "nearest other cluster", key: "nearest_similarity" },
    ],
    difference: "separation",
    percentileKey: "separation_percentile",
    bands: [
      "Clean cluster separation",
      "Unclear cluster separation",
      "Poor cluster separation",
    ],
    note:
      "The median cosine over EVERY pair of genes in this locus, against the highest such median " +
      "against another locus; separation is the difference. Both are over whole sets — no member " +
      "stands in for the locus.",
    // ⛔ Not "centroid": the number is the mean over every PAIR, computed exactly from the centroid (the
    // sum of the members' unit vectors) so that no pair is enumerated. A cosine TO the centroid is a
    // different, larger number — |ū| against (n|ū|² − 1)/(n − 1).
    mean: {
      tab: "mean pair",
      note:
        "The mean cosine over every pair of this locus's embedded genes, against the highest such mean " +
        "against another locus; separation is the difference. Both are exact averages over whole sets, " +
        "computed from each locus's centroid — not a similarity to the centroid, which would read higher. " +
        "Genes without an embedding are not counted; the footer says how many carry one.",
    },
  },
  {
    id: "weak",
    tab: "weakest member",
    rows: [
      { label: "weakest → own", key: "weak_own_similarity" },
      { label: "weakest → other", key: "weak_other_similarity" },
    ],
    difference: "margin",
    percentileKey: "weak_margin_percentile",
    bands: [
      "Weakest member clearly belongs",
      "Weakest member is marginal",
      "Weakest member is closer to another locus",
    ],
    note:
      "The single member least attached to this locus: its nearest neighbour inside the locus, " +
      "against that same gene's nearest gene anywhere else. The point is invariant where a set " +
      "average moves with the set — and the clustering joins points, not averages.",
  },
  {
    id: "own",
    tab: "nearest neighbour",
    rows: [{ label: "nearest gene is own", key: "own_neighbour_fraction" }],
    difference: null,
    percentileKey: "own_fraction_percentile",
    bands: [
      "Every member's nearest gene is in this locus",
      "Some members' nearest gene is elsewhere",
      "Most members' nearest gene is elsewhere",
    ],
    note:
      "The share of this locus's members whose nearest gene in the whole species is another member " +
      "of it. Singletons are excluded: a single gene has no own neighbour to be nearest to.",
  },
];

/**
 * One value marked on the floor strip's axis.
 *
 * ⚠ **`isSecondRow` is the view-row index, not a position among the ticks that survived**: with a
 * missing first row the rival would otherwise be drawn as the subject — filled rather than dotted.
 *
 * ⛔ It lives here rather than in `FloorStrip.vue` because an SFC `<script setup>` block may not
 * export a type at all, and a strip whose tick shape is unnamed is one every caller re-guesses.
 */
export interface FloorTick {
  readonly label: string;
  readonly value: number;
  readonly isSecondRow: boolean;
}

export function viewById(id: SimilarityViewId): SimilarityView {
  // ⛔ Non-null by construction rather than by luck: the union above and the table are written
  // together, and a lookup that could return `undefined` would put that `undefined` on the page.
  return SIMILARITY_VIEWS.find((view) => view.id === id) ?? (SIMILARITY_VIEWS[0] as SimilarityView);
}

/**
 * The estimator the card speaks for: one if every representation agrees, `mixed` if they do not.
 *
 * ⚠ `mixed` is reachable by no served catalogue — each build runs one producer for both — and is handled
 * so that a future one gets a sentence that is true rather than one borrowed from either side.
 */
export type CardEstimator = PairEstimator | "mixed";

export function cardEstimator(baselines: readonly SimilarityBaseline[]): CardEstimator {
  const seen = new Set(baselines.map((baseline) => baseline.estimator));
  if (seen.size > 1) return "mixed";
  return seen.has("mean") ? "mean" : "median";
}

const MIXED_TAB = "all pairs";
const MIXED_NOTE =
  "An average cosine over every pair of genes in this locus, against the highest such average against " +
  "another locus; separation is the difference. On this build one representation's numbers are means " +
  "and the other's medians.";

export function viewTab(view: SimilarityView, estimator: CardEstimator): string {
  if (!view.mean || estimator === "median") return view.tab;
  return estimator === "mean" ? view.mean.tab : MIXED_TAB;
}

export function viewNote(view: SimilarityView, estimator: CardEstimator): string {
  if (!view.mean || estimator === "median") return view.note;
  return estimator === "mean" ? view.mean.note : MIXED_NOTE;
}

export type SimilarityClass = "win" | "warn" | "bad";

export interface SimilarityVerdict {
  /** The quantity the view is ranked on: the difference where it has one, the value where it does not. */
  readonly value: number;
  /** The precomputed midrank, 0..1, over MEASURABLE loci only. */
  readonly percentile: number;
  readonly verdictClass: SimilarityClass;
  readonly label: string;
}

/**
 * ⛔ **Whether this representation has anything to say in this view**, and the one field it must
 * NOT be decided on is `nearest_similarity`: that question is well posed for a single gene and is
 * populated for a singleton, while `within_similarity`, both weak numbers and the own fraction are
 * all null there. Testing the first row is what makes a singleton fall through to its sentence.
 */
export function isMeasurable(
  similarity: LocusSimilarity | null,
  view: SimilarityView,
): similarity is LocusSimilarity {
  return similarity !== null && similarity[(view.rows[0] as SimilarityRow).key] !== null;
}

/**
 * The quantity a view is RANKED on.
 *
 * ⛔ A view with a second row and no value for it has **no score** — never a silent zero, which
 * would put a locus that was not measured at the same place as one whose two sides are equal.
 */
export function viewValue(similarity: LocusSimilarity | null, view: SimilarityView): number | null {
  if (similarity === null) return null;
  const first = similarity[(view.rows[0] as SimilarityRow).key];
  if (first === null) return null;
  if (view.difference === null) return first;
  const second = similarity[(view.rows[1] as SimilarityRow).key];
  return second === null ? null : first - second;
}

/**
 * Three bands, one parallel construction — see the table above for the cuts and why they sit there.
 *
 * ⛔ Returns `null` where the value or the midrank is missing: *not measurable*, which the card must
 * say in words and must never render as `0.000`. A singleton has no pair inside its locus at all;
 * that is an absence, not a bad score.
 */
export function viewVerdict(
  similarity: LocusSimilarity | null,
  view: SimilarityView,
): SimilarityVerdict | null {
  const value = viewValue(similarity, view);
  const percentile = similarity === null ? null : similarity[view.percentileKey];
  if (value === null || percentile === null) return null;
  const level =
    view.difference !== null
      ? value < 0
        ? 2
        : percentile < 0.05
          ? 1
          : 0
      : value < 0.5
        ? 2
        : value < 1
          ? 1
          : 0;
  return {
    value,
    percentile,
    verdictClass: (["win", "warn", "bad"] as const)[level] as SimilarityClass,
    label: view.bands[level] as string,
  };
}

/**
 * A row's value, formatted the way that row may be read.
 *
 * ⛔ **The own fraction is never formatted at 0 dp**, and `lib/formatting`'s `sharePercent` is
 * deliberately not used for it: that keeps a decimal only below 10 %, so 0.9999 would print
 * **"100%"** — the one reading this number must never give, because below 1 means some member's
 * nearest gene in the whole species belongs to another locus. At 1.0 exactly "100%" is the truth;
 * anything that merely rounds there gets `<100%` instead.
 */
export function formatRowValue(value: number, key: SimilarityRowKey): string {
  if (key !== "own_neighbour_fraction") return value.toFixed(3);
  if (value >= 1) return "100%";
  if (Math.round(1000 * value) >= 1000) return "<100%";
  return `${(100 * value).toFixed(1)}%`;
}
