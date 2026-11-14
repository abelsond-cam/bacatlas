import { describe, expect, it } from "vitest";

import type { LocusSimilarity, PairEstimator, SimilarityBaseline } from "@/api/types";

import {
  cardEstimator,
  formatRowValue,
  isMeasurable,
  SIMILARITY_VIEWS,
  viewById,
  viewNote,
  viewTab,
  viewValue,
  viewVerdict,
  type SimilarityView,
} from "./similarityViews";

const MEDIAN = viewById("median");
const WEAK = viewById("weak");
const OWN = viewById("own");

function similarity(overrides: Partial<LocusSimilarity> = {}): LocusSimilarity {
  return {
    within_similarity: 0.94,
    nearest_similarity: 0.61,
    weak_own_similarity: 0.78,
    weak_other_similarity: 0.74,
    own_neighbour_fraction: 1,
    separation_percentile: 0.62,
    weak_margin_percentile: 0.55,
    own_fraction_percentile: 0.53,
    nearest_loci: [],
    ...overrides,
  };
}

describe("⭐ three views, in the order the strip offers them", () => {
  it("opens on the median pair, and it is FIRST", () => {
    // The default is a real choice about which question the card asks: the reading over whole sets.
    // Both of the others are secondary readings of the same locus.
    expect(SIMILARITY_VIEWS.map((view) => view.id)).toEqual(["median", "weak", "own"]);
    expect(viewById("median")).toBe(SIMILARITY_VIEWS[0]);
  });

  it("⚠ the id is what the store keeps, so a lookup never depends on a position", () => {
    for (const view of SIMILARITY_VIEWS) expect(viewById(view.id).id).toBe(view.id);
  });

  it("⛔ only the two COSINE views carry a difference — the own fraction is a share", () => {
    expect([MEDIAN.difference, WEAK.difference, OWN.difference]).toEqual([
      "separation",
      "margin",
      null,
    ]);
  });
});

describe("⛔ the difference is subtracted in the CLIENT, never served", () => {
  it("is the gap between a locus's own members and its nearest rival", () => {
    expect(viewValue(similarity(), MEDIAN)).toBeCloseTo(0.33, 10);
    expect(viewValue(similarity(), WEAK)).toBeCloseTo(0.04, 10);
  });

  it("⛔ has NO score where the second row is missing — never a silent zero", () => {
    // A zero would put a locus that was not measured at the same place as one whose two sides are
    // exactly equal, which is a real and very different reading.
    expect(viewValue(similarity({ nearest_similarity: null }), MEDIAN)).toBeNull();
    expect(viewValue(similarity({ weak_other_similarity: null }), WEAK)).toBeNull();
  });

  it("a view with no difference is ranked on its own value", () => {
    expect(viewValue(similarity({ own_neighbour_fraction: 0.97 }), OWN)).toBeCloseTo(0.97, 10);
  });
});

describe("⛔ measurability is decided on the FIRST row, never on `nearest_similarity`", () => {
  /**
   * ⚠ The trap this pins: a singleton HAS a `nearest_similarity` — that question is well posed for
   * one gene — while `within_similarity`, both weak numbers and the own fraction are all null. A
   * card testing the nearest would give a singleton a block with one rail and no partner.
   */
  const singleton = similarity({
    within_similarity: null,
    weak_own_similarity: null,
    weak_other_similarity: null,
    own_neighbour_fraction: null,
    separation_percentile: null,
    weak_margin_percentile: null,
    own_fraction_percentile: null,
  });

  it("finds a singleton unmeasurable in every view, although it carries a nearest similarity", () => {
    expect(singleton.nearest_similarity).not.toBeNull();
    for (const view of SIMILARITY_VIEWS) expect(isMeasurable(singleton, view)).toBe(false);
  });

  it("finds an ordinary locus measurable in every view, and a missing row not at all", () => {
    for (const view of SIMILARITY_VIEWS) expect(isMeasurable(similarity(), view)).toBe(true);
    for (const view of SIMILARITY_VIEWS) expect(isMeasurable(null, view)).toBe(false);
  });
});

describe("the three bands, one parallel construction per view", () => {
  it("calls a NEGATIVE difference poor — the nearest rival is closer than its own members", () => {
    const bad = similarity({ within_similarity: 0.6, nearest_similarity: 0.7, separation_percentile: 0.01 });
    expect(viewVerdict(bad, MEDIAN)?.verdictClass).toBe("bad");
    expect(viewVerdict(bad, MEDIAN)?.label).toBe("Poor cluster separation");
  });

  it("⚠ calls the bottom 5 % unclear, and everything above it clean", () => {
    // Cut at p5 because the zero crossing already sits near p2.5, so `Poor` catches the bottom ~2 %
    // alone and `Unclear` adds a thin sliver.
    expect(viewVerdict(similarity({ separation_percentile: 0.04 }), MEDIAN)?.label).toBe(
      "Unclear cluster separation",
    );
    expect(viewVerdict(similarity({ separation_percentile: 0.05 }), MEDIAN)?.label).toBe(
      "Clean cluster separation",
    );
    expect(viewVerdict(similarity({ separation_percentile: 0.9 }), MEDIAN)?.verdictClass).toBe("win");
  });

  it("⭐ bands the weakest member on ITS OWN wording, not on separation's", () => {
    // The two views flag different loci — of the median's worst decile only 22–50 % have a negative
    // margin — so one set of words for both would say the same thing about two different findings.
    const closer = similarity({ weak_own_similarity: 0.41, weak_other_similarity: 0.57 });
    expect(viewVerdict(closer, WEAK)?.label).toBe("Weakest member is closer to another locus");
    expect(viewVerdict(similarity({ weak_margin_percentile: 0.01 }), WEAK)?.label).toBe(
      "Weakest member is marginal",
    );
  });

  it("⛔ cuts the own fraction on the VALUE, because anything below 1.0 is the thing worth seeing", () => {
    // Below 1 means some member's nearest gene in the whole species is not a member of this locus.
    expect(viewVerdict(similarity(), OWN)?.verdictClass).toBe("win");
    expect(viewVerdict(similarity({ own_neighbour_fraction: 0.97 }), OWN)?.verdictClass).toBe("warn");
    expect(viewVerdict(similarity({ own_neighbour_fraction: 0.4 }), OWN)?.verdictClass).toBe("bad");
  });

  it("⛔ returns null for NOT MEASURABLE, which the card must say in words", () => {
    // A singleton has no pair inside its locus at all. That is an absence, not a bad score, and
    // rendering it as 0.000 would put a made-up number on the page.
    expect(viewVerdict(null, MEDIAN)).toBeNull();
    expect(viewVerdict(similarity({ separation_percentile: null }), MEDIAN)).toBeNull();
    expect(viewVerdict(similarity({ within_similarity: null }), MEDIAN)).toBeNull();
  });

  it("⚠ a measured percentile of ZERO is a rank, not an absence", () => {
    expect(viewVerdict(similarity({ separation_percentile: 0 }), MEDIAN)?.verdictClass).toBe("warn");
  });
});

describe("⛔ the own fraction is never rounded to a whole percent", () => {
  it("prints 100% only at exactly 1.0", () => {
    expect(formatRowValue(1, "own_neighbour_fraction")).toBe("100%");
  });

  it("⛔⛔ says `<100%` where 0 dp would have said 100% — the one reading it must never give", () => {
    // 0.9999 of 143 members is one gene whose nearest neighbour in the whole species belongs to
    // another locus. `pct(v, 0)` renders that as "100%", which is the opposite claim.
    expect(formatRowValue(0.9999, "own_neighbour_fraction")).toBe("<100%");
    expect(formatRowValue(0.99951, "own_neighbour_fraction")).toBe("<100%");
  });

  it("keeps one decimal everywhere else, where `sharePercent` would drop it above 10 %", () => {
    expect(formatRowValue(0.951049, "own_neighbour_fraction")).toBe("95.1%");
    expect(formatRowValue(0.986842, "own_neighbour_fraction")).toBe("98.7%");
    expect(formatRowValue(0.994, "own_neighbour_fraction")).toBe("99.4%");
  });

  it("⚠ prints every COSINE row at 3 dp, including a measured 1.0", () => {
    // ESM's within-cluster median is 1.0 on real loci, and it is a measurement rather than a
    // rounding artefact: a card that dropped it would drop the rail entirely.
    for (const view of [MEDIAN, WEAK] satisfies SimilarityView[]) {
      for (const row of view.rows) expect(formatRowValue(1, row.key)).toBe("1.000");
    }
    expect(formatRowValue(0.814065, "within_similarity")).toBe("0.814");
  });
});

describe("⛔ the pair view names the statistic the build actually holds", () => {
  function baseline(representation: "esm" | "bacformer", estimator: PairEstimator): SimilarityBaseline {
    return {
      representation, form: "raw", floor_median: 0.5, floor_p25: 0.4, floor_p75: 0.6, floor_p99: 0.9,
      measurable_locus_count: 10, neighbour_knn_k: 200, estimator, unmeasured_views: [],
    };
  }

  it("is a median unless every representation says mean, and `mixed` when they disagree", () => {
    expect(cardEstimator([baseline("esm", "median"), baseline("bacformer", "median")])).toBe("median");
    expect(cardEstimator([baseline("esm", "mean"), baseline("bacformer", "mean")])).toBe("mean");
    expect(cardEstimator([baseline("esm", "mean"), baseline("bacformer", "median")])).toBe("mixed");
    // No baseline at all is the full-audit sentence — there is nothing that could say mean.
    expect(cardEstimator([])).toBe("median");
  });

  it("changes the pair view's tab and note, and leaves the views that read the same either way alone", () => {
    expect(viewTab(MEDIAN, "median")).toBe("median pair");
    expect(viewTab(MEDIAN, "mean")).toBe("mean pair");
    expect(viewNote(MEDIAN, "median")).toBe(MEDIAN.note);
    expect(viewNote(MEDIAN, "mean")).toMatch(/^The mean cosine over every pair of this locus's embedded genes/);
    expect(viewNote(MEDIAN, "mean")).not.toMatch(/median/i);
    for (const view of [WEAK, OWN]) {
      for (const estimator of ["median", "mean", "mixed"] as const) {
        expect(viewTab(view, estimator)).toBe(view.tab);
        expect(viewNote(view, estimator)).toBe(view.note);
      }
    }
  });

  it("⚠ no other view's note claims a median — the weakest member's reads true on either build", () => {
    expect(WEAK.note).not.toMatch(/median/i);
    expect(OWN.note).not.toMatch(/median/i);
  });
});
