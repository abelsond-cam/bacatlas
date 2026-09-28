/**
 * @vitest-environment jsdom
 *
 * ⛔ **jsdom computes no layout and the harness loads no stylesheet.** Every assertion about size
 * below can only see an INLINE style, which is why the track sets its widths inline and why moving
 * one into CSS would leave these tests green while measuring nothing. That is the rule
 * `app.js:1198` states, carried forward.
 */
import { mount } from "@vue/test-utils";
import { createPinia, setActivePinia } from "pinia";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { asCatalogueKey } from "@/api/types";
import type {
  Arrangement,
  IntergenicGap,
  Locus,
  LocusDetailResponse,
  NeighbourDisplayRow,
  OffsetMarginal,
} from "@/api/types";
import { SIGNED_OFFSETS, type NeighbourSlot } from "@/lib/slotSpaces";
import { NARROW_PX, NO_LENGTH_PX } from "@/lib/trackGeometry";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { useTrackDisplayStore } from "@/stores/trackDisplayStore";

import GeneTrack from "./GeneTrack.vue";
import IntergenicSlot from "./IntergenicSlot.vue";
import NeighbourGeneSlot from "./NeighbourGeneSlot.vue";

// ⚠ A PARTIAL mock, like every other suite that mocks this module. `anchorQuery` must stay real:
// it is the one place the anchor's kind becomes a query parameter, and a stub would hide a projected
// genome being fetched as `anchor=` — which 404s and reads on the page as "this genome has nothing here".
// ⛔ This was a TOTAL mock until 2026-09-28, which made every other export throw on access. The hover
// test below then exploded inside `cache.prefetch` on its first line, every run — and still PASSED,
// because it only asserts that `prefetchNeighbour` was CALLED. The whole body of that prefetch was dead
// under this suite. The unhandled rejection was the only thing saying so.
const fetchLocus = vi.hoisted(() => vi.fn(async () => ({ ok: false, kind: "network", detail: "x" })));
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  fetchLocus,
}));

const NEIGHBOUR: NeighbourDisplayRow = {
  label: "77",
  catalogue_ordinal: 77,
  best_product: null,
  display_name: "rfaL",
  display_name_source: "bakta_symbol",
  genome_count: 90,
  median_gene_length_nt: 1230,
  prevalence_band: "core",
};

function slot(overrides: Partial<NeighbourSlot> = {}): NeighbourSlot {
  return { signed_offset: -1, locus: "77", absence_reason: null, same_strand: true, ...overrides };
}

function marginal(overrides: Partial<OffsetMarginal> = {}): OffsetMarginal {
  return {
    signed_offset: -1,
    observed_member_count: 100,
    observed_not_listed: 10,
    members_without_an_observation: 3,
    occupants: [
      { rank: 0, locus: "77", gene_count: 70, same_strand_gene_count: 70 },
      { rank: 1, locus: "78", gene_count: 20, same_strand_gene_count: 2 },
    ],
    ...overrides,
  };
}

function mountSlot(props: Partial<InstanceType<typeof NeighbourGeneSlot>["$props"]> = {}) {
  return mount(NeighbourGeneSlot, {
    props: {
      slot: slot(),
      neighbour: NEIGHBOUR,
      labelOffset: -1,
      marginal: marginal(),
      neighboursByLabel: new Map([[NEIGHBOUR.label, NEIGHBOUR]]),
      collectionGenomeCount: 100,
      isPopoverOpen: false,
      ...props,
    },
  });
}

beforeEach(() => setActivePinia(createPinia()));

describe("⛔ the inline width, which is the only version any test can see", () => {
  it("is set on the SLOT, from the neighbour's median gene length", () => {
    // The bar and the offset label are the block's siblings and must line up under it exactly, so
    // the slot carries the width and not the block.
    const wrapper = mountSlot();
    expect(wrapper.attributes("style")).toContain("width: 123px");
    expect(wrapper.find(".block").attributes("style") ?? "").not.toContain("width");
  });

  it("⛔ a slot with NO occupant takes the fallback width, never a scaled one", () => {
    const wrapper = mountSlot({
      slot: slot({ locus: null, absence_reason: "contig_end", same_strand: null }),
      neighbour: null,
    });
    expect(wrapper.attributes("style")).toContain(`width: ${NO_LENGTH_PX}px`);
  });

  it("marks a narrow slot, and hides what no longer fits inside it", () => {
    const narrow = { ...NEIGHBOUR, median_gene_length_nt: 300 };
    const wrapper = mountSlot({ neighbour: narrow });
    expect(wrapper.classes()).toContain("narrow");
    // Below 64 px a block cannot hold its locus id and product; they move to the popover.
    expect(wrapper.find(".lid").exists()).toBe(false);
    expect(wrapper.find(".gsub").exists()).toBe(false);
  });

  it("⛔ carries the neighbour's PRODUCT under its id, not its name a second time", () => {
    // The name is already above the block; `app.js:1609` puts `bestProduct` in the block, and that
    // is what tells a reader what the neighbour does.
    const wrapper = mountSlot({ neighbour: { ...NEIGHBOUR, display_name: "yfcU", best_product: "fimbrial usher protein YfcU" } });
    expect(wrapper.find(".gsub").text()).toBe("fimbrial usher protein YfcU");
    expect(wrapper.find(".gname").text()).toBe("yfcU");
  });

  it("⚠ and does NOT mark one of exactly 64 px, because the threshold is a strict less-than", () => {
    const wrapper = mountSlot({ neighbour: { ...NEIGHBOUR, median_gene_length_nt: NARROW_PX * 10 } });
    expect(wrapper.classes()).not.toContain("narrow");
    expect(wrapper.find(".lid").exists()).toBe(true);
  });
});

describe("⛔ data-dir is the strand relation AS SHOWN", () => {
  it("points forward for a co-oriented neighbour and back for an opposed one", () => {
    expect(mountSlot().attributes("data-dir")).toBe("fwd");
    expect(mountSlot({ slot: slot({ same_strand: false }) }).attributes("data-dir")).toBe("rev");
  });

  it("⛔ a contig end keeps the FORWARD lane and carries no direction on its block", () => {
    // An absence must not be put in a direction lane: that would assert a transcription direction
    // nothing observed.
    const wrapper = mountSlot({
      slot: slot({ locus: null, absence_reason: "contig_end", same_strand: null }),
      neighbour: null,
    });
    expect(wrapper.attributes("data-dir")).toBe("fwd");
    expect(wrapper.find(".block").attributes("data-dir")).toBeUndefined();
    expect(wrapper.find(".block").attributes("disabled")).toBeDefined();
  });

  it("emits the SHOWN strand on a walk, which is what decides the next frame", async () => {
    const wrapper = mountSlot({ slot: slot({ same_strand: false }) });
    await wrapper.find(".block").trigger("click");
    expect(wrapper.emitted("walk")).toEqual([["77", false]]);
  });

  it("a slot with no occupant cannot be walked from at all", async () => {
    const wrapper = mountSlot({
      slot: slot({ locus: null, absence_reason: "contig_end", same_strand: null }),
      neighbour: null,
    });
    await wrapper.find(".block").trigger("click");
    expect(wrapper.emitted("walk")).toBeUndefined();
  });
});

describe("the block's two colour axes stay independent", () => {
  it("⚠ an occupant outside the exported top-N is UNKNOWN, not a zero share", () => {
    // The ramp would paint a zero as the most divergent purple it has, so it comes off the ramp.
    const wrapper = mountSlot({ slot: slot({ locus: "999" }), neighbour: { ...NEIGHBOUR, label: "999" } });
    expect(wrapper.find(".block").classes()).toContain("unknown");
  });

  it("the hue is the NAMED gene's own share, and the alpha is the locus's coverage", () => {
    const wrapper = mountSlot();
    const style = wrapper.find(".fill").attributes("style") ?? "";
    expect(style).toContain("--c: 0.700"); // 70 of 100 observed members
    expect(style).toContain("--p: 0.900"); // the locus is in 90 of 100 genomes
  });
});

describe("the marginal bar", () => {
  it("draws the display cut as its OWN segment, not merged with the occupants", () => {
    const wrapper = mountSlot();
    const segments = wrapper.findAll(".marg i");
    expect(segments).toHaveLength(3);
    expect(segments[2]!.classes()).toContain("oth");
    expect(segments[2]!.attributes("style")).toContain("width: 10%");
  });

  it("⚠ a position with nothing observed gets a plain bar and NO button", () => {
    // There is nothing to open, and a control that could only ever answer "nothing here" is worse
    // than none.
    const wrapper = mountSlot({ marginal: marginal({ observed_member_count: 0, occupants: [] }) });
    expect(wrapper.find(".marg-hit").exists()).toBe(false);
    expect(wrapper.find(".marg").exists()).toBe(true);
  });

  it("opens the popover for the position, and says so to a screen reader", async () => {
    const wrapper = mountSlot();
    await wrapper.find(".marg-hit").trigger("click");
    expect(wrapper.emitted("togglePopover")).toHaveLength(1);
    expect(wrapper.find(".marg-hit").attributes("aria-label")).toContain("A-1");
  });
});

describe("the intergenic block is TWO objects, not one", () => {
  function gap(overrides: Partial<IntergenicGap> = {}): IntergenicGap {
    return {
      flanking_loci: ["1", "2"],
      observed_genome_count: 90,
      median_signed_length_nt: 120,
      quartile1_signed_length_nt: 118,
      quartile3_signed_length_nt: 121,
      minimum_signed_length_nt: 118,
      maximum_signed_length_nt: 121,
      length_variance_score: 0.25,
      modal_length_nt: 120,
      distinct_named_feature_count: 0,
      every_genome_agrees: false,
      ...overrides,
    };
  }

  it("a positive median is a REGION, scaled", () => {
    const wrapper = mount(IntergenicSlot, { props: { gap: gap() } });
    expect(wrapper.classes()).toContain("gap");
    expect(wrapper.attributes("style")).toContain("width: 12px");
  });

  it("⛔ a negative median is an OVERLAP, drawn as a seam and never as a small region", () => {
    // 18.8 % of adjacent pairs overlap. Clamping at zero made "abuts exactly" and "overlaps by 190
    // bases" the same number, and a seam is not a region at all.
    const wrapper = mount(IntergenicSlot, { props: { gap: gap({ median_signed_length_nt: -4 }) } });
    expect(wrapper.classes()).toContain("joint");
    expect(wrapper.find(".block").classes()).toContain("seam");
    expect(wrapper.find(".block").attributes("data-tip")).toContain("4 bases of overlap");
  });

  it("⛔ a MEASURED ZERO variance still renders as measured — the majority case", () => {
    // `v-if="varianceScore"` is false for a measured zero, and white on the track has to mean
    // "identical in every genome", never "small". 80–90 % of gaps vary not at all.
    const wrapper = mount(IntergenicSlot, {
      props: { gap: gap({ length_variance_score: 0, every_genome_agrees: true }) },
    });
    expect(wrapper.find(".block").classes()).not.toContain("unmeasured");
    expect(wrapper.find(".block").attributes("style")).toContain("--v: 0.000");
    expect(wrapper.find(".block").attributes("data-tip")).toContain("identical in every genome");
  });

  it("⛔ tints through the page's display CURVE, not with the raw score", () => {
    // `app.js:1406` sets `--v` to `gapVar(k)` = sqrt(score / 0.5). The raw score here is 0.125, which
    // unmapped is a 12.5 % wash; through the curve it is 50 %. Passing the raw score leaves nine
    // varying regions in ten looking identical to the ones that never vary.
    const wrapper = mount(IntergenicSlot, {
      props: { gap: gap({ length_variance_score: 0.125, every_genome_agrees: false }) },
    });
    expect(wrapper.find(".block").attributes("style")).toContain("--v: 0.500");
  });

  it("⛔ a NULL variance is not measured, and says so differently", () => {
    const wrapper = mount(IntergenicSlot, { props: { gap: gap({ length_variance_score: null }) } });
    expect(wrapper.find(".block").classes()).toContain("unmeasured");
  });

  it("⛔ a null MEDIAN draws nothing at all, never a zero-length region", () => {
    // Drawing one would assert the two genes abut exactly — a measurement nobody made.
    const wrapper = mount(IntergenicSlot, { props: { gap: gap({ median_signed_length_nt: null }) } });
    expect(wrapper.find(".slot").exists()).toBe(false);
  });

  it("and an absent gap draws nothing", () => {
    expect(mount(IntergenicSlot, { props: { gap: undefined } }).find(".slot").exists()).toBe(false);
  });
});

describe("the whole track", () => {
  function arrangement(overrides: Partial<Arrangement> = {}): Arrangement {
    return {
      rank: 0,
      gene_count: 80,
      genome_count: 80,
      is_recorded_reverse_complement: false,
      slots: SIGNED_OFFSETS.map((signed_offset, position) => ({
        signed_offset,
        locus: `n${position}`,
        absence_reason: null,
        same_strand: position % 2 === 0,
      })),
      ...overrides,
    };
  }

  function detail(): LocusDetailResponse {
    return {
      locus: {
        label: "focal",
        display_name: "wzi",
        best_product: "outer membrane protein",
        median_gene_length_nt: 1500,
        gene_count: 100,
        genome_count: 97,
      } as Locus,
      annotations: {},
      uniref50_families: [],
      arrangements: {
        listed: [arrangement()],
        total: 12,
        arrangements_not_listed: 11,
        members_in_arrangements_not_listed: 15,
        members_without_a_neighbourhood: 5,
        membership_is_complete: true,
      },
      anchor: { is_anchored: false, arrangement_ranks: [], kind: null, projected_copies: [] },
      offsets: SIGNED_OFFSETS.map((signed_offset) => marginal({ signed_offset })),
      intergenic_gaps: [],
      pfam_reference: {},
      neighbour_display_rows: SIGNED_OFFSETS.map((_unused, position) => ({
        ...NEIGHBOUR,
        label: `n${position}`,
      })),
      resolved_neighbour_count: 10,
    };
  }

  function mountTrack() {
    const navigation = useLocusNavigationStore();
    navigation.setCatalogue(asCatalogueKey("ecoli"));
    navigation.route = { label: "focal", direction: "forward" };
    navigation.view = { status: "ready", value: detail() };
    return mount(GeneTrack, {
      props: { detail: detail(), collectionGenomeCount: 100 },
    });
  }

  it("draws ten gene slots plus one focal", () => {
    const wrapper = mountTrack();
    // ⛔ Counted by EXCLUDING the non-gene kinds by name. The ±5 window is five GENES either side,
    // and a count that swept up intergenic blocks would stop being able to notice that moving.
    const geneSlots = wrapper
      .findAll(".slot")
      .filter((node) => !node.classes("gap") && !node.classes("joint"));
    expect(geneSlots).toHaveLength(11);
    expect(wrapper.findAll(".slot.focal")).toHaveLength(1);
  });

  it("puts the focal gene in the middle, after five neighbours", () => {
    const wrapper = mountTrack();
    const geneSlots = wrapper
      .findAll(".slot")
      .filter((node) => !node.classes("gap") && !node.classes("joint"));
    expect(geneSlots[5]!.classes()).toContain("focal");
  });

  it("⛔ draws the two focal remainders as TWO segments", () => {
    // "In an arrangement past the display cut" and "no coordinates for the gene, so no window" are
    // different sentences, and merging them tells a reader the first has no coordinates.
    const wrapper = mountTrack();
    const bar = wrapper.find(".slot.focal .marg");
    expect(bar.find("i.oth").exists()).toBe(true);
    expect(bar.find("i.nowin").exists()).toBe(true);
    expect(bar.find("i.nowin").attributes("data-tip")).toContain("no recorded neighbourhood");
    expect(bar.find("i.oth").attributes("data-tip")).toContain("past the display cut");
  });

  it("labels the columns A-5 … A+5 in order", () => {
    const wrapper = mountTrack();
    const labels = wrapper
      .findAll(".slot:not(.focal):not(.gap):not(.joint) .marg-hit")
      .map((node) => (node.attributes("aria-label") ?? "").match(/A[+-]\d/)?.[0]);
    expect(labels).toEqual(["A-5", "A-4", "A-3", "A-2", "A-1", "A+1", "A+2", "A+3", "A+4", "A+5"]);
  });

  it("⭐ prefetches on hover, so the walk feels synchronous", async () => {
    const wrapper = mountTrack();
    const track = useTrackDisplayStore();
    const prefetch = vi.spyOn(track, "prefetchNeighbour");
    fetchLocus.mockClear();
    await wrapper.find(".slot:not(.focal):not(.gap):not(.joint)").trigger("pointerenter");
    expect(prefetch).toHaveBeenCalledWith("ecoli", "n0");
    // ⛔ The assertion above says the hover handler FIRED. It passed for five days while the prefetch
    // it fires threw on its own first line — the store reads `anchorQuery`, which a total mock had
    // replaced with a throwing proxy. Nothing caught it, because `prefetchNeighbour` discards the
    // promise with `void` by design. So assert the prefetch REACHED the network too: this is the line
    // that would have gone red, and the whole body of `locusDetailCacheStore.prefetch` depends on it.
    expect(fetchLocus).toHaveBeenCalledWith("ecoli", "n0", {});
  });
});
