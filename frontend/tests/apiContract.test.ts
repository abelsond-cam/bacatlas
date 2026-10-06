/**
 * @vitest-environment jsdom
 *
 * **The API/front-end boundary, checked against real bytes.**
 *
 * ⭐ Everything in `src/api/types.ts` is hand-transcribed from `serialisers/locus_serialiser.py`,
 * and a transcription error there is invisible: TypeScript happily believes a field exists that the
 * server never sends, and the component renders `undefined` as an empty string. So the fixtures here
 * are **real responses from the Flask test client, for BOTH species** (`api_responses_{species}.json`),
 * and the real components are mounted against them. Every suite runs once per catalogue.
 *
 * The loci are chosen for what they exercise, by a query that says so (see the dump script):
 * - `ordinary` — every arrangement fits, both members-remainders are zero, and ESM's within-cluster
 *   similarity is a MEASURED 1.0 (it was a measured ZERO distance until the medoid geometry went).
 * - `over_cap` — more arrangements than the API lists: on ecoli 84 arrangements, 8 listed, 86 members
 *   past the cap and 7 with no window. Before the remainder split this locus told a reader that 93
 *   member genes had no coordinates.
 * - `no_window` — the most members with no recorded neighbourhood at all.
 * - `function_rich` — EC and KEGG and all three GO namespaces.
 *
 * ⚠ Regenerate with `backend/scripts/dump_api_locus_fixture.py` whenever the serialiser changes. A
 * fixture that silently goes stale is worse than none, so the shape checks below are exhaustive
 * over the keys the client reads rather than a spot check.
 *
 * ⛔ **This suite has already earned itself.** Mounting the real component against real bytes found
 * that `_neighbour_display_rows` looked catalogue ordinals up in the *locus id* set — two key
 * spaces over the same small integers, merged one screen below the comment warning about exactly
 * that. It mostly worked, because an arrangement occupant is usually also a marginal mode; the
 * ones that are not rendered as blank, unwalkable blocks. On the three ecoli loci it was first
 * recorded on, 37 of 190 slots.
 */
import { asCatalogueKey } from "@/api/types";
import { mount } from "@vue/test-utils";
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";

import { createPinia, setActivePinia } from "pinia";
import { beforeEach, describe, expect, it, vi } from "vitest";

import type {
  AnnotationEntry,
  AuditResidualsResponse,
  GenomeListResponse,
  SearchResponse,
  SpeciesCatalogueResponse,
  FunctionResponse,
  GeneOntologyNamespace,
  GeneSequenceResponse,
  LocusDetailResponse,
} from "@/api/types";
import { SLOT_COUNT } from "@/lib/slotSpaces";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";

import ArrangementPopover from "@/components/popover/ArrangementPopover.vue";
import FunctionTab from "@/components/function/FunctionTab.vue";
import SequenceTab from "@/components/sequence/SequenceTab.vue";
import EmbeddingSimilarityCard from "@/components/locusCard/EmbeddingSimilarityCard.vue";
import LocusHeadline from "@/components/locusCard/LocusHeadline.vue";
import SequenceDiversityCard from "@/components/locusCard/SequenceDiversityCard.vue";
import { PREVALENCE_BANDS } from "@/lib/prevalence";
import ArrangementSwitcher from "@/components/track/ArrangementSwitcher.vue";
import GeneTrack from "@/components/track/GeneTrack.vue";
import PangenomeCensus from "@/components/census/PangenomeCensus.vue";
import SiteFooter from "@/components/layout/SiteFooter.vue";
import NavigatingView from "@/components/views/NavigatingView.vue";

// ⚠ Only the locus FETCH is stubbed, and every other export is the real one: a stubbed helper would
// assert that the fixture agrees with the stub.
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  fetchLocus: vi.fn(async () => ({ ok: false, kind: "network", detail: "not used" })),
}));

// ⚠ Resolved from the project root, NOT from `import.meta.url`. This suite runs under jsdom,
// where `import.meta.url` is an http URL rather than a file one, so `fileURLToPath` yields an
// absolute path rooted at `/` and the read fails with a bare ENOENT that reads as a missing
// fixture. Vitest runs with the project root as cwd.
const FIXTURES = resolve(process.cwd(), "tests/fixtures");
/** ⭐ BOTH species. Every suite below runs once per catalogue, on bytes recorded from each. */
const SPECIES_KEYS = ["ecoli", "kp"] as const;
type SpeciesKey = (typeof SPECIES_KEYS)[number];
const fixturePath = (speciesKey: SpeciesKey) => resolve(FIXTURES, `api_responses_${speciesKey}.json`);

interface Recorded {
  readonly recorded_from: string;
  readonly species_key: string;
  readonly loci: Readonly<
    Record<
      string,
      {
        readonly label: string;
        readonly response: LocusDetailResponse;
        /** ⭐ The Function tab's own endpoint, recorded for the SAME locus in the same pass. */
        readonly function: FunctionResponse;
      }
    >
  >;
  /** ⭐ The Sequence tab's own endpoint: a minus-strand gene, and a genome with no gene here. */
  readonly sequences: Readonly<
    Record<
      string,
      {
        readonly sample_id: string;
        readonly locus_label: string;
        readonly response: GeneSequenceResponse;
      }
    >
  >;
  /**
   * ⭐ The species response, recorded in the same pass. The catalogue map's **viewport** comes from
   * this endpoint and its **positions** from the other, so the pair is the only thing that can show
   * the two disagreeing — and a disagreement there is a picture that still looks like a picture.
   */
  readonly species: SpeciesCatalogueResponse;
  /** ⭐ The footer's two residual lists, from their own endpoint. */
  readonly residuals: AuditResidualsResponse;
  /** One real search, so the rows the dropdown draws are real rows. */
  readonly search: { readonly ligase: SearchResponse };
  /** The anchor picker's list, cut at 5 so the `truncated` path is real bytes too. */
  readonly genomes: GenomeListResponse;
}

// ⭐ `function_rich` carries EC *and* KEGG *and* all three GO namespaces — the parts a typical
// locus does not exercise at all. It is in the generic loops too, so every shape assertion gains it.
const CASES = ["ordinary", "over_cap", "no_window", "function_rich", "inferred_function"] as const;

beforeEach(() => setActivePinia(createPinia()));

describe.each(SPECIES_KEYS)("%s", (speciesKey) => {
  const FIXTURE = fixturePath(speciesKey);
  const recorded: Recorded = JSON.parse(readFileSync(FIXTURE, "utf8")) as Recorded;

  function detailFor(name: (typeof CASES)[number]): LocusDetailResponse {
    const entry = recorded.loci[name];
    if (entry === undefined) throw new Error(`fixture has no '${name}' locus`);
    return entry.response;
  }

  describe("the fixture itself", () => {
    it("⛔ exists and says where it came from, rather than being skipped", () => {
      expect(existsSync(FIXTURE)).toBe(true);
      expect(recorded.recorded_from).toBe("the API test client");
      // ⛔ And WHICH catalogue: a copy of one species' file under the other's name would otherwise
      // run every suite twice on the same bytes and report two species checked.
      expect(recorded.species_key).toBe(speciesKey);
      expect(Object.keys(recorded.loci).sort()).toEqual([...CASES].sort());
    });

    it("⭐ actually exercises the case it exists for", () => {
      // Without this the suite would keep passing after a regenerate that happened to pick three
      // ordinary loci, and the discriminating cases would be gone with nothing to say so.
      const overCap = detailFor("over_cap").arrangements;
      expect(overCap.total).toBeGreaterThan(overCap.listed.length);
      expect(overCap.members_in_arrangements_not_listed).toBeGreaterThan(0);
      const noWindow = detailFor("no_window").arrangements;
      expect(noWindow.members_without_a_neighbourhood).toBeGreaterThan(0);
    });
  });

  describe.each(CASES)("%s — the real response satisfies the contract the client reads", (name) => {
    const detail = detailFor(name);

    it("carries every top-level block the client destructures", () => {
      expect(Object.keys(detail).sort()).toEqual(
        [
          "anchor",
          "annotations",
          "arrangements",
          "intergenic_gaps",
          "locus",
          "neighbour_display_rows",
          "offsets",
          "pfam_reference",
          "resolved_neighbour_count",
          "uniref50_families",
        ].sort(),
      );
    });

    it("⛔ has exactly ten offsets, in recorded order", () => {
      expect(detail.offsets).toHaveLength(SLOT_COUNT);
      expect(detail.offsets.map((offset) => offset.signed_offset)).toEqual([
        -5, -4, -3, -2, -1, 1, 2, 3, 4, 5,
      ]);
    });

    it("⛔ every arrangement has exactly ten slots, each null-plus-a-reason or resolved", () => {
      for (const arrangement of detail.arrangements.listed) {
        expect(arrangement.slots).toHaveLength(SLOT_COUNT);
        for (const slot of arrangement.slots) {
          // Never a bare -1: the packed form is where "−1 means five different things" lives.
          expect(typeof slot.signed_offset).toBe("number");
          if (slot.locus === null) {
            expect(slot.absence_reason).not.toBeNull();
            expect(slot.same_strand).toBeNull();
          } else {
            expect(slot.absence_reason).toBeNull();
            expect(typeof slot.same_strand).toBe("boolean");
          }
        }
      }
    });

    it("⛔ names all five remainders, even where they are zero", () => {
      const arrangements = detail.arrangements;
      for (const field of [
        "arrangements_not_listed",
        "members_in_arrangements_not_listed",
        "members_without_a_neighbourhood",
      ] as const) {
        expect(typeof arrangements[field]).toBe("number");
      }
      for (const offset of detail.offsets) {
        expect(typeof offset.observed_not_listed).toBe("number");
        expect(typeof offset.members_without_an_observation).toBe("number");
      }
    });

    it("⚠ the remainders ADD UP, so none of them is quietly absorbing another", () => {
      const arrangements = detail.arrangements;
      const listed = arrangements.listed.reduce((total, one) => total + one.gene_count, 0);
      expect(
        listed +
          arrangements.members_in_arrangements_not_listed +
          arrangements.members_without_a_neighbourhood,
      ).toBe(detail.locus.gene_count);
      expect(arrangements.listed.length + arrangements.arrangements_not_listed).toBe(
        arrangements.total,
      );
    });

    it("resolves every arrangement slot's locus into neighbour_display_rows", () => {
      // ⭐ The fan-out, answered in this one response. A slot naming a locus with no display row
      // would render as a blank block that is nonetheless clickable.
      const known = new Set(detail.neighbour_display_rows.map((row) => row.label));
      const missing = new Set<string>();
      for (const arrangement of detail.arrangements.listed) {
        for (const slot of arrangement.slots) {
          if (slot.locus !== null && !known.has(slot.locus)) missing.add(slot.locus);
        }
      }
      expect([...missing]).toEqual([]);
    });

    it("⭐ mounts the real track component without warnings", () => {
      const warn = vi.spyOn(console, "warn").mockImplementation(() => undefined);
      const error = vi.spyOn(console, "error").mockImplementation(() => undefined);

      const navigation = useLocusNavigationStore();
      navigation.setCatalogue(asCatalogueKey(speciesKey));
      navigation.route = { label: detail.locus.label, direction: "forward" };
      navigation.view = { status: "ready", value: detail };

      const wrapper = mount(GeneTrack, {
        props: { detail, collectionGenomeCount: 123 },
      });

      const geneSlots = wrapper
        .findAll(".slot")
        .filter((node) => !node.classes("gap") && !node.classes("joint"));
      expect(geneSlots).toHaveLength(SLOT_COUNT + 1);
      expect(wrapper.findAll(".slot.focal")).toHaveLength(1);
      // Every gene slot carries an inline width — the one thing a stylesheet-free harness can see.
      for (const node of geneSlots) {
        expect(node.attributes("style") ?? "").toMatch(/width:\s*\d+px/);
      }
      // ⚠ Vue reports a missing required prop or a bad type as a console warning, not a throw, so a
      // shape mismatch between the hand-written types and the real response would otherwise render
      // `undefined` silently and pass every assertion above.
      expect(error).not.toHaveBeenCalled();
      expect(warn).not.toHaveBeenCalled();
      warn.mockRestore();
      error.mockRestore();
    });
  });

  describe("⭐ the remainder split, on the locus that proves it matters", () => {
    it("does NOT tell a reader that members past the cap have no coordinates", () => {
      const arrangements = detailFor("over_cap").arrangements;
      const listed = arrangements.listed.reduce((total, one) => total + one.gene_count, 0);
      // What the old subtraction would have reported as "no coordinates for the gene, so no window".
      const wouldHaveClaimed = detailFor("over_cap").locus.gene_count - listed;
      expect(arrangements.members_without_a_neighbourhood).toBeLessThan(wouldHaveClaimed);
      expect(
        wouldHaveClaimed - arrangements.members_without_a_neighbourhood,
      ).toBe(arrangements.members_in_arrangements_not_listed);
    });
  });

  describe("⭐ the A0 card, mounted on the same real bytes", () => {
    function mountCard(kind: (typeof CASES)[number]) {
      const detail = detailFor(kind);
      return mount(ArrangementPopover, {
        props: {
          locus: detail.locus,
          arrangements: detail.arrangements.listed,
          total: detail.arrangements.total,
          selectedRank: detail.arrangements.listed[0]?.rank ?? null,
          anchorRanks: detail.anchor.arrangement_ranks,
          membersWithoutANeighbourhood: detail.arrangements.members_without_a_neighbourhood,
          walkDirection: "forward" as const,
          arrangementsNotShown: detail.arrangements.arrangements_not_listed,
          loadStatus: "idle" as const,
        },
      });
    }

    it("⛔ derives the display cut to exactly what the SERVER named it", () => {
      // The card cannot use `members_in_arrangements_not_listed` directly — that number is computed
      // against the capped list and stops being true the moment a page arrives — so it derives the
      // same quantity from `gene_count − members_without_a_neighbourhood − Σ listed`. This asserts
      // the two agree on real bytes, which is the only thing that proves the identity holds.
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const card = mountCard(kind);
        const tail = card.find(".pop-tail");
        const expected = detail.arrangements.members_in_arrangements_not_listed;
        if (expected > 0) {
          expect(tail.text()).toContain(`${expected} member genes sit in arrangements not listed here`);
        } else {
          // ⚠ Read the whole card, not the tail: with both remainders zero there IS no tail element,
          // and `.text()` on an absent wrapper throws rather than returning "" — which would have
          // made this branch pass without looking at anything.
          expect(card.text()).not.toContain("sit in arrangements not listed here");
        }
      }
    });

    it("⛔ still names the members with no window as their own, separate sentence", () => {
      const detail = detailFor("no_window");
      expect(detail.arrangements.members_without_a_neighbourhood).toBeGreaterThan(0);
      expect(mountCard("no_window").find(".pop-tail").text()).toContain(
        `${detail.arrangements.members_without_a_neighbourhood} member genes have no recorded neighbourhood`,
      );
    });

    it("⚠ does not claim completeness at a locus where 76 arrangements are missing", () => {
      const detail = detailFor("over_cap");
      expect(detail.arrangements.total).toBeGreaterThan(detail.arrangements.listed.length);
      const lede = mountCard("over_cap").find(".pop-lede").text();
      expect(lede).not.toContain("every one of them");
      expect(lede).toContain(`${detail.arrangements.listed.length} of them listed below`);
      expect(lede).toContain(`sits in ${detail.arrangements.total}`);
    });

    it("says every one IS listed where every one is, and prints no tail at all", () => {
      const detail = detailFor("ordinary");
      expect(detail.arrangements.arrangements_not_listed).toBe(0);
      expect(detail.arrangements.members_in_arrangements_not_listed).toBe(0);
      expect(detail.arrangements.members_without_a_neighbourhood).toBe(0);
      const card = mountCard("ordinary");
      expect(card.find(".pop-lede").text()).toContain("every one of them listed below");
      expect(card.find(".pop-tail").exists()).toBe(false);
      expect(card.find(".arr-more").exists()).toBe(false);
    });

    it("⛔ lists one row per arrangement the response carried, and none of them selectable", () => {
      for (const kind of CASES) {
        const card = mountCard(kind);
        expect(card.findAll(".alt")).toHaveLength(detailFor(kind).arrangements.listed.length);
        expect(card.findAll("button.alt")).toHaveLength(0);
      }
    });
  });

  describe("⭐ the switcher, on the same real bytes", () => {
    function mountSwitcher(kind: (typeof CASES)[number], anchor: string | null = null) {
      const detail = detailFor(kind);
      return mount(ArrangementSwitcher, {
        props: {
          locus: detail.locus,
          arrangements: detail.arrangements.listed,
          total: detail.arrangements.total,
          selectedIndex: 0,
          anchorRanks: detail.anchor.arrangement_ranks,
          anchorGenomeName: anchor,
          membershipIsComplete: detail.arrangements.membership_is_complete,
          // ⚠ These fixtures are MODELLED genomes read from the real payload — never placed
          // ones — so the switcher must take the member's sentence, not the placed one.
          projectedCopies: null,
          walkDirection: "forward" as const,
        },
      });
    }

    it("⛔ the server sends `membership_is_complete`, and the three loci disagree about it", () => {
      // A field every response carried the same value for would be untested by this fixture, so
      // assert the cases actually differ — otherwise the branch below is never taken either way.
      const flags = CASES.map((kind) => detailFor(kind).arrangements.membership_is_complete);
      expect(flags.every((flag) => typeof flag === "boolean")).toBe(true);
      expect(new Set(flags).size).toBe(2);
    });

    it("⛔ picks the anchor sentence from THAT field, not from the gene remainders", () => {
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const text = mountSwitcher(kind, "SAMEA_TEST").find(".arr-anchored").text();
        // Every fixture locus has an empty anchor rank list, so this is the muted branch throughout.
        expect(detail.anchor.arrangement_ranks).toHaveLength(0);
        if (detail.arrangements.membership_is_complete) {
          expect(text).toContain("has no gene at this locus");
        } else {
          expect(text).toContain("has no recorded neighbourhood at this locus");
        }
      }
    });

    it("⚠ `no_window` is the locus that proves the gene remainder cannot answer it", () => {
      // 64 members have no window here — and `over_cap` has 7, an order of magnitude fewer, yet both
      // are incomplete. A rule keyed on the gene count would have to pick a threshold; there isn't one.
      const noWindow = detailFor("no_window").arrangements;
      const overCap = detailFor("over_cap").arrangements;
      expect(noWindow.members_without_a_neighbourhood).toBeGreaterThan(
        overCap.members_without_a_neighbourhood,
      );
      expect(noWindow.membership_is_complete).toBe(false);
      expect(overCap.membership_is_complete).toBe(false);
    });

    it("counts the total in the sentence and the drawn rows in the coverage", () => {
      const detail = detailFor("over_cap");
      const head = mountSwitcher("over_cap").find(".arr-head").text();
      expect(head).toContain(`sits in ${detail.arrangements.total} neighbourhoods`);
      expect(head).toContain(`The ${detail.arrangements.listed.length} here account for`);
      expect(head).not.toContain("100% of its genes");
    });

    it("offers exactly the rows the track holds slots for, and no more", () => {
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const buttons = mountSwitcher(kind).findAll(".arr-opt");
        // The single-arrangement branch draws none; every other locus draws one per listed row.
        const expected = detail.arrangements.total <= 1 ? 0 : detail.arrangements.listed.length;
        expect(buttons).toHaveLength(expected);
      }
    });
  });

  describe("⭐ the locus card, on real bytes", () => {
    /**
     * ⭐ The floors and the measurable-locus counts come from the SPECIES response recorded in the
     * same pass, never from a stand-in. They used to be a hand-written `PROJECTIONS` constant, and a
     * constant cannot show the two endpoints disagreeing — which is the one thing a recorded pair is
     * for. ⚠ `similarity_baselines` is per representation and the two floors are an order of
     * magnitude apart, so a card handed the wrong one draws every tick in the wrong place.
     */

    it("⛔ every prevalence band on the wire is in the union the client switches on", () => {
      // `rare` was absent from `PrevalenceBand` while 30 % of loci carried it. TypeScript cannot see
      // that — the value arrives as a string — so this is the assertion that can.
      const seen = new Set<string>();
      for (const kind of CASES) {
        const detail = detailFor(kind);
        seen.add(detail.locus.prevalence_band);
        for (const row of detail.neighbour_display_rows) seen.add(row.prevalence_band);
      }
      expect(seen.size).toBeGreaterThan(1);
      for (const band of seen) {
        expect(PREVALENCE_BANDS).toContain(band);
      }
    });

    it("mounts the headline against every fixture locus without a warning", () => {
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const head = mount(LocusHeadline, {
          props: {
            locus: detail.locus,
            collectionGenomeCount: 100,
            separationMeasurableLocusCount: 12_104,
          },
        });
        expect(head.find(".locus-name").text()).toBe(detail.locus.display_name);
        expect(head.findAll(".tile")).toHaveLength(6);
        // ⛔ No tile may render `undefined` or `NaN`: both are what a missing field looks like once
        // it has been through `toFixed` or a percentage.
        for (const tile of head.findAll(".tile .v")) {
          expect(tile.text()).not.toContain("undefined");
          expect(tile.text()).not.toContain("NaN");
        }
      }
    });

    it("⭐ resolves every Pfam accession the card chips — no bare accessions survive", () => {
      // The reference is built from the accessions the response itself mentions, so a bare chip means
      // the join failed. This is the assertion that the version-stripping actually works on real data.
      let chipped = 0;
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const card = mount(SequenceDiversityCard, {
          props: {
            locus: detail.locus,
            families: detail.uniref50_families,
            symbols: detail.annotations["gene_symbol"] ?? [],
            pfamReference: detail.pfam_reference,
            listedArchitectureCount: (detail.annotations["pfam_architecture"] ?? []).length,
          },
        });
        for (const chip of card.findAll(".chip.pfam")) {
          chipped += 1;
          expect(chip.text()).not.toMatch(/^PF\d{5}$/);
          expect(chip.attributes("href")).toMatch(/^https:\/\/www\.ebi\.ac\.uk\/interpro\/entry\//);
        }
      }
      expect(chipped).toBeGreaterThan(0);
    });

    it("⛔ takes every family share over the LOCUS SIZE, so they need not sum to 100 %", () => {
      const detail = detailFor("over_cap");
      const card = mount(SequenceDiversityCard, {
        props: {
          locus: detail.locus,
          families: detail.uniref50_families,
          symbols: detail.annotations["gene_symbol"] ?? [],
          pfamReference: detail.pfam_reference,
          listedArchitectureCount: (detail.annotations["pfam_architecture"] ?? []).length,
        },
      });
      const listed = detail.uniref50_families.reduce((total, one) => total + one.gene_count, 0);
      expect(listed).toBeLessThanOrEqual(detail.locus.gene_count);
      const rows = card.findAll("tbody tr");
      expect(rows).toHaveLength(detail.uniref50_families.length);
    });

    it("mounts the similarity card and never prints a distance where a similarity belongs", () => {
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const card = mount(EmbeddingSimilarityCard, {
          props: {
            similarity: detail.locus.similarity,
            baselines: recorded.species.similarity_baselines,
            neighbours: detail.neighbour_display_rows,
            geneCount: detail.locus.gene_count,
            view: "median" as const,
          },
        });
        if (!card.find(".card").exists()) continue;
        // ⛔ Every rail is the SIMILARITY the column stores, printed as itself. Asserted as the
        // exact identity rather than as a plausible range: these numbers used to be DISTANCES that
        // the client converted with `1 − d`, and on real bytes a `nearest other cluster` sits at
        // 0.371, so any threshold that would catch a surviving subtraction also rejects real data.
        const expected = (["bacformer", "esm"] as const).flatMap((representation) => {
          const similarity = detail.locus.similarity[representation];
          if (similarity === null) return [];
          return [similarity.within_similarity, similarity.nearest_similarity]
            .filter((value): value is number => value !== null)
            .map((value) => value.toFixed(3));
        });
        expect(card.findAll(".pair:not(.sep) .val").map((node) => node.text())).toEqual(expected);
      }
    });

    it("⚠ renders a MEASURED 1.0 as a similarity of 1.000, not as a missing rail", () => {
      // ESM's within-cluster median over every gene pair really is 1.0 on real loci — the members
      // are that close — and a card testing truthiness rather than `!== null` would be fine, but
      // one testing `< 1` or rounding to a percentage would drop or flatten the rail entirely.
      const withOne = CASES.map(detailFor).find(
        (detail) => detail.locus.similarity.esm?.within_similarity === 1,
      );
      expect(withOne).toBeDefined();
      const card = mount(EmbeddingSimilarityCard, {
        props: {
          similarity: withOne!.locus.similarity,
          baselines: recorded.species.similarity_baselines,
          neighbours: withOne!.neighbour_display_rows,
          geneCount: withOne!.locus.gene_count,
          view: "median" as const,
        },
      });
      expect(card.findAll(".pair:not(.sep) .val").map((node) => node.text())).toContain("1.000");
    });

    it("⭐ every view renders on real bytes, and none of them prints undefined or NaN", () => {
      // The three views are not three pictures of one thing — only 18–34 % of flagged loci are
      // flagged by all three — so a view that is never mounted against real bytes is a view whose
      // fields were never checked to exist.
      for (const view of ["median", "weak", "own"] as const) {
        for (const kind of CASES) {
          const detail = detailFor(kind);
          const card = mount(EmbeddingSimilarityCard, {
            props: {
              similarity: detail.locus.similarity,
              baselines: recorded.species.similarity_baselines,
              neighbours: detail.neighbour_display_rows,
              geneCount: detail.locus.gene_count,
              view,
            },
          });
          expect(card.text()).not.toContain("undefined");
          expect(card.text()).not.toContain("NaN");
        }
      }
    });

    it("⛔⛔ prints NO percentile row where the own fraction is exactly 1.0", () => {
      // 88.5 % of loci sit at 1.0, so a midrank inside that tie block reads "p53" — *better than
      // half the catalogue* — when it means *tied with nearly all of it*. The fixture must exercise
      // both sides or this assertion is checking nothing.
      let saturated = 0;
      let below = 0;
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const card = mount(EmbeddingSimilarityCard, {
          props: {
            similarity: detail.locus.similarity,
            baselines: recorded.species.similarity_baselines,
            neighbours: detail.neighbour_display_rows,
            geneCount: detail.locus.gene_count,
            view: "own" as const,
          },
        });
        const ranks = card.findAll(".pair.sep").length;
        const shares = (["bacformer", "esm"] as const)
          .map((representation) => detail.locus.similarity[representation]?.own_neighbour_fraction ?? null)
          .filter((share): share is number => share !== null);
        saturated += shares.filter((share) => share >= 1).length;
        below += shares.filter((share) => share < 1).length;
        expect(ranks).toBe(shares.filter((share) => share < 1).length);
        // ⛔ …and a share below 1 never rounds up to "100%" on the row above.
        for (const value of card.findAll(".pair:not(.sep) .val").map((node) => node.text())) {
          if (value === "100%") continue;
          expect(value).toMatch(/^(<100%|\d{1,2}\.\d%)$/);
        }
      }
      expect(saturated).toBeGreaterThan(0);
      expect(below).toBeGreaterThan(0);
    });
  });

  describe("⭐ the five nearest loci, on real bytes", () => {
    it("⛔ resolves every nearest locus — the set the fan-out nearly missed entirely", () => {
      // `catalogue_ordinal` is a THIRD address space alongside arrangement slot ordinals and
      // marginal occupant locus ids, and the retired `nearest_locus_ordinals` was not collected at
      // all: the list would have had a similarity with no name beside it.
      let checked = 0;
      for (const kind of CASES) {
        const detail = detailFor(kind);
        const byOrdinal = new Set(detail.neighbour_display_rows.map((row) => row.catalogue_ordinal));
        for (const representation of ["bacformer", "esm"] as const) {
          for (const nearest of detail.locus.similarity[representation]?.nearest_loci ?? []) {
            checked += 1;
            expect(byOrdinal.has(nearest.catalogue_ordinal)).toBe(true);
          }
        }
      }
      expect(checked).toBeGreaterThan(0);
    });

    it("⚠ the list is RAGGED on real data — a rank is never a position in an array", () => {
      // If every locus carried five, a client indexing by rank would pass every test here and then
      // name the wrong locus on the first short shortlist it met.
      const lengths = CASES.flatMap((kind) =>
        (["bacformer", "esm"] as const).map(
          (representation) => detailFor(kind).locus.similarity[representation]?.nearest_loci.length ?? 0,
        ),
      );
      expect(lengths.some((length) => length > 0 && length < 5)).toBe(true);
      for (const kind of CASES) {
        for (const representation of ["bacformer", "esm"] as const) {
          const nearest = detailFor(kind).locus.similarity[representation]?.nearest_loci ?? [];
          // ⚠ 1-based and strictly ranked, which is what a client may rely on and nothing else.
          expect(nearest.map((one) => one.rank)).toEqual(nearest.map((_, index) => index + 1));
        }
      }
    });

    it("⚠ the two representations really do name different loci on real data", () => {
      // The reason both are stacked rather than tabbed. If they agreed here, every claim about the
      // two disagreeing would be proving nothing.
      let differing = 0;
      for (const kind of CASES) {
        const similarity = detailFor(kind).locus.similarity;
        const context = new Set(
          (similarity.bacformer?.nearest_loci ?? []).map((one) => one.catalogue_ordinal),
        );
        const sequence = new Set(
          (similarity.esm?.nearest_loci ?? []).map((one) => one.catalogue_ordinal),
        );
        if ([...context].some((ordinal) => !sequence.has(ordinal))) differing += 1;
      }
      expect(differing).toBeGreaterThan(0);
    });

    it("⭐ rank 1's similarity IS `nearest_similarity`, which is what lets the list sit under it", () => {
      // The ingest asserts this too; asserting it on the wire is what proves the two travelled
      // together. A disagreement can only mean two different runs.
      let checked = 0;
      for (const kind of CASES) {
        for (const representation of ["bacformer", "esm"] as const) {
          const similarity = detailFor(kind).locus.similarity[representation];
          const first = similarity?.nearest_loci[0];
          if (similarity == null || first === undefined) continue;
          checked += 1;
          expect(first.cross_similarity).toBeCloseTo(similarity.nearest_similarity as number, 6);
        }
      }
      expect(checked).toBeGreaterThan(0);
    });
  });

  describe("⭐ the Function tab, on real bytes from its own endpoint", () => {
    const NAMESPACES = [
      "molecular_function",
      "biological_process",
      "cellular_component",
    ] as const satisfies readonly GeneOntologyNamespace[];

    function functionFor(kind: (typeof CASES)[number]): FunctionResponse {
      return recorded.loci[kind]!.function;
    }

    it("⛔⛔ names each GO namespace the way the COVERAGE block does — the API spoke two dialects", () => {
      // Found by building this tab: `gene_ontology_namespace` came back as the stored 0/1/2 while
      // `coverage.go_annotated_gene_count` was keyed by name, in the same response. A client grouping
      // entries by an integer renders three GO cards with EMPTY term lists under coverage lines that
      // promise otherwise — and no existing test could see it, because nothing read the field.
      let seen = 0;
      for (const kind of CASES) {
        const block = functionFor(kind);
        expect(Object.keys(block.coverage.go_annotated_gene_count).sort()).toEqual([...NAMESPACES].sort());
        for (const entry of block.annotations.gene_ontology_slim ?? []) {
          expect(NAMESPACES).toContain(entry.gene_ontology_namespace);
          seen += 1;
        }
      }
      expect(seen).toBeGreaterThan(0);
    });

    it("⛔⛔ the response carries NO GO verdict — retired, on real bytes", () => {
      // David, 2026-10-05. `worst_relation` is Pfam's domain-architecture comparator; GO passed it
      // no ontology, and `goslim_metagenomics` keeps a term and its own parent as siblings, so
      // `gumC` scored "classes differ" for plasma membrane vs membrane. ⚠ The recorded bytes are
      // the only thing that can show the route still emitting it, which is why this is here and
      // not only in the unit tests.
      for (const kind of CASES) {
        expect("go_verdicts" in (functionFor(kind) as Record<string, unknown>)).toBe(false);
      }
    });

    it("⛔ no coverage count can exceed the locus, and each matches the locus response", () => {
      // Two endpoints describing one locus. Nothing but the pair can show them disagreeing.
      for (const kind of CASES) {
        const block = functionFor(kind);
        const size = block.coverage.gene_count;
        expect(size).toBe(recorded.loci[kind]!.response.locus.gene_count);
        for (const count of [
          block.coverage.cog_annotated_gene_count,
          block.coverage.ec_annotated_gene_count,
          block.coverage.kegg_annotated_gene_count,
          ...NAMESPACES.map((namespace) => block.coverage.go_annotated_gene_count[namespace]),
        ]) {
          expect(count).toBeGreaterThanOrEqual(0);
          expect(count).toBeLessThanOrEqual(size);
        }
      }
    });

    it("⛔ a namespace with CLASSES has the coverage to support them", () => {
      // ⚠ Replaces a cross-check between the stored verdict and the stored count. The verdict is
      // retired, but the pair it guarded is still worth checking from the other side: a namespace
      // that lists classes must have genes carrying them, or the table stands over nothing.
      let checked = 0;
      for (const kind of CASES) {
        const block = functionFor(kind);
        for (const namespace of NAMESPACES) {
          const annotated = block.coverage.go_annotated_gene_count[namespace];
          const classes = (block.annotations.gene_ontology_slim ?? []).filter(
            (entry) => entry.gene_ontology_namespace === namespace,
          );
          if (classes.length > 0) expect(annotated).toBeGreaterThan(0);
          for (const entry of classes) expect(entry.gene_count).toBeLessThanOrEqual(annotated);
          checked += 1;
        }
      }
      expect(checked).toBe(CASES.length * NAMESPACES.length);
    });

    it("⭐⭐ a donor whose GO claim is a NAMESPACE ROOT suggests nothing — on real bytes", () => {
      // ⛔ The clearest case the GO closure produces, and it turned up in the recorded bytes by
      // itself. The kp `over_cap` locus has a 0.986 neighbour — a strong tier where GO agreement is
      // 91.8 % — whose own top GO term is `GO:0008150`, the BIOLOGICAL PROCESS ROOT. It is an
      // ancestor of every biological-process term, so it says only "this protein does something".
      // Before the claim was closed and the roots removed, the page would have offered it as this
      // node's suggested function.
      //
      // ⚠ Empty means "states nothing comparable", which is neither agreement nor disagreement — so
      // the refusal is correct even though the tier and the rate would both have allowed a call.
      const go = recorded.loci.over_cap?.function.inference.vocabularies.find(
        (v) => v.annotation_kind === "gene_ontology_slim",
      );
      if (!go?.candidate) return; // ecoli's `over_cap` has no GO neighbour to refuse
      if (go.candidate.donor.term !== "GO:0008150") return;
      expect(go.candidate.level).toBeNull();
      expect(go.candidate.value).toBeNull();
      expect(go.candidate.tier).toBe("0.98-0.99");
    });

    it("⭐ EC and KEGG arrive NAMED, on real bytes", () => {
      // ⚠ This assertion is the inverse of the one it replaces. Until 2026-10-05 it read "a KEGG
      // row is present and NEVER named" and the EC row carried `name: null` too — the first on a
      // licence that is real but was never resolved either way, the second on one that does not
      // exist (ExPASy ENZYME is CC BY 4.0, not CC BY-ND). Both are vendored now and joined
      // server-side, so a regression that dropped the join would show here as a null.
      const ec = CASES.flatMap((kind) => functionFor(kind).annotations.ec_number ?? []);
      const kegg = CASES.flatMap((kind) => functionFor(kind).annotations.kegg_orthology ?? []);
      expect(ec.length).toBeGreaterThan(0);
      expect(kegg.length).toBeGreaterThan(0);
      for (const entry of [...ec, ...kegg]) {
        expect(entry.name).not.toBeNull();
        expect(entry.name!.length).toBeGreaterThan(0);
        // ⛔ A name, not the code echoed back — which is what a half-wired join would produce.
        expect(entry.name).not.toBe(entry.term);
      }
    });

    it("renders the rich locus with all three namespaces, EC and KEGG", () => {
      const block = functionFor("function_rich");
      const tab = mount(FunctionTab, {
        props: {
          displayName: recorded.loci.function_rich!.response.locus.display_name,
          locusLabel: recorded.loci.function_rich!.label,
          block,
          status: "ready" as const,
        },
      });
      const headings = tab.findAll(".sub-head").map((node) => node.text());
      expect(headings).toContain("GO — molecular function");
      expect(headings).toContain("GO — biological process");
      expect(headings).toContain("GO — cellular component");
      expect(headings).toContain("EC and KEGG");
      // ⛔ and the GO tables are POPULATED — the whole point of the dialect bug above
      const populated = tab.findAll(".card").filter((card) => card.find("tbody tr").exists());
      expect(populated.length).toBeGreaterThanOrEqual(4);
    });

    it("⭐ the inferred-function card renders on real bytes, with its MEASURED rate and the walk", () => {
      // ⛔ The case is useless unless it still exercises what it was chosen for, so assert that
      // first: no COG of its own, and a donor found past rank 1. A regenerate that quietly picked an
      // annotated locus would otherwise leave this test passing while testing nothing.
      const block = functionFor("inferred_function");
      const cog = block.inference.vocabularies.find(
        (entry) => entry.annotation_kind === "cog_orthogroup",
      )!;
      expect(cog.own).toBeNull();
      expect(cog.candidate).not.toBeNull();
      expect(cog.candidate!.rank).toBeGreaterThanOrEqual(2);
      expect(cog.walk.some((step) => !step.carries_annotation)).toBe(true);

      const tab = mount(FunctionTab, {
        props: {
          displayName: recorded.loci.inferred_function!.response.locus.display_name,
          locusLabel: recorded.loci.inferred_function!.label,
          block,
          status: "ready" as const,
        },
      });
      const card = tab.find(".infer");
      expect(card.exists()).toBe(true);
      expect(card.find(".infer-tag").text()).toBe("inferred");
      // the measured rate, with its n — both, or the rate means nothing
      const cell = cog.candidate!.calibration!;
      expect(card.find(".infer-rate").text()).toContain(
        `${(cell.agreement! * 100).toFixed(1)}%`,
      );
      expect(card.find(".infer-rate").text()).toContain(cell.pairs.toLocaleString());
      // and the ranks that carried nothing are named, not skipped past
      expect(card.text()).toContain("carried no COG");
      // ⛔ the ladder shipped with the response, and the row the suggestion was read from
      const used = card.findAll(".infer-row-used");
      expect(used).toHaveLength(1);
      expect(used[0]!.text()).toContain(cell.tier);
    });

    it("⛔ the calibration shipped with every response is self-consistent, and never a bare rate", () => {
      for (const kind of CASES) {
        const block = functionFor(kind);
        for (const ladder of Object.values(block.calibration)) {
          expect(ladder.cells.length).toBeGreaterThan(0);
          for (const cell of ladder.cells) {
            expect(cell.pairs).toBeGreaterThan(0);
            if (cell.pairs < ladder.min_pairs) {
              // ⛔ below the floor there is a COUNT and no rate — a rate off a handful of pairs
              // would be read as a measurement
              expect(cell.agreement).toBeNull();
              expect(cell.lift).toBeNull();
            } else {
              expect(cell.agreement).not.toBeNull();
              expect(cell.interval_low!).toBeLessThanOrEqual(cell.agreement!);
              expect(cell.interval_high!).toBeGreaterThanOrEqual(cell.agreement!);
              expect(cell.lift!).toBeGreaterThan(1);
            }
          }
        }
      }
    });

    it("⚠ every GO row carries an accession AND a readable class name", () => {
      // `GO:0016020` alone is unreadable; "membrane" alone is unlookupable.
      const entries = CASES.flatMap(
        (kind) => (functionFor(kind).annotations.gene_ontology_slim ?? []) as readonly AnnotationEntry[],
      );
      expect(entries.length).toBeGreaterThan(0);
      for (const entry of entries) {
        expect(entry.term).toMatch(/^GO:\d{7}$/);
        expect(entry.name).not.toBeNull();
      }
    });
  });

  describe("⭐ the Sequence tab, on real bases from a real GFF", () => {
    const minus = recorded.sequences.minus_strand!;
    const absent = recorded.sequences.no_gene_here!;

    function mountSequence(record: typeof minus) {
      return mount(SequenceTab, {
        props: {
          displayName: "x",
          locusLabel: record.locus_label,
          sampleId: record.sample_id,
          response: record.response,
          status: "ready" as const,
        },
      });
    }

    it("⛔⛔ a MINUS-strand gene's upstream flank comes from the HIGHER coordinates", () => {
      // The bug of record, on real bases. `load_meta_flanks` orients on strand and `app.js:4376-4386`
      // states the convention; slicing `start - flank` unconditionally returns the DOWNSTREAM flank
      // for about half of all genes and renders as a perfectly plausible 100 bases of DNA.
      const gene = minus.response.genes[0]!;
      expect(gene.strand).toBe("-");
      expect(gene.upstream_flank_span![0]).toBeGreaterThan(gene.end_position);
      expect(gene.downstream_flank_span![1]).toBeLessThan(gene.start_position);
    });

    it("⚠ and the page SAYS the flanks were reverse-complemented", () => {
      // Mislabelling a correct sequence is worse than showing the wrong one: the reader acts on it.
      const provenances = mountSequence(minus)
        .findAll(".seq-note")
        .map((node) => node.text());
      expect(provenances[0]).toContain("reverse-complemented");
    });

    it("⛔ the protein is what the span implies, stop codon included", () => {
      const gene = minus.response.genes[0]!;
      expect(gene.coding_sequence.length).toBe(gene.end_position - gene.start_position + 1);
      expect(gene.coding_sequence.length).toBe(gene.length_nt);
      // The span INCLUDES the stop codon, so the protein is `length / 3 - 1` residues.
      expect(gene.protein_sequence.length).toBe(gene.length_nt / 3 - 1);
      expect(gene.protein_sequence).not.toContain("*");
    });

    it("⭐ the sliced GC reproduces the column ingest wrote", () => {
      // Two computations of one number, from the same coordinates at different times.
      const gene = minus.response.genes[0]!;
      expect(gene.gc_percent).toBeCloseTo(gene.stored_gc_percent!, 9);
    });

    it("⛔ the bases are DNA and nothing else", () => {
      const gene = minus.response.genes[0]!;
      for (const sequence of [
        gene.coding_sequence,
        gene.upstream_flank_sequence,
        gene.downstream_flank_sequence,
      ]) {
        expect(sequence).toMatch(/^[ACGTN]*$/);
      }
      expect(gene.protein_sequence).toMatch(/^[ACDEFGHIKLMNPQRSTVWY*]+$/);
    });

    it("⛔⛔ a genome with NO gene here renders as an answer, not as a failure", () => {
      expect(absent.response.genes).toEqual([]);
      const tab = mountSequence(absent);
      expect(tab.text()).toContain("has no gene at this locus");
      expect(tab.find(".pop-error").exists()).toBe(false);
    });

    it("⚠ a flank span describes exactly the string that came with it", () => {
      for (const record of [minus]) {
        for (const gene of record.response.genes) {
          if (gene.upstream_flank_span !== null) {
            const [from, to] = gene.upstream_flank_span;
            expect(to - from + 1).toBe(gene.upstream_flank_sequence.length);
          } else {
            expect(gene.upstream_flank_sequence).toBe("");
          }
        }
      }
    });
  });

  describe("⭐ the page shell, on real bytes", () => {
    it("⛔ the GENE census partitions every modelled gene, and the per-genome line says so", () => {
      const catalogue = recorded.species;
      const genes = Object.values(catalogue.prevalence_gene_census).reduce((total, value) => total + value, 0);
      expect(genes).toBe(catalogue.pangenome.gene_count);
      const text = mount(PangenomeCensus, { props: { catalogue } }).text();
      expect(text).toContain(`${Math.round(genes / catalogue.pangenome.genome_count).toLocaleString()} genes`);
      expect(text).not.toContain("NaN");
    });

    it("⛔ every example chip is NAMED — never a bare number", () => {
      const view = mount(NavigatingView, { props: { examples: recorded.species.example_locus_rows } });
      const chips = view.findAll(".chip");
      expect(chips.length).toBe(recorded.species.example_loci.length);
      expect(chips.length).toBeGreaterThan(0);
      for (const chip of chips) expect(chip.text()).not.toMatch(/^\d+$/);
    });

    it("⛔ the footer's residual counts agree with the lists their own endpoint returns", () => {
      const headline = recorded.species.audit_headline;
      const count = (key: string) => (typeof headline[key] === "number" ? (headline[key] as number) : 0);
      expect(recorded.residuals.grouped_on_context_alone).toHaveLength(
        count("synteny_only_n_clusters") + count("no_homology_n_clusters"),
      );
      expect(recorded.residuals.pfam_conflicts).toHaveLength(count("pfam_conflict_n_clusters"));
      const footer = mount(SiteFooter, { props: { catalogue: recorded.species } });
      const summaries = footer.findAll("summary").map((node) => node.text());
      expect(summaries).toContain(
        `Every locus grouped on context alone (${recorded.residuals.grouped_on_context_alone.length})`,
      );
      expect(summaries).toContain(`Every Pfam conflict (${recorded.residuals.pfam_conflicts.length})`);
    });

    it("⛔⛔ the footer carries the KEGG licence flag on EVERY page", () => {
      // This is a GATE on deployment, not a courtesy. KEGG requires an academic service provider
      // licence from anyone offering a service, and this build names KO descriptions — so the
      // unanswered question has to travel with the page. David, 2026-10-05: *"Just flag it on site.
      // Will look at it later. This is a prototype only!"* A change that drops this paragraph while
      // the names are still served is the one that quietly turns a prototype into a publication.
      const footer = mount(SiteFooter, { props: { catalogue: recorded.species } }).text();
      expect(footer).toContain("Internal build — not for publication");
      expect(footer).toContain("academic service provider licence");
      expect(footer).toContain("not yet confirmed");
    });

    it("⚠ the footer quotes the audit headline and prints no unformatted number", () => {
      const footer = mount(SiteFooter, { props: { catalogue: recorded.species } }).text();
      expect(footer).toContain("clusters graded");
      expect(footer).not.toContain("undefined");
      expect(footer).not.toContain("NaN");
    });

    it("⭐ the footer dates BOTH the catalogue and the bundle, so 'is this current?' is answerable", () => {
      // David, 2026-10-03: "the numbers otherwise change and I can't find the current version
      // without asking each time." The compose stack serves a BUILT IMAGE against its own database,
      // so "the app is up" and "the app is current" are different facts and an image from last week
      // looks identical to a fresh one. Only the page can tell them apart, and only if it says both.
      const withIngest = {
        ...recorded.species,
        pangenome: { ...recorded.species.pangenome, built_at: null, ingested_at: "2026-09-25T08:35:06.477248+00:00" },
      };
      const footer = mount(SiteFooter, { props: { catalogue: withIngest } }).text();
      // ⛔ `built_at` is NULL on every loaded catalogue, so the row that used to read "built —" has
      // to fall back to the timestamp that IS populated — and say which one it is showing.
      expect(footer).toContain("catalogue loaded");
      expect(footer).toContain("2026-09-25 08:35 UTC");
      expect(footer).not.toContain("built —");
      // …and when the bundle itself was built, which no catalogue field can supply.
      expect(footer).toContain("page built");
      expect(footer).toContain("2026-01-01 00:00 UTC");
    });

    it("⚠ prefers `built_at` where it exists — it is the stronger claim than when we loaded it", () => {
      const withBuilt = {
        ...recorded.species,
        pangenome: { ...recorded.species.pangenome, built_at: "2026-08-25", ingested_at: "2026-09-25T08:35:06Z" },
      };
      const footer = mount(SiteFooter, { props: { catalogue: withBuilt } }).text();
      expect(footer).toContain("catalogue built");
      expect(footer).toContain("2026-08-25");
      expect(footer).not.toContain("catalogue loaded");
    });

    it("⛔ every search row carries a product field, and every band is one the client knows", () => {
      const hits = recorded.search.ligase.hits;
      expect(hits.length).toBeGreaterThan(5);
      for (const hit of hits) {
        expect("best_product" in hit).toBe(true);
        expect(PREVALENCE_BANDS).toContain(hit.prevalence_band);
      }
      // A mid-string product match is the point of substring search; the fixture must exercise it.
      expect(hits.some((hit) => (hit.best_product ?? "").toLowerCase().includes("ligase"))).toBe(true);
    });

    it("⛔ the anchor list is cut honestly, and a genome's arranged loci never exceed its loci", () => {
      const list = recorded.genomes;
      expect(list.genomes).toHaveLength(5);
      expect(list.truncated).toBe(true);
      expect(list.matched_genome_count).toBe(list.genome_count);
      expect(list.genome_count).toBe(recorded.species.pangenome.genome_count);
      for (const genome of list.genomes) {
        // Loci where it sits in a recorded neighbourhood are a SUBSET of loci where it has a gene.
        expect(genome.arrangement_locus_count).toBeLessThanOrEqual(genome.locus_count);
        // …and per LOCUS: a genome at ρ > 1 counted per gene would exceed the catalogue.
        expect(genome.locus_count).toBeLessThanOrEqual(recorded.species.pangenome.locus_count);
      }
      const ordinals = list.genomes.map((genome) => genome.collection_genome_ordinal);
      expect(ordinals).toEqual([...ordinals].sort((a, b) => a - b));
    });

    it("⚠ a residual row's family and architecture counts are counts or null, never -1", () => {
      for (const row of [...recorded.residuals.grouped_on_context_alone, ...recorded.residuals.pfam_conflicts]) {
        expect(row.pfam_architecture_count === null || row.pfam_architecture_count >= 0).toBe(true);
        expect(PREVALENCE_BANDS).toContain(row.prevalence_band);
      }
    });
  });
});
