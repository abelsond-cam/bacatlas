/**
 * @vitest-environment jsdom
 *
 * ⛔ **Regressions from the adversarial review of the page shell** (2026-09-18) — eleven defects, each
 * rendering plausibly, each reproduced before it was fixed. Every test here failed on the code as
 * first written and names the behaviour the reader saw.
 */
import { asCatalogueKey } from "@/api/types";
import { flushPromises, mount } from "@vue/test-utils";
import { createPinia, setActivePinia } from "pinia";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { nextTick } from "vue";

import { failure, success } from "@/api/result";
import App from "@/App.vue";
import AnchorGenomeBox from "@/components/anchor/AnchorGenomeBox.vue";
import LocusSearchBox from "@/components/search/LocusSearchBox.vue";
import { useAnchorGenomeStore } from "@/stores/anchorGenomeStore";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { useViewTabStore } from "@/stores/viewTabStore";

const api = vi.hoisted(() => ({
  fetchSpeciesList: vi.fn(),
  fetchCatalogues: vi.fn(),
  fetchSpeciesCatalogue: vi.fn(),
  fetchLocus: vi.fn(),
  fetchArrangementPage: vi.fn(),
  fetchLocusFunction: vi.fn(),
  fetchGeneSequence: vi.fn(),
  searchLoci: vi.fn(),
  fetchGenomes: vi.fn(),
  fetchAuditResiduals: vi.fn(),
  fetchProjectedGenomes: vi.fn(),
}));
// ⚠ `anchorQuery` is deliberately NOT stubbed: it is the one place the anchor's kind becomes a
// query parameter, and a stub would stop this suite noticing if a projected genome were ever
// fetched as `anchor=` — which 404s, and reads on the page as "this genome has nothing here".
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  ...api,
}));

function detail(label: string) {
  return {
    locus: { label, display_name: `name${label}` },
    anchor: { is_anchored: false, arrangement_ranks: [] },
    arrangements: { listed: [], total: 0 },
    offsets: [],
    intergenic_gaps: [],
    neighbour_display_rows: [],
  } as never;
}

/** The shell's children are not under test here — only App's own address handling. */
const stubs = {
  HoverTip: true,
  SiteHeader: true,
  SpeciesStrip: true,
  TrackPanel: true,
  ViewTabs: true,
  FunctionView: true,
  LocusEvidenceView: true,
  NavigatingView: true,
  SequenceView: true,
  SiteFooter: true,
};

const SPECIES = [
  // ⚠ `catalogue_key` is what a bare `?species=` resolves TO, and the address is pinned to it.
  { key: "ecoli", published: true, scientific_name: "Escherichia coli", catalogue_key: "ecoli-nuna4" },
  { key: "kp", published: true, scientific_name: "Klebsiella pneumoniae", catalogue_key: "kp-nuna4" },
];

const CATALOGUES = [
  { key: "ecoli-nuna4", species: { key: "ecoli", scientific_name: "Escherichia coli" },
    model: { key: "nuna4", label: null, step_count: 4, exclusivity_form: "damped_exclusion" },
    is_default: true, genome_count: 100, gene_count: 1, locus_count: 1, run_id: "r4" },
  { key: "ecoli-nuna5", species: { key: "ecoli", scientific_name: "Escherichia coli" },
    model: { key: "nuna5", label: null, step_count: 5, exclusivity_form: "exclusion" },
    is_default: false, genome_count: 100, gene_count: 1, locus_count: 1, run_id: "r5" },
  { key: "kp-nuna4", species: { key: "kp", scientific_name: "Klebsiella pneumoniae" },
    model: { key: "nuna4", label: null, step_count: 4, exclusivity_form: "damped_exclusion" },
    is_default: true, genome_count: 100, gene_count: 1, locus_count: 1, run_id: "k4" },
];

beforeEach(() => {
  setActivePinia(createPinia());
  for (const mocked of Object.values(api)) mocked.mockReset();
  api.fetchSpeciesList.mockResolvedValue(success({ species: SPECIES }));
  api.fetchCatalogues.mockResolvedValue(success({ catalogues: CATALOGUES }));
  api.fetchSpeciesCatalogue.mockImplementation(async (key: string) =>
    success({
      // ⚠ The shell names the catalogue it IS, which is what the address is pinned to.
      species: { key: key.split("-")[0], scientific_name: key },
      pangenome: { genome_count: 100, catalogue_key: key },
      landing_locus: "2811",
      example_locus_rows: [],
    }),
  );
  api.fetchLocus.mockImplementation(async (_species: string, label: string) => success(detail(label)));
  // ⛔ App's watcher fires `void own.loadPlacedGenomes(…)` as soon as the species list resolves, and the
  // promise is deliberately discarded there — so an unarmed mock returns `undefined`, `result.ok` throws,
  // and the only trace is an unhandled rejection with no assertion to catch it. Six of them, one per
  // `mount(App)` test, red in CI since 2026-09-23.
  // ⚠ Empty `genomes` keeps these tests' meaning: they assert on the hash, the history and which
  // catalogue was fetched, never on placed genomes.
  api.fetchProjectedGenomes.mockResolvedValue(success({ species_key: "ecoli", genomes: [], caveat: "" }));
});
afterEach(() => vi.useRealTimers());

async function settle(): Promise<void> {
  for (let round = 0; round < 6; round += 1) {
    await flushPromises();
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
}

describe("⛔ the address", () => {
  it("landing on the bare URL REPLACES its entry, so the first Back leaves rather than doing nothing", async () => {
    window.history.replaceState(null, "", "/?species=ecoli");
    const before = window.history.length;
    mount(App, { global: { stubs } });
    await settle();
    expect(window.location.hash).toBe("#2811");
    expect(useLocusNavigationStore().route?.label).toBe("2811");
    expect(window.history.length - before).toBe(0);
  });

  it("an unknown ?species= is SAID, and no other catalogue is opened in its place", async () => {
    window.history.replaceState(null, "", "/?species=nosuch#1098");
    const page = mount(App, { global: { stubs } });
    await settle();
    expect(api.fetchSpeciesCatalogue).not.toHaveBeenCalled();
    expect(api.fetchLocus).not.toHaveBeenCalled();
    expect(page.text()).toContain("No species “nosuch” is published on this server");
  });

  it("⚠ a wrong-CASE species is forgiven, and the address PINNED to the catalogue it resolved to", async () => {
    window.history.replaceState(null, "", "/?species=KP#1098");
    mount(App, { global: { stubs } });
    await settle();
    expect(api.fetchSpeciesCatalogue).toHaveBeenCalledWith("kp-nuna4");
    // ⭐ The address now names the CLUSTERING, not a default that can move underneath it. Without
    // this, `?species=kp#1098` opens a different gene the day the default moves — locus labels are
    // model-private — and it renders perfectly (David, 2026-09-26).
    const address = new URL(window.location.href).searchParams;
    expect(address.get("catalogue")).toBe("kp-nuna4");
    expect(address.get("species")).toBeNull();
    expect(api.fetchLocus.mock.calls[0]?.slice(0, 2)).toEqual(["kp-nuna4", "1098"]);
  });

  it("⛔ a ?catalogue= pins that model, and is NOT checked against the offered menu", async () => {
    // A catalogue can be loaded but staged — that is how a model is reviewed before anyone is shown
    // it — so gating the load on `GET /catalogues` would make a deliberate link unopenable.
    window.history.replaceState(null, "", "/?catalogue=ecoli-nuna5#1098");
    mount(App, { global: { stubs } });
    await settle();
    expect(api.fetchSpeciesCatalogue).toHaveBeenCalledWith("ecoli-nuna5");
    expect(api.fetchLocus.mock.calls[0]?.slice(0, 2)).toEqual(["ecoli-nuna5", "1098"]);
  });

  it("⚠ ?catalogue= WINS over a disagreeing ?species=, and the loser is dropped", async () => {
    window.history.replaceState(null, "", "/?species=kp&catalogue=ecoli-nuna5");
    mount(App, { global: { stubs } });
    await settle();
    expect(api.fetchSpeciesCatalogue).toHaveBeenCalledWith("ecoli-nuna5");
    expect(api.fetchSpeciesCatalogue).not.toHaveBeenCalledWith("kp-nuna4");
  });
});

describe("⛔ the `r` probe", () => {
  it("that lands AFTER the reader moved on changes nothing — it no longer yanks them back", async () => {
    window.history.replaceState(null, "", "/?species=ecoli#100");
    let releaseProbe: (value: unknown) => void = () => undefined;
    api.fetchLocus.mockImplementation(async (_species: string, label: string) =>
      label === "200r" ? new Promise((resolve) => (releaseProbe = resolve)) : success(detail(label)),
    );
    mount(App, { global: { stubs } });
    await settle();
    const navigation = useLocusNavigationStore();
    window.location.hash = "#200r";
    await settle();
    await navigation.navigateTo("400");
    await settle();
    releaseProbe(failure("not_found", "no locus '200r'", 404));
    await settle();
    expect(navigation.route?.label).toBe("400");
    expect(window.location.hash).toBe("#400");
  });

  it("asks about a whole-string 404 ONCE — Back onto a reversed entry does not re-probe", async () => {
    api.fetchLocus.mockImplementation(async (_species: string, label: string) =>
      label === "12r" ? failure("not_found", "no", 404) : success(detail(label)),
    );
    const navigation = useLocusNavigationStore();
    navigation.setCatalogue(asCatalogueKey("ecoli-nuna4"));
    await navigation.openHash("#12r");
    await navigation.navigateTo("40");
    await navigation.openHash("#12r");
    expect(api.fetchLocus.mock.calls.filter((call) => call[1] === "12r")).toHaveLength(1);
    expect(navigation.route).toEqual({ label: "12", direction: "reversed" });
  });

  it("⛔ that FAILED is re-asked by Try again, instead of reporting the link's locus as missing", async () => {
    const navigation = useLocusNavigationStore();
    navigation.setCatalogue(asCatalogueKey("ecoli-nuna4"));
    api.fetchLocus.mockResolvedValueOnce(failure("network", "offline"));
    api.fetchLocus.mockResolvedValueOnce(failure("network", "offline"));
    await navigation.openHash("#12r");
    expect(navigation.view.status).toBe("failed");
    // The network is back: `12r` does not exist, `12` does — the link meant 12, reversed.
    api.fetchLocus.mockImplementation(async (_species: string, label: string) =>
      label === "12r" ? failure("not_found", "no", 404) : success(detail(label)),
    );
    await navigation.retry();
    expect(navigation.route).toEqual({ label: "12", direction: "reversed" });
    expect(navigation.view.status).toBe("ready");
  });
});

describe("⛔ the trail and the tabs", () => {
  it("a locus that does not exist leaves no dead crumb", async () => {
    const navigation = useLocusNavigationStore();
    navigation.setCatalogue(asCatalogueKey("ecoli-nuna4"));
    await navigation.navigateTo("1");
    api.fetchLocus.mockResolvedValueOnce(failure("not_found", "no locus", 404));
    await navigation.navigateTo("99999");
    expect(navigation.trail).toEqual(["1"]);
  });

  it("the FIRST successful draw after a dead link brings the reader to the evidence", async () => {
    const navigation = useLocusNavigationStore();
    const tabs = useViewTabStore();
    navigation.setCatalogue(asCatalogueKey("ecoli-nuna4"));
    api.fetchLocus.mockResolvedValueOnce(failure("not_found", "no locus", 404));
    await navigation.navigateTo("99999");
    tabs.showView("navigating");
    await navigation.navigateTo("48");
    await nextTick();
    expect(tabs.view).toBe("locus");
  });
});

function hit(label: string, name: string) {
  return { label, display_name: name, best_product: null, gene_count: 1, genome_count: 1, prevalence_band: "core", rank_band: 1 };
}

describe("⛔ Enter acts on what is in the box NOW", () => {
  it("search: typing on and pressing Enter before the new answer lands chooses for the NEW text", async () => {
    vi.useFakeTimers();
    api.searchLoci.mockImplementation(async (_species: string, query: string) =>
      query === "rfa"
        ? success({ query: "rfa", mode: "substring", truncated: false, hits: [hit("10", "rfaC"), hit("48", "rfaL")] })
        : success({ query, mode: "substring", truncated: false, hits: [hit("48", "rfaL")] }),
    );
    const box = mount(LocusSearchBox, { props: { speciesKey: "ecoli", collectionGenomeCount: 100 } });
    await box.find("input").setValue("rfa");
    await vi.advanceTimersByTimeAsync(200);
    await flushPromises();
    await box.find("input").setValue("rfaL");
    await box.find("input").trigger("keydown", { key: "Enter" });
    await flushPromises();
    expect(box.emitted("go")).toEqual([["48"]]);
  });

  it("anchor: Enter over 'no genome matches' does NOT clear the reader's anchor", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setAnchor("SAMEA1");
    api.fetchGenomes.mockImplementation(async (_species: string, query: string) =>
      success({ species_key: "ecoli", query, genome_count: 100, matched_genome_count: 0, truncated: false, genomes: [] }),
    );
    const box = mount(AnchorGenomeBox, { props: { speciesKey: "ecoli", placement: "track" } });
    await box.find("input").trigger("focus");
    await box.find("input").setValue("SAMEA9999");
    await box.find("input").trigger("input");
    await flushPromises();
    await box.find("input").trigger("keydown", { key: "Enter" });
    expect(anchor.sampleId).toBe("SAMEA1");
  });
});
