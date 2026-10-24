/**
 * @vitest-environment jsdom
 */
/**
 * The anchor picker's TWO LISTS — the tabs, the census behind them, and the two empty states.
 *
 * ⭐ **Why two lists at all.** The cohort is two populations with different properties: ~4,000
 * complete long-read NCBI genomes, short enough to browse, and ~75,000 AllTheBacteria short-read
 * assemblies, which need the search. One combobox over both serves neither — the complete genomes
 * are the ones a reader anchors deliberately and would never find by typing, and they are ~5 % of
 * the list.
 *
 * ⛔ **The assertion this file exists for is the pair of empty states.** *"This catalogue holds no
 * complete genomes"* and *"no genome matches what you typed"* are different claims, and only the
 * first is true of every catalogue loaded today. A reader shown the second keeps typing at a list
 * that can never answer. Everything else here is scaffolding around that distinction.
 *
 * ⚠ jsdom sees no layout, so nothing here can assert that the tabs are legible or that the sticky
 * row stays put while the list scrolls — that needs a screenshot beside the page, which is the
 * standing rule for this repo's UI work.
 */
import { flushPromises, mount } from "@vue/test-utils";
import { createPinia, setActivePinia } from "pinia";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { success } from "@/api/result";
import { asCatalogueKey, type GenomeSourceGroup } from "@/api/types";
import AnchorGenomeBox from "@/components/anchor/AnchorGenomeBox.vue";
import { useAnchorGenomeStore } from "@/stores/anchorGenomeStore";

const api = vi.hoisted(() => ({ fetchGenomes: vi.fn() }));
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  ...api,
}));

const CATALOGUE = asCatalogueKey("ecoli-nuna4");

/** One row, in the shape the endpoint sends it. */
function row(sampleId: string, group: GenomeSourceGroup) {
  return {
    sample_id: sampleId,
    collection_genome_ordinal: 0,
    locus_count: 4866,
    arrangement_locus_count: 4800,
    source_group: group,
  };
}

/**
 * The listing endpoint, answering per group. ⚠ It echoes `source_group` back exactly as the server
 * does — the component relies on that echo to discard a late answer to a superseded tab.
 */
function servingGenomes(byGroup: Partial<Record<GenomeSourceGroup | "all", string[]>>) {
  api.fetchGenomes.mockImplementation(
    async (
      _catalogue: string,
      query: string,
      options: { sourceGroup?: GenomeSourceGroup | null } = {},
    ) => {
      const group = options.sourceGroup ?? null;
      const pool = (group === null ? byGroup.all : byGroup[group]) ?? [];
      const hits = pool.filter((id) => id.toUpperCase().includes(query.toUpperCase()));
      return success({
        species_key: "ecoli",
        query,
        genome_count: 100,
        matched_genome_count: hits.length,
        truncated: false,
        source_group: group,
        genomes: hits.map((id) => row(id, group ?? "atb")),
      });
    },
  );
}

async function openBox() {
  const box = mount(AnchorGenomeBox, {
    props: { catalogueKey: CATALOGUE, placement: "track" as const },
  });
  await box.find("input").trigger("focus");
  await flushPromises();
  return box;
}

/**
 * Each tab as `"<name> <count>"`.
 *
 * ⚠ Composed from the two spans rather than read off `tab.text()`: the gap between them is flex
 * `gap`, which is layout and not text, so `textContent` returns `"Complete (NCBI)0"` with no
 * separator. Asserting against that string would pin a concatenation artefact as if it were the
 * label.
 */
function tabLabels(box: ReturnType<typeof mount>) {
  return box
    .findAll(".anchor-group")
    .map((tab) => `${tab.find(".anchor-group-name").text()} ${tab.find(".anchor-group-n").text()}`);
}

beforeEach(() => {
  setActivePinia(createPinia());
  api.fetchGenomes.mockReset();
});

// ── the tabs ───────────────────────────────────────────────────────────────────────────────────
describe("the group tabs", () => {
  it("offer BOTH promised lists even when one is empty — this is every catalogue today", async () => {
    // ⛔⛔ The policy that changed after looking at the rendered page. Hiding the empty
    // complete-genome tab obeyed the published page's "do not offer a control that cannot answer"
    // rule — but the census removes that rule's precondition: the control now KNOWS the catalogue
    // holds none and can say so. Hidden, a reader cannot tell "none ingested yet" from "this site
    // does not do complete genomes", which is the distinction the census exists for. And on all
    // four catalogues loaded today the group is empty, so hidden means nobody ever learns it exists.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });
    servingGenomes({ atb: ["SAMEA103923484"], ncbi_complete: [] });

    const box = await openBox();

    expect(tabLabels(box)).toEqual(["Complete (NCBI) 0", "AllTheBacteria 100"]);
    // ⚠ And it opens on the non-empty one, so the first thing the reader sees is a usable list.
    expect(anchor.sourceGroup).toBe("atb");
    expect(box.text()).toContain("SAMEA103923484");
  });

  it("do NOT offer the draft group when it is empty — it is the leftover bucket, not a promise", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });
    servingGenomes({ atb: ["SAMEA103923484"] });

    const box = await openBox();

    expect(tabLabels(box).join(" ")).not.toContain("Draft");
  });

  it("appear once a second group has members, COMPLETE GENOMES FIRST", async () => {
    // ⛔ Not alphabetical and not census order: the complete genomes are what a reader anchors
    // deliberately, so they lead. `atb` leading would bury the short list under the long one.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });
    servingGenomes({ ncbi_complete: ["GCF_000240185.1"], atb: ["SAMEA103923484"] });

    const box = await openBox();

    expect(box.findAll(".anchor-group-name").map((n) => n.text())).toEqual([
      "Complete (NCBI)",
      "AllTheBacteria",
    ]);
  });

  it("show what the catalogue HOLDS, not what matched the query", async () => {
    // ⛔ A count that moved as the reader typed would stop being the thing that says whether a
    // group is worth opening at all.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });
    servingGenomes({ ncbi_complete: ["GCF_000240185.1"], atb: ["SAMEA103923484"] });

    const box = await openBox();
    const before = box.findAll(".anchor-group-n").map((span) => span.text());

    await box.find("input").setValue("ZZZZ-matches-nothing");
    await box.find("input").trigger("input");
    await flushPromises();

    expect(box.findAll(".anchor-group-n").map((span) => span.text())).toEqual(before);
    expect(before).toEqual(["4,017", "75,039"]);
  });

  it("omit a group with no members rather than offering an empty tab", async () => {
    // ⚠ The published page's own rule, applied one level down: a control that cannot tell "not
    // here" from "no idea" must not be offered. An empty tab makes the same face at both.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });
    servingGenomes({ ncbi_complete: ["GCF_000240185.1"], atb: ["SAMEA103923484"] });

    const box = await openBox();

    expect(tabLabels(box).join(" ")).not.toContain("Draft");
  });
});

// ── which list opens ───────────────────────────────────────────────────────────────────────────
describe("the opening list", () => {
  it("is the first NON-EMPTY group, never a hardcoded default", async () => {
    // ⛔ Opening on "Complete genomes" would show every reader of every catalogue loaded today an
    // empty list as the first thing the control does.
    const anchor = useAnchorGenomeStore();
    anchor.setSourceGroupCounts({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });

    expect(anchor.sourceGroup).toBe("atb");
  });

  it("is the complete-genome list once that group is populated", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });

    expect(anchor.sourceGroup).toBe("ncbi_complete");
  });

  it("is UNFILTERED when the census is absent — the pre-census behaviour, not a crash", async () => {
    // ⛔⛔ Code and data deploy by separate paths, so this bundle can be talking to a server that
    // predates `genome_source_groups`, or reading a shell cached with `max-age=86400` before it
    // existed. The first version indexed the census directly and took the whole page down.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts(undefined);
    servingGenomes({ all: ["SAMEA103923484"] });

    const box = await openBox();

    expect(anchor.sourceGroup).toBeNull();
    expect(tabLabels(box)).toEqual([]);
    // ⛔ And NO filter is sent: a guessed group would be answered with an empty list, and the
    // picker would look as though it had silently stopped listing the reader's genomes.
    expect(api.fetchGenomes.mock.calls[0]?.[2]).toMatchObject({ sourceGroup: null });
    expect(box.text()).toContain("SAMEA103923484");
  });

  it("survives a PARTIAL census", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setSourceGroupCounts({ atb: 100 } as Record<GenomeSourceGroup, number>);

    expect(anchor.sourceGroupCounts).toEqual({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });
    expect(anchor.sourceGroup).toBe("atb");
  });
});

// ── switching ──────────────────────────────────────────────────────────────────────────────────
describe("switching list", () => {
  beforeEach(() => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });
    servingGenomes({ ncbi_complete: ["GCF_000240185.1"], atb: ["SAMEA103923484"] });
  });

  it("re-queries the other group and KEEPS the typed text", async () => {
    // ⚠ Keeping the query is the point of the tabs: a reader who typed something and switches is
    // asking the same question of the other population.
    const box = await openBox();
    await box.find("input").setValue("SAM");
    await box.find("input").trigger("input");
    await flushPromises();

    await box.findAll(".anchor-group")[1]!.trigger("click");
    await flushPromises();

    const last = api.fetchGenomes.mock.calls.at(-1);
    expect(last?.[1]).toBe("SAM");
    expect(last?.[2]).toMatchObject({ sourceGroup: "atb" });
    expect((box.find("input").element as HTMLInputElement).value).toBe("SAM");
  });

  it("does NOT clear the anchor — browsing the lists is not destructive", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAnchor("GCF_000240185.1");
    const box = await openBox();

    await box.findAll(".anchor-group")[1]!.trigger("click");
    await flushPromises();

    expect(anchor.sampleId).toBe("GCF_000240185.1");
  });

  it("ignores a click on the list already shown", async () => {
    const box = await openBox();
    const before = api.fetchGenomes.mock.calls.length;

    await box.findAll(".anchor-group")[0]!.trigger("click");
    await flushPromises();

    expect(api.fetchGenomes.mock.calls.length).toBe(before);
  });

  it("Enter cannot anchor a genome from the list the reader has LEFT", async () => {
    // ⛔ The cursor indexes into the answer, so a stale answer plus a moved tab anchors a genome
    // from a list the reader is no longer looking at — and the accession looks plausible.
    const box = await openBox();
    await box.find("input").trigger("keydown", { key: "ArrowDown" });

    // Switch tab but never let the new answer arrive.
    api.fetchGenomes.mockImplementation(() => new Promise(() => {}));
    await box.findAll(".anchor-group")[1]!.trigger("click");
    await box.find("input").trigger("keydown", { key: "Enter" });

    const anchor = useAnchorGenomeStore();
    expect(anchor.sampleId).toBeNull();
  });
});

// ── ⭐ the two empty states, which is what this file is for ─────────────────────────────────────
describe("an empty list says WHY it is empty", () => {
  it('says the catalogue HOLDS NONE when the group is empty for this catalogue', async () => {
    // ⛔⛔ The sentence that stops a reader typing. True of the complete-genome list on all four
    // catalogues loaded today.
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });
    anchor.setSourceGroup("ncbi_complete");
    servingGenomes({ atb: ["SAMEA103923484"] });

    const box = await openBox();

    // ⚠ Written per group, not composed from the label: composing produced "This catalogue holds
    // no complete (ncbi)." — a sentence that reads like a bug. And it says *yet*, which is true:
    // the 4,017 complete genomes are not ingested, not absent in principle.
    expect(box.find(".hit-none").text()).toBe("No complete NCBI genomes are in this catalogue yet.");
  });

  it("says NOTHING MATCHED when the group has members but the query excludes them", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts({ ncbi_complete: 0, ncbi_draft: 0, atb: 100 });
    servingGenomes({ atb: ["SAMEA103923484"] });

    const box = await openBox();
    await box.find("input").setValue("ZZZZ");
    await box.find("input").trigger("input");
    await flushPromises();

    expect(box.find(".hit-none").text()).toBe(
      "No genome in AllTheBacteria matches “ZZZZ”.",
    );
  });

  it("falls back to the ORIGINAL sentence when there is no census to scope it", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setSourceGroupCounts(null);
    servingGenomes({ all: [] });

    const box = await openBox();
    await box.find("input").setValue("ZZZZ");
    await box.find("input").trigger("input");
    await flushPromises();

    expect(box.find(".hit-none").text()).toBe(
      "No genome in this catalogue matches “ZZZZ”.",
    );
  });
});

// ── the rows ───────────────────────────────────────────────────────────────────────────────────
it("a chosen genome still becomes the anchor, whichever list it came from", async () => {
  const anchor = useAnchorGenomeStore();
  anchor.setAvailability(true);
  anchor.setSourceGroupCounts({ ncbi_complete: 4017, ncbi_draft: 0, atb: 75039 });
  servingGenomes({ ncbi_complete: ["GCF_000240185.1"], atb: ["SAMEA103923484"] });

  const box = await openBox();
  await box.findAll(".anchor-hit:not(.clear)")[0]!.trigger("click");

  expect(anchor.sampleId).toBe("GCF_000240185.1");
  expect(anchor.kind).toBe("catalogue");
});
