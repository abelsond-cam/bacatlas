/**
 * @vitest-environment jsdom
 *
 * The page shell's own components — the parts assembling the page added, each pinned on the claim it
 * exists to get right rather than on its markup.
 */
import { readFileSync, readdirSync } from "node:fs";
import { resolve } from "node:path";

import { flushPromises, mount } from "@vue/test-utils";
import { createPinia, setActivePinia } from "pinia";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { failure, success } from "@/api/result";
import { asCatalogueKey } from "@/api/types";
import type { SpeciesCatalogueResponse } from "@/api/types";
import AnchorGenomeBox from "@/components/anchor/AnchorGenomeBox.vue";
import PangenomeCensus from "@/components/census/PangenomeCensus.vue";
import AuditResidualLists from "@/components/footer/AuditResidualLists.vue";
import LocusTrail from "@/components/navigation/LocusTrail.vue";
import LocusSearchBox from "@/components/search/LocusSearchBox.vue";
import OwnGenomesBox from "@/components/anchor/OwnGenomesBox.vue";
import OwnGenomesDialog from "@/components/ownGenomes/OwnGenomesDialog.vue";
import TrackPanel from "@/components/layout/TrackPanel.vue";
import SiteHeader from "@/components/layout/SiteHeader.vue";
import ViewTabs from "@/components/layout/ViewTabs.vue";
import { useAnchorGenomeStore } from "@/stores/anchorGenomeStore";
import { useOwnGenomesStore } from "@/stores/ownGenomesStore";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { VIEW_TABS } from "@/stores/viewTabStore";

const searchLoci = vi.hoisted(() => vi.fn());
const fetchGenomes = vi.hoisted(() => vi.fn());
const fetchAuditResiduals = vi.hoisted(() => vi.fn());
const fetchLocus = vi.hoisted(() => vi.fn());
// ⚠ A PARTIAL mock: `anchorQuery` must stay real, because it is the one place the anchor's kind
// becomes a query parameter and a stub would hide a projected genome being fetched as `anchor=`.
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  searchLoci,
  fetchGenomes,
  fetchAuditResiduals,
  fetchLocus,
}));

beforeEach(() => {
  setActivePinia(createPinia());
  searchLoci.mockReset();
  fetchGenomes.mockReset();
  fetchAuditResiduals.mockReset();
});

function catalogue(): SpeciesCatalogueResponse {
  return {
    pangenome: { genome_count: 100, gene_count: 489_146, locus_count: 17_531 },
    prevalence_census: { core: 3_117, soft_core: 268, shell: 2_964, cloud: 5_724, rare: 5_458 },
    prevalence_gene_census: { core: 310_000, soft_core: 27_500, shell: 60_000, cloud: 86_146, rare: 5_500 },
  } as unknown as SpeciesCatalogueResponse;
}

describe("the census — two partitions of one catalogue", () => {
  it("⭐ its per-genome parts SUM to its per-genome whole, because every gene is in one band", () => {
    const text = mount(PangenomeCensus, { props: { catalogue: catalogue() } }).text();
    // 489,146 / 100 = 4,891 genes; core (310,000 + 27,500) / 100 = 3,375; accessory 1,461; singletons 55
    expect(text).toContain("4,891 genes");
    expect(text).toContain("3,375 core");
    expect(text).toContain("1,461 accessory");
    expect(text).toContain("55 singletons");
    expect(3_375 + 1_461 + 55).toBe(4_891);
  });

  it("⭐ is TWO lines: genomes and genes → the loci they were modelled into, then per genome", () => {
    // David, 2026-09-22. The spaces between the parts are flex gaps, not text, so the line is read
    // as pieces IN ORDER rather than as one string.
    const lines = mount(PangenomeCensus, { props: { catalogue: catalogue() } }).findAll(".pg-line");
    expect(lines).toHaveLength(2);
    const first = lines[0]!.text();
    const pieces = [
      "100 genomes",
      "489,146 genes modelled",
      "→",
      "17,531 loci,",
      "including:",
      "3,385 core,",
      "8,688 accessory,",
      "5,458 singletons",
    ];
    const positions = pieces.map((piece) => first.indexOf(piece));
    expect(positions.every((position) => position >= 0)).toBe(true);
    expect([...positions].sort((a, b) => a - b)).toEqual(positions);
    // ⚠ No comma after the LAST band: "5,458 singletons," would read as a list cut short.
    expect(first.endsWith("5,458 singletons")).toBe(true);
    expect(lines[1]!.text()).toMatch(/^per genome:/);
  });

  it("⚠ calls the rare band `singletons` on BOTH lines — one band, one name", () => {
    const lines = mount(PangenomeCensus, { props: { catalogue: catalogue() } }).findAll(".pg-line");
    expect(lines[0]!.text()).toContain("5,458 singletons");
    expect(lines[1]!.text()).toContain("singletons");
    expect(lines.map((line) => line.text()).join(" ")).not.toContain("rare");
  });
});

describe("the view tabs", () => {
  it("⭐ put the three views OF THIS LOCUS first and the site's documentation last", () => {
    // David, 2026-09-22 — "Navigating BacAtlas" had sat second, between the evidence and the sequence.
    const labels = mount(ViewTabs, { props: { view: "locus" } }).findAll(".view-tab").map((tab) => tab.text());
    expect(labels).toEqual(["Syntelog Loci", "Sequence", "Function", "Navigating BacAtlas"]);
  });
});

describe("the breadcrumb", () => {
  it("draws nothing for a single step, and names every crumb it can", () => {
    expect(mount(LocusTrail, { props: { trail: ["1"], displayNames: new Map() } }).text()).toBe("");
    const trail = mount(LocusTrail, {
      props: { trail: ["1", "2"], displayNames: new Map([["1", "wzi"]]) },
    });
    expect(trail.find("button").text()).toBe("wzi");
    // The current locus is not a button: going to where you are is not a step.
    expect(trail.findAll("button")).toHaveLength(1);
    expect(trail.text()).toContain("·2");
  });
});

describe("⛔ the search says which of three things happened", () => {
  async function typed(query: string) {
    vi.useFakeTimers();
    const box = mount(LocusSearchBox, { props: { catalogueKey: asCatalogueKey("ecoli"), collectionGenomeCount: 100 } });
    await box.find("input").setValue(query);
    await vi.advanceTimersByTimeAsync(200);
    vi.useRealTimers();
    await flushPromises();
    return box;
  }

  it("a FAILED search reads as a failure, never as 'no locus matches'", async () => {
    searchLoci.mockResolvedValueOnce(failure("server", "the server answered 500", 500));
    const box = await typed("ligase");
    expect(box.find(".pop-error").text()).toContain("did not complete");
    expect(box.text()).not.toContain("No locus matches");
  });

  it("an EMPTY answer says so, with the naming caveat", async () => {
    searchLoci.mockResolvedValueOnce(success({ query: "waaL", mode: "substring", truncated: false, hits: [] }));
    const box = await typed("waaL");
    expect(box.text()).toContain("No locus matches “waaL”");
    expect(box.text()).toContain("filed as rfaL, not waaL");
  });

  it("⚠ a PREFIX match says it searched less, and Enter goes to the cursor's hit", async () => {
    searchLoci.mockResolvedValueOnce(
      success({
        query: "rf",
        mode: "prefix",
        truncated: false,
        hits: [
          { label: "48", display_name: "rfaL", best_product: "O-antigen ligase", gene_count: 100, genome_count: 100, prevalence_band: "core", rank_band: 1 },
          { label: "49", display_name: "rfaJ", best_product: null, gene_count: 90, genome_count: 90, prevalence_band: "soft_core", rank_band: 1 },
        ],
      }),
    );
    const box = await typed("rf");
    expect(box.text()).toContain("match the start of a name only");
    expect(box.findAll(".hit-desc").map((node) => node.text())).toEqual(["O-antigen ligase", "—"]);
    await box.find("input").trigger("keydown", { key: "ArrowDown" });
    await box.find("input").trigger("keydown", { key: "Enter" });
    expect(box.emitted("go")).toEqual([["49"]]);
  });
});

describe("⛔ the anchor box", () => {
  it("is REMOVED, not disabled, when the catalogue cannot answer an anchor", () => {
    const box = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    expect(box.find("input").exists()).toBe(false);
  });

  it("⭐ offers the OFF switch first, and choosing a genome sets the ONE anchor", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    fetchGenomes.mockResolvedValue(
      success({
        species_key: "ecoli", query: "", genome_count: 100, matched_genome_count: 100, truncated: true,
        genomes: [{ sample_id: "SAMEA1", collection_genome_ordinal: 0, locus_count: 4_321, arrangement_locus_count: 4_300 }],
      }),
    );
    const track = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    const sequence = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "sequence" } });
    await track.find("input").trigger("focus");
    await flushPromises();
    const rows = track.findAll(".anchor-hit");
    expect(rows[0]!.text()).toBe("most common arrangement");
    expect(rows[1]!.text()).toContain("4,321 loci");
    // ⚠ A cut list says how much it left out.
    expect(track.text()).toContain("1 of 100 matching genomes");
    await rows[1]!.trigger("click");
    expect(anchor.sampleId).toBe("SAMEA1");
    // ⭐ Two mountings of ONE state: the other box shows it without being told.
    expect((sequence.find("input").element as HTMLInputElement).value).toBe("SAMEA1");
  });

  it("⛔ a genome list that did not load SAYS so, rather than looking like no genomes", async () => {
    useAnchorGenomeStore().setAvailability(true);
    fetchGenomes.mockResolvedValue(failure("not_found", "the server answered 404", 404));
    const box = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    await box.find("input").trigger("focus");
    await flushPromises();
    expect(box.find(".pop-error").text()).toContain("did not load");
    expect(box.text()).not.toContain("No genome in this catalogue matches");
  });
});

describe("⛔ the residual lists — two, never one, and fetched only when opened", () => {
  const AUDIT = { synteny_only_n_clusters: 8, no_homology_n_clusters: 2, pfam_conflict_n_clusters: 42 };

  it("states each list's LENGTH from the audit headline before fetching anything", () => {
    const lists = mount(AuditResidualLists, { props: { catalogueKey: asCatalogueKey("ecoli"), audit: AUDIT } });
    const summaries = lists.findAll("summary").map((node) => node.text());
    expect(summaries).toEqual(["Every locus grouped on context alone (10)", "Every Pfam conflict (42)"]);
    expect(fetchAuditResiduals).not.toHaveBeenCalled();
  });

  it("fetches on open and lists each row with the evidence to judge it", async () => {
    fetchAuditResiduals.mockResolvedValueOnce(
      success({
        grouped_on_context_alone: [
          {
            label: "4976", display_name: "ydcD", prevalence_band: "shell", gene_count: 32, uniref50_family_count: 3,
            // ⛔ SIMILARITIES, not the distances this row used to carry: 0.98 is what the footer
            // prints, and a surviving `1 − d` in the component would print 0.02.
            pfam_architecture_count: null, syntenic_a5: 0.35, esm_within_similarity: 0.98, esm_nearest_similarity: 0.98,
          },
        ],
        pfam_conflicts: [],
      }),
    );
    const lists = mount(AuditResidualLists, { props: { catalogueKey: asCatalogueKey("ecoli"), audit: AUDIT } });
    const first = lists.find("details");
    (first.element as HTMLDetailsElement).open = true;
    await first.trigger("toggle");
    await flushPromises();
    const row = lists.find(".goto-row");
    expect(row.text()).toContain("ydcD");
    // ⚠ Pfam could not judge it: no "N architectures" at all, rather than "0 architectures".
    expect(row.find(".goto-ev").text()).toBe("shell · 32 genes · 3 UniRef50 · A5 0.35 · ESM 0.98/0.98");
    await row.trigger("click");
    expect(lists.emitted("go")).toEqual([["4976"]]);
  });
});

describe("⛔ a locus the address names but the catalogue does not have", () => {
  it("is named in the reader's terms, and the starting locus is OFFERED, not substituted", async () => {
    const navigation = useLocusNavigationStore();
    navigation.setCatalogue(asCatalogueKey("ecoli"));
    fetchLocus.mockResolvedValueOnce(failure("not_found", "no locus '999999' in pangenome 1. ⚠ Node labels are TEXT", 404));
    await navigation.navigateTo("999999");
    const panel = mount(TrackPanel, { props: { catalogueKey: asCatalogueKey("ecoli"), speciesKey: "ecoli", collectionGenomeCount: 100, landingLocus: "2811" } });
    const notice = panel.find(".track-error").text();
    expect(notice).toContain("There is no locus “999999” in this catalogue.");
    // The server's sentence is for whoever debugs a label, not for a reader following an old link.
    expect(notice).not.toContain("pangenome 1");
    fetchLocus.mockResolvedValueOnce(failure("network", "offline"));
    await panel.find(".track-error button").trigger("click");
    expect(fetchLocus).toHaveBeenLastCalledWith("ecoli", "2811", expect.anything());
  });
});

/** The two published species, as the dialog receives them. */
const SPECIES = [
  // ⚠ `catalogue_key` is the species' DEFAULT catalogue — what the dialog ADDRESSES. The answers
  // are still filed under `key`, the organism, which is what lets it say "modelled here".
  { key: "ecoli", scientific_name: "Escherichia coli", published: true, catalogue_key: "ecoli-nuna4" },
  { key: "kp", scientific_name: "Klebsiella pneumoniae", published: true, catalogue_key: "kp-nuna4" },
] as unknown as Parameters<typeof mountDialog>[0]["species"];

function genomeList(speciesKey: string, sampleIds: readonly string[]) {
  return success({
    species_key: speciesKey,
    query: "",
    genome_count: sampleIds.length,
    matched_genome_count: sampleIds.length,
    truncated: false,
    genomes: sampleIds.map((sample_id, index) => ({
      sample_id,
      collection_genome_ordinal: index,
      locus_count: 4_000,
      arrangement_locus_count: 3_900,
    })),
  });
}

function mountDialog(props: { speciesKey: string | null; species: readonly unknown[] }) {
  return mount(OwnGenomesDialog, { props: props as never });
}

/** Open the dialog at the add step and let the catalogue check run. */
async function openAddStep() {
  const own = useOwnGenomesStore();
  own.open("account");
  const dialog = mountDialog({ speciesKey: "ecoli", species: SPECIES });
  own.showStep("add");
  await flushPromises();
  return { own, dialog };
}

describe("⭐ View your own genomes — the second box in the gutter", () => {
  it("the ⚓ is drawn even when nothing is anchored: it names the control, it does not claim a genome", () => {
    // David, 2026-09-22. The published rule hid it, because "an anchor on an empty box would be a claim
    // about a genome nobody has chosen". Grey makes no claim; the box still colours only when anchored.
    useAnchorGenomeStore().setAvailability(true);
    const box = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    expect(box.find(".arr-anchor-mark").exists()).toBe(true);
    expect(box.find(".arr-anchor").classes()).not.toContain("on");
  });

  it("opens the dialog on the accounts step", async () => {
    const box = mount(OwnGenomesBox);
    await box.find(".own-genomes-box").trigger("click");
    const own = useOwnGenomesStore();
    expect(own.isOpen).toBe(true);
    expect(own.step).toBe("account");
  });

  it("⛔ says on the FIRST screen what is not built, and that nothing is saved", () => {
    useOwnGenomesStore().open("account");
    const dialog = mountDialog({ speciesKey: "ecoli", species: SPECIES });
    expect(dialog.text()).toContain("Accounts and login will be required — to be completed");
    expect(dialog.text()).toContain("forgotten when you reload");
    // The pipeline is described as what the service does, for collaborators to react to.
    expect(dialog.text()).toContain("BakRep");
    expect(dialog.text()).toContain("Bacformer");
    expect(dialog.find(".own-button[disabled]").text()).toContain("Sign in");
  });

  it("sorts a reader's file into: ours, the other species', not-an-accession, and a repeat", async () => {
    // ⚠ Called with the CATALOGUE key now — the roster is fetched from a catalogue, and filed
    // under the organism. Keying this mock on the bare species silently returned kp's list for
    // both, which is how the "modelled elsewhere" row came out as "to be completed".
    fetchGenomes.mockImplementation((catalogueKey: string) =>
      Promise.resolve(
        catalogueKey === "ecoli-nuna4"
          ? genomeList("ecoli", ["SAMEA103923484"])
          : genomeList("kp", ["SAMN03892119"]),
      ),
    );
    const { own, dialog } = await openAddStep();
    own.setFile("mine.txt", "# mine\nSAMEA103923484\nSAMN03892119\nSAMN05374479\nGCA_012642945.1\nsamn05374479\n");
    await flushPromises();

    const rows = dialog.findAll(".own-lines li");
    expect(rows.map((row) => row.classes()[0])).toEqual([
      "modelled-here",
      "modelled-elsewhere",
      "to-be-completed",
      "not-an-accession",
      "repeat",
    ]);
    expect(rows[1]!.text()).toContain("Klebsiella pneumoniae");
    // ⛔ BioSample only: an assembly accession is refused, never converted behind the reader's back.
    expect(rows[3]!.text()).toContain("not a BioSample");
    // ⚠ The reader's own line comes back in the case they typed it.
    expect(rows[4]!.text()).toContain("samn05374479");
    // ⚠ Only genomes that are actually PLACED are offered — this catalogue has none in this test,
    // and the button says what happens next rather than "Show 0 genomes", which would read as a
    // fault in the reader's file rather than as the part of the service still to be built.
    expect(dialog.text()).toContain("What happens next");
  });

  it("⛔ a modelled-genome list that FAILED to load is not an empty one", async () => {
    fetchGenomes.mockResolvedValue(failure("server", "the server answered 500", 500));
    const { own, dialog } = await openAddStep();
    own.setFile("mine.txt", "SAMN05374479\n");
    await flushPromises();
    expect(dialog.find(".pop-error").text()).toContain("did not load");
    // Nothing is claimed about the accession: it is not called new, and it is not called modelled.
    expect(dialog.find(".own-lines li").classes()).toContain("unchecked");
    expect(dialog.text()).toContain("not checked");
  });
});

describe("⭐ a genome PLACED on the model, anchored beside it", () => {
  it("⛔ does not appear behind the ⚓ — that box means 'one of the modelled genomes'", async () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setProjectedAnchor("SAMN05374479");
    const box = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    const own = mount(OwnGenomesBox, { props: { speciesKey: "ecoli" } });
    await flushPromises();

    // The anchor box is back to its empty state: a placed genome is not one of the modelled 100.
    expect(box.find(".arr-anchor").classes()).not.toContain("on");
    expect((box.find(".arr-anchor input").element as HTMLInputElement).value).not.toContain(
      "SAMN05374479",
    );
    // …and the OTHER box carries it, with the word that says what it is.
    expect(own.find(".own-genomes-box").classes()).toContain("on");
    expect(own.text()).toContain("SAMN05374479");
    expect(own.text()).toContain("placed");
  });

  it("⛔ releases the other anchor, because the track draws one genome", () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    anchor.setAnchor("SAMEA103923484");
    expect(anchor.kind).toBe("catalogue");
    anchor.setProjectedAnchor("SAMN05374479");
    expect(anchor.sampleId).toBe("SAMN05374479");
    expect(anchor.kind).toBe("projected");
    anchor.setAnchor("SAMEA103923484");
    expect(anchor.kind).toBe("catalogue");
  });

  it("⚠ an unavailable catalogue anchor does not drop a PLACED one — different table, different question", () => {
    const anchor = useAnchorGenomeStore();
    anchor.setProjectedAnchor("SAMN05374479");
    anchor.setAvailability(false);
    expect(anchor.sampleId).toBe("SAMN05374479");
  });
});

describe("⭐ one header over both boxes: “Anchor to:”", () => {
  it("names the question once, and the boxes are its two answers (David, 2026-09-23)", () => {
    const anchor = useAnchorGenomeStore();
    anchor.setAvailability(true);
    const box = mount(AnchorGenomeBox, { props: { catalogueKey: asCatalogueKey("ecoli"), placement: "track" } });
    const own = mount(OwnGenomesBox, { props: { speciesKey: "ecoli" } });

    // ⚠ The two labels are the two ANSWERS, so neither repeats the question. A box reading
    // "Anchor to a genome" under a header reading "Anchor to:" says it twice.
    expect(box.find("input").attributes("placeholder")).toBe("A genome in BacAtlas");
    expect(own.text()).toContain("Add my own genome");
    expect(box.find("input").attributes("placeholder")).not.toContain("Anchor");
    expect(own.text()).not.toContain("Anchor");
  });

  it("⛔ the two boxes are the SAME WIDTH — they are a pair, and one narrower reads as subordinate", () => {
    // `.arr-anchor` carries `align-self: center` from when it was itself the grid item; inside the
    // wrapper that made it shrink to its intrinsic width (203px against its sibling's 260px) and
    // clipped the placeholder. jsdom computes no layout, so this pins the RULE that fixes it.
    // ⚠ From the project root, not `import.meta.url`: this suite runs under jsdom, where
    // `import.meta.url` is an http URL and resolving against it fails outright.
    const css = readFileSync(resolve(process.cwd(), "src/styles/app.css"), "utf8");
    expect(css).toMatch(/\.arr-gutter \.arr-anchor \{\s*align-self: stretch;/);
  });
});


describe("⛔ the site's own name, in every place a reader meets it", () => {
  /**
   * nuna's `test_the_page_carries_the_site_name_in_both_places` has guarded the static page since the
   * site went bacformer.org → BacMapper.org → back → Syntitude.org → BacAtlas.org. The app had no
   * equivalent, and a rename that reaches the header but not a tab label, a `<title>` or the ship's
   * `aria-label` shows up only somewhere nobody looks.
   */
  it("names the site in the wordmark, with the masthead beside it on ONE line", () => {
    const header = mount(SiteHeader, { props: { catalogueKey: asCatalogueKey("ecoli"), collectionGenomeCount: 100 } });
    // The wordmark splits over three spans, so read the text rather than the markup.
    expect(header.find(".mark").text()).toBe("BacAtlas.org — Navigate your Microbe");
  });

  it("names it in the documentation tab", () => {
    expect(VIEW_TABS.map((tab) => tab.label)).toContain("Navigating BacAtlas");
  });

  /**
   * ⭐ The real guard: the old name survives NOWHERE — not in a component, not in the ship's
   * `aria-label`, not in the two published pages at the repo root. Those pages are generated by nuna's
   * `render_page`, and on 2026-09-24 they were patched in place rather than re-rendered (the local
   * payload is schema 14 and the live pages are schema 9, so a re-render would publish an unverified
   * Sequence tab). This is what keeps that patch honest until the next real render.
   *
   * ⚠ It checks the capitalised NAME only — deliberately. The environment variables, the database,
   * the Python package and the repo name were renamed separately on 2026-09-26; matching them here
   * would fail for a reason that has nothing to do with what a reader sees.
   * ⛔ The needle below IS the old name. Renaming it makes this guard match every file, and the
   * obvious fix for that — broadening the exclusion — turns the guard into a permanent no-op.
   */
  it("⭐ the OLD name survives nowhere — not in the app, not in the two published pages", () => {
    const root = resolve(process.cwd(), "..");
    const files: string[] = [resolve(process.cwd(), "index.html")];
    const walk = (dir: string): void => {
      for (const entry of readdirSync(dir, { withFileTypes: true })) {
        const full = resolve(dir, entry.name);
        if (entry.isDirectory()) walk(full);
        else files.push(full);
      }
    };
    walk(resolve(process.cwd(), "src"));
    for (const page of ["ecoli.html", "kp.html", "index.html"]) files.push(resolve(root, page));

    // ⛔ Coverage before the verdict: a walk that found nothing would report "no occurrences" and look
    // exactly like a clean tree. Name the count, and name the pages by hand — they are the whole point.
    expect(files.length).toBeGreaterThan(100);
    for (const page of ["ecoli.html", "kp.html"]) {
      expect(readFileSync(resolve(root, page), "utf8")).toContain("BacAtlas<span>.org</span>");
    }

    // ⚠ This file is excluded BY NAME, not by a clever needle: its own docblock recounts the rename,
    // and a test that cannot say what it is guarding is worse than one that names its exception.
    const self = "frontend/src/components/layout/shell.test.ts";
    const offenders = files
      .filter((file) => readFileSync(file, "utf8").includes("Syntitude"))
      .map((file) => file.slice(root.length + 1))
      .filter((file) => file !== self);
    expect(offenders).toEqual([]);
  });
});
