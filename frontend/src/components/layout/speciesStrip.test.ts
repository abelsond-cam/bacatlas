/**
 * @vitest-environment jsdom
 *
 * The cohort picker — the control the browser's axis change is FOR, and which had no test at all.
 *
 * ⛔ **`SpeciesStrip.vue` was mounted by nothing.** It is stubbed in `tests/pageShellBehaviour.
 * test.ts`, the only suite that mounts `App`, and imported nowhere else — so `modelName()`, the one
 * live reader of the picker's wire block, had zero coverage, and the three edits that turned it
 * into a cohort picker would have been unverifiable by all 735 frontend tests.
 *
 * ⚠ And this file lives under `src/`, deliberately. `frontend/tsconfig.json` includes only
 * `src/**`, while vitest runs `src/**` AND `tests/**` — so everything under `frontend/tests/` is
 * RUN but never type-checked (measured: `vue-tsc --listFilesOnly` compiles 418 files, 0 of them
 * under `tests/`). A fixture written there is pinned to no type at all, however it is annotated.
 * Colocated here, the `CatalogueEntry` literals below are checked by `npm run build`.
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it } from "vitest";

import { asCatalogueKey } from "@/api/types";
import type { CatalogueEntry } from "@/api/types";
import SpeciesStrip from "@/components/layout/SpeciesStrip.vue";
import type { SpeciesEntry } from "@/stores/speciesCatalogueStore";

/** ⚠ Typed, and therefore checked — the point of this file living under `src/` rather than
 * `tests/`. The first draft omitted four fields and `vue-tsc` said so. */
const SPECIES: readonly SpeciesEntry[] = [
  {
    key: "kp", scientific_name: "Klebsiella pneumoniae", ncbi_taxonomy_id: 573, published: true,
    genome_count: 100, gene_count: 532_851, locus_count: 15_670,
    catalogue_key: asCatalogueKey("kp-probe100"),
  },
  {
    key: "ecoli", scientific_name: "Escherichia coli", ncbi_taxonomy_id: 562, published: true,
    genome_count: 100, gene_count: 489_146, locus_count: 17_531,
    catalogue_key: asCatalogueKey("ecoli-probe100"),
  },
];

function catalogue(
  key: string,
  label: string,
  cohortGenomes: number,
  isDefault = false,
): CatalogueEntry {
  return {
    key: asCatalogueKey(key),
    species: { key: "kp", scientific_name: "Klebsiella pneumoniae" },
    cohort: { key: key.split("-")[1] ?? key, label, genome_count: cohortGenomes },
    model_label: "kp_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
    is_default: isDefault,
    genome_count: cohortGenomes,
    gene_count: 532_851,
    locus_count: 15_670,
    run_id: "kp_bacformer_clever_exploded_preclusterkp98pm3b-3b0.5-excl_k100_g100_res0.1_seed0",
  };
}

const KP_AMR = catalogue("kp-probe100", "Kp AMR", 100, true);
const KP_FULL = catalogue("kp-71k", "Kp", 71_435);

function strip(catalogues: CatalogueEntry[]) {
  return mount(SpeciesStrip, {
    props: {
      species: SPECIES,
      speciesKey: "kp",
      catalogueKey: catalogues[0]?.key ?? null,
      catalogues,
      catalogue: null,
    },
    global: { stubs: { PangenomeCensus: true } },
  });
}

describe("the option names the COHORT and states its size", () => {
  it("composes the label and the count, and never reads a count out of the label", () => {
    const text = strip([KP_AMR, KP_FULL]).text();
    expect(text).toContain("Kp AMR — 100 genomes");
    // ⛔ 71,435 with a thousands separator, from `cohort.genome_count` — NOT the "71,500" David
    // wrote when naming the cohorts. The label carries the words and the roster carries the
    // number, so the two can never disagree; stored verbatim his string would have rendered
    // "Kp 71,500 — 71,435 genomes".
    expect(text).toContain("Kp — 71,435 genomes");
    expect(text).not.toContain("71,500");
  });

  it("marks the default, so a reader can tell which a bare species link resolves to", () => {
    expect(strip([KP_AMR, KP_FULL]).text()).toContain("Kp AMR — 100 genomes (default)");
  });

  it("⛔ says Cohort, not Model — the legend IS the claim the axis change makes", () => {
    const text = strip([KP_AMR, KP_FULL]).text();
    expect(text).toContain("Cohort:");
    expect(text).not.toContain("Model:");
  });
});

describe("⛔ the cohort is shown even when there is only one, reversing the old rule", () => {
  it("names the single cohort as static text rather than hiding the control", () => {
    // The control used to be hidden below two catalogues, on the sound argument that a choice of
    // one is not a choice. But it is now the ONLY place the cohort is named, and after the
    // re-ingest E. coli has exactly one — so hidden, that reader would never learn which genomes
    // they are looking at.
    const wrapper = strip([KP_AMR]);
    expect(wrapper.text()).toContain("Kp AMR — 100 genomes");
    expect(wrapper.find("select#cohort").exists()).toBe(false);
    expect(wrapper.find(".cohort-sole").exists()).toBe(true);
  });

  it("offers a select, not static text, as soon as there is a second", () => {
    const wrapper = strip([KP_AMR, KP_FULL]);
    expect(wrapper.find("select#cohort").exists()).toBe(true);
    expect(wrapper.find(".cohort-sole").exists()).toBe(false);
    expect(wrapper.findAll("select#cohort option")).toHaveLength(2);
  });

  it("shows neither when the species has no catalogue at all", () => {
    const wrapper = strip([]);
    expect(wrapper.find("select#cohort").exists()).toBe(false);
    expect(wrapper.find(".cohort-sole").exists()).toBe(false);
  });
});

describe("choosing a cohort is a navigation", () => {
  it("emits the catalogue key of the cohort picked, not the species", async () => {
    const wrapper = strip([KP_AMR, KP_FULL]);
    await wrapper.find("select#cohort").setValue("kp-71k");
    expect(wrapper.emitted("selectCatalogue")).toEqual([["kp-71k"]]);
  });

  it("⛔ the species control emits that species' DEFAULT catalogue key, never its species key", async () => {
    // A bare species key in `?catalogue=` resolves to the default with a 200 and is never re-pinned,
    // so the address would then name something the picker cannot match to an option.
    const wrapper = strip([KP_AMR, KP_FULL]);
    await wrapper.find("select#species").setValue("ecoli");
    expect(wrapper.emitted("selectCatalogue")).toEqual([["ecoli-probe100"]]);
  });
});
