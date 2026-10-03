/**
 * @vitest-environment jsdom
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it } from "vitest";

import type { AnnotationEntry, FunctionResponse, GoVerdict } from "@/api/types";

import CogCard from "./CogCard.vue";
import EnzymeAndKeggCard from "./EnzymeAndKeggCard.vue";
import FunctionTab from "./FunctionTab.vue";
import GeneOntologyCard from "./GeneOntologyCard.vue";

function entry(term: string, geneCount: number, overrides: Partial<AnnotationEntry> = {}): AnnotationEntry {
  return { rank: 0, term, name: `${term} name`, gene_count: geneCount, ...overrides };
}

function coverage(overrides: Partial<FunctionResponse["coverage"]> = {}): FunctionResponse["coverage"] {
  return {
    gene_count: 100,
    cog_annotated_gene_count: 60,
    cog_distinct_id_count: 1,
    modal_cog_categories: ["N"],
    ec_annotated_gene_count: 0,
    kegg_annotated_gene_count: 0,
    go_annotated_gene_count: {
      molecular_function: 40,
      biological_process: 40,
      cellular_component: 0,
    },
    ...overrides,
  };
}

function block(overrides: Partial<FunctionResponse> = {}): FunctionResponse {
  return {
    annotations: {
      cog_orthogroup: [entry("COG1132", 60)],
      gene_ontology_slim: [
        entry("GO:0016874", 40, { gene_ontology_namespace: "molecular_function" }),
        entry("GO:0006412", 40, { gene_ontology_namespace: "biological_process" }),
      ],
    },
    coverage: coverage(),
    go_verdicts: {
      molecular_function: "single",
      biological_process: "nested",
      cellular_component: "no_coverage",
    },
    ...overrides,
  };
}

function mountTab(props: Partial<InstanceType<typeof FunctionTab>["$props"]> = {}) {
  return mount(FunctionTab, {
    props: {
      displayName: "lysS",
      locusLabel: "17373",
      block: block(),
      status: "ready" as const,
      ...props,
    },
  });
}

// ── COG ────────────────────────────────────────────────────────────────────────────────────────
describe("⭐ partial COG coverage is HEADROOM, not a caveat", () => {
  it("says the unlabelled genes are unlabelled rather than different", () => {
    // David, 2026-08-20. COG is assigned by best hit, not by a profile, so a gene without one is
    // much weaker evidence of absence than a Pfam miss — which really is a profile failing to match.
    const card = mount(CogCard, { props: { coverage: coverage(), entries: [entry("COG1132", 60)] } });
    expect(card.text()).toContain("60 of 100 genes carry a COG assignment");
    expect(card.text()).toContain(
      "the remaining 40 are unlabelled rather than different — the gap this locus could fill, not a gap in it",
    );
  });

  it("⛔ does NOT predict a label for them", () => {
    // Propagating a COG from the annotated members is annotation transfer: it needs a transfer
    // rule, a confidence measure, a decision about discordant loci and a marking that can never be
    // mistaken for a Bakta call. None of that is shipped, so neither is a predicted label.
    const card = mount(CogCard, { props: { coverage: coverage(), entries: [entry("COG1132", 60)] } });
    const rows = card.findAll("tbody tr");
    expect(rows).toHaveLength(1);
    expect(card.text()).not.toContain("predicted");
  });

  it("omits the headroom sentence when every gene is annotated", () => {
    const card = mount(CogCard, {
      props: { coverage: coverage({ cog_annotated_gene_count: 100 }), entries: [entry("COG1132", 100)] },
    });
    expect(card.text()).not.toContain("unlabelled rather than different");
  });

  it("⛔ says what NO coverage means, and shows no chips at all", () => {
    const card = mount(CogCard, {
      props: {
        coverage: coverage({ cog_annotated_gene_count: 0, cog_distinct_id_count: 0, modal_cog_categories: null }),
        entries: [],
      },
    });
    expect(card.text()).toContain("no coverage here, so it can neither support nor contradict");
    expect(card.findAll(".chip")).toHaveLength(0);
    expect(card.text()).not.toContain("unlabelled rather than different");
  });
});

describe("⛔ the distinct-group count is a COUNT, never a relation", () => {
  it("calls two groups ordinary rather than a fault", () => {
    // More than one orthologous group is an ordinary consequence of grouping above the family
    // level, which is what this method does.
    const card = mount(CogCard, {
      props: { coverage: coverage({ cog_distinct_id_count: 2 }), entries: [entry("COG1132", 40), entry("COG0642", 20)] },
    });
    expect(card.text()).toContain("2 orthologous groups");
    expect(card.text()).toContain("an ordinary consequence of grouping above the family level");
    expect(card.find(".chip.neutral").exists()).toBe(true);
  });

  it("marks a single group as the quiet win, with no explanatory note", () => {
    const card = mount(CogCard, { props: { coverage: coverage(), entries: [entry("COG1132", 60)] } });
    expect(card.text()).toContain("one orthologous group");
    expect(card.text()).not.toContain("ordinary consequence");
    expect(card.find(".chip.win").exists()).toBe(true);
  });
});

describe("⚠ a COG category is one letter OR SEVERAL", () => {
  it("names each letter of a multi-letter category", () => {
    const card = mount(CogCard, {
      props: { coverage: coverage({ modal_cog_categories: ["EP"] }), entries: [entry("COG1132", 60)] },
    });
    expect(card.text()).toContain(
      "EP — Amino acid transport and metabolism · Inorganic ion transport and metabolism",
    );
  });

  it("shows no category chip where the locus has none", () => {
    const card = mount(CogCard, {
      props: { coverage: coverage({ modal_cog_categories: null }), entries: [entry("COG1132", 60)] },
    });
    expect(card.findAll(".chip")).toHaveLength(1);
  });
});

// ── GO ─────────────────────────────────────────────────────────────────────────────────────────
describe("⛔ GO: coverage before the verdict, and `no_coverage` gets no chip", () => {
  function mountGo(verdict: GoVerdict | null, annotated = 40) {
    return mount(GeneOntologyCard, {
      props: {
        namespace: "molecular_function" as const,
        annotatedGeneCount: annotated,
        geneCount: 100,
        verdict,
        entries: [entry("GO:0016874", 40, { name: "ligase activity" })],
      },
    });
  }

  it("shows a chip for every verdict on the ladder", () => {
    expect(mountGo("single").find(".chip").text()).toBe("one class");
    expect(mountGo("nested").find(".chip").text()).toBe("partial annotation");
    expect(mountGo("disjoint").find(".chip").text()).toBe("classes differ");
  });

  it("⛔ shows NO chip for no_coverage — it is the absence of a verdict, not one", () => {
    const card = mountGo("no_coverage", 0);
    expect(card.find(".chip").exists()).toBe(false);
    expect(card.text()).toContain("no coverage here, so it can neither support nor contradict");
  });

  it("names the namespace in both the heading and the coverage line", () => {
    const card = mountGo("single");
    expect(card.find(".sub-head").text()).toBe("GO — molecular function");
    expect(card.text()).toContain("40 of 100 genes carry a molecular function term");
  });

  it("⚠ shows the class NAME and the accession, because neither alone is usable", () => {
    // `GO:0016020` alone is unreadable; "membrane" alone is unlookupable.
    const card = mountGo("single");
    expect(card.find("tbody").text()).toContain("ligase activity");
    expect(card.find("tbody .acc-link").text()).toBe("GO:0016874");
  });
});

describe("⛔ all three namespaces, or one line — never three empty cards", () => {
  it("renders three cards when the locus has GO anywhere", () => {
    const tab = mountTab();
    const headings = tab.findAll(".sub-head").map((node) => node.text());
    expect(headings).toContain("GO — molecular function");
    expect(headings).toContain("GO — biological process");
    // ⭐ cellular_component has NO coverage and is still stated: `no coverage` and `the members
    // disagree` are different findings, and only the second is evidence.
    expect(headings).toContain("GO — cellular component");
  });

  it("collapses to ONE line when the locus has no GO at all", () => {
    const tab = mountTab({
      block: block({
        annotations: { cog_orthogroup: [entry("COG1132", 60)] },
        coverage: coverage({
          go_annotated_gene_count: {
            molecular_function: 0,
            biological_process: 0,
            cellular_component: 0,
          },
        }),
      }),
    });
    const headings = tab.findAll(".sub-head").map((node) => node.text());
    expect(headings.filter((heading) => heading.startsWith("GO"))).toEqual(["GO"]);
    expect(tab.text()).toContain("None of these 100 genes carries a GO term");
  });

  it("⛔ files each term under ITS OWN namespace, not by position", () => {
    // They arrive in one list carrying their namespace. Reading them positionally would file a
    // cellular-component class under molecular function, with three cards that all look right.
    const tab = mountTab();
    const cards = tab.findAllComponents(GeneOntologyCard);
    const molecular = cards.find((card) => card.props("namespace") === "molecular_function")!;
    const biological = cards.find((card) => card.props("namespace") === "biological_process")!;
    expect(molecular.props("entries").map((e: AnnotationEntry) => e.term)).toEqual(["GO:0016874"]);
    expect(biological.props("entries").map((e: AnnotationEntry) => e.term)).toEqual(["GO:0006412"]);
  });

  it("⚠ gives each namespace its OWN verdict, never the first one three times", () => {
    const cards = mountTab().findAllComponents(GeneOntologyCard);
    expect(cards.map((card) => card.props("verdict"))).toEqual(["single", "nested", "no_coverage"]);
  });
});

// ── EC and KEGG ────────────────────────────────────────────────────────────────────────────────
describe("⛔ a KEGG KO is linked and never named", () => {
  it("shows the accession, links it, and prints no description", () => {
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30, kegg_annotated_gene_count: 20 }),
        enzymeEntries: [entry("6.1.1.18", 30, { name: null })],
        keggEntries: [entry("K01886", 20, { name: null })],
      },
    });
    const kegg = card.findAll(".kv").at(-1)!;
    expect(kegg.find(".acc-link").text()).toBe("K01886");
    expect(kegg.find(".acc-link").attributes("href")).toBe("https://www.genome.jp/entry/K01886");
    expect(kegg.text()).not.toContain("name");
  });

  it("⚠ SEPARATES the count from the accession", () => {
    // Run together it reads as one token: `K15540` beside a count of 77 came back from a reader as
    // the question "what is K1554077?"
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ kegg_annotated_gene_count: 77 }),
        enzymeEntries: [],
        keggEntries: [entry("K15540", 77, { name: null })],
      },
    });
    expect(card.find(".alt-n").text()).toBe("×77");
    expect(card.find(".v").text()).not.toContain("K1554077");
  });

  it("states coverage per vocabulary, against the locus", () => {
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30, kegg_annotated_gene_count: 20 }),
        enzymeEntries: [entry("6.1.1.18", 30, { name: null })],
        keggEntries: [entry("K01886", 20, { name: null })],
      },
    });
    expect(card.findAll(".cov").map((node) => node.text())).toEqual([
      "30 of 100 annotated",
      "20 of 100 annotated",
    ]);
  });

  it("⚠ is ABSENT entirely when neither vocabulary has anything", () => {
    // Two headings over nothing say less than no card at all.
    const card = mount(EnzymeAndKeggCard, {
      props: { coverage: coverage(), enzymeEntries: [], keggEntries: [] },
    });
    expect(card.find(".card").exists()).toBe(false);
  });

  it("shows only the row that HAS something", () => {
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30 }),
        enzymeEntries: [entry("6.1.1.18", 30, { name: null })],
        keggEntries: [],
      },
    });
    expect(card.findAll(".kv")).toHaveLength(1);
    expect(card.find(".k").text()).toBe("EC");
  });
});

// ── the tab itself ─────────────────────────────────────────────────────────────────────────────
describe("⛔ three loading states, and a failure never reads as 'nothing here'", () => {
  it("says it is loading rather than showing an empty panel", () => {
    const tab = mountTab({ block: null, status: "pending" });
    expect(tab.text()).toContain("Loading the function annotation…");
    expect(tab.findAll(".card")).toHaveLength(0);
  });

  it("⛔ renders a FAILURE as a failure, with the server's own sentence and a retry", async () => {
    // `app.js:4604` — a panel that fails silently is one a reader will read as "this locus has no
    // function annotation", which on this tab is itself a finding, and a false one.
    const tab = mountTab({ block: null, status: "failed", failureDetail: "the server said 503" });
    expect(tab.find(".pop-error").text()).toContain("the server said 503");
    expect(tab.text()).not.toContain("Loading the function annotation");
    await tab.find(".pop-retry").trigger("click");
    expect(tab.emitted("retry")).toHaveLength(1);
  });

  it("heads the panel with the locus it is about", () => {
    const head = mountTab().find(".func-head");
    expect(head.find("h2").text()).toBe("Function");
    expect(head.find(".lid").text()).toBe("lysS ·17373");
  });

  it("⛔ draws the EC row — 4,544 loci carried an EC number and NOT ONE of them showed it", () => {
    // The tab asked `annotations` for `enzyme_commission`; the backend's `AnnotationKind` enum has
    // always said `ec_number`. `annotations` was `Record<string, …>`, so the wrong key type-checked,
    // rendered, and silently dropped the whole row — and this fixture had no EC entry, so 668 green
    // tests never touched it. gumC (kp-nuna4 locus 722) carries EC on 16 of its 100 genes and showed
    // no EC/KEGG card at all. The compile-time guard is `ANNOTATION_KINDS`; this is the behaviour.
    const tab = mountTab({
      block: block({
        annotations: {
          cog_orthogroup: [entry("COG1132", 60)],
          ec_number: [entry("2.7.10.-", 16)],
        },
        coverage: coverage({ ec_annotated_gene_count: 16 }),
      }),
    });
    expect(tab.text()).toContain("2.7.10.-");
    expect(tab.text()).toContain("16 of 100 annotated");
  });

  it("⚠ says these are not an independent source", () => {
    // They are the same annotations the gene names come from, so a reader taking the COG agreement
    // as corroboration of the merge is double-counting one source.
    expect(mountTab().text()).toContain("not an independent source");
    expect(mountTab().text()).toContain("folded onto the metagenomics GO slim");
  });
});
