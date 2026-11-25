/**
 * @vitest-environment jsdom
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it } from "vitest";

import type {
  AnnotationEntry,
  FunctionResponse,
  OwnSupport,
  PredictionsBlock,
  PropagationRate,
  ServedPrediction,
  VocabularyInference,
} from "@/api/types";
import { INFERENCE_KINDS } from "@/api/types";

import CogCard from "./CogCard.vue";
import EnzymeAndKeggCard from "./EnzymeAndKeggCard.vue";
import FunctionTab from "./FunctionTab.vue";
import GeneOntologyCard from "./GeneOntologyCard.vue";
import PredictedFunctionCard from "./PredictedFunctionCard.vue";

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

/**
 * ⛔ The own-call block and the predictions are REQUIRED on a `FunctionResponse`, so a fixture must carry
 * them. Empty by default — no own call, and predictions NOT LOADED — because a fixture that silently
 * shipped a prediction would make every tab under test look like a predicted one.
 */
function emptyInference(): Pick<FunctionResponse, "inference" | "predictions"> {
  return {
    inference: {
      vocabularies: INFERENCE_KINDS.map((annotation_kind) => ({ annotation_kind, own: null })),
    },
    predictions: noPredictions(),
  };
}

function noPredictions(state: "not_loaded" | "none_produced" = "not_loaded"): PredictionsBlock {
  return { state, annotation_model_version: null, annotation_model_label: null, serve_cut: null, calls: [] };
}

/**
 * ⛔ **Required on every `own` and every donor, so a fixture must carry it** — and `core` with a
 * real rate by default, because that is the published catalogues' ordinary case. The
 * not-measurable branch is the one a test has to ask for, which is the right way round: a fixture
 * that defaulted to `rate: null` would make every card under test look unmeasured.
 */
function propagation(overrides: Partial<PropagationRate> = {}): PropagationRate {
  return {
    annotation_kind: "cog_orthogroup",
    prevalence_band: "core",
    checkable_node_count: 3_253,
    unanimous_node_count: 3_248,
    one_gene_node_count: 43,
    rate: 0.9985,
    interval_low: 0.9967,
    interval_high: 0.9993,
    min_nodes: 30,
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
    ...emptyInference(),
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
describe("⭐ partial COG coverage is what the node's call COVERS", () => {
  it("says the call reaches the whole node, at the measured rate", () => {
    // ⛔ This asserted *"the remaining 40 are unlabelled rather than different — the gap this locus
    // could fill, not a gap in it"* (David, 2026-08-20). That sentence is about the ANNOTATOR's
    // method; the reader's question is what this locus is. David, 2026-10-05: "It can just be COG
    // is assigned to the whole node from the hit within it with measured accuracy > 99.5 %. That
    // is it." Rewritten rather than deleted, because it is the sentence being replaced.
    const card = mount(CogCard, { props: { propagation: propagation(), coverage: coverage(), entries: [entry("COG1132", 60)] } });
    expect(card.text()).toContain("60 of 100 genes carry a COG assignment");
    expect(card.text()).toContain("COG is assigned to the whole node from the hit within it");
    expect(card.text()).toContain("99.85 % of 3,253 comparable core syntelogues agree");
    expect(card.text()).not.toContain("unlabelled rather than different");
  });

  it("⛔ says nothing about propagation where the rate is not measurable", () => {
    // A `rare` band has no comparable node, so there is no rate — and a sentence quoting one would
    // be inventing it. The coverage line still stands on its own.
    const card = mount(CogCard, {
      props: { propagation: null, coverage: coverage(), entries: [entry("COG1132", 60)] },
    });
    expect(card.text()).toContain("60 of 100 genes carry a COG assignment");
    expect(card.text()).not.toContain("assigned to the whole node");
  });

  it("⛔ does NOT predict a label for them", () => {
    // Propagating a COG from the annotated members is annotation transfer: it needs a transfer
    // rule, a confidence measure, a decision about discordant loci and a marking that can never be
    // mistaken for a Bakta call. None of that is shipped, so neither is a predicted label.
    const card = mount(CogCard, { props: { propagation: propagation(), coverage: coverage(), entries: [entry("COG1132", 60)] } });
    const rows = card.findAll("tbody tr");
    expect(rows).toHaveLength(1);
    expect(card.text()).not.toContain("predicted");
  });

  it("omits it when every gene is annotated — there is nowhere left to propagate", () => {
    const card = mount(CogCard, {
      props: { propagation: propagation(), coverage: coverage({ cog_annotated_gene_count: 100 }), entries: [entry("COG1132", 100)] },
    });
    expect(card.text()).not.toContain("assigned to the whole node");
  });

  it("⛔ says what NO coverage means, and shows no chips at all", () => {
    const card = mount(CogCard, {
      props: { propagation: propagation(),
        coverage: coverage({ cog_annotated_gene_count: 0, cog_distinct_id_count: 0, modal_cog_categories: null }),
        entries: [],
      },
    });
    expect(card.text()).toContain("no coverage here, so it can neither support nor contradict");
    expect(card.findAll(".chip")).toHaveLength(0);
    expect(card.text()).not.toContain("assigned to the whole node");
  });
});

describe("⛔ the distinct-group count is a COUNT, never a relation", () => {
  it("calls two groups ordinary rather than a fault", () => {
    // More than one orthologous group is an ordinary consequence of grouping above the family
    // level, which is what this method does.
    const card = mount(CogCard, {
      props: { propagation: propagation(), coverage: coverage({ cog_distinct_id_count: 2 }), entries: [entry("COG1132", 40), entry("COG0642", 20)] },
    });
    expect(card.text()).toContain("2 orthologous groups");
    expect(card.text()).toContain("an ordinary consequence of grouping above the family level");
    expect(card.find(".chip.neutral").exists()).toBe(true);
  });

  it("marks a single group as the quiet win, with no explanatory note", () => {
    const card = mount(CogCard, { props: { propagation: propagation(), coverage: coverage(), entries: [entry("COG1132", 60)] } });
    expect(card.text()).toContain("one orthologous group");
    expect(card.text()).not.toContain("ordinary consequence");
    expect(card.find(".chip.win").exists()).toBe(true);
  });

  it("⛔⛔ claims NOTHING where one gene carries the COG — one voter is not a vote", () => {
    // `nunique` over a single row is 1, so the card printed a `win`-toned "one orthologous group"
    // off a single voter, on 518 ecoli / 498 kp loci. David, 2026-10-05: "The one orthologous group
    // is silly. There is one gene." The rest of the codebase already refuses this claim.
    const card = mount(CogCard, {
      props: {
        propagation: propagation(),
        coverage: coverage({ cog_annotated_gene_count: 1, cog_distinct_id_count: 1 }),
        entries: [entry("COG3206", 1)],
      },
    });
    expect(card.text()).not.toContain("orthologous group");
    expect(card.find(".chip.win").exists()).toBe(false);
    // ⭐ But the category — what the locus IS — is still shown, and now as the lead chip.
    expect(card.find(".chip.lead").exists()).toBe(true);
  });
});

describe("⚠ a COG category is one letter OR SEVERAL", () => {
  it("names each letter of a multi-letter category", () => {
    const card = mount(CogCard, {
      props: { propagation: propagation(), coverage: coverage({ modal_cog_categories: ["EP"] }), entries: [entry("COG1132", 60)] },
    });
    expect(card.text()).toContain(
      "EP — Amino acid transport and metabolism · Inorganic ion transport and metabolism",
    );
  });

  it("shows no category chip where the locus has none", () => {
    const card = mount(CogCard, {
      props: { propagation: propagation(), coverage: coverage({ modal_cog_categories: null }), entries: [entry("COG1132", 60)] },
    });
    expect(card.findAll(".chip")).toHaveLength(1);
    expect(card.find(".chip.lead").exists()).toBe(false);
  });

  it("⭐ gives the category the LEAD tone, not the quietest one the sheet has", () => {
    // David, 2026-10-05: "This is the most important thing and yet it is grey."
    const card = mount(CogCard, {
      props: { propagation: propagation(), coverage: coverage(), entries: [entry("COG1132", 60)] },
    });
    const category = card.find(".chip.lead");
    expect(category.exists()).toBe(true);
    expect(category.text()).toContain("Cell motility");
    expect(category.classes()).not.toContain("neutral");
  });
});

// ── GO ─────────────────────────────────────────────────────────────────────────────────────────
describe("⛔⛔ GO: coverage, the classes, and NO verdict chip", () => {
  function mountGo(annotated = 40, entries = [entry("GO:0016874", 40, { name: "ligase activity" })]) {
    return mount(GeneOntologyCard, {
      props: {
        namespace: "molecular_function" as const,
        annotatedGeneCount: annotated,
        geneCount: 100,
        entries,
      },
    });
  }

  it("⛔ renders NO chip at all — the verdict is retired", () => {
    // David, 2026-10-05. The chip said "classes differ" for 16 genes saying plasma membrane and one
    // saying membrane — parent and child, compared by `worst_relation`, which has no ontology. It
    // also said "one class" over two classes, because `classify_sets` returns `single` for one
    // distinct SET. ⚠ Asserted rather than deleted: this is what stops it coming back.
    expect(mountGo().find(".chip").exists()).toBe(false);
    for (const word of ["one class", "partial annotation", "classes overlap", "classes differ"]) {
      expect(mountGo().text()).not.toContain(word);
    }
  });

  it("still leads with coverage, against the LOCUS size", () => {
    expect(mountGo().text()).toContain("40 of 100");
    expect(mountGo(0, []).text()).toContain(
      "no coverage here, so it can neither support nor contradict",
    );
  });

  it("⚠ does not imply the class list is complete — it is capped at 4 per namespace", () => {
    // 688 ecoli / 269 kp loci hit `TOP_GO = 4` with nothing on the page saying so, which is how
    // `fcl` showed four identical-looking rows beside a chip claiming the classes differed.
    expect(mountGo().text()).toContain("The commonest classes here");
    expect(mountGo(0, []).text()).not.toContain("The commonest classes here");
  });

  it("⭐ shows EVERY class — the thing the chip used to contradict", () => {
    // `gumC`'s molecular function is two classes on every annotated gene and the chip said
    // "one class". The table was right all along; the chip was the wrong half.
    const card = mountGo(16, [
      entry("GO:0016301", 16, { name: "kinase activity" }),
      entry("GO:0043167", 16, { name: "ion binding" }),
    ]);
    expect(card.findAll("tbody tr")).toHaveLength(2);
    expect(card.text()).toContain("kinase activity");
    expect(card.text()).toContain("ion binding");
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
    // ⚠ The guarantee being held is that the COVERAGE side collapses to one line rather than three
    // empty namespace cards — so filter to the cards outside `.infer`, the predicted-function card.
    const coverageHeadings = tab
      .findAll(".card:not(.infer) .sub-head")
      .map((node) => node.text());
    expect(coverageHeadings.filter((heading) => heading.startsWith("GO"))).toEqual(["GO"]);
    expect(tab.text()).toContain("None of these 100 genes carries a GO term");
    // ⚠ And the predicted card, which here predicts nothing, carries NO "predicted" tag. The tag marks a
    // prediction that was made, not a card that could have made one.
    expect(tab.find(".infer").text()).toContain("Annotation-transfer predictions are not loaded for this catalogue");
    expect(tab.find(".infer .infer-tag").exists()).toBe(false);
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

  it("⚠ gives each namespace its OWN entries, never the first one three times", () => {
    // ⛔ They arrive in ONE `gene_ontology_slim` list carrying their own namespace; reading them
    // positionally would file a cellular-component class under molecular function, with three cards
    // that all look right. (This pinned the per-namespace VERDICT until 2026-10-05. The verdict is
    // retired; the per-namespace split it was really exercising is not.)
    const cards = mountTab().findAllComponents(GeneOntologyCard);
    expect(cards.map((card) => card.props("namespace"))).toEqual([
      "molecular_function",
      "biological_process",
      "cellular_component",
    ]);
    expect(cards.map((card) => card.props("entries").length)).toEqual([1, 1, 0]);
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

  it("⭐ glosses the EC code in words — and the KEGG KO NOT, which is the licence", () => {
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30, kegg_annotated_gene_count: 20 }),
        enzymeEntries: [entry("2.7.10.-", 30, { name: null })],
        keggEntries: [entry("K01886", 20, { name: null })],
      },
    });
    const [ec, kegg] = card.findAll(".kv");
    expect(ec!.find(".ec-gloss").text()).toBe("transferase · 4th level not stated");
    // ⛔⛔ The licence decision, as a gate: KEGG's terms permit linking but not redistributing, and
    // ~880 KO descriptions in a published page is redistribution. A future "fix" that glosses this
    // row the way the one above it is glossed fails here, which is the point.
    expect(kegg!.find(".ec-gloss").exists()).toBe(false);
    expect(kegg!.text()).toBe("KEGG KOK01886×2020 of 100 annotated");
  });

  it("⚠ says the two rows are different KINDS of claim", () => {
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30 }),
        enzymeEntries: [entry("6.1.1.18", 30, { name: null })],
        keggEntries: [],
      },
    });
    expect(card.text()).toContain("An EC number classifies the reaction, on four levels");
    expect(card.text()).toContain("a KEGG KO identifies a gene family");
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

  it("⭐ STATES the absent row rather than dropping it", () => {
    // ⛔ `gumC` carries an EC and no KEGG, and the card used to simply omit the KEGG line — leaving
    // a reader unable to tell "no KEGG here" from "KEGG not shown". That is the question David
    // asked of this card: *"we don't need to infer KEGG do we??? We have it for our nodes??"*
    const card = mount(EnzymeAndKeggCard, {
      props: {
        coverage: coverage({ ec_annotated_gene_count: 30 }),
        enzymeEntries: [entry("6.1.1.18", 30, { name: null })],
        keggEntries: [],
      },
    });
    expect(card.findAll(".kv")).toHaveLength(2);
    expect(card.findAll(".k").map((k) => k.text())).toEqual(["EC", "KEGG KO"]);
    expect(card.text()).toContain("No KEGG orthology on any of this node's genes");
    // ⚠ And WHY it is absent — measured over all 1,436,421 genes in the two catalogues, not asserted.
    expect(card.text()).toContain("8.2 % of genes against COG's 67.7 %");
    // The absent row carries no coverage count: "0 of 100 annotated" beside a sentence saying
    // exactly that is the same fact twice.
    expect(card.findAll(".cov")).toHaveLength(1);
  });

  it("⛔ is absent ENTIRELY where neither vocabulary has anything", () => {
    // Two headings over two absences say less than no card at all, and the GO card already states
    // the locus-wide "nothing here" case.
    const card = mount(EnzymeAndKeggCard, {
      props: { coverage: coverage({}), enzymeEntries: [], keggEntries: [] },
    });
    expect(card.find(".card").exists()).toBe(false);
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

  it("⛔ does NOT claim the slim resolves depth — it names the case where it does not", () => {
    // The footnote used to say the fold exists "because a term and its own child are annotated at
    // different depths rather than in disagreement", which reads as a claim that the fold HANDLES
    // that. It does not: `plasma membrane`, `membrane` and `outer membrane` are all three classes
    // in `goslim_metagenomics`, which is exactly how `gumC` came to read "classes differ" for 16
    // genes saying plasma membrane and one saying membrane. Until the measurement is ancestry-aware
    // the page says so itself, so a reader is not told a guarantee the data does not support.
    const text = mountTab().text();
    expect(text).toContain("selection, not a hierarchy");
    expect(text).toContain("plasma membrane");
    expect(text).toMatch(/does not yet account for that/);
    expect(text).not.toContain("rather than in disagreement");
  });
});

// ── the predicted-function card: annotation-transfer v0.1, which replaced the ladder ────────────
function inference(overrides: Partial<VocabularyInference> = {}): VocabularyInference {
  return { annotation_kind: "cog_orthogroup", own: null, ...overrides };
}

function withInference(entry: VocabularyInference): FunctionResponse {
  const base = emptyInference();
  return block({
    annotations: {},
    inference: {
      vocabularies: base.inference.vocabularies.map((row) =>
        row.annotation_kind === entry.annotation_kind ? entry : row,
      ),
    },
  });
}

function call(overrides: Partial<ServedPrediction> = {}): ServedPrediction {
  return {
    annotation_kind: "gene_ontology_slim",
    gene_ontology_namespace: "cellular_component",
    terms: [
      { value: "GO:0005886", name: "plasma membrane" },
      { value: "GO:0016020", name: "membrane" },
    ],
    estimated_accuracy: 0.9206349206349206,
    esm_medoid_similarity: 0.952113,
    bacformer_medoid_similarity: 0.772893,
    annotation_model_version: "at-v0.1",
    donor: {
      node_label: "3565",
      display_name: "ompA",
      gene_count: 72,
      catalogue_key: "kp-probe100",
      species_name: "Klebsiella pneumoniae",
      same_catalogue: true,
    },
    ...overrides,
  };
}

function predictions(calls: ServedPrediction[]): PredictionsBlock {
  return {
    state: "loaded",
    annotation_model_version: "at-v0.1",
    annotation_model_label: "annotation-transfer v0.1",
    serve_cut: 0.8,
    calls,
  };
}

function mountCard(block: PredictionsBlock) {
  return mount(PredictedFunctionCard, { props: { predictions: block } });
}

describe("⭐ the predicted-function card — annotation-transfer v0.1", () => {
  it("⭐ shows each prediction with ITS OWN estimated accuracy, and no band", () => {
    const card = mountCard(
      predictions([
        call(),
        call({ gene_ontology_namespace: "biological_process", estimated_accuracy: 0.8367,
               terms: [{ value: "GO:0006259", name: "DNA metabolic process" },
                       { value: "GO:0008152", name: "metabolic process" }] }),
      ]),
    );  // prettier-ignore
    const text = card.text();
    expect(text).toContain("GO cellular component:");
    expect(text).toContain("Estimated accuracy 92 %");
    expect(text).toContain("Estimated accuracy 83 %");
    // ⛔ No band words anywhere: David, 2026-10-28.
    expect(text).not.toMatch(/\b(high|medium|low) confidence\b/i);
    expect(text).not.toMatch(/0\.(80|85|90|95)\s*[–-]\s*0\./);
    expect(card.find(".infer-tag").text()).toBe("predicted");
  });

  it("⚠ leads with the stated GO class and lists what it implies as 'within', not as a second finding", () => {
    const card = mountCard(predictions([call()]));
    expect(card.find(".infer-value strong").text()).toBe("plasma membrane");
    expect(card.find(".predict-broader").text()).toContain("within membrane");
  });

  it("⭐ names every COG letter, as peers", () => {
    const card = mountCard(
      predictions([
        call({ annotation_kind: "cog_functional_category", gene_ontology_namespace: null,
               terms: [{ value: "E", name: null }, { value: "P", name: null }] }),
      ]),
    );  // prettier-ignore
    expect(card.text()).toContain("COG category:");
    expect(card.find(".infer-value strong").text()).toBe(
      "E — Amino acid transport and metabolism · P — Inorganic ion transport and metabolism",
    );
    expect(card.find(".predict-broader").exists()).toBe(false);
  });

  it("⛔ rounds a promise DOWN, never up — and prints the model's 1.0 as the model's number", () => {
    const text = mountCard(
      predictions([call({ estimated_accuracy: 0.995 }), call({ estimated_accuracy: 0.8 }),
                   call({ estimated_accuracy: 1.0 })]),
    ).text();  // prettier-ignore
    expect(text).toContain("Estimated accuracy 99 %");
    expect(text).toContain("Estimated accuracy 80 %");
    expect(text).toContain("Estimated accuracy 100 %");
  });

  it("walks to a donor in THIS catalogue, so a reader can go and judge it", async () => {
    const card = mountCard(predictions([call()]));
    expect(card.find("a.infer-walk").exists()).toBe(false);
    await card.find("button.infer-walk").trigger("click");
    expect(card.emitted("walk")).toEqual([["3565"]]);
    expect(card.text()).toContain("locus 3565, 72 genes");
    expect(card.text()).toContain("ESM 0.952 · Bacformer 0.773");
  });

  it("⛔ LINKS to a donor in ANOTHER catalogue by its key — labels are catalogue-private", async () => {
    const card = mountCard(
      predictions([
        call({ donor: { node_label: "3565", display_name: "ompA", gene_count: 72, catalogue_key: "ecoli-probe100",
                        species_name: "Escherichia coli", same_catalogue: false } }),
      ]),
    );  // prettier-ignore
    const link = card.find("a.infer-walk");
    expect(link.attributes("href")).toBe("?catalogue=ecoli-probe100#3565");
    expect(card.find("button.infer-walk").exists()).toBe(false);
    expect(card.text()).toContain("Escherichia coli");
    await link.trigger("click");
    expect(card.emitted("walk")).toBeUndefined();
  });

  // ⛔⛔ The five answers, each its own sentence (review of 9cfa481: three shared one). Each test asserts its
  // own sentence AND the absence of every other state's, so mapping one state onto another's words dies.
  const SENTENCES = {
    declined: "makes no prediction for this locus",
    none_produced: "The annotation model produced no predictions for this catalogue",
    not_loaded: "Annotation-transfer predictions are not loaded for this catalogue",
    unknown_model: "which this version of BacAtlas cannot show yet",
  } as const;
  function onlySentence(text: string, which: keyof typeof SENTENCES) {
    for (const [name, sentence] of Object.entries(SENTENCES)) {
      if (name === which) expect(text).toContain(sentence);
      else expect(text).not.toContain(sentence);
    }
  }

  it("⛔ a locus the model DECLINED names the model and its bar — and claims nothing about its genes", () => {
    const card = mountCard(predictions([]));
    onlySentence(card.text(), "declined");
    expect(card.text()).toContain("annotation-transfer v0.1 makes no prediction for this locus");
    expect(card.text()).toContain("estimated accuracy is 80 % or more");
    // ⛔ The retired claim: false for 426 served GO predictions on the probes.
    expect(card.text()).not.toContain("own genes");
    expect(card.find(".infer-tag").exists()).toBe(false);
  });

  it("⛔ 'none produced' is a statement about the MODEL", () => {
    const card = mountCard(noPredictions("none_produced"));
    onlySentence(card.text(), "none_produced");
    expect(card.text()).toContain("so no locus in it carries a predicted function");
  });

  it("⛔ 'not loaded' is a statement about THIS DATABASE, and says it is not one about the model", () => {
    const card = mountCard(noPredictions("not_loaded"));
    onlySentence(card.text(), "not_loaded");
    expect(card.text()).toContain("says nothing about whether the model has predictions for it");
  });

  it("⛔ a version this build does not know FAILS CLOSED, naming it — never an empty bar", () => {
    const card = mountCard({
      state: "unknown_model",
      annotation_model_version: "at-v0.2",
      annotation_model_label: "at-v0.2",
      serve_cut: null,
      calls: [],
    });
    onlySentence(card.text(), "unknown_model");
    expect(card.text()).toContain("come from at-v0.2");
    expect(card.text()).not.toMatch(/accuracy of\s+or more/);
    // ⛔ and even a calls list arriving beside a null bar is not rendered without the bar it cleared
    const guarded = mountCard({ ...predictions([call()]), serve_cut: null });
    expect(guarded.find(".infer-value").exists()).toBe(false);
    onlySentence(guarded.text(), "unknown_model");
  });

  it("⛔ never writes a bare 'v0.1' — nuna's clustering ladder has one too", () => {
    for (const block of [predictions([call()]), predictions([])]) {
      const text = mountCard(block).text();
      expect(text).toContain("annotation-transfer v0.1");
      expect(text.replaceAll("annotation-transfer v0.1", "")).not.toContain("v0.1");
    }
  });
});

describe("the predicted card's place on the tab", () => {
  it("⛔ a node with its own call in every vocabulary is offered no prediction here", () => {
    // ⚠ BOTH vocabularies, deliberately. With only COG given an `own` the EC card still renders —
    // which is correct, and would make a bare `.infer` assertion pass for the wrong reason.
    const own = {
      term: "COG1132",
      name: "x",
      gene_count: 60,
      annotated_gene_count: 60,
      member_gene_count: 100,
      checkable: true,
      propagation: propagation(),
    };
    const tab = mountTab({
      block: block({
        annotations: {},
        inference: {
          vocabularies: INFERENCE_KINDS.map((annotation_kind) => ({ annotation_kind, own })),
        },
      }),
    });
    // The predicted card still renders — here saying the catalogue has none — and predicts nothing.
    expect(tab.find(".infer .infer-value").exists()).toBe(false);
    // ⭐ But the propagation note DOES appear — the claim is made wherever a node carries a call,
    // not only where a single gene carries it. ⛔ Once per vocabulary EXCEPT COG, whose own card
    // states the rate in David's wording; listing it here as well printed the same figure in the
    // same words a few lines above that card.
    expect(tab.findAll(".infer-spread li")).toHaveLength(INFERENCE_KINDS.length - 1);
    expect(tab.find(".infer-spread").text()).not.toContain("COG — carried by");
    // … and the figure is still on the page exactly once, on the COG card itself.
    expect(tab.text()).toContain("COG is assigned to the whole node from the hit within it");
  });

  it("⭐ renders the predictions the response carries, AFTER every card of real annotation", () => {
    const tab = mountTab({ block: block({ predictions: predictions([call()]) }) });
    const cards = tab.findAll(".card");
    expect(cards.at(-1)!.classes()).toContain("predict");
    expect(tab.find(".predict").text()).toContain("plasma membrane");
  });

  it("passes a same-catalogue donor walk up to the view", async () => {
    const tab = mountTab({ block: block({ predictions: predictions([call()]) }) });
    await tab.find(".predict button.infer-walk").trigger("click");
    expect(tab.emitted("walk")).toEqual([["3565"]]);
  });
});

/*
 * ⭐⭐ The block this tab was rewritten for. It previously said a call resting on one gene meant the
 * members *"cannot be checked against each other"* — true, and the inverse of the method's central
 * claim, which is that carrying by one is enough. The tests below pin the new claim in all three
 * states; the old assertion ("agrees with itself") is replaced rather than deleted, because it is
 * the sentence being retired and a deletion would leave nothing saying so.
 */
describe("⭐ the call propagates across the whole syntelogue", () => {
  // ⛔ EC, not COG. COG states its rate on its OWN card now — David's wording, which the summary
  // block was repeating a few lines above it — so it is deliberately absent from this list. The
  // mechanics under test (the band named once, unmeasured vs 0 %, too-few-nodes, a node that
  // propagates nowhere) are the same `propagationRateSentence` either way.
  function propagationNote(own: Partial<OwnSupport>) {
    return mountTab({
      block: withInference(
        inference({
          annotation_kind: "ec_number",
          own: {
            term: "2.7.10.-",
            name: "Protein-tyrosine kinases",
            gene_count: 1,
            annotated_gene_count: 1,
            member_gene_count: 100,
            checkable: false,
            propagation: propagation(),
            ...own,
          },
        }),
      ),
    }).find(".infer-spread");
  }

  it("⛔ says a call on ONE gene APPLIES to all 100 — the claim, not a caveat", () => {
    const note = propagationNote({});
    expect(note.exists()).toBe(true);
    expect(note.text()).toContain("Each of these applies to the whole syntelogue.");
    expect(note.text()).toContain("that propagation is the method, not a gap in it");
    expect(note.text()).toContain("EC number — carried by 1 of 100 genes");
    expect(note.text()).toContain("99.85 % of 3,253 comparable syntelogues agree");
    // ⚠ The band is named ONCE, in the lead, and not on every line — the first render repeated it
    // three times on one screen.
    expect(note.text()).toContain("this catalogue's core syntelogues agree where two or more");
    // ⛔ The retired sentence, named so this test fails if it ever comes back.
    expect(note.text()).not.toContain("cannot be checked");
    expect(note.text()).not.toContain("agrees with itself");
  });

  it("⭐ says the SAME thing where seventeen genes carry it — only the count differs", () => {
    const note = propagationNote({ annotated_gene_count: 17, gene_count: 17, checkable: true });
    expect(note.text()).toContain("EC number — carried by 17 of 100 genes");
    expect(note.text()).toContain("99.85 % of 3,253 comparable syntelogues agree");
    expect(note.text()).toContain("Each of these applies to the whole syntelogue.");
  });

  it("⛔ a band with NO comparable node reads as unmeasured, never as 0 %", () => {
    const note = propagationNote({
      member_gene_count: 1,
      propagation: propagation({
        prevalence_band: "rare",
        checkable_node_count: 0,
        unanimous_node_count: 0,
        one_gene_node_count: 367,
        rate: null,
        interval_low: null,
        interval_high: null,
      }),
    });
    expect(note.text()).toContain("carried by its single gene");
    expect(note.text()).toContain("no other member to carry it to");
    expect(note.text()).toContain("no comparable syntelogue in this catalogue, so no rate");
    expect(note.text()).toContain("unmeasured, which is not the same as unreliable");
    expect(note.text()).not.toContain("0.00 %");
    expect(note.text()).not.toContain("0 %");
    // ⛔ And no sentence explaining a figure that is not there.
    expect(note.text()).not.toContain("The figure beside each");
  });

  it("⚠ too few comparable nodes says HOW few, and still quotes no rate", () => {
    const note = propagationNote({
      propagation: propagation({
        prevalence_band: "soft_core",
        checkable_node_count: 12,
        rate: null,
        interval_low: null,
        interval_high: null,
      }),
    });
    expect(note.text()).toContain("only 12 comparable syntelogues");
    expect(note.text()).toContain("fewer than the 30 needed for a rate");
    expect(note.text()).not.toContain("The figure beside each");
  });

  it("⚠ a node every gene of which carries the call propagates NOWHERE, and says so", () => {
    const note = propagationNote({ annotated_gene_count: 100, gene_count: 100, checkable: true });
    expect(note.text()).toContain("carried by every one of its 100 genes");
    expect(note.text()).toContain("nothing left to propagate to");
    expect(note.text()).not.toContain("applies to all 100");
  });
});
