/**
 * @vitest-environment jsdom
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it } from "vitest";

import type {
  AnnotationEntry,
  CalibrationCell,
  CalibrationLadder,
  FunctionResponse,
  InferenceKind,
  OwnSupport,
  PropagationRate,
  VocabularyInference,
  WalkStep,
} from "@/api/types";
import { INFERENCE_KINDS } from "@/api/types";

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

/**
 * ⛔ The inference block is REQUIRED on a `FunctionResponse`, so a fixture must carry it. Empty by
 * default: a node that needs no suggestion is the ordinary case, and a fixture that silently
 * shipped a candidate would make every card under test look like an inferred one.
 */
function emptyInference(): Pick<FunctionResponse, "inference" | "calibration"> {
  return {
    inference: {
      representation: "esm",
      vocabularies: INFERENCE_KINDS.map((annotation_kind) => ({
        annotation_kind,
        own: null,
        walk: [],
        candidate: null,
      })),
    },
    // ⛔ Derived from INFERENCE_KINDS, never a hand-written pair: adding a vocabulary must not
    // require editing fixtures, and a Record missing a key fails the build rather than the test.
    calibration: Object.fromEntries(
      INFERENCE_KINDS.map((kind) => [kind, ladder(kind)]),
    ) as FunctionResponse["calibration"],
  };
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

function ladder(annotation_kind: InferenceKind): CalibrationLadder {
  return {
    annotation_kind,
    representation: "esm",
    annotated_locus_count: 0,
    locus_count: 0,
    min_pairs: 30,
    calling_floor: 0.8,
    cells: [],
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
    // ⚠ TWO GO headings is correct, and the distinction is the point: one card reports what the
    // locus's own genes say (and here, that none of them say anything), the other offers what a
    // neighbour suggests. The guarantee being held is that the COVERAGE side collapses to one line
    // rather than three empty namespace cards — so filter to the cards outside `.infer`.
    const coverageHeadings = tab
      .findAll(".card:not(.infer) .sub-head")
      .map((node) => node.text());
    expect(coverageHeadings.filter((heading) => heading.startsWith("GO"))).toEqual(["GO"]);
    expect(tab.text()).toContain("None of these 100 genes carries a GO term");
    // ⚠ And the other GO heading is the inference card's, which here suggests nothing — so it
    // carries NO "inferred" tag. The tag marks a suggestion that was made, not a card that could
    // have made one; announcing one over a paragraph withholding it is the bug this guards.
    expect(tab.find(".infer").text()).toContain("nothing near enough to infer from");
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

// ── the inferred-function card: the three states a fixture of real bytes cannot all reach ────────
function cell(overrides: Partial<CalibrationCell> = {}): CalibrationCell {
  return {
    tier: "0.98-0.99",
    facet: null,
    level: 2,
    level_label: "COG orthogroup",
    pairs: 770,
    agreement: 0.908,
    interval_low: 0.886,
    interval_high: 0.926,
    chance: 0.0024,
    lift: 381.9,
    ...overrides,
  };
}

function step(rank: number, cosine: number, carries: boolean): WalkStep {
  return {
    rank,
    cosine,
    tier: "0.98-0.99",
    node_label: `n${rank}`,
    catalogue_ordinal: rank,
    display_name: `node ${rank}`,
    carries_annotation: carries,
  };
}

function inference(overrides: Partial<VocabularyInference> = {}): VocabularyInference {
  return { annotation_kind: "cog_orthogroup", own: null, walk: [], candidate: null, ...overrides };
}

function withInference(entry: VocabularyInference): FunctionResponse {
  const base = emptyInference();
  return block({
    annotations: {},
    inference: {
      representation: "esm",
      vocabularies: base.inference.vocabularies.map((row) =>
        row.annotation_kind === entry.annotation_kind ? entry : row,
      ),
    },
    calibration: {
      ...(Object.fromEntries(
        INFERENCE_KINDS.map((kind) => [kind, ladder(kind)]),
      ) as FunctionResponse["calibration"]),
      cog_orthogroup: {
        ...ladder("cog_orthogroup"),
        cells: [cell(), cell({ tier: "0.90-0.94", agreement: 0.336, pairs: 8525, lift: 4.9 })],
      },
    },
  });
}

describe("the inferred-function card", () => {
  it("⭐ states the suggestion, the depth, the measured rate AND its n", () => {
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.999, false), step(2, 0.985, true)],
          candidate: {
            rank: 2,
            cosine: 0.985,
            tier: "0.98-0.99",
            donor: {
              node_label: "n2",
              catalogue_ordinal: 2,
              display_name: "wza",
              term: "COG0450",
              name: "peroxiredoxin",
              gene_count: 97,
              annotated_gene_count: 97,
              member_gene_count: 99,
              checkable: true,
              propagation: propagation(),
            },
            offers: [
              { facet: null, level: 2, suggested: true, value: ["COG0450"], calibration: cell() },
            ],
            suggested_count: 1,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.find(".infer-value").text()).toContain("COG0450");
    const rate = card.find(".infer-rate").text();
    expect(rate).toContain("90.8%");
    expect(rate).toContain("770");
    // ⛔ The LIFT left this sentence (David, 2026-10-05). "agree 24.3 % of the time, which is 145×
    // more often than a random annotated node" is two true numbers telling opposite stories, and
    // the reader carries away the second. Lift stays in the ladder table below, where comparing a
    // 0.2 % baseline with a 16 % one is the whole point.
    expect(rate).not.toContain("×");
    expect(rate).not.toContain("more often than a random");
    expect(card.find(".infer-table").text()).toContain("382×");
    // ⛔ the rank that carried nothing is named — a card showing only the donor would imply rank 1
    expect(card.text()).toContain("rank 1 at 0.999");
  });

  it("⛔⛔ refuses a CLOSE donor on the measured rate, and says which bar it missed", () => {
    // David, 2026-10-05, on `gumC`: the page suggested a KEGG orthology from `wza`, a rank-5
    // neighbour at cosine 0.9293, where kp KEGG agreement is 24.3 % over 481 pairs — wrong about
    // three times in four. ⛔ "Too remote to call" would be false here: 0.929 is close. It is the
    // measurement that refuses, and the reader is owed the number that did it.
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.9293, true)],
          candidate: {
            rank: 1,
            cosine: 0.9293,
            tier: "0.90-0.94",
            donor: {
              node_label: "3511",
              catalogue_ordinal: 3511,
              display_name: "wza",
              term: "K01991",
              name: null,
              gene_count: 96,
              annotated_gene_count: 96,
              member_gene_count: 99,
              checkable: true,
              propagation: propagation(),
            },
            // ⛔ The refused rung is OFFERED, carrying the cell that refused it. An empty list would
            // mean "the donor states nothing readable", which is a different and far weaker finding.
            offers: [
              {
                facet: null,
                level: 1,
                suggested: false,
                value: ["K01991"],
                calibration: cell({
                  tier: "0.90-0.94",
                  level: 1,
                  level_label: "KEGG orthology",
                  agreement: 0.243,
                  pairs: 481,
                }),
              },
            ],
            suggested_count: 0,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.find(".infer-value").exists()).toBe(false);
    expect(card.text()).toContain("agree on their KEGG orthology only 24.3% of the time");
    expect(card.text()).toContain("481 measured pairs");
    expect(card.text()).toContain("below the 80% this page needs");
    expect(card.text()).toContain("KEGG orthology \u2014 not suggested.");
    expect(card.text()).not.toContain("too remote to call");
    // ⭐ and the donor is still named, so a reader can go and judge it
    expect(card.text()).toContain("wza");
    // ⛔ but the card must not announce a suggestion it is withholding, nor argue the donor's case
    expect(card.text()).toContain("— nothing suggested");
    expect(card.text()).not.toContain("inferred from a neighbour");
    expect(card.find(".infer-tag").exists()).toBe(false);
    expect(card.text()).not.toContain("counts as a donor at all");
  });

  it("⛔ says NOTHING is suggested where the only donor is too remote to call", () => {
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.84, true)],
          candidate: {
            rank: 1,
            cosine: 0.84,
            tier: "< 0.90",
            donor: {
              node_label: "n1",
              catalogue_ordinal: 1,
              display_name: null,
              term: "COG0450",
              name: null,
              gene_count: 4,
              annotated_gene_count: 4,
              member_gene_count: 4,
              checkable: true,
              propagation: propagation(),
            },
            offers: [],
            suggested_count: 0,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.text()).toContain("too remote to call");
    expect(card.find(".infer-value").exists()).toBe(false);
  });

  it("⭐ states the DARK SET as a finding where no neighbour carries the vocabulary", () => {
    const tab = mountTab({
      block: withInference(
        inference({ walk: [1, 2, 3, 4, 5].map((rank) => step(rank, 0.97, false)) }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.text()).toContain("None of the 5 nearest nodes");
    expect(card.find(".infer-value").exists()).toBe(false);
    // ⚠ Three of the four vocabularies take "a" and the fourth takes "an", so a hardcoded article
    // rendered "carries a EC number" on exactly one card.
    const everyCard = tab.findAll(".infer").map((one) => one.text());
    expect(everyCard.some((one) => one.includes("carries an EC number either"))).toBe(true);
    expect(everyCard.some((one) => one.includes("carries a COG either"))).toBe(true);
    expect(everyCard.join(" ")).not.toContain("a EC number");
  });

  it("⭐ says the DONOR's own call covers the donor's whole node, at the donor's rate", () => {
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.985, true)],
          candidate: {
            rank: 1,
            cosine: 0.985,
            tier: "0.98-0.99",
            donor: {
              node_label: "n1",
              catalogue_ordinal: 1,
              display_name: null,
              term: "COG4733",
              name: null,
              gene_count: 1,
              annotated_gene_count: 1,
              member_gene_count: 120,
              checkable: false,
              propagation: propagation(),
            },
            offers: [
              { facet: null, level: 2, suggested: true, value: ["COG4733"], calibration: cell() },
            ],
            suggested_count: 1,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.text()).toContain("carried by 1 of its 120 genes and covers that whole node");
    expect(card.text()).toContain("which is why a neighbour resting on one annotated gene is still a donor");
    // ⚠ That clause is for the ONE-gene donor only; a donor resting on 57 of 100 gets the general
    // rule instead, because the striking case read as a non-sequitur beside it.
    expect(card.text()).not.toContain("the rule by which it counts as a donor at all");
    // ⚠ The donor card DOES name the band: it quotes the donor's, and has no lead to carry it.
    expect(card.text()).toContain("99.85 % of 3,253 comparable core syntelogues agree.");
    // ⛔ The retired claim, named so it cannot come back: a one-gene donor is not a weaker donor,
    // which is the whole reason a donor "carries" a vocabulary by the any-gene rule.
    expect(card.text()).not.toContain("never checked against itself");
  });

  it("⛔ a donor EVERY gene of which carries the call says nothing about propagation", () => {
    // ⚠ There is no propagation on that side to report, and a sentence saying so on every card
    // would be noise — this is the common case (a 100-of-100 donor).
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.985, true)],
          candidate: {
            rank: 1,
            cosine: 0.985,
            tier: "0.98-0.99",
            donor: {
              node_label: "n1",
              catalogue_ordinal: 1,
              display_name: null,
              term: "COG4733",
              name: null,
              gene_count: 100,
              annotated_gene_count: 100,
              member_gene_count: 100,
              checkable: true,
              propagation: propagation(),
            },
            offers: [
              { facet: null, level: 2, suggested: true, value: ["COG4733"], calibration: cell() },
            ],
            suggested_count: 1,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.text()).toContain("COG4733");
    expect(card.text()).not.toContain("covers that whole node");
    expect(card.text()).not.toContain("comparable core syntelogues agree");
  });

  it("⚠ a donor carried by MANY of its genes gets the general rule, not the one-gene case", () => {
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.985, true)],
          candidate: {
            rank: 1,
            cosine: 0.985,
            tier: "0.98-0.99",
            donor: {
              node_label: "n1",
              catalogue_ordinal: 1,
              display_name: null,
              term: "COG4733",
              name: null,
              gene_count: 57,
              annotated_gene_count: 57,
              member_gene_count: 100,
              checkable: true,
              propagation: propagation(),
            },
            offers: [
              { facet: null, level: 2, suggested: true, value: ["COG4733"], calibration: cell() },
            ],
            suggested_count: 1,
          },
        }),
      ),
    });
    const card = tab.find(".infer");
    expect(card.text()).toContain("carried by 57 of its 100 genes and covers that whole node");
    expect(card.text()).toContain("the rule by which it counts as a donor at all");
    expect(card.text()).not.toContain("resting on one annotated gene");
  });

  it("⛔ a node with its own call is offered no suggestion at all", () => {
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
          representation: "esm",
          vocabularies: INFERENCE_KINDS.map((annotation_kind) => ({
            annotation_kind,
            own,
            walk: [],
            candidate: null,
          })),
        },
      }),
    });
    expect(tab.find(".infer").exists()).toBe(false);
    // ⭐ But the propagation note DOES appear — the claim is made wherever a node carries a call,
    // not only where a single gene carries it. ⛔ Once per vocabulary EXCEPT COG, whose own card
    // states the rate in David's wording; listing it here as well printed the same figure in the
    // same words a few lines above that card.
    expect(tab.findAll(".infer-spread li")).toHaveLength(INFERENCE_KINDS.length - 1);
    expect(tab.find(".infer-spread").text()).not.toContain("COG — carried by");
    // … and the figure is still on the page exactly once, on the COG card itself.
    expect(tab.text()).toContain("COG is assigned to the whole node from the hit within it");
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

  it("walks to the donor, so a reader can go and judge it", async () => {
    const tab = mountTab({
      block: withInference(
        inference({
          walk: [step(1, 0.985, true)],
          candidate: {
            rank: 1,
            cosine: 0.985,
            tier: "0.98-0.99",
            donor: {
              node_label: "4242",
              catalogue_ordinal: 7,
              display_name: "wza",
              term: "COG0450",
              name: null,
              gene_count: 9,
              annotated_gene_count: 9,
              member_gene_count: 9,
              checkable: true,
              propagation: propagation(),
            },
            offers: [
              { facet: null, level: 2, suggested: true, value: ["COG0450"], calibration: cell() },
            ],
            suggested_count: 1,
          },
        }),
      ),
    });
    await tab.find(".infer-walk").trigger("click");
    expect(tab.emitted("walk")).toEqual([["4242"]]);
  });
});
