import { describe, expect, it } from "vitest";

import {
  COG_CATEGORIES,
  ENZYME_CLASSES,
  GENE_ONTOLOGY_NAMESPACES,
  cogCategoryNames,
  coverageParts,
  coverageSentence,
  enzymeClassSummary,
  keggOrthologyUrl,
} from "./functionVocabulary";

describe("⛔⛔ the GO verdict ladder is RETIRED, and this is why", () => {
  it("exports no verdict ladder and no verdict function", async () => {
    // ⛔ David, 2026-10-05. The ladder was `worst_relation` — Pfam's domain-architecture comparator
    // — applied to GO slim classes. Pfam passes a CLAN into it and families never nest; GO passed
    // `null` and `goslim_metagenomics` keeps `membrane`, `plasma membrane` and `outer membrane` as
    // siblings, plus the three namespace ROOTS. So a term and its own parent scored `disjoint`,
    // the comparator's condemning verdict: `gumC` read "classes differ" for 16 genes saying plasma
    // membrane and one saying membrane, and 11 of the 43 `disjoint` verdicts across the two
    // catalogues involve a root — a term that is an ancestor of everything.
    //
    // ⚠ Asserted rather than deleted. A deletion leaves nothing saying the ladder must not come
    // back, and the next reader who sees an unused `GoVerdict` in the schema would restore it.
    const module = await import("./functionVocabulary");
    expect("GENE_ONTOLOGY_VERDICTS" in module).toBe(false);
    expect("geneOntologyVerdict" in module).toBe(false);
  });
});

describe("COG categories", () => {
  it("⚠ names EVERY letter of a multi-letter category", () => {
    // `EP`, `KT`, `NUW`, `EHJQ` are all real. "Category W" is unreadable; the name is the point.
    expect(cogCategoryNames(["EP"])).toEqual([
      "Amino acid transport and metabolism",
      "Inorganic ion transport and metabolism",
    ]);
    expect(cogCategoryNames(["NUW"])).toHaveLength(3);
  });

  it("⛔ passes an unknown letter through as ITSELF rather than dropping it", () => {
    // A vanished letter turns "EÖ" into a category the locus does not have.
    expect(cogCategoryNames(["EÖ"])).toEqual([
      "Amino acid transport and metabolism",
      "Ö",
    ]);
  });

  it("returns nothing for a locus with no category", () => {
    expect(cogCategoryNames(null)).toEqual([]);
    expect(cogCategoryNames([])).toEqual([]);
  });

  it("covers the 26 single-letter categories", () => {
    expect(Object.keys(COG_CATEGORIES)).toHaveLength(26);
  });
});

describe("⛔ coverage is stated before any verdict, against the LOCUS size", () => {
  it("says what no coverage MEANS rather than printing a zero", () => {
    // A gene the annotator never labelled says nothing either way — it can neither support nor
    // contradict the locus, which is a different statement from "the members disagree".
    expect(coverageSentence(0, 40, "a COG assignment")).toBe(
      "None of these 40 genes carries a COG assignment — no coverage here, so it can neither " +
        "support nor contradict this locus.",
    );
  });

  it("⛔ takes the share over the locus, never over the annotated subset", () => {
    // Against the annotated subset this would read "1 of 1", i.e. 100 %, where one gene in forty
    // carries a label — the single most misleading thing this tab could say.
    expect(coverageSentence(1, 40, "a GO term")).toBe("1 of 40 genes carry a GO term.");
  });
});

describe("the three namespaces", () => {
  it("are in the published page's own order", () => {
    expect(GENE_ONTOLOGY_NAMESPACES).toEqual([
      "molecular_function",
      "biological_process",
      "cellular_component",
    ]);
  });
});

describe("⛔ KEGG is linked and never named", () => {
  it("builds an entry URL — and this module deliberately holds no KO name table", () => {
    // KEGG's terms permit linking freely but not redistributing its content; embedding ~880 KO
    // descriptions in a published page is redistribution. NCBI's COG is a US Government work, which
    // is why `COG_CATEGORIES` ships names and there is no `KEGG_NAMES` beside it. The asymmetry is
    // the licence, not an oversight — the API sends `name: null` for every KO row to match.
    expect(keggOrthologyUrl("K01886")).toBe("https://www.genome.jp/entry/K01886");
  });
});

describe("⚠ the coverage COUNT is emphasised on its own", () => {
  it("keeps the count separable so the markup can bold it", () => {
    // `app.js:2980` wraps it in `<b>` — it is what a reader scans for. A `<b>` dropped here cannot
    // be recovered by any stylesheet later, so the pieces travel separately.
    expect(coverageParts(60, 100, "a COG assignment")).toEqual({
      emphasis: "60 of 100",
      rest: " genes carry a COG assignment.",
    });
  });

  it("⚠ emphasises NOTHING when there is no coverage", () => {
    // "None" is a sentence, not a count; bolding it would make absence the loudest thing on the card.
    expect(coverageParts(0, 40, "a GO term").emphasis).toBeNull();
  });

  it("reassembles into exactly the sentence", () => {
    for (const [annotated, size] of [[0, 40], [1, 40], [40, 40]] as const) {
      const parts = coverageParts(annotated, size, "a GO term");
      expect(`${parts.emphasis ?? ""}${parts.rest}`).toBe(coverageSentence(annotated, size, "a GO term"));
    }
  });
});

describe("⭐ an EC code reads as words, and KEGG deliberately does not", () => {
  it("names the class and says which fields are NOT stated", () => {
    // ⚠ 2,257 of 9,398 values in the published catalogues are dash-padded, so this is the common
    // shape rather than the edge: `gumC`'s own EC is `2.7.10.-`.
    expect(enzymeClassSummary("2.7.10.-")).toBe("transferase · 4th level not stated");
    expect(enzymeClassSummary("3.1.-.-")).toBe("hydrolase · 3rd and 4th levels not stated");
    expect(enzymeClassSummary("6.-.-.-")).toBe("ligase · only the class is stated");
  });

  it("says nothing about depth where all four fields ARE stated", () => {
    expect(enzymeClassSummary("6.1.1.18")).toBe("ligase");
    expect(enzymeClassSummary("7.1.1.2")).toBe("translocase");
  });

  it("⛔ a value is a SET: it splits on commas and takes the SHALLOWEST depth", () => {
    // The same two shapes that make `split_part` return a code neither side holds.
    expect(enzymeClassSummary("1.6.5.9,7.1.1.-")).toBe(
      "oxidoreductase · translocase · 4th level not stated",
    );
    expect(enzymeClassSummary("2.7.1.1,2.7.10.1")).toBe("transferase");
  });

  it("returns null rather than inventing a class it does not have", () => {
    expect(enzymeClassSummary("-.-.-.-")).toBeNull();
    expect(enzymeClassSummary("")).toBeNull();
    expect(enzymeClassSummary("9.1.1.1")).toBeNull();
  });

  it("⛔ has SEVEN classes and no eighth — it is the numbering, not a database", () => {
    // Vendoring ExPASy's ~8,000 names would be a CC BY-ND derivative; seven generic words are not.
    expect(Object.keys(ENZYME_CLASSES)).toEqual(["1", "2", "3", "4", "5", "6", "7"]);
    for (const name of Object.values(ENZYME_CLASSES)) expect(name.split(" ")).toHaveLength(1);
  });
});
