/**
 * The Function tab's vocabularies — COG categories, the GO verdict ladder, and where each
 * accession is looked up.
 *
 * ⭐ **Every one of these is a *link*, and that is a licensing decision as much as a design one.**
 * KEGG's terms permit linking freely but not redistributing its content, so a KO id is shown and
 * **never named** — embedding ~880 KO descriptions in a page is redistribution. NCBI's COG is a US
 * Government work, so its names ship. The asymmetry is deliberate; do not "fix" it by adding KEGG
 * names to the reference tables.
 */

import type { GeneOntologyNamespace, GoVerdict, OwnSupport, PropagationRate } from "@/api/types";
import { prevalenceBandLabel } from "./prevalence";
import type { Verdict } from "./evidenceVocabulary";

/** The three namespaces, in the order the published page lays them out. */
export const GENE_ONTOLOGY_NAMESPACES = [
  "molecular_function",
  "biological_process",
  "cellular_component",
] as const satisfies readonly GeneOntologyNamespace[];

export const GENE_ONTOLOGY_NAMESPACE_LABEL: Readonly<Record<GeneOntologyNamespace, string>> = {
  molecular_function: "molecular function",
  biological_process: "biological process",
  cellular_component: "cellular component",
};

/**
 * ⭐ **The same ladder Pfam uses, because it means the same thing** — but worded for an ontology:
 * "architecture" is wrong for a GO term, so `single` reads *one class*, not *one architecture*.
 *
 * ⛔ `no_coverage` is a VALUE and is deliberately absent from this table: it is not a verdict, it is
 * the absence of one, and giving it a chip would put "no coverage" in the same visual position as
 * "classes differ". The coverage line above it already says so, in words, first.
 */
/*
 * ⚠ Every `note` is empty, and that is faithful rather than unfinished: the published page renders
 * a GO verdict as a chip and nothing else (`app.js:3113`). The Pfam ladder carries notes because a
 * Pfam verdict is the §6.2 conflict argument; a GO verdict sits above a coverage line that has
 * already said the thing a note would say.
 */
export const GENE_ONTOLOGY_VERDICTS: Readonly<Record<string, Verdict>> = {
  single: { tone: "win", label: "one class", note: "" },
  same_domains: { tone: "win", label: "one class", note: "" },
  nested: { tone: "neutral", label: "partial annotation", note: "" },
  overlapping: { tone: "warn", label: "classes overlap", note: "" },
  disjoint: { tone: "bad", label: "classes differ", note: "" },
};

/**
 * The chip for a namespace's verdict, or `null` where there is nothing to say.
 *
 * ⛔ Returns `null` for `no_coverage` **and** for an unrecognised value, rather than inventing a
 * neutral chip: a verdict this page does not understand must not be rendered as one it does.
 */
export function geneOntologyVerdict(verdict: GoVerdict | null): Verdict | null {
  if (verdict === null || verdict === "no_coverage") return null;
  return GENE_ONTOLOGY_VERDICTS[verdict] ?? null;
}

/**
 * The 26 COG functional categories.
 *
 * ⚠ **A locus's category is one letter OR SEVERAL** — `EP`, `KT`, `NUW`, `EHJQ` are all real — so
 * every letter is named separately and an unknown letter passes through as itself rather than
 * vanishing. "Category W" is unreadable; "W — Extracellular structures" is the reason to show it.
 */
export const COG_CATEGORIES: Readonly<Record<string, string>> = {
  J: "Translation, ribosomal structure and biogenesis",
  A: "RNA processing and modification",
  K: "Transcription",
  L: "Replication, recombination and repair",
  B: "Chromatin structure and dynamics",
  D: "Cell cycle control, cell division, chromosome partitioning",
  Y: "Nuclear structure",
  V: "Defense mechanisms",
  T: "Signal transduction mechanisms",
  M: "Cell wall/membrane/envelope biogenesis",
  N: "Cell motility",
  Z: "Cytoskeleton",
  W: "Extracellular structures",
  U: "Intracellular trafficking, secretion, and vesicular transport",
  O: "Posttranslational modification, protein turnover, chaperones",
  X: "Mobilome: prophages, transposons",
  C: "Energy production and conversion",
  G: "Carbohydrate transport and metabolism",
  E: "Amino acid transport and metabolism",
  F: "Nucleotide transport and metabolism",
  H: "Coenzyme transport and metabolism",
  I: "Lipid transport and metabolism",
  P: "Inorganic ion transport and metabolism",
  Q: "Secondary metabolites biosynthesis, transport and catabolism",
  R: "General function prediction only",
  S: "Function unknown",
};

/** `"EHJQ"` → the four names, in order, with unknown letters kept as themselves. */
export function cogCategoryNames(categories: readonly string[] | null): string[] {
  if (categories === null) return [];
  return categories.flatMap((group) =>
    [...group].map((letter) => COG_CATEGORIES[letter] ?? letter),
  );
}

/**
 * ⭐ **The seven top-level EC classes — and ONLY the seven.**
 *
 * `term_name` is populated for 100 % of this catalogue's COG and GO entries and **0 %** of its EC
 * and KEGG ones, so an EC code reaches the page as `2.7.10.-` and nothing else. Naming the first
 * field costs seven generic words and turns that into *"transferase"*, which is the difference
 * between a reader recognising the locus and not.
 *
 * ⛔ **The deeper fields stay as digits, and that is a LICENCE decision.** ExPASy's ENZYME is
 * CC BY-ND; folding its ~8,000 names onto our four rungs is plausibly a derivative, so it is not
 * vendored. Seven class names are not a database — they are the definition of the numbering.
 */
export const ENZYME_CLASSES: Readonly<Record<string, string>> = {
  "1": "oxidoreductase",
  "2": "transferase",
  "3": "hydrolase",
  "4": "lyase",
  "5": "isomerase",
  "6": "ligase",
  "7": "translocase",
};

/** How many of an EC code's four fields are actually stated — `3.1.-.-` states two. */
function statedDepth(code: string): number {
  const fields = code.split(".");
  let depth = 0;
  for (const field of fields.slice(0, 4)) {
    if (!/^\d+$/.test(field)) break;
    depth += 1;
  }
  return depth;
}

/**
 * ⭐ An EC value in words: its class, and which of its four fields are not stated.
 *
 * ⛔ **A value is a SET, not a string.** 1,034 of 9,398 are comma-joined lists and 2,257 carry `-`
 * placeholders, so this splits on commas and reads each code's own depth — the same two shapes that
 * make `split_part` return a code neither side holds. `null` where no field is a class digit, which
 * is the honest answer for a value this page does not understand.
 */
export function enzymeClassSummary(value: string): string | null {
  const codes = value
    .split(",")
    .map((one) => one.trim())
    .filter(Boolean);
  const names = [...new Set(codes.map((code) => ENZYME_CLASSES[code.split(".")[0] ?? ""]))].filter(
    (name): name is string => name !== undefined,
  );
  if (names.length === 0) return null;
  const depth = Math.min(...codes.map(statedDepth));
  const missing =
    depth >= 4
      ? null
      : depth === 3
        ? "4th level not stated"
        : depth === 2
          ? "3rd and 4th levels not stated"
          : "only the class is stated";
  return missing === null ? names.join(" · ") : `${names.join(" · ")} · ${missing}`;
}

export function cogEntryUrl(accession: string): string {
  return `https://www.ncbi.nlm.nih.gov/research/cog/cog/${encodeURIComponent(accession)}/`;
}

export function geneOntologyTermUrl(accession: string): string {
  return `https://amigo.geneontology.org/amigo/term/${encodeURIComponent(accession)}`;
}

export function enzymeCommissionUrl(accession: string): string {
  return `https://enzyme.expasy.org/EC/${encodeURIComponent(accession)}`;
}

/** ⚠ Linked, never named — see the module docstring. */
export function keggOrthologyUrl(accession: string): string {
  return `https://www.genome.jp/entry/${encodeURIComponent(accession)}`;
}

/**
 * ⛔ **Coverage, stated BEFORE any verdict, and against the LOCUS size.**
 *
 * A share taken against the annotated subset reads 100 % where one gene in forty carries a label,
 * which is the single most misleading thing this tab could say. And *no coverage* is its own
 * sentence: a gene the annotator never labelled says nothing either way, so it can neither support
 * nor contradict the locus — which is a different statement from "the members disagree".
 */
export function coverageParts(
  annotated: number,
  geneCount: number,
  what: string,
): { emphasis: string | null; rest: string } {
  if (annotated === 0) {
    return {
      // ⚠ Nothing to emphasise: the published page bolds the COUNT, and "none" is a sentence rather
      // than a count. Bolding it would make absence the loudest thing on the card.
      emphasis: null,
      rest:
        `None of these ${geneCount} genes carries ${what} — no coverage here, so it can neither ` +
        "support nor contradict this locus.",
    };
  }
  // ⚠ The count is emphasised on its own (`app.js:2980` wraps it in `<b>`), because it is what a
  // reader scans for. Kept as a separate piece rather than folded into the sentence so the markup
  // can carry it — a `<b>` that is dropped here cannot be recovered by any stylesheet later.
  return { emphasis: `${annotated} of ${geneCount}`, rest: ` genes carry ${what}.` };
}

export function coverageSentence(annotated: number, geneCount: number, what: string): string {
  const parts = coverageParts(annotated, geneCount, what);
  return `${parts.emphasis ?? ""}${parts.rest}`;
}

/**
 * ⭐ **How many of a node's genes carry the call, and what follows for the ones that do not.**
 *
 * ⛔ **Three cases, and the third is not a weaker version of the others.** A node every one of whose
 * genes already carries the call propagates nothing — there is nowhere to propagate to — and a
 * one-member node has no other gene at all. "Applies to all 1" there would be true and useless,
 * which is how a page starts sounding like it is reciting rather than reporting.
 *
 * ⚠ **Deliberately terse: the WARRANT is stated once, in the block's lead**, not on every line. The
 * first render repeated *"that propagation is the method, not a gap in it"* three times on one
 * screen, which reads as insistence rather than explanation.
 */
export function propagationCarriage(own: OwnSupport): string {
  const members = own.member_gene_count.toLocaleString();
  if (own.member_gene_count <= 1) {
    return "carried by its single gene — a one-genome syntelogue has no other member to carry it to";
  }
  if (own.annotated_gene_count >= own.member_gene_count) {
    return `carried by every one of its ${members} genes — nothing left to propagate to`;
  }
  return `carried by ${own.annotated_gene_count.toLocaleString()} of ${members} genes`;
}

/**
 * The measured rate the propagation is quoted at, or **why there is none**.
 *
 * ⛔ The absence of a rate is a sentence, never a `0 %`. `rare` has zero comparable nodes in both
 * published catalogues while holding most of the one-gene ones, so this branch is the common case
 * there rather than an edge: *unmeasured* and *unreliable* are different claims and the page may
 * not let one read as the other.
 *
 * ⚠ `namesBand` is false where the band has already been named once for the whole block, and true
 * on the inferred-function card, which quotes the DONOR's band with no lead to carry it.
 */
export function propagationRateSentence(
  propagation: PropagationRate,
  { namesBand = true }: { namesBand?: boolean } = {},
): string {
  const band = namesBand ? `${prevalenceBandLabel(propagation.prevalence_band)} ` : "";
  const comparable = propagation.checkable_node_count.toLocaleString();
  if (propagation.rate === null) {
    return propagation.checkable_node_count === 0
      ? `no comparable ${band}syntelogue in this catalogue, so no rate — unmeasured, which is not ` +
          "the same as unreliable"
      : `only ${comparable} comparable ${band}syntelogues, fewer than the ${propagation.min_nodes} ` +
          "needed for a rate — unmeasured, which is not the same as unreliable";
  }
  return (
    `${(propagation.rate * 100).toFixed(2)} % of ${comparable} comparable ` +
    `${band}syntelogues agree`
  );
}
