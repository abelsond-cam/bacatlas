/**
 * @vitest-environment jsdom
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it } from "vitest";

import type {
  AnnotationEntry,
  Locus,
  LocusSimilarity,
  NeighbourDisplayRow,
  PfamFamilyReference,
  Representation,
  SimilarityBaseline,
  UnirefFamily,
} from "@/api/types";
import type { SimilarityViewId } from "@/lib/similarityViews";

import EmbeddingSimilarityCard from "./EmbeddingSimilarityCard.vue";
import LocusHeadline from "./LocusHeadline.vue";
import SequenceDiversityCard from "./SequenceDiversityCard.vue";

/**
 * The two random GENE-pair floors, an order of magnitude apart — which is the whole reason they are
 * per representation, and the reason the card draws them at all: the same 0.41 reads oppositely.
 *
 * ⛔ These are `similarity_baselines[rep].floor_*` and NOT the retired medoid null. On the published
 * ecoli/Bacformer catalogue those sit at 0.0587 and 0.0651 — close enough to look interchangeable
 * and not be. Read from the shipped catalogue, median then p25/p75/p99.
 */
const BACFORMER_FLOOR = [0.0587, 0.0273, 0.0938, 0.2252] as const;
const ESM_FLOOR = [0.7417, 0.6014, 0.823, 0.9309] as const;

function similarity(overrides: Partial<LocusSimilarity> = {}): LocusSimilarity {
  return {
    within_similarity: 0.94,
    nearest_similarity: 0.61,
    weak_own_similarity: 0.78,
    weak_other_similarity: 0.74,
    own_neighbour_fraction: 1,
    separation_percentile: 0.62,
    weak_margin_percentile: 0.55,
    own_fraction_percentile: 0.53,
    nearest_loci: [],
    ...overrides,
  };
}

function locus(overrides: Partial<Locus> = {}): Locus {
  return {
    label: "4222",
    catalogue_ordinal: 10,
    display_name: "traC",
    display_name_source: "bakta_symbol",
    display_name_source_accession: null,
    best_product: "conjugal transfer protein TraC",
    bakta_gene_symbol: "traC",
    gene_count: 100,
    genome_count: 97,
    named_gene_count: 80,
    prevalence_band: "shell",
    median_gene_length_nt: 2400,
    gene_length_interquartile_range_nt: 120,
    uniref50: {
      family_count: 3,
      major_family_count: 2,
      labelled_gene_count: 95,
      impurity: 0.12,
      coverage: 0.95,
    },
    pfam: { annotated_gene_count: 90, architecture_count: 2, concordance_class: "nested" },
    evidence: {
      syntenic_a5: 0.91,
      collapse_tier: "mmseq@0.98",
      collapse_bucket: null,
      resolved_threshold: 0.98,
      resolved_threshold_is_capped_at_50_members: false,
    },
    similarity: { esm: similarity(), bacformer: similarity() },
    interest_score: 0.5,
    ...overrides,
  } as Locus;
}

// ── the headline ───────────────────────────────────────────────────────────────────────────────
function mountHeadline(overrides: Partial<Locus> = {}, collectionGenomeCount = 100) {
  return mount(LocusHeadline, {
    props: { locus: locus(overrides), collectionGenomeCount },
  });
}

function tileValues(wrapper: ReturnType<typeof mountHeadline>): string[] {
  return wrapper.findAll(".tile .v").map((node) => node.text());
}

function tileCaptions(wrapper: ReturnType<typeof mountHeadline>): string[] {
  return wrapper.findAll(".tile .k").map((node) => node.text());
}

describe("⭐ the headline — what it is, then how good it is", () => {
  it("names the locus, its band and its id", () => {
    const head = mountHeadline();
    expect(head.find(".locus-name").text()).toBe("traC");
    expect(head.find(".band").text()).toBe("shell");
    expect(head.find(".locus-id").text()).toBe("locus 4222");
  });

  it("⚠ shades the band from its POSITION, so the ramp stays monotonic", () => {
    expect(mountHeadline({ prevalence_band: "core" }).find(".band").attributes("style")).toContain(
      "--b: 1.00",
    );
    expect(mountHeadline({ prevalence_band: "rare" }).find(".band").attributes("style")).toContain(
      "--b: 0.00",
    );
  });

  it("⛔ renders `rare`, the band that was missing from the contract for 30 % of loci", () => {
    expect(mountHeadline({ prevalence_band: "rare" }).find(".band").text()).toBe("rare");
  });

  it("reads `soft_core` without its underscore", () => {
    expect(mountHeadline({ prevalence_band: "soft_core" }).find(".band").text()).toBe("soft core");
  });

  it("shows six tiles: three saying WHAT, three saying HOW GOOD", () => {
    expect(tileValues(mountHeadline())).toEqual([
      "97%", // 97 of 100 genomes
      "1.03", // 100 genes / 97 genomes
      "98%", // ⭐ sequence identity — `resolved_threshold`, promoted out of the prose below
      "0.91", // synteny A5
      "0.940", // within cluster · Bacformer — the median over every gene pair, shown as itself
      "0.940", // within cluster · ESM
    ]);
  });

  it("⛔ tiles `resolved_threshold` and NEVER `seqid_coverage`", () => {
    // `seqid_coverage` is populated on 0 of 17,531 ecoli-nuna4 loci, so a tile reading it would print
    // a dash on every locus in both catalogues. The identity tile must follow the threshold.
    const head = mountHeadline({
      evidence: { ...locus().evidence, resolved_threshold: 0.6 },
    });
    expect(tileValues(head)[2]).toBe("60%");
  });

  it("⚠ states the 50-member cap in the identity tile's title, because 35 % of loci sit at it", () => {
    const capped = mountHeadline({
      evidence: { ...locus().evidence, resolved_threshold_is_capped_at_50_members: true },
    });
    const title = capped.findAll(".tile")[2]!.attributes("title")!;
    expect(title).toContain("first 50 member genes");
    expect(title).toContain("not all 100");
    // …and it is NOT claimed where it does not apply.
    expect(mountHeadline().findAll(".tile")[2]!.attributes("title")).not.toContain("first 50");
  });

  it("⛔ tells the two NO-THRESHOLD cases apart — they are opposite findings, not one dash", () => {
    // Measured on ecoli-nuna4: 5,896 loci carry no threshold — 5,427 singletons (no pair to align)
    // and 469 multi-member loci MMseqs never grouped at any rung. One tooltip for both would merge
    // "does not arise" with "alignment found nothing".
    const singleton = mountHeadline({
      gene_count: 1,
      evidence: { ...locus().evidence, resolved_threshold: null },
    });
    expect(tileValues(singleton)[2]).toBe("—");
    expect(singleton.findAll(".tile")[2]!.attributes("title")).toContain("no pair to align");

    const unaligned = mountHeadline({
      evidence: { ...locus().evidence, resolved_threshold: null, collapse_tier: "esm_homology" },
    });
    const title = unaligned.findAll(".tile")[2]!.attributes("title")!;
    expect(title).toContain("never grouped these members at any threshold");
    expect(title).toContain("esm_homology");
  });

  it("⛔ prints `—` where a number was never measured, never 0", () => {
    const head = mountHeadline({
      evidence: { ...locus().evidence, syntenic_a5: null },
      similarity: {
        // ⛔ Two different absences, both of which must read `—`: ESM has no row at all, while
        // Bacformer has one whose within-cluster median is null — a singleton, measured and found
        // to have no pair. Neither is a zero.
        esm: null,
        bacformer: similarity({ within_similarity: null, separation_percentile: null }),
      },
    });
    expect(tileValues(head)).toEqual(["97%", "1.03", "98%", "—", "—", "—"]);
  });

  it("⛔ says `one member` for a singleton's cohesion, which is neither 0 nor a dash", () => {
    // `within_similarity` is null for EXACTLY the singletons — 5,427 of 17,531 ecoli-nuna4 loci and
    // 4,398 of 15,670 kp-nuna4 — and for no multi-member locus. A 0.0 would claim the members are
    // unrelated; a dash would say the question was asked and missed.
    const head = mountHeadline({
      gene_count: 1,
      similarity: {
        esm: similarity({ within_similarity: null }),
        bacformer: similarity({ within_similarity: null }),
      },
    });
    expect(tileValues(head).slice(4, 6)).toEqual(["one member", "one member"]);
    // ⚠ And it is set as prose, not in the 19px number face, or it overflows three columns.
    expect(head.findAll(".tile .v")[4]!.classes()).toContain("v-note");
  });

  it("⛔ shows the within-cluster median ITSELF, never a rescaling against the floor", () => {
    // These two tiles were `cohesion()`: `(intra − null_mean) / (1 − null_mean)`, which turned one
    // stored 0.70 into "68%" on Bacformer and "15%" on ESM — a rescaling of a distance to ONE gene,
    // captioned as if it were the locus. The label and the number have moved together, and the
    // scale a reader needs is on the card's floor strip, where it is drawn rather than folded in.
    const head = mountHeadline({
      similarity: {
        esm: similarity({ within_similarity: 0.7 }),
        bacformer: similarity({ within_similarity: 0.7 }),
      },
    });
    expect(tileValues(head).slice(4, 6)).toEqual(["0.700", "0.700"]);
  });

  it("⛔ no longer tiles cluster separation — the card below is where it lives now", () => {
    // It held the sixth slot as a percentile until 2026-10-03. `EmbeddingSimilarityCard` already drew
    // it as a `.pair.sep` row against each representation's OWN measurable count, so the tile was the
    // duplicate — and the slot went to sequence identity, which was in no tile at all.
    expect(tileCaptions(mountHeadline())).not.toContain("cluster separation");
    expect(tileCaptions(mountHeadline())).toEqual([
      "of genomes",
      "copies per genome",
      "sequence identity",
      "synteny A5",
      "within cluster · Bacformer",
      "within cluster · ESM",
    ]);
  });
});

describe("⭐ the name says what it RESTS ON — the denominator, beside the name", () => {
  it("⛔ gives gumC's case the fraction, which the page used to say nothing about", () => {
    // David, 2026-10-03: "Only 2 / 100 gumC nodes have Bakta names. Yet this is very quiet on the
    // current page." kp-nuna4 locus 722 — 100 members, named_gene_count 2 — and `bakta_symbol`
    // returned null, so the line was absent entirely.
    const head = mountHeadline({ gene_count: 100, named_gene_count: 2, display_name: "gumC" });
    expect(head.find(".inferred-note").text()).toBe(
      "Name inferred from 2 of 100 named genes — the other 98 carry no name",
    );
  });

  it("says so where every member gene carries a name — 7,132 of 8,245 ecoli loci", () => {
    // "No note" is too weak a signal to carry 86 % of the named loci: silence would then mean both
    // "all 100 agree" and "2 of 100", which is the ambiguity this line exists to remove.
    expect(mountHeadline({ gene_count: 100, named_gene_count: 100 }).find(".inferred-note").text()).toBe(
      "Name carried by all 100 member genes",
    );
  });

  it("reads `gene` in the singular where exactly one member is named", () => {
    expect(mountHeadline({ gene_count: 40, named_gene_count: 1 }).find(".inferred-note").text()).toBe(
      "Name inferred from 1 of 40 named gene — the other 39 carry no name",
    );
  });

  it("⛔ states the ZERO for the three sources measured to always have it", () => {
    // `label`, `product` and `pfam_architecture` have named_gene_count = 0 in EVERY case across
    // ecoli-nuna4, ecoli-nuna5, kp-nuna4 and kp-nuna5 — so each can state it flatly rather than hedge.
    expect(mountHeadline({ display_name_source: "product", named_gene_count: 0 }).find(".inferred-note").text()).toBe(
      "No member gene carries a name — inferred from the Bakta product",
    );
    const pfam = mountHeadline({
      display_name_source: "pfam_architecture",
      display_name_source_accession: "PF00126",
      named_gene_count: 0,
    });
    expect(pfam.find(".inferred-note").text()).toBe(
      "No member gene carries a name — inferred from PF00126",
    );
  });

  it("⛔ breaks the silence on the LARGEST group — the 8,741 loci named only by their id", () => {
    // These got no line at all, which is how the biggest group on the page said the least. The name
    // IS the locus id, so there is no name to caveat — but there is a fact to state.
    const head = mountHeadline({ display_name_source: "label", display_name: "4222", named_gene_count: 0 });
    expect(head.find(".inferred-note").text()).toBe(
      "No member gene carries a name — this is the locus id",
    );
  });
});

// ── sequence diversity ─────────────────────────────────────────────────────────────────────────
function family(overrides: Partial<UnirefFamily> = {}): UnirefFamily {
  return {
    rank: 0,
    uniref50_accession: "UniRef50_P76362",
    gene_count: 60,
    modal_product: "LysR family transcriptional regulator",
    modal_architecture: "PF00126,PF03466",
    pfam_annotated_gene_count: 60,
    modal_symbol: "lysR",
    distinct_symbol_count: 1,
    named_gene_count: 60,
    ...overrides,
  };
}

const PFAM_REFERENCE: Record<string, PfamFamilyReference> = {
  PF00126: {
    short_name: "HTH_1",
    description: "Bacterial regulatory helix-turn-helix protein",
    interpro_accession: "IPR000847",
    interpro_name: "HTH transcriptional regulator LysR",
    clan_accession: "CL0123",
    clan_name: "HTH",
  },
  PF03466: {
    short_name: "LysR_substrate",
    description: "LysR substrate binding domain",
    interpro_accession: "",
    clan_accession: "",
    clan_name: "",
    interpro_name: "",
  },
};

function symbol(term: string, count: number, rank = 0): AnnotationEntry {
  return { rank, term, name: null, gene_count: count };
}

function mountDiversity(
  overrides: {
    locus?: Partial<Locus>;
    families?: readonly UnirefFamily[];
    symbols?: readonly AnnotationEntry[];
    listedArchitectureCount?: number;
  } = {},
) {
  return mount(SequenceDiversityCard, {
    props: {
      locus: locus(overrides.locus ?? {}),
      families: overrides.families ?? [family(), family({ rank: 1, uniref50_accession: "UniRef50_Q664A4", gene_count: 30 })],
      symbols: overrides.symbols ?? [symbol("lysR", 80)],
      pfamReference: PFAM_REFERENCE,
      listedArchitectureCount: overrides.listedArchitectureCount ?? 2,
    },
  });
}

describe("⭐ families and names are ONE table", () => {
  it("claims with the MAJOR-family count, not every straggler", () => {
    // The table lists every family; claiming with that number would overstate the split at every
    // locus with a long tail.
    expect(mountDiversity().find(".lede").text()).toContain("UniRef50 files these genes under 2 families");
  });

  it("says so plainly where the locus and the sequence family agree", () => {
    const card = mountDiversity({
      locus: { uniref50: { ...locus().uniref50, family_count: 1, major_family_count: 1 } },
      families: [family()],
    });
    expect(card.find(".lede").text()).toBe(
      "One UniRef50 family — this locus and the sequence family agree.",
    );
  });

  it("⛔ takes every share over the LOCUS SIZE, not the sum of the listed families", () => {
    // 60 + 30 of 100 member genes. Summing the rows would print 67 % / 33 % and claim the table is
    // the whole locus.
    expect(mountDiversity().findAll("td.n:not(.share)").map((n) => n.text())).toEqual(["60%", "30%"]);
  });

  it("links each accession to UniProt and chips its architecture by NAME", () => {
    const card = mountDiversity();
    expect(card.find("a.acc-link").attributes("href")).toBe(
      "https://www.uniprot.org/uniref/UniRef50_P76362",
    );
    expect(card.findAll(".chip.pfam").map((n) => n.text())).toEqual([
      "HTH_1",
      "LysR_substrate",
      "HTH_1",
      "LysR_substrate",
    ]);
  });

  it("⭐ links a chip to INTERPRO where there is one, and to Pfam where there is not", () => {
    const chips = mountDiversity().findAll(".chip.pfam");
    expect(chips[0]!.attributes("href")).toBe(
      "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000847/",
    );
    expect(chips[1]!.attributes("href")).toBe("https://www.ebi.ac.uk/interpro/entry/pfam/PF03466/");
  });

  it("⚠ shows an UNRESOLVED accession as itself rather than dropping the chip", () => {
    // A dropped chip would say the locus has no such domain, which is a different, false claim.
    const card = mountDiversity({ families: [family({ modal_architecture: "PF99999" })] });
    expect(card.find(".chip.pfam").text()).toBe("PF99999");
  });

  it("⚠ marks a family split across gene names with a COUNT, never a list", () => {
    const card = mountDiversity({
      families: [family({ distinct_symbol_count: 3 })],
    });
    const marker = card.find(".fam-sym .alt-n");
    expect(marker.text()).toBe("+2");
    expect(marker.attributes("title")).toContain("3 different gene names among this family's 60 genes");
  });

  it("⛔ shows NO `+N` marker where the count was never measured", () => {
    // `null` is not measured and `1` is measured agreement — two different facts, and neither is a
    // disagreement. A `> 1` test on a null reads as false by accident rather than by decision.
    expect(mountDiversity({ families: [family({ distinct_symbol_count: null })] })
      .find(".fam-sym .alt-n").exists()).toBe(false);
    expect(mountDiversity({ families: [family({ distinct_symbol_count: 1 })] })
      .find(".fam-sym .alt-n").exists()).toBe(false);
  });

  it("⚠ a family with a measured coverage but NO architecture still earns the Pfam column", () => {
    // Its cell reads "no domain" plus its coverage, which is a statement. Dropping the column would
    // turn "0 of 60 annotated" into silence.
    const card = mountDiversity({
      families: [family({ modal_architecture: null, pfam_annotated_gene_count: 0 })],
    });
    expect(card.findAll("th").map((n) => n.text())).toContain("Pfam");
    expect(card.find(".fam-pf .none").text()).toBe("no domain");
    expect(card.find(".fam-cov").text()).toBe("0/60 annotated");
  });

  it("omits the optional columns it has not earned", () => {
    const card = mountDiversity({
      families: [
        family({ modal_symbol: null, modal_architecture: null, pfam_annotated_gene_count: null }),
      ],
    });
    expect(card.findAll("th").map((n) => n.text())).toEqual([
      "UniRef50",
      "modal Bakta product",
      "genes",
      "share",
    ]);
  });
});

describe("the naming line", () => {
  it("names how many genes the locus is naming", () => {
    expect(mountDiversity().findAll(".lede")[1]!.text()).toContain(
      "20 of 100 genes here carry no gene name. This locus names them traC",
    );
  });

  it("says when there is no symbol at all", () => {
    expect(mountDiversity({ locus: { named_gene_count: 0 } }).findAll(".lede")[1]!.text()).toBe(
      "None of these 100 genes carries a gene name — a coherent locus with no symbol at all.",
    );
  });

  it("says when every gene already has one", () => {
    expect(mountDiversity({ locus: { named_gene_count: 100 } }).findAll(".lede")[1]!.text()).toBe(
      "Every one of these 100 genes already carries a gene name.",
    );
  });

  it("counts the names to reconcile only when there is more than one", () => {
    expect(mountDiversity({ symbols: [symbol("a", 5), symbol("b", 3)] }).text()).toContain(
      "2 different gene names are in use across this one locus, the names to reconcile",
    );
    expect(mountDiversity().text()).not.toContain("names to reconcile");
  });
});

describe("⛔ the two verdicts, both READ and neither re-derived", () => {
  it("chips the collapse tier and the Pfam class side by side", () => {
    const chips = mountDiversity().findAll(".chip-row .chip");
    expect(chips.map((n) => n.text())).toEqual([
      "holds together at 98% identity",
      "partial annotation",
    ]);
  });

  it("⚠ calls `nested` partial annotation and says so under the chip", () => {
    expect(mountDiversity().find(".verdict-note").text()).toContain(
      "Absent evidence, not conflicting evidence",
    );
  });

  it("⛔ shows NO Pfam chip where no gene carries a domain — a verdict needs coverage", () => {
    const card = mountDiversity({
      locus: { pfam: { annotated_gene_count: 0, architecture_count: 0, concordance_class: "single" } },
    });
    expect(card.findAll(".chip-row .chip").map((n) => n.text())).toEqual([
      "holds together at 98% identity",
    ]);
    expect(card.find(".verdict-note").exists()).toBe(false);
  });
});

describe("⛔ Pfam coverage is stated AS COVERAGE", () => {
  it("says the unannotated remainder is absent evidence, not a competing architecture", () => {
    expect(mountDiversity().findAll(".muted").at(-1)!.text()).toBe(
      "90 of 100 genes carry a Pfam-A domain; the other 10 carry none at all — absent annotation, " +
        "not a competing architecture · 2 architectures across the locus.",
    );
  });

  it("says Pfam can neither support nor contradict where it has no coverage", () => {
    const card = mountDiversity({
      locus: { pfam: { annotated_gene_count: 0, architecture_count: 0, concordance_class: null } },
    });
    expect(card.findAll(".muted").at(-1)!.text()).toBe(
      "None of these 100 genes carries a Pfam-A domain — Pfam has no coverage here, so it can " +
        "neither support nor contradict this locus.",
    );
  });

  it("⚠ names how many architectures are LISTED when the locus has more", () => {
    const card = mountDiversity({
      locus: { pfam: { annotated_gene_count: 100, architecture_count: 9, concordance_class: "nested" } },
      listedArchitectureCount: 5,
    });
    expect(card.findAll(".muted").at(-1)!.text()).toBe(
      "All 100 genes carry a Pfam-A domain · 9 architectures across the locus (5 listed).",
    );
  });

  it("⛔ shows a family's own coverage only when it is SHORT, and a measured zero counts", () => {
    // `v-if` on the number would hide a family where NO gene is annotated — the loudest case.
    const none = mountDiversity({
      families: [family({ pfam_annotated_gene_count: 0, modal_architecture: null })],
    });
    expect(none.find(".fam-cov").text()).toBe("0/60 annotated");
    expect(mountDiversity().find(".fam-cov").exists()).toBe(false);
  });
});

// ── embedding similarity ───────────────────────────────────────────────────────────────────────
/** The two floors, an order of magnitude apart — which is the whole reason they are per rep. */
function baseline(representation: Representation): SimilarityBaseline {
  const floor = representation === "esm" ? ESM_FLOOR : BACFORMER_FLOOR;
  return {
    representation,
    form: "raw",
    floor_median: floor[0],
    floor_p25: floor[1],
    floor_p75: floor[2],
    floor_p99: floor[3],
    measurable_locus_count: 12_104,
    neighbour_knn_k: 200,
  };
}

const BASELINES = [baseline("bacformer"), baseline("esm")];

/** Two real fan-out rows, so a nearest-locus ordinal has something to resolve through. */
const NEIGHBOURS: NeighbourDisplayRow[] = [
  {
    label: "1065",
    catalogue_ordinal: 1065,
    display_name: "rfaL",
    best_product: "O-antigen ligase",
    display_name_source: "bakta_symbol",
    genome_count: 90,
    median_gene_length_nt: 1230,
    prevalence_band: "core",
  },
  {
    label: "2404",
    catalogue_ordinal: 2404,
    display_name: "waaL",
    best_product: null,
    display_name_source: "bakta_symbol",
    genome_count: 20,
    median_gene_length_nt: 900,
    prevalence_band: "shell",
  },
];

function mountSimilarity(
  overrides: {
    similarity?: Partial<Record<Representation, LocusSimilarity | null>>;
    geneCount?: number;
    baselines?: readonly SimilarityBaseline[];
    view?: SimilarityViewId;
    neighbours?: readonly NeighbourDisplayRow[];
  } = {},
) {
  return mount(EmbeddingSimilarityCard, {
    props: {
      similarity: {
        bacformer:
          overrides.similarity && "bacformer" in overrides.similarity
            ? (overrides.similarity.bacformer ?? null)
            : similarity(),
        esm:
          overrides.similarity && "esm" in overrides.similarity
            ? (overrides.similarity.esm ?? null)
            : similarity(),
      },
      baselines: overrides.baselines ?? BASELINES,
      neighbours: overrides.neighbours ?? NEIGHBOURS,
      geneCount: overrides.geneCount ?? 100,
      view: overrides.view ?? "median",
    },
  });
}

function rowValues(card: ReturnType<typeof mountSimilarity>): string[] {
  return card.findAll(".pair:not(.sep) .val").map((node) => node.text());
}

describe("⛔ embedding similarity is consistency, never corroboration", () => {
  it("says so, in the card's own footer", () => {
    expect(mountSimilarity().find(".sim-note").text()).toContain(
      "reported as consistency and diagnosis — never as corroboration",
    );
  });

  it("⚠ leads with Bacformer, the context axis the track is built on", () => {
    expect(mountSimilarity().findAll(".sub-head").map((n) => n.text())).toEqual(["Bacformer", "ESM"]);
  });

  it("⚠ stacks BOTH representations inside the view — they are not tabs", () => {
    // The two disagree about which loci are weak (of 1,059 ecoli loci Bacformer flags, ESM rescues
    // 817), and a tab would hide exactly that disagreement behind a click.
    expect(rowValues(mountSimilarity())).toEqual(["0.940", "0.610", "0.940", "0.610"]);
  });

  it("⛔ prints a SIMILARITY, with nothing left to convert from a distance", () => {
    // The fixture's 0.940 is the stored number. A surviving `1 − d` would print 0.060.
    expect(rowValues(mountSimilarity())).not.toContain("0.060");
  });

  it("names the denominator of the percentile, both halves of the sentence", () => {
    expect(mountSimilarity().findAll(".sep-pct")[0]!.text()).toBe("p62 of 12,104 loci");
  });

  it("⛔ subtracts the difference HERE, so it cannot drift from the two numbers beside it", () => {
    // `separation` is not served: 0.940 − 0.610 printed from the same two numbers the rows show.
    const card = mountSimilarity();
    expect(card.findAll(".pair.sep .lab").map((n) => n.text())).toEqual(["separation", "separation"]);
    expect(card.findAll(".pair.sep .val")[0]!.text()).toBe("+0.330");
  });
});

describe("⭐ three views, and they REPLACE each other rather than stacking", () => {
  it("offers exactly three, with the median pair on by default", () => {
    const strip = mountSimilarity().findAll(".sim-view");
    expect(strip.map((n) => n.text())).toEqual(["median pair", "weakest member", "nearest neighbour"]);
    expect(strip.filter((n) => n.classes("on")).map((n) => n.text())).toEqual(["median pair"]);
  });

  it("asks for a view rather than switching itself — the choice is the reader's and survives a walk", () => {
    const card = mountSimilarity();
    card.findAll(".sim-view")[1]!.trigger("click");
    expect(card.emitted("selectView")).toEqual([["weak"]]);
  });

  it("⛔ the weakest-member view REPLACES the median rows, never adds to them", () => {
    const card = mountSimilarity({ view: "weak" });
    expect(card.findAll(".pair:not(.sep) .lab").map((n) => n.text())).toEqual([
      "weakest → own",
      "weakest → other",
      "weakest → own",
      "weakest → other",
    ]);
    expect(rowValues(card)).toEqual(["0.780", "0.740", "0.780", "0.740"]);
    // …and it is ranked on ITS OWN midrank, not on separation's.
    expect(card.findAll(".pair.sep .lab")[0]!.text()).toBe("margin");
    expect(card.findAll(".sep-pct")[0]!.text()).toBe("p55 of 12,104 loci");
    expect(card.findAll(".pair.sep .val")[0]!.text()).toBe("+0.040");
  });

  it("⛔ the nearest-neighbour view shows ONE row, because it is one share", () => {
    const card = mountSimilarity({ view: "own" });
    expect(card.findAll(".pair:not(.sep) .lab").map((n) => n.text())).toEqual([
      "nearest gene is own",
      "nearest gene is own",
    ]);
  });
});

describe("⛔⛔ the own fraction is a SHARE, and below 1.0 is the entire point of it", () => {
  it("⛔ never renders at 0 dp: 0.9999 must not print as 100%", () => {
    // 0.9999 of a large locus is a real member whose nearest gene in the whole species belongs to
    // another locus. `pct(v, 0)` renders that as "100%", which is the opposite claim.
    const card = mountSimilarity({
      view: "own",
      similarity: {
        bacformer: similarity({ own_neighbour_fraction: 0.9999, own_fraction_percentile: 0.03 }),
        esm: null,
      },
    });
    expect(rowValues(card)).toEqual(["<100%"]);
  });

  it("keeps one decimal where a whole percent would round two different loci together", () => {
    const card = mountSimilarity({
      geneCount: 143,
      view: "own",
      similarity: {
        bacformer: similarity({ own_neighbour_fraction: 0.951049, own_fraction_percentile: 0.09 }),
        esm: null,
      },
    });
    expect(rowValues(card)).toEqual(["95.1%"]);
    // ⛔ The rank row does not repeat the share — it prints what the share IMPLIES and does not say.
    expect(card.find(".pair.sep .val").text()).toBe("7 out");
    expect(card.find(".pair.sep .lab").text()).toBe("rank");
  });

  it("⛔⛔ a share of EXACTLY 1.0 gets NO percentile row at all", () => {
    // 88.5 % of loci sit at 1.0, so the midrank inside that tie block prints "p53" — reading as
    // *better than half the catalogue* when it means *tied with nearly all of it*.
    const card = mountSimilarity({ view: "own" });
    expect(rowValues(card)).toEqual(["100%", "100%"]);
    expect(card.findAll(".pair.sep")).toHaveLength(0);
  });

  it("…and shows it again the moment the share drops below 1.0", () => {
    const card = mountSimilarity({
      view: "own",
      similarity: {
        bacformer: similarity({ own_neighbour_fraction: 0.97, own_fraction_percentile: 0.08 }),
        esm: null,
      },
    });
    expect(card.findAll(".sep-pct")[0]!.text()).toBe("p8 of 12,104 loci");
  });
});

describe("⛔ not measurable is a SENTENCE, never a zero", () => {
  /**
   * ⚠ A singleton keeps its `nearest_similarity` — that question is well posed for one gene — and
   * loses every other number. 5,427 of *E. coli*'s 17,531 loci are single genes.
   */
  function singleton(): LocusSimilarity {
    return similarity({
      within_similarity: null,
      weak_own_similarity: null,
      weak_other_similarity: null,
      own_neighbour_fraction: null,
      separation_percentile: null,
      weak_margin_percentile: null,
      own_fraction_percentile: null,
    });
  }

  it("⛔ a single gene is told WHY, rather than getting a strip with nothing between it and the caption", () => {
    const card = mountSimilarity({
      geneCount: 1,
      similarity: { bacformer: singleton(), esm: singleton() },
    });
    expect(card.find(".sim-body .muted").text()).toContain(
      "A single gene has no pair inside its locus",
    );
    // ⛔ And no rails, no floor strip and no rank row — an absence is not a measurement of zero.
    expect(card.findAll(".pair")).toHaveLength(0);
    expect(card.find(".nullstrip").exists()).toBe(false);
  });

  it("drops the reason where the locus is not a singleton", () => {
    const card = mountSimilarity({
      geneCount: 40,
      similarity: { bacformer: singleton(), esm: singleton() },
    });
    expect(card.find(".sim-body .muted").text()).toBe("Not measurable for this locus.");
  });

  it("⛔ says NOT MEASURABLE in words where the midrank is missing, and `—` for the difference", () => {
    const card = mountSimilarity({
      similarity: {
        bacformer: similarity({ separation_percentile: null }),
        esm: similarity({ separation_percentile: null }),
      },
    });
    expect(card.findAll(".sep-pct")[0]!.text()).toBe("not measurable");
    expect(card.findAll(".pair.sep .val")[0]!.text()).toBe("—");
  });

  it("renders no card at all where neither representation has a row", () => {
    expect(
      mountSimilarity({ similarity: { bacformer: null, esm: null } }).find(".card").exists(),
    ).toBe(false);
  });

  it("⚠ still renders the card where only the OTHER representation is measurable", () => {
    // The two disagree about which loci are weak; losing the card would lose ESM's answer with it.
    const card = mountSimilarity({ similarity: { bacformer: null } });
    expect(card.findAll(".sub-head").map((n) => n.text())).toEqual(["ESM"]);
  });
});

describe("⭐ the floor strip is what makes a raw cosine readable", () => {
  it("draws the p25–p75 box, its whisker to p99 and the median line", () => {
    const strip = mountSimilarity().findAll(".nullstrip")[0]!;
    // ⚠ jsdom normalises `2.7%` but not the width, so both edges are asserted.
    expect(strip.find(".fbox i.q").attributes("style")).toContain("left: 2.7%");
    expect(strip.find(".fbox i.w").attributes("style")).toContain("left: 9.4%");
    expect(strip.find(".fbox i.m").attributes("style")).toContain("left: 5.9%");
  });

  it("⛔ draws a BOX, not a density — the floor is summarised by quartiles, not binned", () => {
    // A smooth hump would be a picture making a claim the artifact does not support.
    expect(mountSimilarity().find(".nullstrip .nd").exists()).toBe(false);
  });

  it("names the floor, which differs by an order of magnitude between the two", () => {
    // ⛔ The random GENE-pair floor, and NOT the retired medoid null: on the published ecoli
    // catalogue those sit at 0.0587 and 0.0651 — close enough to look interchangeable and not be.
    expect(mountSimilarity().findAll(".nsfoot span:first-child").map((n) => n.text())).toEqual([
      "random gene pairs 0.059 (p25–p75 0.03–0.09)",
      "random gene pairs 0.742 (p25–p75 0.60–0.82)",
    ]);
  });

  it("⚠ places a tick to sub-percent precision — two loci 0.5 % apart must not coincide", () => {
    const strip = mountSimilarity({
      similarity: {
        bacformer: similarity({ within_similarity: 0.935, nearest_similarity: 0.612 }),
      },
    }).findAll(".nullstrip")[0]!;
    expect(strip.find(".nt.intra").attributes("style")).toContain("left: 93.5%");
    expect(strip.find(".nt.inter").attributes("style")).toContain("left: 61.2%");
  });

  it("⛔ draws NO axis on the own view — a share of members is not a cosine", () => {
    // Drawing it against the random gene-pair box would invite reading one as the other.
    expect(mountSimilarity({ view: "own" }).find(".nullstrip").exists()).toBe(false);
    expect(mountSimilarity({ view: "weak" }).findAll(".nullstrip")).toHaveLength(2);
  });

  it("⛔ draws no strip where the floor was never measured", () => {
    const card = mountSimilarity({ baselines: [] });
    expect(card.find(".nullstrip").exists()).toBe(false);
    // …and the rails and the difference row are still there: an absent floor is not an absent
    // measurement. Its denominator is what goes, and the sentence says so.
    expect(card.findAll(".pair:not(.sep)")).toHaveLength(4);
    expect(card.findAll(".sep-pct")[0]!.text()).toBe("p62 of the measurable loci");
  });
});

describe("⭐ the five nearest other loci — the evidence for `nearest other cluster`", () => {
  const NEAREST = [
    { rank: 1, catalogue_ordinal: 1065, cross_similarity: 0.371015 },
    { rank: 2, catalogue_ordinal: 2404, cross_similarity: 0.364022 },
  ];

  function withNearest(view: SimilarityViewId = "median") {
    return mountSimilarity({
      view,
      similarity: {
        bacformer: similarity({ nearest_loci: NEAREST }),
        esm: null,
      },
    });
  }

  it("⛔⛔ resolves each catalogue ordinal through neighbour_display_rows", () => {
    // A slot code carries `catalogue_ordinal * 2 + strand` and an occupant carries a surrogate
    // locus id: both are small integers over the same range, so the wrong index names one locus
    // where another belongs, on a list that still looks entirely right (`2b99bb4`).
    expect(withNearest().findAll(".sim-row .nm").map((n) => n.text())).toEqual(["rfaL", "waaL"]);
  });

  it("shows the product beside the name, because the locus NUMBER says nothing", () => {
    const rows = withNearest().findAll(".sim-row");
    expect(rows[0]!.find(".desc").text()).toBe("O-antigen ligase");
    expect(rows[1]!.find(".desc").exists()).toBe(false);
    expect(rows.map((row) => row.find(".cos").text())).toEqual(["0.371", "0.364"]);
  });

  it("walks to a neighbour by its LABEL, not by the ordinal it was addressed with", () => {
    const card = withNearest();
    card.findAll(".sim-row")[1]!.trigger("click");
    expect(card.emitted("walk")).toEqual([["2404"]]);
  });

  it("⛔ only on the median view, whose row it is the evidence for", () => {
    expect(withNearest("weak").find(".sim-key").exists()).toBe(false);
    expect(withNearest("own").find(".sim-key").exists()).toBe(false);
  });

  it("⚠ drops an unresolvable ordinal rather than drawing a blank, walkable row", () => {
    const card = mountSimilarity({
      neighbours: [NEIGHBOURS[0]!],
      similarity: { bacformer: similarity({ nearest_loci: NEAREST }), esm: null },
    });
    expect(card.findAll(".sim-row")).toHaveLength(1);
  });

  it("⚠ shows nothing at all where the shortlist is empty — it is RAGGED, never padded", () => {
    expect(mountSimilarity().find(".sim-key").exists()).toBe(false);
  });
});

// ── the vote that names the locus ─────────────────────────────────────────────────────
describe("⭐ the counts the NAME came from", () => {
  /** kp locus 3992 as it really is: 97 genes, 7 named, `mviN` 5 against `rfbX` 2. */
  function mviN() {
    return mountDiversity({
      locus: { gene_count: 97, named_gene_count: 7, display_name: "mviN", bakta_gene_symbol: "mviN" },
      symbols: [symbol("mviN", 5), symbol("rfbX", 2, 1)],
    });
  }

  it("⭐ shows each NAME's gene count — the numbers the locus was named on", () => {
    // ⛔ The bug this answers: the only counts beside a gene name were the FAMILY gene counts in
    // the table below, so this locus read as "9 genes are rfbX and 4 are mviN" off two family rows.
    const vote = mviN().find(".sym-vote");
    expect(vote.text()).toContain("mviN 5");
    expect(vote.text()).toContain("rfbX 2");
  });

  it("⛔ counts the UNNAMED genes as their own remainder — no name is not a rival name", () => {
    expect(mviN().find(".sym-vote").text()).toContain("90 unnamed");
  });

  it("marks the winner, and only the winner", () => {
    const won = mviN().findAll(".sym-vote-item.won");
    expect(won).toHaveLength(1);
    expect(won[0]!.text()).toContain("mviN");
  });

  it("⛔ keeps the two remainders apart: outside the top-N cut, and named at all", () => {
    // 80 named, 60 of them listed — so 20 carry a name this card does not show, and 20 more carry
    // none. One combined "other 40" would let missing annotation read as a rival name.
    const vote = mountDiversity({
      locus: { gene_count: 100, named_gene_count: 80 },
      symbols: [symbol("lysR", 40), symbol("ybeF", 20, 1)],
    }).find(".sym-vote");
    expect(vote.text()).toContain("20 under names not listed");
    expect(vote.text()).toContain("20 unnamed");
  });

  it("⛔ says NOTHING where one name covers every gene — the headline already said it", () => {
    expect(
      mountDiversity({
        locus: { gene_count: 100, named_gene_count: 100 },
        symbols: [symbol("lysR", 100)],
      }).find(".sym-vote").exists(),
    ).toBe(false);
  });

  it("⚠ names a TIE as a tie — 29 ecoli and 12 kp loci were settled by the alphabet", () => {
    const card = mountDiversity({
      locus: { gene_count: 100, named_gene_count: 60, display_name: "aaaA", bakta_gene_symbol: "aaaA" },
      symbols: [symbol("aaaA", 30), symbol("zzzZ", 30, 1)],
    });
    expect(card.find(".sym-vote-note").text()).toContain("alphabetical order alone");
  });

  it("⛔ marks no winner where the name did NOT come from this vote", () => {
    // An inferred name is a display fallback. Marking a symbol as its source would dress a Pfam
    // short name as a Bakta annotation — the one thing the inferred-name rule must never do.
    const card = mountDiversity({
      locus: { gene_count: 97, named_gene_count: 7, display_name_source: "pfam_architecture" },
      symbols: [symbol("mviN", 5), symbol("rfbX", 2, 1)],
    });
    expect(card.findAll(".sym-vote-item.won")).toHaveLength(0);
    expect(card.find(".sym-vote-note").exists()).toBe(false);
  });
});

// ── the family symbol column's denominator ────────────────────────────────────────────
describe("⛔ a family's modal gene name travels with its coverage", () => {
  it("says how many of the family's genes are NAMED where that differs from its size", () => {
    // kp 3992's 9-gene family names TWO genes `rfbX`. It reported `distinct_symbol_count = 1` — one
    // distinct name, no "+N" marker — and read as nine genes agreeing.
    const card = mountDiversity({
      families: [family({ gene_count: 9, named_gene_count: 2, modal_symbol: "rfbX", distinct_symbol_count: 1 })],
    });
    expect(card.find(".fam-sym .fam-cov").text()).toBe("2/9 named");
    expect(card.find(".fam-sym .alt-n").exists()).toBe(false);
  });

  it("stays silent where every gene in the family carries a name", () => {
    expect(
      mountDiversity({ families: [family({ gene_count: 60, named_gene_count: 60 })] })
        .find(".fam-sym .fam-cov").exists(),
    ).toBe(false);
  });

  it("⛔ `null` is NOT MEASURED and shows nothing; a measured ZERO is the loudest case", () => {
    expect(
      mountDiversity({ families: [family({ named_gene_count: null })] })
        .find(".fam-sym .fam-cov").exists(),
    ).toBe(false);
    // 0 of 60 named, with a modal symbol that cannot exist — the row shows its coverage rather
    // than a bare name, which is what a `v-if` on the NUMBER would have hidden.
    expect(
      mountDiversity({ families: [family({ named_gene_count: 0 })] })
        .find(".fam-sym .fam-cov").text(),
    ).toBe("0/60 named");
  });

  it("quotes the denominator in the +N marker only where it says something", () => {
    const split = mountDiversity({
      families: [family({ gene_count: 9, named_gene_count: 2, distinct_symbol_count: 2 })],
    });
    expect(split.find(".fam-sym .alt-n").attributes("title")).toContain("2 named genes of 9");
    const whole = mountDiversity({ families: [family({ distinct_symbol_count: 3 })] });
    expect(whole.find(".fam-sym .alt-n").attributes("title")).toContain("this family's 60 genes");
  });
});
