/**
 * @vitest-environment jsdom
 */
import { mount } from "@vue/test-utils";
import { describe, expect, it, vi } from "vitest";

import type { GeneSequence, GeneSequenceResponse } from "@/api/types";

import CopySequenceButton from "./CopySequenceButton.vue";
import GeneSequenceCard from "./GeneSequenceCard.vue";
import SequenceBlock from "./SequenceBlock.vue";
import SequenceTab from "./SequenceTab.vue";

function gene(overrides: Partial<GeneSequence> = {}): GeneSequence {
  return {
    copy_ordinal: 1,
    copy_count: 1,
    flat_index: 3581,
    contig_name: "contig00083",
    seqid: "SAMEA1.contig00083",
    start_position: 1488,
    end_position: 1976,
    strand: "+",
    strand_is_observed: true,
    gff_phase: 0,
    is_five_prime_partial: false,
    length_nt: 489,
    protein_length_aa: 162,
    gene_symbol: "yfcQ",
    product: "fimbrial-like protein",
    locus_tag: "ONFFFB_19145",
    coding_sequence: "ATG" + "AAA".repeat(161) + "TAG",
    protein_sequence: "M" + "K".repeat(161),
    upstream_flank_sequence: "G".repeat(100),
    downstream_flank_sequence: "C".repeat(100),
    upstream_flank_is_truncated_by_contig_end: false,
    downstream_flank_is_truncated_by_contig_end: false,
    upstream_flank_span: [1388, 1487],
    downstream_flank_span: [1977, 2076],
    gc_percent: 46.4,
    ...overrides,
  };
}

function response(genes: readonly GeneSequence[]): GeneSequenceResponse {
  return { genome: { sample_id: "SAMEA1" }, locus: { label: "2811" }, genes };
}

function mountTab(props: Partial<InstanceType<typeof SequenceTab>["$props"]> = {}) {
  return mount(SequenceTab, {
    props: {
      displayName: "yfcQ",
      locusLabel: "2811",
      sampleId: "SAMEA1",
      response: response([gene()]),
      status: "ready" as const,
      ...props,
    },
  });
}

describe("⛔ the flanking blocks are NOT the gene, and the card says so", () => {
  it("carries the warning in its own paragraph, above the blocks", () => {
    // A guide designed against flanking DNA believing it to be coding is a wasted experiment, so
    // this cannot rely on the block headings alone.
    const card = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    expect(card.find(".seq-warn").text()).toContain("are NOT part of the gene");
    expect(card.find(".seq-warn").text()).toContain("only the middle block is the coding sequence");
  });

  it("⚠ repeats it in every flank HEADING too", () => {
    const headings = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } })
      .findAllComponents(SequenceBlock)
      .map((block) => block.props("heading"));
    expect(headings[0]).toBe("100 bases upstream — not the gene");
    expect(headings[2]).toBe("100 bases downstream — not the gene");
  });
});

describe("⛔ a flank names the contig coordinates it came from", () => {
  it("says reverse-complemented on a MINUS-strand gene, and not on a plus one", () => {
    // On a minus gene the upstream flank sits at HIGHER contig coordinates. The server does the
    // orienting; the card must report it rather than implying the bases are as the contig reads.
    const minus = mount(GeneSequenceCard, {
      props: {
        gene: gene({ strand: "-", upstream_flank_span: [1977, 2076], downstream_flank_span: [1388, 1487] }),
        flankLength: 100,
      },
    });
    const provenances = minus.findAllComponents(SequenceBlock).map((b) => b.props("provenance"));
    expect(provenances[0]).toBe("Contig 1,977–2,076, reverse-complemented");
    expect(provenances[2]).toBe("Contig 1,388–1,487, reverse-complemented");

    const plus = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    expect(plus.findAllComponents(SequenceBlock)[0]!.props("provenance")).toBe("Contig 1,388–1,487");
  });

  it("⚠ prints NOTHING rather than a span for an empty flank", () => {
    const card = mount(GeneSequenceCard, {
      props: {
        gene: gene({ upstream_flank_sequence: "", upstream_flank_span: null, upstream_flank_is_truncated_by_contig_end: true }),
        flankLength: 100,
      },
    });
    expect(card.findAllComponents(SequenceBlock)[0]!.props("provenance")).toBe("");
  });

  it("says the coding block includes the stop codon", () => {
    const card = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    expect(card.findAllComponents(SequenceBlock)[1]!.props("provenance")).toContain(
      "Includes the stop codon",
    );
  });
});

describe("⚠ a gene at a contig end says its flank is short", () => {
  it("names which end, and says it is ordinary", () => {
    // Draft assemblies: mean contig 14.4 kb, ~14 genes. A short flank that does not say it is short
    // reads as a complete one.
    const card = mount(GeneSequenceCard, {
      props: {
        gene: gene({ upstream_flank_is_truncated_by_contig_end: true, upstream_flank_sequence: "G".repeat(12) }),
        flankLength: 100,
      },
    });
    const edge = card.find(".seq-edge").text();
    expect(edge).toContain("within 100 bases of the end of contig00083");
    expect(edge).toContain("upstream flank is short");
    expect(edge).toContain("common rather than exceptional");
  });

  it("pluralises when BOTH ends run off", () => {
    const card = mount(GeneSequenceCard, {
      props: {
        gene: gene({
          upstream_flank_is_truncated_by_contig_end: true,
          downstream_flank_is_truncated_by_contig_end: true,
        }),
        flankLength: 100,
      },
    });
    expect(card.find(".seq-edge").text()).toContain("upstream and downstream flanks are short");
  });

  it("says nothing at all when neither is short", () => {
    const card = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    expect(card.find(".seq-edge").exists()).toBe(false);
  });
});

describe("the coordinates the card states", () => {
  it("⛔ names the CONTIG, never an index", () => {
    const card = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    const rows = card.findAll(".seq-stats dd").map((node) => node.text());
    expect(rows[0]).toBe("contig00083");
  });

  it("states the span as 1-based inclusive, and the direction in words", () => {
    const card = mount(GeneSequenceCard, { props: { gene: gene({ strand: "-" }), flankLength: 100 } });
    const text = card.find(".seq-stats").text();
    expect(text).toContain("1,488–1,976 (1-based, inclusive)");
    expect(text).toContain("− (reverse)");
    expect(text).toContain("read as the reverse complement of the contig");
  });

  it("⚠ notes a 5′-partial CDS, because its first residue is not a start codon", () => {
    const card = mount(GeneSequenceCard, {
      props: { gene: gene({ is_five_prime_partial: true }), flankLength: 100 },
    });
    expect(card.find(".seq-stats").text()).toContain("does not begin at a start codon");
  });

  it("⛔ says when the strand was ASSUMED — silence there would read as an observation", () => {
    const observed = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    expect(observed.find(".seq-stats").text()).not.toContain("assumed forward");

    const assumed = mount(GeneSequenceCard, {
      props: { gene: gene({ strand_is_observed: false }), flankLength: 100 },
    });
    expect(assumed.find(".seq-stats").text()).toContain(
      "the strand was not recorded for this genome and has been assumed forward",
    );
  });
});

describe("⭐ rho > 1 shows every copy, numbered", () => {
  it("heads each card with `copy n of m`", () => {
    const tab = mountTab({
      response: response([
        gene({ copy_ordinal: 1, copy_count: 2, flat_index: 10 }),
        gene({ copy_ordinal: 2, copy_count: 2, flat_index: 11 }),
      ]),
    });
    expect(tab.findAll(".seq-copy-h").map((node) => node.text())).toEqual([
      "Copy 1 of 2 at this locus",
      "Copy 2 of 2 at this locus",
    ]);
    expect(tab.text()).toContain("2 copies at this locus");
  });

  it("⚠ shows NO copy heading for a single gene — `copy 1 of 1` is noise", () => {
    expect(mountTab().find(".seq-copy-h").exists()).toBe(false);
  });
});

describe("⛔ four states, and no two of them look alike", () => {
  it("asks for a genome rather than choosing one", () => {
    // A sequence belongs to one genome; "some member's" is not an answer to "what is this gene".
    const tab = mountTab({ sampleId: null, response: null, status: "idle" });
    expect(tab.text()).toContain("Anchor a genome to read its DNA here");
    expect(tab.findAllComponents(GeneSequenceCard)).toHaveLength(0);
  });

  it("says it is reading rather than showing an empty panel", () => {
    const tab = mountTab({ response: null, status: "pending" });
    expect(tab.text()).toContain("Reading the annotation file…");
  });

  it("⛔ renders a FAILURE as a failure, with the server's sentence and a retry", async () => {
    const tab = mountTab({ response: null, status: "failed", failureDetail: "contig not in the file" });
    expect(tab.find(".pop-error").text()).toContain("contig not in the file");
    await tab.find(".pop-retry").trigger("click");
    expect(tab.emitted("retry")).toHaveLength(1);
  });

  it("⛔⛔ 'no gene here' is an ANSWER and reads nothing like a failure", () => {
    // The whole distinction the tab turns on. Both sentences exist and they are different
    // sentences: one is a fact about this genome, the other is a fact about this request.
    const tab = mountTab({ response: response([]) });
    expect(tab.text()).toContain("SAMEA1 has no gene at this locus — it is not one of its members");
    expect(tab.find(".pop-error").exists()).toBe(false);
    expect(tab.text()).not.toContain("did not load");
  });
});

describe("⚠ copying says whether it worked", () => {
  it("reports a failure rather than doing nothing visible", async () => {
    // Clipboard access is refused on an insecure origin and absent in jsdom; a button that silently
    // does nothing sends a reader to the bench with an empty clipboard and no idea.
    const block = mount(SequenceBlock, {
      props: { kind: "cds", heading: "h", provenance: "p", sequence: "ATG", what: "the gene" },
    });
    await block.find(".seq-copy").trigger("click");
    await block.vm.$nextTick();
    expect(block.find(".seq-copy").text()).toBe("copy failed");
  });

  it("says `copied` when the clipboard takes it", async () => {
    const writeText = vi.fn(async () => undefined);
    vi.stubGlobal("navigator", { clipboard: { writeText } });
    const block = mount(SequenceBlock, {
      props: { kind: "cds", heading: "h", provenance: "p", sequence: "ATGAAA", what: "the gene" },
    });
    await block.find(".seq-copy").trigger("click");
    await block.vm.$nextTick();
    expect(writeText).toHaveBeenCalledWith("ATGAAA");
    expect(block.find(".seq-copy").text()).toBe("copied");
    vi.unstubAllGlobals();
  });

  it("⭐ offers the gene WITH its flanks, labelled as not the coding sequence alone", async () => {
    // It is the block a guide designer usually wants and the one most easily mistaken for the gene.
    // A copy button and its arithmetic (`app.js::seqGene`), not a fourth block of bases already on
    // screen above it.
    const card = mount(GeneSequenceCard, { props: { gene: gene(), flankLength: 100 } });
    const both = card.find(".seq-both");
    const copy = both.findComponent(CopySequenceButton);
    expect(copy.props("text")).toBe("G".repeat(100) + gene().coding_sequence + "C".repeat(100));
    expect(copy.props("what")).toContain("not the coding sequence alone");
    expect(copy.text()).toBe("copy gene + 100 bp each side");
    expect(both.find(".seq-both-note").text()).toContain(
      `${(200 + gene().coding_sequence.length).toLocaleString()} bases: 100 + ${gene().coding_sequence.length} + 100.`,
    );
    // ⛔ and the four blocks are the four things, not five.
    expect(card.findAllComponents(SequenceBlock)).toHaveLength(4);
  });

  it("⭐ numbers each block from 1 WITHIN itself, sixty to a line, and copies the bases alone", () => {
    const block = mount(SequenceBlock, {
      props: { kind: "cds", heading: "h", provenance: "p", sequence: "A".repeat(125), what: "w" },
    });
    // `textContent`, not `.text()`: the latter trims, and the leading pad IS the alignment.
    const lines = (block.find(".seq-pre").element.textContent ?? "").split("\n");
    expect(lines).toHaveLength(3);
    expect(lines[0]).toBe(`  1  ${"A".repeat(60)}`);
    expect(lines[2]).toBe(`121  ${"A".repeat(5)}`);
    expect(block.findComponent(CopySequenceButton).props("text")).toBe("A".repeat(125));
  });

  it("⚠ an empty flank is a FACT about the assembly, and says so", () => {
    const block = mount(SequenceBlock, {
      props: { kind: "up", heading: "h", provenance: "", sequence: "", what: "w" },
    });
    expect(block.find(".seq-none").text()).toBe("None in the assembly.");
    expect(block.find(".seq-copy").exists()).toBe(false);
  });
});
