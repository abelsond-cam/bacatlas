import { describe, expect, it } from "vitest";

import type { Locus, PfamFamilyReference } from "@/api/types";

import {
  collapseTierStatement,
  cosineLowerBound,
  dominantDomainsText,
  pfamAccessionsIn,
  pfamArchitectureFacts,
  pfamEntryUrl,
} from "./evidenceVocabulary";

function family(shortName: string): PfamFamilyReference {
  return {
    short_name: shortName,
    description: "",
    interpro_accession: "",
    interpro_name: "",
    clan_accession: "",
    clan_name: "",
  };
}

const REFERENCE = { PF14052: family("Wzi"), PF00126: family("HTH_1"), PF03466: family("LysR_substrate"), PF09999: family("") };

/** A locus with only the fields these sentences read — and DIFFERENT values for each, so none can stand in for another. */
function locus(
  evidence: Partial<Locus["evidence"]> = {},
  pfam: Partial<Locus["pfam"]> = {},
  geneCount = 100,
): Locus {
  return {
    gene_count: geneCount,
    pfam: {
      annotated_gene_count: 88,
      architecture_count: 3,
      dominant_domains: "PF00126,PF03466",
      dominant_domain_gene_count: 82,
      other_architecture_gene_count: 6,
      architecture_shared: false,
      ...pfam,
    },
    evidence: {
      syntenic_a5: 0.9,
      collapse_tier: "mmseq@0.98",
      collapse_bucket: null,
      resolved_threshold: 0.98,
      resolved_threshold_is_capped_at_50_members: true,
      mmseqs_group_share: 0.98,
      esm_member_similarity_p05: 0.9126,
      esm_member_similarity_median: 0.9671,
      ...evidence,
    },
  } as unknown as Locus;
}

/** ⛔ Every sentence the page must never say again (David, 2026-11-26). */
const RETIRED_WORDING = /no homology|over-merge|grouped on context alone|conflicting|contested|evidence against/i;

describe("⭐ one sentence per tier, built from the numbers (David, 2026-11-26)", () => {
  it("MMseqs: quotes the share AND the identity rung", () => {
    const statement = collapseTierStatement(locus(), REFERENCE)!;
    expect(statement.label).toBe("groups at 98% identity");
    expect(statement.sentence).toBe(
      "98% of sampled members (a uniform sample of up to 50) group at 98% identity by MMseqs.",
    );
  });

  it("⚠ the rung and the share are different numbers, and each lands in its own place", () => {
    const statement = collapseTierStatement(
      locus({ collapse_tier: "mmseq@0.4", resolved_threshold: 0.4, mmseqs_group_share: 1 }),
      REFERENCE,
    )!;
    expect(statement.sentence).toBe(
      "100% of sampled members (a uniform sample of up to 50) group at 40% identity by MMseqs.",
    );
  });

  it("⛔ an older catalogue's rung, with no share, states only the rung — true under either definition", () => {
    const statement = collapseTierStatement(locus({ mmseqs_group_share: null }), REFERENCE)!;
    expect(statement.sentence).toBe("MMseqs groups these members at 98% identity.");
  });

  it("alignable_not_grouped: aligns at the permissive floor, groups at no identity", () => {
    const statement = collapseTierStatement(
      locus({ collapse_tier: "alignable_not_grouped", resolved_threshold: null, mmseqs_group_share: null }),
      REFERENCE,
    )!;
    expect(statement.sentence).toBe(
      "Cannot be grouped by MMseqs at any sequence identity down to 30%, although at least 95% of sampled " +
        "members join at its permissive floor (E ≤ 1e-3, coverage ≥ 50%).",
    );
  });

  it("pfam_shared: the share of ALL genes, and the domain set named by the reference — kp-71k 4008, Wzi", () => {
    const statement = collapseTierStatement(
      locus(
        { collapse_tier: "pfam_shared", resolved_threshold: null, mmseqs_group_share: null },
        { dominant_domains: "PF14052", dominant_domain_gene_count: 69_825, annotated_gene_count: 69_884,
          other_architecture_gene_count: 3, architecture_shared: true },
        69_884,
      ),
      REFERENCE,
    )!;
    // ⛔ 69,825 / 69,884 = 0.99916 — "100%" would claim every gene, so it keeps a decimal.
    expect(statement.sentence).toBe(
      "Cannot be grouped by MMseqs at any sequence identity; 99.9% of its genes share one Pfam architecture: " +
        "PF14052 (Wzi).",
    );
    expect(statement.sentence).not.toMatch(RETIRED_WORDING);
  });

  it("esm_only: no MMseqs, no shared architecture, and the ESM bound with its median", () => {
    const statement = collapseTierStatement(
      locus({ collapse_tier: "esm_only", resolved_threshold: null, mmseqs_group_share: null }),
      REFERENCE,
    )!;
    expect(statement.label).toBe("grouped by neither MMseqs nor Pfam");
    expect(statement.sentence).toBe(
      "Cannot be grouped by MMseqs at any sequence identity, and no more than 95% of its genes share one " +
        "Pfam architecture. ESM: 95% of members have cosine ≥ 0.912 to the locus medoid (median 0.967).",
    );
  });

  it("⚠ esm_only with no ESM medoid SAYS so — never a cosine of 0", () => {
    const statement = collapseTierStatement(
      locus({ collapse_tier: "esm_only", esm_member_similarity_p05: null, esm_member_similarity_median: null }),
      REFERENCE,
    )!;
    expect(statement.sentence).toContain("ESM: no member → medoid cosine was measured for this locus.");
  });

  it("⛔ carries an UNKNOWN tier through neutrally, as its raw token, with no sentence", () => {
    for (const token of ["some_future_tier", "synteny_only", "no_homology", "pfam_not_alignable", "esm_homology"]) {
      expect(collapseTierStatement(locus({ collapse_tier: token }), REFERENCE)).toEqual({ label: token, sentence: "" });
    }
  });

  it("shows nothing where there is no tier (a singleton)", () => {
    expect(collapseTierStatement(locus({ collapse_tier: null }), REFERENCE)).toBeNull();
    expect(collapseTierStatement(locus({ collapse_tier: "" }), REFERENCE)).toBeNull();
  });

  it("⛔ no tier's sentence or label says anything the 2026-11-26 ruling retired", () => {
    for (const tier of ["mmseq@0.3", "alignable_not_grouped", "pfam_shared", "esm_only"]) {
      const statement = collapseTierStatement(locus({ collapse_tier: tier }), REFERENCE)!;
      expect(`${statement.label} ${statement.sentence}`).not.toMatch(RETIRED_WORDING);
    }
  });
});

describe("⭐ the Pfam block is numbers only", () => {
  it("dominant share · other architectures · no Pfam — each over ALL the genes", () => {
    expect(pfamArchitectureFacts(locus(), REFERENCE)).toBe(
      "Pfam: 82% of genes carry PF00126 (HTH_1) + PF03466 (LysR_substrate) · 6% other architectures · " +
        "12% no Pfam domain.",
    );
  });

  it("leaves out a part that is exactly zero", () => {
    const facts = pfamArchitectureFacts(
      locus({}, { annotated_gene_count: 100, dominant_domain_gene_count: 100, other_architecture_gene_count: 0 }),
      REFERENCE,
    );
    expect(facts).toBe("Pfam: 100% of genes carry PF00126 (HTH_1) + PF03466 (LysR_substrate).");
  });

  it("⛔ keeps a decimal where whole percent would cross the strict 95 % line — 0.9496 is not 95%", () => {
    const facts = pfamArchitectureFacts(
      locus({}, { annotated_gene_count: 10_000, dominant_domain_gene_count: 9_496, other_architecture_gene_count: 504 }, 10_000),
      REFERENCE,
    );
    expect(facts).toContain("94.96% of genes carry");
  });

  it("⛔ is null on an older catalogue (never measured) and where no gene carries a domain", () => {
    expect(pfamArchitectureFacts(locus({}, { dominant_domains: null, dominant_domain_gene_count: null }), REFERENCE)).toBeNull();
    expect(
      pfamArchitectureFacts(
        locus({}, { annotated_gene_count: 0, dominant_domains: "", dominant_domain_gene_count: 0, other_architecture_gene_count: 0 }),
        REFERENCE,
      ),
    ).toBeNull();
  });

  it("names a family where the reference can, and leaves an empty short name as the bare accession", () => {
    expect(dominantDomainsText("PF09999,PF14052", REFERENCE)).toBe("PF09999 + PF14052 (Wzi)");
    expect(dominantDomainsText("PF77777", REFERENCE)).toBe("PF77777");
  });
});

describe("⛔ an ESM bound is FLOORED, never rounded up past what the data reaches", () => {
  it("floors to 3 dp, and survives a binary representation just below the stored value", () => {
    expect(cosineLowerBound(0.9126)).toBe("0.912");
    expect(cosineLowerBound(0.912)).toBe("0.912");
    expect(cosineLowerBound(0.9999)).toBe("0.999");
    expect(cosineLowerBound(-0.1234)).toBe("-0.124");
  });
});

describe("⛔ a Pfam accession is cut at the dot before it is looked up", () => {
  it("strips the version, which is how the reference is keyed", () => {
    expect(pfamAccessionsIn("PF00126.29,PF03466")).toEqual(["PF00126", "PF03466"]);
  });

  it("handles one, none, and stray whitespace", () => {
    expect(pfamAccessionsIn("PF01381")).toEqual(["PF01381"]);
    expect(pfamAccessionsIn(" PF01381 , PF00126 ")).toEqual(["PF01381", "PF00126"]);
    expect(pfamAccessionsIn(null)).toEqual([]);
    expect(pfamAccessionsIn("")).toEqual([]);
  });
});

describe("⭐ a Pfam chip links to InterPro where there is one", () => {
  it("prefers the integrated record", () => {
    expect(pfamEntryUrl("PF00126", "IPR000847")).toBe(
      "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000847/",
    );
  });

  it("⚠ falls back to the Pfam entry where the family has no integrated entry", () => {
    // "" is an ABSENCE here, and building a URL out of it would 404 on every clanless family.
    expect(pfamEntryUrl("PF01327", "")).toBe(
      "https://www.ebi.ac.uk/interpro/entry/pfam/PF01327/",
    );
  });
});
