/**
 * What the locus card says about how a locus holds together — **the audit's numbers, never its verdicts.**
 *
 * ⭐ **David's ruling of 2026-11-26.** At 71,435 genomes the audit's per-locus verdicts were decided by
 * stray genes: kp-71k locus 4008 (Wzi — PF14052 on 69,825 of 69,884 genes) was shown as *"grouped on
 * context alone — Unalignable, no Pfam, no ESM homology … the over-merge the audit exists to count"* and
 * *"conflicting architecture … evidence AGAINST this merge"*. Those were audit bugs. The page now never
 * says "no homology", "over-merge", "grouped on context alone", "conflicting" or "contested"; it states
 * three measured facts instead, one sentence per tier of the redesigned ladder:
 *
 * 1. **MMseqs** — a locus *groups at t % identity* when at least 95 % of its sampled members (a uniform
 *    sample of up to 50) fall in one connected component of MMseqs hits at identity ≥ t.
 * 2. **Pfam** — *one Pfam architecture is shared* when MORE THAN 95 % of ALL the locus's genes carry the
 *    same domain set (order and repeat count ignored). Under that, numbers only — no verdict word.
 * 3. **ESM** — when neither groups it, the 5th percentile of member → medoid ESM cosine, median beside it.
 *
 * The ladder, one tier per non-singleton locus, in order: `mmseq@{t}` → `alignable_not_grouped` →
 * `pfam_shared` → `esm_only`. ⛔ An UNKNOWN token — including the retired `synteny_only`, `no_homology`,
 * `pfam_not_alignable` and `esm_homology` an older catalogue still carries — renders as itself, neutrally,
 * with no sentence: the token is a real record, and no wording here describes it any more.
 *
 * ⛔ **SUPERSEDED (David, 2026-11-26).** This file used to carry a "SETTLED, do not re-raise (David,
 * 2026-09-13)" note keeping the published page's `synteny_only` wording — *"the over-merge the audit exists
 * to count"*. The 2026-11-26 ruling replaces it: that wording, the tier vocabulary it belonged to and the
 * Pfam class chips (`PFAM_VERDICTS`) are gone from the card.
 */
import type { Locus, PfamFamilyReference } from "@/api/types";
import { percentKeepingSide } from "@/lib/formatting";

/**
 * ⛔ **The thresholds the sentences state, as literals** — and the ingest refuses any catalogue whose audit
 * summary records other values (`catalogue_frames.PAGE_STATED_AUDIT_THRESHOLDS`), so these words cannot go
 * on stating a number the audit did not use.
 */
export const MMSEQS_GROUP_MIN_SHARE = 0.95;
export const PFAM_SHARED_MIN_SHARE = 0.95;
export const ESM_QUOTED_QUANTILE = 0.05;

/** How a locus holds together, as the card states it: a short neutral chip and one sentence. */
export interface TierStatement {
  /** The chip's words — a description, never a grade. An unknown tier's chip is its raw token. */
  readonly label: string;
  /** The one sentence, built from the numbers. Empty for an unknown token, which has no wording. */
  readonly sentence: string;
}

/** `mmseq@0.98` — the tier that names the identity at which the members group. */
const MMSEQS_TIER = /^mmseq@([0-9.]+)$/;

const NO_MMSEQS = "Cannot be grouped by MMseqs at any sequence identity";

/** A whole-percent identity rung: `0.98` → `98%`. Rungs are whole percents by construction. */
function rung(identity: number): string {
  return `${Math.round(identity * 100)}%`;
}

/**
 * A cosine quoted as a LOWER bound, floored rather than rounded to 3 dp.
 *
 * ⛔ "95 % of members have cosine ≥ X" is a bound: rounding 0.9126 up to 0.913 would state a bound the data
 * does not reach. The epsilon keeps a value stored as 0.912 (= 911.9999… × 10⁻³ in binary) from flooring
 * to 0.911.
 */
export function cosineLowerBound(value: number): string {
  return (Math.floor(value * 1000 + 1e-9) / 1000).toFixed(3);
}

/**
 * The dominant domain set, named where the Pfam reference can: `PF00126 (HTH_1) + PF03466 (LysR_substrate)`.
 *
 * ⚠ An empty `short_name` is what the vendored table HAD, not a missing key, so the accession stands alone.
 */
export function dominantDomainsText(
  domains: string,
  reference: Readonly<Record<string, PfamFamilyReference>>,
): string {
  return pfamAccessionsIn(domains)
    .map((accession) => {
      const name = reference[accession]?.short_name ?? "";
      return name === "" ? accession : `${accession} (${name})`;
    })
    .join(" + ");
}

/**
 * The dominant domain set's share of ALL the locus's genes — `null` where the audit never measured it.
 *
 * ⛔ Over `gene_count`, never over the annotated genes: a set on 40 of 40 annotated genes in a locus of 100
 * is on 40 % of the locus, and a share over the annotated subset would print 100 %.
 */
export function dominantDomainShare(locus: Locus): number | null {
  const count = locus.pfam.dominant_domain_gene_count;
  return count === null || locus.gene_count === 0 ? null : count / locus.gene_count;
}

/**
 * The tier as the card states it, or `null` where there is no tier (a singleton).
 *
 * ⚠ `mmseq@…` is parametric, so it is parsed rather than looked up. On a catalogue audited before
 * 2026-11-26 the share is `null` and the rung's DEFINITION was different (the 50 longest members, not a
 * sample at 95 %), so the sentence then states only the rung, which is true under either definition.
 */
export function collapseTierStatement(
  locus: Locus,
  reference: Readonly<Record<string, PfamFamilyReference>>,
): TierStatement | null {
  const tier = locus.evidence.collapse_tier;
  if (tier === null || tier === "") return null;

  const mmseqs = MMSEQS_TIER.exec(tier);
  if (mmseqs !== null) {
    const identity = rung(Number.parseFloat(mmseqs[1] as string));
    const share = locus.evidence.mmseqs_group_share;
    return {
      label: `groups at ${identity} identity`,
      sentence:
        share === null
          ? `MMseqs groups these members at ${identity} identity.`
          : `${percentKeepingSide(share, MMSEQS_GROUP_MIN_SHARE)} of sampled members (a uniform sample of up ` +
            `to 50) group at ${identity} identity by MMseqs.`,
    };
  }

  switch (tier) {
    case "alignable_not_grouped":
      return {
        label: "aligns · groups at no identity",
        sentence:
          `${NO_MMSEQS} down to 30%, although at least ${rung(MMSEQS_GROUP_MIN_SHARE)} of sampled members ` +
          "join at its permissive floor (E ≤ 1e-3, coverage ≥ 50%).",
      };
    case "pfam_shared": {
      const share = dominantDomainShare(locus);
      const domains = locus.pfam.dominant_domains;
      // ⚠ After a colon, not in parentheses: a named family is already `PF14052 (Wzi)`, and nesting it in a
      // second pair reads as a typo.
      const named = domains ? `: ${dominantDomainsText(domains, reference)}` : "";
      return {
        label: "one Pfam architecture",
        sentence:
          share === null
            ? `${NO_MMSEQS}; more than ${rung(PFAM_SHARED_MIN_SHARE)} of its genes share one Pfam ` +
              `architecture${named}.`
            : `${NO_MMSEQS}; ${percentKeepingSide(share, PFAM_SHARED_MIN_SHARE)} of its genes share one ` +
              `Pfam architecture${named}.`,
      };
    }
    case "esm_only": {
      const p05 = locus.evidence.esm_member_similarity_p05;
      const median = locus.evidence.esm_member_similarity_median;
      // ⚠ "no more than", not "fewer than": `pfam_shared` is a STRICT `> 95 %`, so 19 of 20 genes (exactly
      // 95 %) is not shared — and "fewer than 95 %" would be false of it.
      const head =
        `${NO_MMSEQS}, and no more than ${rung(PFAM_SHARED_MIN_SHARE)} of its genes share one Pfam ` +
        "architecture.";
      const esm =
        p05 === null
          ? " ESM: no member → medoid cosine was measured for this locus."
          : ` ESM: ${rung(1 - ESM_QUOTED_QUANTILE)} of members have cosine ≥ ${cosineLowerBound(p05)} to ` +
            `the locus medoid${median === null ? "" : ` (median ${median.toFixed(3)})`}.`;
      return { label: "grouped by neither MMseqs nor Pfam", sentence: head + esm };
    }
    default:
      // ⛔ Carried through rather than dropped: the token is a real record. No sentence — none fits it.
      return { label: tier, sentence: "" };
  }
}

/**
 * The Pfam block's numbers — the dominant set's share, the other architectures', and no Pfam at all — each
 * over ALL the locus's genes. `null` where the audit never measured them (a catalogue audited before
 * 2026-11-26), and for a locus where no gene carries a domain, which the card says in words instead.
 *
 * ⛔ **No verdict word, whatever the share** (David, 2026-11-26): "82% carry PF00126 · 6% other
 * architectures · 12% no Pfam" is the whole statement. A part that is exactly zero is left out.
 */
export function pfamArchitectureFacts(
  locus: Locus,
  reference: Readonly<Record<string, PfamFamilyReference>>,
): string | null {
  const { dominant_domains: domains, dominant_domain_gene_count: dominant } = locus.pfam;
  const other = locus.pfam.other_architecture_gene_count;
  const annotated = locus.pfam.annotated_gene_count;
  const size = locus.gene_count;
  if (domains === null || dominant === null || other === null || annotated === 0 || size === 0) return null;
  const parts = [
    `${percentKeepingSide(dominant / size, PFAM_SHARED_MIN_SHARE)} of genes carry ${dominantDomainsText(domains, reference)}`,
  ];
  if (other > 0) parts.push(`${percentKeepingSide(other / size)} other architectures`);
  const unannotated = size - annotated;
  if (unannotated > 0) parts.push(`${percentKeepingSide(unannotated / size)} no Pfam domain`);
  return `Pfam: ${parts.join(" · ")}.`;
}

/**
 * `PF00126.29,PF03466` → `["PF00126", "PF03466"]`.
 *
 * ⛔ The version suffix is cut before lookup, because the reference is keyed version-stripped. An
 * accession looked up as written silently misses and falls back to a bare chip — which reads as
 * *"this family has no name"* rather than as a failed join.
 */
export function pfamAccessionsIn(architecture: string | null): string[] {
  if (architecture === null || architecture === "") return [];
  const out: string[] = [];
  for (const raw of architecture.split(",")) {
    const accession = (raw.split(".")[0] ?? "").trim();
    if (accession !== "") out.push(accession);
  }
  return out;
}

/**
 * Where a Pfam chip links.
 *
 * ⭐ The **InterPro** entry is preferred over the Pfam one: it is the integrated record, and the
 * page a reader following a domain actually wants. Families with no integrated entry fall back to
 * the Pfam entry — which is why `""` here is an absence to test for, never a URL fragment.
 */
export function pfamEntryUrl(accession: string, interproAccession: string): string {
  return interproAccession === ""
    ? `https://www.ebi.ac.uk/interpro/entry/pfam/${encodeURIComponent(accession)}/`
    : `https://www.ebi.ac.uk/interpro/entry/InterPro/${encodeURIComponent(interproAccession)}/`;
}
