/**
 * The two verdict vocabularies the locus card chips: the **collapse tier** (how the members hold
 * together by sequence) and the **Pfam concordance class** (what Pfam makes of them).
 *
 * ⛔ **Both are READ, never re-derived.** `app.js:3208` records that the page once re-derived the
 * Pfam verdict by a different rule — every architecture against the first, where the audit's
 * `worst_relation` compares all PAIRS and drops non-maximal sets — and the two disagreed on 2 of
 * 22,624 loci. *"A page that quotes the report must not be able to contradict it, ever."* What is
 * here is the **vocabulary**: what each class is called and how loudly to say it. The class itself
 * arrives in the response.
 *
 * ⚠ Only `disjoint` is the §6.2 conflict. `nested` is the signature of PARTIAL annotation, which is
 * absent evidence and must never be dressed as disagreement.
 */

/** How loudly a chip speaks. `neutral` is a real value — it is not "no verdict". */
export type VerdictTone = "win" | "neutral" | "warn" | "bad";

export interface Verdict {
  readonly tone: VerdictTone;
  /** The chip's own words. */
  readonly label: string;
  /** The sentence under it. Empty where the chip says everything there is to say. */
  readonly note: string;
}

/**
 * ⚠ **`synteny_only`'s wording is the published page's, and stays that way — SETTLED, do not
 * re-raise (David, 2026-09-13).**
 *
 * Its note calls the tier *"the over-merge the audit exists to count"*, which reads against
 * `CLAUDE.md`'s later position that *"`synteny only` names the evidence, NOT a mistake"* (David,
 * 2026-08-18). Asked directly, David ruled: **not relevant to the backend rebuild, not a porting
 * decision, and not a major discrepancy — leave it.** The rebuild's job here is to carry the page
 * across faithfully; re-wording the science it states is a separate piece of work on a separate
 * day. This note exists so the next reader who spots the tension finds the answer rather than
 * re-opening it.
 */
export const COLLAPSE_TIERS: Readonly<Record<string, Verdict>> = {
  pfam_not_alignable: {
    tone: "win",
    label: "never aligns · shared Pfam architecture",
    note:
      "MMseqs cannot align these members at all, yet they share a Pfam domain architecture — " +
      "homologous, not an over-merge.",
  },
  esm_homology: {
    tone: "win",
    label: "never aligns · ESM homology",
    note:
      "No alignment and no shared Pfam, but ESM (sequence-only, and blind to the context signal " +
      "that built this locus) finds them homologous — the short proteins alignment cannot bridge.",
  },
  synteny_only: {
    tone: "bad",
    label: "grouped on context alone",
    note:
      "Unalignable, no Pfam, no ESM homology: grouped on genomic context only. This is the " +
      "over-merge the audit exists to count.",
  },
  /**
   * ⚠ **RETIRED as a category** — `CLAUDE.md`: *"`no_homology` is retired — it meant not measured,
   * not nothing found."* Kept in the vocabulary because catalogues graded before the retirement
   * still carry the token, and a tier with no entry would render as a bare unstyled string.
   */
  no_homology: {
    tone: "bad",
    label: "no evidence at all",
    note: "Unalignable with no supporting evidence of any kind — the worst residual.",
  },
  alignable_not_grouped: {
    tone: "warn",
    label: "aligns, but never as one group",
    note:
      "Members align pairwise at the permissive floor but never resolve into a single group on " +
      "the identity ladder.",
  },
};

/** `mmseq@0.98` — the tier that names the identity at which the members first group. */
const MMSEQS_TIER = /^mmseq@([0-9.]+)$/;

/**
 * A collapse tier as a chip.
 *
 * ⚠ The `mmseq@…` tiers are **parametric**, not a fixed list, so they are parsed rather than looked
 * up — and the cut at 0.5 is where the sentence changes from *"ordinary sequence homology"* to a
 * threshold far below any pangenome tool's floor (Panaroo 70 %, Roary 95 %).
 *
 * ⛔ An unknown tier returns a `neutral` chip carrying the raw token rather than `null`: the token
 * is a real measurement, and dropping it would say the evidence was never looked at.
 */
export function collapseTierVerdict(tier: string | null): Verdict | null {
  if (tier === null || tier === "") return null;
  const known = COLLAPSE_TIERS[tier];
  if (known !== undefined) return known;

  const match = MMSEQS_TIER.exec(tier);
  if (match !== null) {
    const identity = Number.parseFloat(match[1] as string);
    const percent = Math.round(identity * 100);
    return identity <= 0.5
      ? {
          tone: "warn",
          label: `holds together at ${percent}% identity`,
          note:
            `Sequence alignment only joins these members once the identity threshold drops to ` +
            `${percent}% — far below any pangenome tool's floor (Panaroo 70%, Roary 95%).`,
        }
      : {
          tone: "neutral",
          label: `holds together at ${percent}% identity`,
          note:
            `Ordinary sequence homology: the members group on the identity ladder at ${percent}%.`,
        };
  }
  return { tone: "neutral", label: tier, note: "" };
}

/**
 * What each Pfam concordance class means, and how loudly to say it.
 *
 * ⛔ `no_coverage` is deliberately absent from this map. Fewer than one annotated gene is neither
 * agreement nor disagreement, and giving it a chip would put a verdict on screen where none was
 * reached — the card states coverage in words instead.
 */
export const PFAM_VERDICTS: Readonly<Record<string, Verdict>> = {
  single: {
    tone: "win",
    label: "one architecture",
    note:
      "Every annotated gene here carries the same domain architecture — nothing discordant to " +
      "weigh against the merge.",
  },
  same_domains: {
    tone: "win",
    label: "one architecture",
    note:
      "The annotated genes carry the same domains, differing only in repeat count or order — the " +
      "same protein, not a different one.",
  },
  nested: {
    tone: "neutral",
    label: "partial annotation",
    note:
      "The architectures are nested: some members carry a domain the others do not, which is a " +
      "weaker Pfam hit rather than a competing architecture. Absent evidence, not conflicting " +
      "evidence.",
  },
  overlapping: {
    tone: "warn",
    label: "architectures overlap",
    note:
      "The architectures share a domain but neither contains the other. Weigh this against the " +
      "sequence and context evidence.",
  },
  // ⛔ The reader-facing note cited `nuna_structural_syntology.md §6.2` until 2026-10-03. That file
  // exists and §6.2 ("Pfam domain architecture — structural evidence") is the argument behind this
  // chip — but it lives in the PRIVATE nuna repo, so to everyone reading BacAtlas the citation named
  // a document they cannot open. Provenance belongs in the source, where it can be followed; the note
  // keeps the claim and drops the dead reference.
  disjoint: {
    tone: "bad",
    label: "conflicting architecture",
    note:
      "The minor architecture shares NO domain family or clan with the dominant one — different " +
      "folds at the same locus. A discordant architecture is evidence AGAINST this merge, and this " +
      "locus should be treated as contested until the sequence and context evidence is weighed.",
  },
};

/**
 * The Pfam verdict chip, or `null` where none may be shown.
 *
 * ⛔ **Two separate reasons to show nothing, and both must hold before a chip appears**: the audit
 * reached no class, or no gene here carries a Pfam-A domain at all. A verdict issued without
 * coverage behind it is the exact failure `no_coverage` exists to prevent — see T7.
 */
export function pfamVerdict(
  concordanceClass: string | null,
  annotatedGeneCount: number | null,
): Verdict | null {
  if (concordanceClass === null || annotatedGeneCount === null || annotatedGeneCount === 0) {
    return null;
  }
  return PFAM_VERDICTS[concordanceClass] ?? null;
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
