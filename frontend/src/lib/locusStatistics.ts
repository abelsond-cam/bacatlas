/**
 * The small numbers on the locus card that are *computed* rather than read — and the rules for what
 * each one may claim.
 *
 * ⛔ **Almost everything that used to live here has gone, and each departure is a fact about the
 * data rather than a tidy-up.** `similarityFromDistance` converted a stored cosine DISTANCE into a
 * similarity; the API now serves similarities everywhere it once served distances — on the locus
 * card and on the footer's residual rows alike — so a surviving `1 − d` would silently print the
 * complement of a real number. `cohesion` rescaled a member's distance to ONE gene (its medoid)
 * against the random-pair baseline and captioned the result "cohesion"; the label and the number
 * have now moved together, and the within-cluster median over every gene pair is shown directly, in
 * the same units the card below it shows. `separationVerdict` banded exactly one of what are now
 * three views, so it lives with the view table in `lib/similarityViews`.
 *
 * ⛔ Three numbers were re-derived on the published page and must not be re-derived here either:
 * every percentile is a precomputed midrank over measurable loci, the Pfam facts are the audit's own
 * (`pfam.dominant_domain_gene_count` and its kin — since 2026-11-26 numbers only, the clan verdict having
 * left the wire), and the display name and best product are ingest columns.
 */

/**
 * `p12` — a percentile as the card writes it.
 *
 * ⚠ Floored at 1, never 0: `p0` reads as "no data" beside every other tile that uses `—` for
 * exactly that, and the bottom of a midrank is a real rank rather than an absence.
 */
export function percentileLabel(percentile: number): string {
  return `p${Math.max(1, Math.round(percentile * 100))}`;
}

/** A signed number the card prints with its sign, e.g. `+0.412` / `-0.008`. */
export function signedFixed(value: number, digits = 3): string {
  return `${value >= 0 ? "+" : ""}${value.toFixed(digits)}`;
}

/** Copies per genome — `null` where the locus has no genomes to divide by. */
export function copiesPerGenome(geneCount: number, genomeCount: number): number | null {
  return genomeCount > 0 ? geneCount / genomeCount : null;
}
