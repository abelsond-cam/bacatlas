/**
 * The page's shared number formatting. One rule per number, in one place, so two panels showing
 * the same quantity cannot round it differently.
 */

/**
 * A share of a whole, as a percentage.
 *
 * ⚠ **A deliberate, single divergence from the published page**, written down because an
 * undocumented one is how two "authoritative" renderings come to disagree. `app.js`'s `pct` rounds
 * to whole percent everywhere it is called (`app.js:858`), so an occupant carried by 2 of 500
 * member genes prints as **"2 · 0%"** — a real, listed row reporting zero. Below 10 % this keeps
 * one decimal, so that row reads "2 · 0.4%" instead. Nothing above 10 % changes, which is every
 * number a parity diff of the shipped catalogues actually compares.
 */
export function sharePercent(share: number): string {
  // ⛔ The mirror case, at the top (2026-11-14): at 71,435 genomes a core locus in 71,136 of them
  // rounded to "100% of genomes" beside a count that says 299 do not carry it. A share short of the
  // whole keeps one decimal where whole-percent rounding would print 100, and never reads "100.0".
  if (share < 1 && Math.round(100 * share) >= 100) {
    const tenths = (100 * share).toFixed(1);
    return `${tenths === "100.0" ? "99.9" : tenths}%`;
  }
  return `${(100 * share).toFixed(share < 0.1 ? 1 : 0)}%`;
}

/**
 * A share as a percentage that never lands on the wrong side of what the page says about it.
 *
 * ⛔ **Whole percent, unless rounding would make the number contradict its own sentence.** The audit's
 * Pfam test is a strict `> 95 %`: 0.9496 printed as "95%" beside "no more than 95%" reads as a tie, and
 * 0.9504 printed as "95%" beside "more than 95%" reads as false. Decimals are added — up to four — until
 * the printed value sits on the same side of `threshold` as the true one. The same guard keeps a share short
 * of the whole from printing "100%" (69,825 of 69,884 genes → "99.9%") and a share above nothing from
 * printing "0%".
 */
export function percentKeepingSide(share: number, threshold?: number): string {
  if (share <= 0 || share >= 1) return `${Math.round(share * 100)}%`;
  const side = (value: number): number =>
    threshold === undefined ? 0 : value > threshold ? 1 : value < threshold ? -1 : 0;
  for (let digits = 0; digits < 4; digits += 1) {
    const text = (100 * share).toFixed(digits);
    const printed = Number.parseFloat(text) / 100;
    if (printed > 0 && printed < 1 && side(printed) === side(share)) return `${text}%`;
  }
  return `${(100 * share).toFixed(4)}%`;
}

/** `n` with the singular or plural of `noun`, e.g. `1 neighbourhood` / `3 neighbourhoods`. */
/**
 * `"a"` or `"an"`, for a noun the caller does not control.
 *
 * ⚠ **The vocabulary labels are the reason this exists**: three of the four take "a" (COG, GO,
 * KEGG) and the fourth takes "an" (EC number), so a hardcoded article renders *"carries a EC
 * number"* on exactly one tab in exactly one place. The letter test is enough for these four —
 * it is not a general English rule, and a noun whose spelling and sound disagree ("a UniRef50",
 * "an hour") needs the article passed in rather than guessed.
 */
export function indefiniteArticle(noun: string): string {
  return /^[aeiou]/i.test(noun) ? "an" : "a";
}

export function pluralise(count: number, noun: string, plural = `${noun}s`): string {
  return `${count} ${count === 1 ? noun : plural}`;
}
