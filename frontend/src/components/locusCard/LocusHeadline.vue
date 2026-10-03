<script setup lang="ts">
/**
 * What this locus IS, and how good it is — the card the `--anchor` accent ties to the focal block
 * on the track, so it is visible at a glance which gene is being described.
 *
 * ⭐ **Two rows of three tiles: the first says WHAT, the second says HOW GOOD.** Shares, not
 * fractions — the denominator is the same on every locus, so printing it six times a page adds
 * nothing a reader has to hold. Member genes, named genes and the family count are all stated
 * better lower down and are deliberately absent here.
 *
 * ⛔ **Every tile is now shown in its own units, and none is a rank.** `cluster separation` used to
 * hold the sixth slot as a percentile; it moved down into `EmbeddingSimilarityCard`, which already
 * draws it as a `.pair.sep` row against each representation's OWN measurable count — so nothing was
 * lost, and the slot went to the number David asked for. The two within-cluster tiles stay in their
 * own units because almost every locus in these catalogues is cohesive, so ranking cohesion against
 * its peers says the opposite of the truth. See `lib/similarityViews`.
 *
 * ⭐ **The sixth tile is SEQUENCE IDENTITY** (David, 2026-10-03): *"this is effectively a minimum for
 * the identity across all members of the set … but currently buried in text."* It was — the
 * `holds together at 60% identity` chip on the card below. It is `evidence.resolved_threshold`.
 *
 * ⛔ **NOT `seqid_coverage`.** That column is populated on 0 of 17,531 ecoli-nuna4 loci; tiling it
 * would print a dash on every locus in both catalogues.
 */
import { computed } from "vue";

import type { Locus, Representation } from "@/api/types";
import { sharePercent } from "@/lib/formatting";
import { copiesPerGenome } from "@/lib/locusStatistics";
import { prevalenceBandLabel, prevalenceBandShade } from "@/lib/prevalence";

const props = defineProps<{
  locus: Locus;
  /** How many genomes the whole collection has — the denominator for "of genomes". */
  collectionGenomeCount: number;
}>();

/**
 * ⛔ These two tiles used to be `cohesion()`: `(intra − null_mean) / (1 − null_mean)`, a medoid
 * similarity rescaled to "distance from random towards perfect". The label said cohesion and the
 * number was a rescaling of a distance to ONE member, so the label and the number have now moved
 * together — this is the within-cluster median over every gene pair, shown directly and in the units
 * the card below shows. No rescaling, and nothing to misread as a percentile. The reasoning
 * `cohesion` carried survives where it belongs: on the SEPARATION tile, where the variation lives.
 */
function withinSimilarity(representation: Representation): number | null {
  return props.locus.similarity[representation]?.within_similarity ?? null;
}

/**
 * ⛔ **`null` here is a MEASUREMENT, and it is only ever one thing.** `within_similarity` is null for
 * exactly the singletons — 5,427 of 17,531 ecoli-nuna4 loci and 4,398 of 15,670 kp-nuna4 loci — and
 * for **no** multi-member locus (measured against both catalogues, 2026-10-03). A `0.0` would claim
 * their members are unrelated and a dash would say the question was asked and missed; `one member`
 * says the question does not arise. `schemas/locus_browser.schema.md:244` makes the same point.
 *
 * A dash is still right where a whole representation is absent from a multi-member locus, which is
 * not measured rather than not applicable.
 */
function withinTile(representation: Representation): { value: string; valueClass: string | null } {
  const value = withinSimilarity(representation);
  if (value !== null) return { value: value.toFixed(3), valueClass: null };
  if (props.locus.gene_count === 1) return { value: "one member", valueClass: "v-note" };
  return { value: "—", valueClass: null };
}

/**
 * ⚠ **Two different absences, and the tooltip must say which.** Measured on ecoli-nuna4: 5,896 loci
 * carry no `resolved_threshold`, of which **5,427 are singletons** (one gene — there is no pair to
 * align) and **469 are multi-member loci MMseqs never grouped at any rung of the ladder**, held
 * together by other evidence instead — `esm_homology` 290, `alignable_not_grouped` 113,
 * `pfam_not_alignable` 56, `synteny_only` 10. Printing one dash for both would merge "does not
 * arise" with "alignment found nothing", which are opposite findings.
 *
 * ⛔ **The 50-member cap is stated, never implied.** `resolved_threshold` is computed over at most
 * `FAMILY_MEMBER_CAP = 50` members, and 4,092 of 11,635 ecoli-nuna4 loci that have one sit at that
 * cap — 35 %. A threshold read as if it covered all 100 members would be wrong more than a third of
 * the time, so the flag the serialiser already ships goes in the title.
 */
const identityTitle = computed(() => {
  const evidence = props.locus.evidence;
  if (evidence.resolved_threshold === null) {
    return props.locus.gene_count === 1
      ? "One member gene, so there is no pair to align — an identity threshold does not arise here."
      : "Sequence alignment never grouped these members at any threshold on the ladder" +
          `${evidence.collapse_tier ? `; they are held together by ${evidence.collapse_tier}` : ""}.`;
  }
  const capped = evidence.resolved_threshold_is_capped_at_50_members;
  return (
    `The members group into one family once the identity threshold falls to ` +
    `${sharePercent(evidence.resolved_threshold)} — effectively a floor for identity across the set.` +
    (capped
      ? ` ⚠ Measured over the first 50 member genes, not all ${props.locus.gene_count.toLocaleString()}.`
      : "")
  );
});

/** ⛔ `—` where a number was never measured. Never `0`, which is a measurement. */
function orDash(value: number | null, format: (value: number) => string): string {
  return value === null ? "—" : format(value);
}

const copies = computed(() => copiesPerGenome(props.locus.gene_count, props.locus.genome_count));

/** ⚠ Named, because a bare array literal infers a UNION of six shapes and the template can then only
 *  reach the keys every member happens to carry. */
type Tile = {
  key: string;
  value: string;
  caption: string;
  valueClass?: string | null;
  title?: string;
};

const tiles = computed<Tile[]>(() => [
  {
    key: "prevalence",
    value: orDash(
      props.collectionGenomeCount > 0 ? props.locus.genome_count / props.collectionGenomeCount : null,
      (share) => sharePercent(share),
    ),
    caption: "of genomes",
  },
  { key: "copies", value: orDash(copies.value, (rho) => rho.toFixed(2)), caption: "copies per genome" },
  {
    key: "identity",
    value: orDash(props.locus.evidence.resolved_threshold, (threshold) => sharePercent(threshold)),
    caption: "sequence identity",
    title: identityTitle.value,
  },
  {
    key: "a5",
    value: orDash(props.locus.evidence.syntenic_a5, (a5) => a5.toFixed(2)),
    caption: "synteny A5",
  },
  { key: "bacformer-within", ...withinTile("bacformer"), caption: "within cluster · Bacformer" },
  { key: "esm-within", ...withinTile("esm"), caption: "within cluster · ESM" },
]);

/**
 * ⚠ The caveat sits BESIDE the name, not in the typography. A name rendered in a different weight
 * to say "inferred" just looks broken; a sentence says what it means.
 *
 * ⭐ **Every locus now states what its name RESTS ON, because the silence was the complaint**
 * (David, 2026-10-03): *"some nodes like gumC are largely being inferred. Only 2 / 100 gumC nodes
 * have Bakta names. Yet this is very quiet on the current page. Hard to see."* `gumC` is kp-nuna4
 * locus 722 — 100 member genes, `named_gene_count` 2 (`gumC` ×1, `wzc` ×1, 98 unnamed) — and this
 * line said nothing at all, because `bakta_symbol` returned `null`. It is the same argument
 * `locus_serialiser` already makes per UniRef50 family — *"the symbol column's DENOMINATOR, and it
 * travels with the modal symbol always"* — one level up, on the name the page leads with.
 *
 * ⛔ **Measured over all four catalogues, not assumed (2026-10-03).** `label`, `product` and
 * `pfam_architecture` have `named_gene_count = 0` in EVERY case — 8,741 / 175 / 370 loci in
 * ecoli-nuna4, with no exception in ecoli-nuna5, kp-nuna4 or kp-nuna5 either — so those three can
 * state flatly that no member gene carries a name. `bakta_symbol` is never 0: 1,113 of 8,245 are
 * partly named and 7,132 are fully named. Both halves get a line. "No note" is too weak a signal to
 * carry the 86 % that are fully named, and it is what left the 8,741 `label` loci — the largest
 * group in the catalogue — saying nothing.
 */
const nameProvenance = computed(() => {
  const named = props.locus.named_gene_count;
  const total = props.locus.gene_count;
  switch (props.locus.display_name_source) {
    case "bakta_symbol":
      if (total === 0) return null;
      return named < total
        ? `Name inferred from ${named.toLocaleString()} of ${total.toLocaleString()} named ` +
            `${named === 1 ? "gene" : "genes"} — the other ${(total - named).toLocaleString()} carry no name`
        : `Name carried by all ${total.toLocaleString()} member genes`;
    case "product":
      return "No member gene carries a name — inferred from the Bakta product";
    case "pfam_architecture":
      // The accession the name was read out of — `soleArch`'s single architecture.
      return `No member gene carries a name — inferred from ${props.locus.display_name_source_accession ?? "Pfam"}`;
    case "label":
      return "No member gene carries a name — this is the locus id";
    default:
      return null;
  }
});
</script>

<template>
  <div class="card focus-card">
    <div class="locus-head">
      <h1 class="locus-name">{{ locus.display_name }}</h1>
      <span
        class="band"
        :style="{ '--b': prevalenceBandShade(locus.prevalence_band).toFixed(2) }"
      >{{ prevalenceBandLabel(locus.prevalence_band) }}</span>
      <span class="locus-id">locus {{ locus.label }}</span>
    </div>

    <p v-if="nameProvenance" class="inferred-note">{{ nameProvenance }}</p>
    <p v-if="locus.best_product" class="lede">{{ locus.best_product }}</p>

    <div class="tiles">
      <div v-for="tile in tiles" :key="tile.key" class="tile" :title="tile.title">
        <div class="v" :class="tile.valueClass">{{ tile.value }}</div>
        <div class="k">{{ tile.caption }}</div>
      </div>
    </div>
  </div>
</template>
