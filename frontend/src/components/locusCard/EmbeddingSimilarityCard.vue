<script setup lang="ts">
/**
 * **Embedding Similarity — consistency, never corroboration.**
 *
 * ⛔ Embedding distance is *the clustering's own objective*, so this card is reported as consistency
 * and diagnosis and must never be read as independent support. The independent evidence is above:
 * UniRef50, Pfam and genomic context. The published card says so in its own footer and that sentence
 * is carried here verbatim.
 *
 * ⛔ **This replaced the medoid `Embedding geometry` card, and the measurement changed with it.**
 * The old card reduced a locus to ONE member — the gene nearest its centroid — and measured that
 * point. Every number here is a median over whole gene SETS, or is anchored on the single member
 * least attached to its locus. The old construction was biased optimistic: the all-pairs median sits
 * 0.009–0.017 below the medoid one, and separation < 0 was understated at 2.5 % against a true 4.2 %.
 *
 * ⭐ **Three views, and they REPLACE each other rather than stacking** (David, 2026-09-24). The
 * choice is the reader's and survives a walk, so it lives in `similarityViewStore` rather than being
 * rebuilt from the locus — and never in the URL.
 *
 * ⚠ **Both representations are STACKED inside whichever view is showing — they are not tabs.** The
 * two disagree about WHICH loci are weak (of the 1,059 ecoli loci Bacformer flags, ESM rescues 817),
 * and a tab would hide exactly that disagreement behind a click.
 *
 * ⚠ **Not measurable is a sentence, not a zero.** A singleton has no pair inside its locus at all —
 * 5,427 of *E. coli*'s 17,531 loci — and rendering that as `0.000` puts a made-up number beside real
 * ones.
 */
import { computed } from "vue";

import type {
  LocusSimilarity,
  NeighbourDisplayRow,
  Representation,
  SimilarityBaseline,
} from "@/api/types";
import { percentileLabel, signedFixed } from "@/lib/locusStatistics";
import {
  cardEstimator,
  formatRowValue,
  isMeasurable,
  SIMILARITY_VIEWS,
  viewById,
  viewNote,
  viewTab,
  viewValue,
  viewVerdict,
  type FloorTick,
  type SimilarityViewId,
} from "@/lib/similarityViews";

import FloorStrip from "./FloorStrip.vue";
import NearestLociCard from "./NearestLociCard.vue";

const props = defineProps<{
  similarity: Readonly<Record<Representation, LocusSimilarity | null>>;
  /** The random gene-pair floor per representation — without it a cosine on this card has no scale. */
  baselines: readonly SimilarityBaseline[];
  /** ⛔ The fan-out rows the nearest loci resolve their catalogue ordinals through. */
  neighbours: readonly NeighbourDisplayRow[];
  /** The locus's member-gene count: how many members point OUT, and why a single gene has no pair. */
  geneCount: number;
  view: SimilarityViewId;
}>();

const emit = defineEmits<{
  walk: [locusLabel: string];
  selectView: [view: SimilarityViewId];
}>();

/** ⚠ Bacformer first: it is the context axis the track is built on. */
const ORDER: readonly { readonly representation: Representation; readonly heading: string }[] = [
  { representation: "bacformer", heading: "Bacformer" },
  { representation: "esm", heading: "ESM" },
];

const activeView = computed(() => viewById(props.view));
/** ⛔ Median or mean, from the catalogue's own baselines — the card never assumes which. */
const estimator = computed(() => cardEstimator(props.baselines));

const baselineFor = computed(() => {
  const index = new Map<string, SimilarityBaseline>();
  for (const baseline of props.baselines) index.set(baseline.representation, baseline);
  return index;
});

/**
 * ⛔ Whether the card exists at all, tested across EVERY row key in both representations rather than
 * on the view in hand: a locus measurable only in the weakest-member view would otherwise lose its
 * card the moment the reader opened on the median.
 */
const hasAnySimilarity = computed(() =>
  ORDER.some(({ representation }) => {
    const similarity = props.similarity[representation];
    if (similarity === null) return false;
    return SIMILARITY_VIEWS.some((view) => view.rows.some((row) => similarity[row.key] !== null));
  }),
);

const blocks = computed(() => {
  const view = activeView.value;
  return ORDER.flatMap(({ representation, heading }) => {
    const similarity = props.similarity[representation] ?? null;
    // ⛔ Tested on the FIRST row, never on `nearest_similarity`: that question is well posed for a
    // single gene and is populated for a singleton, so testing it would give a singleton a block
    // with one rail and no partner — the exact thing the empty sentence below exists to replace.
    if (!isMeasurable(similarity, view)) return [];

    const rows = view.rows.flatMap((row, index) => {
      const value = similarity[row.key];
      if (value === null) return [];
      return [{ key: row.key, label: row.label, value, text: formatRowValue(value, row.key), index }];
    });
    // ⛔ The tick class is the VIEW-row index, not a position in the surviving rows: with a missing
    // first row the rival would otherwise be drawn as the subject.
    const ticks: FloorTick[] = rows.map((row) => ({
      label: row.label,
      value: row.value,
      isSecondRow: row.index > 0,
    }));

    const verdict = viewVerdict(similarity, view);
    const value = viewValue(similarity, view);
    const baseline = baselineFor.value.get(representation) ?? null;

    return [
      {
        representation,
        heading,
        rows,
        ticks,
        verdict,
        baseline,
        nearestLoci: similarity.nearest_loci,
        /**
         * ⛔ **A share of exactly 1.0 gets NO percentile row.** 88.5 % of loci sit there, so the
         * midrank inside that tie block prints "p53" — reading as *better than half the catalogue*
         * when it means *tied with nearly all of it*. The share is already on the row above; below
         * 1.0 the whole tie block is above this locus and the rank means what it looks like.
         *
         * ⚠ Gated on the VALUE and not on the verdict, so a 1.0 whose midrank is missing is silent
         * for the same reason rather than falling through to "not measurable".
         */
        showsRankRow: view.difference !== null || value === null || value < 1,
      },
    ];
  });
});

/** ⚠ The denominator of "p12 of 12,104 loci", per representation — ESM's need not equal Bacformer's. */
function rankText(block: (typeof blocks.value)[number]): string {
  if (block.verdict === null) return "not measurable";
  const measurable = block.baseline?.measurable_locus_count ?? null;
  const over = measurable === null ? "the measurable loci" : `${measurable.toLocaleString()} loci`;
  return `${percentileLabel(block.verdict.percentile)} of ${over}`;
}

/**
 * The difference where the view has one; where it does not, the members pointing OUT of the locus.
 *
 * ⛔ Not the share again — the row above already prints it. ⚠ Rounded from the share and the member
 * count rather than served: at the largest locus (143 members) 4 dp resolves a single gene.
 */
function rankValue(block: (typeof blocks.value)[number]): string {
  if (block.verdict === null) return "—";
  if (activeView.value.difference !== null) return signedFixed(block.verdict.value);
  return `${Math.round((1 - block.verdict.value) * props.geneCount)} out`;
}

/**
 * ⚠ **5,427 of E. coli's 17,531 loci are single genes.** Without this the card renders its view
 * strip, its caption and nothing between them — which reads as a bug rather than as an answer.
 */
const emptyReason = computed(() =>
  props.geneCount === 1
    ? "A single gene has no pair inside its locus and no neighbour of its own, so none of these is " +
      "measurable for it."
    : "Not measurable for this locus.",
);

const note = computed(
  () =>
    `${viewNote(activeView.value, estimator.value)} Embedding distance is the clustering's own objective, so this is ` +
    "reported as consistency and diagnosis — never as corroboration. The independent evidence is " +
    "above: UniRef50, Pfam and genomic context.",
);

/** A similarity as a rail width, clamped into the same 0..1 axis the floor strip uses. */
function railPercent(value: number): string {
  return `${(100 * Math.max(0, Math.min(1, value))).toFixed(1)}%`;
}
</script>

<template>
  <div v-if="hasAnySimilarity" class="card sim-card">
    <h2>Embedding Similarity</h2>

    <!-- ⛔ Its OWN class, never a strip's class reused. Two button strips on one page must not share
         a selector, or a stylesheet rule — or a test — can mean one and hit the other. -->
    <div class="sim-views">
      <button
        v-for="option in SIMILARITY_VIEWS"
        :key="option.id"
        type="button"
        class="sim-view"
        :class="{ on: option.id === view }"
        :title="viewNote(option, estimator)"
        @click="emit('selectView', option.id)"
      >{{ viewTab(option, estimator) }}</button>
    </div>

    <div class="sim-body">
      <template v-for="block in blocks" :key="block.representation">
        <h3 class="sub-head">{{ block.heading }}</h3>

        <div v-for="row in block.rows" :key="row.key" class="pair" :class="{ inter: row.index > 0 }">
          <div class="lab">{{ row.label }}</div>
          <!-- inline width: jsdom computes no layout, so a stylesheet version is untestable -->
          <div class="rail"><i :style="{ width: railPercent(row.value) }" /></div>
          <div class="val">{{ row.text }}</div>
        </div>

        <!-- ⛔ Only a view whose rows are COSINES gets the axis. The own fraction is a share of
             members, and drawing it against a random gene-pair box would invite reading one as the
             other — the exact confusion the floor strip's own comment exists to prevent. -->
        <FloorStrip
          v-if="activeView.difference !== null && block.baseline"
          :floor-median="block.baseline.floor_median"
          :floor-p25="block.baseline.floor_p25"
          :floor-p75="block.baseline.floor_p75"
          :floor-p99="block.baseline.floor_p99"
          :ticks="block.ticks"
        />

        <div v-if="block.showsRankRow" class="pair sep" :title="block.verdict?.label">
          <div class="lab">{{ activeView.difference ?? "rank" }}</div>
          <div class="sep-pct">{{ rankText(block) }}</div>
          <div class="val">{{ rankValue(block) }}</div>
        </div>

        <NearestLociCard
          v-if="activeView.id === 'median'"
          :nearest-loci="block.nearestLoci"
          :neighbours="neighbours"
          @walk="emit('walk', $event)"
        />
      </template>

      <p v-if="!blocks.length" class="muted">{{ emptyReason }}</p>
    </div>

    <!-- its own class as well as `muted`: the body can hold a muted line of its own, and a selector
         that could match either would read one for the other -->
    <p class="muted sim-note">{{ note }}</p>
  </div>
</template>
