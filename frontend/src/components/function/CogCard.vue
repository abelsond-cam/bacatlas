<script setup lang="ts">
/**
 * **COG** — the orthologous groups these genes were assigned, and how many of them there are.
 *
 * ⭐ **Partial COG coverage is HEADROOM, not a caveat (David, 2026-08-20).** COG is assigned by best
 * hit rather than by a profile HMM, so a gene without one is *unlabelled*, not *different*. That is
 * far weaker evidence of absence than a Pfam miss, which really is a profile failing to match — and
 * the card says so, because the alternative reads as a hole in the locus when it is a hole in the
 * annotation.
 *
 * ⛔ **This card still states only what the genes SAY — it does not predict.** Propagating a COG to
 * the unannotated members is *annotation transfer*, and it now exists, but as its own card
 * (`InferredFunctionCard.vue`) with the four things it needed before it could: a transfer rule (the
 * similarity ladder David chose, 2026-10-03), a confidence measure (the agreement rate measured on
 * the loaded catalogue), a decision about discordant and unverifiable loci (a call resting on one
 * gene is marked as unchecked, and the measured rate at which it still covers the node is stated
 * above the cards), and a marking that cannot be mistaken for a Bakta call. ⚠ **Keep
 * that separation.** Nothing inferred belongs on this card; the two sit side by side so a reader can
 * see which is which, which is the whole reason the transfer was held back until it could be.
 *
 * ⛔ **The distinct-id count is a COUNT, never a relation.** More than one orthologous group in a
 * locus is an ordinary consequence of grouping above the family level, which is what this method
 * does — it is read beside the sequence and context evidence, not as a fault.
 */
import { computed } from "vue";

import type { AnnotationEntry, FunctionResponse, PropagationRate } from "@/api/types";
import { pluralise } from "@/lib/formatting";
import {
  cogCategoryNames,
  cogEntryUrl,
  coverageParts,
  propagationRateSentence,
} from "@/lib/functionVocabulary";

import CountTable from "../shared/CountTable.vue";

const props = defineProps<{
  coverage: FunctionResponse["coverage"];
  entries: readonly AnnotationEntry[];
  /** The band base rate this node's COG propagates at — `null` where the node states none. */
  propagation: PropagationRate | null;
}>();

const annotated = computed(() => props.coverage.cog_annotated_gene_count);
const geneCount = computed(() => props.coverage.gene_count);

const coverage = computed(() =>
  coverageParts(annotated.value, geneCount.value, "a COG assignment"),
);

/**
 * ⭐ **What the COG on this node actually covers.** Present only where there is something to cover.
 *
 * ⛔ This replaced *"COG is assigned by best hit, not by a profile, so the remaining 99 are
 * unlabelled rather than different — the gap this locus could fill, not a gap in it"* (David,
 * 2026-10-05: *"Is awfully worded. It can just be COG is assigned to the whole node from the hit
 * within it with measured accuracy > 99.5 %. That is it."*). The old sentence was about the
 * annotator's method; the reader's question is what this locus is, and the answer is that the call
 * covers every member at a rate we have measured.
 */
const propagates = computed(() => {
  if (annotated.value === 0 || geneCount.value - annotated.value <= 0) return null;
  if (props.propagation === null) return null;
  return (
    "COG is assigned to the whole node from the hit within it — " +
    `${propagationRateSentence(props.propagation)}.`
  );
});

const categoryChip = computed(() => {
  const categories = props.coverage.modal_cog_categories;
  if (categories === null || categories.length === 0) return null;
  // ⚠ A category is one letter OR SEVERAL (`EP`, `KT`, `NUW`, `EHJQ` are all real), so each letter
  // is named. "Category W" is unreadable; "W — Extracellular structures" is the point of showing it.
  return { code: categories.join(""), names: cogCategoryNames(categories).join(" · ") };
});

/**
 * ⛔ **`null` where fewer than TWO genes carry a COG — one voter is not a vote.**
 *
 * `cog_distinct_id_count` is `nunique` over the annotated members, so on a node where one gene
 * carries a COG it is 1 by construction and the card was printing a `win`-toned *"one orthologous
 * group"* — unanimity claimed from a single voter, on **518 ecoli / 498 kp** loci. The rest of the
 * codebase already refuses that claim (`function_inference_service`'s `checkable`); this card did
 * not. David, 2026-10-05: *"The one orthologous group is silly. There is one gene."*
 */
const groupChip = computed(() => {
  const distinct = props.coverage.cog_distinct_id_count;
  if (annotated.value < 2) return null;
  return {
    tone: distinct <= 1 ? "win" : "neutral",
    label: distinct <= 1 ? "one orthologous group" : pluralise(distinct, "orthologous group"),
    isPlural: distinct > 1,
  };
});

const rows = computed(() =>
  props.entries.map((entry) => ({
    key: entry.term,
    count: entry.gene_count,
    term: entry.term,
    name: entry.name,
  })),
);
</script>

<template>
  <div class="card">
    <h3 class="sub-head">COG</h3>
    <!-- ⛔ Coverage first, always, and against the LOCUS size — a share against the annotated subset
         would read 100 % where one gene in forty carries a label. -->
    <p class="muted cover"><b v-if="coverage.emphasis">{{ coverage.emphasis }}</b>{{ coverage.rest }}</p>
    <p v-if="propagates" class="muted cover">{{ propagates }}</p>

    <template v-if="annotated > 0">
      <div class="chip-row">
        <!-- ⭐ The category is what this locus IS, and it was styled `neutral` — grey, the
             quietest tone the sheet has. David, 2026-10-05: "This is the most important thing and
             yet it is grey." -->
        <span
          v-if="categoryChip"
          class="chip lead"
          :title="`COG functional category ${categoryChip.code}`"
        >{{ categoryChip.code }} — {{ categoryChip.names }}</span>
        <span v-if="groupChip" class="chip" :class="groupChip.tone">{{ groupChip.label }}</span>
      </div>

      <CountTable v-if="rows.length" :headings="['COG', 'what it is']" :rows="rows" :total="geneCount">
        <template #cells="{ row }">
          <td class="acc">
            <a class="acc-link" :href="cogEntryUrl(row.term)" target="_blank" rel="noopener">{{ row.term }}</a>
          </td>
          <!-- ⚠ An em dash, not an empty cell: the vocabulary having no name is a fact about the
               vocabulary, and a blank reads as a rendering failure. -->
          <td>{{ row.name ?? "—" }}</td>
        </template>
      </CountTable>

      <p v-if="groupChip?.isPlural" class="muted">
        More than one orthologous group is an ordinary consequence of grouping above the family
        level, which is what this method does — read it beside the sequence and context evidence
        rather than as a fault.
      </p>
    </template>
  </div>
</template>
