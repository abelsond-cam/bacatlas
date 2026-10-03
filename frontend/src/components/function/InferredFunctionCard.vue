<script setup lang="ts">
/**
 * **Inferred function** — what can be said about a node that Bakta never annotated, and how sure.
 *
 * ⛔ **Separate from every other card on this tab, and it must stay separate.** The cards above print
 * what this node's own genes SAY; this one prints what a neighbour SUGGESTS. A reader left to tell
 * those apart by context will read a suggestion as an annotation, which is the one failure this card
 * is designed around — so it carries its own heading, its own "suggested" wording, and never a bare
 * term where the other cards put a real one.
 *
 * ⭐ **Every rate shown here was measured on the loaded catalogue**, by
 * `instruments/annotation_transfer.py`, which is also the code
 * `nuna/docs/model_evaluation/COG_function_inference.md` reports from. Not a constant lifted out of
 * the write-up: that would keep reading plausibly long after the catalogue moved underneath it.
 *
 * ⛔ **Three states, and only one of them is a suggestion.**
 *  1. a donor was found and its similarity earns a depth → the suggestion, with its measured rate;
 *  2. a donor was found but sits below 0.90, or states nothing at the quotable depth → the
 *     neighbour is named and NOTHING is suggested ("too remote to call");
 *  3. no neighbour carries this vocabulary at all → the dark-set statement, which is a FINDING about
 *     where this node sits and not a failure to load.
 */
import { computed } from "vue";

import type { CalibrationLadder, VocabularyInference } from "@/api/types";

const props = defineProps<{
  inference: VocabularyInference;
  ladder: CalibrationLadder;
  /** How the vocabulary is named to a reader — "COG", not `cog_orthogroup`. */
  vocabularyLabel: string;
}>();

const emit = defineEmits<{ walk: [label: string] }>();

/**
 * ⚠ Three decimals, matching `NearestLociCard` and `LocusHeadline` — the SAME quantity is on the
 * similarity card, and two panels rounding one cosine differently is the thing `lib/formatting`
 * exists to prevent. Local rather than shared only because those two also inline it; moving all
 * three behind one helper is a separate change.
 */
function cosine(value: number | null): string {
  return value === null ? "—" : value.toFixed(3);
}

const candidate = computed(() => props.inference.candidate);

/** The ranks tried before the donor — the part a card showing only the donor would hide. */
const emptySteps = computed(() => {
  const found = candidate.value;
  return props.inference.walk.filter(
    (step) => !step.carries_annotation && (found === null || step.rank < found.rank),
  );
});

/**
 * ⛔ The one sentence this card exists to get right. It names the depth, the rate, the pair count
 * and the lift — in that order, because the rate is meaningless without the n and the n is
 * meaningless without knowing what was compared.
 */
const confidenceSentence = computed(() => {
  const found = candidate.value;
  if (found === null || found.level === null) return null;
  const cell = found.calibration;
  if (cell === null || cell.agreement === null) {
    return `Too few measured pairs at this similarity to state how often a transfer of this kind is right${
      cell === null ? "" : ` (only ${cell.pairs.toLocaleString()})`
    } — treat it as a neighbour, not a suggestion.`;
  }
  const percent = (cell.agreement * 100).toFixed(1);
  const lift = cell.lift === null ? null : Math.round(cell.lift).toLocaleString();
  return (
    `Nodes this similar agree on their ${cell.level_label} ${percent}% of the time ` +
    `(${cell.pairs.toLocaleString()} measured pairs in this catalogue)` +
    (lift === null ? "." : `, which is ${lift}× more often than a random annotated node.`)
  );
});

/** ⚠ The donor's OWN support, by the same rule the node's own call is judged by. */
const donorCaveat = computed(() => {
  const found = candidate.value;
  if (found === null || found.donor.checkable) return null;
  return (
    `⚠ The neighbour's own ${props.vocabularyLabel} rests on ` +
    `${found.donor.annotated_gene_count} of its ${found.donor.member_gene_count.toLocaleString()} ` +
    `genes, so it was never checked against itself either.`
  );
});

const measuredCells = computed(() => props.ladder.cells.filter((cell) => cell.agreement !== null));
</script>

<template>
  <div class="card infer">
    <!--
      ⛔ The heading must not repeat the card above it. The real COG card already says "None of these
      102 genes carries a COG assignment"; a second card headed "none assigned here" reads as the
      same finding twice, which is what the first render of this looked like. It says instead where
      the suggestion COMES FROM, which is the one thing the card above cannot say.
    -->
    <h3 class="sub-head">
      {{ vocabularyLabel }} —
      {{ candidate === null ? "nothing near enough to infer from" : "inferred from a neighbour" }}
      <span class="infer-tag">inferred</span>
    </h3>

    <!-- 3 · nothing within reach. A finding, stated as one. -->
    <p v-if="candidate === null" class="muted cover">
      None of the {{ inference.walk.length }} nearest nodes by ESM similarity carries a
      {{ vocabularyLabel }} either, so there is nothing to transfer. That is common and it is
      informative: unannotated nodes sit next to unannotated nodes far more often than chance, which
      is why widening the neighbour list is not the missing piece.
    </p>

    <template v-else>
      <!-- 1 · the suggestion -->
      <p v-if="candidate.level !== null" class="infer-value">
        <span class="infer-lead">Suggested {{ candidate.calibration?.level_label ?? "value" }}:</span>
        <strong>{{ candidate.value?.join(", ") }}</strong>
      </p>
      <!-- 2 · a donor too remote, or too shallow, to quote -->
      <p v-else class="muted cover">
        The nearest node carrying a {{ vocabularyLabel }} is at
        {{ cosine(candidate.cosine) }} similarity — too remote to call, so nothing is suggested.
      </p>

      <p class="infer-from">
        from
        <button type="button" class="infer-walk" @click="emit('walk', candidate.donor.node_label)">
          {{ candidate.donor.display_name ?? candidate.donor.node_label }}
        </button>
        <template v-if="candidate.donor.term">
          ({{ candidate.donor.term
          }}<template v-if="candidate.donor.name">, {{ candidate.donor.name }}</template>) on
          {{ candidate.donor.gene_count }} of
          {{ candidate.donor.member_gene_count.toLocaleString() }} of its genes
        </template>
        · rank {{ candidate.rank }} neighbour at {{ cosine(candidate.cosine) }}
      </p>

      <p v-if="confidenceSentence" class="infer-rate">{{ confidenceSentence }}</p>
      <p v-if="donorCaveat" class="muted cover">{{ donorCaveat }}</p>

      <!--
        ⛔ The ranks that carried nothing. For *E. coli* COG only 1,679 of 4,143 transfers come from
        rank 1, so a card that printed the donor alone would imply the nearest node supplied it.
      -->
      <p v-if="emptySteps.length" class="muted cover">
        Nearer {{ emptySteps.length === 1 ? "neighbour" : "neighbours" }} carried no
        {{ vocabularyLabel }}:
        {{ emptySteps.map((step) => `rank ${step.rank} at ${cosine(step.cosine)}`).join(", ") }}.
      </p>
    </template>

    <!--
      ⭐ The whole ladder, not just the cell that applies — David's explicit ask, so the reader can
      see the rate above in the company of the rates it is not.
    -->
    <details v-if="measuredCells.length" class="infer-ladder">
      <summary>How this was calibrated, on this catalogue</summary>
      <!--
        ⚠ Its own scroll container. At 100 % width the five columns spread into a sparse strip that
        reads as unrelated numbers, and at phone width the last one clips off the edge — neither
        visible to a test that computes no layout. The table hugs its content and the container
        scrolls instead.
      -->
      <div class="infer-scroll">
        <table class="infer-table">
          <thead>
            <tr>
              <th scope="col">ESM similarity</th>
              <th scope="col">quoted as</th>
              <th scope="col">agree</th>
              <th scope="col">pairs</th>
              <th scope="col">vs chance</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="cell in measuredCells"
              :key="`${cell.tier}-${cell.level}`"
              :class="{ 'infer-row-used': cell.tier === candidate?.tier && cell.level === candidate?.level }"
            >
              <td>{{ cell.tier }}</td>
              <td>{{ cell.level_label }}</td>
              <td class="num">{{ ((cell.agreement ?? 0) * 100).toFixed(1) }}%</td>
              <td class="num">{{ cell.pairs.toLocaleString() }}</td>
              <td class="num">{{ cell.lift === null ? "—" : `${Math.round(cell.lift)}×` }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p class="muted cover">
        Measured over every stored ESM neighbour edge in this catalogue, at every rank — not rank 1
        alone, because most transfers come from further out. A row is shown only where at least
        {{ ladder.min_pairs }} pairs could be compared.
      </p>
    </details>
  </div>
</template>
