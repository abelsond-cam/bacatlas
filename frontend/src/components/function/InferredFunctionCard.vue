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
 * `nuna/docs/model_evaluation/function_inference.md` reports from. Not a constant lifted out of
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
import { indefiniteArticle } from "@/lib/formatting";
import { propagationRateSentence } from "@/lib/functionVocabulary";

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
 * ⭐⭐ **One line per rung the donor can be read at, each with its own verdict.**
 *
 * ⛔ Until 2026-10-06 this card showed ONE suggestion per vocabulary, which for GO meant averaging
 * three unrelated claims. Measured at 0.94-0.96, cellular component agrees 84.2 % / 86.2 % of the
 * time and molecular function 60.3 % / 70.5 % — so a single pooled 82 % was letting the weak one
 * through on the strength of the other two. They are separate questions and now get separate lines.
 *
 * ⛔ **A refused rung still gets a line**, carrying the rate that refused it. A reader who sees
 * nothing cannot tell "not considered" from "considered and declined", and the second is the whole
 * point of having a floor.
 */
const offerLines = computed(() => {
  const found = candidate.value;
  if (found === null) return [];
  const floor = Math.round(props.ladder.calling_floor * 100);
  return found.offers.map((offer) => {
    const cell = offer.calibration;
    const label = cell?.level_label ?? props.vocabularyLabel;
    if (cell === null || cell.agreement === null) {
      return {
        key: `${offer.facet ?? ""}-${offer.level}`,
        suggested: false,
        label,
        value: offer.value.join(", "),
        sentence:
          "Too few measured pairs at this similarity to state how often a transfer of this kind " +
          `is right${cell === null ? "" : ` (only ${cell.pairs.toLocaleString()})`} — treat it as ` +
          "a neighbour, not a suggestion.",
      };
    }
    const percent = (cell.agreement * 100).toFixed(1);
    const pairs = cell.pairs.toLocaleString();
    // ⛔ **The lift is deliberately absent** (David, 2026-10-05). "agree 24.3 % of the time, which
    // is 145× more often than a random annotated node" is two true numbers telling opposite
    // stories, and the second is the one a reader carries away. Lift lives in the ladder below,
    // where the vocabularies are compared with each other.
    return {
      key: `${offer.facet ?? ""}-${offer.level}`,
      suggested: offer.suggested,
      label,
      value: offer.value.join(", "),
      sentence: offer.suggested
        ? `Nodes this similar agree on their ${label} ${percent}% of the time (${pairs} measured pairs in this catalogue).`
        : `Nodes this similar agree on their ${label} only ${percent}% of the time (${pairs} measured pairs) — below the ${floor}% this page needs, so this one is not suggested.`,
    };
  });
});

/**
 * ⭐ **Why nothing is suggested, when a perfectly close neighbour was found.**
 *
 * ⛔ Two different refusals wore one sentence until 2026-10-05. *"Too remote to call"* is true of a
 * neighbour below 0.90; it is false of `wza` at **0.929**, which is close, carries a KEGG, and is
 * refused because nodes that similar agree only **24.3 %** of the time. Saying "too remote" there
 * blames the geometry for what the measurement decided, and hides the number the reader needs.
 *
 * ⚠ Only for the case where the donor states nothing READABLE at all — where it does, each rung
 * explains itself in `offerLines`.
 */
const refusalSentence = computed(() => {
  const found = candidate.value;
  if (found === null || found.offers.length > 0) return null;
  return (
    `The nearest node carrying ${indefiniteArticle(props.vocabularyLabel)} ` +
    `${props.vocabularyLabel} is at ${cosine(found.cosine)} similarity — too remote to call, so ` +
    "nothing is suggested."
  );
});

/**
 * ⭐ **The donor's own node, judged by the same rule as the focal one.**
 *
 * ⛔ This used to fire only where the donor rested on ONE gene, and said it *"was never checked
 * against itself either"* — the same retired claim the Function tab led with, in the place it does
 * the most damage: it implied a one-gene donor was a weaker donor, when a donor carries a
 * vocabulary by the **any-gene** rule precisely because the call covers its whole node. Now it
 * fires wherever the donor's own call is itself propagated, and states the rate that justifies it.
 *
 * ⚠ `null` where every one of the donor's genes carries the call: there is no propagation on that
 * side to report, and a sentence saying so on every card would be noise.
 */
const donorPropagation = computed(() => {
  const found = candidate.value;
  // ⛔ Silent where nothing is suggested. Explaining that the donor is a legitimate donor, directly
  // under a paragraph refusing to take anything from it, reads as the card arguing with itself.
  if (found === null || found.suggested_count === 0) return null;
  const donor = found.donor;
  if (donor.annotated_gene_count >= donor.member_gene_count) return null;
  const rate = propagationRateSentence(donor.propagation);
  // ⚠ Two endings, because the general rule and the striking case are different sentences. "Why a
  // neighbour resting on one annotated gene is still a donor" is a non-sequitur beside a donor
  // whose call is carried by 57 of its 100 genes — it is the one-gene donor that needs saying.
  const rule =
    donor.annotated_gene_count === 1
      ? "which is why a neighbour resting on one annotated gene is still a donor"
      : "which is the rule by which it counts as a donor at all";
  return (
    `The neighbour's own ${props.vocabularyLabel} is carried by ` +
    `${donor.annotated_gene_count.toLocaleString()} of its ` +
    `${donor.member_gene_count.toLocaleString()} genes and covers that whole node — ${rule}. ` +
    // ⚠ The helper returns a CLAUSE, because its other caller lists several after a dash. Here it
    // ends a sentence of its own, so it is capitalised and stopped at the point of use.
    `${rate[0]!.toUpperCase()}${rate.slice(1)}.`
  );
});

/** Which of the card's three outcomes this is — see the template. */
const heading = computed(() => {
  if (candidate.value === null) return "nothing near enough to infer from";
  if (candidate.value.suggested_count === 0) return "nothing suggested";
  return "inferred from a neighbour";
});

const measuredCells = computed(() => props.ladder.cells.filter((cell) => cell.agreement !== null));

/**
 * Which rows of the ladder below this card actually judged THIS suggestion — up to three for GO.
 * ⚠ Keyed on the facet too: a key of (tier, level) alone would light all three GO rows of a tier
 * when only one of them was used.
 */
const usedCells = computed(
  () =>
    new Set(
      (candidate.value?.offers ?? []).map(
        (offer) => `${candidate.value!.tier}-${offer.facet ?? ""}-${offer.level}`,
      ),
    ),
);
</script>

<template>
  <div class="card infer">
    <!--
      ⛔ The heading must not repeat the card above it. The real COG card already says "None of these
      102 genes carries a COG assignment"; a second card headed "none assigned here" reads as the
      same finding twice, which is what the first render of this looked like. It says instead where
      the suggestion COMES FROM, which is the one thing the card above cannot say.
    -->
    <!--
      ⛔ THREE headings, because the card now has three outcomes and two of them suggest nothing.
      It read "inferred from a neighbour … [inferred]" over a paragraph explaining that nothing was
      inferred — the floor had refused the call. A card must not announce a suggestion it is about
      to withhold, and the "inferred" tag belongs only where something was.
    -->
    <h3 class="sub-head">
      {{ vocabularyLabel }} — {{ heading }}
      <span v-if="candidate !== null && candidate.suggested_count > 0" class="infer-tag">inferred</span>
    </h3>

    <!-- 3 · nothing within reach. A finding, stated as one. -->
    <p v-if="candidate === null" class="muted cover">
      None of the {{ inference.walk.length }} nearest nodes by ESM similarity carries
      {{ indefiniteArticle(vocabularyLabel) }} {{ vocabularyLabel }} either, so there is nothing to
      transfer. That is common and it is
      informative: unannotated nodes sit next to unannotated nodes far more often than chance, which
      is why widening the neighbour list is not the missing piece.
    </p>

    <template v-else>
      <!--
        1 · one line per rung the donor can be read at. ⛔ GO has three and they do not stand or
        fall together, so a refused rung is SHOWN with the rate that refused it rather than omitted
        — a reader who sees nothing cannot tell "not considered" from "considered and declined".
      -->
      <template v-for="line in offerLines" :key="line.key">
        <p v-if="line.suggested" class="infer-value">
          <span class="infer-lead">Suggested {{ line.label }}:</span>
          <strong>{{ line.value }}</strong>
        </p>
        <p v-else class="muted cover">
          <strong>{{ line.label }} — not suggested.</strong> {{ line.sentence }}
        </p>
        <p v-if="line.suggested" class="infer-rate">{{ line.sentence }}</p>
      </template>
      <!-- 2 · the donor states nothing readable at all: a different refusal, in its own words -->
      <p v-if="refusalSentence" class="muted cover">{{ refusalSentence }}</p>

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

      <p v-if="donorPropagation" class="muted cover">{{ donorPropagation }}</p>

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
              :key="`${cell.tier}-${cell.facet ?? ''}-${cell.level}`"
              :class="{ 'infer-row-used': usedCells.has(`${cell.tier}-${cell.facet ?? ''}-${cell.level}`) }"
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
