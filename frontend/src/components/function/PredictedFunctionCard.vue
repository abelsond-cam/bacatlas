<script setup lang="ts">
/**
 * **Predicted function** — what annotation-transfer suggests for this locus, each with its own accuracy.
 *
 * ⭐ It replaced the 5-nearest-loci ladder (David, 2026-11-14). annotation-transfer v0.1 copies the
 * classes of ONE similar locus — the donor, whose medoid gene is the nearest in both ESM and Bacformer
 * space by a fitted score — and serves the copy only where its calibrated accuracy is 0.80 or more.
 *
 * ⛔ **Separate from every other card on this tab, and it must stay separate.** The cards above print
 * what this locus's own genes SAY; this one prints what the model PREDICTS. So it keeps the dashed
 * `.infer` frame and a "predicted" tag, and its terms sit on the tinted strip, never in the headline
 * weight the real terms use.
 *
 * ⛔ **Five answers, one sentence each, and none may read like another** (`PredictionsBlock.state`; the
 * review of 9cfa481 found three of them sharing one sentence):
 *  1. `loaded`, with calls → one line per claim, each with ITS OWN estimated accuracy (no bands — David,
 *     2026-10-28) and its donor, linked;
 *  2. `loaded`, no calls → the model DECLINED this locus;
 *  3. `none_produced` → the model produced none for this catalogue: a finding about the model;
 *  4. `not_loaded` → none are loaded into this database: says nothing about the model;
 *  5. `unknown_model` → loaded from a version this build cannot show, so none are shown.
 *
 * ⚠ The declined sentence claims nothing about the locus's own annotation. It used to say the model
 * predicts only where the locus's own genes state nothing — false for 426 of 35,588 served GO predictions
 * on the probes, where a namespace root won the producer's tie-break (`served_locus_annotation`'s module
 * docstring). Only what is true of every served prediction is said.
 */
import { computed } from "vue";

import type { PredictionsBlock, ServedPrediction } from "@/api/types";
import {
  COG_CATEGORIES,
  GENE_ONTOLOGY_NAMESPACE_LABEL,
  estimatedAccuracyPercent,
} from "@/lib/functionVocabulary";

const props = defineProps<{ predictions: PredictionsBlock }>();

const emit = defineEmits<{ walk: [label: string] }>();

/** How the claim is named to a reader. ⛔ A COG CATEGORY, never "COG" alone: the card above is accessions. */
function claimLabel(call: ServedPrediction): string {
  if (call.annotation_kind === "cog_functional_category") return "COG category";
  return call.gene_ontology_namespace === null
    ? "GO"
    : `GO ${GENE_ONTOLOGY_NAMESPACE_LABEL[call.gene_ontology_namespace]}`;
}

/** A COG letter in words; a GO class by its name, its accession where the reference has none. */
function termText(call: ServedPrediction, index: number): string {
  const term = call.terms[index]!;
  if (call.annotation_kind === "cog_functional_category") {
    const name = COG_CATEGORIES[term.value];
    return name === undefined ? term.value : `${term.value} — ${name}`;
  }
  return term.name ?? term.value;
}

/**
 * ⭐ The headline and what it implies, kept apart.
 *
 * For GO the stored set is a CLOSED lineage — the donor's class first, then every slim class above it —
 * and the accuracy is scored on that whole set. Leading with the first and listing the rest as "within"
 * shows what the donor states AND what the number covers; printing them as equals ("plasma membrane,
 * membrane") would read as two findings. COG letters are peers, so every one is a headline.
 */
const lines = computed(() =>
  props.predictions.calls.map((call, index) => {
    const isGeneOntology = call.annotation_kind === "gene_ontology_slim";
    const headline = isGeneOntology
      ? termText(call, 0)
      : call.terms.map((_term, termIndex) => termText(call, termIndex)).join(" · ");
    const broader = isGeneOntology
      ? call.terms.slice(1).map((_term, termIndex) => termText(call, termIndex + 1))
      : [];
    return {
      key: `${call.annotation_kind}-${call.gene_ontology_namespace ?? ""}-${index}`,
      label: claimLabel(call),
      headline,
      broader,
      accuracy: estimatedAccuracyPercent(call.estimated_accuracy),
      donor: call.donor,
      // ⛔ Locus labels are catalogue-private: a donor in ANOTHER catalogue is reached by its key, and
      // a walk inside this page would land on a different locus wearing the same number.
      href: call.donor.same_catalogue
        ? null
        : `?${new URLSearchParams({ catalogue: call.donor.catalogue_key }).toString()}#${call.donor.node_label}`,
      similarity: `ESM ${call.esm_medoid_similarity.toFixed(3)} · Bacformer ${call.bacformer_medoid_similarity.toFixed(3)}`,
    };
  }),
);

/** "80 %" — the bar, in the page's own percent form. */
const bar = computed(() =>
  props.predictions.serve_cut === null ? null : estimatedAccuracyPercent(props.predictions.serve_cut),
);
</script>

<template>
  <div class="card infer predict">
    <h3 class="sub-head">
      Predicted function
      <span v-if="lines.length" class="infer-tag">predicted</span>
    </h3>

    <!-- 4 · not loaded HERE. A statement about this database, never about the model. -->
    <p v-if="predictions.state === 'not_loaded'" class="muted cover predict-state">
      Annotation-transfer predictions are not loaded for this catalogue, so none can be shown. That says
      nothing about whether the model has predictions for it.
    </p>

    <!-- 3 · the model produced none for this catalogue. A finding about the model. -->
    <p v-else-if="predictions.state === 'none_produced'" class="muted cover predict-state">
      The annotation model produced no predictions for this catalogue, so no locus in it carries a
      predicted function.
    </p>

    <!-- 5 · a version this build cannot show: fail closed, and say which. -->
    <p v-else-if="predictions.state === 'unknown_model' || bar === null" class="muted cover predict-state">
      This catalogue's predictions come from {{ predictions.annotation_model_label ?? "an annotation model" }},
      which this version of BacAtlas cannot show yet, so none are shown.
    </p>

    <!-- 2 · loaded, and the model declined this locus. -->
    <p v-else-if="!lines.length" class="muted cover predict-state">
      {{ predictions.annotation_model_label }} makes no prediction for this locus. A prediction is shown
      only where its estimated accuracy is {{ bar }} or more.
    </p>

    <!-- 1 · the predictions, one line per claim, each with its own number. -->
    <template v-else>
      <p class="muted cover predict-model">
        From {{ predictions.annotation_model_label }}: each copied from the one most similar
        annotated locus, and shown only at an estimated accuracy of {{ bar }} or more.
      </p>
      <div v-for="line in lines" :key="line.key" class="predict-call">
        <p class="infer-value">
          <span class="infer-lead">{{ line.label }}:</span>
          <strong>{{ line.headline }}</strong>
          <span v-if="line.broader.length" class="predict-broader"> — within {{ line.broader.join(", ") }}</span>
        </p>
        <p class="infer-rate">
          Estimated accuracy <strong class="predict-accuracy">{{ line.accuracy }}</strong>
        </p>
        <p class="infer-from">
          from
          <button
            v-if="line.href === null"
            type="button"
            class="infer-walk"
            @click="emit('walk', line.donor.node_label)"
          >{{ line.donor.display_name }}</button>
          <a v-else class="infer-walk" :href="line.href">{{ line.donor.display_name }}</a>
          (locus {{ line.donor.node_label }}, {{ line.donor.gene_count.toLocaleString() }}
          {{ line.donor.gene_count === 1 ? "gene" : "genes" }}<template v-if="!line.donor.same_catalogue">,
            <em>{{ line.donor.species_name }}</em></template>)
          · medoid similarity {{ line.similarity }}
        </p>
      </div>
      <p class="muted cover predict-meaning">
        Estimated accuracy is the model's calibrated chance that this locus carries at least one of the
        classes shown — measured on loci whose own annotation was held out, never on this one. A GO
        prediction counts as right if it matches at its broadest class.
      </p>
    </template>
  </div>
</template>
