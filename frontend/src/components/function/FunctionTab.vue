<script setup lang="ts">
/**
 * **The Function tab** — COG, GO and EC as Bakta assigned them.
 *
 * ⚠ **Not an independent source.** These are the same annotations the gene names come from, so they
 * are not evidence the clustering did not already see. The closing note says so on the page, because
 * a reader who takes the COG agreement as corroboration of the merge is double-counting one source.
 *
 * ⛔ **Three loading states, never two.** `pending` is not `failed`, and neither is "this locus has
 * no function annotation". `app.js:4604` — *"a panel that fails silently is one a reader will read
 * as 'this genome has nothing here', which is a different claim and a false one."* Here the false
 * claim would be that a locus is unannotated, which on this tab is itself a finding.
 */
import { computed } from "vue";

import type {
  AnnotationEntry,
  AnnotationKind,
  FunctionResponse,
  InferenceKind,
  OwnSupport,
  VocabularyInference,
} from "@/api/types";
import { INFERENCE_KINDS } from "@/api/types";
import {
  GENE_ONTOLOGY_NAMESPACES,
  coverageSentence,
  propagationCarriage,
  propagationRateSentence,
} from "@/lib/functionVocabulary";
import { prevalenceBandLabel } from "@/lib/prevalence";

import CogCard from "./CogCard.vue";
import EnzymeAndKeggCard from "./EnzymeAndKeggCard.vue";
import GeneOntologyCard from "./GeneOntologyCard.vue";
import InferredFunctionCard from "./InferredFunctionCard.vue";

/** How each inferable vocabulary is named to a reader — never the enum value. */
const VOCABULARY_LABELS: Readonly<Record<InferenceKind, string>> = {
  cog_orthogroup: "COG",
  gene_ontology_slim: "GO",
  kegg_orthology: "KEGG",
  ec_number: "EC number",
};

const props = defineProps<{
  /** The locus this panel is about — its name and label head the tab. */
  displayName: string;
  locusLabel: string;
  block: FunctionResponse | null;
  status: "idle" | "pending" | "ready" | "failed";
  /** The server's own sentence when the request failed. */
  failureDetail?: string | null;
}>();

const emit = defineEmits<{ retry: []; walk: [label: string] }>();

const inferenceByKind = computed(() => {
  const byKind = new Map<InferenceKind, VocabularyInference>();
  for (const entry of props.block?.inference.vocabularies ?? []) byKind.set(entry.annotation_kind, entry);
  return byKind;
});

/** The vocabularies this node has none of — the ones a neighbour can be asked about. */
const inferable = computed(() =>
  INFERENCE_KINDS.map((kind) => inferenceByKind.value.get(kind)).filter(
    (entry): entry is VocabularyInference => entry !== undefined && entry.own === null,
  ),
);

/**
 * ⭐ **One line per vocabulary this node states: that the call covers the whole syntelogue, and the
 * measured rate it is quoted at.**
 *
 * ⛔ **On EVERY node that carries a call, not only the one-gene ones.** This replaced a paragraph
 * that appeared only where a single gene carried the vocabulary and said its members *"cannot be
 * checked against each other"* — true, and the inverse of what the method claims. The claim and the
 * number are the same whether 1 or 17 genes carry it, because the rate answers *is the call right
 * for the members that carry nothing*, and that question does not change with the number that do.
 * Only the first clause differs (David, 2026-10-04).
 *
 * ⛔ No similarity appears here, and that is measured rather than overlooked: within-node unanimity
 * is 99.8-100 % at every band of both ESM and Bacformer, so a floor would be a gate doing no work.
 */
const propagationLines = computed(() =>
  INFERENCE_KINDS.map((kind) => ({ kind, own: inferenceByKind.value.get(kind)?.own ?? null }))
    .filter((row): row is { kind: InferenceKind; own: OwnSupport } => row.own !== null)
    .map(({ kind, own }) => ({
      kind,
      label: VOCABULARY_LABELS[kind],
      carriage: propagationCarriage(own),
      // ⚠ Band-free: every line of one node shares that node's band, so the lead names it once.
      rate: propagationRateSentence(own.propagation, { namesBand: false }),
    })),
);

/**
 * The lead's closing clause — what the figures beside the lines ARE.
 *
 * ⛔ `null` where no vocabulary has a rate at all, which is the ordinary state of a `rare` node:
 * a sentence explaining a figure that is not there would be the page talking about itself. The
 * lines then say, each for itself, why there is none.
 */
const propagationRateLead = computed(() => {
  const band = propagationLines.value.length
    ? (inferenceByKind.value.get(propagationLines.value[0]!.kind)!.own!.propagation.prevalence_band)
    : null;
  const measured = INFERENCE_KINDS.map((kind) => inferenceByKind.value.get(kind)?.own ?? null).some(
    (own) => own !== null && own.propagation.rate !== null,
  );
  if (band === null || !measured) return null;
  return (
    `The figure beside each is how often this catalogue's ${prevalenceBandLabel(band)} ` +
    "syntelogues agree where two or more of their genes do carry one."
  );
});

function entriesOf(kind: AnnotationKind): readonly AnnotationEntry[] {
  return props.block?.annotations[kind] ?? [];
}

/**
 * ⛔ Split by namespace HERE rather than trusting the order they arrived in. They come back in one
 * `gene_ontology_slim` list carrying their own namespace, and reading them positionally would file
 * a cellular-component class under molecular function — with three cards that all look right.
 */
const geneOntologyEntries = computed(() => {
  const byNamespace = new Map<string, AnnotationEntry[]>();
  for (const entry of entriesOf("gene_ontology_slim")) {
    const namespace = entry.gene_ontology_namespace;
    if (namespace === undefined) continue;
    byNamespace.set(namespace, [...(byNamespace.get(namespace) ?? []), entry]);
  }
  return byNamespace;
});

/**
 * ⭐ Whether this locus has ANY GO at all, across the three namespaces.
 *
 * All three cards, or one line — never three empty cards saying the same thing three times. But the
 * test is the annotated-gene COUNTS, not the term lists: a namespace can have coverage and still
 * list nothing, and it is the coverage that decides whether there is anything to report.
 */
const hasAnyGeneOntology = computed(() =>
  props.block === null
    ? false
    : GENE_ONTOLOGY_NAMESPACES.some(
        (namespace) => props.block!.coverage.go_annotated_gene_count[namespace] > 0,
      ),
);

const noGeneOntologySentence = computed(() =>
  props.block === null
    ? ""
    : coverageSentence(0, props.block.coverage.gene_count, "a GO term"),
);
</script>

<template>
  <!-- ⚠ No `.wrap` of its own: the page's view panel already is one, and a second nests the gutter. -->
  <div>
    <div class="func-head">
      <h2>Function</h2>
      <span class="lid">{{ displayName }} ·{{ locusLabel }}</span>
    </div>

    <!-- ⛔ Failure, waiting and "nothing here" are three sentences, and only the third is a finding. -->
    <p v-if="status === 'failed'" class="pop-error" role="alert">
      {{ failureDetail ?? "the function annotation did not load" }}
      <button type="button" class="pop-retry" @click="emit('retry')">Try again</button>
    </p>
    <p v-else-if="status !== 'ready' || block === null" class="muted">Loading the function annotation…</p>

    <template v-else>
      <!--
        ⛔ BEFORE any card, because it is the claim the cards below are instances of rather than a
        footnote to them. A reader who scrolls no further should still have read it.
      -->
      <section v-if="propagationLines.length" class="infer-spread">
        <p class="infer-spread-lead">
          <strong>Each of these applies to the whole syntelogue.</strong> A syntelogue is one gene
          seen in many genomes, so a function carried by some of its members is carried by all of
          them — that propagation is the method, not a gap in it.
          <template v-if="propagationRateLead">{{ propagationRateLead }}</template>
        </p>
        <ul>
          <li v-for="line in propagationLines" :key="`propagates-${line.kind}`">
            <strong>{{ line.label }}</strong> — {{ line.carriage }} · {{ line.rate }}
          </li>
        </ul>
      </section>

      <CogCard
        :coverage="block.coverage"
        :entries="entriesOf('cog_orthogroup')"
        :propagation="inferenceByKind.get('cog_orthogroup')?.own?.propagation ?? null"
      />

      <div v-if="!hasAnyGeneOntology" class="card">
        <h3 class="sub-head">GO</h3>
        <p class="muted cover">{{ noGeneOntologySentence }}</p>
      </div>
      <GeneOntologyCard
        v-for="namespace in hasAnyGeneOntology ? GENE_ONTOLOGY_NAMESPACES : []"
        :key="namespace"
        :namespace="namespace"
        :annotated-gene-count="block.coverage.go_annotated_gene_count[namespace]"
        :gene-count="block.coverage.gene_count"
        :entries="geneOntologyEntries.get(namespace) ?? []"
      />

      <EnzymeAndKeggCard
        :coverage="block.coverage"
        :enzyme-entries="entriesOf('ec_number')"
        :kegg-entries="entriesOf('kegg_orthology')"
      />

      <!-- ⛔ AFTER every card of real annotation, never among them. What a neighbour SUGGESTS. -->
      <InferredFunctionCard
        v-for="entry in inferable"
        :key="`infer-${entry.annotation_kind}`"
        :inference="entry"
        :ladder="block.calibration[entry.annotation_kind]"
        :vocabulary-label="VOCABULARY_LABELS[entry.annotation_kind]"
        @walk="(label) => emit('walk', label)"
      />

      <!--
        ⛔ This was ONE five-line paragraph doing four jobs — provenance, independence, the GO-slim
        method and the coverage principle — and a reader looking for any one of them had to read all
        four (David, 2026-10-03: "The text on current page is a bit hard to absorb. A cleaner layout
        needed."). Split and labelled: the same four claims, each findable on its own.
      -->
      <dl class="func-notes">
        <div>
          <dt>Source</dt>
          <dd>
            COG, GO and EC as Bakta assigned them, from the same annotation the gene names come
            from — not an independent source, and not evidence the clustering did not already see.
          </dd>
        </div>
        <div>
          <dt>GO slim</dt>
          <dd>
            GO terms are folded onto the metagenomics GO slim — 111 classes — so members are
            compared on one vocabulary rather than on 38,000 terms of varying depth.
            <strong>The slim is a selection, not a hierarchy:</strong> <em>plasma membrane</em> and
            <em>membrane</em> are both classes in it, so two members can carry different classes and
            still be saying compatible things. Agreement here does not yet account for that.
          </dd>
        </div>
        <div>
          <dt>Coverage</dt>
          <dd>
            Stated first on every card above: a gene Bakta never annotated says nothing either way,
            which is a different finding from the members disagreeing.
          </dd>
        </div>
        <div>
          <dt>Density</dt>
          <!--
            ⚠ Measured on ecoli-nuna4, loci of ≥10 members that carry the axis at all (2026-10-03):
            COG median 1.000 (n=4,031, 75.9 % on every gene) · EC median 1.000 (n=2,114, 84.8 %) ·
            Pfam median 1.000 (n=6,304, 96.6 %) · KEGG Q1 0.091, median 0.737, only 36.5 % on every
            gene (n=908). KEGG is the one axis where a reader must look at the count rather than
            assume it, so it is the one named here.
          -->
          <dd>
            Where COG, EC or Pfam appear at all, the median locus carries them on <em>every</em>
            member gene — so a partial count there is worth noticing. KEGG is the exception and is
            routinely partial: of the loci carrying any KEGG, a quarter carry it on under a tenth of
            their genes. Read the count on the card, never the terms alone.
          </dd>
        </div>
      </dl>
    </template>
  </div>
</template>
