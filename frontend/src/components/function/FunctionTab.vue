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

import type { AnnotationEntry, AnnotationKind, FunctionResponse } from "@/api/types";
import { GENE_ONTOLOGY_NAMESPACES, coverageSentence } from "@/lib/functionVocabulary";

import CogCard from "./CogCard.vue";
import EnzymeAndKeggCard from "./EnzymeAndKeggCard.vue";
import GeneOntologyCard from "./GeneOntologyCard.vue";

const props = defineProps<{
  /** The locus this panel is about — its name and label head the tab. */
  displayName: string;
  locusLabel: string;
  block: FunctionResponse | null;
  status: "idle" | "pending" | "ready" | "failed";
  /** The server's own sentence when the request failed. */
  failureDetail?: string | null;
}>();

const emit = defineEmits<{ retry: [] }>();

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
      <CogCard :coverage="block.coverage" :entries="entriesOf('cog_orthogroup')" />

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
        :verdict="block.go_verdicts[namespace]"
        :entries="geneOntologyEntries.get(namespace) ?? []"
      />

      <EnzymeAndKeggCard
        :coverage="block.coverage"
        :enzyme-entries="entriesOf('ec_number')"
        :kegg-entries="entriesOf('kegg_orthology')"
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
            GO terms are folded onto the metagenomics GO slim before the members are compared,
            because a term and its own child are annotated at different depths rather than in
            disagreement.
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
