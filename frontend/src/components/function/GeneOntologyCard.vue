<script setup lang="ts">
/**
 * **GO, one card per namespace** — molecular function, biological process, cellular component.
 *
 * ⛔ **All three, or none.** A namespace with nothing to show still *states* that it has nothing,
 * because `no coverage` and `the members disagree` are different findings and only the second is
 * evidence. A locus with no GO at all is the one case that collapses to a single line, because three
 * identical empty cards say the same thing three times.
 *
 * ⭐ **The terms are GO SLIM classes, not raw terms**, folded onto the metagenomics slim because a
 * term and its own child are annotated at different *depths* rather than in disagreement.
 *
 * ⛔⛔ **The verdict chip is RETIRED (David, 2026-10-05), and the paragraph above says exactly why
 * it had to be.** The principle was right and the fold did not deliver it: `goslim_metagenomics`
 * keeps `membrane`, `plasma membrane` and `outer membrane` as *siblings*, and the three namespace
 * ROOTS besides, so a term and its own parent survived the fold and were then compared by
 * `worst_relation` — Pfam's domain-architecture comparator, which has no ontology and for which
 * `disjoint` is the condemning verdict. `gumC` read *"classes differ"* for 16 genes saying plasma
 * membrane and one saying membrane; 11 of the 43 `disjoint` verdicts across both catalogues involve
 * a root, a term that is an ancestor of everything. A chip that is wrong on its loudest case is
 * worse than no chip. ⚠ What replaces it is not silence but a **correct measurement** — the
 * propagation rate above these cards — and that work is where the ontology now belongs.
 *
 * ⛔ **Coverage first, every time**, and it is still the thing the card leads with.
 */
import { computed } from "vue";

import type { AnnotationEntry, GeneOntologyNamespace } from "@/api/types";
import {
  GENE_ONTOLOGY_NAMESPACE_LABEL,
  coverageParts,
  geneOntologyTermUrl,
} from "@/lib/functionVocabulary";

import CountTable from "../shared/CountTable.vue";

const props = defineProps<{
  namespace: GeneOntologyNamespace;
  annotatedGeneCount: number;
  geneCount: number;
  entries: readonly AnnotationEntry[];
}>();

const label = computed(() => GENE_ONTOLOGY_NAMESPACE_LABEL[props.namespace]);
const coverage = computed(() =>
  coverageParts(props.annotatedGeneCount, props.geneCount, `a ${label.value} term`),
);

const rows = computed(() =>
  props.entries.map((entry) => ({
    key: entry.term,
    count: entry.gene_count,
    term: entry.term,
    // ⚠ The class NAME is the readable half; the accession is what a reader follows. Both, always —
    // `GO:0016020` alone is unreadable and "membrane" alone is unlookupable.
    name: entry.name ?? entry.term,
  })),
);
</script>

<template>
  <div class="card">
    <h3 class="sub-head">GO — {{ label }}</h3>
    <p class="muted cover"><b v-if="coverage.emphasis">{{ coverage.emphasis }}</b>{{ coverage.rest }}</p>
    <!--
      ⚠ **"commonest", because the list is CAPPED and the card cannot yet say at what.** The export
      keeps `TOP_GO = 4` classes per namespace and 688 ecoli / 269 kp loci hit that cap, with
      nothing on the page saying classes were dropped — which is how `fcl` came to show four
      identical-looking rows beside a chip claiming the classes differed. Stage 2 carries the total
      so this can read "4 of 7"; until then the card must not imply the list is complete.
    -->
    <p v-if="rows.length" class="muted cover">
      The commonest classes here, by how many member genes carry each.
    </p>
    <CountTable v-if="rows.length" :headings="['class', 'GO']" :rows="rows" :total="geneCount">
      <template #cells="{ row }">
        <td>{{ row.name }}</td>
        <td class="acc">
          <a class="acc-link" :href="geneOntologyTermUrl(row.term)" target="_blank" rel="noopener" :title="`${row.term} on AmiGO`">{{ row.term }}</a>
        </td>
      </template>
    </CountTable>
  </div>
</template>
