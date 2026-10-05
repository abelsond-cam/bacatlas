<script setup lang="ts">
/**
 * **EC and KEGG** — thin enough that they share one compact card rather than each getting a table.
 *
 * ⭐ **Both are NAMED now, and the asymmetry that used to be here was two separate mistakes.**
 * EC showed seven generic words because a note in this codebase recorded ExPASy ENZYME as CC BY-ND;
 * it is **CC BY 4.0**, so the real names are vendored and `2.7.10.-` reads *Protein-tyrosine
 * kinases*. KEGG showed nothing at all on a licence that is genuinely unresolved — and still is.
 *
 * ⛔⛔ **KEGG's licence is NOT settled.** Its terms require an academic service provider licence
 * from anyone *offering services*, and serving a KO description is offering one. Named here on
 * David's decision of 2026-10-05 for an internal prototype, with the question carried on the page
 * rather than quietly resolved. **Nothing ships to BacAtlas.org until the subscription is
 * confirmed.** The published static pages still link out and name nothing.
 *
 * ⚠ **The count is SEPARATED from the accession.** Run together it reads as one token: `K15540`
 * beside a count of 77 came back from a reader as the question *"what is K1554077?"*
 *
 * ⚠ **The two rows are different KINDS of claim, and sharing a card invites confusing them.** An EC
 * number classifies a *reaction* on four hierarchical levels; a KEGG KO identifies a *gene family*.
 * Neither implies the other, and the card says so once rather than leaving the heading to imply it.
 *
 * ⭐ **An absent row is STATED, not omitted.** `gumC` carries an EC and no KEGG, and the card simply
 * dropped the KEGG line — leaving a reader unable to tell *"no KEGG here"* from *"KEGG not shown"*,
 * which is the question David asked of this card (*"we don't need to infer KEGG do we??? We have it
 * for our nodes??"*). The answer is that KEGG is the thinnest vocabulary Bakta emits — **8.2 % of
 * the 1,436,421 genes in the two catalogues, against COG's 67.7 %** — and that is said where the
 * missing row would have been.
 */
import { computed } from "vue";

import type { AnnotationEntry, FunctionResponse } from "@/api/types";
import { enzymeClassSummary, enzymeCommissionUrl, keggOrthologyUrl } from "@/lib/functionVocabulary";

const props = defineProps<{
  coverage: FunctionResponse["coverage"];
  enzymeEntries: readonly AnnotationEntry[];
  keggEntries: readonly AnnotationEntry[];
}>();

/**
 * ⛔ **The served name first, the seven-word gloss only as a fallback.** The reference names 1,409
 * of the 1,410 EC codes these catalogues carry; the gloss exists for the one it does not.
 */
function enzymeGloss(entry: AnnotationEntry): string | null {
  return entry.name ?? enzymeClassSummary(entry.term);
}

const lines = computed(() =>
  [
    {
      key: "ec",
      label: "EC",
      entries: props.enzymeEntries,
      annotated: props.coverage.ec_annotated_gene_count,
      url: enzymeCommissionUrl,
      title: (term: string) => `EC ${term} on ExPASy`,
      gloss: enzymeGloss,
      // ⚠ Said once per row, where the absence is, rather than as a footnote three cards away.
      absent: "No EC number on any of this node's genes.",
    },
    {
      key: "kegg",
      label: "KEGG KO",
      entries: props.keggEntries,
      annotated: props.coverage.kegg_annotated_gene_count,
      url: keggOrthologyUrl,
      title: (term: string) => `${term} on KEGG`,
      // ⚠ No fallback: a KO id encodes nothing a reader can decode, so an unnamed one stays bare.
      gloss: (entry: AnnotationEntry) => entry.name,
      absent:
        "No KEGG orthology on any of this node's genes — KEGG is the thinnest label Bakta " +
        "emits, on 8.2 % of genes against COG's 67.7 %.",
    },
  ].filter((line) => line.entries.length > 0 || props.coverage.gene_count > 0),
);

const geneCount = computed(() => props.coverage.gene_count);
const anyEntries = computed(
  () => props.enzymeEntries.length > 0 || props.keggEntries.length > 0,
);
</script>

<template>
  <!-- ⚠ The whole card is absent when neither vocabulary has anything, rather than present and
       empty: two headings over nothing say less than no card at all. Where ONE is present the other
       states its absence, because that is a finding and a dropped row is not. -->
  <div v-if="anyEntries" class="card">
    <h3 class="sub-head">EC and KEGG</h3>
    <div v-for="line in lines" :key="line.key" class="kv">
      <span class="k">{{ line.label }}</span>
      <span class="v">
        <template v-if="line.entries.length">
          <template v-for="(entry, index) in line.entries" :key="entry.term">
            <template v-if="index"> · </template>
            <a
              class="acc-link"
              :href="line.url(entry.term)"
              target="_blank"
              rel="noopener"
              :title="line.title(entry.term)"
            >{{ entry.term }}</a>
            <span class="alt-n">×{{ entry.gene_count }}</span>
            <span v-if="line.gloss(entry)" class="ec-gloss">{{ line.gloss(entry) }}</span>
          </template>
        </template>
        <span v-else class="muted">{{ line.absent }}</span>
      </span>
      <span v-if="line.entries.length" class="cov">
        {{ line.annotated }} of {{ geneCount }} annotated
      </span>
    </div>
    <p class="muted cover">
      An EC number classifies the reaction, on four levels; a KEGG KO identifies a gene family.
      They are different vocabularies and neither implies the other — a locus can carry one, both or
      neither, and this card is absent where it carries neither.
    </p>
    <!-- ⭐ The CC BY condition, discharged. It is the whole of ExPASy's licence, so it belongs on
         the page that uses the names and not only in a module docstring. -->
    <p class="muted cover">
      EC names from
      <a href="https://enzyme.expasy.org/" target="_blank" rel="noopener">ExPASy ENZYME</a>, SIB
      Swiss Institute of Bioinformatics (CC BY 4.0). KO names from
      <a href="https://www.kegg.jp/" target="_blank" rel="noopener">KEGG</a>, Kanehisa Laboratories
      — <strong>licence unresolved, internal prototype only.</strong>
    </p>
  </div>
</template>
