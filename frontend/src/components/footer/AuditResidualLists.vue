<script setup lang="ts">
/**
 * Two lists of what the audit MEASURED, each locus one click away — the model's own residuals, on its own
 * page.
 *
 * ⭐ **What they list (David, 2026-11-26).** These were *"Every locus grouped on context alone"* and
 * *"Every Pfam conflict"*. At 71,435 genomes both were decided by stray genes — kp-71k locus 4008, Wzi on
 * 69,825 of 69,884 genes, sat in both — and David ruled them audit bugs. The lists now name facts, never a
 * verdict:
 *
 * 1. **Loci grouped by neither MMseqs nor a shared Pfam architecture** — the ladder's last rung, `esm_only`.
 * 2. **Loci where at most 95 % of genes share one Pfam architecture** — among loci with members and with
 *    Pfam. ⚠ "at most", not "fewer than": the audit's test is a strict `> 95 %`, so a locus at exactly
 *    95 % (19 of 20 genes) is listed, and "fewer than" would be false of it.
 *
 * ⛔ **Two lists, never one, and they overlap.** A locus in both appears in both.
 *
 * ⚠ **Fetched when a list is first opened.** The COUNTS come from the shell (`residual_locus_counts`),
 * counted from the same clauses the list endpoint selects with — so a closed list still says how long it
 * is, and the number cannot disagree with the rows it opens to.
 */
import { computed, ref } from "vue";

import { fetchAuditResiduals } from "@/api/client";
import type { Failure } from "@/api/result";
import type { AuditResidualsResponse, CatalogueKey, ResidualLocusCounts, ResidualLocusRow } from "@/api/types";
import { cosineLowerBound } from "@/lib/evidenceVocabulary";
import { percentKeepingSide } from "@/lib/formatting";
import { prevalenceBandLabel } from "@/lib/prevalence";

const props = defineProps<{
  catalogueKey: CatalogueKey;
  counts: ResidualLocusCounts;
}>();

const emit = defineEmits<{ go: [locusLabel: string] }>();

const lists = ref<AuditResidualsResponse | null>(null);
const failed = ref<Failure | null>(null);
const isLoading = ref(false);

// ⚠ The shell's count until the rows arrive, then the rows' OWN count — a summary must never promise ten
// over a list of nine, even if the two were somehow read a deploy apart.
const esmOnlyCount = computed(() => lists.value?.esm_only.length ?? props.counts.esm_only);
const notSharedCount = computed(
  () => lists.value?.pfam_architecture_not_shared.length ?? props.counts.pfam_architecture_not_shared,
);

async function load(event: Event): Promise<void> {
  if (!(event.target as HTMLDetailsElement).open || lists.value !== null || isLoading.value) return;
  isLoading.value = true;
  const result = await fetchAuditResiduals(props.catalogueKey);
  isLoading.value = false;
  if (result.ok) {
    lists.value = result.value;
    failed.value = null;
  } else {
    failed.value = result;
  }
}

/** `render_page._goto_rows`: enough evidence beside the name to judge before clicking. */
function evidence(row: ResidualLocusRow): string {
  const bits = [prevalenceBandLabel(row.prevalence_band), `${row.gene_count.toLocaleString()} genes`];
  // ⚠ Omitted when not measured — `?? 0` would have printed a measurement nobody made.
  if (row.uniref50_family_count !== null) bits.push(`${row.uniref50_family_count} UniRef50`);
  if (row.pfam_architecture_count !== null) bits.push(`${row.pfam_architecture_count} architectures`);
  // ⭐ The two facts the lists are ABOUT, so each row states its own reason. The share is over ALL genes.
  if (row.pfam_dominant_domain_gene_count !== null && row.gene_count > 0) {
    bits.push(`${percentKeepingSide(row.pfam_dominant_domain_gene_count / row.gene_count, 0.95)} on one Pfam set`);
  }
  if (row.esm_member_similarity_p05 !== null) {
    bits.push(`95% ESM ≥ ${cosineLowerBound(row.esm_member_similarity_p05)}`);
  }
  if (row.syntenic_a5 !== null) bits.push(`A5 ${row.syntenic_a5.toFixed(2)}`);
  // ⛔ SIMILARITIES on the wire now, and there is nothing left to convert. These were medoid
  // DISTANCES and this line turned each into `1 − d`; they are the same set-to-set numbers the
  // locus card shows, so a surviving subtraction would print the complement of a real similarity —
  // 0.98 as 0.02 — on the one list whose whole job is to be judged gene by gene.
  const within = row.esm_within_similarity;
  const nearest = row.esm_nearest_similarity;
  if (within !== null && nearest !== null) bits.push(`ESM ${within.toFixed(2)}/${nearest.toFixed(2)}`);
  return bits.join(" · ");
}
</script>

<template>
  <details v-if="esmOnlyCount > 0" class="foot-detail" @toggle="load">
    <summary>Loci grouped by neither MMseqs nor a shared Pfam architecture ({{ esmOnlyCount.toLocaleString() }})</summary>
    <p class="muted">
      MMseqs groups none of these at any sequence identity, and no Pfam domain set is carried by more than
      95% of any one's genes. The locus card quotes each one's ESM similarity to its own medoid instead. Each
      is one click away.
    </p>
    <p v-if="failed !== null" class="pop-error" role="alert">The list did not load — {{ failed.detail }}.</p>
    <p v-else-if="lists === null" class="muted">Loading…</p>
    <div v-else class="goto-list">
      <button
        v-for="row in lists.esm_only"
        :key="row.label"
        type="button"
        class="goto-row"
        @click="emit('go', row.label)"
      >
        <b>{{ row.display_name }}</b> <span class="lid">&middot;{{ row.label }}</span>
        <span class="goto-ev">{{ evidence(row) }}</span>
      </button>
    </div>
  </details>
  <details v-if="notSharedCount > 0" class="foot-detail" @toggle="load">
    <summary>Loci where at most 95% of genes share one Pfam architecture ({{ notSharedCount.toLocaleString() }})</summary>
    <p class="muted">
      Loci with more than one gene and at least one Pfam domain, where the commonest domain set — order and
      repeat count ignored — is carried by 95% of the genes or fewer. Genes with no Pfam domain count in the
      denominator. Each is one click away.
    </p>
    <p v-if="failed !== null" class="pop-error" role="alert">The list did not load — {{ failed.detail }}.</p>
    <p v-else-if="lists === null" class="muted">Loading…</p>
    <div v-else class="goto-list">
      <button
        v-for="row in lists.pfam_architecture_not_shared"
        :key="row.label"
        type="button"
        class="goto-row"
        @click="emit('go', row.label)"
      >
        <b>{{ row.display_name }}</b> <span class="lid">&middot;{{ row.label }}</span>
        <span class="goto-ev">{{ evidence(row) }}</span>
      </button>
    </div>
  </details>
</template>
