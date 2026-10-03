<script setup lang="ts">
/**
 * "Syntelog Loci" — the home view: what this locus is, and the evidence for it being one.
 *
 * Two columns, as the published card drew them (`app.js::renderCard`): the ARGUMENT on the left —
 * the headline and the sequence-diversity card — and the REFERENCE column beside it.
 *
 * ⛔ That column held two cards until 2026-09-24: the embedding geometry, and a neighbourhood map
 * immediately below the numbers it was a picture of. Both went together, because they were ONE
 * construction — the map was an MDS of six loci's MEDOIDS and each ring was that locus's members'
 * median distance to its own medoid, so a set-to-set similarity has nothing for it to draw. The five
 * nearest loci survive: they are now ranked by a median over gene pairs and are listed on the
 * similarity card itself, beside the row they are the evidence for.
 */
import { storeToRefs } from "pinia";
import { computed } from "vue";

import type { LocusDetailResponse, SpeciesCatalogueResponse } from "@/api/types";
import EmbeddingSimilarityCard from "@/components/locusCard/EmbeddingSimilarityCard.vue";
import LocusHeadline from "@/components/locusCard/LocusHeadline.vue";
import SequenceDiversityCard from "@/components/locusCard/SequenceDiversityCard.vue";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { useSimilarityViewStore } from "@/stores/similarityViewStore";

const props = defineProps<{
  detail: LocusDetailResponse;
  catalogue: SpeciesCatalogueResponse;
}>();

const navigation = useLocusNavigationStore();
const similarityView = useSimilarityViewStore();
const { view } = storeToRefs(similarityView);

const baselines = computed(() => props.catalogue.similarity_baselines);
/*
 * ⛔ `separationCount` lived here to give the headline's separation TILE its "of N loci" denominator,
 * reading Bacformer's measurable count because the tile read Bacformer. The tile moved into
 * `EmbeddingSimilarityCard` on 2026-10-03, and that card already derives each representation's OWN
 * measurable count from the same `baselines` — so the prop went with the tile rather than being
 * threaded through a component that no longer shows it.
 */
</script>

<template>
  <div class="work">
    <div>
      <div id="card">
        <LocusHeadline
          :locus="detail.locus"
          :collection-genome-count="catalogue.pangenome.genome_count"
        />
        <SequenceDiversityCard
          :locus="detail.locus"
          :families="detail.uniref50_families"
          :symbols="detail.annotations['gene_symbol'] ?? []"
          :pfam-reference="detail.pfam_reference"
          :listed-architecture-count="(detail.annotations['pfam_architecture'] ?? []).length"
        />
      </div>
    </div>
    <aside id="side" class="side">
      <EmbeddingSimilarityCard
        :similarity="detail.locus.similarity"
        :baselines="baselines"
        :neighbours="detail.neighbour_display_rows"
        :gene-count="detail.locus.gene_count"
        :view="view"
        @walk="navigation.navigateTo($event)"
        @select-view="similarityView.selectView($event)"
      />
    </aside>
  </div>
</template>
