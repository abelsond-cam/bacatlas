<script setup lang="ts">
/**
 * The Function tab — COG, GO and EC/KEGG across every member gene, fetched when the tab opens.
 *
 * ⚠ It also carries the inferred-function card, whose donor is **walkable**: a suggestion is only as
 * good as the node it came from, so the reader can go and look at that node. The walk goes through
 * the same navigation store every other walkable block uses, so it lands in history identically.
 */
import { storeToRefs } from "pinia";

import type { LocusDetailResponse } from "@/api/types";
import FunctionTab from "@/components/function/FunctionTab.vue";
import { useFunctionBlockStore } from "@/stores/functionBlockStore";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";

defineProps<{ detail: LocusDetailResponse }>();

const functionBlock = useFunctionBlockStore();
const navigation = useLocusNavigationStore();
const { block, status, lastFailure } = storeToRefs(functionBlock);
</script>

<template>
  <FunctionTab
    :display-name="detail.locus.display_name"
    :locus-label="detail.locus.label"
    :block="block"
    :status="status"
    :failure-detail="lastFailure?.detail ?? null"
    @retry="functionBlock.load()"
    @walk="navigation.navigateTo($event)"
  />
</template>
