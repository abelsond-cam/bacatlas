<script setup lang="ts">
/**
 * The species picker, the COHORT picker beside it, and the pangenome census. (The mark moved up to
 * the instrument bar on 2026-09-22 — `SiteHeader.vue`.)
 *
 * ⛔ **THE SECOND CONTROL CHOOSES A COHORT, NOT A MODEL (David, 2026-10-27).** It names the genomes
 * — `Kp AMR — 100 genomes` — and the model is one free-text label the footer prints. The option
 * text is composed here from `cohort.label` and `cohort.genome_count`: the label NEVER carries the
 * count, because `genome_count` is passed, never inferred, and a label that spelled its own number
 * could contradict the roster it names.
 *
 * ⭐ **A switch is a navigation, not a filter** — this page holds one catalogue at a time, so
 * choosing another species, or another cohort of the same species, leaves this one entirely (the
 * published picker set `location.href`).
 *
 * ⛔ **Both controls emit the same event, because both choose the same thing.** A catalogue is one
 * species × one cohort, and the page is only ever showing one of them. Picking a species means
 * "that species' default catalogue"; picking a cohort means "this species, that cohort". Two emits
 * would have made the parent decide which control had spoken last, and that is a state machine
 * nobody needs.
 */
import { computed } from "vue";

import type { CatalogueEntry, CatalogueKey, SpeciesCatalogueResponse } from "@/api/types";
import PangenomeCensus from "@/components/census/PangenomeCensus.vue";
import type { SpeciesEntry } from "@/stores/speciesCatalogueStore";

const props = defineProps<{
  species: readonly SpeciesEntry[];
  /** ⚠ The organism, for the left-hand control — NOT the catalogue key. */
  speciesKey: string | null;
  /** The catalogue actually being shown, which is what the cohort control marks as selected. */
  catalogueKey: CatalogueKey | null;
  /** This species' catalogues. Length one for a species with a single cohort. */
  catalogues: readonly CatalogueEntry[];
  catalogue: SpeciesCatalogueResponse | null;
}>();

const emit = defineEmits<{ selectCatalogue: [catalogueKey: CatalogueKey] }>();

function onSpeciesChange(event: Event): void {
  const chosen = (event.target as HTMLSelectElement).value;
  const entry = props.species.find((candidate) => candidate.key === chosen);
  // ⚠ A species with nothing published has no catalogue to go to, so the control does not pretend.
  if (entry?.catalogue_key) emit("selectCatalogue", entry.catalogue_key);
}

function onCohortChange(event: Event): void {
  const chosen = (event.target as HTMLSelectElement).value;
  if (chosen) emit("selectCatalogue", chosen as CatalogueKey);
}

/**
 * `E. coli AMR — 100 genomes`. The label and the count are composed HERE and never stored together.
 *
 * ⛔ **No optional chaining, and its absence is the point.** The function this replaced read
 * `entry.model?.key ?? entry.key`, so a response with no model block rendered the catalogue key as
 * the option text — plausible output from missing data, which is the hardest kind of wrong to
 * notice. `cohort` is NOT NULL in the database and non-optional in `CatalogueEntry`, so a missing
 * one is a bug that should throw in the one place a human is looking.
 */
function cohortName(entry: CatalogueEntry): string {
  return `${entry.cohort.label} — ${entry.cohort.genome_count.toLocaleString()} genomes`;
}

/** The one cohort a single-cohort species has, or `null`. ⚠ A computed rather than `catalogues[0]`
 * in the template: `noUncheckedIndexedAccess` types that as possibly-undefined, and the honest
 * narrowing is this one, not a `!`. */
const soleCohort = computed<CatalogueEntry | null>(() =>
  props.catalogues.length === 1 ? (props.catalogues[0] ?? null) : null,
);
</script>

<template>
  <section class="species-strip">
    <div class="wrap species-in">
      <label class="species-pick">
        <span class="legend">Choose a microbe to navigate:</span>
        <select id="species" aria-label="Species" :value="speciesKey ?? ''" @change="onSpeciesChange">
          <option v-for="entry in species" :key="entry.key" :value="entry.key">
            {{ entry.scientific_name }}
          </option>
        </select>
      </label>
      <!-- ⛔ **SHOWN EVEN WHEN THERE IS ONE COHORT, which reverses the old rule and has to.** The
           control used to be hidden below two catalogues, on the sound argument that a choice of
           one is not a choice. But it is now the ONLY place the cohort is named, and the cohort is
           what the whole page is about — hidden, a reader of a single-cohort species would never
           learn which genomes they are looking at. One cohort renders as static text rather than a
           select, because a dropdown with one option still says "there is a choice here". -->
      <label v-if="catalogues.length > 1" class="species-pick cohort-pick">
        <span class="legend">Cohort:</span>
        <select id="cohort" aria-label="Genome cohort" :value="catalogueKey ?? ''" @change="onCohortChange">
          <option v-for="entry in catalogues" :key="entry.key" :value="entry.key">
            {{ cohortName(entry) }}{{ entry.is_default ? " (default)" : "" }}
          </option>
        </select>
      </label>
      <div v-else-if="soleCohort" class="species-pick cohort-pick">
        <span class="legend">Cohort:</span>
        <span class="cohort-sole">{{ cohortName(soleCohort) }}</span>
      </div>
      <div class="pangenome">
        <span class="legend">Pangenome</span>
        <PangenomeCensus v-if="catalogue !== null" :catalogue="catalogue" />
      </div>
    </div>
  </section>
</template>
