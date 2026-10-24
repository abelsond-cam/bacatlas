/**
 * Which species the page is showing, and everything about it that arrives before a locus.
 *
 * ⭐ **A species switch is a navigation, not a filter** (`app.js`, the species picker). One page holds
 * one catalogue: the loci, the anchor genomes and every cache key belong to it, and nothing carries
 * across. So choosing a species resets the navigation store rather than re-filtering it.
 *
 * ⚠ **Three loading states here too, never conflated** — the same rule as the locus view. A species
 * list that failed to load and a deployment with no published species are different claims, and the
 * picker must not make the same face at both.
 */

import { defineStore } from "pinia";
import { computed, ref, shallowRef } from "vue";

import { fetchCatalogues, fetchSpeciesCatalogue, fetchSpeciesList } from "@/api/client";
import type { Failure } from "@/api/result";
import type {
  CatalogueEntry,
  CatalogueKey,
  SpeciesCatalogueResponse,
  SpeciesListResponse,
} from "@/api/types";

import { useAnchorGenomeStore } from "./anchorGenomeStore";
import { useLocusNavigationStore } from "./locusNavigationStore";

export type CatalogueLoad<T> =
  | { readonly status: "idle" }
  | { readonly status: "pending" }
  | { readonly status: "ready"; readonly value: T }
  | { readonly status: "failed"; readonly failure: Failure };

export type SpeciesEntry = SpeciesListResponse["species"][number];

export const useSpeciesCatalogueStore = defineStore("speciesCatalogue", () => {
  const navigation = useLocusNavigationStore();
  const anchor = useAnchorGenomeStore();

  const speciesList = shallowRef<CatalogueLoad<readonly SpeciesEntry[]>>({ status: "idle" });
  const catalogues = shallowRef<CatalogueLoad<readonly CatalogueEntry[]>>({ status: "idle" });
  const catalogue = shallowRef<CatalogueLoad<SpeciesCatalogueResponse>>({ status: "idle" });
  /**
   * ⛔ The CATALOGUE being shown — `ecoli-nuna5`, not `ecoli`.
   *
   * It was the species key while a species had one catalogue. Everything that ADDRESSES the
   * server keys on this; anything that asks "which organism is this" must use
   * `current.species.key` instead, because the two stopped being the same string.
   */
  const catalogueKey = ref<CatalogueKey | null>(null);

  /** The loaded catalogue, or `null` — for components that render nothing until it exists. */
  const current = computed<SpeciesCatalogueResponse | null>(() =>
    catalogue.value.status === "ready" ? catalogue.value.value : null,
  );

  /** Every catalogue on offer, or an empty list until they arrive. */
  const offeredCatalogues = computed<readonly CatalogueEntry[]>(() =>
    catalogues.value.status === "ready" ? catalogues.value.value : [],
  );

  /**
   * The catalogues of one species — what the model picker offers beside it.
   *
   * ⚠ Empty until the list loads, and of length one for a species with a single model. The strip
   * hides itself in that case rather than showing a control with nothing to choose.
   */
  const cataloguesForCurrentSpecies = computed<readonly CatalogueEntry[]>(() => {
    const species = current.value?.species.key;
    return species ? offeredCatalogues.value.filter((entry) => entry.species.key === species) : [];
  });

  async function loadCatalogues(): Promise<void> {
    catalogues.value = { status: "pending" };
    const result = await fetchCatalogues();
    catalogues.value = result.ok
      ? { status: "ready", value: result.value.catalogues }
      : { status: "failed", failure: result };
  }

  /** Only the species that serve something are offered as a destination. */
  const publishedSpecies = computed<readonly SpeciesEntry[]>(() =>
    speciesList.value.status === "ready"
      ? speciesList.value.value.filter((entry) => entry.published)
      : [],
  );

  async function loadSpeciesList(): Promise<void> {
    speciesList.value = { status: "pending" };
    const result = await fetchSpeciesList();
    speciesList.value = result.ok
      ? { status: "ready", value: result.value.species }
      : { status: "failed", failure: result };
  }

  /**
   * Show one species. Resets the navigation (and with it the trail, the cache and the anchor)
   * BEFORE the catalogue arrives, so nothing from the previous species can be drawn under the new
   * one's name while the request is in flight.
   */
  async function selectCatalogue(nextKey: CatalogueKey): Promise<void> {
    // ⛔ The guard keys on the CATALOGUE, not the species. Keyed on the species it made switching
    // model within one species a silent no-op — the picker would move and the page would not.
    if (catalogueKey.value === nextKey && catalogue.value.status !== "failed") return;
    catalogueKey.value = nextKey;
    navigation.setCatalogue(nextKey);
    catalogue.value = { status: "pending" };
    const result = await fetchSpeciesCatalogue(nextKey);
    // ⚠ The reader may have picked another catalogue while this was in flight.
    if (catalogueKey.value !== nextKey) return;
    if (!result.ok) {
      catalogue.value = { status: "failed", failure: result };
      anchor.setAvailability(false);
      return;
    }
    catalogue.value = { status: "ready", value: result.value };
    // A published pangenome always carries gene membership in the database, so an anchor can always
    // be answered — the published page withheld the box only on payloads that shipped no `arr.gid`.
    anchor.setAvailability(result.value.pangenome.genome_count > 0);
    // ⭐ The per-group census, so the picker knows which lists this catalogue can offer BEFORE the
    // reader opens it. ⚠ Read from the shell rather than counted from a listing: a count derived
    // from a filtered list could only ever say "nothing matched your query".
    anchor.setSourceGroupCounts(result.value.pangenome.genome_source_groups);
  }

  return {
    speciesList,
    catalogues,
    catalogue,
    catalogueKey,
    current,
    publishedSpecies,
    offeredCatalogues,
    cataloguesForCurrentSpecies,
    loadSpeciesList,
    loadCatalogues,
    selectCatalogue,
  };
});
