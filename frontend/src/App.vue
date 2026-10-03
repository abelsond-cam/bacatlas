<script setup lang="ts">
import { asCatalogueKey } from "@/api/types";
import type { CatalogueKey } from "@/api/types";
/**
 * The page: header and search, the species strip, the track, four views of one locus, the footer.
 *
 * ⭐ **The URL carries two things and no more.** The species is a query parameter (`?species=kp`),
 * because switching species is a full navigation to a different catalogue — the published site made
 * it a different page. The locus is the published page's own HASH (`#2811`, `#2811r`), so the
 * browser's Back walks the trail and a saved link keeps working. The anchor, the tab, the popover
 * and the map representation are not places, and stay out (`lib/locusHashRoute.ts`).
 *
 * ⚠ **Hash ⇄ route sync goes one way at a time.** A click changes the route, and the route writes
 * the hash (a history entry, as `app.js::go` made one). Back changes the hash, and the hash asks the
 * navigation store — which, finding the route already there, does nothing. Comparing DECODED forms
 * is what stops the two echoing each other forever over a percent-encoded label.
 */
import { storeToRefs } from "pinia";
import { computed, onBeforeUnmount, onMounted, ref, watch } from "vue";

import HoverTip from "@/components/layout/HoverTip.vue";
import OwnGenomesDialog from "@/components/ownGenomes/OwnGenomesDialog.vue";
import SiteFooter from "@/components/layout/SiteFooter.vue";
import SiteHeader from "@/components/layout/SiteHeader.vue";
import SpeciesStrip from "@/components/layout/SpeciesStrip.vue";
import TrackPanel from "@/components/layout/TrackPanel.vue";
import ViewTabs from "@/components/layout/ViewTabs.vue";
import FunctionView from "@/components/views/FunctionView.vue";
import LocusEvidenceView from "@/components/views/LocusEvidenceView.vue";
import NavigatingView from "@/components/views/NavigatingView.vue";
import SequenceView from "@/components/views/SequenceView.vue";
import { formatLocusHash } from "@/lib/locusHashRoute";
import { FORWARD } from "@/lib/walkDirection";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { useOwnGenomesStore } from "@/stores/ownGenomesStore";
import { useSpeciesCatalogueStore } from "@/stores/speciesCatalogueStore";
import { useViewTabStore } from "@/stores/viewTabStore";

const species = useSpeciesCatalogueStore();
const own = useOwnGenomesStore();
const navigation = useLocusNavigationStore();
const tabs = useViewTabStore();
const { speciesList, catalogue, current, publishedSpecies, catalogueKey, cataloguesForCurrentSpecies } =
  storeToRefs(species);
const { drawable, hash } = storeToRefs(navigation);
const { view } = storeToRefs(tabs);

const SPECIES_PARAMETER = "species";
/**
 * ⭐ The parameter that PINS a clustering. `?catalogue=ecoli-nuna5` names one exactly;
 * `?species=ecoli` follows whatever that species currently serves.
 *
 * ⚠ Both are read, and `catalogue` wins — but a followed `?species=` is REWRITTEN to the key it
 * resolved to, so a link copied from the address bar always names the clustering the sender was
 * looking at. Without that, `?species=ecoli#2811` silently opens a DIFFERENT GENE the day the
 * default moves, because locus labels are model-private — and it renders perfectly (David,
 * 2026-09-26).
 */
const CATALOGUE_PARAMETER = "catalogue";

function decoded(text: string): string {
  try {
    return decodeURIComponent(text.replace(/^#/, ""));
  } catch {
    return text.replace(/^#/, "");
  }
}

/** Jump to a locus from search, a crumb, a chip or the footer — no direction of travel, so forward. */
function go(locusLabel: string): void {
  void navigation.navigateTo(locusLabel);
}

/** Open whatever the address bar names, or the catalogue's landing locus if it names nothing. */
async function openFromAddress(): Promise<void> {
  const opened = await navigation.openHash(window.location.hash);
  const landing = current.value?.landing_locus ?? null;
  if (opened || landing === null) return;
  // ⛔ REPLACE the bare entry rather than pushing a second one: with a push, the reader's first Back
  // returned to the hash-less address, which re-opened the same landing locus — a Back that did
  // nothing, on every visit to the home URL.
  const url = new URL(window.location.href);
  url.hash = formatLocusHash({ label: landing, direction: FORWARD });
  window.history.replaceState(window.history.state, "", url);
  await navigation.navigateTo(landing);
}

function onHashChange(): void {
  if (current.value === null) return;
  if (decoded(window.location.hash) === decoded(hash.value)) return;
  void openFromAddress();
}

watch(hash, (next) => {
  if (!next || decoded(next) === decoded(window.location.hash)) return;
  window.location.hash = next;
});

/**
 * ⭐ A species switch LEAVES this catalogue: a new address, a fresh page, nothing carried across —
 * exactly as the published picker did. The trail, the cache and the anchor all belong to one
 * catalogue, and a reload is the one reset that cannot forget any of them.
 */
function selectCatalogue(nextKey: CatalogueKey): void {
  if (nextKey === catalogueKey.value) return;
  const url = new URL(window.location.href);
  url.searchParams.set(CATALOGUE_PARAMETER, nextKey);
  // ⚠ The bare species parameter is DROPPED, so no address can hold two answers that disagree.
  url.searchParams.delete(SPECIES_PARAMETER);
  // ⛔ The hash goes too, and for a model switch that is not housekeeping. Locus labels are
  // MODEL-PRIVATE: `#2811` is a different gene in nuna5 than in nuna4, so carrying it across would
  // land the reader on a real locus that is not the one they were reading, with nothing saying so.
  // The new catalogue's own landing locus opens instead.
  url.hash = "";
  window.location.assign(url.toString());
}

/**
 * ⭐ The placed genomes are fetched as soon as the species list is known, not when the dialog opens.
 *
 * ⚠ The gutter box has to KNOW whether anything is placed before a reader clicks it: it opens on the
 * summary when there is something to show and on the accounts page when there is not, and it draws
 * its placed state from the same list. Fetching only inside the dialog left the box unable to answer
 * either question until after the reader had already been sent to the wrong screen. Two small
 * responses, once per page load, against a cache that is immutable for a build.
 */
watch(
  publishedSpecies,
  (entries) => {
    if (entries.length === 0) return;
    void own.loadPlacedGenomes(
      entries
        .filter((entry) => entry.catalogue_key !== null)
        .map((entry) => ({
          key: entry.key,
          catalogueKey: entry.catalogue_key as CatalogueKey,
          scientificName: entry.scientific_name,
        })),
    );
  },
  { immediate: true },
);

const pageTitle = computed(() =>
  current.value ? `${current.value.species.scientific_name} · BacAtlas` : "BacAtlas",
);
watch(pageTitle, (title) => {
  document.title = title;
});

/**
 * A `?species=` naming nothing this server publishes. ⛔ It is SAID, never swapped for another
 * catalogue: locus labels are small integers that exist in both (locus 1098 is `yohK` in ecoli and
 * `fimA` in kp), so a silent substitution would open the link's locus number in the wrong organism,
 * under a URL still naming the one the reader asked for.
 */
const unknownSpecies = ref<string | null>(null);

onMounted(async () => {
  window.addEventListener("hashchange", onHashChange);
  await Promise.all([species.loadSpeciesList(), species.loadCatalogues()]);
  const address = new URL(window.location.href).searchParams;
  const pinned = address.get(CATALOGUE_PARAMETER);
  let chosen: CatalogueKey | undefined;

  if (pinned !== null) {
    // ⛔ NOT validated against `GET /catalogues`. A catalogue can be loaded but not offered — that
    // is how a model is staged for review before anyone is shown it — so gating on the menu would
    // make a deliberate link unopenable. Send it and let the server answer; a typo lands as a named
    // catalogue failure below, which is honest, rather than as "no such species".
    chosen = asCatalogueKey(pinned);
  } else {
    const requested = address.get(SPECIES_PARAMETER);
    const entry =
      requested === null
        ? publishedSpecies.value[0]
        : // Case is forgiven — `KP` is plainly `kp`.
          publishedSpecies.value.find((candidate) => candidate.key === requested.toLowerCase());
    if (entry === undefined) {
      if (requested !== null && speciesList.value.status === "ready") unknownSpecies.value = requested;
      return;
    }
    if (entry.catalogue_key === null) return;
    chosen = entry.catalogue_key;
    // ⭐ PIN THE ADDRESS to what the species actually resolved to (David, 2026-09-26). A link
    // copied from here then names a clustering rather than a default that can move underneath it.
    const url = new URL(window.location.href);
    url.searchParams.set(CATALOGUE_PARAMETER, chosen);
    url.searchParams.delete(SPECIES_PARAMETER);
    window.history.replaceState(window.history.state, "", url);
  }

  if (chosen === undefined) return;
  await species.selectCatalogue(chosen);
  if (current.value !== null) await openFromAddress();
});

onBeforeUnmount(() => window.removeEventListener("hashchange", onHashChange));
</script>

<template>
  <HoverTip />
  <!-- ⚠ Mounted ONCE, here, not in the gutter: a <dialog> opened with `showModal()` belongs to the
       page, and the box that opens it is redrawn on every walk. -->
  <!-- ⛔ The SPECIES key, not the catalogue: this dialog files answers by organism. Handed
       `ecoli-nuna5` it would miss every lookup and report genomes modelled on this very page
       as "modelled elsewhere". -->
  <OwnGenomesDialog :species-key="current?.species.key ?? null" :species="publishedSpecies" />
  <SiteHeader
    :catalogue-key="catalogueKey"
    :collection-genome-count="current?.pangenome.genome_count ?? null"
    @go="go"
  />
  <SpeciesStrip
    :species="publishedSpecies"
    :species-key="current?.species.key ?? null"
    :catalogue-key="catalogueKey"
    :catalogues="cataloguesForCurrentSpecies"
    :catalogue="current"
    @select-catalogue="selectCatalogue"
  />

  <!-- ⛔ Each failure SAYS which thing failed; none of them renders as an empty catalogue. -->
  <section v-if="speciesList.status === 'failed'" class="wrap load-failure">
    <p class="pop-error" role="alert">The list of species did not load — {{ speciesList.failure.detail }}.</p>
  </section>
  <section v-else-if="unknownSpecies !== null" class="wrap load-failure">
    <p class="pop-error" role="alert">
      No species “{{ unknownSpecies }}” is published on this server — choose one above.
    </p>
  </section>
  <section v-else-if="speciesList.status === 'ready' && publishedSpecies.length === 0" class="wrap load-failure">
    <p class="lede">No species is published on this server yet.</p>
  </section>
  <section v-else-if="catalogue.status === 'failed'" class="wrap load-failure">
    <p class="pop-error" role="alert">This catalogue did not load — {{ catalogue.failure.detail }}.</p>
  </section>

  <template v-if="current !== null">
    <TrackPanel
      :catalogue-key="catalogueKey"
      :species-key="current.species.key"
      :collection-genome-count="current.pangenome.genome_count"
      :landing-locus="current.landing_locus"
    />

    <!-- Four views of one locus. They sit OUTSIDE the track panel because the track is common to all
         four — you navigate on it whichever view you are reading. -->
    <main class="views">
      <div class="wrap">
        <ViewTabs :view="view" @show="tabs.showView($event)" />
      </div>
      <!-- In the strip's order (`VIEW_TABS`): the three views of this locus, then the documentation. -->
      <section id="view-node" class="view-panel" role="tabpanel" aria-label="Syntelog Loci" :hidden="view !== 'locus'">
        <div class="wrap">
          <LocusEvidenceView v-if="drawable !== null" :detail="drawable" :catalogue="current" />
        </div>
      </section>
      <section id="view-sequence" class="view-panel" role="tabpanel" aria-label="Sequence" :hidden="view !== 'sequence'">
        <div class="wrap">
          <SequenceView
            v-if="drawable !== null"
            :detail="drawable"
            :catalogue-key="current.pangenome.catalogue_key"
          />
        </div>
      </section>
      <section id="view-function" class="view-panel" role="tabpanel" aria-label="Function" :hidden="view !== 'function'">
        <div class="wrap">
          <FunctionView v-if="drawable !== null" :detail="drawable" />
        </div>
      </section>
      <section
        id="view-navigating"
        class="view-panel"
        role="tabpanel"
        aria-label="Navigating BacAtlas"
        :hidden="view !== 'navigating'"
      >
        <div class="wrap">
          <NavigatingView :examples="current.example_locus_rows" @go="go" />
        </div>
      </section>
    </main>

    <SiteFooter :catalogue="current" @go="go" />
  </template>
</template>
