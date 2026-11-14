<script setup lang="ts">
/**
 * The structural locus browser: the breadcrumb, the track, its two popovers, and the arrangement
 * switcher anchored under the focal gene.
 *
 * ⛔ **Three loading states, drawn three ways** (`locusNavigationStore`): `pending` has nothing to
 * show and draws a placeholder in the track's own place, so the page does not jump; `refreshing`
 * keeps the previous track up and dims it only after 150 ms; `failed` says so in words BESIDE the
 * last good track, never by clearing it — an empty track and an unanswered request must never look
 * alike.
 *
 * ⚠ **A popover that outlives what it points at is worse than no popover** (`app.js::placePop`). It
 * is measured from the live position of the bar that opened it, on open, on every track scroll and on
 * resize — the track scrolls sideways under a panel that does not — and it closes on Escape and on
 * any click outside it.
 */
import { storeToRefs } from "pinia";
import { computed, nextTick, onBeforeUnmount, onMounted, ref, watch } from "vue";

import type { CatalogueKey, NeighbourDisplayRow } from "@/api/types";
import AnchorGenomeBox from "@/components/anchor/AnchorGenomeBox.vue";
import OwnGenomesBox from "@/components/anchor/OwnGenomesBox.vue";
import LocusTrail from "@/components/navigation/LocusTrail.vue";
import ArrangementPopover from "@/components/popover/ArrangementPopover.vue";
import OffsetPopover from "@/components/popover/OffsetPopover.vue";
import ArrangementFan from "@/components/track/ArrangementFan.vue";
import ArrangementSwitcher from "@/components/track/ArrangementSwitcher.vue";
import BandKey from "@/components/track/BandKey.vue";
import GeneTrack from "@/components/track/GeneTrack.vue";
import { sharePercent } from "@/lib/formatting";
import { placePopover } from "@/lib/popoverPlacement";
import type { WalkDirection } from "@/lib/walkDirection";
import { useAnchorGenomeStore } from "@/stores/anchorGenomeStore";
import { useArrangementBrowserStore } from "@/stores/arrangementBrowserStore";
import { useLocusNavigationStore } from "@/stores/locusNavigationStore";
import { useTrackDisplayStore } from "@/stores/trackDisplayStore";

const props = defineProps<{
  /**
   * ⛔ TWO keys, because this panel feeds one child that ADDRESSES and one that IDENTIFIES.
   *
   * `catalogueKey` reaches `AnchorGenomeBox`, which fetches from the server. `speciesKey` reaches
   * `OwnGenomesBox`, whose store files answers under the organism — hand that one a catalogue key
   * and every genome modelled on this very page is reported as "modelled elsewhere".
   */
  catalogueKey: CatalogueKey | null;
  speciesKey: string | null;
  collectionGenomeCount: number;
  /** Where to send a reader whose address named no locus — offered, never taken for them. */
  landingLocus: string | null;
}>();

const navigation = useLocusNavigationStore();
const track = useTrackDisplayStore();
const browser = useArrangementBrowserStore();
const { sampleId: anchorSampleId, kind: anchorKind } = storeToRefs(useAnchorGenomeStore());
const { view, drawable, trail, displayNames, walkDirection } = storeToRefs(navigation);

/**
 * ⭐ The locus as a whole, at the apex of the neighbourhood fan: "n genes · x % of genomes" (David,
 * 2026-11-14). The options below divide the genomes that CARRY the locus; at 71,435 genomes a reader
 * needs the whole before the parts. ⚠ The share is of the collection, the same denominator the gene
 * card's tile uses, so the two can never disagree.
 */
const apexLabel = computed(() => {
  const locus = drawable.value?.locus;
  if (!locus || props.collectionGenomeCount <= 0) return null;
  const share = sharePercent(locus.genome_count / props.collectionGenomeCount);
  return `${locus.gene_count.toLocaleString()} genes · ${share} of genomes`;
});
const {
  selectedArrangementIndex,
  drawnArrangement,
  displayMirror,
  openPopoverSlot,
  openPopoverHeading,
  openPopoverMarginal,
  openPopoverDrawnLocus,
  isFocalPopoverOpen,
} = storeToRefs(track);

const isDimmed = computed(() => view.value.status === "refreshing" && view.value.isDimmed);
const failure = computed(() => (view.value.status === "failed" ? view.value.failure : null));

/*
 * ⚠ `topNeighbourCount` is DERIVED, and exactly: the export keeps the top N occupants per position,
 * so wherever something was left out (`observed_not_listed > 0`) the list that remains IS N long.
 * The published page read N from `meta.top_neighbours`, an export setting the database does not
 * hold; where nothing was cut there is no tail sentence to put N in.
 */
const neighboursByLabel = computed(() => {
  const index = new Map<string, NeighbourDisplayRow>();
  for (const row of drawable.value?.neighbour_display_rows ?? []) index.set(row.label, row);
  return index;
});

// ── the popover: placed from the live bar, dismissed from outside ─────────────────────────────
const popover = ref<HTMLDivElement | null>(null);
const frame = ref<HTMLDivElement | null>(null);
const popoverStyle = ref<Record<string, string>>({});
const isPopoverOpen = computed(() => openPopoverSlot.value !== null || isFocalPopoverOpen.value);

// ── the fan: redrawn from live positions, like the popover ────────────────────────────────────
const arrangements = ref<HTMLDivElement | null>(null);
const fan = ref<InstanceType<typeof ArrangementFan> | null>(null);

function redrawFan(): void {
  fan.value?.redraw(arrangements.value, frame.value);
}

watch([drawable, selectedArrangementIndex, walkDirection], async () => {
  await nextTick();
  // Twice: once for the new row, once after the track has centred its focal gene under it.
  redrawFan();
  requestAnimationFrame(redrawFan);
});

function onTrackScroll(): void {
  placeOpenPopover();
  redrawFan();
}

function onResize(): void {
  placeOpenPopover();
  redrawFan();
}

function placeOpenPopover(): void {
  const element = popover.value;
  const host = element?.offsetParent as HTMLElement | null | undefined;
  if (!element || !host || !isPopoverOpen.value) return;
  const anchor =
    frame.value?.querySelector(".slot.sel .marg-hit") ?? frame.value?.querySelector(".slot.sel .block");
  if (!anchor) {
    track.closePopover();
    return;
  }
  const position = placePopover(
    anchor.getBoundingClientRect(),
    { ...host.getBoundingClientRect().toJSON(), clientWidth: host.clientWidth },
    element.offsetWidth,
  );
  popoverStyle.value = { left: `${position.leftPx}px`, top: `${position.topPx}px` };
}

watch([openPopoverSlot, isFocalPopoverOpen, drawable], async () => {
  await nextTick();
  placeOpenPopover();
});

function onDocumentClick(event: MouseEvent): void {
  if (!isPopoverOpen.value || popover.value === null) return;
  // ⚠ The PATH, not `contains(target)`: a click on "Try again" re-renders the card before this
  // handler runs, removing the button from the DOM — `contains` then says "outside" and the card
  // closed on the reader mid-retry. The path is fixed at dispatch time.
  if (event.composedPath().includes(popover.value)) return;
  track.closePopover();
}

function onDocumentKeydown(event: KeyboardEvent): void {
  if (event.key === "Escape" && isPopoverOpen.value) track.closePopover();
}

let scroller: Element | null = null;
onMounted(() => {
  document.addEventListener("click", onDocumentClick);
  document.addEventListener("keydown", onDocumentKeydown);
  window.addEventListener("resize", onResize);
});
// The scroller is rendered by `GeneTrack` and replaced when the track first appears, so it is found
// after each draw rather than once.
watch(
  drawable,
  async () => {
    await nextTick();
    const next = frame.value?.querySelector(".track-scroll") ?? null;
    if (next === scroller) return;
    scroller?.removeEventListener("scroll", onTrackScroll);
    scroller = next;
    scroller?.addEventListener("scroll", onTrackScroll);
  },
  { immediate: true },
);
onBeforeUnmount(() => {
  document.removeEventListener("click", onDocumentClick);
  document.removeEventListener("keydown", onDocumentKeydown);
  window.removeEventListener("resize", onResize);
  scroller?.removeEventListener("scroll", onTrackScroll);
});

function walk(locus: string, direction: WalkDirection): void {
  track.closePopover();
  void navigation.navigateTo(locus, direction);
}

function jump(locus: string): void {
  void navigation.navigateTo(locus);
}
</script>

<template>
  <!-- ⭐ `is-projected` redefines `--anchor` for the whole panel, so every ⚓ mark, band and
       popover rule takes the placed colour without knowing this state exists. -->
  <section class="track-panel" :class="{ 'is-projected': anchorKind === 'projected' }">
    <div class="wrap">
      <!-- Heading and breadcrumb share ONE row (David, 2026-08-25): stacked they spent two lines on
           ~20 px of text. `.trail:not(:empty)` draws the divider, so it appears with the trail. -->
      <div class="panel-bar">
        <h2 class="panel-head">Structural locus browser</h2>
        <LocusTrail :trail="trail" :display-names="displayNames" @go="jump" />
      </div>

      <div ref="frame" class="track-frame">
        <!-- ⚠ Idle is NOT loading: nothing was asked for, and "Loading…" would wait forever. -->
        <div v-if="view.status === 'pending'" class="track-scroll track-waiting">
          <p class="muted">Loading the locus…</p>
        </div>
        <div v-else-if="view.status === 'idle'" class="track-scroll track-waiting">
          <p class="muted">Search a gene name above to open a locus.</p>
        </div>
        <GeneTrack
          v-else-if="drawable !== null"
          :detail="drawable"
          :collection-genome-count="props.collectionGenomeCount"
          :is-dimmed="isDimmed"
        />
      </div>

      <!-- ⛔ A failure is SAID, beside whatever is still drawn — never an empty track.
           ⚠ A missing locus is named in the reader's terms, not the server's: its 404 sentence is written
           for whoever is debugging a label, not for someone who followed an old link. And the published
           page's answer — silently opening the landing locus instead — is not repeated: a reader who
           asked for one locus and got another, with nothing said, has been told something false. The
           landing locus is OFFERED instead. -->
      <div v-if="failure !== null" class="pop-error track-error" role="alert">
        <span v-if="failure.kind === 'not_found'">
          There is no locus “{{ navigation.route?.label }}” in this catalogue.
          <template v-if="drawable !== null">The track below is still the previous locus.</template>
        </span>
        <span v-else>
          The locus did not load — {{ failure.detail.replace(/\.$/, "") }}.
          <template v-if="drawable !== null">The track below is still the previous locus.</template>
        </span>
        <button v-if="failure.kind !== 'not_found'" type="button" class="pop-retry" @click="navigation.retry()">
          Try again
        </button>
        <button
          v-else-if="drawable === null && landingLocus !== null"
          type="button"
          class="pop-retry"
          @click="jump(landingLocus)"
        >
          Open this catalogue's starting locus
        </button>
      </div>

      <div
        id="pop"
        ref="popover"
        class="pop"
        role="dialog"
        aria-label="Position detail"
        :hidden="!isPopoverOpen || drawable === null"
        :style="popoverStyle"
      >
        <template v-if="drawable !== null">
          <OffsetPopover
            v-if="openPopoverSlot !== null && openPopoverMarginal !== null && openPopoverHeading !== null"
            :marginal="openPopoverMarginal"
            :labelled-offset="openPopoverHeading.labelledOffset"
            :recorded-offset="openPopoverHeading.recordedOffset"
            :display-mirror="displayMirror"
            :drawn-locus="openPopoverDrawnLocus"
            :focal-gene-count="drawable.locus.gene_count"
            :neighbours-by-label="neighboursByLabel"
            :top-neighbour-count="openPopoverMarginal.observed_not_listed > 0 ? openPopoverMarginal.occupants.length : null"
            @walk="walk"
            @close="track.closePopover()"
          />
          <ArrangementPopover
            v-else-if="isFocalPopoverOpen"
            :locus="drawable.locus"
            :arrangements="browser.rows"
            :total="browser.total"
            :selected-rank="drawnArrangement?.rank ?? null"
            :anchor-ranks="drawable.anchor.arrangement_ranks"
            :members-without-a-neighbourhood="drawable.arrangements.members_without_a_neighbourhood"
            :walk-direction="walkDirection"
            :arrangements-not-shown="browser.arrangementsNotShown"
            :load-status="browser.status"
            :load-failure-detail="browser.lastFailure?.detail ?? null"
            @load-more="browser.loadMore()"
            @close="track.closePopover()"
          />
        </template>
      </div>

      <!-- The switcher, ANCHORED UNDER THE FOCAL GENE (David, 2026-08-23): the neighbourhoods are
           choices for that one locus. The anchor box sits in its left gutter and the channel key in
           its right, so the options stay centred under the focal gene whichever is occupied. -->
      <div v-if="drawable !== null" ref="arrangements" class="arr-wrap">
        <ArrangementFan ref="fan" :apex-label="apexLabel" />
        <!-- ⚠ ONE grid cell, two boxes: the anchor box has to stay centred on the option row
             (David, 2026-08-25), so the grid places this wrapper and the own-genomes box hangs out
             of flow beneath it — `.arr-gutter`, in the ADDED section of app.css. -->
        <div class="arr-gutter">
          <!-- ⭐ One header over both boxes (David, 2026-09-23). They are two answers to the same
               question — which genome the neighbourhood is drawn for — and naming that question once
               is what makes them read as a choice rather than as two unrelated controls.
               ⚠ `.legend` is the class "CHOOSE A MICROBE TO NAVIGATE:" already uses, so the two
               headers are the same kind of thing on the page rather than merely similar. -->
          <span class="legend arr-gutter-head">Anchor to:</span>
          <AnchorGenomeBox :catalogue-key="catalogueKey" placement="track" />
          <OwnGenomesBox :species-key="speciesKey" />
        </div>
        <ArrangementSwitcher
          :locus="drawable.locus"
          :arrangements="drawable.arrangements.listed"
          :total="drawable.arrangements.total"
          :selected-index="selectedArrangementIndex"
          :anchor-ranks="drawable.anchor.arrangement_ranks"
          :anchor-genome-name="anchorSampleId"
          :membership-is-complete="drawable.arrangements.membership_is_complete"
          :projected-copies="drawable.anchor.kind === 'projected' ? drawable.anchor.projected_copies : null"
          :walk-direction="walkDirection"
          @select="track.selectArrangement($event)"
        />
        <BandKey />
      </div>
    </div>
  </section>
</template>
