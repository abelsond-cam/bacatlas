<script setup lang="ts">
/**
 * Anchor the page to one genome: the track draws ITS neighbourhood at every locus instead of the
 * commonest one, and the Sequence tab shows its bases.
 *
 * ⭐ **One anchor, MOUNTED TWICE** — the switcher's gutter and the top of the Sequence tab. On the
 * published page `setAnchor` walked every mounted box writing the value through, because *"the boxes
 * are two views of one state, not two states that are kept in step"*. Here that sentence is the
 * implementation: both instances read one store and there is nothing to keep in step.
 *
 * ⭐ **The search idiom, deliberately** — input → filter → a list of buttons, ArrowUp/Down/Enter,
 * `mousedown` prevented so blur cannot beat the click, a short close delay on blur. The cursor runs
 * from −1, the "most common arrangement" row, which is ALWAYS first and never hidden behind a ×: it
 * is the off switch, and an off switch a reader has to hunt for is one they conclude does not exist.
 *
 * ⛔ **REMOVED, not disabled, when the catalogue cannot answer an anchor** — the rule the published
 * page set: a control that cannot tell "this genome is not here" from "the page has no idea" must
 * not be offered, because the box would make the same face at both.
 *
 * ⚠ **The genomes are filtered on the server now.** The published page held all of them in
 * `meta.genomes`; at 80,000 genomes that is an array nobody should ship, so a list that was cut says
 * so and says how much it left out.
 *
 * ⭐ **TWO LISTS, because the cohort is two populations.** ~4,000 complete long-read NCBI genomes,
 * short enough to browse, and ~75,000 AllTheBacteria short-read assemblies, which need the search.
 * One combobox over both serves neither: the complete genomes are the ones a reader anchors
 * deliberately and would never find by typing, and they are 5 % of the list. So the dropdown opens
 * with a row of group tabs and the search narrows WITHIN the chosen group.
 *
 * ⛔ **A group with no members is not offered as a tab, and its absence is SAID rather than shown
 * as an empty list** — the same rule this control already applies to itself. On every catalogue
 * loaded today the complete-genome group is empty (all 280 probe genomes are BioSample-keyed ATB
 * assemblies, measured), so this is the ordinary case and not an edge one: the reader is told the
 * catalogue holds none, which is a different sentence from "none match what you typed" and is the
 * one that stops them typing.
 */
import { storeToRefs } from "pinia";
import { computed, ref, watch } from "vue";

import { fetchGenomes } from "@/api/client";
import type { Failure } from "@/api/result";
import type { CatalogueKey, GenomeListResponse, GenomeSourceGroup } from "@/api/types";
import {
  SOURCE_GROUP_DESCRIPTIONS,
  SOURCE_GROUP_EMPTY_SENTENCES,
  SOURCE_GROUP_LABELS,
  useAnchorGenomeStore,
} from "@/stores/anchorGenomeStore";

const props = defineProps<{
  catalogueKey: CatalogueKey | null;
  /** Which mounting this is — it only decides the element ids, never the state. */
  placement: "track" | "sequence";
}>();

/** How many genomes one answer lists before saying it was cut. */
const LIST_LIMIT = 200;
const CLOSE_AFTER_BLUR_MS = 120;

const anchor = useAnchorGenomeStore();
const {
  sampleId: anchoredSampleId,
  isAvailable,
  kind,
  sourceGroup,
  sourceGroupCounts,
  offeredSourceGroups,
  hasSourceGroupChoice,
} = storeToRefs(anchor);

/**
 * ⛔ **This box shows only a CATALOGUE anchor.** There is one anchor and both boxes read it, but a
 * genome placed on the model afterwards belongs to the other box — showing its accession here would
 * put it behind the ⚓, which means "one of the modelled genomes", about a genome that is not one.
 * So a projected anchor reads as no anchor here, and the box goes back to its grey empty state.
 */
const sampleId = computed(() => (kind.value === "catalogue" ? anchoredSampleId.value : null));
const isAnchored = computed(() => sampleId.value !== null);

const ids = computed(() =>
  props.placement === "track"
    ? { wrap: "anchor-wrap", box: "anchor-genome", results: "anchor-results" }
    : { wrap: "seq-anchor-wrap", box: "seq-anchor-genome", results: "seq-anchor-results" },
);

const box = ref<HTMLInputElement | null>(null);
const text = ref(sampleId.value ?? "");
const isOpen = ref(false);
const answer = ref<GenomeListResponse | null>(null);
const failed = ref<Failure | null>(null);
/** −1 is the "most common arrangement" row. */
const cursor = ref(-1);
let issued = 0;

// ⭐ Either box setting the anchor rewrites BOTH — they read the same store.
watch(sampleId, (next) => {
  text.value = next ?? "";
});

async function search(query: string): Promise<void> {
  const species = props.catalogueKey;
  if (species === null) return;
  const token = ++issued;
  // ⚠ The group is read at ISSUE time and sent with the request, so a late answer is matched
  // against the filter it was asked for — not against whichever tab is open when it lands.
  const asked = sourceGroup.value;
  const result = await fetchGenomes(species, query.trim(), {
    limit: LIST_LIMIT,
    sourceGroup: asked,
  });
  if (token !== issued) return;
  if (result.ok) {
    answer.value = result.value;
    failed.value = null;
  } else {
    answer.value = null;
    failed.value = result;
  }
}

/**
 * ⛔ **Switching tab re-queries and resets the cursor, and keeps the typed text.** Keeping the
 * query is the point of the tabs: a reader who typed `SAMEA22` and switches group is asking the
 * same question of the other population. Resetting the cursor is not cosmetic — it indexes into
 * `answer.genomes`, so a stale index can point past the end of the new list or, worse, at a
 * different genome than the one highlighted.
 */
function chooseSourceGroup(group: GenomeSourceGroup): void {
  if (group === sourceGroup.value) return;
  anchor.setSourceGroup(group);
  cursor.value = -1;
  answer.value = null;
  void search(text.value === sampleId.value ? "" : text.value);
}

/** ⚠ The count of what the catalogue HOLDS in this group, never of what matched — see the store. */
function heldInGroup(group: GenomeSourceGroup): number {
  return sourceGroupCounts.value[group] ?? 0;
}

/** The sentence shown when the chosen group has nothing to list, and WHY it has nothing. */
const emptyGroupReason = computed(() => {
  if (answer.value === null) return null;
  if (answer.value.genomes.length > 0) return null;
  const group = sourceGroup.value;
  // ⚠ Unfiltered — the pre-census behaviour, and the original sentence unchanged.
  if (group === null) return `No genome in this catalogue matches \u201c${answer.value.query}\u201d.`;
  // ⛔ The two cases the whole census exists to separate: the catalogue holds none of these at all
  // (true of the complete-genome list on every catalogue loaded today), against none matching what
  // was typed. A reader shown the second would keep typing at a list that can never answer.
  if (heldInGroup(group) === 0) return SOURCE_GROUP_EMPTY_SENTENCES[group];
  return `No genome in ${SOURCE_GROUP_LABELS[group]} matches \u201c${answer.value.query}\u201d.`;
});

/** The group's name, or nothing where the list is unfiltered. */
const sourceGroupLabel = computed(() =>
  sourceGroup.value === null ? null : SOURCE_GROUP_LABELS[sourceGroup.value],
);

function open(): void {
  isOpen.value = true;
  cursor.value = -1;
  // On focus the box holds the current anchor; listing only that one genome would hide the others,
  // so an unedited box lists everything.
  void search(text.value === sampleId.value ? "" : text.value);
}

function close(): void {
  isOpen.value = false;
  cursor.value = -1;
}

function choose(next: string | null): void {
  anchor.setAnchor(next);
  text.value = next ?? "";
  close();
  box.value?.blur();
}

function onInput(): void {
  isOpen.value = true;
  cursor.value = -1;
  void search(text.value);
}

function onKeydown(event: KeyboardEvent): void {
  if (event.key === "Escape") {
    close();
    box.value?.blur();
    return;
  }
  const hits = answer.value?.genomes ?? [];
  // ⛔ No hits, or an answer to an earlier query: every key below does nothing, as on the published
  // page (`if (!A.hits.length) return;`). Otherwise Enter on "No genome matches" hit the −1 row and
  // CLEARED the reader's anchor — the one thing they were not asking for.
  const answeredQuery = answer.value?.query ?? null;
  const typed = text.value === sampleId.value ? "" : text.value.trim();
  // ⛔ The answer must match the query AND the group it is being navigated in. Guarding the query
  // alone, Enter on a list still showing the previous group's genomes would anchor a genome from a
  // list the reader had already left — and the accession would look perfectly plausible.
  const answeredGroup = answer.value?.source_group ?? null;
  if (hits.length === 0 || answeredQuery !== typed || answeredGroup !== sourceGroup.value) return;
  if (event.key === "ArrowDown") {
    cursor.value = cursor.value + 1 >= hits.length ? -1 : cursor.value + 1;
    event.preventDefault();
  } else if (event.key === "ArrowUp") {
    cursor.value = cursor.value <= -1 ? hits.length - 1 : cursor.value - 1;
    event.preventDefault();
  } else if (event.key === "Enter") {
    choose(cursor.value < 0 ? null : (hits[cursor.value]?.sample_id ?? null));
    event.preventDefault();
  }
}

function onBlur(): void {
  setTimeout(() => {
    close();
    // A half-typed query that was never chosen is not an anchor; put back the one that is.
    text.value = sampleId.value ?? "";
  }, CLOSE_AFTER_BLUR_MS);
}
</script>

<template>
  <div v-if="isAvailable" :id="ids.wrap" class="arr-anchor" :class="{ on: isAnchored, 'seq-anchor': placement === 'sequence' }">
    <label class="arr-anchor-box" :for="ids.box">
      <input
        :id="ids.box"
        ref="box"
        v-model="text"
        type="search"
        autocomplete="off"
        spellcheck="false"
        placeholder="A genome in BacAtlas"
        :aria-label="placement === 'track' ? 'Anchor the neighbourhood to one genome' : 'Anchor the sequence to one genome'"
        @focus="open"
        @input="onInput"
        @keydown="onKeydown"
        @blur="onBlur"
      />
      <!-- ⚠ U+2693 U+FE0E — the variation selector asks for the TEXT anchor. The bare character can
           draw as a colour emoji, which ignores `color`, and the grey state would never be grey. -->
      <span class="arr-anchor-mark" aria-hidden="true">⚓︎</span>
      <span class="arr-anchor-caret" aria-hidden="true">▾</span>
    </label>
    <!-- ⚠ Clicks stop HERE: the document-level closer dismisses the popover on any outside click,
         and would take this list with it before a click on a row could land. -->
    <!-- ⛔ `role="listbox"` moved INWARD, onto the options alone. The group tabs are
         `role="tab"` controls, and a listbox may contain only options — tabs inside one are
         invalid ARIA, and a screen reader would announce two conflicting widgets as one. -->
    <div :id="ids.results" class="results anchor-results" @click.stop>
      <template v-if="isOpen">
        <div v-if="failed !== null" class="hit-none pop-error" role="alert">
          The genome list did not load — {{ failed.detail }}.
        </div>
        <template v-else-if="answer !== null">
          <!-- ⭐ The group tabs. ⚠ Rendered only where there is a CHOICE: one list needs no tabs,
               and a single tab that cannot be switched is a control that does nothing. -->
          <div v-if="hasSourceGroupChoice" class="anchor-groups" role="tablist">
            <button
              v-for="group in offeredSourceGroups"
              :key="group"
              type="button"
              class="anchor-group"
              :class="{ on: group === sourceGroup, empty: heldInGroup(group) === 0 }"
              role="tab"
              :aria-selected="group === sourceGroup ? 'true' : 'false'"
              :title="SOURCE_GROUP_DESCRIPTIONS[group]"
              @mousedown.prevent
              @click="chooseSourceGroup(group)"
            >
              <!-- ⛔ The name is its OWN element so it can ellipsis alone; as a bare text node it
                   took the whole button's `overflow: hidden` with it and clipped the count away. -->
              <span class="anchor-group-name">{{ SOURCE_GROUP_LABELS[group] }}</span>
              <!-- ⚠ What the catalogue HOLDS, not what matched: the number must not move as the
                   reader types, or it stops being the thing that says a group is worth opening. -->
              <span class="anchor-group-n">{{ heldInGroup(group).toLocaleString() }}</span>
            </button>
          </div>
          <div
            role="listbox"
            :aria-label="
              sourceGroupLabel === null
                ? 'Genomes in this catalogue'
                : `Genomes in ${sourceGroupLabel}`
            "
          >
          <button
            type="button"
            class="anchor-hit clear"
            role="option"
            :aria-selected="cursor === -1 ? 'true' : 'false'"
            @mousedown.prevent
            @click="choose(null)"
          >
            most common arrangement
          </button>
          <!-- ⛔ "holds none" and "none matched" are DIFFERENT sentences — see `emptyGroupReason`.
               The first is true of the complete-genome list on every catalogue loaded today, and a
               reader shown the second would keep typing at a list that can never answer. -->
          <div v-if="emptyGroupReason !== null" class="hit-none">{{ emptyGroupReason }}</div>
          <button
            v-for="(genome, index) in answer.genomes"
            :key="genome.sample_id"
            type="button"
            class="anchor-hit"
            role="option"
            :aria-selected="index === cursor ? 'true' : 'false'"
            @mousedown.prevent
            @click="choose(genome.sample_id)"
          >
            <span>{{ genome.sample_id }}</span>
            <!-- ⚠ `locus_count` — loci where this genome has a GENE. The published page could only
                 count loci where it appears in a recorded neighbourhood (`arrangement_locus_count`),
                 which undercounts wherever a gene got no window. Per LOCUS either way. -->
            <span class="anchor-n">{{ genome.locus_count.toLocaleString() }} loci</span>
          </button>
          </div>
          <div v-if="answer.truncated" class="hit-none alt-desc">
            {{ answer.genomes.length }} of {{ answer.matched_genome_count.toLocaleString() }} matching
            genomes<template v-if="sourceGroupLabel !== null"> in {{ sourceGroupLabel }}</template>
            — type more of the accession to narrow it.
          </div>
        </template>
      </template>
    </div>
  </div>
</template>
