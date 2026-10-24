/**
 * The anchor: one named genome, drawn instead of the commonest arrangement.
 *
 * ⭐ **Session state, deliberately not in the URL.** The hash already carries a trailing `r` for
 * walk direction, and putting the anchor there too would make every anchor change a history entry
 * (`app.js:195-200`). Switching species is a full navigation to a different catalogue, so the
 * anchor drops with it — which is right, because the two BioSample sets are disjoint.
 *
 * ⭐ **One assignment, not a walk over mounted boxes.** `app.js:3950`'s `setAnchor` walks `ANCHORS`
 * writing through every box that is mounted, because *"the boxes are two views of one state, not
 * two states that are kept in step"*. Here that sentence is the implementation: there is one ref,
 * both boxes read it, and there is nothing to keep in step.
 *
 * ⚠ `null` is "no anchor", and the page must be able to hold it with the control absent altogether:
 * a catalogue with no membership offers no anchor box rather than one that could only ever answer
 * "not found".
 */

import { defineStore } from "pinia";
import { computed, ref } from "vue";

import type { GenomeSourceGroup } from "@/api/types";

/**
 * ⛔ **The order the picker offers the two lists in, and it is not alphabetical.** Complete genomes
 * first because they are the ones a reader anchors deliberately — a closed long-read genome is the
 * useful ruler — while the ATB list is the long tail they search. `ncbi_draft` is last because it
 * is not one of the two lists the reader was promised; it exists so that an NCBI-keyed assembly
 * which is not closed appears SOMEWHERE rather than vanishing from a control that claims to list
 * every genome.
 */
/*
 * ⚠ A const TUPLE, not `readonly GenomeSourceGroup[]`: under `noUncheckedIndexedAccess` the array
 * form makes `SOURCE_GROUP_ORDER[0]` possibly `undefined`, and the fallback below would then be
 * able to put `undefined` into `sourceGroup` — which the API client would send as no filter at all.
 * `satisfies` keeps the values checked against the union while preserving the literal types.
 */
export const SOURCE_GROUP_ORDER = ["ncbi_complete", "atb", "ncbi_draft"] as const satisfies readonly GenomeSourceGroup[];

/** What each list is called, and what it is. Shown as a heading, so it says the provenance. */
/* ⚠ Short enough to fit three abreast in the popover — measured by screenshotting it, not guessed.
 * The full sentence lives in `SOURCE_GROUP_DESCRIPTIONS`, which is the tab's `title`. */
export const SOURCE_GROUP_LABELS: Readonly<Record<GenomeSourceGroup, string>> = {
  ncbi_complete: "Complete (NCBI)",
  atb: "AllTheBacteria",
  ncbi_draft: "Draft (NCBI)",
};

/**
 * ⛔ **The two lists the reader is PROMISED, offered whether or not this catalogue has any.**
 *
 * The published page's rule was *a control that cannot tell "this genome is not here" from "the
 * page has no idea" must not be offered* — and hiding an empty group obeyed it. But the census
 * removes that rule's precondition: the control now KNOWS the catalogue holds none, and can say so.
 * Hiding the tab instead leaves a reader unable to tell "this catalogue has no complete genomes
 * yet" from "this site does not do complete genomes", which is the very distinction the census
 * exists to make — and on every catalogue loaded today the complete-genome group is empty, so
 * hiding it means nobody ever learns the list exists.
 *
 * ⚠ `ncbi_draft` is NOT promised, so it appears only when it has members: it is the leftover
 * bucket for an NCBI assembly that is not closed, not a list anyone asked for.
 */
export const ALWAYS_OFFERED_SOURCE_GROUPS: readonly GenomeSourceGroup[] = ["ncbi_complete", "atb"];

/** ⚠ One sentence each, because "complete" and "draft" are claims a reader may not share. */
export const SOURCE_GROUP_DESCRIPTIONS: Readonly<Record<GenomeSourceGroup, string>> = {
  ncbi_complete: "Closed long-read assemblies — chromosome and plasmids circularised.",
  atb: "Short-read draft assemblies, keyed on BioSample.",
  ncbi_draft: "Long-read assemblies from NCBI that are not closed.",
};

/**
 * What an EMPTY list says. ⛔ Written per group rather than composed from the label, because
 * "This catalogue holds no complete (ncbi)." is what composing produced — a sentence that reads
 * like a bug. It also says *yet* for the complete genomes, which is true and is the honest thing
 * to tell a reader: they are not ingested, not absent in principle.
 */
export const SOURCE_GROUP_EMPTY_SENTENCES: Readonly<Record<GenomeSourceGroup, string>> = {
  ncbi_complete: "No complete NCBI genomes are in this catalogue yet.",
  atb: "This catalogue holds no AllTheBacteria assemblies.",
  ncbi_draft: "This catalogue holds no draft NCBI genomes.",
};

/**
 * ⛔ **Which RELATION this genome has to the catalogue, and it is not decoration.** A `catalogue`
 * genome is one of the modelled hundred: its genes are in the model, and it is COUNTED in the
 * arrangements it occupies. A `projected` genome was never clustered: its genes were placed on the
 * loci of their nearest modelled genes afterwards, and it MATCHES an arrangement without being
 * counted in one. The two are fetched by different query parameters, cached under different keys,
 * and described by different sentences — and a page that lost the distinction would tell a reader
 * their genome is part of a model it is not part of.
 */
export type AnchorKind = "catalogue" | "projected";

export const useAnchorGenomeStore = defineStore("anchorGenome", () => {
  /** The anchored genome's `sample_id`, or `null` for none. */
  const sampleId = ref<string | null>(null);
  /** ⚠ Meaningless while `sampleId` is null, and never read then. */
  const kind = ref<AnchorKind>("catalogue");
  /** Whether this catalogue can offer an anchor at all — false removes the control. */
  const isAvailable = ref(false);

  /**
   * ⭐ How many genomes this catalogue holds per source group, from the SHELL.
   *
   * ⛔ **Independent of any query, and that is the whole reason it lives here rather than being
   * read off a listing.** The picker must be able to say *"this catalogue holds no complete NCBI
   * genomes"* — which is true of every catalogue loaded today — as distinct from *"none match what
   * you typed"*. A count taken from a filtered list can only ever say the second.
   */
  const sourceGroupCounts = ref<Readonly<Record<GenomeSourceGroup, number>>>({
    ncbi_complete: 0,
    ncbi_draft: 0,
    atb: 0,
  });

  /**
   * Which list the picker is showing. ⛔ **`null` means UNFILTERED — every genome, no tabs** — and
   * it is the state the control degrades to, not an error state. The API reads an absent
   * `source_group` as "every group", so null is a real, correct value rather than a missing one.
   */
  const sourceGroup = ref<GenomeSourceGroup | null>(null);

  /**
   * Whether a census was supplied at all — ⛔ distinct from one that is all zeros.
   *
   * "The server did not tell us" and "the catalogue holds none" are different facts and lead to
   * different controls: the first degrades to the unfiltered pre-census list with no tabs, the
   * second shows the tabs and says which group is empty. A single all-zero record cannot tell them
   * apart, which is why this is its own flag rather than a `sum(counts) === 0` test.
   */
  const hasCensus = ref(false);

  /**
   * The groups this catalogue can actually offer, in `SOURCE_GROUP_ORDER`.
   *
   * ⚠ **A group with no members is NOT offered as a tab**, on the published page's own rule for the
   * anchor box itself: *a control that cannot tell "this genome is not here" from "the page has no
   * idea" must not be offered.* An empty tab makes the same face at both. What replaces it is a
   * sentence in the list saying the catalogue holds none — see `AnchorGenomeBox`.
   */
  const offeredSourceGroups = computed(() => {
    // ⚠ No census at all means no tabs and no filter — the pre-census behaviour. Offering the two
    // promised groups here would filter a list on counts we do not have.
    if (!hasCensus.value) return [];
    return SOURCE_GROUP_ORDER.filter(
      (group) =>
        sourceGroupCounts.value[group] > 0 || ALWAYS_OFFERED_SOURCE_GROUPS.includes(group),
    );
  });

  /** Whether the reader has a choice to make at all. One list needs no tabs. */
  const hasSourceGroupChoice = computed(() => offeredSourceGroups.value.length > 1);

  const isAnchored = computed(() => sampleId.value !== null);
  const isProjected = computed(() => sampleId.value !== null && kind.value === "projected");

  /** Absolute, like every other setter here: say what it becomes.
   *
   * ⛔ **Setting either kind releases the other**, because there is one anchor. That is the whole
   * reason both live in this store rather than one each: two independent refs would let a reader
   * anchor a modelled genome and one of their own at the same time, and the track can only draw one.
   */
  function setAnchor(nextSampleId: string | null, nextKind: AnchorKind = "catalogue"): void {
    // ⚠ `""` from a cleared input is "no anchor", not a genome named the empty string. Normalising
    // here is what keeps the cache key stable — see `locusCacheKey`.
    sampleId.value = nextSampleId ? nextSampleId : null;
    kind.value = nextSampleId ? nextKind : "catalogue";
  }

  /** The projected half of the same one anchor. */
  function setProjectedAnchor(nextSampleId: string | null): void {
    setAnchor(nextSampleId, "projected");
  }

  function clearAnchor(): void {
    sampleId.value = null;
    kind.value = "catalogue";
  }

  /** Called when a catalogue loads: does it carry the membership an anchor needs?
   *
   * ⚠ Only a CATALOGUE anchor depends on that membership. A projected genome is anchored through
   * its own table, so an unavailable catalogue anchor must not silently drop one.
   */
  function setAvailability(available: boolean): void {
    isAvailable.value = available;
    if (!available && kind.value === "catalogue") sampleId.value = null;
  }

  /**
   * Called when a catalogue loads: the per-group census, and the group to open on.
   *
   * ⛔ **The opening group is the first NON-EMPTY one, never a hardcoded default.** Opening on
   * "Complete genomes (NCBI)" would show every reader of every catalogue loaded today an empty
   * list as the first thing the control does — the group is empty on all four, measured. Falling
   * back to `SOURCE_GROUP_ORDER[0]` when nothing has members keeps the value a valid group rather
   * than `undefined`, which the API would then send as a filter.
   */
  function setSourceGroupCounts(
    counts: Readonly<Partial<Record<GenomeSourceGroup, number>>> | null | undefined,
  ): void {
    // ⛔⛔ **Tolerates an ABSENT census, and that is a production requirement, not test comfort.**
    // Code and data deploy by separate paths, so this bundle can be talking to a server that
    // predates `genome_source_groups`, or reading a shell response cached with `max-age=86400`
    // before it existed. The first version of this indexed `counts` directly and threw
    // `Cannot read properties of undefined`, taking the whole page down — found by
    // `pageShellBehaviour.test.ts`, whose fixtures are exactly that older shape.
    const normalised = {
      ncbi_complete: counts?.ncbi_complete ?? 0,
      ncbi_draft: counts?.ncbi_draft ?? 0,
      atb: counts?.atb ?? 0,
    };
    sourceGroupCounts.value = normalised;
    hasCensus.value = counts != null;
    // ⛔ **No census means NO FILTER, never a guessed group.** Defaulting to a group here would
    // send a filter the server answers with an empty list, and the reader would see a picker that
    // has silently stopped listing their genomes. Unfiltered is the pre-existing behaviour.
    sourceGroup.value = counts == null
      ? null
      : (SOURCE_GROUP_ORDER.find((group) => normalised[group] > 0) ?? SOURCE_GROUP_ORDER[0]);
  }

  /** ⚠ Switching list does NOT clear the anchor: the anchored genome is still anchored, it is
   * simply in the other list. Clearing it would make browsing the groups destructive. */
  function setSourceGroup(next: GenomeSourceGroup | null): void {
    sourceGroup.value = next;
  }

  return {
    sampleId,
    kind,
    isAvailable,
    isAnchored,
    isProjected,
    sourceGroup,
    sourceGroupCounts,
    hasCensus,
    offeredSourceGroups,
    hasSourceGroupChoice,
    setAnchor,
    setProjectedAnchor,
    clearAnchor,
    setAvailability,
    setSourceGroup,
    setSourceGroupCounts,
  };
});
