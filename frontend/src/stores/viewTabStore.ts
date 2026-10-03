/**
 * Which of the four views of one locus is showing.
 *
 * ⭐ **The track sits above the tabs and stays visible in every view**, so a reader on any tab can
 * still walk. Walking is a request to SEE the locus, so it brings the reader back to the evidence —
 * except from Sequence, which is per-locus content too: a reader stepping along the track there is
 * exactly the one who wants the sequence to follow, not to be thrown back a tab on every step
 * (`app.js::show`).
 *
 * ⛔ **Absolute, never a toggle** — the same rule as every open/closed flag in this app. And the tab
 * is NOT in the URL: it is not a place (`lib/locusHashRoute.ts`).
 */

import { defineStore, storeToRefs } from "pinia";
import { ref, watch } from "vue";

import { useFunctionBlockStore } from "./functionBlockStore";
import { useGeneSequenceStore } from "./geneSequenceStore";
import { useLocusNavigationStore } from "./locusNavigationStore";

/**
 * In the order the strip draws them. `locus` is home. ⭐ The three views OF THIS LOCUS come first and
 * the site's own documentation last (David, 2026-09-22) — it had sat second, between the evidence and
 * the sequence.
 */
export const VIEW_TABS = [
  { id: "locus", label: "Syntelog Loci" },
  { id: "sequence", label: "Sequence" },
  // ⭐ `Function`, not `EggNOG` (David, 2026-10-03): *"Change 'EggNOG' to just be a 'Function' tab
  // and it can include the EggNOG annotations."* The tab names what a reader is looking for; EggNOG
  // is one of the sources inside it, and the panel's own heading has read `Function` all along — the
  // strip was the only place still naming the database instead of the question.
  { id: "function", label: "Function" },
  { id: "navigating", label: "Navigating BacAtlas" },
] as const;

export type ViewId = (typeof VIEW_TABS)[number]["id"];

export const useViewTabStore = defineStore("viewTab", () => {
  const navigation = useLocusNavigationStore();
  const functionBlock = useFunctionBlockStore();
  const geneSequence = useGeneSequenceStore();
  const { drawable } = storeToRefs(navigation);

  const view = ref<ViewId>("locus");

  /**
   * Show a view. The two tabs that fetch are told whether they are open, which is what lets a
   * reader who never opens them fetch nothing at all.
   */
  function showView(next: ViewId): void {
    view.value = next;
    functionBlock.setOpen(next === "function");
    geneSequence.setOpen(next === "sequence");
  }

  /**
   * Walking to a new locus brings the reader back to its evidence — unless they are reading its
   * sequence. Keyed on the DRAWN locus, like every per-locus store, so a response that is still in
   * flight does not move the reader off a tab under a locus that is still on screen.
   */
  watch(
    () => drawable.value?.locus.label ?? null,
    (label, previousLabel) => {
      // ⚠ NOT skipped when there was no previous locus: a reader whose first link was dead, who then
      // opened the documentation and clicked a chip, asked to SEE that locus (`app.js::show`:
      // `current !== i` held for -1 too).
      if (label === null || label === previousLabel) return;
      if (view.value !== "locus" && view.value !== "sequence") showView("locus");
    },
  );

  return { view, showView };
});
