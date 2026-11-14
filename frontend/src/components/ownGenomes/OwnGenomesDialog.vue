<script setup lang="ts">
import type { CatalogueKey } from "@/api/types";
/**
 * "View your own genomes" — the accounts placeholder, the add-genomes step, and what was added.
 *
 * ⭐ **A native `<dialog>`, opened with `showModal()`** — it brings the backdrop, the focus trap and
 * Escape with it, and none of that is worth hand-writing. It cannot be a URL: the hash on this page
 * is the locus label (`lib/locusHashRoute.ts`). When accounts exist this becomes a real route.
 *
 * ⛔ **It says what is NOT built, on the first screen a reader meets.** Accounts, saving, the BakRep
 * download and the ESM-C → Bacformer forward pass are all described as what the service does and
 * marked "to be completed" — David wants collaborators to see the whole shape, without being able to
 * mistake the sketch for the service.
 *
 * ⚠ The file is read with `File.text()` in the browser. Nothing is uploaded, nothing is saved, and it
 * is gone on reload — which the first screen also says.
 */
import { computed, ref, watch } from "vue";

import type { SpeciesEntry } from "@/stores/speciesCatalogueStore";
import { useAnchorGenomeStore } from "@/stores/anchorGenomeStore";
import { useOwnGenomesStore, type LineOutcome } from "@/stores/ownGenomesStore";
import { MAXIMUM_LINES } from "@/lib/parseAccessionList";

const props = defineProps<{
  speciesKey: string | null;
  species: readonly SpeciesEntry[];
}>();

const own = useOwnGenomesStore();
const anchor = useAnchorGenomeStore();

/**
 * ⭐ **The GPU figure is measured** (David, 2026-09-22: "I thought ESM and bacformer single forward
 * was around 15 seconds per genome. Can we please check this"). It is — CSD3 job 36053709, one
 * A100, bf16, the pinned ESM-C + Bacformer path over three probe genomes: **median 14.3 s per
 * genome**, of which Bacformer is 0.12–0.47 s and ESM-C is all the rest. Rounded to 15 on purpose:
 * the three spanned 13–23 s (the first pays CUDA warm-up), so a decimal would claim a precision the
 * measurement does not have. ⚠ It is the HuggingFace attention path — the run warns faESM is not
 * installed — so a tuned service is faster than what is quoted here, never slower.
 *
 * ⭐ **And so is the CPU figure** — CSD3 job 36055321, 32 Icelake cores, **fp32 because an Icelake
 * core has no native bf16 matmul**: 527.7 s and 372.2 s for the same two genomes, so ~7.5 min,
 * quoted as 8. ⚠ **The core count is part of the number and is said on screen**: this is a 32-core
 * server, and a 4-core laptop is roughly eight times slower again. (The laptop measurement was
 * abandoned — MPS ran out of memory at 13 GiB on a 16 GB machine — and is not missed: 32 cores is
 * what a CPU deployment would have.) Both numbers and their jobs: nuna's `PROJECT_STATE.md` §6.
 */
const GPU_SECONDS_PER_GENOME = 15;
const CPU_MINUTES_PER_GENOME = 8;
const CPU_CORES_TIMED = 32;

const element = ref<HTMLDialogElement | null>(null);

/** ⚠ Guarded: jsdom has no `showModal`, and a test that mounts this must not blow up on it. */
watch(
  () => own.isOpen,
  (open) => {
    const dialog = element.value;
    if (dialog === null) return;
    if (open && typeof dialog.showModal === "function" && !dialog.open) dialog.showModal();
    if (!open && typeof dialog.close === "function" && dialog.open) dialog.close();
  },
);

/** Escape and the backdrop close the dialog themselves; the store has to hear about it. */
function onClose(): void {
  own.close();
}

const speciesName = computed(
  () => props.species.find((entry) => entry.key === props.speciesKey)?.scientific_name ?? "this catalogue",
);

/**
 * ⛔ Species whose catalogue has NO placed genomes yet, so the demonstration below would walk a reader
 * through a workflow that has nothing behind it (David, 2026-11-14: placing a genome on kp-71k is not
 * built). Remove a key here the day its catalogue carries projected genomes.
 */
const PLACEMENT_PENDING_SPECIES: ReadonlySet<string> = new Set(["kp"]);
const placementPending = computed(
  () => props.speciesKey !== null && PLACEMENT_PENDING_SPECIES.has(props.speciesKey),
);
/** Where the worked example lives: the first published species whose placement is not pending. */
const exampleSpecies = computed(
  () => props.species.find((entry) => entry.catalogue_key !== null && !PLACEMENT_PENDING_SPECIES.has(entry.key)) ?? null,
);

/** The rosters are loaded when the reader reaches the step that needs them, never before. */
watch(
  () => [own.isOpen, own.step] as const,
  ([open, step]) => {
    if (!open || step !== "add") return;
    if (own.isLoadingRosters || Object.keys(own.modelledBySpecies).length > 0) return;
    // ⚠ TWO keys per entry, deliberately: `key` files the answers under the species, and
    // `catalogueKey` decides which catalogue's placements are read. A species with no default
    // serves nothing and is skipped rather than addressed by its species key.
    const roster = props.species
      .filter((entry) => entry.catalogue_key !== null)
      .map((entry) => ({
        key: entry.key,
        catalogueKey: entry.catalogue_key as CatalogueKey,
        scientificName: entry.scientific_name,
      }));
    // ⛔ Both lists, together. "Placed here", "modelled here" and "neither" are three different
    // answers, and with only one list loaded the dialog would give the wrong one of the three.
    void own.loadRosters(roster);
    void own.loadPlacedGenomes(roster);
  },
);

/** The genomes placed on the catalogue now open — what step 3 actually shows. */
const placed = computed(() => own.placedFor(props.speciesKey));

/** Anchoring from here closes the dialog: the reader asked to SEE the genome, not to read about it. */
function anchorPlaced(sampleId: string): void {
  anchor.setProjectedAnchor(sampleId);
  own.close();
}

/**
 * ⚠ The agreement and its denominator, always together. A locus with *m* modelled genes can supply
 * at most min(n, m) of the n checkers, so "10 of 10" and "1 of 1" are both unanimous and the second
 * says far less — printing the numerator alone would make a singleton locus look like a failure.
 */
function share(numerator: number, denominator: number): string {
  return `${numerator.toLocaleString()} of ${denominator.toLocaleString()}`;
}

function percent(part: number, whole: number): string {
  return whole === 0 ? "—" : `${Math.round((100 * part) / whole)}%`;
}

async function onFile(event: Event): Promise<void> {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file === undefined) return;
  own.setFile(file.name, await file.text());
}

const rows = computed(() =>
  own.lines.map((line) => ({ line, outcome: own.outcomeFor(line, props.speciesKey) })),
);

/** ⚠ Counted from the outcomes, never from the line count: a repeat is not a genome. */
const showable = computed(() => rows.value.filter((row) => row.outcome.kind === "placed-here").length);
/** The ones the service would have to fetch and place — the part that is not built. */
const toBeCompleted = computed(
  () => rows.value.filter((row) => row.outcome.kind === "to-be-completed").length,
);

const MARK: Readonly<Record<LineOutcome["kind"], string>> = {
  "not-an-accession": "✗",
  repeat: "–",
  "modelled-here": "–",
  "modelled-elsewhere": "→",
  "placed-here": "✓",
  "placed-elsewhere": "→",
  "to-be-completed": "○",
  unchecked: "?",
};

function sentenceFor(outcome: LineOutcome): string {
  switch (outcome.kind) {
    case "not-an-accession":
      return "not a BioSample — one per line, SAMN…, SAMEA… or SAMD…";
    case "repeat":
      return "the same genome again";
    case "modelled-here":
      return (
        `already one of the 100 genomes ${speciesName.value} was modelled from — ` +
        "anchor to it with “A genome in BacAtlas”"
      );
    case "modelled-elsewhere":
      return `one of the 100 modelled genomes of ${outcome.scientificName} — switch species to anchor it`;
    case "placed-here":
      return (
        `ready — ${outcome.entry.placed_gene_count.toLocaleString()} genes placed on ` +
        `${outcome.entry.distinct_locus_count.toLocaleString()} loci, in advance for this demonstration`
      );
    case "placed-elsewhere":
      return `placed on the ${outcome.scientificName} catalogue — switch species to see it`;
    case "to-be-completed":
      return "a BioSample we do not hold yet — with an account it would be fetched from BakRep and placed (to be completed)";
    case "unchecked":
      return "the modelled genome list did not load, so this was not checked";
  }
}
</script>

<template>
  <dialog ref="element" class="own-dialog" aria-labelledby="own-dialog-title" @close="onClose">
    <div v-if="own.isOpen" class="own-dialog-in">
      <header class="own-dialog-head">
        <h2 id="own-dialog-title">View your own genomes</h2>
        <button type="button" class="own-close" aria-label="Close" @click="own.close()">✕</button>
      </header>

      <!-- ── 0 · not built for this catalogue yet ───────────────────────────────────────────── -->
      <section v-if="placementPending" class="own-step own-pending">
        <p class="lede">Pending implementation.</p>
        <p class="muted">
          View the
          <a v-if="exampleSpecies" :href="`?catalogue=${exampleSpecies.catalogue_key}`">
            <i>{{ exampleSpecies.scientific_name }}</i> pangenome</a>
          <template v-else>E. coli pangenome</template>
          for a precomputed example of the workflow.
        </p>
      </section>

      <!-- ── 1 · accounts, and what the service does ─────────────────────────────────────────── -->
      <section v-else-if="own.step === 'account'" class="own-step">
        <p class="lede">Accounts and login will be required — to be completed.</p>
        <p class="muted">
          With an account, the genomes you add are saved and are here when you come back.
        </p>
        <p class="muted">
          Each genome is fetched from <b>BakRep</b> — the Bakta annotation of its assembly — every
          protein is embedded with <b>ESM-C</b>, the genome is run through <b>Bacformer</b>, and each
          gene is placed on the locus of its nearest modelled gene. That is about
          <b>{{ GPU_SECONDS_PER_GENOME }} s per genome on an A100</b>, or about
          {{ CPU_MINUTES_PER_GENOME }} minutes on {{ CPU_CORES_TIMED }} CPU cores — both measured,
          and almost all of either is ESM-C rather than Bacformer.
        </p>
        <p class="muted own-caveat">
          Until then this is a demonstration: nothing you add is saved, fetched or computed, and it is
          forgotten when you reload.
        </p>
        <div class="own-actions">
          <button type="button" class="own-button" disabled>Sign in — to be completed</button>
          <button type="button" class="own-button own-go" @click="own.showStep('add')">
            Continue with the demo →
          </button>
        </div>
      </section>

      <!-- ── 2 · add genomes from a text file ────────────────────────────────────────────────── -->
      <section v-else-if="own.step === 'add'" class="own-step">
        <p class="lede">Add genomes · {{ speciesName }}</p>
        <p class="muted">
          A text file with <b>one BioSample per line</b> (SAMN…, SAMEA…, SAMD…). Lines starting with
          <code>#</code> are ignored, and the first column of a CSV or TSV is taken. Up to
          {{ MAXIMUM_LINES }} lines. The file is read here in your browser — nothing is uploaded.
        </p>
        <div class="own-actions">
          <label class="own-button own-file">
            Choose a .txt file
            <input type="file" accept=".txt,.tsv,.csv,text/plain" @change="onFile" />
          </label>
          <span v-if="own.fileName" class="muted own-file-name">{{ own.fileName }}</span>
        </div>

        <p v-if="own.isLoadingRosters" class="muted">Checking your accessions against the catalogue…</p>
        <p v-else-if="own.rosterFailure" class="pop-error" role="alert">
          The modelled genome list did not load — {{ own.rosterFailure.detail }}. Your accessions are
          listed, but they could not be checked against the catalogue.
        </p>

        <ul v-if="rows.length" class="own-lines">
          <li v-for="row in rows" :key="`${row.line.lineNumber}-${row.line.text}`" :class="row.outcome.kind">
            <span class="own-line-n">{{ row.line.lineNumber }}</span>
            <span class="own-line-mark" aria-hidden="true">{{ MARK[row.outcome.kind] }}</span>
            <span class="own-line-id mono">{{ row.line.text }}</span>
            <span class="own-line-say">{{ sentenceFor(row.outcome) }}</span>
          </li>
        </ul>
        <p v-else-if="own.fileName" class="muted">
          No accession in {{ own.fileName }} — every line was blank or a comment.
        </p>

        <p v-if="own.list?.truncated" class="muted alt-desc">
          Only the first {{ MAXIMUM_LINES }} lines were read.
        </p>
        <p v-if="own.list?.tooLarge" class="muted alt-desc">
          That file is larger than 64 kB, so only its beginning was read.
        </p>

        <div class="own-actions">
          <button type="button" class="own-button" @click="own.showStep('account')">← Back</button>
          <button type="button" class="own-button own-go" @click="own.showStep('yours')">
            <!-- ⚠ The label names what will actually be shown. With nothing placed on this
                 catalogue, "Show 0 genomes" would read as a failure of the reader's file rather
                 than as the part of the service that is not built yet. -->
            {{
              showable === 0
                ? "What happens next →"
                : showable === 1
                  ? "Show 1 genome →"
                  : `Show ${showable} genomes →`
            }}
          </button>
        </div>
      </section>

      <!-- ── 3 · what was added ──────────────────────────────────────────────────────────────── -->
      <section v-else class="own-step">
        <p class="lede">Your genomes · {{ speciesName }}</p>
        <p class="muted own-banner">
          <b>Demonstration — placed in advance.</b> None of these genomes was clustered by the model.
          Each of its genes was placed on the locus of its <b>nearest modelled gene</b>, and its
          neighbours checked that placement: an approximation of what the model would have done, not
          the model's own clustering. Fetching a genome from BakRep and embedding it on demand is the
          part still to be built.
        </p>

        <ul v-if="placed.length" class="own-placed">
          <li v-for="entry in placed" :key="entry.sample_id">
            <div class="own-placed-head">
              <span class="mono own-placed-id">{{ entry.sample_id }}</span>
              <button type="button" class="own-button own-go" @click="anchorPlaced(entry.sample_id)">
                Show it on the track →
              </button>
            </div>
            <div class="own-placed-facts">
              <span>
                <b>{{ entry.placed_gene_count.toLocaleString() }}</b> genes placed on
                <b>{{ entry.distinct_locus_count.toLocaleString() }}</b> loci
              </span>
              <span>
                <b>{{ percent(entry.window_matched_gene_count, entry.placed_gene_count) }}</b> sit in a
                neighbourhood the 100 already have
              </span>
              <span>
                nearest-gene similarity
                <b>{{ entry.nearest_cosine.median?.toFixed(3) ?? "—" }}</b> median, down to
                {{ entry.nearest_cosine.minimum?.toFixed(2) ?? "—" }}
              </span>
              <span>
                {{ share(entry.contested_gene_count, entry.placed_gene_count) }} genes have a rival
                locus with more support
              </span>
              <span v-if="entry.genes_without_a_neighbourhood > 0" class="alt-desc">
                {{ entry.genes_without_a_neighbourhood.toLocaleString() }} genes are alone on their
                contig, so they have no neighbourhood at all
              </span>
            </div>
          </li>
        </ul>
        <p v-else class="muted">
          No genome has been placed on this catalogue. Add some, or switch species.
        </p>

        <p class="muted alt-desc">
          Every count elsewhere on this page — prevalence, the bands, the census, the neighbourhoods —
          is the 100 modelled genomes, and stays that way.
        </p>
        <div class="own-actions">
          <button type="button" class="own-button" @click="own.showStep('add')">← Add genomes</button>
          <button type="button" class="own-button own-go" @click="own.close()">Close</button>
        </div>
      </section>
    </div>
  </dialog>
</template>
