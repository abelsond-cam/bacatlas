<script setup lang="ts">
/**
 * **Sequence diversity and Bakta annotations — ONE card, deliberately.**
 *
 * ⭐ They were two, and the second largely restated the first: UniRef50 is the sequence fact and the
 * Bakta name is *derived from it*, so a locus's families and its names are one table, not two lists
 * in two boxes. That separation is what hid `ybeF` (David, 2026-08-19) — its families split 63/33
 * and its symbols 63/35, in different cards, and the correspondence between them had to be guessed.
 * The join ships, so a row is a statement about its own genes.
 *
 * The Pfam verdict rides here for the same reason: it is a verdict about *these families*, and it
 * had a whole card whose visible content was the accession chips this table already carries.
 *
 * ⛔ **Both verdicts are READ, never re-derived.** `app.js:3208` — the page once re-derived the Pfam
 * class by a different rule and disagreed with the audit on 2 of 22,624 loci. *"A page that quotes
 * the report must not be able to contradict it, ever."*
 */
import { computed } from "vue";

import type { AnnotationEntry, Locus, PfamFamilyReference, UnirefFamily } from "@/api/types";
import { collapseTierVerdict, pfamAccessionsIn, pfamVerdict } from "@/lib/evidenceVocabulary";
import CountTable from "@/components/shared/CountTable.vue";
import PfamChip from "@/components/shared/PfamChip.vue";

const props = defineProps<{
  locus: Locus;
  families: readonly UnirefFamily[];
  /** Every gene symbol in use here, so the card can say how many names there are to reconcile. */
  symbols: readonly AnnotationEntry[];
  pfamReference: Readonly<Record<string, PfamFamilyReference>>;
  /** How many Pfam architectures are LISTED, which may be fewer than the locus has. */
  listedArchitectureCount: number;
}>();

const size = computed(() => props.locus.gene_count);

/**
 * ⛔ Claim with the MAJOR-family count — families holding ≥10 % of the labelled genes, the audit's
 * own bar. The table below lists every family including single-gene stragglers, and claiming with
 * that number would overstate the split at every locus with a long tail.
 */
const majorFamilyCount = computed(
  () => props.locus.uniref50.major_family_count || props.locus.uniref50.family_count,
);

const naming = computed(() => {
  const named = props.locus.named_gene_count;
  const unnamed = size.value - named;
  if (named === 0) {
    return {
      kind: "none" as const,
      text: `None of these ${size.value} genes carries a gene name — a coherent locus with no symbol at all.`,
    };
  }
  if (unnamed > 0) return { kind: "some" as const, unnamed, total: size.value };
  return { kind: "all" as const, text: `Every one of these ${size.value} genes already carries a gene name.` };
});

/**
 * ⭐ **The vote that NAMES the locus, shown as a vote.** The locus is called its commonest Bakta
 * symbol counted as an exact string, and until this line existed the card showed that number
 * nowhere: the only counts standing beside a gene name were the FAMILY gene counts in the table
 * below, so kp 3992 — named `mviN` on 5 genes against `rfbX`’s 2, with 90 of 97 unnamed — read as
 * *"9 genes are rfbX and 4 are mviN"* off two family rows whose counts are not symbol counts at
 * all. That misreading is what this block answers, and it costs one line of data already served.
 *
 * ⛔ **Two remainders, each its own sentence**, because they are different facts: genes carrying a
 * name outside the top-N cut, and genes carrying no name at all. Folding them into one "other"
 * would let missing annotation read as a rival name.
 *
 * ⚠ **A TIE is named as a tie.** Rank 0 and rank 1 equal means the name was settled by
 * alphabetical order and nothing else — true of 29 *E. coli* and 12 kp loci — and a card that
 * printed the winner silently would be presenting a coin toss as a finding.
 *
 * ⭐ **`corroborated_symbol`: a thin vote the locus's other annotation settled.** The winner is then
 * the NAME, not rank 0 — kp 722 lists `gumC` 1 above `wzc` 1 and is named `wzc` because PF23607
 * names it — and the note says which accession did, so the mark never reads as the count's verdict.
 */
const symbolVote = computed(() => {
  const listedEntries = props.symbols;
  if (!listedEntries.length) return null;
  const listed = listedEntries.reduce((sum, entry) => sum + entry.gene_count, 0);
  const unnamed = Math.max(0, size.value - props.locus.named_gene_count);
  const unlistedNamed = Math.max(0, props.locus.named_gene_count - listed);
  // One name, every gene carrying it: the headline and the sentence above already say this, and
  // repeating it as a one-item vote is noise.
  if (listedEntries.length === 1 && !unnamed && !unlistedNamed) return null;
  const top = listedEntries[0];
  if (top === undefined) return null;
  const runnerUp = listedEntries[1];
  return {
    entries: listedEntries,
    // ⚠ Only when the name came FROM this vote. An inferred name has no winner here to mark.
    winner:
      props.locus.display_name_source === "bakta_symbol"
        ? top.term
        : props.locus.display_name_source === "corroborated_symbol"
          ? props.locus.display_name
          : null,
    corroboratedBy:
      props.locus.display_name_source === "corroborated_symbol"
        ? (props.locus.display_name_source_accession ?? "their other annotation")
        : null,
    isTie: runnerUp !== undefined && runnerUp.gene_count === top.gene_count,
    unlistedNamed,
    unnamed,
  };
});

const tier = computed(() => collapseTierVerdict(props.locus.evidence.collapse_tier));
const verdict = computed(() =>
  pfamVerdict(props.locus.pfam.concordance_class, props.locus.pfam.annotated_gene_count),
);

/** Which optional columns the table earns — an empty column is worse than an absent one. */
const hasSymbols = computed(() => props.families.some((family) => family.modal_symbol));
const hasArchitectures = computed(() => props.families.some(contributesPfamColumn));

const headings = computed(() => {
  const out = ["UniRef50"];
  if (hasSymbols.value) out.push("gene name");
  out.push("modal Bakta product");
  if (hasArchitectures.value) out.push("Pfam");
  return out;
});

const rows = computed(() =>
  props.families.map((family) => ({
    key: family.uniref50_accession,
    count: family.gene_count,
    family,
    accessions: pfamAccessionsIn(family.modal_architecture),
    // ⚠ `null` is NOT MEASURED and shows no marker; `1` is measured agreement and shows none
    // either. Only a measured disagreement earns one, so the two absences collapse here rather
    // than in the template where a `> 1` on a null would read as false by accident.
    extraSymbolCount:
      family.distinct_symbol_count !== null && family.distinct_symbol_count > 1
        ? family.distinct_symbol_count - 1
        : null,
    // ⛔ Tested against null, never truthiness: a family where NO gene is annotated is the loudest
    // case and a `v-if` on the number would hide exactly that one.
    showsOwnCoverage:
      family.pfam_annotated_gene_count !== null &&
      family.pfam_annotated_gene_count < family.gene_count,
    // ⛔ The SAME rule for the gene-name column, which had no denominator at all: a family of 9
    // genes naming 2 drew a bare `rfbX` and read as nine genes agreeing. Tested against null for
    // the reason above — a family where nothing is named is the loudest case of all.
    showsOwnNameCoverage:
      family.named_gene_count !== null && family.named_gene_count < family.gene_count,
  })),
);

function unirefUrl(accession: string): string {
  return `https://www.uniprot.org/uniref/${encodeURIComponent(accession)}`;
}

/**
 * ⚠ `distinct_symbol_count` is a COUNT, not a list. The marker says how many and never invents the
 * others — showing a bare modal would let a family that is itself split across names read as
 * unanimous, which is the same mistake the two-card layout made.
 */
function symbolMarkerTitle(family: UnirefFamily): string {
  // ⚠ The denominator is quoted only where it SAYS something. "60 named genes of 60" is noise;
  // "2 named genes of 9" is the whole point.
  const named =
    family.named_gene_count === null || family.named_gene_count === family.gene_count
      ? `${family.gene_count} genes`
      : `${family.named_gene_count} named genes of ${family.gene_count}`;
  return (
    `${family.distinct_symbol_count} different gene names among this family's ` +
    `${named}; ${family.modal_symbol} is the commonest`
  );
}

/** ⚠ A family earns the Pfam column if it has an architecture OR a measured coverage count. */
function contributesPfamColumn(family: UnirefFamily): boolean {
  return Boolean(family.modal_architecture) || family.pfam_annotated_gene_count !== null;
}

/**
 * ⛔ **Coverage is stated AS COVERAGE.** A gene with no Pfam domain is *missing evidence*, not
 * evidence of a different architecture, and reading the second as the first is the
 * "14 of 99 have different Pfam structures" misreading a bare domain count invites.
 */
const coverage = computed(() => {
  const annotated = props.locus.pfam.annotated_gene_count;
  if (annotated === null) return null;
  if (annotated === 0) {
    return (
      `None of these ${size.value} genes carries a Pfam-A domain — Pfam has no coverage here, ` +
      "so it can neither support nor contradict this locus."
    );
  }
  const architectures = props.locus.pfam.architecture_count;
  const head =
    annotated < size.value
      ? `${annotated} of ${size.value} genes carry a Pfam-A domain; the other ${size.value - annotated} ` +
        "carry none at all — absent annotation, not a competing architecture"
      : `All ${size.value} genes carry a Pfam-A domain`;
  // ⚠ The TRUE architecture count, which the top-N list on screen cannot show — so the card says
  // how many are listed whenever that is fewer.
  const listed =
    architectures > props.listedArchitectureCount ? ` (${props.listedArchitectureCount} listed)` : "";
  const tail =
    architectures > 1
      ? ` · ${architectures} architectures across the locus${listed}`
      : " · one architecture across the locus";
  return `${head}${tail}.`;
});
</script>

<template>
  <div class="card" :class="verdict ? `pf-${verdict.tone}` : null">
    <h2>Sequence diversity and Bakta annotations</h2>

    <p v-if="majorFamilyCount > 1" class="lede">
      UniRef50 files these genes under <b>{{ majorFamilyCount }} families</b> — each holding at least
      a tenth of the labelled genes. Nuna holds them as one locus.
    </p>
    <p v-else-if="families.length >= 1" class="lede">
      One UniRef50 family — this locus and the sequence family agree.
    </p>

    <p class="lede">
      <template v-if="naming.kind === 'some'">
        <b>{{ naming.unnamed }}</b> of {{ naming.total }} genes here carry no gene name. This locus
        names them <b>{{ locus.display_name }}</b>.
      </template>
      <template v-else>{{ naming.text }}</template>
    </p>

    <!-- ⛔ The counts the NAME came from. Every other number on this card is a count of genes in a
         family; these are counts of genes carrying a name, and the two were indistinguishable on
         the page until this block existed. -->
    <p v-if="symbolVote" class="sym-vote">
      <span class="sym-vote-lead">gene names here</span>
      <span
        v-for="entry in symbolVote.entries"
        :key="entry.term"
        class="sym-vote-item"
        :class="{ won: entry.term === symbolVote.winner }"
      >
        <span class="acc">{{ entry.term }}</span> {{ entry.gene_count }}
      </span>
      <span v-if="symbolVote.unlistedNamed" class="sym-vote-item rest"
        >{{ symbolVote.unlistedNamed }} under names not listed</span
      >
      <!-- ⚠ A separate remainder from the one above: no name at all is not a rival name. -->
      <span v-if="symbolVote.unnamed" class="sym-vote-item rest"
        >{{ symbolVote.unnamed }} unnamed</span
      >
    </p>
    <!-- ⚠ A tie is the only case where the winner needs explaining, so it is the only case that
         gets a sentence. The counting rule itself is stated once, in the card's footnote. -->
    <p v-if="symbolVote && symbolVote.winner" class="muted sym-vote-note">
      <template v-if="symbolVote.corroboratedBy"
        >Too few genes are named for these counts to decide, so <b>{{ symbolVote.winner }}</b> names
        the locus as the one of these names {{ symbolVote.corroboratedBy }} also gives.</template
      >
      <template v-else-if="symbolVote.isTie"
        >The two commonest names are level, so <b>{{ symbolVote.winner }}</b> won this locus on
        alphabetical order alone — not on evidence.</template
      >
      <template v-else>The commonest of these names the locus.</template>
    </p>

    <!-- Both verdicts on one row: how the members hold together by sequence, and what Pfam makes
         of them. -->
    <div v-if="tier || verdict" class="chip-row">
      <span v-if="tier" class="chip" :class="tier.tone">{{ tier.label }}</span>
      <span v-if="verdict" class="chip" :class="verdict.tone">{{ verdict.label }}</span>
    </div>
    <p v-if="tier && tier.note" class="muted">{{ tier.note }}</p>
    <p v-if="verdict" class="muted verdict-note">{{ verdict.note }}</p>

    <CountTable v-if="families.length" :headings="headings" :rows="rows" :total="size">
      <template #cells="{ row }">
        <td class="acc">
          <a
            class="acc-link"
            :href="unirefUrl(row.family.uniref50_accession)"
            :title="`${row.family.uniref50_accession} on UniProt`"
            target="_blank"
            rel="noopener noreferrer"
          >{{ row.family.uniref50_accession }}</a>
        </td>

        <td v-if="hasSymbols" class="fam-sym">
          <span class="fam-sym-in">
            <template v-if="row.family.modal_symbol">
              <span class="acc">{{ row.family.modal_symbol }}</span>
              <span
                v-if="row.extraSymbolCount !== null"
                class="alt-n"
                :title="symbolMarkerTitle(row.family)"
              >+{{ row.extraSymbolCount }}</span>
              <!-- ⛔ `v-if` on the FLAG, never on the number — see `showsOwnNameCoverage`. -->
              <span v-if="row.showsOwnNameCoverage" class="fam-cov"
                >{{ row.family.named_gene_count }}/{{ row.family.gene_count }} named</span
              >
            </template>
            <span v-else class="none">—</span>
          </span>
        </td>

        <td>{{ row.family.modal_product || "—" }}</td>

        <td v-if="hasArchitectures" class="fam-pf">
          <span class="fam-pf-in">
            <template v-if="row.accessions.length">
              <PfamChip
                v-for="accession in row.accessions"
                :key="accession"
                :accession="accession"
                :family="pfamReference[accession]"
              />
            </template>
            <span v-else class="none">no domain</span>
            <!-- ⛔ `v-if` on a NUMBER would hide a measured zero — see `showsOwnCoverage`. -->
            <span v-if="row.showsOwnCoverage" class="fam-cov"
              >{{ row.family.pfam_annotated_gene_count }}/{{ row.family.gene_count }} annotated</span
            >
          </span>
        </td>
      </template>
    </CountTable>

    <p v-if="families.length" class="muted">
      A UniRef50 accession carries no name of its own, so each family is labelled by the commonest
      Bakta product and gene name among its members here<template v-if="symbols.length > 1"> —
      {{ symbols.length }} different gene names are in use across this one locus, the names to
      reconcile</template>. ⚠ The count beside each row is that <b>family's</b> gene count, not the
      gene name's: where the two differ the cell says how many of the family's genes are named at
      all. Names are counted exactly, with an allele tag (<span class="acc">ompK36_T333N</span>)
      folded onto the gene it belongs to first, so one gene cannot split its own vote. The
      clustering never saw UniRef50, Pfam, products or genomic context.
    </p>

    <p v-if="coverage" class="muted">{{ coverage }}</p>
  </div>
</template>
