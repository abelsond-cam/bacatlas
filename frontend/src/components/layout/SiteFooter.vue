<script setup lang="ts">
/**
 * The footer: what this catalogue is, how its model was built, what the audit counted against it,
 * and whose data it shows.
 *
 * ⭐ **Every number here is read from the model's own output and its audit — nothing re-derived**,
 * which is the promise its first sentence makes. The audit headline is QUOTED under the keys the
 * report used, in a closed `<details>`: a reader who wants to know what was counted against this
 * model must be able to find it, but a stranger landing here wants a locus, not a scorecard.
 *
 * ⚠ **Attribution is binding**, not a courtesy: GO and UniProt are CC BY 4.0. KEGG appears only as a
 * linked accession — its terms permit linking and not redistribution, which is why no KO is named.
 */
import { computed } from "vue";

import type { SpeciesCatalogueResponse } from "@/api/types";
import AuditResidualLists from "@/components/footer/AuditResidualLists.vue";

const props = defineProps<{ catalogue: SpeciesCatalogueResponse }>();
const emit = defineEmits<{ go: [locusLabel: string] }>();

/**
 * The graded headline in the order a reader needs it (`render_page._HEADLINE_ROWS`): how big, what was
 * found against it, what was rescued. Each is formatted exactly as the published footer formatted it.
 */
const HEADLINE_ROWS: readonly (readonly [string, string, (value: number) => string])[] = [
  ["clusters graded", "n_clusters_total", (value) => value.toLocaleString()],
  ["genes", "n_genes_total", (value) => value.toLocaleString()],
  ["synteny only — no homology by mmseqs", "synteny_only_n_clusters", (value) => `${value.toLocaleString()} clusters`],
  ["synteny only, as a gene rate", "synteny_only_gene_rate", (value) => `${(value * 100).toFixed(4)}%`],
  ["Pfam conflict", "pfam_conflict_n_clusters", (value) => `${value.toLocaleString()} clusters`],
  // NOT "of which": the denominator is different — Pfam judges only loci with ≥ 2 annotated members.
  ["loci Pfam can judge (≥2 annotated members)", "pfam_judgeable_n_clusters", (value) => `${value.toLocaleString()} clusters`],
  ["UniRef50 family split across loci", "split_gene_rate_excl_singletons", (value) => `${(value * 100).toFixed(2)}% of genes`],
  ["rescued by ESM homology", "n_clusters_esm_rescued", (value) => `${value.toLocaleString()} clusters`],
];

const provenance = computed(() => {
  const pangenome = props.catalogue.pangenome;
  const rows: [string, string][] = props.catalogue.provenance_rows.map(([label, value]) => [label, value]);
  rows.push(["model", pangenome.run_id]);
  // ⛔ `built_at` is NULL on every catalogue loaded so far — nothing in `ingest` writes it — so this
  // row read "built —" and the footer could not answer the one question it exists to answer. The
  // catalogue's `ingested_at` IS populated and is the honest version of it: when this stack loaded
  // the data it is serving. `built_at` still wins where it exists, because that is the stronger fact.
  const catalogueDate = pangenome.built_at ?? shortDate(pangenome.ingested_at);
  rows.push([pangenome.built_at ? "catalogue built" : "catalogue loaded", catalogueDate ?? "—"]);
  if (pangenome.git_sha) rows.push(["model code", pangenome.git_sha]);
  // ⭐ And when THIS bundle was built. The compose stack serves a built image, so "the app is up" and
  // "the app is current" are different facts, and an image from last week looks identical to a fresh
  // one. Printing both dates is what makes "is this current?" answerable from the page.
  rows.push(["page built", shortDate(__BUILD_STAMP__) ?? "—"]);
  return rows;
});

/** `2026-10-03T08:35:06.477248+00:00` → `2026-10-03 08:35 UTC`; `null` stays `null`. */
function shortDate(value: string | null | undefined): string | null {
  if (!value) return null;
  const parsed = new Date(value);
  if (Number.isNaN(parsed.getTime())) return value; // ⚠ Say what we were given rather than "Invalid Date".
  return `${parsed.toISOString().slice(0, 16).replace("T", " ")} UTC`;
}

const headline = computed(() =>
  HEADLINE_ROWS.flatMap(([label, key, format]) => {
    const value = props.catalogue.audit_headline[key];
    if (value === undefined || value === null) return [];
    return [[label, typeof value === "number" ? format(value) : String(value)] as const];
  }),
);

const omitted = computed(() => {
  const reasons = Object.values(props.catalogue.pangenome.omitted_sections);
  return reasons.length ? `Not shown on this build: ${reasons.join("; ")}.` : "";
});
</script>

<template>
  <footer class="foot">
    <div class="wrap">
      <p>
        Prototype built from a single <em>{{ catalogue.species.scientific_name }}</em> clustering. Every
        number on this page is read from that model's own output and its accessory-fidelity audit — nothing
        is re-derived for display.
      </p>
      <!--
        ⛔⛔ **The licence flag, site-wide, and it is a GATE on deployment rather than a courtesy.**
        David, 2026-10-05: *"Don't worry about KEGG licence! Just flag it on site. Will look at it
        later. This is a prototype only! On my own device right now!"* KEGG's terms require an
        academic service provider licence from anyone offering services, and serving a KO
        description is offering one. This paragraph is what carries the unanswered question forward
        so it cannot quietly become a public deployment; it is pinned by a test for that reason.
        EC needs no flag — CC BY 4.0 — and is attributed on the card that uses it.
      -->
      <p class="foot-notice">
        <strong>Internal build — not for publication.</strong> It names KEGG orthologies, and KEGG
        requires an academic service provider licence from anyone offering a service. That licence is
        <strong>not yet confirmed</strong>, so this build must not be served publicly until it is.
      </p>
      <dl>
        <div v-for="([label, value], index) in provenance" :key="index">
          <dt>{{ label }}</dt>
          <dd>{{ value }}</dd>
        </div>
      </dl>
      <details v-if="headline.length" class="foot-detail">
        <summary>Accessory-fidelity audit — what was counted against this model</summary>
        <dl>
          <div v-for="[label, value] in headline" :key="label">
            <dt>{{ label }}</dt>
            <dd>{{ value }}</dd>
          </div>
        </dl>
        <p class="muted">
          Read verbatim from the accessory audit. <b>Synteny only</b> means no homology mmseqs, Pfam or ESM
          could find — the members were grouped on genomic context. That names the evidence, not a mistake:
          many are good calls, and the way to tell is the evidence beside each one. Clusters rescued by Pfam
          or by ESM are homologous and are not counted here.
        </p>
      </details>
      <!-- ⛔ The CATALOGUE's key below, not the species'. This block's whole claim is "a model's own
           residuals, on its own page" — handed the species key it would fetch the DEFAULT model's
           residuals and print them under another model's name. -->
      <AuditResidualLists
        :catalogue-key="catalogue.pangenome.catalogue_key"
        :audit="catalogue.audit_headline"
        @go="emit('go', $event)"
      />
      <p v-if="omitted" class="muted">{{ omitted }}</p>
      <p class="muted">
        Reference data: gene and protein annotation by <b>Bakta</b>; protein families from
        <a href="https://www.uniprot.org/">UniProt / UniRef50</a> and
        <a href="https://www.ebi.ac.uk/interpro/">Pfam / InterPro</a>; functional classification from the
        <a href="https://geneontology.org/">Gene Ontology</a> and
        <a href="https://www.ncbi.nlm.nih.gov/research/cog">NCBI COG</a>;
        <a href="https://www.genome.jp/kegg/">KEGG</a> orthology accessions are linked, not reproduced. Gene
        Ontology and UniProt data are used under
        <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>.
      </p>
    </div>
  </footer>
</template>
