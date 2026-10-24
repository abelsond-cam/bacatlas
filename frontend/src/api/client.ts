/**
 * One function per endpoint. **No component builds a URL**, so a route change is one edit here.
 *
 * ⭐ Every response is immutable for a build — writes are ours alone and offline — so the browser's
 * own HTTP cache does the work, and the client adds no invalidation logic because nothing can
 * invalidate. That is the single largest payoff of the read-only scope, and it is why prefetching
 * on `pointerenter` is free rather than a risk.
 */

import { requestJson } from "./request";
import type { Result } from "./result";
import type {
  ArrangementPageResponse,
  AuditResidualsResponse,
  CatalogueKey,
  CataloguesResponse,
  FunctionResponse,
  GenomeListResponse,
  GenomeSourceGroup,
  GeneSequenceResponse,
  LocusDetailResponse,
  ProjectedGenomesResponse,
  SearchResponse,
  SpeciesCatalogueResponse,
  SpeciesListResponse,
} from "./types";

/**
 * ⛔ A locus label goes into a PATH segment, so it is encoded here rather than interpolated. The
 * API routes it with `<path:locus_label>`, which happily accepts a `/` and would otherwise turn one
 * locus into a 404 on a route nobody wrote.
 */
function labelSegment(locusLabel: string): string {
  return encodeURIComponent(locusLabel);
}

export function fetchSpeciesList(signal?: AbortSignal): Promise<Result<SpeciesListResponse>> {
  return requestJson<SpeciesListResponse>("species", signal ? { signal } : {});
}

/**
 * Every catalogue the server OFFERS — the model picker's source.
 *
 * ⚠ One row per catalogue, so a species with two models appears twice. A catalogue that is loaded
 * but staged is absent here and still loads by key, so never gate a fetch on this list.
 */
export function fetchCatalogues(signal?: AbortSignal): Promise<Result<CataloguesResponse>> {
  return requestJson<CataloguesResponse>("catalogues", signal ? { signal } : {});
}

export function fetchSpeciesCatalogue(
  catalogueKey: CatalogueKey,
  signal?: AbortSignal,
): Promise<Result<SpeciesCatalogueResponse>> {
  return requestJson<SpeciesCatalogueResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}`,
    signal ? { signal } : {},
  );
}

/**
 * ⭐ The hot path. One round trip, and the popover is then **offline**: all ten offset slots, the
 * drawn arrangements, the gaps and both representations' similarity are in this one response, so
 * opening a popover, switching arrangement, flipping walk direction and switching the similarity
 * card's view all fetch nothing.
 *
 * `anchorSampleId` names a genome; given one, the arrangement that genome carries is included even
 * if it sits past the display cap — otherwise the reader is told in words that their genome sits in
 * #37 and has no button to go back to it.
 */
export function fetchLocus(
  catalogueKey: CatalogueKey,
  locusLabel: string,
  options: { anchorSampleId?: string; projectedSampleId?: string; signal?: AbortSignal } = {},
): Promise<Result<LocusDetailResponse>> {
  // ⛔ `anchor` and `projected` are different parameters because they are different relations: a
  // modelled genome is COUNTED in an arrangement, a projected one MATCHES it. The API refuses both
  // at once with a named 400 rather than picking one silently, so this never sends both.
  const query = options.anchorSampleId
    ? { anchor: options.anchorSampleId }
    : options.projectedSampleId
      ? { projected: options.projectedSampleId }
      : undefined;
  return requestJson<LocusDetailResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/loci/${labelSegment(locusLabel)}`,
    {
      ...(options.signal ? { signal: options.signal } : {}),
      ...(query ? { query } : {}),
    },
  );
}

/**
 * The one place the anchor's KIND becomes a query parameter — `{}` for no anchor.
 *
 * ⚠ Everything that fetches a locus goes through this, so a new call site cannot forget that a
 * projected genome is not fetched with `anchor=`. Sending it that way returns a 404 (it is in no
 * collection), which renders as "this genome has nothing here" — a claim, and a false one.
 */
export function anchorQuery(
  sampleId: string | null,
  kind: "catalogue" | "projected",
): { anchorSampleId?: string; projectedSampleId?: string } {
  if (!sampleId) return {};
  return kind === "projected" ? { projectedSampleId: sampleId } : { anchorSampleId: sampleId };
}

/** Genomes placed on this catalogue after the model was built — never members of it. */
export function fetchProjectedGenomes(
  catalogueKey: CatalogueKey,
  signal?: AbortSignal,
): Promise<Result<ProjectedGenomesResponse>> {
  return requestJson<ProjectedGenomesResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/projected-genomes`,
    { ...(signal ? { signal } : {}) },
  );
}

/** Arrangements past the display cut — the full scroller, paged. */
export function fetchArrangementPage(
  catalogueKey: CatalogueKey,
  locusLabel: string,
  offset: number,
  signal?: AbortSignal,
): Promise<Result<ArrangementPageResponse>> {
  return requestJson<ArrangementPageResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/loci/${labelSegment(locusLabel)}/arrangements`,
    { query: { offset }, ...(signal ? { signal } : {}) },
  );
}

/** The EggNOG tab — fetched on tab open, not on every walk. */
export function fetchLocusFunction(
  catalogueKey: CatalogueKey,
  locusLabel: string,
  signal?: AbortSignal,
): Promise<Result<FunctionResponse>> {
  return requestJson<FunctionResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/loci/${labelSegment(locusLabel)}/function`,
    signal ? { signal } : {},
  );
}

/**
 * The Sequence tab — one genome's gene(s) at one locus, sliced from the original GFF on the server.
 *
 * ⚠ Fetched on tab open and per genome, never on a walk: it is the only endpoint here that opens a
 * file, and nobody wants a genome's bases on every one of forty steps.
 */
export function fetchGeneSequence(
  catalogueKey: CatalogueKey,
  sampleId: string,
  locusLabel: string,
  signal?: AbortSignal,
): Promise<Result<GeneSequenceResponse>> {
  return requestJson<GeneSequenceResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/genomes/${encodeURIComponent(sampleId)}/loci/` +
      `${labelSegment(locusLabel)}/sequence`,
    signal ? { signal } : {},
  );
}

/**
 * Substring search with the page's exact semantics — `pg_trgm` accelerating `ILIKE '%q%'`, so
 * `ligase` still finds *O-antigen ligase RfaL* mid-string. Replaces the resident 3.0 MB haystack.
 */
export function searchLoci(
  catalogueKey: CatalogueKey,
  query: string,
  options: { limit?: number; signal?: AbortSignal } = {},
): Promise<Result<SearchResponse>> {
  return requestJson<SearchResponse>(`catalogues/${encodeURIComponent(catalogueKey)}/search`, {
    query: { q: query, limit: options.limit },
    ...(options.signal ? { signal: options.signal } : {}),
  });
}

/**
 * The genomes an anchor can name, filtered server-side — the published page held all of them in
 * `meta.genomes` and filtered in the browser, which at 80,000 genomes is an array nobody should ship.
 * Case-insensitive substring on the sample id, the same semantics as `app.js::anchorSearch`.
 */
export function fetchGenomes(
  catalogueKey: CatalogueKey,
  query: string,
  options: { limit?: number; sourceGroup?: GenomeSourceGroup | null; signal?: AbortSignal } = {},
): Promise<Result<GenomeListResponse>> {
  return requestJson<GenomeListResponse>(`catalogues/${encodeURIComponent(catalogueKey)}/genomes`, {
    // ⚠ `source_group` is omitted when null rather than sent empty: the server reads an absent
    // parameter as "every group", and an unknown VALUE as a named 400. Sending `undefined` here
    // keeps those two cases apart instead of relying on the empty-string carve-out.
    query: { q: query, limit: options.limit, source_group: options.sourceGroup ?? undefined },
    ...(options.signal ? { signal: options.signal } : {}),
  });
}

/** The footer's two residual lists — fetched when the reader opens one, never with the page. */
export function fetchAuditResiduals(
  catalogueKey: CatalogueKey,
  signal?: AbortSignal,
): Promise<Result<AuditResidualsResponse>> {
  return requestJson<AuditResidualsResponse>(
    `catalogues/${encodeURIComponent(catalogueKey)}/audit/residual-loci`,
    signal ? { signal } : {},
  );
}
