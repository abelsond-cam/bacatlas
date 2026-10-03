import { createPinia, setActivePinia } from "pinia";
import { beforeEach, describe, expect, it, vi } from "vitest";
import { nextTick } from "vue";

import { failure, success } from "@/api/result";
import type { CalibrationLadder, InferenceKind } from "@/api/types";
import { INFERENCE_KINDS, asCatalogueKey } from "@/api/types";
import type { FunctionResponse, LocusDetailResponse } from "@/api/types";

import { useFunctionBlockStore } from "./functionBlockStore";
import { useLocusNavigationStore } from "./locusNavigationStore";

const fetchLocus = vi.hoisted(() => vi.fn());
const fetchLocusFunction = vi.hoisted(() => vi.fn());
// ⚠ A PARTIAL mock: `anchorQuery` must stay real, because it is the one place the anchor's
// kind becomes a query parameter and a stub would hide a projected genome being fetched as
// `anchor=` — which 404s and reads on the page as "this genome has nothing here".
vi.mock("@/api/client", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/api/client")>()),
  fetchLocus, fetchLocusFunction,
}));

function locusDetail(label: string): LocusDetailResponse {
  return {
    locus: { label, gene_count: 100 } as LocusDetailResponse["locus"],
    annotations: {},
    uniref50_families: [],
    arrangements: {
      listed: [],
      total: 0,
      arrangements_not_listed: 0,
      members_in_arrangements_not_listed: 0,
      members_without_a_neighbourhood: 0,
      membership_is_complete: true,
    },
    anchor: { is_anchored: false, arrangement_ranks: [], kind: null, projected_copies: [] },
    offsets: [],
    intergenic_gaps: [],
    pfam_reference: {},
    neighbour_display_rows: [],
    resolved_neighbour_count: 0,
  };
}

/**
 * ⛔ The inference block is REQUIRED on a `FunctionResponse`, so a fixture must carry it. Empty by
 * default: a node that needs no suggestion is the ordinary case, and a fixture that silently
 * shipped a candidate would make every card under test look like an inferred one.
 */
function emptyInference(): Pick<FunctionResponse, "inference" | "calibration"> {
  return {
    inference: {
      representation: "esm",
      vocabularies: INFERENCE_KINDS.map((annotation_kind) => ({
        annotation_kind,
        own: null,
        walk: [],
        candidate: null,
      })),
    },
    calibration: {
      cog_orthogroup: ladder("cog_orthogroup"),
      ec_number: ladder("ec_number"),
    },
  };
}

function ladder(annotation_kind: InferenceKind): CalibrationLadder {
  return {
    annotation_kind,
    representation: "esm",
    annotated_locus_count: 0,
    locus_count: 0,
    min_pairs: 30,
    cells: [],
  };
}

function functionBlock(cogCount: number): FunctionResponse {
  return {
    annotations: { cog_orthogroup: [{ rank: 0, term: "COG1132", name: "x", gene_count: cogCount }] },
    ...emptyInference(),
    coverage: {
      gene_count: 100,
      cog_annotated_gene_count: cogCount,
      cog_distinct_id_count: 1,
      modal_cog_categories: ["N"],
      ec_annotated_gene_count: 0,
      kegg_annotated_gene_count: 0,
      go_annotated_gene_count: {
        molecular_function: 0,
        biological_process: 0,
        cellular_component: 0,
      },
    },
    go_verdicts: {
      molecular_function: "no_coverage",
      biological_process: "no_coverage",
      cellular_component: "no_coverage",
    },
  };
}

async function drawLocus(label: string) {
  const navigation = useLocusNavigationStore();
  navigation.setCatalogue(asCatalogueKey("ecoli"));
  fetchLocus.mockResolvedValueOnce(success(locusDetail(label)));
  await navigation.navigateTo(label);
  await nextTick();
}

beforeEach(() => {
  setActivePinia(createPinia());
  fetchLocus.mockReset();
  fetchLocusFunction.mockReset();
});

describe("⭐ fetched on tab OPEN, not on every walk", () => {
  it("asks for nothing until the reader opens the tab", async () => {
    await drawLocus("1");
    expect(fetchLocusFunction).not.toHaveBeenCalled();

    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(60)));
    useFunctionBlockStore().setOpen(true);
    await nextTick();
    expect(fetchLocusFunction).toHaveBeenCalledTimes(1);
  });

  it("⛔ does not re-ask when the tab is closed and opened again on the same locus", async () => {
    // The block is already in hand; a second request would buy nothing and spend a round trip on
    // every tab click.
    await drawLocus("1");
    const store = useFunctionBlockStore();
    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(60)));
    store.setOpen(true);
    await nextTick();
    await nextTick();
    store.setOpen(false);
    store.setOpen(true);
    await nextTick();
    expect(fetchLocusFunction).toHaveBeenCalledTimes(1);
  });

  it("⛔ setOpen is ABSOLUTE, never a toggle", () => {
    const store = useFunctionBlockStore();
    store.setOpen(true);
    store.setOpen(true);
    expect(store.isOpen).toBe(true);
  });
});

describe("⚠ a reader who leaves the tab open and walks sees the NEW locus", () => {
  it("discards the old block and asks again", async () => {
    await drawLocus("1");
    const store = useFunctionBlockStore();
    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(60)));
    store.setOpen(true);
    await nextTick();
    await nextTick();
    expect(store.block?.coverage.cog_annotated_gene_count).toBe(60);

    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(11)));
    await drawLocus("2");
    await nextTick();
    expect(fetchLocusFunction).toHaveBeenCalledTimes(2);
    expect(store.block?.coverage.cog_annotated_gene_count).toBe(11);
  });

  it("⛔ does NOT ask again when the tab is shut", async () => {
    await drawLocus("1");
    const store = useFunctionBlockStore();
    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(60)));
    store.setOpen(true);
    await nextTick();
    await nextTick();
    store.setOpen(false);

    await drawLocus("2");
    await nextTick();
    expect(fetchLocusFunction).toHaveBeenCalledTimes(1);
    // ⛔ And the previous locus's block is GONE rather than left showing under a new locus.
    expect(store.block).toBeNull();
  });
});

describe("⛔ a failure renders as a failure, never as 'this locus has no function'", () => {
  it("keeps the failure and holds no block", async () => {
    await drawLocus("1");
    const store = useFunctionBlockStore();
    fetchLocusFunction.mockResolvedValueOnce(failure("server", "503"));
    store.setOpen(true);
    await nextTick();
    await nextTick();
    expect(store.status).toBe("failed");
    expect(store.block).toBeNull();
    expect(store.lastFailure?.detail).toBe("503");
  });

  it("⚠ a retry after a failure DOES ask again — unlike a second tab click after a success", async () => {
    await drawLocus("1");
    const store = useFunctionBlockStore();
    fetchLocusFunction.mockResolvedValueOnce(failure("network", "offline"));
    store.setOpen(true);
    await nextTick();
    await nextTick();

    fetchLocusFunction.mockResolvedValueOnce(success(functionBlock(60)));
    await store.load();
    expect(fetchLocusFunction).toHaveBeenCalledTimes(2);
    expect(store.status).toBe("ready");
    expect(store.lastFailure).toBeNull();
  });
});
