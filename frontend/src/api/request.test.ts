import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import { API_BASE_URL, BUILD_STAMP, RequestLane, apiUrl, requestJson } from "./request";
import { failure, mapResult, success } from "./result";
import * as resultModule from "./result";

/** A `fetch` stand-in that resolves a JSON body with a chosen status. */
function respondWith(status: number, body: unknown, options: { json?: () => never } = {}) {
  return vi.fn(async () => ({
    ok: status >= 200 && status < 300,
    status,
    json: options.json ?? (async () => body),
  })) as unknown as typeof fetch;
}

const originalFetch = globalThis.fetch;
afterEach(() => {
  globalThis.fetch = originalFetch;
  vi.restoreAllMocks();
});

describe("the base URL", () => {
  it("⛔ is same-origin and relative by default — never an absolute URL in source", () => {
    // The institutional host comes first, probably under a subpath; the public origin later. Both
    // are build-time config, and this is the assertion that keeps a literal out of the bundle.
    expect(API_BASE_URL.startsWith("http")).toBe(false);
    expect(API_BASE_URL).toBe("/api/v1");
  });

  it("⚠ treats a SET-but-EMPTY base as unset, rather than sending requests to the origin root", async () => {
    // `VITE_API_BASE_URL=` is how a Dockerfile ARG or a CI matrix "leaves it out". With `??` the empty
    // string survived and every request went to `/species/…`.
    vi.stubEnv("VITE_API_BASE_URL", "");
    vi.resetModules();
    const reloaded = await import("./request");
    expect(reloaded.API_BASE_URL).toBe("/api/v1");
    vi.unstubAllEnvs();
    vi.resetModules();
  });

  it("joins without doubling or dropping a slash", () => {
    expect(apiUrl("species")).toBe("/api/v1/species");
    expect(apiUrl("/species")).toBe("/api/v1/species");
  });
});

describe("requestJson names every failure mode rather than collapsing them", () => {
  it("returns the parsed body on 200", async () => {
    globalThis.fetch = respondWith(200, { hello: "world" });
    const result = await requestJson<{ hello: string }>("thing");
    expect(result).toEqual(success({ hello: "world" }));
  });

  it("⛔ a 404 is `not_found` and carries the SERVER's sentence, not one we invented", async () => {
    // The API answers a 404 naming which of the species, the pangenome or the locus was missing.
    globalThis.fetch = respondWith(404, { error: "not_found", detail: "no locus '9999'" });
    const result = await requestJson("thing");
    expect(result).toEqual(failure("not_found", "no locus '9999'", 404));
  });

  it("a 500 is `server`, which is a different sentence to a reader than a 404", async () => {
    globalThis.fetch = respondWith(500, {});
    const result = await requestJson("thing");
    expect(result.ok).toBe(false);
    if (!result.ok) {
      expect(result.kind).toBe("server");
      expect(result.status).toBe(500);
    }
  });

  it("an unreachable server is `network` — retrying may work, which 500 does not imply", async () => {
    globalThis.fetch = vi.fn(async () => {
      throw new TypeError("Failed to fetch");
    }) as unknown as typeof fetch;
    const result = await requestJson("thing");
    expect(result.ok).toBe(false);
    if (!result.ok) expect(result.kind).toBe("network");
  });

  it("⚠ a body that is not JSON is `malformed` — OUR bug, and it says so", async () => {
    globalThis.fetch = vi.fn(async () => ({
      ok: true,
      status: 200,
      json: async () => {
        throw new SyntaxError("Unexpected token <");
      },
    })) as unknown as typeof fetch;
    const result = await requestJson("thing");
    expect(result.ok).toBe(false);
    if (!result.ok) expect(result.kind).toBe("malformed");
  });

  it("⚠ an HTML ERROR page is its STATUS, not a parse error — a missing route reads as a 404", async () => {
    // Flask answers an unknown route with an HTML 404. Reported as `malformed`, a reader saw
    // "Unexpected token '<'" in the anchor list and nothing said it was a 404 at all.
    globalThis.fetch = vi.fn(async () => ({
      ok: false,
      status: 404,
      json: async () => {
        throw new SyntaxError("Unexpected token '<', \"<!doctype \"... is not valid JSON");
      },
    })) as unknown as typeof fetch;
    const result = await requestJson("thing");
    expect(result.ok).toBe(false);
    if (!result.ok) {
      expect(result.kind).toBe("not_found");
      expect(result.detail).toBe("the server answered 404");
    }
  });

  it("puts query parameters on the URL and omits the undefined ones", async () => {
    const stub = respondWith(200, {});
    globalThis.fetch = stub;
    await requestJson("thing", { query: { offset: 50, limit: undefined } });
    const [url] = (stub as unknown as ReturnType<typeof vi.fn>).mock.calls[0]!;
    expect(url).toBe(`/api/v1/thing?offset=50&v=${BUILD_STAMP}`);
  });

  it("⛔ carries the BUILD STAMP on every request, even with no query at all", async () => {
    // ⛔⛔ The regression this exists for is a page that is silently one deploy behind. Catalogue
    // responses are cached `max-age=86400` with an ETag that was the pangenome id alone, so a
    // release that changed a response's SHAPE left every warm browser replaying the old body and
    // never revalidating. New bundle, day-old JSON, no error anywhere — the 2026-10-25 anchor tabs
    // simply did not appear. The stamp makes each deploy's requests a different URI, so the stale
    // entry is bypassed rather than trusted.
    const stub = respondWith(200, {});
    globalThis.fetch = stub;
    await requestJson("species");
    const [url] = (stub as unknown as ReturnType<typeof vi.fn>).mock.calls[0]!;
    expect(url).toBe(`/api/v1/species?v=${BUILD_STAMP}`);
    expect(BUILD_STAMP).toBeTruthy();
  });
});

describe("⛔ the Result type has no escape hatch", () => {
  it("exports no unwrap, no getOrElse and no default", () => {
    // `T | null` is the shape that permits the bug, because null is exactly as falsy as an empty
    // list. A `getOrElse` would be that bug in one function call, so the surface is asserted.
    expect(Object.keys(resultModule).sort()).toEqual([
      "SupersededError",
      "failure",
      "mapResult",
      "success",
    ]);
  });

  it("mapResult leaves a failure alone", () => {
    const bad = failure("network", "nope");
    expect(mapResult(bad, () => "transformed")).toBe(bad);
    expect(mapResult(success(2), (value) => value * 3)).toEqual(success(6));
  });
});

describe("a request lane", () => {
  let lane: RequestLane;
  beforeEach(() => {
    lane = new RequestLane("test");
  });

  it("returns the result of a single request", async () => {
    const outcome = await lane.run(async () => success("value"));
    expect(outcome).toEqual({ superseded: false, result: success("value") });
  });

  it("⭐ a superseded request is neither an error NOR an empty result", async () => {
    // A reader who walks on twice quickly supersedes their own first request. Rendering an error
    // for that, or rendering nothing, are both wrong: the right behaviour is to change nothing,
    // and `superseded` is what lets a component say so in one visible early return.
    let releaseFirst: (value: unknown) => void = () => {};
    const first = lane.run(async () => {
      await new Promise((resolve) => {
        releaseFirst = resolve;
      });
      return success("stale");
    });
    const second = lane.run(async () => success("fresh"));
    releaseFirst(null);

    expect(await first).toEqual({ superseded: true });
    expect(await second).toEqual({ superseded: false, result: success("fresh") });
  });

  it("aborts the previous request when a new one starts", async () => {
    const signals: AbortSignal[] = [];
    const held = lane.run(async (signal) => {
      signals.push(signal);
      await new Promise((resolve) => setTimeout(resolve, 5));
      return success("first");
    });
    await lane.run(async (signal) => {
      signals.push(signal);
      return success("second");
    });
    await held;
    expect(signals[0]!.aborted).toBe(true);
    expect(signals[1]!.aborted).toBe(false);
  });

  it("⚠ checks the TOKEN as well as the abort, because an abort is a request and not a guarantee", async () => {
    // A response can already have resolved when the abort fires. The token is the guard that is
    // actually sufficient; the abort is the one that saves the bytes.
    let releaseFirst: (value: unknown) => void = () => {};
    const first = lane.run(async (signal) => {
      await new Promise((resolve) => {
        releaseFirst = resolve;
      });
      // Resolves successfully DESPITE having been aborted.
      expect(signal.aborted).toBe(true);
      return success("stale");
    });
    await lane.run(async () => success("fresh"));
    releaseFirst(null);
    expect(await first).toEqual({ superseded: true });
  });

  it("an AbortError from fetch reads as superseded, not as a failure the reader caused", async () => {
    const aborted = new Error("aborted");
    aborted.name = "AbortError";
    const outcome = await lane.run(async () => {
      throw aborted;
    });
    expect(outcome).toEqual({ superseded: true });
  });

  it("a thrown error that is NOT an abort is a real failure", async () => {
    const outcome = await lane.run(async () => {
      throw new TypeError("boom");
    });
    expect(outcome.superseded).toBe(false);
    if (!outcome.superseded) {
      expect(outcome.result.ok).toBe(false);
      if (!outcome.result.ok) expect(outcome.result.kind).toBe("network");
    }
  });

  it("cancel abandons what is in flight and leaves the view alone", async () => {
    let release: (value: unknown) => void = () => {};
    const held = lane.run(async () => {
      await new Promise((resolve) => {
        release = resolve;
      });
      return success("stale");
    });
    lane.cancel();
    release(null);
    expect(await held).toEqual({ superseded: true });
    expect(lane.isBusy).toBe(false);
  });

  it("reports busy only while a request is in flight", async () => {
    expect(lane.isBusy).toBe(false);
    let release: (value: unknown) => void = () => {};
    const held = lane.run(async () => {
      await new Promise((resolve) => {
        release = resolve;
      });
      return success("done");
    });
    expect(lane.isBusy).toBe(true);
    release(null);
    await held;
    expect(lane.isBusy).toBe(false);
  });
});

describe("⛔ the lanes are a closed list, and adding one is a visible edit", () => {
  it("names every one of them", async () => {
    // Three are the hot path's own rule — navigating, the function tab, the sequence tab. The
    // fourth, `arrangements`, pages the A0 card past the display cut: not on the hot path, and on
    // its own lane because a lane is a CANCELLATION domain, so sharing one would abort the reader's
    // walk to fetch a scroller. The argument is beside `LANES`; this test is what makes a fifth
    // impossible to add quietly.
    const { LANES } = await import("./request");
    expect(Object.keys(LANES).sort()).toEqual([
      "arrangements",
      "function",
      "navigation",
      "sequence",
    ]);
  });

  it("⛔ keeps them independent — one lane's cancel does not touch another's", async () => {
    // The whole reason `arrangements` is not folded into `navigation`: `run` aborts whatever its
    // own lane is holding, so a shared lane would make paging the card cancel the walk.
    const { LANES } = await import("./request");
    let release: (value: unknown) => void = () => {};
    const walking = LANES.navigation.run(
      () => new Promise<never>((resolve) => { release = resolve as never; }),
    );
    LANES.arrangements.cancel();
    expect(LANES.navigation.isBusy).toBe(true);
    release(success(1));
    await walking;
  });
});
