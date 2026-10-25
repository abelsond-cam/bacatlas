/**
 * The **one** async boundary, and the stale-response guard generalised from the page's `seqTok`.
 *
 * `app.js:4118` states the invariant this replaces: *"app.js is otherwise wholly synchronous by
 * choice… so the asynchrony is confined to `withSeq` and never leaks upward."* It cannot stay that
 * way once the catalogue lives in a database, so the discipline moves here instead:
 *
 * ⭐ **Exactly three things fetch**, and each gets its own lane: navigating to a locus, opening the
 * function tab, opening the sequence tab. Everything else the reader can do — switching arrangement,
 * flipping walk direction, opening a popover, changing map representation — is answered out of the
 * response already in hand. If a fourth lane appears, something has been moved off the hot path
 * that belongs on it.
 *
 * ⛔ **A superseded response is not a failure and not an empty result — it is a third thing**, and
 * {@link LaneOutcome} makes a component say so. A reader who walks on twice quickly supersedes their
 * own first request; rendering an error for that, or rendering nothing, are both wrong. The right
 * behaviour is to change nothing at all, and the early return that does it is visible in the code.
 */

import { type Failure, type Result, failure, success } from "./result";

/**
 * ⛔ Same origin by default. A literal origin in source is what makes the next deploy a rewrite.
 *
 * ⚠ `||`, not `??`: a build run with `VITE_API_BASE_URL=` (set, but empty — the natural way to "leave
 * it unset" in a Dockerfile ARG or a CI matrix) must fall back too. With `??` the empty string
 * survived, every request went to `/species/…` at the origin root, and each one 404'd.
 */
export const API_BASE_URL: string = import.meta.env.VITE_API_BASE_URL || "/api/v1";

/** Join the base with a path, tolerating a trailing slash on either side. */
export function apiUrl(path: string): string {
  const base = API_BASE_URL.replace(/\/$/, "");
  return `${base}/${path.replace(/^\//, "")}`;
}

/**
 * ⛔⛔ **The deploy stamp every API request carries, and why a cached page is otherwise a LIE.**
 *
 * Catalogue responses are served `Cache-Control: public, max-age=86400` with
 * `ETag: "{pangenome_id}"`, on the stated reasoning that *for one URI the body changes only when
 * the pangenome does*. That premise holds for the DATA and fails for the SHAPE: the 2026-10-25
 * deploy added `genome_source_groups` to the shell without the pangenome changing, so every
 * browser that had loaded the site in the previous 24 h kept replaying the old body — no new
 * field, therefore no anchor tabs — and `max-age` meant it never even asked. New code reading a
 * day-old response is the worst version of this: both halves are individually correct and the page
 * is silently a deploy behind, with nothing to see in any log.
 *
 * Appending the build stamp makes each deploy's requests a **different URI**, so a stale entry is
 * bypassed rather than revalidated — no hard reload, no reduced caching, and the existing ETag
 * reasoning is untouched because the query string is already part of the target URI (RFC 9110
 * §8.8.3), which is exactly how `?q=` is already treated.
 *
 * ⚠ It must be on EVERY request, not just the shell: any response whose shape a future deploy
 * changes has the same problem, and remembering which is the failure mode this replaces.
 */
/*
 * ⚠ Digits only. `__BUILD_STAMP__` is an ISO timestamp, and its colons percent-encode into
 * `%3A` on every request — correct but unreadable in a log and longer than it needs to be. The
 * digits change on every build just as well, which is the whole requirement. The human-readable
 * form stays in the footer, which uses `__BUILD_STAMP__` directly.
 */
export const BUILD_STAMP: string = __BUILD_STAMP__.replace(/\D/g, "");

/**
 * One GET, parsed, as a {@link Result}.
 *
 * ⚠ Every failure mode is *named* rather than collapsed into a boolean: a 404 is the server
 * answering, and a reader can act on it; a network error is not, and retrying may work.
 */
export async function requestJson<T>(
  path: string,
  options: { signal?: AbortSignal; query?: Record<string, string | number | undefined> } = {},
): Promise<Result<T>> {
  const url = new URL(apiUrl(path), globalThis.location?.origin ?? "http://localhost");
  for (const [key, value] of Object.entries(options.query ?? {})) {
    if (value !== undefined) url.searchParams.set(key, String(value));
  }
  // ⛔ See BUILD_STAMP: this is what stops a browser replaying the previous deploy's response shape.
  url.searchParams.set("v", BUILD_STAMP);
  // Relative when it can be — the URL object above exists to build the query string safely, not to
  // pin an origin into the request.
  const target = url.origin === (globalThis.location?.origin ?? "http://localhost")
    ? `${url.pathname}${url.search}`
    : url.toString();

  let response: Response;
  try {
    // `exactOptionalPropertyTypes` refuses `signal: undefined` here, and it is right to: an
    // absent key and a key holding `undefined` are different things to `RequestInit`.
    response = await fetch(target, {
      headers: { Accept: "application/json" },
      ...(options.signal ? { signal: options.signal } : {}),
    });
  } catch (error) {
    if (error instanceof Error && error.name === "AbortError") throw error;
    return failure("network", describe(error, "the request did not reach the server"));
  }

  let body: unknown;
  try {
    body = await response.json();
  } catch (error) {
    if (error instanceof Error && error.name === "AbortError") throw error;
    // ⚠ A non-JSON ERROR page is the server answering with a status, not a contract breach: a route
    // that does not exist gets the framework's HTML 404, and reporting that as "Unexpected token '<'"
    // tells a reader nothing and hides that it was a 404 at all. Only a SUCCESS that is not JSON is
    // `malformed`, because only that one is our bug.
    if (!response.ok) {
      return failure(
        response.status === 404 ? "not_found" : "server",
        `the server answered ${response.status}`,
        response.status,
      );
    }
    return failure(
      "malformed",
      describe(error, "the server's answer was not JSON"),
      response.status,
    );
  }

  if (response.ok) return success(body as T);

  // The API answers a 404 with `{error, detail}` naming which of the species, the pangenome or the
  // locus was missing — so show the server's own sentence rather than inventing one.
  const detail =
    isRecord(body) && typeof body["detail"] === "string"
      ? body["detail"]
      : `the server answered ${response.status}`;
  return failure(response.status === 404 ? "not_found" : "server", detail, response.status);
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null;
}

function describe(error: unknown, fallback: string): string {
  return error instanceof Error && error.message ? error.message : fallback;
}

/**
 * What a lane hands back. Superseded is its own case, so it cannot be confused with either an
 * error or an empty answer.
 */
export type LaneOutcome<T> =
  | { readonly superseded: true }
  | { readonly superseded: false; readonly result: Result<T> };

const SUPERSEDED: LaneOutcome<never> = { superseded: true };

/**
 * A single in-flight request, with the one that came before it cancelled.
 *
 * ⚠ The token is checked **after** the await as well as relying on the abort: an abort is a request
 * to stop, not a guarantee the response has not already resolved. Both guards are needed, and the
 * token is the one that is actually sufficient.
 */
export class RequestLane {
  readonly name: string;
  private token = 0;
  private controller: AbortController | null = null;

  constructor(name: string) {
    this.name = name;
  }

  /** Whether a request is in flight. Drives the `refreshing` state, never the `pending` one. */
  get isBusy(): boolean {
    return this.controller !== null;
  }

  async run<T>(work: (signal: AbortSignal) => Promise<Result<T>>): Promise<LaneOutcome<T>> {
    this.controller?.abort();
    const controller = new AbortController();
    this.controller = controller;
    const issued = ++this.token;
    try {
      const result = await work(controller.signal);
      if (issued !== this.token) return SUPERSEDED;
      return { superseded: false, result };
    } catch (error) {
      if (issued !== this.token) return SUPERSEDED;
      if (error instanceof Error && error.name === "AbortError") return SUPERSEDED;
      return {
        superseded: false,
        result: failure("network", describe(error, "the request failed")) satisfies Failure,
      };
    } finally {
      if (issued === this.token) this.controller = null;
    }
  }

  /** Abandon whatever is in flight — leaving the view as it is, not as an error. */
  cancel(): void {
    this.controller?.abort();
    this.controller = null;
    this.token++;
  }
}

/**
 * ⛔ **The lanes, named here rather than created ad hoc so that adding one is a visible edit to
 * this list and has to be argued for.** Three of them are the hot path's own rule — navigating,
 * opening the EggNOG tab, opening the Sequence tab — and the fourth is argued below.
 *
 * ⭐ **`arrangements` — the argument, made once, so the next lane has to clear the same bar.**
 * The three-lane rule exists to keep the *hot path* one round trip: everything a reader does to
 * the drawn locus — switching arrangement, flipping the frame, opening a popover, changing map
 * representation — is answered from the response already in hand. Paging past the arrangement cut
 * is not on that path. It is a deliberate second act by a reader who has opened the A0 card, and
 * it is the fourth fetch the design named from the start.
 *
 * It cannot ride either of the other lanes, and the reason is supersession rather than tidiness:
 * `RequestLane.run` **cancels whatever that lane is holding**. On `navigation`, asking for the next
 * page of arrangements would abort the reader's own walk; on `function`, opening the EggNOG tab
 * while the scroller loads would abort one of two independent surfaces. Lanes are not a namespace —
 * they are a cancellation domain, and these are three domains.
 *
 * ⚠ Navigating away must still discard an arrangement page in flight, and the lane cannot do that
 * on its own. `arrangementBrowserStore` cancels this lane when the drawn locus changes, and
 * re-checks the label after the await, because a Vue watcher runs after the microtask that resolves
 * the fetch.
 */
export const LANES = {
  navigation: new RequestLane("navigation"),
  function: new RequestLane("function"),
  sequence: new RequestLane("sequence"),
  arrangements: new RequestLane("arrangements"),
} as const;
