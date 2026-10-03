/// <reference types="vite/client" />

interface ImportMetaEnv {
  /**
   * Where the API lives. ⛔ Defaults to the **same origin** at `/api/v1` — never an absolute URL in
   * source. The institutional host comes first, probably under a subpath; the public origin later.
   */
  readonly VITE_API_BASE_URL?: string;
  /** Where the app itself is served from. Consumed by `vite.config.ts`, not by application code. */
  readonly VITE_PUBLIC_BASE?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

/**
 * When this bundle was built — substituted by `vite.config.ts`'s `define`, so it is a literal in the
 * output rather than a value read at run time. ⛔ Declared here because a `define` is invisible to
 * TypeScript otherwise: the identifier compiles, ships, and throws `ReferenceError` in the browser if
 * the `define` is ever dropped. With this declaration, dropping it is a build failure.
 */
declare const __BUILD_STAMP__: string;
