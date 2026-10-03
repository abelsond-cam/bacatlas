import { fileURLToPath, URL } from "node:url";

import vue from "@vitejs/plugin-vue";
import { defineConfig } from "vitest/config";

export default defineConfig({
  // ⚠ The Vue plugin is needed HERE as well as in `vite.config.ts`: this config replaces that one
  // rather than extending it, and without the plugin a `.vue` import fails as "invalid JS syntax"
  // — which reads as a broken component rather than a missing build step.
  //
  // ⛔ And the `define` has to be repeated for the SAME reason. `__BUILD_STAMP__` is substituted by
  // `vite.config.ts`, which this file replaces, so without it the footer threw `ReferenceError:
  // __BUILD_STAMP__ is not defined` in four suites. A fixed value, not a wall clock: a test must not
  // depend on the second it ran in, and pinning it here means a dropped `define` in the real build
  // is still caught — by the build, which is where that failure belongs.
  define: { __BUILD_STAMP__: JSON.stringify("2026-01-01T00:00:00.000Z") },
  plugins: [vue()],
  resolve: { alias: { "@": fileURLToPath(new URL("./src", import.meta.url)) } },
  test: {
    globals: true,
    // `lib/`, `api/` and `stores/` are pure and need no DOM. The component suites declare `jsdom`
    // per-file with a docblock, so a pure test never pays for an environment it does not use.
    environment: "node",
    include: ["src/**/*.test.ts", "tests/**/*.test.ts"],
  },
});
