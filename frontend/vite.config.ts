import { fileURLToPath, URL } from "node:url";

import vue from "@vitejs/plugin-vue";
import { defineConfig } from "vite";

// ⛔ NO absolute origin and NO hardcoded base path anywhere in the app. The service is
// institution-only first — probably served under a subpath — and public later, from a different
// origin. Both come from build-time config so that neither is a code change:
//   VITE_API_BASE_URL   where the API lives   (default: same origin, `/api/v1`)
//   VITE_PUBLIC_BASE    where the app is served from (default: `/`)
// A literal `https://…` or `/bacatlas/…` compiled into a component is exactly the thing that
// makes the second deployment a rewrite instead of an environment variable.
// ⭐ **When this bundle was built, compiled into it.** David, 2026-10-03: *"the numbers otherwise
// change and I can't find the current version without asking each time."* The compose stack serves a
// BUILT IMAGE, so "the app is up" and "the app is current" are different facts and only the page can
// tell them apart — an image from last week serves perfectly and looks identical. The footer prints
// this beside the catalogue's own load time.
// ⚠ Overridable, so a reproducible build can pass a fixed value instead of a wall clock.
const BUILD_STAMP = process.env.VITE_BUILD_STAMP ?? new Date().toISOString();

export default defineConfig(({ mode }) => ({
  base: process.env.VITE_PUBLIC_BASE ?? "/",
  define: { __BUILD_STAMP__: JSON.stringify(BUILD_STAMP) },
  plugins: [vue()],
  resolve: {
    alias: { "@": fileURLToPath(new URL("./src", import.meta.url)) },
  },
  // Dev only. In every deployed configuration the API is reached at `VITE_API_BASE_URL`, and this
  // proxy does not exist. Spread rather than `proxy: undefined`, which `exactOptionalPropertyTypes`
  // rightly refuses: an absent key and a key holding `undefined` are different things.
  ...(mode === "development"
    ? {
        server: {
          proxy: {
            // ⚠ 5001, not Flask's usual 5000: on macOS the AirPlay receiver owns 5000 and answers a
            // proxied request with a 403 that looks exactly like the API refusing it.
            "/api": { target: process.env.BACATLAS_DEV_API_TARGET ?? "http://localhost:5001", changeOrigin: true },
          },
        },
      }
    : {}),
}));
