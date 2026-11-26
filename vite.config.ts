import { defineConfig } from "vitest/config";

export default defineConfig({
  // Relative base so the built site works from any path, including a subdirectory
  // or a plain file:// copy on a lab share.
  base: "./",
  build: {
    outDir: "dist",
    target: "es2022",
    assetsInlineLimit: 0,
  },
  test: {
    globals: true,
    environment: "node",
    include: ["tests/**/*.test.ts"],
  },
});
