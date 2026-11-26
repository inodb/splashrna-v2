/**
 * End-to-end parity with the original server.
 *
 * The 2017 CherryPy instance was still running while this rewrite was built, so this
 * compares the complete Entrez pipeline -- isoform intersection, tiling across block
 * boundaries, scoring, in-gene de-duplication, off-target filtering and ranking --
 * against what that server actually returned.
 *
 * The fixture was captured with apa and cdsOnly off, so the server tiled the
 * intersection directly rather than using its precomputed cache, which is the path
 * this implementation takes.
 *
 * Requires both derived indexes; skipped when they have not been built.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, it } from "vitest";
import { StaticGeneIndex } from "../src/data/genes.js";
import { StaticOffTargetIndex } from "../src/data/offtarget.js";
import { designGuides } from "../src/engine/design.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import legacy from "./fixtures/legacy-genes.json" with { type: "json" };

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const BUILT =
  existsSync(`${PUBLIC}/genes/manifest.json`) &&
  existsSync(`${PUBLIC}/offtarget/manifest.json`);

/** Serve public/ over a stubbed fetch, honouring Range. */
function installFileFetch() {
  const original = globalThis.fetch;
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const name = String(input).replace(/^\.\//, "");
    const path = `${PUBLIC}/${name}`;
    if (!existsSync(path)) return new Response(null, { status: 404 });

    const body = readFileSync(path);
    const range = (init?.headers as Record<string, string> | undefined)?.Range;
    const match = range ? /bytes=(\d+)-(\d+)/.exec(range) : null;
    if (match) {
      return new Response(body.subarray(Number(match[1]), Number(match[2]) + 1), {
        status: 206,
      });
    }
    if (name.endsWith(".json")) return new Response(body.toString("utf8"), { status: 200 });
    return new Response(body, { status: 200 });
  }) as typeof fetch;
  return () => {
    globalThis.fetch = original;
  };
}

let restore: (() => void) | null = null;
afterEach(() => {
  restore?.();
  restore = null;
});

describe.skipIf(!BUILT)("parity with the 2017 server", () => {
  for (const testCase of legacy.cases) {
    it(`returns the same ten shRNAs, in order, for Entrez gene ${testCase.entrezId}`, async () => {
      restore = installFileFetch();

      const scorer = new SplashScorer(
        decodeModel(
          JSON.parse(
            readFileSync(`${PUBLIC}/splashrna-model.json`, "utf8"),
          ) as RawModelBundle,
        ),
      );
      const genes = await StaticGeneIndex.load("./genes");
      const offTarget = await StaticOffTargetIndex.load("./offtarget");

      const gene = await genes.get(testCase.entrezId);
      expect(gene, `gene ${testCase.entrezId} must be in the index`).not.toBeNull();

      const result = await designGuides(
        { label: `Entrez ${testCase.entrezId}`, sequence: gene!.sequence },
        scorer,
        {
          predictionsPerGene: legacy.options.predictionsPerGene,
          removeRestrictionSites: legacy.options.removeRestrictionSites,
          allowInGeneDuplicates: legacy.options.allowInGeneDuplicates,
          allowOffTargetDuplicates: legacy.options.allowOffTargetDuplicates,
          // The original filtered every tile before truncating.
          offTargetMode: "exhaustive",
        },
        offTarget,
      );

      expect(result.guides.map((g) => g.guide)).toEqual(testCase.guides);
    });
  }
});
