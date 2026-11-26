/**
 * End-to-end scoring tests against the original server.
 *
 * `scores.json` holds per-model scores for 44 guides; `tiles.json` holds every tile of
 * a real 2 kb transcript (Entrez 100129520), scored by the Shogun models in 2017-era
 * form. Agreement on the *ranking* of that transcript is the acceptance criterion for
 * this rewrite: it is what a user actually sees.
 */

import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { beforeAll, describe, expect, it } from "vitest";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import scores from "./fixtures/scores.json" with { type: "json" };
import tiles from "./fixtures/tiles.json" with { type: "json" };

/**
 * Shogun passes kernel values through a float32 cache (`KERNELCACHE_ELEM`), so the
 * original scores carry roughly four float32 ulps of noise. We reproduce them to well
 * inside that, which means this implementation is slightly *more* precise than the one
 * it replaces.
 */
const TOLERANCE = 1e-6;

let scorer: SplashScorer;

beforeAll(() => {
  const path = fileURLToPath(new URL("../public/splashrna-model.json", import.meta.url));
  const bundle = JSON.parse(readFileSync(path, "utf8")) as RawModelBundle;
  scorer = new SplashScorer(decodeModel(bundle));
});

describe("model bundle", () => {
  it("decodes both classifiers", () => {
    expect(scorer.citation).toContain("Pelossof");
  });
});

describe("individual guide scores", () => {
  it("reproduces the miR-30 classifier", () => {
    scores.tests.forEach((guide, i) => {
      expect(scorer.score(guide).mir30).toBeCloseTo(scores.mir30_scores[i], 6);
    });
  });

  it("reproduces the miR-E classifier", () => {
    scores.tests.forEach((guide, i) => {
      expect(scorer.score(guide).mirE).toBeCloseTo(scores.mirE_scores[i], 6);
    });
  });

  it("rejects guides of the wrong length", () => {
    expect(() => scorer.score("ACGT")).toThrow(/22 nt/);
  });
});

describe("whole-transcript agreement", () => {
  let computed: { mir30: number; mirE: number; splash: number }[];

  beforeAll(() => {
    computed = scorer.scoreAll(tiles.tiles);
  });

  it("matches every tile of a real transcript", () => {
    let worst = 0;
    for (let i = 0; i < tiles.tiles.length; i++) {
      worst = Math.max(
        worst,
        Math.abs(computed[i].mir30 - tiles.mir30[i]),
        Math.abs(computed[i].mirE - tiles.mirE[i]),
      );
    }
    expect(worst).toBeLessThan(TOLERANCE);
  });

  it("produces an identical ranking", () => {
    const reference = tiles.tiles
      .map((guide, i) => ({
        guide,
        splash: scorer.combine(tiles.mir30[i], tiles.mirE[i]),
      }))
      .sort((a, b) => b.splash - a.splash);

    const ours = tiles.tiles
      .map((guide, i) => ({ guide, splash: computed[i].splash }))
      .sort((a, b) => b.splash - a.splash);

    for (const n of [5, 10, 25, 100]) {
      expect(ours.slice(0, n).map((g) => g.guide)).toEqual(
        reference.slice(0, n).map((g) => g.guide),
      );
    }
  });
});
