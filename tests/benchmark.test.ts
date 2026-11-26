/**
 * Benchmarks behind docs/BENCHMARK.md.
 *
 * Not a test — it asserts almost nothing. It exists so the numbers in that document
 * can be reproduced rather than taken on trust:
 *
 *     SPLASH_BENCH=1 npx vitest run tests/benchmark.test.ts
 *
 * Skipped by default so it never slows CI.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { beforeAll, describe, it } from "vitest";
import { StaticGeneIndex, selectTranscript } from "../src/data/genes.js";
import { StaticOffTargetIndex } from "../src/data/offtarget.js";
import { designGuides } from "../src/engine/design.js";
import { CollapsedWd } from "../src/engine/collapsed.js";
import { SEQUENCE_LENGTH, wdKernel } from "../src/engine/kernels.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import { encode, tileTranscript } from "../src/engine/sequence.js";

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const ENABLED = process.env.SPLASH_BENCH === "1";
const INDEXES = existsSync(`${PUBLIC}/genes/manifest.json`);

function installFileFetch() {
  const original = globalThis.fetch;
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const name = String(input).replace(/^\.\//, "");
    const path = `${PUBLIC}/${name}`;
    if (!existsSync(path)) return new Response(null, { status: 404 });
    const body = readFileSync(path);
    const range = (init?.headers as Record<string, string> | undefined)?.Range;
    const m = range ? /bytes=(\d+)-(\d+)/.exec(range) : null;
    if (m) return new Response(body.subarray(Number(m[1]), Number(m[2]) + 1), { status: 206 });
    if (name.endsWith(".json")) return new Response(body.toString("utf8"), { status: 200 });
    return new Response(body, { status: 200 });
  }) as typeof fetch;
  return () => {
    globalThis.fetch = original;
  };
}

/** Median of `runs` timings, in milliseconds. */
async function timed(runs: number, fn: () => unknown | Promise<unknown>): Promise<number> {
  const samples: number[] = [];
  for (let i = 0; i < runs; i++) {
    const t0 = performance.now();
    await fn();
    samples.push(performance.now() - t0);
  }
  samples.sort((a, b) => a - b);
  return samples[Math.floor(samples.length / 2)];
}

const row = (label: string, value: string) => console.log(`  ${label.padEnd(38)} ${value}`);

describe.skipIf(!ENABLED)("benchmarks", () => {
  let scorer: SplashScorer;
  let bundle: RawModelBundle;

  beforeAll(() => {
    installFileFetch();
    bundle = JSON.parse(readFileSync(`${PUBLIC}/splashrna-model.json`, "utf8")) as RawModelBundle;
    scorer = new SplashScorer(decodeModel(bundle));
  });

  it("startup", async () => {
    console.log("\nSTARTUP");
    row("parse + decode model bundle", `${(await timed(5, () => decodeModel(bundle))).toFixed(0)} ms`);
    row("model bundle on disk", `${(readFileSync(`${PUBLIC}/splashrna-model.json`).length / 1024).toFixed(0)} KB`);
  });

  it("scoring throughput", async () => {
    console.log("\nSCORING");
    let seed = 1;
    const guide = Array.from({ length: SEQUENCE_LENGTH }, () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return "ACGT"[(seed >>> 16) & 3];
    }).join("");

    const one = await timed(200, () => scorer.score(guide));
    row("one guide, both classifiers", `${(one * 1000).toFixed(0)} µs`);
    row("throughput", `${(1000 / one).toFixed(0)} guides/s`);

    // Distinct pseudo-random guides. Rotating one sequence would be unrealistically
    // cache- and branch-friendly and would overstate throughput several fold.
    const guides = Array.from({ length: 2000 }, () =>
      Array.from({ length: SEQUENCE_LENGTH }, () => {
        seed = (seed * 1103515245 + 12345) & 0x7fffffff;
        return "ACGT"[(seed >>> 16) & 3];
      }).join(""),
    );
    const many = await timed(3, () => scorer.scoreAll(guides));
    row("2,000 distinct guides (a 2 kb transcript)", `${(many / 1000).toFixed(2)} s`);
    row("  per guide", `${((many / 2000) * 1000).toFixed(0)} µs`);
  });

  it("the collapsed-table optimisation", async () => {
    console.log("\nCOLLAPSED TABLES vs NAIVE LOOP");
    const model = decodeModel(bundle).mir30;
    const query = encode("TAGTGTTGCAGAATCTACGCAA");

    const fast = await timed(50, () => model.wd.apply(query));
    const naive = await timed(3, () => {
      let total = 0;
      for (let i = 0; i < model.count; i++) {
        total += model.coefWd[i] * wdKernel(model.supportVectors, i * SEQUENCE_LENGTH, query, 0);
      }
      return total;
    });
    row("naive loop over 11,585 vectors", `${(naive * 1000).toFixed(0)} µs`);
    row("collapsed lookup tables", `${(fast * 1000).toFixed(0)} µs`);
    row("speedup", `${(naive / fast).toFixed(0)}x`);

    const build = await timed(3, () => new CollapsedWd(model.supportVectors, model.coefWd, model.count));
    row("one-off table construction", `${build.toFixed(0)} ms`);
  });

  it.skipIf(!INDEXES)("end to end, by gene ID", async () => {
    console.log("\nEND TO END (gene ID -> ranked shRNAs)");
    const genes = await StaticGeneIndex.load("./genes");
    const offTarget = await StaticOffTargetIndex.load("./offtarget");

    const load = await timed(3, () => StaticGeneIndex.load("./genes"));
    row("gene index load (1.27 MB)", `${load.toFixed(0)} ms`);

    for (const id of [7157, 672]) {
      const gene = (await genes.get(id))!;
      const choice = selectTranscript(gene, { polyACut: true });
      const tiles = tileTranscript(choice.sequence).length;

      const noFilter = await timed(3, () =>
        designGuides({ label: "x", sequence: choice.sequence }, scorer, { predictionsPerGene: 10 }),
      );
      const filtered = await timed(3, () =>
        designGuides(
          { label: "x", sequence: choice.sequence },
          scorer,
          { predictionsPerGene: 10, offTargetMode: "exhaustive" },
          offTarget,
        ),
      );
      row(`gene ${id} (${tiles.toLocaleString()} tiles)`, `${(noFilter / 1000).toFixed(2)} s`);
      row(`  + exhaustive off-target filtering`, `${(filtered / 1000).toFixed(2)} s`);
    }
  });
});
