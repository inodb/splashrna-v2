/**
 * Parity for the polyA cut and CDS-only options.
 *
 * The 2017 server trimmed each transcript at its annotated cleavage site before
 * designing, and had that ON BY DEFAULT. The earlier parity fixtures were captured
 * with it off, so they prove nothing about this path — these do.
 *
 * Captured across 8 genes at n=15 (which forces the real computation) and n=10 (which
 * goes through the server's precomputed cache). The two agreed on every gene, so both
 * are pinned: if a future change broke only the cached path's assumptions, that would
 * still show up here.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterEach, beforeAll, describe, expect, it } from "vitest";
import { StaticGeneIndex, selectTranscript } from "../src/data/genes.js";
import { StaticOffTargetIndex } from "../src/data/offtarget.js";
import { designGuides } from "../src/engine/design.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import legacy from "./fixtures/legacy-apa.json" with { type: "json" };

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const BUILT =
  existsSync(`${PUBLIC}/genes/manifest.json`) && existsSync(`${PUBLIC}/offtarget/manifest.json`);

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

let restore: (() => void) | null = null;
afterEach(() => {
  restore?.();
  restore = null;
});

let scorer: SplashScorer;
beforeAll(() => {
  scorer = new SplashScorer(
    decodeModel(JSON.parse(readFileSync(`${PUBLIC}/splashrna-model.json`, "utf8")) as RawModelBundle),
  );
});

describe.skipIf(!BUILT)("polyA cut and CDS-only parity with the 2017 server", () => {
  for (const testCase of legacy.cases) {
    const label = testCase.options.cdsOnly ? "CDS only" : "polyA cut";
    it(`${label}, n=${testCase.predictionsPerGene}, gene ${testCase.entrezId}`, async () => {
      restore = installFileFetch();
      const genes = await StaticGeneIndex.load("./genes");
      const offTarget = await StaticOffTargetIndex.load("./offtarget");

      const gene = await genes.get(testCase.entrezId);
      expect(gene, `gene ${testCase.entrezId} must be in the index`).not.toBeNull();

      const choice = selectTranscript(gene!, testCase.options);
      expect(choice.sequence.length, "selected transcript must be non-empty").toBeGreaterThan(0);

      const result = await designGuides(
        { label: `Entrez ${testCase.entrezId}`, sequence: choice.sequence },
        scorer,
        {
          predictionsPerGene: testCase.predictionsPerGene,
          removeRestrictionSites: false,
          allowInGeneDuplicates: false,
          allowOffTargetDuplicates: false,
          offTargetMode: "exhaustive",
        },
        offTarget,
      );

      expect(result.guides.map((g) => g.guide)).toEqual(testCase.guides);
    });
  }
});

describe("transcript selection rules", () => {
  const gene = {
    entrezId: 1,
    strand: "+" as const,
    sequence: "ACGTACGTACGTACGTACGTACGT",
    cds: "TTTTGGGG",
    polyACutIndex: 12,
  };

  it("trims at the cut when one is known", () => {
    const c = selectTranscript(gene, { polyACut: true });
    expect(c.sequence).toBe(gene.sequence.slice(0, 12));
    expect(c.source).toBe("polyA-trimmed");
    expect(c.warnings).toEqual([]);
  });

  it("leaves the transcript untouched when the cut is switched off", () => {
    const c = selectTranscript(gene, { polyACut: false });
    expect(c.sequence).toBe(gene.sequence);
    expect(c.source).toBe("intersection");
  });

  it("substitutes the CDS when no polyA annotation exists, and says so", () => {
    const c = selectTranscript({ ...gene, polyACutIndex: -1 }, { polyACut: true });
    expect(c.sequence).toBe(gene.cds);
    expect(c.source).toBe("cds");
    expect(c.warnings.join(" ")).toMatch(/CDS intersection was substituted/);
  });

  it("falls back to the full transcript when there is no CDS either", () => {
    const c = selectTranscript({ ...gene, polyACutIndex: -1, cds: "" }, { polyACut: true });
    expect(c.sequence).toBe(gene.sequence);
    expect(c.warnings.join(" ")).toMatch(/no intersection for CDS regions/);
  });

  it("uses the CDS directly when CDS-only is set", () => {
    const c = selectTranscript(gene, { cdsOnly: true });
    expect(c.sequence).toBe(gene.cds);
  });

  it("reports rather than guesses when CDS-only is set but no CDS exists", () => {
    const c = selectTranscript({ ...gene, cds: "" }, { cdsOnly: true });
    expect(c.sequence).toBe("");
    expect(c.warnings.join(" ")).toMatch(/no CDS intersection/);
  });
});
