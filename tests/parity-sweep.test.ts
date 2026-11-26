/**
 * Comprehensive parity sweep against the 2017 server.
 *
 * The other parity tests cover a handful of hand-picked genes. This one covers a
 * stratified sample drawn across the whole feature space — genes whose polyA site
 * resolves and genes that fall back to the CDS or to the untrimmed transcript, short
 * and long transcripts, both strands — crossed with every transcript-selection and
 * filtering option.
 *
 * Every case is a response the live legacy server actually produced, captured while it
 * was still running. It cannot be regenerated once that host is gone, which is why the
 * fixture is committed.
 *
 * The full sweep takes a few minutes, so by default a deterministic subset runs.
 * Set SPLASH_FULL_PARITY=1 to run all of it.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterAll, beforeAll, describe, expect, it } from "vitest";
import { StaticGeneIndex, selectTranscript } from "../src/data/genes.js";
import { StaticOffTargetIndex } from "../src/data/offtarget.js";
import { designGuides, type DesignOptions } from "../src/engine/design.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const FIXTURE = fileURLToPath(new URL("./fixtures/parity-sweep.json", import.meta.url));
const READY =
  existsSync(FIXTURE) &&
  existsSync(`${PUBLIC}/genes/manifest.json`) &&
  existsSync(`${PUBLIC}/offtarget/manifest.json`);

const FULL = process.env.SPLASH_FULL_PARITY === "1";
const SUBSET = 120;

interface Case {
  tag: string;
  entrezId: number;
  options: { n: number; apa: string; cds: string; re_: string; ing: string; outg: string };
  /** Guide sequences in table order, duplicates preserved. */
  guides: string[];
  /** The combined SplashRNA score the 2017 server printed for each. */
  scores: number[];
}

/**
 * Shogun passed kernel values through a float32 cache, so the original scores carry
 * a few float32 ulps of noise. Anything inside this is agreement.
 */
const SCORE_TOLERANCE = 1e-6;

/** Serve public/ over fetch, honouring Range exactly as a static host would. */
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

describe.skipIf(!READY)("comprehensive parity sweep", () => {
  let restore: () => void;
  let scorer: SplashScorer;
  let genes: StaticGeneIndex;
  let offTarget: StaticOffTargetIndex;
  let cases: Case[];

  beforeAll(async () => {
    restore = installFileFetch();
    scorer = new SplashScorer(
      decodeModel(JSON.parse(readFileSync(`${PUBLIC}/splashrna-model.json`, "utf8")) as RawModelBundle),
    );
    genes = await StaticGeneIndex.load("./genes");
    offTarget = await StaticOffTargetIndex.load("./offtarget");

    const all = JSON.parse(readFileSync(FIXTURE, "utf8")) as Case[];
    // Deterministic stride, so the subset still spans every variant and stratum.
    cases = FULL ? all : all.filter((_, i) => i % Math.ceil(all.length / SUBSET) === 0);
  });

  afterAll(() => restore?.());

  it(`reproduces the 2017 server on every captured case`, async () => {
    const mismatches: string[] = [];
    const byTag = new Map<string, { total: number; ok: number }>();
    let worst = 0;

    for (const c of cases) {
      const stat = byTag.get(c.tag) ?? { total: 0, ok: 0 };
      stat.total++;

      const gene = await genes.get(c.entrezId);
      if (!gene) {
        mismatches.push(`${c.tag}/${c.entrezId}: gene missing from the index`);
        byTag.set(c.tag, stat);
        continue;
      }

      const choice = selectTranscript(gene, {
        polyACut: c.options.apa === "on",
        cdsOnly: c.options.cds === "on",
      });

      const options: Partial<DesignOptions> = {
        predictionsPerGene: c.options.n,
        removeRestrictionSites: c.options.re_ === "on",
        allowInGeneDuplicates: c.options.ing === "on",
        allowOffTargetDuplicates: c.options.outg === "on",
        // The original filtered every tile before truncating.
        offTargetMode: "exhaustive",
      };

      const designed = choice.sequence
        ? (
            await designGuides(
              { label: String(c.entrezId), sequence: choice.sequence },
              scorer,
              options,
              offTarget,
            )
          ).guides
        : [];

      const ours = designed.map((g) => g.guide);
      const guidesAgree = JSON.stringify(ours) === JSON.stringify(c.guides);
      const worstScore = guidesAgree
        ? Math.max(0, ...designed.map((g, i) => Math.abs(g.splash - c.scores[i])))
        : Infinity;

      if (guidesAgree && worstScore <= SCORE_TOLERANCE) {
        stat.ok++;
        worst = Math.max(worst, worstScore);
      } else if (mismatches.length < 10) {
        mismatches.push(
          guidesAgree
            ? `${c.tag}/${c.entrezId}: guides agree but score differs by ${worstScore.toExponential(2)}`
            : `${c.tag}/${c.entrezId}: expected ${c.guides.length} guides ` +
              `[${c.guides.slice(0, 2).join(",")}…], got ${ours.length} [${ours.slice(0, 2).join(",")}…]`,
        );
      }
      byTag.set(c.tag, stat);
    }

    const lines = [...byTag.entries()]
      .sort()
      .map(([tag, s]) => `  ${tag.padEnd(9)} ${s.ok}/${s.total}`);
    console.log(
      `\nparity by option variant (${FULL ? "full sweep" : `subset of ${cases.length}`}):\n` +
        lines.join("\n") +
        `\n  largest score difference: ${worst.toExponential(2)}`,
    );

    expect(mismatches, mismatches.slice(0, 10).join("\n")).toEqual([]);
  }, 900_000);
});
