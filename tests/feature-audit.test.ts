/**
 * Feature audit: every option must demonstrably change the output.
 *
 * The other test files check correctness — that scores and rankings match the 2017
 * server. This one checks *presence*: that each switch the UI exposes is actually
 * wired through and does something. A silently inert option is worse than a missing
 * one, because the user believes it applied.
 *
 * It deliberately asserts on observable behaviour (which guides come back) rather than
 * on internals.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { beforeAll, describe, expect, it } from "vitest";
import { StaticGeneIndex } from "../src/data/genes.js";
import { StaticOffTargetIndex } from "../src/data/offtarget.js";
import { designGuides, type OffTargetIndex } from "../src/engine/design.js";
import { toCsv } from "../src/engine/export.js";
import { classify, parseFasta } from "../src/engine/fasta.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import { reverseComplement, restrictionSites } from "../src/engine/sequence.js";
import tiles from "./fixtures/tiles.json" with { type: "json" };

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const INDEXES_BUILT =
  existsSync(`${PUBLIC}/genes/manifest.json`) && existsSync(`${PUBLIC}/offtarget/manifest.json`);

function realTranscript(count: number): string {
  const guides = tiles.tiles.slice(0, count);
  let sequence = reverseComplement(guides[0]);
  for (let i = 1; i < guides.length; i++) sequence += reverseComplement(guides[i]).slice(-1);
  return sequence;
}

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

let scorer: SplashScorer;
const sequence = realTranscript(700);

beforeAll(() => {
  scorer = new SplashScorer(
    decodeModel(JSON.parse(readFileSync(`${PUBLIC}/splashrna-model.json`, "utf8")) as RawModelBundle),
  );
});

describe("input types", () => {
  it("accepts all three forms the original accepted", () => {
    expect(classify("ACGTACGTACGTACGTACGTAC")).toBe("sequence");
    expect(classify("NM_001127233.1")).toBe("refseq");
    expect(classify("71458")).toBe("entrez");
  });

  it("parses multi-record FASTA like the original", () => {
    const entries = parseFasta("> a\n71458\n> b\nNM_001127233.1\n> c\nACGTACGTACGTACGTACGTAC");
    expect(entries.map((e) => e.kind)).toEqual(["entrez", "refseq", "sequence"]);
  });
});

describe("option: predictions per gene", () => {
  it("returns exactly the requested count", async () => {
    for (const n of [1, 5, 25]) {
      const r = await designGuides({ label: "x", sequence }, scorer, { predictionsPerGene: n });
      expect(r.guides, `requested ${n}`).toHaveLength(n);
    }
  });
});

describe("option: remove restriction sites", () => {
  it("on, no surviving guide carries an EcoRI or XhoI site", async () => {
    const r = await designGuides({ label: "x", sequence }, scorer, {
      removeRestrictionSites: true,
      predictionsPerGene: 100,
    });
    for (const g of r.guides) expect(restrictionSites(g.guide)).toEqual([]);
    expect(r.discarded.restrictionSites).toBeGreaterThan(0);
  });

  it("off, guides with sites are kept and flagged rather than dropped", async () => {
    const r = await designGuides({ label: "x", sequence }, scorer, {
      removeRestrictionSites: false,
      predictionsPerGene: 400,
    });
    expect(r.discarded.restrictionSites).toBe(0);
    const flagged = r.guides.filter((g) => restrictionSites(g.guide).length > 0);
    expect(flagged.length).toBeGreaterThan(0);
    for (const g of flagged) {
      expect(g.warnings.join(" ")).toMatch(/EcoRI|XhoI/);
    }
  });
});

describe("option: in-gene duplicates", () => {
  // A transcript with a deliberately repeated block, so some 22-mers occur twice.
  const repeated = sequence.slice(0, 300) + sequence.slice(0, 300);

  it("off, repeated guides are removed", async () => {
    const r = await designGuides({ label: "x", sequence: repeated }, scorer, {
      allowInGeneDuplicates: false,
      predictionsPerGene: 50,
    });
    expect(r.discarded.inGeneDuplicates).toBeGreaterThan(0);
    expect(new Set(r.guides.map((g) => g.guide)).size).toBe(r.guides.length);
  });

  it("on, repeated guides are kept and warned about", async () => {
    const r = await designGuides({ label: "x", sequence: repeated }, scorer, {
      allowInGeneDuplicates: true,
      predictionsPerGene: 50,
    });
    expect(r.discarded.inGeneDuplicates).toBe(0);
    expect(r.guides.some((g) => /more than once/i.test(g.warnings.join(" ")))).toBe(true);
  });
});

describe("option: off-target duplicates", () => {
  const stub = (hit: string): OffTargetIndex => ({
    async lookup(guides) {
      const map = new Map();
      for (const g of guides) {
        map.set(g, g === hit ? { humanGenes: 2, mouseGenes: 1 } : { humanGenes: 0, mouseGenes: 0 });
      }
      return map;
    },
  });

  it("off, matching guides are removed; on, they are kept and flagged", async () => {
    const top = (await designGuides({ label: "x", sequence }, scorer, { predictionsPerGene: 1 }))
      .guides[0].guide;

    const removed = await designGuides({ label: "x", sequence }, scorer,
      { predictionsPerGene: 5, allowOffTargetDuplicates: false }, stub(top));
    expect(removed.guides.map((g) => g.guide)).not.toContain(top);
    expect(removed.discarded.offTarget).toBeGreaterThan(0);

    const kept = await designGuides({ label: "x", sequence }, scorer,
      { predictionsPerGene: 5, allowOffTargetDuplicates: true }, stub(top));
    expect(kept.guides[0].guide).toBe(top);
    expect(kept.guides[0].warnings.join(" ")).toMatch(/2 human genes and 1 mouse gene/);
  });
});

describe("option: off-target mode", () => {
  it("exhaustive checks every tile; top-candidates checks only the best", async () => {
    const seen: number[] = [];
    const counting: OffTargetIndex = {
      async lookup(guides) {
        seen.push(guides.length);
        return new Map(guides.map((g) => [g, { humanGenes: 0, mouseGenes: 0 }]));
      },
    };

    seen.length = 0;
    await designGuides({ label: "x", sequence }, scorer,
      { predictionsPerGene: 5, offTargetMode: "top-candidates", offTargetCandidates: 50 }, counting);
    const top = seen[0];

    seen.length = 0;
    await designGuides({ label: "x", sequence }, scorer,
      { predictionsPerGene: 5, offTargetMode: "exhaustive" }, counting);
    const all = seen[0];

    expect(top).toBe(50);
    expect(all).toBeGreaterThan(top * 5);
  });
});

describe("output: CSV export", () => {
  it("carries every column a user needs, with one row per guide", async () => {
    const r = await designGuides({ label: "gene1", sequence }, scorer, { predictionsPerGene: 3 });
    const csv = toCsv([r]);
    const header = csv.split("\n").find((l) => l.startsWith("Feature"))!;
    for (const col of [
      "shRNA.name", "Position", "Antisense.Guide.Sequence", "SplashRNA",
      "miR30.score", "miRE.score", "97mer.construct", "Warnings",
      "Human.22mer.match.gene.count", "Mouse.22mer.match.gene.count",
      "Human.match.entrezIDs", "Mouse.match.entrezIDs",
    ]) {
      expect(header, `missing column ${col}`).toContain(col);
    }
    expect(csv).toContain("Pelossof");
    expect(csv.split("\n").filter((l) => l.startsWith("gene1,"))).toHaveLength(3);
  });
});

describe.skipIf(!INDEXES_BUILT)("data-backed features", () => {
  it("Entrez IDs resolve to the isoform intersection", async () => {
    const restore = installFileFetch();
    try {
      const genes = await StaticGeneIndex.load("./genes");
      const gene = await genes.get(7157);
      expect(gene).not.toBeNull();
      expect(gene!.sequence.length).toBeGreaterThan(1000);
      expect(gene!.strand).toMatch(/^[+-]$/);
    } finally {
      restore();
    }
  });

  it("off-target filtering consults both genomes", async () => {
    const restore = installFileFetch();
    try {
      const index = await StaticOffTargetIndex.load("./offtarget");
      const hits = await index.lookup([reverseComplement("A".repeat(22))]);
      const hit = [...hits.values()][0];
      expect(hit).toHaveProperty("humanGenes");
      expect(hit).toHaveProperty("mouseGenes");
    } finally {
      restore();
    }
  });
});
