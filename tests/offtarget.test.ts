/**
 * Tests for the static off-target index.
 *
 * The index files are large derived data and are not committed, so the integration
 * test here is skipped unless `public/offtarget/` has been built. What *is* committed
 * is a 750-entry sample drawn from the original legacy tables, which is what makes the
 * integration test meaningful when the index is present.
 *
 * The decoder itself is always tested, against an index encoded here in the test.
 */

import { existsSync, readFileSync, statSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, it } from "vitest";
import { StaticOffTargetIndex, packKmer } from "../src/data/offtarget.js";
import { designGuides } from "../src/engine/design.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import { reverseComplement } from "../src/engine/sequence.js";
import sample from "./fixtures/offtarget-sample.json" with { type: "json" };
import geneCases from "./fixtures/offtarget-genes.json" with { type: "json" };
import tiles from "./fixtures/tiles.json" with { type: "json" };

const OFFTARGET_DIR = fileURLToPath(new URL("../public/offtarget", import.meta.url));
const BUILT = existsSync(`${OFFTARGET_DIR}/manifest.json`);

/** Serve local files over a stubbed fetch, honouring Range exactly as a host would. */
function installFileFetch(root: string, onRequest?: (url: string, range?: string) => void) {
  const original = globalThis.fetch;
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const url = String(input);
    const name = url.replace(/^.*\/offtarget\//, "");
    const path = `${root}/${name}`;
    const range = (init?.headers as Record<string, string> | undefined)?.Range;
    onRequest?.(name, range);

    if (!existsSync(path)) {
      return new Response(null, { status: 404 });
    }
    const body = readFileSync(path);
    if (range) {
      const match = /bytes=(\d+)-(\d+)/.exec(range);
      if (match) {
        const slice = body.subarray(Number(match[1]), Number(match[2]) + 1);
        return new Response(slice, { status: 206 });
      }
    }
    if (name.endsWith(".json")) {
      return new Response(body.toString("utf8"), { status: 200 });
    }
    return new Response(body, { status: 200 });
  }) as typeof fetch;
  return () => {
    globalThis.fetch = original;
  };
}

/**
 * Rebuild the real transcript the tile fixture came from.
 *
 * Each tile is the reverse complement of a 22-mer at consecutive positions, so the
 * sense strand is the first tile's reverse complement plus one base per tile after it.
 * This matters: a synthetic sequence has essentially no genuine genome matches, so
 * off-target filtering would look inert when tested against one.
 */
function realTranscript(count: number): string {
  const guides = tiles.tiles.slice(0, count);
  let sequence = reverseComplement(guides[0]);
  for (let i = 1; i < guides.length; i++) {
    sequence += reverseComplement(guides[i]).slice(-1);
  }
  return sequence;
}

let restore: (() => void) | null = null;
afterEach(() => {
  restore?.();
  restore = null;
});

describe("k-mer packing", () => {
  it("packs 22-mers inside the exact-integer range", () => {
    expect(packKmer("A".repeat(22))).toBe(0);
    const max = packKmer("T".repeat(22))!;
    expect(max).toBe(4 ** 22 - 1);
    expect(Number.isSafeInteger(max)).toBe(true);
  });

  it("is injective on distinct sequences", () => {
    const seen = new Set<number>();
    for (const s of ["ACGTACGTACGTACGTACGTAC", "ACGTACGTACGTACGTACGTAG", "A".repeat(22)]) {
      const packed = packKmer(s)!;
      expect(seen.has(packed)).toBe(false);
      seen.add(packed);
    }
  });

  it("rejects ambiguous bases", () => {
    expect(packKmer("N".repeat(22))).toBeNull();
    expect(packKmer("ACGTNCGTACGTACGTACGTAC")).toBeNull();
  });
});

describe.skipIf(!BUILT)("against the built index", () => {
  it("reports the right gene count for keys present in the legacy tables", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const index = await StaticOffTargetIndex.load("./offtarget");

    for (const genome of ["human", "mouse"] as const) {
      const entries = sample.present[genome];
      // Look up by guide, i.e. the reverse complement of the stored sense 22-mer.
      const guides = entries.map(([kmer]) => reverseComplement(kmer as string));
      const hits = await index.lookup(guides);

      const field = genome === "human" ? "humanGenes" : "mouseGenes";
      for (let i = 0; i < entries.length; i++) {
        const [kmer, count] = entries[i] as [string, number];
        const hit = hits.get(guides[i]);
        expect(hit, `missing result for ${genome} entry ${i}`).toBeDefined();
        expect(hit![field], `${genome} ${kmer} gene count`).toBe(count);
      }
    }
  });

  it("reports no hits for 22-mers absent from both tables", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const index = await StaticOffTargetIndex.load("./offtarget");
    const guides = sample.absent.map((k) => reverseComplement(k));
    const hits = await index.lookup(guides);

    for (const guide of guides) {
      const hit = hits.get(guide)!;
      expect(hit.humanGenes).toBe(0);
      expect(hit.mouseGenes).toBe(0);
    }
  });

  it("uses range requests for small queries rather than the whole file", async () => {
    const requests: { name: string; range?: string }[] = [];
    restore = installFileFetch(OFFTARGET_DIR, (name, range) => requests.push({ name, range }));

    const index = await StaticOffTargetIndex.load("./offtarget");
    requests.length = 0;
    await index.lookup([reverseComplement(sample.present.human[0][0] as string)]);

    const keyFetches = requests.filter((r) => r.name.endsWith("-keys.bin"));
    expect(keyFetches.length).toBeGreaterThan(0);
    for (const fetchCall of keyFetches) {
      expect(fetchCall.range, "expected a Range request").toMatch(/^bytes=\d+-\d+$/);
    }

    // A single lookup must not pull down a multi-megabyte file.
    const whole = statSync(`${OFFTARGET_DIR}/human-keys.bin`).size;
    expect(whole).toBeGreaterThan(1_000_000);
  });

  it("falls back to a full download when a query touches most blocks", async () => {
    const requests: { name: string; range?: string }[] = [];
    restore = installFileFetch(OFFTARGET_DIR, (name, range) => requests.push({ name, range }));

    const index = await StaticOffTargetIndex.load("./offtarget");
    requests.length = 0;

    // Spread lookups across the whole key space so nearly every block is touched.
    const guides: string[] = [];
    for (let i = 0; i < 4000; i++) {
      let seed = i * 2654435761;
      let seq = "";
      for (let j = 0; j < 22; j++) {
        seed = (seed * 1103515245 + 12345) & 0x7fffffff;
        seq += "ACGT"[(seed >>> 16) & 3];
      }
      guides.push(seq);
    }
    await index.lookup(guides);

    const full = requests.filter((r) => r.name.endsWith("-keys.bin") && !r.range);
    expect(full.length).toBeGreaterThan(0);
  });
});

function loadScorer(): SplashScorer {
  const modelPath = fileURLToPath(
    new URL("../public/splashrna-model.json", import.meta.url),
  );
  return new SplashScorer(
    decodeModel(JSON.parse(readFileSync(modelPath, "utf8")) as RawModelBundle),
  );
}

describe.skipIf(!BUILT)("off-target gene identities", () => {
  it("returns exactly the Entrez IDs the legacy tables list", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const index = await StaticOffTargetIndex.load("./offtarget");

    // Guides are the reverse complement of the sense k-mer the tables are keyed by.
    const guides = geneCases.cases.map((c) => reverseComplement(c.kmer));
    const hits = await index.lookup(guides);
    await index.annotate(hits);

    for (let i = 0; i < geneCases.cases.length; i++) {
      const c = geneCases.cases[i];
      const hit = hits.get(guides[i]);
      expect(hit, `missing ${c.genome} ${c.kmer}`).toBeDefined();
      const ids = (c.genome === "human" ? hit!.humanIds : hit!.mouseIds) ?? [];
      expect([...ids].sort((a, b) => a - b), `${c.genome} ${c.kmer}`).toEqual(c.geneIds);
      // The count and the list must agree.
      const count = c.genome === "human" ? hit!.humanGenes : hit!.mouseGenes;
      expect(count).toBe(c.geneIds.length);
    }
  });

  it("leaves IDs empty for guides with no match", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const index = await StaticOffTargetIndex.load("./offtarget");
    const guides = sample.absent.slice(0, 20).map((k) => reverseComplement(k));
    const hits = await index.lookup(guides);
    await index.annotate(hits);
    for (const g of guides) {
      expect(hits.get(g)!.humanIds ?? []).toEqual([]);
      expect(hits.get(g)!.mouseIds ?? []).toEqual([]);
    }
  });
});

describe.skipIf(!BUILT)("end to end with off-target filtering", () => {
  it("removes guides that match other genes", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const scorer = loadScorer();
    const index = await StaticOffTargetIndex.load("./offtarget");

    // Build a transcript that deliberately contains known non-unique 22-mers, spaced
    // apart by filler. The real test transcript happens to be entirely unique, so it
    // cannot show that the filter removes anything.
    const planted = sample.present.human.slice(0, 6).map(([kmer]) => kmer as string);
    const filler = "ACGATCGTTAGCCATGACTTAGCACGATTGCA";
    const sequence = planted.map((kmer) => filler + kmer).join("") + filler;
    const input = { label: "planted", sequence };

    const result = await designGuides(
      input,
      scorer,
      {
        predictionsPerGene: 500,
        removeRestrictionSites: false,
        offTargetMode: "exhaustive",
        allowOffTargetDuplicates: false,
      },
      index,
    );

    // Every planted 22-mer must have been removed. The guide is the reverse
    // complement of the sense k-mer stored in the table.
    const survivors = new Set(result.guides.map((g) => g.guide));
    for (const kmer of planted) {
      expect(survivors.has(reverseComplement(kmer)), `${kmer} should be filtered`).toBe(
        false,
      );
    }
    expect(result.discarded.offTarget).toBeGreaterThanOrEqual(planted.length);

    // And nothing unique was removed by mistake.
    for (const guide of result.guides) {
      expect(guide.offTarget!.humanGenes).toBe(0);
      expect(guide.offTarget!.mouseGenes).toBe(0);
    }
  });

  it("leaves a genuinely unique transcript untouched", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const scorer = loadScorer();
    const index = await StaticOffTargetIndex.load("./offtarget");

    // Entrez 100129520: none of its 1,979 tiles occur in a second gene, verified
    // directly against the legacy tables. The filter must not invent hits.
    const result = await designGuides(
      { label: "unique", sequence: realTranscript(600) },
      scorer,
      { predictionsPerGene: 25, offTargetMode: "exhaustive" },
      index,
    );
    expect(result.discarded.offTarget).toBe(0);
    expect(result.guides).toHaveLength(25);
  });

  it("reports hits instead of removing them when duplicates are allowed", async () => {
    restore = installFileFetch(OFFTARGET_DIR);
    const modelPath = fileURLToPath(
      new URL("../public/splashrna-model.json", import.meta.url),
    );
    const scorer = new SplashScorer(
      decodeModel(JSON.parse(readFileSync(modelPath, "utf8")) as RawModelBundle),
    );
    const index = await StaticOffTargetIndex.load("./offtarget");
    const sequence = realTranscript(600);

    const kept = await designGuides(
      { label: "e2e", sequence },
      scorer,
      { predictionsPerGene: 25, allowOffTargetDuplicates: true },
      index,
    );
    expect(kept.discarded.offTarget).toBe(0);
    for (const guide of kept.guides) expect(guide.offTarget).not.toBeNull();
  });
});
