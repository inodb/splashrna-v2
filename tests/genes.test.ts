/**
 * Tests for the static gene index and for `*`-aware tiling.
 *
 * As with the off-target index, the 58 MB blob is derived data and is not committed,
 * so the integration tests skip unless `public/genes/` has been built. The committed
 * fixture carries a SHA-256 of each sampled gene's sequence, so a decode error cannot
 * slip through.
 */

import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, it } from "vitest";
import { StaticGeneIndex } from "../src/data/genes.js";
import { tileTranscript } from "../src/engine/sequence.js";
import sample from "./fixtures/genes-sample.json" with { type: "json" };

const GENES_DIR = fileURLToPath(new URL("../public/genes", import.meta.url));
const BUILT = existsSync(`${GENES_DIR}/manifest.json`);

function installFileFetch(root: string, onRequest?: (name: string, range?: string) => void) {
  const original = globalThis.fetch;
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const name = String(input).replace(/^.*\/genes\//, "");
    const path = `${root}/${name}`;
    const range = (init?.headers as Record<string, string> | undefined)?.Range;
    onRequest?.(name, range);
    if (!existsSync(path)) return new Response(null, { status: 404 });

    const body = readFileSync(path);
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

describe("tiling across intersection boundaries", () => {
  it("numbers positions ignoring the boundary markers", () => {
    // Two 22-nt blocks joined by a single '*'. Positions must run 1,2 then 23,
    // i.e. the '*' occupies no coordinate.
    const left = "ACGTACGTACGTACGTACGTAC";
    const right = "TTTTGGGGCCCCAAAATTTTGG";
    const tiles = tileTranscript(`${left}*${right}`);
    expect(tiles[0].position).toBe(1);
    expect(tiles.at(-1)!.position).toBe(right.length + 1);
  });

  it("emits no window that spans a boundary", () => {
    const seq = `${"A".repeat(30)}*${"C".repeat(30)}`;
    const tiles = tileTranscript(seq);
    for (const tile of tiles) {
      expect(tile.guide).not.toContain("*");
      // A window spanning the join would mix both homopolymers.
      expect(/^T+$|^G+$/.test(tile.guide)).toBe(true);
    }
  });

  it("is unchanged for ordinary sequences", () => {
    const seq = "ACGT".repeat(10);
    const tiles = tileTranscript(seq);
    expect(tiles).toHaveLength(seq.length - 22 + 1);
    expect(tiles.map((t) => t.position)).toEqual(
      tiles.map((_, i) => i + 1),
    );
  });
});

describe.skipIf(!BUILT)("static gene index", () => {
  it("decodes sampled genes byte-for-byte", async () => {
    restore = installFileFetch(GENES_DIR);
    const index = await StaticGeneIndex.load("./genes");
    expect(index.count).toBe(39760);

    for (const expected of sample.genes) {
      const gene = await index.get(expected.entrezId);
      expect(gene, `gene ${expected.entrezId} missing`).not.toBeNull();
      expect(gene!.strand).toBe(expected.strand);
      expect(gene!.sequence).toHaveLength(expected.length);
      expect(gene!.sequence.slice(0, 40)).toBe(expected.head);
      expect(gene!.sequence.slice(-40)).toBe(expected.tail);
      expect(createHash("sha256").update(gene!.sequence).digest("hex")).toBe(
        expected.sha256,
      );
    }
  });

  it("preserves boundary markers and the single N in the table", async () => {
    restore = installFileFetch(GENES_DIR);
    const index = await StaticGeneIndex.load("./genes");

    for (const expected of sample.genes) {
      const gene = (await index.get(expected.entrezId))!;
      expect([...gene.sequence].filter((c) => c === "*")).toHaveLength(expected.stars);
      expect([...gene.sequence].filter((c) => c === "N")).toHaveLength(expected.ns);
    }
    // Gene 73047 holds the only N in the whole table; it must still be retrievable.
    const withN = (await index.get(73047))!;
    expect(withN.sequence).toContain("N");
  });

  it("returns null for genes that are not in the table", async () => {
    restore = installFileFetch(GENES_DIR);
    const index = await StaticGeneIndex.load("./genes");
    for (const id of sample.absent) {
      expect(await index.get(id)).toBeNull();
    }
  });

  it("reads a gene with one range request, not the whole blob", async () => {
    const requests: { name: string; range?: string }[] = [];
    restore = installFileFetch(GENES_DIR, (name, range) => requests.push({ name, range }));
    const index = await StaticGeneIndex.load("./genes");
    requests.length = 0;

    await index.get(71458);
    const blobFetches = requests.filter((r) => r.name === "genes.bin");
    expect(blobFetches).toHaveLength(1);
    expect(blobFetches[0].range).toMatch(/^bytes=\d+-\d+$/);
  });

  it("tiles a real intersection transcript into scorable guides", async () => {
    restore = installFileFetch(GENES_DIR);
    const index = await StaticGeneIndex.load("./genes");
    const gene = (await index.get(71458))!;
    const tiles = tileTranscript(gene.sequence);

    expect(tiles.length).toBeGreaterThan(100);
    for (const tile of tiles) {
      expect(tile.guide).toHaveLength(22);
      expect(tile.guide).not.toContain("*");
      expect(tile.position).toBeGreaterThan(0);
    }
    // Positions must be strictly increasing and within the '*'-free length.
    const stripped = gene.sequence.replace(/\*/g, "").length;
    expect(tiles.at(-1)!.position).toBeLessThanOrEqual(stripped);
    for (let i = 1; i < tiles.length; i++) {
      expect(tiles[i].position).toBeGreaterThan(tiles[i - 1].position);
    }
  });
});
