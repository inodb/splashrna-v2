/**
 * Annotation parsing and the transcript-to-genomic mapping behind the diagram.
 *
 * The mapping is the part that can be silently wrong: a guide tick drawn at the wrong
 * coordinate looks perfectly plausible. So it is checked against the genomic positions
 * the 2017 server actually printed, on both strands.
 */

import { existsSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, it } from "vitest";
import {
  StaticAnnotationIndex,
  parseIsoforms,
  parseRanges,
  transcriptCoordinates,
} from "../src/data/annotations.js";
import positions from "./fixtures/legacy-positions.json" with { type: "json" };

const PUBLIC = fileURLToPath(new URL("../public", import.meta.url));
const BUILT = existsSync(`${PUBLIC}/annotations/manifest.json`);

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

describe("parsing the legacy annotation encoding", () => {
  const text =
    "chr17,-:NM_000546,cds,7676521,7676594;cds,7676382,7676403;utr5,7687377,7687550" +
    ":NM_001126112,cds,7676521,7676594;utr3,7668402,7669608";

  it("reads chromosome, strand and one entry per isoform", () => {
    const { chrom, strand, isoforms } = parseIsoforms(text);
    expect(chrom).toBe("chr17");
    expect(strand).toBe("-");
    expect(isoforms.map((i) => i.refseq)).toEqual(["NM_000546", "NM_001126112"]);
  });

  it("keeps the accession off the exon list and classifies each block", () => {
    const { isoforms } = parseIsoforms(text);
    expect(isoforms[0].exons).toHaveLength(3);
    expect(isoforms[0].exons.map((e) => e.kind).sort()).toEqual(["cds", "cds", "utr5"]);
    // Exons come back in ascending genomic order regardless of input order.
    const starts = isoforms[0].exons.map((e) => e.start);
    expect([...starts].sort((a, b) => a - b)).toEqual(starts);
  });

  it("survives malformed input rather than throwing", () => {
    expect(parseIsoforms("").isoforms).toEqual([]);
    expect(parseIsoforms("chr1,+").isoforms).toEqual([]);
    expect(parseIsoforms("chr1,+:NM_1,cds").isoforms).toEqual([]);
  });

  it("parses intersection ranges, including the asterisk marker", () => {
    expect(parseRanges("100,200*;300,400")).toEqual([
      { start: 100, end: 200 },
      { start: 300, end: 400 },
    ]);
    expect(parseRanges("")).toEqual([]);
    expect(parseRanges("100,.")).toEqual([]);
  });

  it("strips DOUBLE asterisks, which occur in the real table", () => {
    // 12048303,12048545** is genuine data. Stripping only the first asterisk leaves a
    // non-numeric end, silently dropping the block and shifting every later
    // coordinate — which put guide ticks in the wrong place until this was caught.
    expect(parseRanges("12048303,12048545**")).toEqual([
      { start: 12048303, end: 12048545 },
    ]);
    expect(parseRanges("1,2**;3,4*;5,6")).toHaveLength(3);
  });

  it("drops ranges it cannot make sense of instead of emitting NaN", () => {
    expect(parseRanges("abc,def")).toEqual([]);
    expect(parseRanges("500,100")).toEqual([]);
  });
});

describe("transcript coordinates", () => {
  const ranges = [
    { start: 10, end: 12 },
    { start: 20, end: 21 },
  ];

  it("runs along the genome on the + strand", () => {
    expect(transcriptCoordinates(ranges, "+")).toEqual([10, 11, 12, 20, 21]);
  });

  it("runs the other way on the - strand", () => {
    expect(transcriptCoordinates(ranges, "-")).toEqual([21, 20, 12, 11, 10]);
  });
});

describe.skipIf(!BUILT)("against the 2017 server's reported positions", () => {
  it("places every guide at the coordinate the server printed", async () => {
    restore = installFileFetch();
    const index = await StaticAnnotationIndex.load("./annotations");

    const byGene = new Map<number, number[]>();
    let checked = 0;

    for (const c of positions.cases) {
      let coordinates = byGene.get(c.entrezId);
      if (!coordinates) {
        const annotation = await index.get(c.entrezId);
        expect(annotation, `gene ${c.entrezId} must be in the index`).not.toBeNull();
        coordinates = transcriptCoordinates(annotation!.intersection, annotation!.strand);
        byGene.set(c.entrezId, coordinates);
      }
      expect(
        coordinates[c.offset],
        `gene ${c.entrezId}, guide ${c.guide} at offset ${c.offset}`,
      ).toBe(c.genomic);
      checked++;
    }

    expect(checked).toBe(positions.cases.length);
  });

  it("covers both strands", async () => {
    restore = installFileFetch();
    const index = await StaticAnnotationIndex.load("./annotations");
    const strands = new Set<string>();
    for (const id of new Set(positions.cases.map((c) => c.entrezId))) {
      const a = await index.get(id);
      if (a) strands.add(a.strand);
    }
    expect([...strands].sort()).toEqual(["+", "-"]);
  });

  it("returns null for genes that are not annotated", async () => {
    restore = installFileFetch();
    const index = await StaticAnnotationIndex.load("./annotations");
    expect(await index.get(999999991)).toBeNull();
  });
});
