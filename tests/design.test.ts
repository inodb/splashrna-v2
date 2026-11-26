/** Tests for the design pipeline, sequence helpers, FASTA parsing and CSV export. */

import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { beforeAll, describe, expect, it } from "vitest";
import { designGuides, type OffTargetIndex } from "../src/engine/design.js";
import { toCsv } from "../src/engine/export.js";
import { classify, parseFasta } from "../src/engine/fasta.js";
import { decodeModel, type RawModelBundle } from "../src/engine/model.js";
import { SplashScorer } from "../src/engine/scorer.js";
import {
  mir30Construct,
  restrictionSites,
  reverseComplement,
  tileTranscript,
} from "../src/engine/sequence.js";
import tiles from "./fixtures/tiles.json" with { type: "json" };

describe("sequence helpers", () => {
  it("reverse complements", () => {
    expect(reverseComplement("ATGC")).toBe("GCAT");
    expect(reverseComplement("AAAA")).toBe("TTTT");
  });

  it("tiles into 22-mers with 1-based positions", () => {
    const seq = "A".repeat(30);
    const result = tileTranscript(seq);
    expect(result).toHaveLength(30 - 22 + 1);
    expect(result[0].position).toBe(1);
    expect(result[0].guide).toHaveLength(22);
  });

  it("builds a 97-mer miR-30 construct", () => {
    const guide = "TTGGATCTTTCTTCAGACTGTG";
    const construct = mir30Construct(guide);
    expect(construct).toHaveLength(97);
    expect(construct.startsWith("TGCTGTTGACAGTGAGCG")).toBe(true);
    expect(construct.endsWith("TGCCTACTGCCTCGGA")).toBe(true);
    expect(construct).toContain(guide);
  });

  it("introduces the passenger mismatch at position 1", () => {
    // The passenger strand's first base is replaced: C/G becomes A, otherwise C.
    const guide = "AAAAAAAAAAAAAAAAAAAAAA"; // passenger starts with T -> replaced by C
    expect(mir30Construct(guide)).toContain("TGCTGTTGACAGTGAGCGC");
  });

  it("detects restriction sites in the construct", () => {
    expect(restrictionSites("TTGGATCTTTCTTCAGACTGTG")).toEqual([]);
    const withEcoRI = reverseComplement("GAATTC") + "A".repeat(16);
    expect(restrictionSites(withEcoRI).length + restrictionSites("GAATTC" + "A".repeat(16)).length)
      .toBeGreaterThan(0);
  });
});

describe("FASTA parsing", () => {
  it("classifies entry types", () => {
    expect(classify("NM_001127233.1")).toBe("refseq");
    expect(classify("71458")).toBe("entrez");
    expect(classify("ACGTACGT")).toBe("sequence");
    expect(classify("not a sequence!")).toBe("unknown");
  });

  it("reads multi-record input and skips comments", () => {
    const entries = parseFasta(
      "# a comment\n> one\n71458\n> two\nACGTACGTAA\nCCGG\n",
    );
    expect(entries).toHaveLength(2);
    expect(entries[0]).toMatchObject({ label: "one", value: "71458", kind: "entrez" });
    expect(entries[1].value).toBe("ACGTACGTAACCGG");
  });

  it("accepts a bare sequence with no header", () => {
    const entries = parseFasta("ACGTACGTACGT");
    expect(entries).toHaveLength(1);
    expect(entries[0].kind).toBe("sequence");
  });
});

describe("design pipeline", () => {
  let scorer: SplashScorer;
  // A real transcript, reconstructed from the fixture's first tile positions.
  const sequence = reverseComplement(tiles.tiles.slice(0, 400).map((t) => t[0]).join(""));

  beforeAll(() => {
    const path = fileURLToPath(new URL("../public/splashrna-model.json", import.meta.url));
    scorer = new SplashScorer(
      decodeModel(JSON.parse(readFileSync(path, "utf8")) as RawModelBundle),
    );
  });

  it("returns the requested number of guides, best first", async () => {
    const result = await designGuides({ label: "test", sequence }, scorer, {
      predictionsPerGene: 5,
    });
    expect(result.guides).toHaveLength(5);
    for (let i = 1; i < result.guides.length; i++) {
      expect(result.guides[i - 1].splash).toBeGreaterThanOrEqual(result.guides[i].splash);
    }
  });

  it("refuses sequences shorter than one tile", async () => {
    const result = await designGuides({ label: "short", sequence: "ACGT" }, scorer);
    expect(result.guides).toHaveLength(0);
    expect(result.warnings[0]).toMatch(/at least 22 nt/i);
  });

  it("drops guides with restriction sites when asked", async () => {
    const withSites = await designGuides({ label: "x", sequence }, scorer, {
      removeRestrictionSites: true,
      predictionsPerGene: 50,
    });
    for (const guide of withSites.guides) {
      expect(restrictionSites(guide.guide)).toEqual([]);
    }
  });

  it("applies an off-target index and reports hits", async () => {
    const first = (
      await designGuides({ label: "x", sequence }, scorer, { predictionsPerGene: 1 })
    ).guides[0].guide;

    const index: OffTargetIndex = {
      async lookup(guides) {
        const map = new Map();
        for (const g of guides) {
          map.set(g, g === first ? { human: ["SOMEGENE"], mouse: [] } : { human: [], mouse: [] });
        }
        return map;
      },
    };

    const filtered = await designGuides(
      { label: "x", sequence },
      scorer,
      { predictionsPerGene: 3, allowOffTargetDuplicates: false },
      index,
    );
    expect(filtered.guides.map((g) => g.guide)).not.toContain(first);

    const kept = await designGuides(
      { label: "x", sequence },
      scorer,
      { predictionsPerGene: 3, allowOffTargetDuplicates: true },
      index,
    );
    expect(kept.guides[0].guide).toBe(first);
    expect(kept.guides[0].offTargetHuman).toEqual(["SOMEGENE"]);
  });

  it("exports CSV with a header, citation and one row per guide", async () => {
    const result = await designGuides({ label: "gene1", sequence }, scorer, {
      predictionsPerGene: 3,
    });
    const csv = toCsv([result]);
    const lines = csv.split("\n");
    expect(csv).toContain("Pelossof");
    expect(lines.some((l) => l.startsWith("Feature,"))).toBe(true);
    expect(lines.filter((l) => l.startsWith("gene1,"))).toHaveLength(3);
  });
});
