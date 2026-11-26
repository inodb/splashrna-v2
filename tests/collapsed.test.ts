/**
 * The collapsed weighted-degree tables must be exactly equivalent to looping the
 * kernel over every support vector. This test is the guard on that optimisation:
 * if the table construction ever drifts, this fails before the slower end-to-end
 * fixtures do, and points straight at the cause.
 */

import { describe, expect, it } from "vitest";
import { CollapsedWd } from "../src/engine/collapsed.js";
import { SEQUENCE_LENGTH, wdKernel } from "../src/engine/kernels.js";
import { encode } from "../src/engine/sequence.js";

/** Deterministic pseudo-random sequences, so failures are reproducible. */
function makeSequences(n: number, seed = 42): string[] {
  let state = seed;
  const next = () => {
    state = (state * 1103515245 + 12345) & 0x7fffffff;
    return state / 0x7fffffff;
  };
  return Array.from({ length: n }, () =>
    Array.from({ length: SEQUENCE_LENGTH }, () => "ACGT"[Math.floor(next() * 4)]).join(""),
  );
}

function naive(
  svs: Uint8Array,
  coef: Float64Array,
  count: number,
  query: Uint8Array,
): number {
  let total = 0;
  for (let i = 0; i < count; i++) {
    total += coef[i] * wdKernel(svs, i * SEQUENCE_LENGTH, query, 0);
  }
  return total;
}

describe("collapsed weighted-degree tables", () => {
  const supportStrings = makeSequences(300, 7);
  const queries = makeSequences(40, 99);
  const svs = encode(supportStrings.join(""));
  const coef = Float64Array.from(supportStrings, (_, i) => Math.sin(i) * 0.03);
  const collapsed = new CollapsedWd(svs, coef, supportStrings.length);

  it("equals the per-support-vector loop", () => {
    for (const q of queries) {
      const encoded = encode(q);
      expect(collapsed.apply(encoded)).toBeCloseTo(
        naive(svs, coef, supportStrings.length, encoded),
        10,
      );
    }
  });

  it("handles long matches, where the sparse degree-8 index is exercised", () => {
    // A query identical to a support vector maximises run lengths, so every degree
    // from 1 to 10 contributes and both the dense and sparse paths are used.
    for (const seq of supportStrings.slice(0, 10)) {
      const encoded = encode(seq);
      expect(collapsed.apply(encoded)).toBeCloseTo(
        naive(svs, coef, supportStrings.length, encoded),
        10,
      );
    }
  });

  it("handles homopolymers, the worst case for run lengths", () => {
    for (const seq of ["A".repeat(22), "ACACACACACACACACACACAC", "G".repeat(22)]) {
      const encoded = encode(seq);
      expect(collapsed.apply(encoded)).toBeCloseTo(
        naive(svs, coef, supportStrings.length, encoded),
        10,
      );
    }
  });
});
