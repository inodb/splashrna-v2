/**
 * Kernel-level tests against values captured from the original Shogun 2.1.0 build.
 *
 * These are the tests that would catch a silent change in the scoring maths, so they
 * compare against recorded ground truth rather than against our own implementation.
 */

import { describe, expect, it } from "vitest";
import {
  DEGREE,
  SPECTRUM_WEIGHTS,
  WD_WEIGHTS,
  spectrumCounts,
  spectrumDot,
  wdKernel,
} from "../src/engine/kernels.js";
import { encode } from "../src/engine/sequence.js";
import probe from "./fixtures/spectrum-probe.json" with { type: "json" };
import scores from "./fixtures/scores.json" with { type: "json" };

const wd = (a: string, b: string) => {
  const ea = encode(a);
  const eb = encode(b);
  return wdKernel(ea, 0, eb, 0, a.length);
};

describe("weighted-degree position kernel", () => {
  it("uses Shogun's degree weighting", () => {
    // w_k = 2(d - k + 1) / (d(d + 1)); the series must sum to 1.
    expect(WD_WEIGHTS).toHaveLength(DEGREE);
    expect(WD_WEIGHTS.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 12);
    expect(WD_WEIGHTS[0]).toBeCloseTo(2 / 11, 12);
  });

  it("gives the analytic value for a non-repetitive self-comparison", () => {
    // With no shift contribution, K(x, x) = sum_k w_k (L - k + 1) = 19 for L = 22.
    expect(wd("ACGTACGTACGTACGTACGTAC", "ACGTACGTACGTACGTACGTAC")).toBeCloseTo(19, 10);
  });

  it("matches Shogun on the recorded pair matrix", () => {
    const { pairs_a, pairs_b, raw_wd } = scores;
    for (let i = 0; i < pairs_a.length; i++) {
      for (let j = 0; j < pairs_b.length; j++) {
        expect(wd(pairs_a[i], pairs_b[j])).toBeCloseTo(raw_wd[i][j], 10);
      }
    }
  });

  it("is symmetric", () => {
    expect(wd("TAAAAATCGAGGCTGTTCAGTT", "TAAAAATCGAGGCTGTTCAGTA")).toBeCloseTo(
      wd("TAAAAATCGAGGCTGTTCAGTA", "TAAAAATCGAGGCTGTTCAGTT"),
      12,
    );
  });
});

describe("weighted spectrum kernel", () => {
  const spec = (a: string, b: string) => {
    const ca = spectrumCounts(encode(a), 0, a.length);
    const cb = spectrumCounts(encode(b), 0, b.length);
    return spectrumDot(ca, cb);
  };

  it("uses suffix-mask weights of 1/2, 1/3, 1/6", () => {
    expect(SPECTRUM_WEIGHTS).toEqual([1 / 2, 1 / 3, 1 / 6]);
    expect(SPECTRUM_WEIGHTS.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 12);
  });

  it("matches Shogun across the full probe matrix", () => {
    const { probes, M } = probe;
    for (let i = 0; i < probes.length; i++) {
      for (let j = 0; j < probes.length; j++) {
        expect(spec(probes[i], probes[j])).toBeCloseTo(M[i][j], 10);
      }
    }
  });

  it("masks suffixes rather than prefixes", () => {
    // AAAAAAA vs AAAAAAC is 20 under suffix masking and 24.167 under prefix masking.
    // This single value is what pins down the orientation of Shogun's reverse encoding.
    expect(spec("AAAAAAA", "AAAAAAC")).toBeCloseTo(20, 10);
  });

  it("matches the recorded split-half kernels", () => {
    const { pairs_a, pairs_b, raw_spec1, raw_spec2 } = scores;
    for (let i = 0; i < pairs_a.length; i++) {
      for (let j = 0; j < pairs_b.length; j++) {
        expect(spec(pairs_a[i].slice(0, 15), pairs_b[j].slice(0, 15))).toBeCloseTo(
          raw_spec1[i][j],
          10,
        );
        expect(spec(pairs_a[i].slice(15), pairs_b[j].slice(15))).toBeCloseTo(
          raw_spec2[i][j],
          10,
        );
      }
    }
  });
});
