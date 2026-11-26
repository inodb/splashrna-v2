/**
 * The two string kernels behind the SplashRNA classifiers.
 *
 * These reproduce Shogun 2.1.0's `WeightedDegreePositionStringKernel` and
 * `WeightedCommWordStringKernel` exactly; see docs/MODEL-ARCHAEOLOGY.md for how the
 * parameters were recovered and tools/reference.py for the Python mirror of this code.
 *
 * Nucleotides are represented as base indices 0..3 (A, C, G, T).
 */

export const DEGREE = 10;
export const SEQUENCE_LENGTH = 22;
export const SPLIT = 15;
export const SPECTRUM_ORDER = 3;

/** w_k = 2(d - k + 1) / (d(d + 1)) for k = 1..d */
export const WD_WEIGHTS: number[] = Array.from(
  { length: DEGREE },
  (_, i) => (2 * (DEGREE - i)) / (DEGREE * (DEGREE + 1)),
);

/**
 * WD_CUMULATIVE[k] is the total contribution of a match run of length k, i.e. the
 * sum of the first k degree weights. Indexing by run length turns the inner loop over
 * degrees into a single table lookup.
 */
export const WD_CUMULATIVE: Float64Array = (() => {
  const table = new Float64Array(DEGREE + 1);
  let sum = 0;
  for (let k = 1; k <= DEGREE; k++) {
    sum += WD_WEIGHTS[k - 1];
    table[k] = sum;
  }
  return table;
})();

/** Suffix-mask weights for spectrum orders 1, 2 and 3. */
export const SPECTRUM_WEIGHTS = [1 / 2, 1 / 3, 1 / 6] as const;

/** Shift weight beta_s = 1 / (s + 1); the models use a single shift of 1. */
const BETA_1 = 1 / 2;

/**
 * Unnormalized weighted-degree-position kernel between two sequences of equal length,
 * with a shift of one position in each direction.
 *
 * The three alignments (direct, `a` shifted forward, `b` shifted forward) are walked in
 * a single backward pass. Going right to left lets each position reuse the run length
 * computed for its successor, and holding the previous bases in locals keeps this to
 * two array reads per position instead of six. This function dominates the runtime of
 * the whole application, so the slightly awkward shape is deliberate.
 */
export function wdKernel(
  a: Uint8Array,
  offA: number,
  b: Uint8Array,
  offB: number,
  len: number = SEQUENCE_LENGTH,
): number {
  const last = len - 1;

  // Trailing position: only the direct alignment has a comparable pair there.
  let prevA = a[offA + last];
  let prevB = b[offB + last];
  let runDirect = prevA === prevB ? 1 : 0;
  let direct = WD_CUMULATIVE[runDirect];
  let shifted = 0;
  let runAShift = 0;
  let runBShift = 0;

  for (let j = last - 1; j >= 0; j--) {
    const curA = a[offA + j];
    const curB = b[offB + j];

    if (curA === curB) {
      runDirect = runDirect < DEGREE ? runDirect + 1 : DEGREE;
      direct += WD_CUMULATIVE[runDirect];
    } else {
      runDirect = 0;
    }

    // a[j + 1] against b[j]
    if (prevA === curB) {
      runAShift = runAShift < DEGREE ? runAShift + 1 : DEGREE;
      shifted += WD_CUMULATIVE[runAShift];
    } else {
      runAShift = 0;
    }

    // a[j] against b[j + 1]
    if (curA === prevB) {
      runBShift = runBShift < DEGREE ? runBShift + 1 : DEGREE;
      shifted += WD_CUMULATIVE[runBShift];
    } else {
      runBShift = 0;
    }

    prevA = curA;
    prevB = curB;
  }

  return direct + BETA_1 * shifted;
}

/**
 * Suffix-masked k-mer counts for a sequence slice, laid out as a single 84-element
 * vector: 4 one-mers, then 16 two-mers, then 64 three-mers.
 *
 * Shogun builds its order-3 words with `reverse=true`, which makes the degree masking
 * retain the LAST k characters of each word. Getting this backwards is the single
 * easiest way to break compatibility, so it is asserted by the spectrum-probe test.
 */
export const SPECTRUM_DIMS = 4 + 16 + 64;

export function spectrumCounts(
  seq: Uint8Array,
  offset: number,
  len: number,
  out?: Float64Array,
): Float64Array {
  const counts = out ?? new Float64Array(SPECTRUM_DIMS);
  if (out) counts.fill(0);

  const words = len - SPECTRUM_ORDER + 1;
  for (let p = 0; p < words; p++) {
    const c0 = seq[offset + p];
    const c1 = seq[offset + p + 1];
    const c2 = seq[offset + p + 2];
    counts[c2] += 1; // suffix of length 1
    counts[4 + (c1 * 4 + c2)] += 1; // suffix of length 2
    counts[20 + (c0 * 16 + c1 * 4 + c2)] += 1; // full word
  }
  return counts;
}

/** Weighted inner product of two spectrum count vectors. */
export function spectrumDot(a: Float64Array, b: Float64Array): number {
  let s1 = 0;
  let s2 = 0;
  let s3 = 0;
  for (let i = 0; i < 4; i++) s1 += a[i] * b[i];
  for (let i = 4; i < 20; i++) s2 += a[i] * b[i];
  for (let i = 20; i < SPECTRUM_DIMS; i++) s3 += a[i] * b[i];
  return (
    SPECTRUM_WEIGHTS[0] * s1 + SPECTRUM_WEIGHTS[1] * s2 + SPECTRUM_WEIGHTS[2] * s3
  );
}
