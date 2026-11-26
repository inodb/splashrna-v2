/**
 * Collapsed form of the weighted-degree kernel sum.
 *
 * Scoring needs only the weighted total over support vectors,
 *
 *     S(x) = sum_i coef_i * K(sv_i, x)
 *
 * and never the individual kernel values. Expanding the run-length formulation back
 * into per-degree match indicators,
 *
 *     K(sv, x) = sum_j sum_k w_k * [ sv[j..j+k] == x[j..j+k] ]
 *
 * lets the sum over support vectors move inside:
 *
 *     S(x) = sum_j sum_k T[j][k][ code(x, j, k) ]        where
 *     T[j][k][c] = w_k * sum_i coef_i * [ code(sv_i, j, k) == c ]
 *
 * T is built once at load time, so scoring a guide costs a few hundred table lookups
 * instead of a pass over every support vector. For the miR-30 model that is a few
 * hundred lookups in place of 11,585 kernel evaluations.
 *
 * Both shifted alignments reuse the same table, because the table for "support vector
 * k-mer at position j+1" is just T[j+1].
 *
 * Degrees 1..7 are dense arrays (4^k entries). Degrees 8..10 are rare enough that they
 * are handled by looking up an 8-mer index and checking the handful of candidates
 * directly, which keeps memory to a few megabytes instead of a few hundred.
 */

import { DEGREE, SEQUENCE_LENGTH, WD_WEIGHTS } from "./kernels.js";

const DENSE_MAX_DEGREE = 7;
const INDEX_DEGREE = 8;
const BETA_1 = 1 / 2;

/** Base-4 code of the k-mer starting at `pos`. */
function kmerCode(seq: Uint8Array, offset: number, pos: number, k: number): number {
  let code = 0;
  for (let t = 0; t < k; t++) code = code * 4 + seq[offset + pos + t];
  return code;
}

export class CollapsedWd {
  /** Dense tables indexed by [position * (DENSE_MAX_DEGREE + 1) + degree]. */
  private readonly dense: (Float64Array | null)[];
  /** Per position, a map from 8-mer code to the support vectors carrying it. */
  private readonly index8: Map<number, Int32Array>[];

  constructor(
    private readonly supportVectors: Uint8Array,
    private readonly coef: Float64Array,
    private readonly count: number,
  ) {
    const stride = DENSE_MAX_DEGREE + 1;
    this.dense = new Array(SEQUENCE_LENGTH * stride).fill(null);
    for (let j = 0; j < SEQUENCE_LENGTH; j++) {
      for (let k = 1; k <= DENSE_MAX_DEGREE && j + k <= SEQUENCE_LENGTH; k++) {
        this.dense[j * stride + k] = new Float64Array(4 ** k);
      }
    }

    const builders: Map<number, number[]>[] = Array.from(
      { length: SEQUENCE_LENGTH },
      () => new Map<number, number[]>(),
    );

    for (let i = 0; i < count; i++) {
      const base = i * SEQUENCE_LENGTH;
      const weight = coef[i];
      for (let j = 0; j < SEQUENCE_LENGTH; j++) {
        const maxK = Math.min(INDEX_DEGREE, SEQUENCE_LENGTH - j);
        let code = 0;
        for (let k = 1; k <= maxK; k++) {
          code = code * 4 + supportVectors[base + j + k - 1];
          if (k <= DENSE_MAX_DEGREE) {
            // Fold the degree weight in now so scoring is pure addition.
            this.dense[j * stride + k]![code] += weight * WD_WEIGHTS[k - 1];
          } else {
            const bucket = builders[j].get(code);
            if (bucket) bucket.push(i);
            else builders[j].set(code, [i]);
          }
        }
      }
    }

    this.index8 = builders.map((m) => {
      const out = new Map<number, Int32Array>();
      for (const [code, list] of m) out.set(code, Int32Array.from(list));
      return out;
    });
  }

  /** Dense contribution for a support-vector position and a query position. */
  private denseAt(svPos: number, q: Uint8Array, qPos: number, maxK: number): number {
    const stride = DENSE_MAX_DEGREE + 1;
    const limit = Math.min(DENSE_MAX_DEGREE, maxK);
    let total = 0;
    let code = 0;
    for (let k = 1; k <= limit; k++) {
      code = code * 4 + q[qPos + k - 1];
      total += this.dense[svPos * stride + k]![code];
    }
    return total;
  }

  /**
   * Contribution of degrees 8..10, for which we check candidate support vectors
   * directly. A query shares an 8-mer with only a handful of them.
   */
  private longAt(svPos: number, q: Uint8Array, qPos: number, maxK: number): number {
    if (maxK < INDEX_DEGREE) return 0;
    const candidates = this.index8[svPos].get(kmerCode(q, 0, qPos, INDEX_DEGREE));
    if (!candidates) return 0;

    let total = 0;
    for (let n = 0; n < candidates.length; n++) {
      const i = candidates[n];
      const base = i * SEQUENCE_LENGTH;
      let add = WD_WEIGHTS[INDEX_DEGREE - 1];
      for (let k = INDEX_DEGREE + 1; k <= Math.min(DEGREE, maxK); k++) {
        if (this.supportVectors[base + svPos + k - 1] !== q[qPos + k - 1]) break;
        add += WD_WEIGHTS[k - 1];
      }
      total += this.coef[i] * add;
    }
    return total;
  }

  private at(svPos: number, q: Uint8Array, qPos: number, maxK: number): number {
    return this.denseAt(svPos, q, qPos, maxK) + this.longAt(svPos, q, qPos, maxK);
  }

  /** sum_i coef_i * K(sv_i, query), identical to looping wdKernel over every vector. */
  apply(q: Uint8Array): number {
    let direct = 0;
    for (let j = 0; j < SEQUENCE_LENGTH; j++) {
      direct += this.at(j, q, j, SEQUENCE_LENGTH - j);
    }

    let shifted = 0;
    for (let j = 0; j < SEQUENCE_LENGTH - 1; j++) {
      // A k-mer match in either shifted alignment needs j + 1 + k <= SEQUENCE_LENGTH.
      const maxK = SEQUENCE_LENGTH - 1 - j;
      if (maxK < 1) continue;
      shifted += this.at(j + 1, q, j, maxK); // support vector shifted forward
      shifted += this.at(j, q, j + 1, maxK); // query shifted forward
    }

    return direct + BETA_1 * shifted;
  }

  get supportVectorCount(): number {
    return this.count;
  }
}
