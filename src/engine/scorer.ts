/**
 * SplashRNA scoring.
 *
 * The published algorithm is sequential: a miR-30 classifier gates, and the miR-E
 * classifier only contributes for candidates that clear the threshold.
 *
 *   Pelossof R, Fairchild L, et al. Prediction of potent shRNAs with a sequential
 *   classification algorithm. Nat Biotechnol. 2017;35(4):350-353.
 */

import {
  SEQUENCE_LENGTH,
  SPECTRUM_DIMS,
  SPLIT,
  spectrumCounts,
  spectrumDot,
  wdKernel,
} from "./kernels.js";
import type { SplashModel, SubModel } from "./model.js";
import { encode } from "./sequence.js";

/** Scratch buffers, reused across every tile so scoring allocates nothing per call. */
interface Scratch {
  query: Uint8Array;
  left: Float64Array;
  right: Float64Array;
}

function makeScratch(): Scratch {
  return {
    query: new Uint8Array(SEQUENCE_LENGTH),
    left: new Float64Array(SPECTRUM_DIMS),
    right: new Float64Array(SPECTRUM_DIMS),
  };
}

function scoreSubModel(model: SubModel, scratch: Scratch): number {
  const { query, left, right } = scratch;

  // Weighted-degree term, via the collapsed lookup tables. The self-normalizer is
  // still an ordinary kernel evaluation, since it involves only the query.
  const normWd = Math.sqrt(wdKernel(query, 0, query, 0, SEQUENCE_LENGTH));
  const wdSum = model.wd.apply(query);

  // Spectrum terms over the two halves of the guide.
  spectrumCounts(query, 0, SPLIT, left);
  spectrumCounts(query, SPLIT, SEQUENCE_LENGTH - SPLIT, right);
  const normLeft = Math.sqrt(spectrumDot(left, left));
  const normRight = Math.sqrt(spectrumDot(right, right));

  return (
    wdSum / normWd +
    spectrumDot(model.spectrumLeft, left) / normLeft +
    spectrumDot(model.spectrumRight, right) / normRight +
    model.bias
  );
}

export interface GuideScore {
  mir30: number;
  mirE: number;
  splash: number;
}

export class SplashScorer {
  private readonly scratch = makeScratch();

  constructor(private readonly model: SplashModel) {}

  get citation(): string {
    return this.model.citation;
  }

  /** Score a single 22-nucleotide antisense guide. */
  score(guide: string): GuideScore {
    if (guide.length !== SEQUENCE_LENGTH) {
      throw new Error(`Guide must be ${SEQUENCE_LENGTH} nt, got ${guide.length}`);
    }
    encode(guide, this.scratch.query);

    const mir30 = scoreSubModel(this.model.mir30, this.scratch);
    const mirE = scoreSubModel(this.model.mirE, this.scratch);
    return { mir30, mirE, splash: this.combine(mir30, mirE) };
  }

  /**
   * Score many guides, optionally reporting progress.
   *
   * `onProgress` is called with a fraction in [0, 1]; the web UI uses it to drive a
   * progress bar, since a long transcript takes a noticeable moment.
   */
  scoreAll(guides: string[], onProgress?: (fraction: number) => void): GuideScore[] {
    const out: GuideScore[] = new Array(guides.length);
    const step = Math.max(1, Math.floor(guides.length / 100));
    for (let i = 0; i < guides.length; i++) {
      out[i] = this.score(guides[i]);
      if (onProgress && i % step === 0) onProgress(i / guides.length);
    }
    onProgress?.(1);
    return out;
  }

  /** The sequential combination rule from the paper. */
  combine(mir30: number, mirE: number): number {
    const { threshold, weight } = this.model.combine;
    return mir30 > threshold ? mir30 * weight + mirE * (1 - weight) : mir30 * weight;
  }
}
