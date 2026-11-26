/** Loading and decoding of the SplashRNA model bundle produced by tools/build_model.py. */

import { CollapsedWd } from "./collapsed.js";
import { SEQUENCE_LENGTH, SPECTRUM_DIMS } from "./kernels.js";
import { encode } from "./sequence.js";

export interface RawSubModel {
  bias: number;
  numSupportVectors: number;
  supportVectors: string;
  coefWd: string;
  spectrumLeft: number[];
  spectrumRight: number[];
}

export interface RawModelBundle {
  formatVersion: number;
  source: string;
  citation: string;
  kernel: {
    sequenceLength: number;
    split: number;
    wdDegree: number;
    wdShift: number;
    spectrumOrder: number;
    spectrumWeights: number[];
  };
  combine: { threshold: number; weight: number };
  models: { mir30: RawSubModel; mirE: RawSubModel };
}

/** Decoded, ready-to-score form of one sub-classifier. */
export interface SubModel {
  bias: number;
  count: number;
  /** Base indices for every support vector, concatenated, 22 per vector. */
  supportVectors: Uint8Array;
  /** alpha_i / ||sv_i|| for the weighted-degree kernel. */
  coefWd: Float64Array;
  /** Pre-collapsed lookup tables for the weighted-degree kernel sum. */
  wd: CollapsedWd;
  /** Collapsed spectrum weights for nucleotides 1..15 and 16..22. */
  spectrumLeft: Float64Array;
  spectrumRight: Float64Array;
}

export interface SplashModel {
  citation: string;
  source: string;
  combine: { threshold: number; weight: number };
  mir30: SubModel;
  mirE: SubModel;
}

function decodeBase64Float64(encoded: string): Float64Array {
  const binary = atob(encoded);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return new Float64Array(bytes.buffer);
}

function decodeSubModel(raw: RawSubModel): SubModel {
  const count = raw.numSupportVectors;
  const expected = count * SEQUENCE_LENGTH;
  if (raw.supportVectors.length !== expected) {
    throw new Error(
      `Support vector block is ${raw.supportVectors.length} characters, expected ${expected}`,
    );
  }

  const supportVectors = encode(raw.supportVectors);
  const coefWd = decodeBase64Float64(raw.coefWd);
  if (coefWd.length !== count) {
    throw new Error(`Got ${coefWd.length} coefficients for ${count} support vectors`);
  }
  for (const vec of [raw.spectrumLeft, raw.spectrumRight]) {
    if (vec.length !== SPECTRUM_DIMS) {
      throw new Error(`Spectrum weights must have ${SPECTRUM_DIMS} entries`);
    }
  }

  return {
    bias: raw.bias,
    count,
    supportVectors,
    coefWd,
    wd: new CollapsedWd(supportVectors, coefWd, count),
    spectrumLeft: Float64Array.from(raw.spectrumLeft),
    spectrumRight: Float64Array.from(raw.spectrumRight),
  };
}

export function decodeModel(raw: RawModelBundle): SplashModel {
  if (raw.formatVersion !== 1) {
    throw new Error(`Unsupported model format version ${raw.formatVersion}`);
  }
  return {
    citation: raw.citation,
    source: raw.source,
    combine: raw.combine,
    mir30: decodeSubModel(raw.models.mir30),
    mirE: decodeSubModel(raw.models.mirE),
  };
}

/** Fetch and decode the model bundle. */
export async function loadModel(url = "splashrna-model.json"): Promise<SplashModel> {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Could not load the SplashRNA model (HTTP ${response.status})`);
  }
  return decodeModel((await response.json()) as RawModelBundle);
}
