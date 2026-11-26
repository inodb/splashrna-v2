/** Nucleotide encoding, reverse complement, tiling and the miR-30 / miR-E constructs. */

import { SEQUENCE_LENGTH } from "./kernels.js";

export const BASES = "ACGT";

/** Base index, or -1 for anything that is not an unambiguous ACGT nucleotide. */
const BASE_INDEX = new Int8Array(128).fill(-1);
for (let i = 0; i < BASES.length; i++) {
  BASE_INDEX[BASES.charCodeAt(i)] = i;
  BASE_INDEX[BASES.toLowerCase().charCodeAt(i)] = i;
}

const COMPLEMENT: Record<string, string> = {
  A: "T", C: "G", G: "C", T: "A", N: "N",
  a: "T", c: "G", g: "C", t: "A", n: "N",
};

export function isPureACGT(seq: string): boolean {
  for (let i = 0; i < seq.length; i++) {
    if (BASE_INDEX[seq.charCodeAt(i)] < 0) return false;
  }
  return seq.length > 0;
}

export function encode(seq: string, out?: Uint8Array): Uint8Array {
  const arr = out ?? new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const code = BASE_INDEX[seq.charCodeAt(i)];
    if (code < 0) throw new Error(`Cannot encode non-ACGT character at position ${i}`);
    arr[i] = code;
  }
  return arr;
}

export function reverseComplement(seq: string): string {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    out += COMPLEMENT[seq[i]] ?? "N";
  }
  return out;
}

export interface Tile {
  /** Antisense guide: the reverse complement of the sense 22-mer. */
  guide: string;
  /** 1-based start of the sense 22-mer within the supplied transcript. */
  position: number;
}

/**
 * Tile a transcript into overlapping antisense 22-mers.
 *
 * Matches the original server's behaviour: it walks every start position and takes the
 * reverse complement, so the guide strand is what gets scored.
 */
export function tileTranscript(sequence: string): Tile[] {
  const seq = sequence.toUpperCase().replace(/\s/g, "");
  const tiles: Tile[] = [];
  for (let i = 0; i + SEQUENCE_LENGTH <= seq.length; i++) {
    tiles.push({
      guide: reverseComplement(seq.slice(i, i + SEQUENCE_LENGTH)),
      position: i + 1,
    });
  }
  return tiles;
}

// --- Hairpin constructs --------------------------------------------------------------

const MIR30_START = "TGCTGTTGACAGTGAGCG";
const MIR30_LOOP = "TAGTGAAGCCACAGATGTA";
const MIR30_END = "TGCCTACTGCCTCGGA";

/**
 * Build the 97-mer miR-30 construct for a guide sequence.
 *
 * A deliberate mismatch is introduced at the first passenger position, exactly as the
 * original server did: a C/G there becomes A, anything else becomes C.
 */
export function mir30Construct(guide: string): string {
  const passenger = reverseComplement(guide);
  const first = passenger[0] === "C" || passenger[0] === "G" ? "A" : "C";
  return MIR30_START + first + passenger.slice(1) + MIR30_LOOP + guide + MIR30_END;
}

/** Restriction sites that make a construct difficult to clone. */
export const RESTRICTION_SITES: Record<string, string> = {
  EcoRI: "GAATTC",
  XhoI: "CTCGAG",
};

/** Which restriction sites (if any) appear in the 97-mer for this guide. */
export function restrictionSites(guide: string): string[] {
  const construct = mir30Construct(guide);
  return Object.entries(RESTRICTION_SITES)
    .filter(([, site]) => construct.includes(site))
    .map(([name]) => name);
}
