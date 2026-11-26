/**
 * The shRNA design pipeline: tile a transcript, score every tile, filter, rank.
 *
 * This is the business logic that used to live in `Splash.calculate()` on the old
 * CherryPy server, restated so each step is independently testable.
 */

import { SEQUENCE_LENGTH } from "./kernels.js";
import type { SplashScorer } from "./scorer.js";
import {
  isPureACGT,
  mir30Construct,
  restrictionSites,
  tileTranscript,
} from "./sequence.js";

/**
 * How thoroughly to apply the off-target filter.
 *
 * The original server looked up every tile in the genome-wide match tables before
 * truncating to the requested number of results. That is exact but expensive, because
 * it touches thousands of tiles to return five.
 *
 * `"top-candidates"` ranks first and only looks up the best `offTargetCandidates`
 * tiles. In practice this returns the same list, because a tile would have to be
 * knocked out of the top N and replaced by one ranked below the candidate cut-off to
 * differ. `"exhaustive"` reproduces the original ordering exactly.
 */
export type OffTargetMode = "top-candidates" | "exhaustive";

export interface OffTargetHits {
  /**
   * Number of genes whose transcripts contain this 22-mer. The legacy tables only
   * record 22-mers occurring in two or more genes, so any non-zero count means the
   * guide is not unique.
   */
  humanGenes: number;
  mouseGenes: number;
  /**
   * Entrez gene IDs behind those counts, filled in only for guides that are displayed.
   * Absent when the index was not asked for them.
   */
  humanIds?: number[];
  mouseIds?: number[];
}

export const NO_OFF_TARGET_HITS: OffTargetHits = { humanGenes: 0, mouseGenes: 0 };

export function hasOffTargetHits(hits: OffTargetHits | null): boolean {
  return hits !== null && hits.humanGenes + hits.mouseGenes > 0;
}

/** Source of genome-wide 22-mer match data. */
export interface OffTargetIndex {
  lookup(guides: string[]): Promise<Map<string, OffTargetHits>>;
  /**
   * Optional second pass filling in Entrez IDs. Separate from `lookup` because counts
   * are needed for every candidate while identities are only needed for the few shown.
   */
  annotate?(hits: Map<string, OffTargetHits>): Promise<void>;
}

export interface DesignOptions {
  /** How many shRNAs to return per input sequence. */
  predictionsPerGene: number;
  /** Drop guides whose 97-mer contains an EcoRI or XhoI site. */
  removeRestrictionSites: boolean;
  /** Keep guides that occur more than once within the same transcript. */
  allowInGeneDuplicates: boolean;
  /** Keep guides that exactly match a different human or mouse transcript. */
  allowOffTargetDuplicates: boolean;
  offTargetMode: OffTargetMode;
  /** Tiles to off-target check in "top-candidates" mode. */
  offTargetCandidates: number;
}

export const DEFAULT_OPTIONS: DesignOptions = {
  predictionsPerGene: 5,
  removeRestrictionSites: true,
  allowInGeneDuplicates: false,
  allowOffTargetDuplicates: false,
  offTargetMode: "top-candidates",
  offTargetCandidates: 200,
};

export interface DesignedGuide {
  guide: string;
  position: number;
  mir30: number;
  mirE: number;
  splash: number;
  construct: string;
  warnings: string[];
  /** Null when no off-target index was supplied, so "unknown" is distinct from "none". */
  offTarget: OffTargetHits | null;
}

export interface DesignResult {
  label: string;
  sequenceLength: number;
  guides: DesignedGuide[];
  warnings: string[];
  /** Tiles discarded at each stage, for transparency in the UI. */
  discarded: {
    ambiguousBases: number;
    restrictionSites: number;
    inGeneDuplicates: number;
    offTarget: number;
  };
}

export interface DesignInput {
  label: string;
  sequence: string;
}

const SCORE_WARNING_THRESHOLD = 1.0;
const SHORT_TRANSCRIPT_LENGTH = 1000;

export async function designGuides(
  input: DesignInput,
  scorer: SplashScorer,
  options: Partial<DesignOptions> = {},
  offTargetIndex?: OffTargetIndex,
  onProgress?: (fraction: number) => void,
): Promise<DesignResult> {
  const opts = { ...DEFAULT_OPTIONS, ...options };
  const sequence = input.sequence.toUpperCase().replace(/\s/g, "");
  const discarded = {
    ambiguousBases: 0,
    restrictionSites: 0,
    inGeneDuplicates: 0,
    offTarget: 0,
  };
  const warnings: string[] = [];

  if (sequence.length < SEQUENCE_LENGTH) {
    return {
      label: input.label,
      sequenceLength: sequence.length,
      guides: [],
      warnings: [
        `Sequence is only ${sequence.length} nt. At least ${SEQUENCE_LENGTH} nt are needed to design an shRNA.`,
      ],
      discarded,
    };
  }

  // 1. Tile, dropping anything containing N or other ambiguous bases.
  let tiles = tileTranscript(sequence);
  const total = tiles.length;
  tiles = tiles.filter((t) => isPureACGT(t.guide));
  discarded.ambiguousBases = total - tiles.length;
  if (discarded.ambiguousBases > 0) {
    warnings.push(
      `${discarded.ambiguousBases} tile(s) containing ambiguous bases were removed.`,
    );
  }

  // 2. Restriction sites, if the user wants them gone.
  if (opts.removeRestrictionSites) {
    const kept = tiles.filter((t) => restrictionSites(t.guide).length === 0);
    discarded.restrictionSites = tiles.length - kept.length;
    tiles = kept;
  }

  // 3. In-gene duplicates.
  const occurrences = new Map<string, number>();
  for (const t of tiles) occurrences.set(t.guide, (occurrences.get(t.guide) ?? 0) + 1);
  if (!opts.allowInGeneDuplicates) {
    const kept = tiles.filter((t) => occurrences.get(t.guide) === 1);
    discarded.inGeneDuplicates = tiles.length - kept.length;
    tiles = kept;
  }

  if (tiles.length === 0) {
    warnings.push("Every candidate was removed by the selected filters.");
    return { label: input.label, sequenceLength: sequence.length, guides: [], warnings, discarded };
  }

  // 4. Score and rank.
  const scores = scorer.scoreAll(tiles.map((t) => t.guide), onProgress);
  let ranked: DesignedGuide[] = tiles.map((tile, i) => ({
    guide: tile.guide,
    position: tile.position,
    mir30: scores[i].mir30,
    mirE: scores[i].mirE,
    splash: scores[i].splash,
    construct: mir30Construct(tile.guide),
    warnings: [
      ...(occurrences.get(tile.guide)! > 1
        ? ["Occurs more than once in this transcript."]
        : []),
      ...restrictionSites(tile.guide).map((s) => `${s} site in the 97-mer construct.`),
    ],
    offTarget: null,
  }));
  ranked.sort((a, b) => b.splash - a.splash);

  // 5. Off-target filtering, over either the whole ranked list or just the top slice.
  if (offTargetIndex) {
    const pool =
      opts.offTargetMode === "exhaustive"
        ? ranked
        : ranked.slice(0, Math.max(opts.offTargetCandidates, opts.predictionsPerGene));
    const hits = await offTargetIndex.lookup(pool.map((g) => g.guide));
    for (const guide of pool) {
      guide.offTarget = hits.get(guide.guide) ?? NO_OFF_TARGET_HITS;
      if (hasOffTargetHits(guide.offTarget)) {
        const { humanGenes, mouseGenes } = guide.offTarget;
        const parts = [
          humanGenes > 0 ? `${humanGenes} human gene${humanGenes === 1 ? "" : "s"}` : "",
          mouseGenes > 0 ? `${mouseGenes} mouse gene${mouseGenes === 1 ? "" : "s"}` : "",
        ].filter(Boolean);
        guide.warnings.push(`Exact 22-mer match in ${parts.join(" and ")}.`);
      }
    }
    if (!opts.allowOffTargetDuplicates) {
      const before = ranked.length;
      ranked = ranked.filter((g) => !hasOffTargetHits(g.offTarget));
      discarded.offTarget = before - ranked.length;
    }
  }

  const guides = ranked.slice(0, opts.predictionsPerGene);

  // Name the off-target genes, but only for the handful being returned.
  if (offTargetIndex?.annotate) {
    const shown = new Map<string, OffTargetHits>();
    for (const g of guides) if (g.offTarget) shown.set(g.guide, g.offTarget);
    if (shown.size > 0) await offTargetIndex.annotate(shown);
  }

  if (guides.length > 0 && guides[0].splash < SCORE_WARNING_THRESHOLD) {
    warnings.push(
      `The best predicted score is below ${SCORE_WARNING_THRESHOLD}. You may need to test more shRNAs from this list to find a potent one.`,
    );
    if (sequence.length < SHORT_TRANSCRIPT_LENGTH) {
      warnings.push(
        "The transcript is shorter than 1 kb. Choosing a single isoform may improve results.",
      );
    }
  }

  return { label: input.label, sequenceLength: sequence.length, guides, warnings, discarded };
}
