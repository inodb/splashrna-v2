/**
 * Static off-target index.
 *
 * Answers "does this 22-mer occur in more than one gene?" against the human and mouse
 * tables, with no server. The index is a block-structured file built by
 * tools/build_offtarget.py; see that file for the byte layout.
 *
 * Only the block directories (~21 KB) are fetched up front. A lookup binary-searches
 * the directory in memory and then fetches just the block it needs with a Range
 * request. When a query touches a large fraction of the file -- which happens in
 * exhaustive mode -- the whole file is fetched instead, since that is cheaper than
 * hundreds of ranges.
 *
 * A 22-mer packs into 44 bits, which is well inside the 53 bits a JavaScript number
 * represents exactly, so everything here uses plain numbers. That does mean arithmetic
 * rather than bitwise operators, because JavaScript's bitwise operators truncate to 32
 * bits and would silently corrupt these values.
 */

import type { OffTargetHits, OffTargetIndex } from "../engine/design.js";
import { reverseComplement } from "../engine/sequence.js";

const BASE_INDEX: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

/** Pack an ACGT 22-mer into a 44-bit integer, or null if it has ambiguous bases. */
export function packKmer(seq: string): number | null {
  let value = 0;
  for (let i = 0; i < seq.length; i++) {
    const base = BASE_INDEX[seq[i]];
    if (base === undefined) return null;
    value = value * 4 + base;
  }
  return value;
}

interface GenomeManifest {
  assembly: string;
  count: number;
  blocks: number;
  keysBytes: number;
  dirBytes: number;
  keys: string;
  directory: string;
  genes: string;
  geneIds: string;
  geneCount: number;
  geneRefs: number;
  dirRecordBytes: number;
}

interface Manifest {
  formatVersion: number;
  blockSize: number;
  genomes: Record<string, GenomeManifest>;
}

/**
 * One entry of a decoded block: the key, how many genes carry it, and where that
 * gene list starts. The start is derived by running the counts forward from the
 * block's own starting index, so no per-key offset table has to be stored.
 */
interface BlockEntry {
  key: number;
  geneCount: number;
  geneStart: number;
}

function decodeBlock(view: Uint8Array, limit: number, geneStart: number): BlockEntry[] {
  const entries: BlockEntry[] = [];
  let pos = 0;
  let geneCursor = geneStart;

  // First key: six big-endian bytes. Multiplication, not shifts -- see the file note.
  let key = 0;
  for (let i = 0; i < 6; i++) key = key * 256 + view[pos++];

  const readVarint = (): number => {
    let result = 0;
    let multiplier = 1;
    for (;;) {
      const byte = view[pos++];
      result += (byte & 0x7f) * multiplier;
      if ((byte & 0x80) === 0) return result;
      multiplier *= 128;
    }
  };

  let count = readVarint();
  entries.push({ key, geneCount: count, geneStart: geneCursor });
  geneCursor += count;

  while (entries.length < limit && pos < view.length) {
    key += readVarint();
    count = readVarint();
    entries.push({ key, geneCount: count, geneStart: geneCursor });
    geneCursor += count;
  }
  return entries;
}

class Genome {
  /** First key of each block, ascending. */
  private readonly firstKeys: Float64Array;
  /** Byte offset of each block, plus a final sentinel of the file length. */
  private readonly offsets: Float64Array;
  /** Index of each block's first gene reference in the genes payload. */
  private readonly geneStarts: Float64Array;
  /** Entrez IDs the uint16 gene references point at; loaded on first use. */
  private geneIds: number[] | null = null;
  private readonly cache = new Map<number, BlockEntry[]>();
  private whole: Uint8Array | null = null;

  constructor(
    private readonly baseUrl: string,
    private readonly manifest: GenomeManifest,
    private readonly blockSize: number,
    directory: ArrayBuffer,
  ) {
    const view = new DataView(directory);
    const blocks = manifest.blocks;
    const stride = manifest.dirRecordBytes;
    this.firstKeys = new Float64Array(blocks);
    this.offsets = new Float64Array(blocks + 1);
    this.geneStarts = new Float64Array(blocks);

    for (let i = 0; i < blocks; i++) {
      const base = i * stride;
      this.offsets[i] = view.getUint32(base, false);
      this.geneStarts[i] = view.getUint32(base + 4, false);
      let key = 0;
      for (let b = 0; b < 6; b++) key = key * 256 + view.getUint8(base + 8 + b);
      this.firstKeys[i] = key;
    }
    this.offsets[blocks] = manifest.keysBytes;
  }

  get assembly(): string {
    return this.manifest.assembly;
  }

  /** Index of the block that would contain `key`, or -1 if it precedes the first. */
  private blockFor(key: number): number {
    let low = 0;
    let high = this.firstKeys.length - 1;
    if (key < this.firstKeys[0]) return -1;
    while (low < high) {
      const mid = (low + high + 1) >> 1;
      if (this.firstKeys[mid] <= key) low = mid;
      else high = mid - 1;
    }
    return low;
  }

  private async fetchBlock(index: number): Promise<BlockEntry[]> {
    const cached = this.cache.get(index);
    if (cached) return cached;

    const start = this.offsets[index];
    const end = this.offsets[index + 1];
    let bytes: Uint8Array;

    if (this.whole) {
      bytes = this.whole.subarray(start, end);
    } else {
      const response = await fetch(`${this.baseUrl}/${this.manifest.keys}`, {
        headers: { Range: `bytes=${start}-${end - 1}` },
      });
      if (!response.ok) {
        throw new Error(`Off-target index fetch failed (HTTP ${response.status})`);
      }
      const buffer = new Uint8Array(await response.arrayBuffer());
      // A host that ignores Range returns 200 and the entire file; keep it.
      if (response.status === 200 && buffer.length === this.manifest.keysBytes) {
        this.whole = buffer;
        bytes = buffer.subarray(start, end);
      } else {
        bytes = buffer;
      }
    }

    const entries = decodeBlock(bytes, this.blockSize, this.geneStarts[index]);
    this.cache.set(index, entries);
    return entries;
  }

  /** Download the whole keys file, for queries that would touch most of it anyway. */
  private async fetchWhole(): Promise<void> {
    if (this.whole) return;
    const response = await fetch(`${this.baseUrl}/${this.manifest.keys}`);
    if (!response.ok) {
      throw new Error(`Off-target index fetch failed (HTTP ${response.status})`);
    }
    this.whole = new Uint8Array(await response.arrayBuffer());
  }

  /** The Entrez ID table, fetched once and only when gene names are first asked for. */
  private async loadGeneIds(): Promise<number[]> {
    if (this.geneIds) return this.geneIds;
    const response = await fetch(`${this.baseUrl}/${this.manifest.geneIds}`);
    if (!response.ok) {
      throw new Error(`Could not load ${this.assembly} gene IDs (HTTP ${response.status})`);
    }
    this.geneIds = (await response.json()) as number[];
    return this.geneIds;
  }

  /**
   * Entrez IDs of the genes carrying each key.
   *
   * Only called for guides that are actually displayed, so this fetches a few dozen
   * bytes per guide rather than the 6 MB payload.
   */
  async genesFor(keys: number[]): Promise<Map<number, number[]>> {
    const wanted: { key: number; entry: BlockEntry }[] = [];
    for (const key of keys) {
      const block = this.blockFor(key);
      if (block < 0) continue;
      const entries = this.cache.get(block);
      const entry = entries?.find((e) => e.key === key);
      if (entry && entry.geneCount > 0) wanted.push({ key, entry });
    }
    if (wanted.length === 0) return new Map();

    const ids = await this.loadGeneIds();
    const out = new Map<number, number[]>();

    await Promise.all(
      wanted.map(async ({ key, entry }) => {
        const start = entry.geneStart * 2;
        const end = start + entry.geneCount * 2;
        const response = await fetch(`${this.baseUrl}/${this.manifest.genes}`, {
          headers: { Range: `bytes=${start}-${end - 1}` },
        });
        if (!response.ok) return;
        let bytes = new Uint8Array(await response.arrayBuffer());
        if (response.status === 200 && bytes.length === this.manifest.geneRefs * 2) {
          bytes = bytes.subarray(start, end);
        }
        const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
        const genes: number[] = [];
        for (let i = 0; i < entry.geneCount; i++) {
          const index = view.getUint16(i * 2, false);
          if (index < ids.length) genes.push(ids[index]);
        }
        out.set(key, genes);
      }),
    );
    return out;
  }

  /**
   * Gene counts for the given packed keys. A key absent from the table occurs in at
   * most one gene and so is not an off-target risk.
   */
  async lookup(keys: number[]): Promise<Map<number, number>> {
    const needed = new Set<number>();
    for (const key of keys) {
      const block = this.blockFor(key);
      if (block >= 0 && !this.cache.has(block)) needed.add(block);
    }

    // Past roughly a quarter of the file, one sequential download beats many ranges.
    if (needed.size > this.manifest.blocks / 4) await this.fetchWhole();

    await Promise.all([...needed].map((index) => this.fetchBlock(index)));

    const found = new Map<number, number>();
    for (const key of keys) {
      const block = this.blockFor(key);
      if (block < 0) continue;
      const entries = this.cache.get(block);
      if (!entries) continue;
      // Blocks hold at most blockSize entries, so a linear scan is fine and avoids
      // another binary search; most lookups miss anyway.
      for (const entry of entries) {
        if (entry.key === key) {
          found.set(key, entry.geneCount);
          break;
        }
        if (entry.key > key) break;
      }
    }
    return found;
  }
}

/**
 * Off-target index backed by static files.
 *
 * The legacy tables are keyed by the *sense* 22-mer, so a guide is looked up by its
 * reverse complement -- matching `getMatches()` in the original server.
 */
export class StaticOffTargetIndex implements OffTargetIndex {
  private constructor(
    private readonly human: Genome,
    private readonly mouse: Genome,
  ) {}

  static async load(baseUrl = "./offtarget"): Promise<StaticOffTargetIndex> {
    const response = await fetch(`${baseUrl}/manifest.json`);
    if (!response.ok) {
      throw new Error(`No off-target index at ${baseUrl} (HTTP ${response.status})`);
    }
    const manifest = (await response.json()) as Manifest;
    if (manifest.formatVersion !== 2) {
      throw new Error(
        `Unsupported off-target index version ${manifest.formatVersion}; rebuild with tools/build_offtarget.py`,
      );
    }

    const genomes = await Promise.all(
      (["human", "mouse"] as const).map(async (name) => {
        const entry = manifest.genomes[name];
        if (!entry) throw new Error(`Off-target index is missing the ${name} genome`);
        const dir = await fetch(`${baseUrl}/${entry.directory}`);
        if (!dir.ok) {
          throw new Error(`Could not load ${name} block directory (HTTP ${dir.status})`);
        }
        return new Genome(baseUrl, entry, manifest.blockSize, await dir.arrayBuffer());
      }),
    );

    return new StaticOffTargetIndex(genomes[0], genomes[1]);
  }

  async lookup(guides: string[]): Promise<Map<string, OffTargetHits>> {
    const packed = new Map<string, number>();
    for (const guide of guides) {
      const key = packKmer(reverseComplement(guide));
      if (key !== null) packed.set(guide, key);
    }
    const keys = [...new Set(packed.values())];

    const [human, mouse] = await Promise.all([
      this.human.lookup(keys),
      this.mouse.lookup(keys),
    ]);

    const result = new Map<string, OffTargetHits>();
    for (const [guide, key] of packed) {
      result.set(guide, {
        humanGenes: human.get(key) ?? 0,
        mouseGenes: mouse.get(key) ?? 0,
      });
    }
    return result;
  }

  /**
   * Fill in the Entrez IDs behind the counts, for a handful of guides.
   *
   * Separate from `lookup` because the counts are needed for every candidate while the
   * IDs are needed only for the few that get shown.
   */
  async annotate(hits: Map<string, OffTargetHits>): Promise<void> {
    const keys = new Map<string, number>();
    for (const guide of hits.keys()) {
      const key = packKmer(reverseComplement(guide));
      if (key !== null) keys.set(guide, key);
    }
    const unique = [...new Set(keys.values())];
    const [human, mouse] = await Promise.all([
      this.human.genesFor(unique),
      this.mouse.genesFor(unique),
    ]);
    for (const [guide, key] of keys) {
      const hit = hits.get(guide);
      if (!hit) continue;
      hit.humanIds = human.get(key) ?? [];
      hit.mouseIds = mouse.get(key) ?? [];
    }
  }
}
