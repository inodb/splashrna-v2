/**
 * Static gene index: Entrez gene ID -> isoform-intersection transcript.
 *
 * The original server designed against the intersection of a gene's RefSeq isoforms --
 * the stretch common to all of them -- rather than any single transcript. This restores
 * that, from static files built by tools/build_genes.py.
 *
 * Only the 636 KB index is fetched up front. The 58 MB sequence blob is never
 * downloaded whole: a gene is read with a single Range request for its few kilobytes.
 */

const SYMBOLS = "ACGT*N";

interface Manifest {
  formatVersion: number;
  count: number;
  blobBytes: number;
  indexBytes: number;
  blob: string;
  index: string;
}

export interface GeneTranscript {
  entrezId: number;
  strand: "+" | "-";
  /** Intersection sequence; `*` marks a break between non-contiguous blocks. */
  sequence: string;
}

/** Record layout in genes.idx: id, offset, length, strand, 3 bytes padding. */
const RECORD_BYTES = 16;

export class StaticGeneIndex {
  private readonly ids: Uint32Array;
  private readonly offsets: Uint32Array;
  private readonly lengths: Uint32Array;
  private readonly strands: Uint8Array;
  private readonly cache = new Map<number, GeneTranscript>();

  private constructor(
    private readonly baseUrl: string,
    private readonly manifest: Manifest,
    index: ArrayBuffer,
  ) {
    const count = index.byteLength / RECORD_BYTES;
    const view = new DataView(index);
    this.ids = new Uint32Array(count);
    this.offsets = new Uint32Array(count);
    this.lengths = new Uint32Array(count);
    this.strands = new Uint8Array(count);

    for (let i = 0; i < count; i++) {
      const base = i * RECORD_BYTES;
      this.ids[i] = view.getUint32(base, true);
      this.offsets[i] = view.getUint32(base + 4, true);
      this.lengths[i] = view.getUint32(base + 8, true);
      this.strands[i] = view.getUint8(base + 12);
    }
  }

  static async load(baseUrl = "./genes"): Promise<StaticGeneIndex> {
    const response = await fetch(`${baseUrl}/manifest.json`);
    if (!response.ok) {
      throw new Error(`No gene index at ${baseUrl} (HTTP ${response.status})`);
    }
    const manifest = (await response.json()) as Manifest;
    if (manifest.formatVersion !== 1) {
      throw new Error(`Unsupported gene index version ${manifest.formatVersion}`);
    }
    const index = await fetch(`${baseUrl}/${manifest.index}`);
    if (!index.ok) {
      throw new Error(`Could not load the gene index (HTTP ${index.status})`);
    }
    return new StaticGeneIndex(baseUrl, manifest, await index.arrayBuffer());
  }

  get count(): number {
    return this.ids.length;
  }

  /** Index of `entrezId` in the sorted id array, or -1. */
  private find(entrezId: number): number {
    let low = 0;
    let high = this.ids.length - 1;
    while (low <= high) {
      const mid = (low + high) >> 1;
      const value = this.ids[mid];
      if (value === entrezId) return mid;
      if (value < entrezId) low = mid + 1;
      else high = mid - 1;
    }
    return -1;
  }

  has(entrezId: number): boolean {
    return this.find(entrezId) >= 0;
  }

  /** Fetch one gene's intersection transcript, or null if the gene is not in the table. */
  async get(entrezId: number): Promise<GeneTranscript | null> {
    const cached = this.cache.get(entrezId);
    if (cached) return cached;

    const slot = this.find(entrezId);
    if (slot < 0) return null;

    const symbols = this.lengths[slot];
    const start = this.offsets[slot];
    const end = start + Math.ceil(symbols / 2);

    const response = await fetch(`${this.baseUrl}/${this.manifest.blob}`, {
      headers: { Range: `bytes=${start}-${end - 1}` },
    });
    if (!response.ok) {
      throw new Error(`Could not read gene ${entrezId} (HTTP ${response.status})`);
    }

    let bytes = new Uint8Array(await response.arrayBuffer());
    // A host that ignores Range hands back the whole blob; take our slice out of it.
    if (response.status === 200 && bytes.length === this.manifest.blobBytes) {
      bytes = bytes.subarray(start, end);
    }

    let sequence = "";
    for (let i = 0; i < symbols; i++) {
      const byte = bytes[i >> 1];
      sequence += SYMBOLS[i % 2 === 0 ? byte >> 4 : byte & 0x0f];
    }

    const transcript: GeneTranscript = {
      entrezId,
      strand: this.strands[slot] === 1 ? "-" : "+",
      sequence,
    };
    this.cache.set(entrezId, transcript);
    return transcript;
  }
}
