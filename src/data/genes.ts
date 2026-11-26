/**
 * Static gene index: Entrez gene ID -> isoform-intersection transcript, CDS
 * intersection, and polyadenylation cut point.
 *
 * The original server designed against the intersection of a gene's RefSeq isoforms --
 * the stretch common to all of them -- rather than any single transcript, and trimmed
 * it at the annotated polyA site before designing. This restores both, from static
 * files built by tools/build_genes.py.
 *
 * Only the 1.27 MB index is fetched up front. The 58 MB transcript blob and 32 MB CDS
 * blob are never downloaded whole: a gene is read with a single Range request for its
 * few kilobytes.
 */

const SYMBOLS = "ACGT*N";

interface Manifest {
  formatVersion: number;
  count: number;
  blobBytes: number;
  cdsBytes: number;
  recordBytes: number;
  blob: string;
  cds: string;
  index: string;
}

export interface GeneTranscript {
  entrezId: number;
  strand: "+" | "-";
  /** Intersection sequence; `*` marks a break between non-contiguous blocks. */
  sequence: string;
  /** CDS-only intersection, empty when the gene has none. */
  cds: string;
  /** Prefix length to keep when trimming at the polyA site, or -1 if unresolved. */
  polyACutIndex: number;
}

/** Which transcript a design run should actually use. */
export interface TranscriptChoice {
  sequence: string;
  warnings: string[];
  source: "intersection" | "polyA-trimmed" | "cds";
}

export interface TranscriptOptions {
  /** Trim at the annotated polyA site. On by default, as in the original server. */
  polyACut: boolean;
  /** Design against the CDS intersection only. */
  cdsOnly: boolean;
}

export const DEFAULT_TRANSCRIPT_OPTIONS: TranscriptOptions = {
  polyACut: true,
  cdsOnly: false,
};

/**
 * Apply the transcript-selection options, reproducing `polyaCut()` from the original
 * including its fallbacks and their wording.
 */
export function selectTranscript(
  gene: GeneTranscript,
  options: Partial<TranscriptOptions> = {},
): TranscriptChoice {
  const { polyACut, cdsOnly } = { ...DEFAULT_TRANSCRIPT_OPTIONS, ...options };

  if (cdsOnly) {
    if (!gene.cds) {
      return {
        sequence: "",
        warnings: [
          `Gene ${gene.entrezId} has no CDS intersection. Clearing "CDS only" may help.`,
        ],
        source: "cds",
      };
    }
    return { sequence: gene.cds, warnings: [], source: "cds" };
  }

  if (!polyACut) {
    return { sequence: gene.sequence, warnings: [], source: "intersection" };
  }

  if (gene.polyACutIndex >= 0) {
    return {
      sequence: gene.sequence.slice(0, gene.polyACutIndex),
      warnings: [],
      source: "polyA-trimmed",
    };
  }

  // No usable polyA annotation: the original substituted the CDS intersection, or
  // fell back to the untrimmed transcript when there was no CDS either.
  if (gene.cds) {
    return {
      sequence: gene.cds,
      warnings: [
        "Could not find annotation for polyadenylation site. CDS intersection was substituted.",
      ],
      source: "cds",
    };
  }
  return {
    sequence: gene.sequence,
    warnings: [
      "Could not find annotation for polyadenylation site and there is no intersection for CDS regions. Predictions on full transcript are shown.",
    ],
    source: "intersection",
  };
}

export class StaticGeneIndex {
  private readonly ids: Uint32Array;
  private readonly seqOffsets: Uint32Array;
  private readonly seqLengths: Uint32Array;
  private readonly cuts: Int32Array;
  private readonly cdsOffsets: Uint32Array;
  private readonly cdsLengths: Uint32Array;
  private readonly strands: Uint8Array;
  private readonly cache = new Map<number, GeneTranscript>();

  private constructor(
    private readonly baseUrl: string,
    private readonly manifest: Manifest,
    index: ArrayBuffer,
  ) {
    const size = manifest.recordBytes;
    const count = index.byteLength / size;
    const view = new DataView(index);

    this.ids = new Uint32Array(count);
    this.seqOffsets = new Uint32Array(count);
    this.seqLengths = new Uint32Array(count);
    this.cuts = new Int32Array(count);
    this.cdsOffsets = new Uint32Array(count);
    this.cdsLengths = new Uint32Array(count);
    this.strands = new Uint8Array(count);

    for (let i = 0; i < count; i++) {
      const at = i * size;
      this.ids[i] = view.getUint32(at, true);
      this.seqOffsets[i] = view.getUint32(at + 4, true);
      this.seqLengths[i] = view.getUint32(at + 8, true);
      this.cuts[i] = view.getInt32(at + 12, true);
      this.cdsOffsets[i] = view.getUint32(at + 16, true);
      this.cdsLengths[i] = view.getUint32(at + 20, true);
      this.strands[i] = view.getUint8(at + 24);
    }
  }

  static async load(baseUrl = "./genes"): Promise<StaticGeneIndex> {
    const response = await fetch(`${baseUrl}/manifest.json`);
    if (!response.ok) {
      throw new Error(`No gene index at ${baseUrl} (HTTP ${response.status})`);
    }
    const manifest = (await response.json()) as Manifest;
    if (manifest.formatVersion !== 2) {
      throw new Error(
        `Unsupported gene index version ${manifest.formatVersion}; rebuild with tools/build_genes.py`,
      );
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

  /** Read one slice of a blob with a Range request and unpack it. */
  private async readSequence(
    file: string,
    totalBytes: number,
    offset: number,
    symbols: number,
  ): Promise<string> {
    if (symbols === 0) return "";
    const end = offset + Math.ceil(symbols / 2);

    const response = await fetch(`${this.baseUrl}/${file}`, {
      headers: { Range: `bytes=${offset}-${end - 1}` },
    });
    if (!response.ok) {
      throw new Error(`Could not read ${file} (HTTP ${response.status})`);
    }
    let bytes = new Uint8Array(await response.arrayBuffer());
    // A host that ignores Range hands back the whole blob; take our slice out of it.
    if (response.status === 200 && bytes.length === totalBytes) {
      bytes = bytes.subarray(offset, end);
    }

    let out = "";
    for (let i = 0; i < symbols; i++) {
      const byte = bytes[i >> 1];
      out += SYMBOLS[i % 2 === 0 ? byte >> 4 : byte & 0x0f];
    }
    return out;
  }

  /** Fetch one gene, or null if it is not in the table. */
  async get(entrezId: number): Promise<GeneTranscript | null> {
    const cached = this.cache.get(entrezId);
    if (cached) return cached;

    const slot = this.find(entrezId);
    if (slot < 0) return null;

    const [sequence, cds] = await Promise.all([
      this.readSequence(
        this.manifest.blob,
        this.manifest.blobBytes,
        this.seqOffsets[slot],
        this.seqLengths[slot],
      ),
      this.readSequence(
        this.manifest.cds,
        this.manifest.cdsBytes,
        this.cdsOffsets[slot],
        this.cdsLengths[slot],
      ),
    ]);

    const transcript: GeneTranscript = {
      entrezId,
      strand: this.strands[slot] === 1 ? "-" : "+",
      sequence,
      cds,
      polyACutIndex: this.cuts[slot],
    };
    this.cache.set(entrezId, transcript);
    return transcript;
  }
}
