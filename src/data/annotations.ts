/**
 * Isoform annotations behind the transcript diagram.
 *
 * These are the 2015 annotations the predictions were computed against, stored in the
 * legacy encoding and parsed here, so the picture agrees with the numbers beside it.
 * See tools/build_annotations.py for why this ships its own data rather than embedding
 * a genome browser pointed at current RefSeq.
 *
 * Only the 0.8 MB index is loaded, and only when a diagram is first opened. Each gene
 * is about 700 bytes, fetched with one Range request.
 */

export type ExonKind = "cds" | "utr5" | "utr3" | "other";

export interface Exon {
  kind: ExonKind;
  start: number;
  end: number;
}

export interface Isoform {
  refseq: string;
  exons: Exon[];
}

export interface GeneAnnotation {
  entrezId: number;
  chrom: string;
  strand: "+" | "-";
  isoforms: Isoform[];
  /** Genomic ranges of the isoform intersection — the region actually designed against. */
  intersection: { start: number; end: number }[];
}

interface Manifest {
  formatVersion: number;
  blobBytes: number;
  recordBytes: number;
  blob: string;
  index: string;
}

const RECORD_FIELDS = 5;

function parseKind(value: string): ExonKind {
  return value === "cds" || value === "utr5" || value === "utr3" ? value : "other";
}

/**
 * Parse `chrom,strand:REFSEQ,annot,start,end;annot,start,end:REFSEQ2,...`
 *
 * The first `;`-separated piece of each isoform carries the accession, the rest do not
 * — an artefact of how the original server concatenated these strings.
 */
export function parseIsoforms(text: string): {
  chrom: string;
  strand: "+" | "-";
  isoforms: Isoform[];
} {
  const parts = text.split(":");
  const [chrom = "", strand = "+"] = (parts[0] ?? "").split(",");
  const isoforms: Isoform[] = [];

  for (const block of parts.slice(1)) {
    if (!block) continue;
    const pieces = block.split(";");
    const first = pieces[0].split(",");
    if (first.length < 4) continue;

    const refseq = first[0];
    const exons: Exon[] = [
      { kind: parseKind(first[1]), start: Number(first[2]), end: Number(first[3]) },
    ];
    for (const piece of pieces.slice(1)) {
      const f = piece.split(",");
      if (f.length < 3) continue;
      exons.push({ kind: parseKind(f[0]), start: Number(f[1]), end: Number(f[2]) });
    }
    exons.sort((a, b) => a.start - b.start);
    isoforms.push({ refseq, exons });
  }

  return { chrom, strand: strand === "-" ? "-" : "+", isoforms };
}

/**
 * Parse `start,end*;start,end` into ascending genomic ranges.
 *
 * Ends carry one *or more* trailing asterisks — `12048303,12048545**` occurs in the
 * real table. Stripping only the first leaves a non-numeric string, which silently
 * drops that whole block from the transcript and shifts every downstream coordinate.
 * The original used `rstrip("*")`, so strip all of them.
 */
export function parseRanges(text: string): { start: number; end: number }[] {
  const out: { start: number; end: number }[] = [];
  for (const item of text.replace(/;$/, "").split(";")) {
    if (!item) continue;
    const [a, b] = item.replace(/:/g, ",").split(",");
    if (b === undefined || b === ".") continue;
    const start = Number(a);
    const end = Number(b.replace(/\*+$/, ""));
    if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) continue;
    out.push({ start, end });
  }
  out.sort((a, b) => a.start - b.start);
  return out;
}

/**
 * Every genomic coordinate of the intersection, in transcript order.
 *
 * Ascending along the genome for a + strand gene; reversed for a - strand gene, since
 * the transcript runs the other way. This mirrors `get_position()` followed by the
 * reverse in `entrezID_prepare()` in the original server.
 */
export function transcriptCoordinates(
  intersection: { start: number; end: number }[],
  strand: "+" | "-",
): number[] {
  const out: number[] = [];
  for (const { start, end } of intersection) {
    for (let p = start; p <= end; p++) out.push(p);
  }
  if (strand === "-") out.reverse();
  return out;
}

export class StaticAnnotationIndex {
  private readonly ids: Uint32Array;
  private readonly offsets: Uint32Array; // [isoOffset, isoLen, rangeOffset, rangeLen] per gene
  private readonly cache = new Map<number, GeneAnnotation>();

  private constructor(
    private readonly baseUrl: string,
    private readonly manifest: Manifest,
    index: ArrayBuffer,
  ) {
    const size = manifest.recordBytes;
    const count = index.byteLength / size;
    const view = new DataView(index);
    this.ids = new Uint32Array(count);
    this.offsets = new Uint32Array(count * 4);
    for (let i = 0; i < count; i++) {
      const at = i * size;
      this.ids[i] = view.getUint32(at, true);
      for (let f = 0; f < RECORD_FIELDS - 1; f++) {
        this.offsets[i * 4 + f] = view.getUint32(at + 4 + f * 4, true);
      }
    }
  }

  static async load(baseUrl = "./annotations"): Promise<StaticAnnotationIndex> {
    const response = await fetch(`${baseUrl}/manifest.json`);
    if (!response.ok) {
      throw new Error(`No annotation index at ${baseUrl} (HTTP ${response.status})`);
    }
    const manifest = (await response.json()) as Manifest;
    if (manifest.formatVersion !== 1) {
      throw new Error(`Unsupported annotation index version ${manifest.formatVersion}`);
    }
    const index = await fetch(`${baseUrl}/${manifest.index}`);
    if (!index.ok) {
      throw new Error(`Could not load the annotation index (HTTP ${index.status})`);
    }
    return new StaticAnnotationIndex(baseUrl, manifest, await index.arrayBuffer());
  }

  private find(entrezId: number): number {
    let low = 0;
    let high = this.ids.length - 1;
    while (low <= high) {
      const mid = (low + high) >> 1;
      if (this.ids[mid] === entrezId) return mid;
      if (this.ids[mid] < entrezId) low = mid + 1;
      else high = mid - 1;
    }
    return -1;
  }

  private async readRange(offset: number, length: number): Promise<string> {
    if (length === 0) return "";
    const response = await fetch(`${this.baseUrl}/${this.manifest.blob}`, {
      headers: { Range: `bytes=${offset}-${offset + length - 1}` },
    });
    if (!response.ok) {
      throw new Error(`Could not read annotations (HTTP ${response.status})`);
    }
    if (response.status === 200) {
      // Host ignored Range and sent the whole blob.
      const whole = await response.text();
      return whole.slice(offset, offset + length);
    }
    return response.text();
  }

  async get(entrezId: number): Promise<GeneAnnotation | null> {
    const cached = this.cache.get(entrezId);
    if (cached) return cached;

    const slot = this.find(entrezId);
    if (slot < 0) return null;

    const base = slot * 4;
    const [isoText, rangeText] = await Promise.all([
      this.readRange(this.offsets[base], this.offsets[base + 1]),
      this.readRange(this.offsets[base + 2], this.offsets[base + 3]),
    ]);

    const { chrom, strand, isoforms } = parseIsoforms(isoText);
    const annotation: GeneAnnotation = {
      entrezId,
      chrom,
      strand,
      isoforms,
      intersection: parseRanges(rangeText),
    };
    this.cache.set(entrezId, annotation);
    return annotation;
  }
}
