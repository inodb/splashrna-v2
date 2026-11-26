/**
 * Sequence retrieval from NCBI E-utilities, straight from the browser.
 *
 * E-utilities send permissive CORS headers, so the static site can call them with no
 * backend of its own. NCBI asks unauthenticated clients to stay under three requests
 * per second and to identify themselves with a `tool` parameter; we send a generic
 * tool name and deliberately send no email address.
 *
 * If you run this for a whole lab and start hitting rate limits, add an API key via
 * `setApiKey()` rather than an email address.
 *
 * https://www.ncbi.nlm.nih.gov/books/NBK25497/
 */

const BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
const TOOL = "splashrna-v2";

let apiKey: string | null = null;

/** Optional NCBI API key, which raises the rate limit from 3/s to 10/s. */
export function setApiKey(key: string | null): void {
  apiKey = key;
}

function url(endpoint: string, params: Record<string, string>): string {
  const query = new URLSearchParams({ tool: TOOL, ...params });
  if (apiKey) query.set("api_key", apiKey);
  return `${BASE}/${endpoint}?${query.toString()}`;
}

// NCBI's guidance is at most three requests per second; serialise and space them out.
let queue: Promise<unknown> = Promise.resolve();
const MIN_INTERVAL_MS = 350;

function throttle<T>(task: () => Promise<T>): Promise<T> {
  const result = queue.then(async () => {
    const started = Date.now();
    try {
      return await task();
    } finally {
      const wait = MIN_INTERVAL_MS - (Date.now() - started);
      if (wait > 0) await new Promise((r) => setTimeout(r, wait));
    }
  });
  queue = result.catch(() => undefined);
  return result as Promise<T>;
}

async function request(endpoint: string, params: Record<string, string>): Promise<string> {
  return throttle(async () => {
    const response = await fetch(url(endpoint, params));
    if (!response.ok) {
      throw new Error(`NCBI returned HTTP ${response.status} for ${endpoint}`);
    }
    return response.text();
  });
}

export interface FetchedSequence {
  accession: string;
  description: string;
  sequence: string;
}

/** Parse a single-record FASTA payload. */
function parseSingleFasta(text: string): { header: string; sequence: string } {
  const lines = text.trim().split(/\r?\n/);
  if (lines.length === 0 || !lines[0].startsWith(">")) {
    throw new Error("NCBI did not return a FASTA record");
  }
  return {
    header: lines[0].slice(1).trim(),
    sequence: lines.slice(1).join("").replace(/\s/g, "").toUpperCase(),
  };
}

/** Fetch a transcript by RefSeq accession, e.g. NM_001127233.1. */
export async function fetchRefSeq(accession: string): Promise<FetchedSequence> {
  const text = await request("efetch.fcgi", {
    db: "nuccore",
    id: accession,
    rettype: "fasta",
    retmode: "text",
  });
  const { header, sequence } = parseSingleFasta(text);
  if (!sequence) throw new Error(`No sequence returned for ${accession}`);
  return { accession, description: header, sequence };
}

export interface GeneTranscripts {
  entrezId: string;
  symbol: string;
  description: string;
  accessions: string[];
}

/**
 * Resolve an Entrez gene ID to its RefSeq mRNA transcripts.
 *
 * NOTE ON SEMANTICS: the original SplashRNA server designed against the *intersection*
 * of a gene's isoforms, using a precomputed table. That table is not part of this
 * phase, so a gene ID resolves to its individual transcripts and the caller picks one.
 * Results for a multi-isoform gene may therefore differ from the 2017 server, and the
 * UI says so rather than hiding the difference.
 */
export async function fetchGeneTranscripts(entrezId: string): Promise<GeneTranscripts> {
  const summary = await request("esummary.fcgi", {
    db: "gene",
    id: entrezId,
    retmode: "json",
  });
  const parsed = JSON.parse(summary) as {
    result?: Record<string, { name?: string; description?: string }>;
  };
  const record = parsed.result?.[entrezId];
  if (!record) throw new Error(`Entrez gene ${entrezId} was not found`);

  const linkXml = await request("elink.fcgi", {
    dbfrom: "gene",
    db: "nuccore",
    id: entrezId,
    linkname: "gene_nuccore_refseqrna",
  });
  const ids = [...linkXml.matchAll(/<Id>(\d+)<\/Id>/g)].map((m) => m[1]).slice(1);

  let accessions: string[] = [];
  if (ids.length > 0) {
    const docs = await request("esummary.fcgi", {
      db: "nuccore",
      id: ids.join(","),
      retmode: "json",
    });
    const parsedDocs = JSON.parse(docs) as {
      result?: Record<string, { caption?: string; accessionversion?: string }>;
    };
    accessions = Object.values(parsedDocs.result ?? {})
      .map((d) => d.accessionversion ?? d.caption ?? "")
      .filter((a) => a.startsWith("NM_"));
  }

  return {
    entrezId,
    symbol: record.name ?? entrezId,
    description: record.description ?? "",
    accessions,
  };
}
