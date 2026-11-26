/** Parsing of the FASTA-style input the SplashRNA form accepts. */

export type EntryKind = "sequence" | "refseq" | "entrez" | "unknown";

export interface FastaEntry {
  label: string;
  value: string;
  kind: EntryKind;
}

const REFSEQ = /^(NM_[0-9.]+)$/i;
const ENTREZ = /^[0-9]+$/;
const SEQUENCE = /^[ACGTN*]+$/i;

export function classify(value: string): EntryKind {
  if (REFSEQ.test(value)) return "refseq";
  if (ENTREZ.test(value)) return "entrez";
  if (SEQUENCE.test(value)) return "sequence";
  return "unknown";
}

/**
 * Parse FASTA-ish input.
 *
 * Each record is a `>` header followed by one or more lines holding either a
 * nucleotide sequence, a RefSeq accession, or an Entrez gene ID. Lines beginning with
 * `#` are comments. A bare sequence with no header is accepted as a single record.
 */
export function parseFasta(text: string): FastaEntry[] {
  const entries: FastaEntry[] = [];
  let label: string | null = null;
  let buffer: string[] = [];

  const flush = () => {
    if (buffer.length === 0) return;
    const value = buffer.join("").trim();
    if (value) {
      entries.push({
        label: label ?? `Sequence ${entries.length + 1}`,
        value,
        kind: classify(value),
      });
    }
    buffer = [];
  };

  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith("#")) continue;
    if (line.startsWith(">")) {
      flush();
      label = line.slice(1).trim() || `Sequence ${entries.length + 1}`;
    } else {
      buffer.push(line.replace(/\s/g, ""));
    }
  }
  flush();
  return entries;
}

export const EXAMPLE_INPUT = `# Entries may be an Entrez gene ID, a RefSeq accession, or a sequence
> entrezID
71458
> refseqID
NM_001127233.1
> Sequence
ATGTAAATAACTTTAATAAAACGTTTATAGCAGTTACACAGAATTTCAATCCTAGTATATAGTACCTAGT
`;
