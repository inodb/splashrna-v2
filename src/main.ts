/** Browser entry point: wires the form to the engine. */

import { fetchGeneTranscripts, fetchRefSeq } from "./data/ncbi.js";
import {
  StaticAnnotationIndex,
  transcriptCoordinates,
} from "./data/annotations.js";
import {
  DEFAULT_TRANSCRIPT_OPTIONS,
  StaticGeneIndex,
  selectTranscript,
  type TranscriptOptions,
} from "./data/genes.js";
import { renderDiagram, type DiagramGuide } from "./ui/diagram.js";
import { StaticOffTargetIndex } from "./data/offtarget.js";
import type { OffTargetIndex } from "./engine/design.js";
import {
  DEFAULT_OPTIONS,
  type DesignOptions,
  type DesignResult,
  designGuides,
} from "./engine/design.js";
import { downloadCsv, toCsv } from "./engine/export.js";
import { EXAMPLE_INPUT, type FastaEntry, parseFasta } from "./engine/fasta.js";
import { loadModel } from "./engine/model.js";
import { SplashScorer } from "./engine/scorer.js";

const $ = <T extends HTMLElement>(id: string): T =>
  document.getElementById(id) as T;

const els = {
  fasta: $<HTMLTextAreaElement>("fasta"),
  count: $<HTMLInputElement>("count"),
  removeRE: $<HTMLInputElement>("removeRE"),
  inGeneDup: $<HTMLInputElement>("inGeneDup"),
  offTargetDup: $<HTMLInputElement>("offTargetDup"),
  offTargetMode: $<HTMLSelectElement>("offTargetMode"),
  polyACut: $<HTMLInputElement>("polyACut"),
  cdsOnly: $<HTMLInputElement>("cdsOnly"),
  run: $<HTMLButtonElement>("run"),
  csv: $<HTMLButtonElement>("csv"),
  example: $<HTMLButtonElement>("example"),
  clear: $<HTMLButtonElement>("clear"),
  status: $<HTMLDivElement>("status"),
  progress: $<HTMLProgressElement>("progress"),
  results: $<HTMLElement>("results"),
  offTargetHint: $<HTMLElement>("offTargetHint"),
  academic: $<HTMLInputElement>("academic"),
  email: $<HTMLInputElement>("email"),
};

/**
 * The academic-use gate, as the original server had it.
 *
 * That server blocked results until the box was ticked and an address matching this
 * same pattern was supplied, then recorded the address. Neither side verified
 * anything, so this was always an attestation rather than enforcement — and a static
 * site has nowhere to record it, so only the attestation survives. Whether it is
 * required at all is a licensing question, not a technical one.
 */
const EMAIL_PATTERN = /[^@]+@[^@]+\.[^@]+/;
const STORAGE_KEY = "splashrna.academic";

function gateSatisfied(): boolean {
  return els.academic.checked && EMAIL_PATTERN.test(els.email.value.trim());
}

/** Remember the attestation locally so it is not retyped on every visit. */
function rememberGate(): void {
  try {
    localStorage.setItem(
      STORAGE_KEY,
      JSON.stringify({ academic: els.academic.checked, email: els.email.value.trim() }),
    );
  } catch {
    // Private browsing, or storage disabled. Nothing here is important enough to warn about.
  }
}

function restoreGate(): void {
  try {
    const saved = localStorage.getItem(STORAGE_KEY);
    if (!saved) return;
    const { academic, email } = JSON.parse(saved) as { academic: boolean; email: string };
    els.academic.checked = Boolean(academic);
    els.email.value = typeof email === "string" ? email : "";
  } catch {
    // Ignore anything unparseable.
  }
}

const plural = (n: number, singular: string, plural = `${singular}s`) =>
  `${n.toLocaleString()} ${n === 1 ? singular : plural}`;

let scorer: SplashScorer | null = null;
let loading: Promise<SplashScorer> | null = null;
let lastResults: DesignResult[] = [];

/**
 * The off-target index is an optional deployment. `undefined` means we have not looked
 * yet; `null` means it is genuinely not deployed, and we stop retrying.
 */
let offTarget: OffTargetIndex | null | undefined;

async function getOffTargetIndex(): Promise<OffTargetIndex | null> {
  if (offTarget !== undefined) return offTarget;
  try {
    offTarget = await StaticOffTargetIndex.load("./offtarget");
  } catch {
    offTarget = null;
  }
  return offTarget;
}

/**
 * Annotations for the transcript diagram. Loaded on first use only — a session that
 * never designs against a gene ID never pays for it.
 */
let annotations: StaticAnnotationIndex | null | undefined;

async function getAnnotationIndex(): Promise<StaticAnnotationIndex | null> {
  if (annotations !== undefined) return annotations;
  try {
    annotations = await StaticAnnotationIndex.load("./annotations");
  } catch {
    annotations = null;
  }
  return annotations;
}

/** Diagrams built during a run, keyed by the result label they belong to. */
const diagrams = new Map<string, SVGSVGElement>();

/** Likewise optional: without it, Entrez IDs fall back to a single RefSeq transcript. */
let geneIndex: StaticGeneIndex | null | undefined;

async function getGeneIndex(): Promise<StaticGeneIndex | null> {
  if (geneIndex !== undefined) return geneIndex;
  try {
    geneIndex = await StaticGeneIndex.load("./genes");
  } catch {
    geneIndex = null;
  }
  return geneIndex;
}

function setStatus(message: string, isError = false): void {
  els.status.textContent = message;
  els.status.classList.toggle("error", isError);
}

/** Load the model once, on first use, so the page paints immediately. */
function getScorer(): Promise<SplashScorer> {
  if (scorer) return Promise.resolve(scorer);
  loading ??= loadModel("./splashrna-model.json").then((model) => {
    scorer = new SplashScorer(model);
    return scorer;
  });
  return loading;
}

function readTranscriptOptions(): TranscriptOptions {
  return {
    ...DEFAULT_TRANSCRIPT_OPTIONS,
    polyACut: els.polyACut.checked,
    cdsOnly: els.cdsOnly.checked,
  };
}

function readOptions(): DesignOptions {
  const count = Number.parseInt(els.count.value, 10);
  return {
    ...DEFAULT_OPTIONS,
    predictionsPerGene: Number.isFinite(count) && count > 0 ? count : 5,
    removeRestrictionSites: els.removeRE.checked,
    allowInGeneDuplicates: els.inGeneDup.checked,
    allowOffTargetDuplicates: els.offTargetDup.checked,
    offTargetMode: els.offTargetMode.value as DesignOptions["offTargetMode"],
  };
}

/** Turn one parsed record into a label and a nucleotide sequence. */
interface ResolvedEntry {
  label: string;
  sequence: string;
  note?: string;
  /** Set when the sequence came from the gene index, so a diagram can be drawn. */
  entrezId?: number;
  /** Guide offsets only map to genomic coordinates for the intersection itself. */
  mappable?: boolean;
}

async function resolveEntry(entry: FastaEntry): Promise<ResolvedEntry> {
  switch (entry.kind) {
    case "sequence":
      return { label: entry.label, sequence: entry.value.replace(/\*/g, "") };

    case "refseq": {
      setStatus(`Fetching ${entry.value} from NCBI…`);
      const record = await fetchRefSeq(entry.value);
      return { label: `${entry.label} (${record.accession})`, sequence: record.sequence };
    }

    case "entrez": {
      // Preferred path: the isoform intersection, as the original server used.
      const genes = await getGeneIndex();
      if (genes) {
        setStatus(`Looking up gene ${entry.value}…`);
        const intersection = await genes.get(Number.parseInt(entry.value, 10));
        if (intersection) {
          // polyA trimming and CDS-only apply to gene IDs, where the intersection and
          // CDS tables exist. A pasted sequence or a RefSeq accession has neither.
          const choice = selectTranscript(intersection, readTranscriptOptions());
          return {
            label: `${entry.label} (Entrez ${entry.value})`,
            sequence: choice.sequence,
            note: choice.warnings[0],
            entrezId: intersection.entrezId,
            // A polyA-trimmed transcript is a prefix of the intersection, so offsets
            // still line up. The CDS is a different sequence with different
            // coordinates, which this index does not carry.
            mappable: choice.source !== "cds",
          };
        }
      }

      // Fallback: a single RefSeq transcript, which is not the same thing.
      setStatus(`Looking up gene ${entry.value} at NCBI…`);
      const gene = await fetchGeneTranscripts(entry.value);
      if (gene.accessions.length === 0) {
        throw new Error(`No RefSeq mRNA transcripts found for Entrez gene ${entry.value}`);
      }
      const accession = gene.accessions[0];
      const record = await fetchRefSeq(accession);
      return {
        label: `${entry.label} — ${gene.symbol} (${accession})`,
        sequence: record.sequence,
        note: `Gene data is not available here, so ${accession} was used${
          gene.accessions.length > 1
            ? ` — note ${gene.symbol} has ${gene.accessions.length} RefSeq transcripts`
            : ""
        }. Designing against a single transcript can give different guides than the isoform intersection.`,
      };
    }

    default:
      throw new Error(
        `"${entry.value.slice(0, 40)}" is not a nucleotide sequence, a RefSeq accession, or an Entrez gene ID`,
      );
  }
}

function renderResults(results: DesignResult[]): void {
  els.results.replaceChildren();

  for (const result of results) {
    const section = document.createElement("article");
    section.className = "result";

    const header = document.createElement("header");
    const title = document.createElement("h3");
    title.textContent = result.label;
    const meta = document.createElement("span");
    meta.className = "meta";
    const dropped =
      result.discarded.ambiguousBases +
      result.discarded.restrictionSites +
      result.discarded.inGeneDuplicates +
      result.discarded.offTarget;
    meta.textContent =
      `${result.sequenceLength.toLocaleString()} nt · ` +
      `${result.guides.length} shown · ` +
      `${plural(dropped, "candidate")} filtered`;
    header.append(title, meta);
    section.append(header);

    if (result.warnings.length > 0) {
      const box = document.createElement("div");
      box.className = "warnings";
      const list = document.createElement("ul");
      for (const warning of result.warnings) {
        const item = document.createElement("li");
        item.textContent = warning;
        list.append(item);
      }
      box.append(list);
      section.append(box);
    }

    if (result.guides.length > 0) {
      const scroll = document.createElement("div");
      scroll.className = "scroll";
      const table = document.createElement("table");
      table.innerHTML =
        "<thead><tr><th>Position</th><th>Antisense guide</th><th>SplashRNA</th>" +
        "<th>miR-30</th><th>miR-E</th><th>Off-target</th>" +
        "<th>97-mer construct</th><th>Notes</th></tr></thead>";
      const body = document.createElement("tbody");

      for (const guide of result.guides) {
        const row = document.createElement("tr");
        const cells: HTMLTableCellElement[] = [];

        const add = (text: string, className?: string) => {
          const cell = document.createElement("td");
          cell.textContent = text;
          if (className) cell.className = className;
          cells.push(cell);
          return cell;
        };

        add(String(guide.position));
        add(guide.guide, "seq");
        add(guide.splash.toFixed(3), guide.splash >= 1 ? "score high" : "score");
        add(guide.mir30.toFixed(3), "score");
        add(guide.mirE.toFixed(3), "score");

        // "—" distinguishes "index not consulted" from a genuine zero.
        if (guide.offTarget === null) {
          add("—");
        } else {
          const { humanGenes, mouseGenes, humanIds, mouseIds } = guide.offTarget;
          if (humanGenes + mouseGenes === 0) {
            add("unique");
          } else {
            const cell = document.createElement("td");
            cell.append(`${humanGenes} human / ${mouseGenes} mouse`);
            const ids = [...(humanIds ?? []), ...(mouseIds ?? [])];
            if (ids.length > 0) {
              const list = document.createElement("div");
              list.className = "offtarget-genes";
              // Cap the visible list; a 22-mer can hit over a hundred genes.
              ids.slice(0, 6).forEach((id, i) => {
                if (i > 0) list.append(" ");
                const a = document.createElement("a");
                a.href = `https://www.ncbi.nlm.nih.gov/gene/${id}`;
                a.target = "_blank";
                a.rel = "noopener noreferrer";
                a.textContent = String(id);
                list.append(a);
              });
              if (ids.length > 6) list.append(` +${ids.length - 6}`);
              cell.append(list);
            }
            cells.push(cell);
          }
        }

        const construct = document.createElement("td");
        const code = document.createElement("code");
        code.className = "construct";
        code.textContent = guide.construct;
        code.title = guide.construct;
        construct.append(code);
        cells.push(construct);

        add(guide.warnings.join(" "));
        row.append(...cells);
        body.append(row);
      }

      table.append(body);
      scroll.append(table);
      section.append(scroll);

      const diagram = diagrams.get(result.label);
      if (diagram) {
        const figure = document.createElement("figure");
        figure.className = "diagram-figure";
        figure.append(diagram);
        const caption = document.createElement("figcaption");
        caption.innerHTML =
          '<span class="key key-cds"></span> coding ' +
          '<span class="key key-utr"></span> UTR ' +
          '<span class="key key-int"></span> designed region ' +
          '<span class="key key-guide"></span> shRNA';
        figure.append(caption);
        section.append(figure);
      }
    }

    els.results.append(section);
  }
}

async function run(): Promise<void> {
  if (!gateSatisfied()) {
    setStatus(
      "Please confirm that you are an academic user by checking the box and entering your email.",
      true,
    );
    (els.academic.checked ? els.email : els.academic).focus();
    return;
  }
  rememberGate();

  const entries = parseFasta(els.fasta.value);
  if (entries.length === 0) {
    setStatus("Add a target first — a sequence, a gene ID, or a RefSeq accession.", true);
    return;
  }

  els.run.disabled = true;
  els.csv.disabled = true;
  els.results.replaceChildren();
  diagrams.clear();
  els.progress.hidden = false;
  els.progress.value = 0;

  try {
    setStatus("Loading the scoring model…");
    const engine = await getScorer();
    const options = readOptions();

    const index = await getOffTargetIndex();
    const missingIndex = index === null;

    const results: DesignResult[] = [];

    for (const [idx, entry] of entries.entries()) {
      const resolved = await resolveEntry(entry);
      const { label, sequence, note } = resolved;
      setStatus(`Scoring ${label} (${idx + 1} of ${entries.length})…`);

      const result = await designGuides(
        { label, sequence },
        engine,
        options,
        index ?? undefined,
        (fraction) => {
          els.progress.value = (idx + fraction) / entries.length;
        },
      );
      if (note) result.warnings.unshift(note);

      // Draw the transcript diagram when we know which gene this was.
      if (resolved.entrezId !== undefined) {
        const index = await getAnnotationIndex();
        const annotation = index ? await index.get(resolved.entrezId) : null;
        if (annotation) {
          const coordinates = resolved.mappable
            ? transcriptCoordinates(annotation.intersection, annotation.strand)
            : [];
          const marks: DiagramGuide[] = result.guides
            .map((g) => ({
              position: g.position,
              genomic: coordinates[g.position] ?? NaN,
              guide: g.guide,
              splash: g.splash,
            }))
            .filter((g) => Number.isFinite(g.genomic));
          const svg = renderDiagram(annotation, marks);
          if (svg) diagrams.set(label, svg);
        }
      }
      if (missingIndex) {
        result.warnings.push(
          "Off-target filtering was skipped — that data is not available here, so these guides may also hit other human or mouse genes.",
        );
      }
      results.push(result);
      renderResults(results);
    }

    lastResults = results;
    els.csv.disabled = results.length === 0;
    const total = results.reduce((n, r) => n + r.guides.length, 0);
    setStatus(`Done — ${total} shRNA${total === 1 ? "" : "s"} across ${results.length} input${results.length === 1 ? "" : "s"}.`);
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  } finally {
    els.run.disabled = false;
    els.progress.hidden = true;
  }
}

els.run.addEventListener("click", () => void run());
els.example.addEventListener("click", () => {
  els.fasta.value = EXAMPLE_INPUT;
  setStatus("");
});
els.clear.addEventListener("click", () => {
  els.fasta.value = "";
  els.results.replaceChildren();
  els.csv.disabled = true;
  setStatus("");
});
els.csv.addEventListener("click", () => {
  if (lastResults.length > 0) downloadCsv(toCsv(lastResults));
});

els.fasta.value = EXAMPLE_INPUT;
restoreGate();

// Say up front whether off-target filtering is actually available on this deployment,
// rather than claiming either way and being wrong.
void getOffTargetIndex().then((index) => {
  els.offTargetHint.textContent = index
    ? "Otherwise a guide is dropped if its 22-mer matches another human or mouse gene exactly."
    : "Off-target data is not available here, so this cannot be checked.";
});
