/** Browser entry point: wires the form to the engine. */

import { fetchGeneTranscripts, fetchRefSeq } from "./data/ncbi.js";
import {
  DEFAULT_TRANSCRIPT_OPTIONS,
  StaticGeneIndex,
  selectTranscript,
  type TranscriptOptions,
} from "./data/genes.js";
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
};

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
async function resolveEntry(
  entry: FastaEntry,
): Promise<{ label: string; sequence: string; note?: string }> {
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
        setStatus(`Loading the isoform intersection for Entrez gene ${entry.value}…`);
        const intersection = await genes.get(Number.parseInt(entry.value, 10));
        if (intersection) {
          // polyA trimming and CDS-only apply to gene IDs, where the intersection and
          // CDS tables exist. A pasted sequence or a RefSeq accession has neither.
          const choice = selectTranscript(intersection, readTranscriptOptions());
          return {
            label: `${entry.label} (Entrez ${entry.value})`,
            sequence: choice.sequence,
            note: choice.warnings[0],
          };
        }
      }

      // Fallback: a single RefSeq transcript, which is not the same thing.
      setStatus(`Looking up Entrez gene ${entry.value}…`);
      const gene = await fetchGeneTranscripts(entry.value);
      if (gene.accessions.length === 0) {
        throw new Error(`No RefSeq mRNA transcripts found for Entrez gene ${entry.value}`);
      }
      const accession = gene.accessions[0];
      const record = await fetchRefSeq(accession);
      return {
        label: `${entry.label} — ${gene.symbol} (${accession})`,
        sequence: record.sequence,
        note: `The isoform-intersection index is not deployed on this host, so ${accession} was used${
          gene.accessions.length > 1
            ? ` (${gene.symbol} has ${gene.accessions.length} RefSeq transcripts)`
            : ""
        }. The original server designed against the intersection of all isoforms, so results may differ.`,
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
          const { humanGenes, mouseGenes } = guide.offTarget;
          add(
            humanGenes + mouseGenes === 0
              ? "unique"
              : `${humanGenes} human / ${mouseGenes} mouse`,
          );
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
    }

    els.results.append(section);
  }
}

async function run(): Promise<void> {
  const entries = parseFasta(els.fasta.value);
  if (entries.length === 0) {
    setStatus("Enter at least one sequence, accession, or gene ID.", true);
    return;
  }

  els.run.disabled = true;
  els.csv.disabled = true;
  els.results.replaceChildren();
  els.progress.hidden = false;
  els.progress.value = 0;

  try {
    setStatus("Loading the SplashRNA model…");
    const engine = await getScorer();
    const options = readOptions();

    setStatus("Checking for the off-target index…");
    const index = await getOffTargetIndex();
    const missingIndex = index === null;

    const results: DesignResult[] = [];

    for (const [idx, entry] of entries.entries()) {
      const { label, sequence, note } = await resolveEntry(entry);
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
      if (missingIndex) {
        result.warnings.push(
          "Off-target filtering was skipped: the index is not deployed on this host. Guides may match other human or mouse transcripts.",
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

// Say up front whether off-target filtering is actually available on this deployment,
// rather than claiming either way and being wrong.
void getOffTargetIndex().then((index) => {
  els.offTargetHint.textContent = index
    ? "Guides with an exact match in another human or mouse gene."
    : "Unavailable: the off-target index is not deployed on this host.";
});
