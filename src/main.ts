/** Browser entry point: wires the form to the engine. */

import { fetchGeneTranscripts, fetchRefSeq } from "./data/ncbi.js";
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
  run: $<HTMLButtonElement>("run"),
  csv: $<HTMLButtonElement>("csv"),
  example: $<HTMLButtonElement>("example"),
  clear: $<HTMLButtonElement>("clear"),
  status: $<HTMLDivElement>("status"),
  progress: $<HTMLProgressElement>("progress"),
  results: $<HTMLElement>("results"),
};

let scorer: SplashScorer | null = null;
let loading: Promise<SplashScorer> | null = null;
let lastResults: DesignResult[] = [];

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
        note:
          gene.accessions.length > 1
            ? `${gene.symbol} has ${gene.accessions.length} RefSeq transcripts; ${accession} was used. The original server designed against the intersection of all isoforms, so results for this gene may differ.`
            : undefined,
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
    meta.textContent = `${result.sequenceLength.toLocaleString()} nt · ${result.guides.length} shown · ${dropped.toLocaleString()} candidates filtered`;
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
        "<th>miR-30</th><th>miR-E</th><th>97-mer construct</th><th>Notes</th></tr></thead>";
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
    const results: DesignResult[] = [];

    for (const [index, entry] of entries.entries()) {
      const { label, sequence, note } = await resolveEntry(entry);
      setStatus(`Scoring ${label} (${index + 1} of ${entries.length})…`);

      const result = await designGuides(
        { label, sequence },
        engine,
        options,
        undefined,
        (fraction) => {
          els.progress.value = (index + fraction) / entries.length;
        },
      );
      if (note) result.warnings.unshift(note);
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
