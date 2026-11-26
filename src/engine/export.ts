/** CSV export of design results, matching the columns the original server produced. */

import type { DesignResult } from "./design.js";

const CITATION =
  "Pelossof R, Fairchild L, et al. (2017). Prediction of potent shRNAs with a sequential classification algorithm. Nat Biotechnol 35(4):350-353.";
const METHODS_NOTE =
  "Please add to your methods: The shRNAs were designed using the SplashRNA algorithm (Pelossof and Fairchild et al., 2017).";

const COLUMNS = [
  "Feature",
  "shRNA.name",
  "Position",
  "Antisense.Guide.Sequence",
  "SplashRNA",
  "miR30.score",
  "miRE.score",
  "97mer.construct",
  "Warnings",
  "Human.22mer.match.genes",
  "Mouse.22mer.match.genes",
];

function escapeCell(value: string | number): string {
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

export function toCsv(results: DesignResult[], generatedAt = new Date()): string {
  const rows: (string | number)[][] = [
    [`Generated ${generatedAt.toISOString().replace("T", " ").slice(0, 19)} UTC`],
    [METHODS_NOTE],
    [`Ref: ${CITATION}`],
    [],
    COLUMNS,
  ];

  for (const result of results) {
    if (result.guides.length === 0) {
      rows.push([result.label, "no predictions", "", "", "", "", "", "", result.warnings.join("; ")]);
      continue;
    }
    result.guides.forEach((guide, i) => {
      rows.push([
        result.label,
        `${result.label}_${guide.position}`,
        guide.position,
        guide.guide,
        guide.splash.toFixed(6),
        guide.mir30.toFixed(6),
        guide.mirE.toFixed(6),
        guide.construct,
        [...guide.warnings, ...(i === 0 ? result.warnings : [])].join("; "),
        guide.offTargetHuman.join(";"),
        guide.offTargetMouse.join(";"),
      ]);
    });
  }

  return rows.map((row) => row.map(escapeCell).join(",")).join("\n");
}

/** Trigger a browser download of the given CSV text. */
export function downloadCsv(csv: string, filename = "splashrna-predictions.csv"): void {
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  const href = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = href;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(href);
}
