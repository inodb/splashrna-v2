/**
 * Transcript diagram: isoforms, the designed region, and where the guides land.
 *
 * The original server drew this with Scribl on a canvas. This is inline SVG instead —
 * Scribl is from 2013, needs jQuery 1.4, and pulls scripts over plain http:// from
 * third-party hosts, which is both broken on an HTTPS page and exactly the kind of
 * dependency this rewrite exists to remove. The data model is simple enough that SVG
 * is less work than reviving it, and it scales, prints and themes for free.
 *
 * Layout, top to bottom:
 *   one lane per RefSeq isoform — CDS as tall blocks, UTRs as short ones
 *   the isoform intersection — the region actually designed against
 *   a tick per predicted shRNA, at its genomic position
 */

import type { GeneAnnotation } from "../data/annotations.js";

const SVG = "http://www.w3.org/2000/svg";

const LEFT = 108; // room for accession labels
const RIGHT = 12;
const WIDTH = 760;
const LANE_HEIGHT = 22;
const CDS_HEIGHT = 11;
const UTR_HEIGHT = 5;
const GUIDE_LANE = 26;
const TOP = 18;

export interface DiagramGuide {
  position: number; // transcript offset, as shown in the results table
  genomic: number;
  guide: string;
  splash: number;
}

function el<K extends keyof SVGElementTagNameMap>(
  name: K,
  attrs: Record<string, string | number>,
): SVGElementTagNameMap[K] {
  const node = document.createElementNS(SVG, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, String(v));
  return node;
}

/**
 * Render the diagram, or return null when there is nothing meaningful to draw —
 * a gene with no isoform annotation, for instance.
 */
export function renderDiagram(
  annotation: GeneAnnotation,
  guides: DiagramGuide[],
): SVGSVGElement | null {
  const spans = [
    ...annotation.isoforms.flatMap((i) => i.exons),
    ...annotation.intersection,
  ];
  if (spans.length === 0) return null;

  const min = Math.min(...spans.map((s) => s.start));
  const max = Math.max(...spans.map((s) => s.end));
  if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) return null;

  const plotWidth = WIDTH - LEFT - RIGHT;
  // Genomic coordinates ascend left to right regardless of strand, matching how the
  // original drew it; the arrow marker carries the direction instead.
  const x = (pos: number) => LEFT + ((pos - min) / (max - min)) * plotWidth;

  const laneCount = annotation.isoforms.length + 1;
  const height = TOP + laneCount * LANE_HEIGHT + GUIDE_LANE + 26;
  const svg = el("svg", {
    xmlns: SVG,
    viewBox: `0 0 ${WIDTH} ${height}`,
    width: "100%",
    height,
    class: "diagram",
    role: "img",
    "aria-label":
      `${annotation.isoforms.length} isoform(s) of Entrez gene ${annotation.entrezId} ` +
      `on ${annotation.chrom} ${annotation.strand} strand, with ${guides.length} predicted shRNA positions`,
  });

  const title = el("title", {});
  title.textContent = `${annotation.chrom}:${min.toLocaleString()}–${max.toLocaleString()} (${annotation.strand} strand)`;
  svg.append(title);

  let y = TOP;

  // --- isoform lanes -----------------------------------------------------------------
  for (const isoform of annotation.isoforms) {
    const first = isoform.exons[0];
    const last = isoform.exons[isoform.exons.length - 1];

    const label = el("text", { x: LEFT - 8, y: y + 4, class: "diagram-label", "text-anchor": "end" });
    label.textContent = isoform.refseq;
    svg.append(label);

    // Intron line spanning the isoform.
    svg.append(
      el("line", {
        x1: x(first.start), y1: y, x2: x(last.end), y2: y, class: "diagram-intron",
      }),
    );

    for (const exon of isoform.exons) {
      const h = exon.kind === "cds" ? CDS_HEIGHT : UTR_HEIGHT;
      svg.append(
        el("rect", {
          x: x(exon.start),
          y: y - h / 2,
          width: Math.max(1, x(exon.end) - x(exon.start)),
          height: h,
          rx: 1,
          class: exon.kind === "cds" ? "diagram-cds" : "diagram-utr",
        }),
      );
    }

    // Direction chevrons on the introns only, as genome browsers draw them. A single
    // arrow at the end of the model collided with the accession labels, and chevrons
    // drawn straight through the exons just added noise.
    const from = x(first.start);
    const to = x(last.end);
    const dir = annotation.strand === "+" ? 1 : -1;
    const blocks = isoform.exons.map((e) => [x(e.start), x(e.end)] as const);
    const overExon = (cx: number) =>
      blocks.some(([a, b]) => cx >= a - 2 && cx <= b + 2);

    for (let cx = from + 30; cx < to - 8; cx += 72) {
      if (overExon(cx)) continue;
      svg.append(
        el("path", {
          d: `M ${cx - 2.5 * dir} ${y - 3} L ${cx + 2.5 * dir} ${y} L ${cx - 2.5 * dir} ${y + 3}`,
          class: "diagram-arrow",
        }),
      );
    }

    y += LANE_HEIGHT;
  }

  // --- the designed region -----------------------------------------------------------
  const intersectionLabel = el("text", {
    x: LEFT - 8, y: y + 4, class: "diagram-label diagram-label-strong", "text-anchor": "end",
  });
  intersectionLabel.textContent = "designed region";
  svg.append(intersectionLabel);

  if (annotation.intersection.length === 0) {
    const none = el("text", { x: LEFT + 4, y: y + 4, class: "diagram-label" });
    none.textContent = "none";
    svg.append(none);
  }
  for (const block of annotation.intersection) {
    svg.append(
      el("rect", {
        x: x(block.start),
        y: y - CDS_HEIGHT / 2,
        width: Math.max(1, x(block.end) - x(block.start)),
        height: CDS_HEIGHT,
        rx: 1,
        class: "diagram-intersection",
      }),
    );
  }
  y += LANE_HEIGHT;

  // --- guides ------------------------------------------------------------------------
  const guideLabel = el("text", {
    x: LEFT - 8, y: y + 8, class: "diagram-label diagram-label-strong", "text-anchor": "end",
  });
  guideLabel.textContent = `shRNAs (${guides.length})`;
  svg.append(guideLabel);

  for (const g of guides) {
    if (!Number.isFinite(g.genomic)) continue;
    const gx = x(g.genomic);
    const tick = el("line", { x1: gx, y1: y, x2: gx, y2: y + 13, class: "diagram-guide" });
    const tip = el("title", {});
    tip.textContent = `${g.guide}\nposition ${g.position}, ${annotation.chrom}:${g.genomic.toLocaleString()}\nSplashRNA ${g.splash.toFixed(3)}`;
    tick.append(tip);
    svg.append(tick);
  }
  y += GUIDE_LANE;

  // --- scale -------------------------------------------------------------------------
  svg.append(el("line", { x1: LEFT, y1: y, x2: LEFT + plotWidth, y2: y, class: "diagram-axis" }));
  for (const [pos, anchor] of [[min, "start"], [max, "end"]] as const) {
    const t = el("text", { x: x(pos), y: y + 13, class: "diagram-tick", "text-anchor": anchor });
    t.textContent = pos.toLocaleString();
    svg.append(t);
  }
  const chrom = el("text", { x: LEFT + plotWidth / 2, y: y + 13, class: "diagram-tick", "text-anchor": "middle" });
  chrom.textContent = `${annotation.chrom} (${annotation.strand})`;
  svg.append(chrom);

  return svg;
}
