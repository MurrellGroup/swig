import { forwardRef, useMemo } from "react";

import {
  normalizedLogoEntries,
  probabilityLogoColor,
  type ProbabilityLogoAlphabet,
  type ProbabilityLogoColumn,
} from "./probability-logo.ts";

interface ProbabilityLogoProps {
  columns: readonly ProbabilityLogoColumn[];
  alphabet: ProbabilityLogoAlphabet;
  title: string;
  description?: string;
  columnWidth?: number;
  stackHeight?: number;
  labelEvery?: number;
}

const SERIF_FONT = "Georgia, 'Times New Roman', Times, serif";

/**
 * Frequency-stack logo. Total stack height is one probability unit for every
 * column; unlike a classical information logo, entropy never rescales it.
 * Nested clipped SVG viewports make every glyph fit its exact probability box.
 */
export const ProbabilityLogo = forwardRef<SVGSVGElement, ProbabilityLogoProps>(function ProbabilityLogo({
  columns,
  alphabet,
  title,
  description = "Every column is a unit-height probability stack; letter height is marginal frequency and is not scaled by entropy.",
  columnWidth = 22,
  stackHeight = 150,
  labelEvery = 5,
}, ref) {
  const normalized = useMemo(() => columns.map((column) => ({
    ...column,
    entries: normalizedLogoEntries(column.entries),
  })), [columns]);
  const left = 38;
  const top = 12;
  const bottom = 34;
  const right = 10;
  const width = Math.max(180, left + right + normalized.length * columnWidth);
  const height = top + stackHeight + bottom;
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  return <div className="probability-logo-viewport">
    <svg ref={ref} className="probability-logo-svg" xmlns="http://www.w3.org/2000/svg" width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
      <title>{title}</title>
      <desc>{description}</desc>
      <rect width={width} height={height} fill="#ffffff" />
      {ticks.map((tick) => {
        const y = top + (1 - tick) * stackHeight;
        return <g key={tick}><line x1={left} x2={width - right} y1={y} y2={y} stroke={tick === 0 ? "#596762" : "#dce2df"} strokeWidth={tick === 0 ? 1 : 0.6} /><text x={left - 6} y={y + 3} textAnchor="end" fill="#61706b" fontFamily={SERIF_FONT} fontSize="8">{tick.toFixed(tick === 0 || tick === 1 ? 0 : 2)}</text></g>;
      })}
      {normalized.map((column, columnIndex) => {
        const x = left + columnIndex * columnWidth + 1;
        const cellWidth = Math.max(1, columnWidth - 2);
        const ordered = [...column.entries].filter((entry) => entry.probability > 0).sort((a, b) => a.probability - b.probability || a.symbol.localeCompare(b.symbol));
        let remaining = stackHeight;
        return <g key={`${column.label}-${columnIndex}`}><title>{column.title ?? `${column.label}: ${column.entries.map((entry) => `${entry.symbol} ${(entry.probability * 100).toFixed(3)}%`).join(", ")}`}</title>{ordered.map((entry, entryIndex) => {
          const glyphHeight = entryIndex === ordered.length - 1 ? remaining : Math.min(remaining, entry.probability * stackHeight);
          const y = top + remaining - glyphHeight;
          remaining -= glyphHeight;
          if (!(glyphHeight > 0)) return null;
          return <svg key={`${entry.symbol}-${entryIndex}`} x={x} y={y} width={cellWidth} height={glyphHeight} viewBox="0 0 100 100" preserveAspectRatio="none" overflow="hidden">
            <text x="50" y="88" textAnchor="middle" fill={entry.color ?? probabilityLogoColor(entry.symbol, alphabet)} fontFamily={SERIF_FONT} fontSize="100" fontWeight="700" textLength="90" lengthAdjust="spacingAndGlyphs">{entry.symbol}</text>
          </svg>;
        })}{(columnIndex % Math.max(1, labelEvery) === 0 || columnIndex === normalized.length - 1) && <text x={x + cellWidth / 2} y={top + stackHeight + 15} textAnchor="middle" fill="#52605b" fontFamily={SERIF_FONT} fontSize="8">{column.label}</text>}</g>;
      })}
      <text x="11" y={top + stackHeight / 2} transform={`rotate(-90 11 ${top + stackHeight / 2})`} textAnchor="middle" fill="#52605b" fontFamily={SERIF_FONT} fontSize="9">Probability</text>
    </svg>
  </div>;
});

export function serializeProbabilityLogoSvg(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  return `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}\n`;
}
