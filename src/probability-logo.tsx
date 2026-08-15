import {
  forwardRef,
  useMemo,
} from "react";

import { logoGlyphRun } from "./logo-glyphs.ts";
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
  /** Left coordinate shared with an aligned annotation surface. */
  leftInset?: number;
  /** Omit the overflow wrapper when a parent owns the shared scroll surface. */
  embedded?: boolean;
}

interface FittedLogoGlyphProps {
  symbol: string;
  x: number;
  y: number;
  width: number;
  height: number;
  fill: string;
}

/**
 * Axis and track labels use an ordinary monospace stack. Probability glyphs
 * themselves are embedded paths and therefore do not depend on browser fonts.
 */
export const LOGO_MONOSPACE_FONT = "'Liberation Mono', 'DejaVu Sans Mono', 'Courier New', monospace";

/** Fit literal contour extrema exactly into a target probability rectangle. */
export function FittedLogoGlyph({ symbol, x, y, width, height, fill }: FittedLogoGlyphProps) {
  const run = logoGlyphRun(symbol);
  if (!run || !(width > 0) || !(height > 0)) return null;
  const scaleX = width / run.width;
  const scaleY = height / run.height;
  const translateX = x - run.xMin * scaleX;
  const translateY = y + run.yMax * scaleY;
  return <g
    fill={fill}
    transform={`matrix(${scaleX} 0 0 ${-scaleY} ${translateX} ${translateY})`}
    aria-hidden="true"
  >{run.paths.map((path, index) => <path key={`${path.x}-${index}`} d={path.d} transform={path.x ? `translate(${path.x} 0)` : undefined} />)}</g>;
}

/**
 * Frequency-stack logo. Total stack height is one probability unit for every
 * column; unlike a classical information logo, entropy never rescales it.
 * Each glyph's painted bounds fill its probability rectangle exactly: there is
 * no inter-column padding and descenders are included rather than clipped.
 */
export const ProbabilityLogo = forwardRef<SVGSVGElement, ProbabilityLogoProps>(function ProbabilityLogo({
  columns,
  alphabet,
  title,
  description = "Every column is a unit-height probability stack; letter height is marginal frequency and is not scaled by entropy.",
  columnWidth = 22,
  stackHeight = 150,
  labelEvery = 5,
  leftInset = 38,
  embedded = false,
}, ref) {
  const normalized = useMemo(() => columns.map((column) => ({
    ...column,
    entries: normalizedLogoEntries(column.entries),
  })), [columns]);
  const bottom = 24;
  const width = Math.max(1, leftInset + normalized.length * columnWidth);
  const height = stackHeight + bottom;
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  const svg = <svg ref={ref} className="probability-logo-svg" xmlns="http://www.w3.org/2000/svg" width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
    <title>{title}</title>
    <desc>{description}</desc>
    <rect width={width} height={height} fill="#ffffff" />
    {ticks.map((tick) => {
      const y = (1 - tick) * stackHeight;
      const labelY = tick === 1 ? 8 : tick === 0 ? stackHeight - 2 : y + 3;
      return <g key={tick}><line x1={leftInset} x2={width} y1={y} y2={y} stroke={tick === 0 ? "#596762" : "#dce2df"} strokeWidth={tick === 0 ? 1 : 0.6} /><text x={leftInset - 6} y={labelY} textAnchor="end" fill="#61706b" fontFamily={LOGO_MONOSPACE_FONT} fontSize="8">{tick.toFixed(tick === 0 || tick === 1 ? 0 : 2)}</text></g>;
    })}
    {normalized.map((column, columnIndex) => {
      const x = leftInset + columnIndex * columnWidth;
      const ordered = [...column.entries].filter((entry) => entry.probability > 0).sort((a, b) => a.probability - b.probability || a.symbol.localeCompare(b.symbol));
      let remaining = stackHeight;
      return <g key={`${column.label}-${columnIndex}`}><title>{column.title ?? `${column.label}: ${column.entries.map((entry) => `${entry.symbol} ${(entry.probability * 100).toFixed(3)}%`).join(", ")}`}</title>{ordered.map((entry, entryIndex) => {
        const glyphHeight = entryIndex === ordered.length - 1 ? remaining : Math.min(remaining, entry.probability * stackHeight);
        const y = remaining - glyphHeight;
        remaining -= glyphHeight;
        if (!(glyphHeight > 0)) return null;
        return <FittedLogoGlyph
          key={`${entry.symbol}-${entryIndex}`}
          symbol={entry.symbol}
          x={x}
          y={y}
          width={columnWidth}
          height={glyphHeight}
          fill={entry.color ?? probabilityLogoColor(entry.symbol, alphabet)}
        />;
      })}{(columnIndex % Math.max(1, labelEvery) === 0 || columnIndex === normalized.length - 1) && <text x={x + columnWidth / 2} y={stackHeight + 15} textAnchor="middle" fill="#52605b" fontFamily={LOGO_MONOSPACE_FONT} fontSize="8">{column.label}</text>}</g>;
    })}
    <text x="11" y={stackHeight / 2} transform={`rotate(-90 11 ${stackHeight / 2})`} textAnchor="middle" fill="#52605b" fontFamily={LOGO_MONOSPACE_FONT} fontSize="9">Probability</text>
  </svg>;
  return embedded ? svg : <div className="probability-logo-viewport">{svg}</div>;
});

export function serializeProbabilityLogoSvg(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  return `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}\n`;
}
