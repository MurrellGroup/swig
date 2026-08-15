import { forwardRef, useMemo } from "react";

import {
  FittedLogoGlyph,
  LOGO_MONOSPACE_FONT,
} from "../probability-logo.tsx";
import { probabilityLogoColor } from "../probability-logo.ts";
import type { PhyloUcaHmmDisplayTrack } from "./hmm-annotation-model.ts";
import type { PhyloUcaSegmentKind } from "./types.ts";

export interface PhyloUcaAnnotationColumnLayout {
  alignmentColumn: number;
  x: number;
  width: number;
}

interface Props {
  tracks: readonly PhyloUcaHmmDisplayTrack[];
  columns: readonly PhyloUcaAnnotationColumnLayout[];
  leftInset: number;
  contentWidth: number;
  title: string;
  labelOffset: number;
}

const CHARACTERS = ["A", "C", "G", "T", "-"] as const;
const TRACK_HEIGHT = 16;

function kindColor(kind: PhyloUcaSegmentKind): string {
  if (kind === "V") return "#4ba996";
  if (kind === "D") return "#d2a12f";
  if (kind === "J") return "#d76755";
  if (kind === "N") return "#87938e";
  return "#6f7a76";
}

function kindBackground(kind: PhyloUcaSegmentKind): string {
  if (kind === "V") return "#e8f6f1";
  if (kind === "D") return "#fff5d8";
  if (kind === "J") return "#ffede8";
  if (kind === "N") return "#eef2f1";
  return "#f4f5f4";
}

function signed(value: number): string {
  return value >= 0 ? `+${value}` : String(value);
}

function trackTitle(track: PhyloUcaHmmDisplayTrack): string {
  const details = [
    `${track.label}; maximum source occupancy ${(track.maximumWeight * 100).toFixed(3)}%`,
    `weighted center column ${track.weightedCenter.toFixed(2)}`,
    `${track.sourceTrackCount} underlying HMM route/register track${track.sourceTrackCount === 1 ? "" : "s"}`,
  ];
  if (track.sourceDOrdinals.length) details.push(`D-use ordinal${track.sourceDOrdinals.length === 1 ? "" : "s"}: ${track.sourceDOrdinals.join(", ")}`);
  if (track.sourceRegistrationOffsets.length) details.push(`alignment register${track.sourceRegistrationOffsets.length === 1 ? "" : "s"}: ${track.sourceRegistrationOffsets.map(signed).join(", ")}`);
  if (track.sourceLabels.length > 1) details.push(`combined rows: ${track.sourceLabels.join("; ")}`);
  return details.join(" · ");
}

function pointTitle(track: PhyloUcaHmmDisplayTrack, alignmentColumn: number, probabilities: readonly number[]): string {
  const total = probabilities.reduce((sum, value) => sum + value, 0);
  const detail = probabilities.map((value, index) => value > 1e-7
    ? `${CHARACTERS[index]} ${(value * 100).toFixed(3)}% overall / ${total > 0 ? (value / total * 100).toFixed(2) : "0.00"}% within track`
    : "").filter(Boolean).join("; ");
  return `${track.label} · alignment column ${alignmentColumn} · source occupancy ${(total * 100).toFixed(3)}%${detail ? ` · nucleotide mass: ${detail}` : ""}`;
}

/**
 * Compact HMM-source rows on the exact same x grid as the posterior logo.
 * Each display row stacks characters to its unnormalized source mass. A
 * template allele normally draws one character per occupied column; when
 * collapsed routes/registers disagree, their characters form a visible
 * mixture. Total glyph height is the posterior track occupancy.
 */
export const PhyloUcaHmmAnnotationTracks = forwardRef<SVGSVGElement, Props>(function PhyloUcaHmmAnnotationTracks({ tracks, columns, leftInset, contentWidth, title, labelOffset }, ref) {
  const layoutByColumn = useMemo(() => new Map(columns.map((column) => [column.alignmentColumn, column])), [columns]);
  const width = Math.max(1, leftInset + contentWidth);
  const height = Math.max(1, tracks.length * TRACK_HEIGHT);
  return <svg ref={ref} className="phylo-uca-annotation-svg" xmlns="http://www.w3.org/2000/svg" width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
    <title>{title}</title>
    <desc>Rows combine equivalent HMM routes and D-alignment registers. Conflicting template registrations and non-templated tracks are nucleotide mixtures. Letter height is posterior source occupancy.</desc>
    <rect width={width} height={height} fill="#ffffff" />
    {tracks.map((track, trackIndex) => {
      const rowY = trackIndex * TRACK_HEIGHT;
      return <g key={track.id}>
        <title>{trackTitle(track)}</title>
        <rect x="0" y={rowY} width={width} height={TRACK_HEIGHT} fill={kindBackground(track.kind)} />
        <line x1="0" x2={width} y1={rowY + TRACK_HEIGHT - 0.25} y2={rowY + TRACK_HEIGHT - 0.25} stroke="#e1e6e3" strokeWidth="0.5" />
        {track.points.flatMap((point) => {
          const layout = layoutByColumn.get(point.alignmentColumn);
          if (!layout) return [];
          const ordered = point.probabilities.map((probability, character) => ({ probability, character })).filter((entry) => entry.probability > 1e-8).sort((left, right) => left.probability - right.probability || left.character - right.character);
          const total = ordered.reduce((sum, entry) => sum + entry.probability, 0);
          let remaining = Math.min(1, total) * TRACK_HEIGHT;
          let stackBottom = rowY + TRACK_HEIGHT;
          const tooltip = pointTitle(track, point.alignmentColumn, point.probabilities);
          return [<g key={`${track.id}-${point.alignmentColumn}`} aria-label={tooltip}><title>{tooltip}</title>{ordered.map((entry, entryIndex) => {
            const glyphHeight = entryIndex === ordered.length - 1 ? remaining : Math.min(remaining, entry.probability * TRACK_HEIGHT);
            const glyphY = stackBottom - glyphHeight;
            remaining -= glyphHeight;
            stackBottom = glyphY;
            if (!(glyphHeight > 0)) return null;
            const symbol = CHARACTERS[entry.character];
            return <FittedLogoGlyph
              key={`${point.alignmentColumn}-${symbol}`}
              symbol={symbol}
              x={leftInset + layout.x}
              y={glyphY}
              width={layout.width}
              height={glyphHeight}
              fill={probabilityLogoColor(symbol, "nucleotide")}
            />;
          })}</g>];
        })}
      </g>;
    })}
    <g className="phylo-uca-track-label-panel" transform={`translate(${Math.max(0, Math.min(contentWidth, labelOffset))} 0)`}>
      {tracks.map((track, trackIndex) => {
        const rowY = trackIndex * TRACK_HEIGHT;
        const shownLabel = track.label.length > 30 ? `${track.label.slice(0, 28)}…` : track.label;
        const tooltip = trackTitle(track);
        return <g key={`label-${track.id}`} aria-label={tooltip}>
          <title>{tooltip}</title>
          <rect x="0" y={rowY} width={leftInset} height={TRACK_HEIGHT} fill={kindBackground(track.kind)} />
          <rect x="0" y={rowY} width="4" height={TRACK_HEIGHT} fill={kindColor(track.kind)} />
          <text x="8" y={rowY + 11} fill="#26332f" fontFamily={LOGO_MONOSPACE_FONT} fontSize="8" fontWeight="700">{shownLabel}</text>
          <text x={leftInset - 7} y={rowY + 11} textAnchor="end" fill="#66736e" fontFamily={LOGO_MONOSPACE_FONT} fontSize="7.5">{(track.maximumWeight * 100).toFixed(track.maximumWeight >= 0.1 ? 1 : 2)}%</text>
          <line x1="0" x2={leftInset} y1={rowY + TRACK_HEIGHT - 0.25} y2={rowY + TRACK_HEIGHT - 0.25} stroke="#d9dfdc" strokeWidth="0.5" />
        </g>;
      })}
      <line x1={leftInset - 0.5} x2={leftInset - 0.5} y1="0" y2={height} stroke="#aebbb5" strokeWidth="1" />
    </g>
  </svg>;
});

export interface PhyloUcaAnnotationSvgViewport {
  x: number;
  y: number;
  width: number;
  height: number;
}

function svgDimensions(svg: SVGSVGElement): { width: number; height: number } {
  return {
    width: Math.max(1, svg.viewBox.baseVal.width || Number(svg.getAttribute("width")) || 1),
    height: Math.max(1, svg.viewBox.baseVal.height || Number(svg.getAttribute("height")) || 1),
  };
}

/** Serialize either the complete track canvas or an exact scrolled crop. */
export function serializePhyloUcaHmmAnnotationSvg(svg: SVGSVGElement, viewport?: PhyloUcaAnnotationSvgViewport): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  const { width: sourceWidth, height: sourceHeight } = svgDimensions(svg);
  const x = viewport ? Math.max(0, Math.min(sourceWidth - 1, viewport.x)) : 0;
  const y = viewport ? Math.max(0, Math.min(sourceHeight - 1, viewport.y)) : 0;
  const width = viewport ? Math.max(1, Math.min(sourceWidth - x, viewport.width)) : sourceWidth;
  const height = viewport ? Math.max(1, Math.min(sourceHeight - y, viewport.height)) : sourceHeight;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  clone.setAttribute("width", String(width));
  clone.setAttribute("height", String(height));
  clone.setAttribute("viewBox", `${x} ${y} ${width} ${height}`);
  const labelPanel = clone.querySelector(".phylo-uca-track-label-panel");
  if (labelPanel) {
    labelPanel.setAttribute("transform", `translate(${viewport ? x : 0} 0)`);
    labelPanel.setAttribute("style", "filter:drop-shadow(2px 0 1px rgba(31,49,43,.18))");
  }
  return `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}\n`;
}

/**
 * Export the HMM source rows and their aligned UCA probability logo as one
 * self-contained SVG. A visible-window export crops the track rows vertically
 * at their inner scroll position and crops both panels to the same horizontal
 * window; the complete corresponding logo remains directly below the tracks.
 */
export function serializePhyloUcaTrackLogoSvg(
  trackSvg: SVGSVGElement,
  logoSvg: SVGSVGElement,
  viewport?: PhyloUcaAnnotationSvgViewport,
): string {
  const trackSize = svgDimensions(trackSvg);
  const logoSize = svgDimensions(logoSvg);
  const sharedWidth = Math.max(1, Math.min(trackSize.width, logoSize.width));
  const x = viewport ? Math.max(0, Math.min(sharedWidth - 1, viewport.x)) : 0;
  const trackY = viewport ? Math.max(0, Math.min(trackSize.height - 1, viewport.y)) : 0;
  const width = viewport ? Math.max(1, Math.min(sharedWidth - x, viewport.width)) : sharedWidth;
  const trackHeight = viewport ? Math.max(1, Math.min(trackSize.height - trackY, viewport.height)) : trackSize.height;
  const height = trackHeight + logoSize.height;

  const namespace = "http://www.w3.org/2000/svg";
  const root = document.createElementNS(namespace, "svg");
  root.setAttribute("xmlns", namespace);
  root.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  root.setAttribute("width", String(width));
  root.setAttribute("height", String(height));
  root.setAttribute("viewBox", `0 0 ${width} ${height}`);
  root.setAttribute("role", "img");
  root.setAttribute("aria-label", "Phylogenetic UCA HMM source tracks and aligned posterior frequency logo");
  const title = document.createElementNS(namespace, "title");
  title.textContent = "Phylogenetic UCA HMM source tracks and posterior frequency logo";
  root.append(title);
  const description = document.createElementNS(namespace, "desc");
  description.textContent = "The upper rows show HMM source occupancy; the aligned UCA probability logo and its numbering and CDR annotations are included below.";
  root.append(description);

  const tracks = trackSvg.cloneNode(true) as SVGSVGElement;
  tracks.setAttribute("x", "0");
  tracks.setAttribute("y", "0");
  tracks.setAttribute("width", String(width));
  tracks.setAttribute("height", String(trackHeight));
  tracks.setAttribute("viewBox", `${x} ${trackY} ${width} ${trackHeight}`);
  const labelPanel = tracks.querySelector(".phylo-uca-track-label-panel");
  if (labelPanel) {
    labelPanel.setAttribute("transform", `translate(${x} 0)`);
    labelPanel.setAttribute("style", "filter:drop-shadow(2px 0 1px rgba(31,49,43,.18))");
  }
  root.append(tracks);

  const logo = logoSvg.cloneNode(true) as SVGSVGElement;
  logo.setAttribute("x", "0");
  logo.setAttribute("y", String(trackHeight));
  logo.setAttribute("width", String(width));
  logo.setAttribute("height", String(logoSize.height));
  logo.setAttribute("viewBox", `${x} 0 ${width} ${logoSize.height}`);
  root.append(logo);

  return `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(root)}\n`;
}
