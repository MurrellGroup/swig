import { useId, useMemo } from "react";

import {
  FittedLogoGlyph,
  LOGO_MONOSPACE_FONT,
} from "../probability-logo.tsx";
import { probabilityLogoColor } from "../probability-logo.ts";
import type { PhyloUcaHmmAnnotationTrack, PhyloUcaSegmentKind } from "./types.ts";

export interface PhyloUcaAnnotationColumnLayout {
  alignmentColumn: number;
  x: number;
  width: number;
}

interface Props {
  tracks: readonly PhyloUcaHmmAnnotationTrack[];
  columns: readonly PhyloUcaAnnotationColumnLayout[];
  leftInset: number;
  contentWidth: number;
  title: string;
}

const CHARACTERS = ["A", "C", "G", "T", "-"] as const;
const TRACK_HEIGHT = 12;

function kindColor(kind: PhyloUcaSegmentKind): string {
  if (kind === "V") return "#4ba996";
  if (kind === "D") return "#d2a12f";
  if (kind === "J") return "#d76755";
  if (kind === "N") return "#87938e";
  return "#6f7a76";
}

function pointTitle(track: PhyloUcaHmmAnnotationTrack, alignmentColumn: number, probabilities: readonly number[]): string {
  const total = probabilities.reduce((sum, value) => sum + value, 0);
  const detail = probabilities.map((value, index) => value > 1e-5 ? `${CHARACTERS[index]} ${(value * 100).toFixed(3)}%` : "").filter(Boolean).join(", ");
  return `${track.label} · column ${alignmentColumn} · source weight ${(total * 100).toFixed(3)}%${detail ? ` · ${detail}` : ""}`;
}

/**
 * Compact HMM-source rows on the exact same x grid as the posterior logo.
 * Pure template rows draw one reference character per occupied column. Mixed
 * N/unresolved rows stack conditional characters to their unnormalized source
 * mass, so total glyph height is the posterior track occupancy.
 */
export function PhyloUcaHmmAnnotationTracks({ tracks, columns, leftInset, contentWidth, title }: Props) {
  const clipId = `phylo-uca-track-label-${useId().replaceAll(":", "")}`;
  const layoutByColumn = useMemo(() => new Map(columns.map((column) => [column.alignmentColumn, column])), [columns]);
  const width = Math.max(1, leftInset + contentWidth);
  const height = Math.max(1, tracks.length * TRACK_HEIGHT);
  return <svg className="phylo-uca-annotation-svg" xmlns="http://www.w3.org/2000/svg" width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
    <title>{title}</title>
    <desc>Allele rows are pure template tracks. N and unresolved rows are nucleotide mixtures. Letter height is posterior source occupancy.</desc>
    <defs><clipPath id={clipId}><rect x="5" y="0" width={Math.max(1, leftInset - 48)} height={height} /></clipPath></defs>
    <rect width={width} height={height} fill="#ffffff" />
    {tracks.map((track, trackIndex) => {
      const rowY = trackIndex * TRACK_HEIGHT;
      return <g key={track.id}>
        <title>{track.label} · maximum column weight {(track.maximumWeight * 100).toFixed(3)}%</title>
        <rect x="0" y={rowY} width="3" height={TRACK_HEIGHT} fill={kindColor(track.kind)} />
        <text x="6" y={rowY + 8.5} clipPath={`url(#${clipId})`} fill="#26332f" fontFamily={LOGO_MONOSPACE_FONT} fontSize="7.5" fontWeight="700">{track.label}</text>
        <text x={leftInset - 4} y={rowY + 8.5} textAnchor="end" fill="#66736e" fontFamily={LOGO_MONOSPACE_FONT} fontSize="7">{(track.maximumWeight * 100).toFixed(track.maximumWeight >= 0.1 ? 1 : 2)}%</text>
        <line x1="0" x2={width} y1={rowY + TRACK_HEIGHT - 0.25} y2={rowY + TRACK_HEIGHT - 0.25} stroke="#e1e6e3" strokeWidth="0.5" />
        {track.points.flatMap((point) => {
          const layout = layoutByColumn.get(point.alignmentColumn);
          if (!layout) return [];
          const ordered = point.probabilities.map((probability, character) => ({ probability, character })).filter((entry) => entry.probability > 1e-8).sort((left, right) => left.probability - right.probability || left.character - right.character);
          const total = ordered.reduce((sum, entry) => sum + entry.probability, 0);
          let remaining = Math.min(1, total) * TRACK_HEIGHT;
          let stackBottom = rowY + TRACK_HEIGHT;
          return [<g key={`${track.id}-${point.alignmentColumn}`}><title>{pointTitle(track, point.alignmentColumn, point.probabilities)}</title>{ordered.map((entry, entryIndex) => {
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
  </svg>;
}
