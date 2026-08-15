import type { PhyloUcaHmmAnnotationPoint, PhyloUcaHmmAnnotationTrack, PhyloUcaSegmentKind } from "./types.ts";

const MASS_EPSILON = 1e-10;

export interface PhyloUcaHmmDisplayTrack extends PhyloUcaHmmAnnotationTrack {
  /** Raw HMM route/register rows combined into this display row. */
  sourceTrackCount: number;
  sourceLabels: string[];
  sourceDOrdinals: number[];
  sourceRegistrationOffsets: number[];
  /** Alignment-column center of mass, used for left-to-right display order. */
  weightedCenter: number;
  /** Occupancy summed over columns; this is an ordering diagnostic, not a probability. */
  integratedWeight: number;
}

function pointWeight(point: PhyloUcaHmmAnnotationPoint): number {
  return point.probabilities.reduce((sum, value) => sum + value, 0);
}

function trackMass(track: PhyloUcaHmmAnnotationTrack): number {
  return track.points.reduce((sum, point) => sum + pointWeight(point), 0);
}

function combinedTrack(
  id: string,
  kind: PhyloUcaSegmentKind,
  label: string,
  source: readonly PhyloUcaHmmAnnotationTrack[],
  call?: string,
): PhyloUcaHmmDisplayTrack {
  const byColumn = new Map<number, Float64Array>();
  for (const track of source) for (const point of track.points) {
    let masses = byColumn.get(point.alignmentColumn);
    if (!masses) {
      masses = new Float64Array(5);
      byColumn.set(point.alignmentColumn, masses);
    }
    for (let character = 0; character < 5; character += 1) masses[character] += point.probabilities[character];
  }
  let maximumWeight = 0;
  let integratedWeight = 0;
  let weightedColumnTotal = 0;
  let pure = true;
  const points = [...byColumn.entries()].sort(([left], [right]) => left - right).flatMap(([alignmentColumn, raw]) => {
    const probabilities = Array.from(raw) as [number, number, number, number, number];
    const total = probabilities.reduce((sum, value) => sum + value, 0);
    if (!(total > MASS_EPSILON)) return [];
    maximumWeight = Math.max(maximumWeight, total);
    integratedWeight += total;
    weightedColumnTotal += alignmentColumn * total;
    if (probabilities.filter((value) => value > MASS_EPSILON).length > 1) pure = false;
    return [{ alignmentColumn, probabilities }];
  });
  const sourceDOrdinals = [...new Set(source.map((track) => track.dOrdinal).filter((value): value is number => value !== undefined))].sort((left, right) => left - right);
  const sourceRegistrationOffsets = [...new Set(source.map((track) => track.registrationOffset).filter((value): value is number => value !== undefined))].sort((left, right) => left - right);
  return {
    id,
    kind,
    label,
    call,
    dOrdinal: sourceDOrdinals.length === 1 ? sourceDOrdinals[0] : undefined,
    pure,
    points,
    maximumWeight,
    sourceTrackCount: source.length,
    sourceLabels: [...new Set(source.map((track) => track.label))],
    sourceDOrdinals,
    sourceRegistrationOffsets,
    weightedCenter: integratedWeight > 0 ? weightedColumnTotal / integratedWeight : 0,
    integratedWeight,
  };
}

function groupTemplateTracks(tracks: readonly PhyloUcaHmmAnnotationTrack[]): PhyloUcaHmmDisplayTrack[] {
  const groups = new Map<string, PhyloUcaHmmAnnotationTrack[]>();
  for (const track of tracks) {
    if (track.kind === "N") continue;
    const key = track.call ? `${track.kind}|${track.call}` : `${track.kind}|${track.id}`;
    const group = groups.get(key);
    if (group) group.push(track); else groups.set(key, [track]);
  }
  return [...groups.entries()].map(([key, group]) => {
    const first = group[0];
    const label = first.call ? `${first.kind} · ${first.call}` : first.label;
    return combinedTrack(`display|${key}`, first.kind, label, group, first.call);
  });
}

interface DBlock {
  ordinal: number;
  center: number;
  integratedWeight: number;
}

function dBlocks(tracks: readonly PhyloUcaHmmAnnotationTrack[]): DBlock[] {
  const byOrdinal = new Map<number, { weightedColumns: number; integratedWeight: number }>();
  for (const track of tracks) {
    if (track.kind !== "D" || track.dOrdinal === undefined) continue;
    const current = byOrdinal.get(track.dOrdinal) ?? { weightedColumns: 0, integratedWeight: 0 };
    for (const point of track.points) {
      const weight = pointWeight(point);
      current.weightedColumns += point.alignmentColumn * weight;
      current.integratedWeight += weight;
    }
    byOrdinal.set(track.dOrdinal, current);
  }
  return [...byOrdinal.entries()].map(([ordinal, value]) => ({
    ordinal,
    center: value.integratedWeight > 0 ? value.weightedColumns / value.integratedWeight : 0,
    integratedWeight: value.integratedWeight,
  })).sort((left, right) => left.center - right.center || right.integratedWeight - left.integratedWeight || left.ordinal - right.ordinal);
}

function nSlot(track: PhyloUcaHmmAnnotationTrack, point: PhyloUcaHmmAnnotationPoint, blocks: readonly DBlock[]): number {
  if (!blocks.length) return 0;
  if (/J-entry boundary/i.test(track.label)) return blocks.length;
  if (/V-trim boundary/i.test(track.label) || track.label === "N0") return 0;
  const match = /^N(\d+)$/.exec(track.label);
  const dUsed = match ? Number(match[1]) : track.dOrdinal;
  if (dUsed !== undefined && /unresolved D/i.test(track.label) === false) {
    const block = blocks.findIndex((value) => value.ordinal === dUsed);
    if (block >= 0) return Math.min(blocks.length, block + 1);
  }
  let slot = 0;
  while (slot < blocks.length && point.alignmentColumn > blocks[slot].center) slot += 1;
  return slot;
}

function spatialNTracks(tracks: readonly PhyloUcaHmmAnnotationTrack[], blocks: readonly DBlock[]): Array<PhyloUcaHmmDisplayTrack | undefined> {
  const sources = tracks.filter((track) => track.kind === "N");
  const slots = Array.from({ length: blocks.length + 1 }, () => new Map<string, PhyloUcaHmmAnnotationTrack>());
  for (const track of sources) for (const point of track.points) {
    const slot = nSlot(track, point, blocks);
    const splitId = `${track.id}|display-slot-${slot}`;
    const split = slots[slot].get(splitId);
    if (split) split.points.push(point);
    else slots[slot].set(splitId, { ...track, id: splitId, points: [point] });
  }
  return slots.map((slot, index) => {
    const source = [...slot.values()];
    return source.length ? combinedTrack(`display|NT${index + 1}`, "N", `NT${index + 1}`, source) : undefined;
  });
}

function dominantDBlock(
  display: PhyloUcaHmmDisplayTrack,
  raw: readonly PhyloUcaHmmAnnotationTrack[],
  blocks: readonly DBlock[],
): number {
  if (!blocks.length) return 0;
  const massByOrdinal = new Map<number, number>();
  for (const track of raw) {
    if (track.kind !== "D" || track.call !== display.call || track.dOrdinal === undefined) continue;
    massByOrdinal.set(track.dOrdinal, (massByOrdinal.get(track.dOrdinal) ?? 0) + trackMass(track));
  }
  let bestIndex = -1;
  let bestMass = -1;
  blocks.forEach((block, index) => {
    const mass = massByOrdinal.get(block.ordinal) ?? 0;
    if (mass > bestMass) { bestMass = mass; bestIndex = index; }
  });
  if (bestMass > 0) return bestIndex;
  let nearest = 0;
  for (let index = 1; index < blocks.length; index += 1) {
    if (Math.abs(blocks[index].center - display.weightedCenter) < Math.abs(blocks[nearest].center - display.weightedCenter)) nearest = index;
  }
  return nearest;
}

function alternativeOrder(left: PhyloUcaHmmDisplayTrack, right: PhyloUcaHmmDisplayTrack): number {
  const center = left.weightedCenter - right.weightedCenter;
  if (Math.abs(center) > 1) return center;
  return right.maximumWeight - left.maximumWeight || right.integratedWeight - left.integratedWeight || left.label.localeCompare(right.label, undefined, { numeric: true });
}

/**
 * Visualization-only aggregation. It never changes the HMM likelihood or UCA
 * posterior: route/register masses already produced by inference are summed by
 * allele, and non-template mass is organized into spatial junction slots.
 */
export function collapseAndOrderHmmAnnotationTracks(
  tracks: readonly PhyloUcaHmmAnnotationTrack[],
  mode: "viterbi" | "marginalized" = "viterbi",
): PhyloUcaHmmDisplayTrack[] {
  const template = groupTemplateTracks(tracks);
  if (mode === "marginalized") {
    const nSources = tracks.filter((track) => track.kind === "N");
    const nt = nSources.length ? combinedTrack("display|NT", "N", "NT · all non-template mass", nSources) : undefined;
    const byKind = (kind: PhyloUcaSegmentKind) => template.filter((track) => track.kind === kind).sort(alternativeOrder);
    const unknown = template.filter((track) => !["V", "D", "J"].includes(track.kind)).sort(alternativeOrder);
    return [...(nt ? [nt] : []), ...byKind("V"), ...byKind("D"), ...byKind("J"), ...unknown];
  }
  const blocks = dBlocks(tracks);
  const nTracks = spatialNTracks(tracks, blocks);
  const v = template.filter((track) => track.kind === "V").sort(alternativeOrder);
  const j = template.filter((track) => track.kind === "J").sort(alternativeOrder);
  const unknown = template.filter((track) => track.kind !== "V" && track.kind !== "D" && track.kind !== "J").sort(alternativeOrder);
  const dByBlock = Array.from({ length: Math.max(1, blocks.length) }, () => [] as PhyloUcaHmmDisplayTrack[]);
  for (const track of template.filter((value) => value.kind === "D")) dByBlock[dominantDBlock(track, tracks, blocks)].push(track);
  dByBlock.forEach((group) => group.sort(alternativeOrder));
  const ordered = [...v];
  if (!blocks.length) {
    ordered.push(...nTracks.filter((track): track is PhyloUcaHmmDisplayTrack => Boolean(track)), ...dByBlock.flat(), ...unknown, ...j);
    return ordered;
  }
  for (let block = 0; block < blocks.length; block += 1) {
    const precedingN = nTracks[block];
    if (precedingN) ordered.push(precedingN);
    ordered.push(...dByBlock[block]);
  }
  const finalN = nTracks[blocks.length];
  if (finalN) ordered.push(finalN);
  ordered.push(...unknown, ...j);
  return ordered;
}
