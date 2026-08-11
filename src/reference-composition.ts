import type { LocusKey, SegmentKey } from "./reference-pack";

const SEGMENTS: SegmentKey[] = ["V", "D", "J", "C"];

export function referenceCellKey(locus: LocusKey, segment: SegmentKey): string {
  return `${locus}:${segment}`;
}

export function segmentAppliesToLocus(locus: LocusKey, segment: SegmentKey): boolean {
  return segment !== "D" || locus === "IGH" || locus === "TRB" || locus === "TRD";
}

export function composeReferenceOverrides(
  loci: LocusKey[],
  cells: Record<string, { text: string }>,
  baseline: (locus: LocusKey, segment: SegmentKey) => string,
): Partial<Record<SegmentKey, string>> {
  const result: Partial<Record<SegmentKey, string>> = {};
  for (const segment of SEGMENTS) {
    if (!loci.some((locus) => cells[referenceCellKey(locus, segment)])) continue;
    result[segment] = loci.map((locus) => cells[referenceCellKey(locus, segment)]?.text
      ?? baseline(locus, segment)).filter(Boolean).join("\n");
  }
  return result;
}
