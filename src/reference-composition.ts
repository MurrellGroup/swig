import type { LocusKey, SegmentKey } from "./reference-pack";
import { filterReferenceFasta } from "./reference-fasta.ts";

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
  exclusions: Record<string, readonly string[]> = {},
): Partial<Record<SegmentKey, string>> {
  const result: Partial<Record<SegmentKey, string>> = {};
  for (const segment of SEGMENTS) {
    if (!loci.some((locus) => {
      const key = referenceCellKey(locus, segment);
      return Boolean(cells[key] || exclusions[key]?.length);
    })) continue;
    result[segment] = loci.map((locus) => {
      const key = referenceCellKey(locus, segment);
      const source = cells[key]?.text ?? baseline(locus, segment);
      return exclusions[key]?.length ? filterReferenceFasta(source, exclusions[key]).fasta : source;
    }).filter(Boolean).join("\n");
  }
  return result;
}
