export type AirrRow = Record<string, string>;
export type AlignmentMode = "nt" | "aa";

interface AlignmentPair {
  query: string;
  reference: string;
}

export interface TrackFeature {
  key: string;
  label: string;
  start: number;
  end: number;
  left: number;
  width: number;
  kind: "region" | "segment";
}

const CODONS: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S",
  TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
  CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

export const REGIONS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"] as const;
export const ALIGNMENT_SEGMENTS = [
  { key: "v", label: "V", color: "var(--segment-v)" },
  { key: "d", label: "D", color: "var(--segment-d)" },
  { key: "j", label: "J", color: "var(--segment-j)" },
  { key: "c", label: "C", color: "var(--segment-c)" },
] as const;

function positiveModulo(value: number, modulus: number): number {
  return ((value % modulus) + modulus) % modulus;
}

function translateAligned(sequence: string): string {
  let protein = "";
  for (let index = 0; index + 2 < sequence.length; index += 3) {
    const codon = sequence.slice(index, index + 3).toUpperCase();
    protein += codon.includes("-") ? "-" : CODONS[codon] ?? "X";
  }
  return protein;
}

function queryOffsetColumn(alignedQuery: string, basesToSkip: number): number {
  if (!basesToSkip) return 0;
  let observed = 0;
  for (let column = 0; column < alignedQuery.length; column += 1) {
    if (alignedQuery[column] !== "-") observed += 1;
    if (observed === basesToSkip) return column + 1;
  }
  return alignedQuery.length;
}

export function biologicalSegmentAlignment(
  nucleotideQuery: string,
  nucleotideReference: string,
  queryStartOneBased: number,
  sequenceFrameOneBased: number,
): AlignmentPair | null {
  if (!nucleotideQuery || !nucleotideReference || !queryStartOneBased || !sequenceFrameOneBased) return null;
  const queryStart = queryStartOneBased - 1;
  const frame = sequenceFrameOneBased - 1;
  const basesToSkip = positiveModulo(frame - queryStart, 3);
  const column = queryOffsetColumn(nucleotideQuery, basesToSkip);
  const query = translateAligned(nucleotideQuery.slice(column));
  const reference = translateAligned(nucleotideReference.slice(column));
  return query || reference ? { query, reference } : null;
}

export function buildTrackFeatures(row: AirrRow): { regions: TrackFeature[]; segments: TrackFeature[] } {
  const length = row.sequence?.length ?? 0;
  if (!length) return { regions: [], segments: [] };
  const feature = (key: string, label: string, start: number, end: number, kind: TrackFeature["kind"]): TrackFeature => ({
    key,
    label,
    start,
    end,
    left: (start - 1) / length * 100,
    width: Math.max(0.5, (end - start + 1) / length * 100),
    kind,
  });
  const regions = REGIONS.flatMap((region) => {
    const start = Number(row[`${region}_start`]);
    const end = Number(row[`${region}_end`]);
    return start > 0 && end >= start ? [feature(region, region.toUpperCase(), start, end, "region")] : [];
  });
  const segments = ALIGNMENT_SEGMENTS.flatMap((segment) => {
    const start = Number(row[`${segment.key}_sequence_start`]);
    const end = Number(row[`${segment.key}_sequence_end`]);
    return row[`${segment.key}_call`] && start > 0 && end >= start
      ? [feature(segment.key, segment.label, start, end, "segment")]
      : [];
  });
  return { regions, segments };
}
