import type { AirrScanRow } from "../result-store.ts";
import type { RefinementInputRow, RefinementSegment } from "./types.ts";

function numeric(value: string): number | null {
  const parsed = Number(value);
  return value !== "" && Number.isFinite(parsed) ? parsed : null;
}

function normalizedIdentity(value: number | null): number | null {
  if (value === null) return null;
  const normalized = value > 1 ? value / 100 : value;
  return normalized >= 0 && normalized <= 1 ? normalized : null;
}

function bestAlternativeIdentity(value: string): number | null {
  let best: number | null = null;
  for (const entry of value.split(";")) {
    const identity = normalizedIdentity(numeric(entry.split("|")[2] ?? ""));
    if (identity !== null && (best === null || identity > best)) best = identity;
  }
  return best;
}

function bestIdentity(primary: string, alternatives: string): number | null {
  const direct = normalizedIdentity(numeric(primary));
  const alternate = bestAlternativeIdentity(alternatives);
  if (direct === null) return alternate;
  if (alternate === null) return direct;
  return Math.max(direct, alternate);
}

export function refinementInputFields(segment: RefinementSegment): string[] {
  const prefix = segment.toLowerCase();
  return [...new Set([
    "sequence_id", "swig_dataset_id", "sample_id", "subject_id", "locus", "duplicate_count",
    "v_identity", "v_alternatives",
    `${prefix}_call`, `${prefix}_score`, `${prefix}_identity`, `${prefix}_alternatives`,
  ])];
}

export function toRefinementInputRow(row: AirrScanRow, segment: RefinementSegment): RefinementInputRow {
  const prefix = segment.toLowerCase();
  const vIdentity = bestIdentity(row.values.v_identity ?? "", row.values.v_alternatives ?? "");
  const segmentIdentity = bestIdentity(row.values[`${prefix}_identity`] ?? "", row.values[`${prefix}_alternatives`] ?? "");
  const closestIdentity = vIdentity ?? segmentIdentity;
  return {
    ordinal: row.ordinal,
    sequenceId: row.values.sequence_id ?? "",
    datasetId: row.values.swig_dataset_id ?? "",
    sampleId: row.values.sample_id ?? "",
    subjectId: row.values.subject_id ?? "",
    locus: row.values.locus ?? "",
    call: row.values[`${prefix}_call`] ?? "",
    score: numeric(row.values[`${prefix}_score`] ?? ""),
    identity: numeric(row.values[`${prefix}_identity`] ?? ""),
    shm: closestIdentity === null ? null : Math.max(0, 1 - closestIdentity),
    alternatives: row.values[`${prefix}_alternatives`] ?? "",
    abundance: Math.max(1, Math.floor(Number(row.values.duplicate_count) || 1)),
  };
}
