import { packSessionVector, unpackSessionVector, type SessionVector } from "../session-state.ts";
import type {
  AlleleRefinementResult,
  RefinementSegment,
  SavedAlleleRefinement,
  SavedSegmentRefinement,
  SegmentRefinementResult,
} from "./types.ts";

function packed<T extends "i32" | "f32">(value: Int32Array | Float32Array, type: T): SessionVector & { type: T } {
  const result = packSessionVector(value);
  if (result.type !== type) throw new Error("Unexpected allele-refinement session vector type.");
  return result as SessionVector & { type: T };
}

function saveSegment(result: SegmentRefinementResult): SavedSegmentRefinement {
  return {
    ...result,
    mapNode: packed(result.mapNode, "i32"),
    mapProbability: packed(result.mapProbability, "f32"),
    posteriorEntropy: packed(result.posteriorEntropy, "f32"),
    localTopNode: packed(result.localTopNode, "i32"),
    localTopProbability: packed(result.localTopProbability, "f32"),
  };
}

function restoreSegment(result: SavedSegmentRefinement): SegmentRefinementResult {
  const mapNode = unpackSessionVector(result.mapNode);
  const mapProbability = unpackSessionVector(result.mapProbability);
  const posteriorEntropy = unpackSessionVector(result.posteriorEntropy);
  const localTopNode = unpackSessionVector(result.localTopNode);
  const localTopProbability = unpackSessionVector(result.localTopProbability);
  if (!(mapNode instanceof Int32Array) || !(mapProbability instanceof Float32Array)
    || !(posteriorEntropy instanceof Float32Array) || !(localTopNode instanceof Int32Array)
    || !(localTopProbability instanceof Float32Array)) {
    throw new Error("The saved repertoire-level allele result has incompatible vector types.");
  }
  return { ...result, mapNode, mapProbability, posteriorEntropy, localTopNode, localTopProbability };
}

export function saveAlleleRefinement(
  result: AlleleRefinementResult,
  applied: boolean,
  applyMinimumPosterior: number,
): SavedAlleleRefinement {
  const segments: SavedAlleleRefinement["segments"] = {};
  for (const segment of ["V", "D", "J"] as RefinementSegment[]) {
    const value = result.segments[segment];
    if (value) segments[segment] = saveSegment(value);
  }
  return { ...result, segments, applied, applyMinimumPosterior };
}

export function restoreAlleleRefinement(saved: SavedAlleleRefinement): AlleleRefinementResult {
  const segments: AlleleRefinementResult["segments"] = {};
  for (const segment of ["V", "D", "J"] as RefinementSegment[]) {
    const value = saved.segments[segment];
    if (value) segments[segment] = restoreSegment(value);
  }
  return {
    version: 1,
    options: saved.options,
    totalRecords: saved.totalRecords,
    activeRecords: saved.activeRecords,
    segments,
    runAt: saved.runAt,
    warnings: saved.warnings,
  };
}
