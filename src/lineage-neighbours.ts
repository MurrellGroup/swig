import { inferLineageGermline, type InferredLineageGermline, type LineageGermlineMethod } from "./lineage-alignment.ts";
import { minHashSketch, sketchSimilarity } from "./post-analysis-core.ts";
import type { AirrDetailRow, AirrResultStore } from "./result-store.ts";
import { datasetScopeValue, type DatasetScope } from "./study-design.ts";

export interface LineageGermlineSketchIndex {
  lineageCount: number;
  sketches: Uint32Array;
  representativeOrdinals: Int32Array;
  loci: string[];
  studyGroups: string[];
  indexedRows: number;
  representedLineages: number;
}

export interface GermlineScreenCandidate {
  lineageId: number;
  sketchSimilarity: number;
  representativeOrdinal: number;
}

export interface GermlineNeighbourScore extends GermlineScreenCandidate {
  /** Original source lineage whose inferred UCA gave the best exact score. */
  sourceLineageId: number;
  germlineIdentity: number;
  sourceGermline: InferredLineageGermline;
  candidateGermline: InferredLineageGermline;
  candidateRowsLoaded: number;
  candidateTotalRows: number;
}

function normalizeAligned(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replaceAll(".", "-").replace(/[^ACGTN-]/g, "N");
}

/** One-row UCA proxy used only for cheap sketch screening, never final scoring. */
function provisionalAncestor(queryValue: string, germlineValue: string): { sequence: string; quality: number } {
  const query = normalizeAligned(queryValue);
  const germline = normalizeAligned(germlineValue);
  const length = Math.max(query.length, germline.length);
  let sequence = "";
  let known = 0;
  let mismatch = 0;
  for (let column = 0; column < length; column += 1) {
    const q = query[column] ?? "-";
    const g = germline[column] ?? "-";
    if (/^[ACGT]$/.test(g)) {
      sequence += g;
      known += 1;
      if (/^[ACGT]$/.test(q) && q !== g) mismatch += 1;
    } else if (/^[ACGT]$/.test(q)) sequence += q;
    else sequence += "N";
  }
  return { sequence, quality: known * 1_000 + length - mismatch * 25 };
}

/**
 * Build one compact eight-word UCA screen signature per assigned lineage in a
 * single bounded-memory AIRR scan. The least-mutated/highest-coverage member
 * supplies the screen sketch; exact hits are subsequently verified against a
 * the user-selected exact inferred-germline method.
 */
export async function buildLineageGermlineSketchIndex(
  store: AirrResultStore,
  assignments: Int32Array,
  lineageCount: number,
  scope: DatasetScope,
  activeMask?: Uint8Array | null,
  onProgress?: (processed: number, total: number) => void,
): Promise<LineageGermlineSketchIndex> {
  if (assignments.length < store.count) throw new Error("Lineage assignments do not cover the AIRR table.");
  const sketches = new Uint32Array((lineageCount + 1) * 8);
  sketches.fill(0xffffffff);
  const representativeOrdinals = new Int32Array(lineageCount + 1);
  representativeOrdinals.fill(-1);
  const quality = new Float64Array(lineageCount + 1);
  quality.fill(Number.NEGATIVE_INFINITY);
  const loci = Array<string>(lineageCount + 1).fill("");
  const studyGroups = Array<string>(lineageCount + 1).fill("");
  let indexedRows = 0;
  await store.scanAirrRows(
    ["sequence_alignment", "sequence", "germline_alignment", "locus", "swig_dataset_id", "sample_id", "subject_id", "swig_cohort"],
    async (rows) => {
      for (const row of rows) {
        const lineageId = assignments[row.ordinal];
        if (!(lineageId > 0) || lineageId > lineageCount) continue;
        indexedRows += 1;
        const provisional = provisionalAncestor(row.values.sequence_alignment || row.values.sequence, row.values.germline_alignment);
        if (provisional.quality <= quality[lineageId]) continue;
        quality[lineageId] = provisional.quality;
        representativeOrdinals[lineageId] = row.ordinal;
        sketches.set(minHashSketch(provisional.sequence), lineageId * 8);
        loci[lineageId] = row.values.locus;
        studyGroups[lineageId] = datasetScopeValue({
          datasetId: row.values.swig_dataset_id,
          sampleId: row.values.sample_id,
          subjectId: row.values.subject_id,
          cohort: row.values.swig_cohort,
        }, scope);
      }
    },
    { batchSize: 2_000, includeMask: activeMask ?? undefined, onProgress },
  );
  let representedLineages = 0;
  for (let lineageId = 1; lineageId <= lineageCount; lineageId += 1) if (representativeOrdinals[lineageId] >= 0) representedLineages += 1;
  return { lineageCount, sketches, representativeOrdinals, loci, studyGroups, indexedRows, representedLineages };
}

interface HeapEntry extends GermlineScreenCandidate {}

function heapLess(left: HeapEntry, right: HeapEntry): boolean {
  return left.sketchSimilarity < right.sketchSimilarity ||
    (left.sketchSimilarity === right.sketchSimilarity && left.lineageId > right.lineageId);
}

function heapPush(heap: HeapEntry[], value: HeapEntry) {
  heap.push(value);
  let index = heap.length - 1;
  while (index > 0) {
    const parent = Math.floor((index - 1) / 2);
    if (!heapLess(heap[index], heap[parent])) break;
    [heap[index], heap[parent]] = [heap[parent], heap[index]];
    index = parent;
  }
}

function heapReplaceRoot(heap: HeapEntry[], value: HeapEntry) {
  heap[0] = value;
  let index = 0;
  while (true) {
    const left = index * 2 + 1;
    const right = left + 1;
    let smallest = index;
    if (left < heap.length && heapLess(heap[left], heap[smallest])) smallest = left;
    if (right < heap.length && heapLess(heap[right], heap[smallest])) smallest = right;
    if (smallest === index) break;
    [heap[index], heap[smallest]] = [heap[smallest], heap[index]];
    index = smallest;
  }
}

export function screenLineageGermlineCandidates(
  index: LineageGermlineSketchIndex,
  source: InferredLineageGermline,
  excludedLineageIds: Iterable<number>,
  locus: string,
  studyGroup: string,
  maximumCandidates = 200,
): GermlineScreenCandidate[] {
  const excluded = new Set(excludedLineageIds);
  const sourceSketch = minHashSketch(source.trimmedUca);
  const heap: HeapEntry[] = [];
  for (let lineageId = 1; lineageId <= index.lineageCount; lineageId += 1) {
    const representativeOrdinal = index.representativeOrdinals[lineageId];
    if (representativeOrdinal < 0 || excluded.has(lineageId)) continue;
    if (locus && index.loci[lineageId] !== locus) continue;
    if (studyGroup && index.studyGroups[lineageId] !== studyGroup) continue;
    const similarity = sketchSimilarity(sourceSketch, index.sketches.subarray(lineageId * 8, lineageId * 8 + 8));
    const candidate = { lineageId, sketchSimilarity: similarity, representativeOrdinal };
    if (heap.length < maximumCandidates) heapPush(heap, candidate);
    else if (similarity > heap[0].sketchSimilarity || (similarity === heap[0].sketchSimilarity && lineageId < heap[0].lineageId)) heapReplaceRoot(heap, candidate);
  }
  return heap.sort((left, right) => right.sketchSimilarity - left.sketchSimilarity || left.lineageId - right.lineageId);
}

/** Banded Levenshtein distance with N as an unknown/wildcard state. */
export function inferredGermlineIdentity(leftValue: string, rightValue: string, minimumIdentity = 0): number {
  const left = normalizeAligned(leftValue).replace(/^-+|-+$/g, "");
  const right = normalizeAligned(rightValue).replace(/^-+|-+$/g, "");
  const denominator = Math.max(left.length, right.length);
  if (!denominator) return 0;
  const maximum = Math.max(Math.abs(left.length - right.length), Math.floor((1 - minimumIdentity) * denominator + 1e-9));
  let previous = new Int32Array(right.length + 1);
  let current = new Int32Array(right.length + 1);
  for (let column = 0; column <= right.length; column += 1) previous[column] = column;
  for (let row = 1; row <= left.length; row += 1) {
    current.fill(maximum + 1);
    current[0] = row;
    const start = Math.max(1, row - maximum);
    const end = Math.min(right.length, row + maximum);
    let rowMinimum = maximum + 1;
    for (let column = start; column <= end; column += 1) {
      const a = left[row - 1];
      const b = right[column - 1];
      const substitution = a === b || a === "N" || b === "N" ? 0 : 1;
      current[column] = Math.min(previous[column] + 1, current[column - 1] + 1, previous[column - 1] + substitution);
      rowMinimum = Math.min(rowMinimum, current[column]);
    }
    if (rowMinimum > maximum) return Math.max(0, 1 - (maximum + 1) / denominator);
    [previous, current] = [current, previous];
  }
  const distance = previous[right.length];
  return Math.max(0, 1 - distance / denominator);
}

export function scoreGermlineCandidate(
  sourceRows: AirrDetailRow[],
  candidateRows: AirrDetailRow[],
  screen: GermlineScreenCandidate,
  candidateTotalRows: number,
  minimumIdentity: number,
  sourceLineageId = 0,
  method: LineageGermlineMethod = "closest",
): GermlineNeighbourScore | null {
  const sourceGermline = inferLineageGermline(sourceRows, method);
  const candidateGermline = inferLineageGermline(candidateRows, method);
  const germlineIdentity = inferredGermlineIdentity(sourceGermline.trimmedUca, candidateGermline.trimmedUca, minimumIdentity);
  if (germlineIdentity + 1e-12 < minimumIdentity) return null;
  return {
    ...screen,
    sourceLineageId,
    germlineIdentity,
    sourceGermline,
    candidateGermline,
    candidateRowsLoaded: candidateRows.length,
    candidateTotalRows,
  };
}
