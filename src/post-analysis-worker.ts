/// <reference lib="webworker" />

import {
  assignLineages,
  DenoiseAccumulator,
  deduplicate,
  expandSingleLinkage,
  minHashSketch,
  normalizeNt,
  queryRecords,
  sequenceFingerprint,
  type DedupKey,
  type DedupResult,
  type DenoiseOptions,
  type ExpansionOptions,
  type LineageOptions,
  type LineageResult,
  type QueryHit,
  type PostAnalysisRecord,
  type QueryOptions,
} from "./post-analysis-core";

interface IngestRow {
  ordinal: number;
  sequence_id: string;
  sequence: string;
  sequence_alignment: string;
  locus: string;
  v_call: string;
  j_call: string;
  cdr3: string;
  cdr3_aa: string;
  productive: string;
  duplicate_count: string;
}

type Request =
  | { id: number; type: "init"; total: number }
  | { id: number; type: "ingest"; rows: IngestRow[] }
  | { id: number; type: "initSketches" }
  | { id: number; type: "ingestSketches"; rows: Array<{ ordinal: number; sequence: string }> }
  | { id: number; type: "dedup"; key: DedupKey }
  | { id: number; type: "denoiseInit"; options: DenoiseOptions }
  | { id: number; type: "denoiseIngest"; rows: Array<{ ordinal: number; sequence: string }> }
  | { id: number; type: "denoiseFinish" }
  | { id: number; type: "applyDedupFilter" }
  | { id: number; type: "setActiveMask"; mask: Uint8Array | null }
  | { id: number; type: "lineages"; options: LineageOptions; useDedup: boolean }
  | { id: number; type: "query"; queries: string[]; options: QueryOptions }
  | { id: number; type: "expand"; seedOrdinals: number[]; options: ExpansionOptions }
  | { id: number; type: "lineageMembers"; lineageId: number; offset: number; limit: number }
  | { id: number; type: "lineageAssignments" }
  | { id: number; type: "dedupMembers"; representative: number; offset: number; limit: number }
  | { id: number; type: "dedupCounts" }
  | { id: number; type: "clear" };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let records: PostAnalysisRecord[] = [];
let expected = 0;
let currentDedup: DedupResult | undefined;
let currentLineages: LineageResult | undefined;
let packedSketches: Uint32Array | undefined;
let currentActiveMask: Uint8Array | undefined;
let denoiseAccumulator: DenoiseAccumulator | undefined;
const interned = new Map<string, string>();

function intern(value: string): string {
  const existing = interned.get(value);
  if (existing !== undefined) return existing;
  interned.set(value, value);
  return value;
}

function compactLineageResult(result: LineageResult) {
  return {
    summaries: result.summaries,
    lineageCount: result.lineageCount,
    sizeHistogram: result.sizeHistogram,
    vUsage: result.vUsage,
    jUsage: result.jUsage,
    assignedRecords: result.assignedRecords,
    unassignedRecords: result.unassignedRecords,
    candidateComparisons: result.candidateComparisons,
    truncatedCandidates: result.truncatedCandidates,
  };
}

function compactDedupResult(result: DedupResult) {
  return {
    mode: result.mode,
    key: result.key,
    algorithm: result.algorithm,
    inputRecords: result.inputRecords,
    inputAbundance: result.inputAbundance,
    uniqueRecords: result.uniqueRecords,
    collapsedRecords: result.collapsedRecords,
    largestGroups: result.largestGroups,
    partitions: result.partitions,
    candidateComparisons: result.candidateComparisons,
    excludedAmbiguous: result.excludedAmbiguous,
    unresolvedRecords: result.unresolvedRecords,
    warnings: result.warnings,
  };
}

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    let result: unknown;
    if (request.type === "init") {
      records = [];
      expected = request.total;
      currentDedup = undefined;
      currentLineages = undefined;
      packedSketches = undefined;
      currentActiveMask = undefined;
      denoiseAccumulator = undefined;
      interned.clear();
      result = { expected };
    } else if (request.type === "ingest") {
      for (const row of request.rows) {
        const sequence = normalizeNt(row.sequence);
        records.push({
          ordinal: row.ordinal,
          sequenceId: row.sequence_id,
          locus: intern(row.locus),
          vCall: intern(row.v_call),
          jCall: intern(row.j_call),
          cdr3Nt: normalizeNt(row.cdr3),
          cdr3Aa: row.cdr3_aa.toUpperCase().replace(/[^A-Z*]/g, ""),
          productive: row.productive === "T" || row.productive.toLowerCase() === "true",
          sequenceFingerprint: sequenceFingerprint(sequence),
          trimmedFingerprint: sequenceFingerprint(normalizeNt(row.sequence_alignment || row.sequence)),
          inputCount: Math.max(1, Math.floor(Number(row.duplicate_count) || 1)),
        });
      }
      result = { indexed: records.length, expected };
    } else if (request.type === "initSketches") {
      packedSketches = new Uint32Array(expected * 8);
      packedSketches.fill(0xffffffff);
      result = { expected };
    } else if (request.type === "ingestSketches") {
      if (!packedSketches) throw new Error("The VDJ sketch index has not been initialized.");
      for (const row of request.rows) packedSketches.set(minHashSketch(row.sequence), row.ordinal * 8);
      result = { sketched: request.rows.length };
    } else if (request.type === "dedup") {
      currentDedup = deduplicate(records, request.key);
      denoiseAccumulator = undefined;
      currentLineages = undefined;
      currentActiveMask = undefined;
      result = compactDedupResult(currentDedup);
    } else if (request.type === "denoiseInit") {
      denoiseAccumulator = new DenoiseAccumulator(records, request.options);
      currentDedup = undefined;
      currentLineages = undefined;
      currentActiveMask = undefined;
      result = { expected };
    } else if (request.type === "denoiseIngest") {
      if (!denoiseAccumulator) throw new Error("Initialize denoising before streaming VDJ sequences.");
      request.rows.forEach((row) => denoiseAccumulator!.add(row.ordinal, row.sequence));
      result = { ingested: request.rows.length };
    } else if (request.type === "denoiseFinish") {
      if (!denoiseAccumulator) throw new Error("Initialize denoising before finalization.");
      currentDedup = denoiseAccumulator.finish();
      denoiseAccumulator = undefined;
      result = compactDedupResult(currentDedup);
    } else if (request.type === "applyDedupFilter") {
      if (!currentDedup) throw new Error("Run deduplication before applying representative filtering.");
      currentActiveMask = Uint8Array.from(currentDedup.counts, (count) => count > 0 ? 1 : 0);
      currentLineages = undefined;
      result = { mask: currentActiveMask, retained: currentDedup.uniqueRecords };
    } else if (request.type === "setActiveMask") {
      if (request.mask && request.mask.length !== records.length) throw new Error("The downstream working-set mask does not match the AIRR record count.");
      currentActiveMask = request.mask ? request.mask.slice() : undefined;
      currentLineages = undefined;
      let retained = records.length;
      if (currentActiveMask) {
        retained = 0;
        for (const value of currentActiveMask) retained += value ? 1 : 0;
      }
      result = { retained };
    } else if (request.type === "lineages") {
      currentLineages = assignLineages(records, request.options, request.useDedup ? currentDedup : undefined, currentActiveMask);
      result = compactLineageResult(currentLineages);
    } else if (request.type === "query") {
      if (request.options.target === "trimmed" && !packedSketches) throw new Error("Build the VDJ sketch index before querying aligned sequences.");
      result = { hits: queryRecords(records, request.queries, request.options, packedSketches, currentActiveMask) satisfies QueryHit[] };
    } else if (request.type === "expand") {
      result = expandSingleLinkage(records, request.seedOrdinals, request.options, currentActiveMask);
    } else if (request.type === "lineageMembers") {
      if (!currentLineages) throw new Error("Run lineage assignment before opening lineage members.");
      const ordinals: number[] = [];
      let total = 0;
      for (let ordinal = 0; ordinal < currentLineages.assignments.length; ordinal += 1) {
        if (currentActiveMask && !currentActiveMask[ordinal]) continue;
        if (currentLineages.assignments[ordinal] !== request.lineageId) continue;
        if (total >= request.offset && ordinals.length < request.limit) ordinals.push(ordinal);
        total += 1;
      }
      result = { ordinals, total };
    } else if (request.type === "lineageAssignments") {
      if (!currentLineages) throw new Error("Run lineage assignment before exporting clone identifiers.");
      result = { assignments: currentLineages.assignments };
    } else if (request.type === "dedupMembers") {
      if (!currentDedup) throw new Error("Run deduplication before opening duplicate members.");
      const ordinals: number[] = [];
      let total = 0;
      for (let ordinal = 0; ordinal < currentDedup.representatives.length; ordinal += 1) {
        if (currentDedup.representatives[ordinal] !== request.representative) continue;
        if (total >= request.offset && ordinals.length < request.limit) ordinals.push(ordinal);
        total += 1;
      }
      result = { ordinals, total };
    } else if (request.type === "dedupCounts") {
      if (!currentDedup) throw new Error("Run deduplication before exporting deduplicated AIRR data.");
      result = { counts: currentDedup.counts };
    } else {
      records = [];
      currentDedup = undefined;
      currentLineages = undefined;
      packedSketches = undefined;
      currentActiveMask = undefined;
      denoiseAccumulator = undefined;
      interned.clear();
      result = { cleared: true };
    }
    worker.postMessage({ id: request.id, result });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
