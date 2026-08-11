/// <reference lib="webworker" />

import {
  assignLineages,
  deduplicate,
  expandSingleLinkage,
  minHashSketch,
  normalizeNt,
  queryRecords,
  sequenceFingerprint,
  type DedupKey,
  type DedupResult,
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
}

type Request =
  | { id: number; type: "init"; total: number }
  | { id: number; type: "ingest"; rows: IngestRow[] }
  | { id: number; type: "initSketches" }
  | { id: number; type: "ingestSketches"; rows: Array<{ ordinal: number; sequence: string }> }
  | { id: number; type: "dedup"; key: DedupKey }
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
    key: result.key,
    inputRecords: result.inputRecords,
    uniqueRecords: result.uniqueRecords,
    collapsedRecords: result.collapsedRecords,
    largestGroups: result.largestGroups,
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
      currentLineages = undefined;
      result = compactDedupResult(currentDedup);
    } else if (request.type === "lineages") {
      currentLineages = assignLineages(records, request.options, request.useDedup ? currentDedup : undefined);
      result = compactLineageResult(currentLineages);
    } else if (request.type === "query") {
      if (request.options.target === "trimmed" && !packedSketches) throw new Error("Build the VDJ sketch index before querying aligned sequences.");
      result = { hits: queryRecords(records, request.queries, request.options, packedSketches) satisfies QueryHit[] };
    } else if (request.type === "expand") {
      result = expandSingleLinkage(records, request.seedOrdinals, request.options);
    } else if (request.type === "lineageMembers") {
      if (!currentLineages) throw new Error("Run lineage assignment before opening lineage members.");
      const ordinals: number[] = [];
      let total = 0;
      for (let ordinal = 0; ordinal < currentLineages.assignments.length; ordinal += 1) {
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
      interned.clear();
      result = { cleared: true };
    }
    worker.postMessage({ id: request.id, result });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
