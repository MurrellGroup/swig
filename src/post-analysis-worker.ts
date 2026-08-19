/// <reference lib="webworker" />

import {
  assignLineages,
  createExactDedupPlan,
  DenoiseAccumulator,
  expandSingleLinkage,
  finishExactDedupPlan,
  findLineageNeighbours,
  minHashSketch,
  queryRecords,
  runDenoisePartitionJob,
  runExactDedupJob,
  type DedupKey,
  type DedupResult,
  type DenoisePartitionJob,
  type DenoisePartitionJobResult,
  type DenoiseOptions,
  type ExactDedupJob,
  type ExactDedupJobResult,
  type ExpansionOptions,
  type LineageOptions,
  type LineageNeighbourOptions,
  type LineageResult,
  type QueryHit,
  type PostAnalysisRecord,
  type QueryOptions,
} from "./post-analysis-core";
import { airrRowToPostAnalysisRecord } from "./post-analysis-record";
import type { DatasetScope } from "./study-design";
import { applyCallOverrides } from "./allele-refinement/apply";
import type { AlleleReassignmentPolicy } from "./allele-refinement/types";

interface IngestRow {
  ordinal: number;
  sequence_id: string;
  swig_dataset_id: string;
  sample_id: string;
  subject_id: string;
  swig_cohort: string;
  swig_timepoint: string;
  swig_compartment: string;
  sequence: string;
  sequence_alignment: string;
  locus: string;
  v_call: string;
  j_call: string;
  c_call: string;
  cdr3: string;
  cdr3_aa: string;
  productive: string;
  duplicate_count: string;
}

type Request =
  | { id: number; type: "init"; total: number; doubleDMask: Uint8Array }
  | { id: number; type: "ingest"; rows: IngestRow[] }
  | { id: number; type: "initSketches" }
  | { id: number; type: "ingestSketches"; rows: Array<{ ordinal: number; sequence: string }> }
  | { id: number; type: "dedup"; key: DedupKey; unresolvedPolicy: "discard" | "retain"; scope: DatasetScope; respectConstantCall: boolean; workers: number }
  | { id: number; type: "denoiseInit"; options: DenoiseOptions }
  | { id: number; type: "denoiseIngest"; rows: Array<{ ordinal: number; sequence: string }> }
  | { id: number; type: "denoiseFinish"; workers: number }
  | { id: number; type: "applyDedupFilter" }
  | { id: number; type: "setActiveMask"; mask: Uint8Array | null }
  | { id: number; type: "activeMask" }
  | { id: number; type: "setCallOverrides"; reassignmentPolicy: AlleleReassignmentPolicy; minimumPosterior: number; v?: { labels: string[]; mapNode: Int32Array; probability: Float32Array }; j?: { labels: string[]; mapNode: Int32Array; probability: Float32Array } }
  | { id: number; type: "lineages"; options: LineageOptions; useDedup: boolean }
  | { id: number; type: "query"; queries: string[]; options: QueryOptions }
  | { id: number; type: "expand"; seedOrdinals: number[]; options: ExpansionOptions }
  | { id: number; type: "lineageMembers"; lineageId: number; offset: number; limit: number }
  | { id: number; type: "lineageMembersMany"; lineageIds: number[]; limitPerLineage: number }
  | { id: number; type: "lineageNeighbours"; options: LineageNeighbourOptions; useDedup: boolean }
  | { id: number; type: "lineageAssignments" }
  | { id: number; type: "dedupMembers"; representative: number; offset: number; limit: number }
  | { id: number; type: "dedupCounts" }
  | { id: number; type: "dedupState" }
  | { id: number; type: "restoreState"; activeMask?: Uint8Array | null; dedup?: { dashboard: Omit<DedupResult, "counts" | "representatives">; counts: Uint32Array; representatives: Int32Array }; lineages?: { dashboard: Omit<LineageResult, "assignments">; assignments: Int32Array } }
  | { id: number; type: "clear" };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let records: PostAnalysisRecord[] = [];
let expected = 0;
let currentDedup: DedupResult | undefined;
let currentLineages: LineageResult | undefined;
let packedSketches: Uint32Array | undefined;
let currentActiveMask: Uint8Array | undefined;
let currentDoubleDMask: Uint8Array | undefined;
let denoiseAccumulator: DenoiseAccumulator | undefined;
const interned = new Map<string, string>();

function intern(value: string): string {
  const existing = interned.get(value);
  if (existing !== undefined) return existing;
  interned.set(value, value);
  return value;
}

type PartitionTask =
  | { kind: "denoise"; job: DenoisePartitionJob }
  | { kind: "exact"; job: ExactDedupJob };

function runPartitionTasks(
  tasks: PartitionTask[],
  requestedWorkers: number,
  onCompleted?: (result: DenoisePartitionJobResult | ExactDedupJobResult) => void,
): Promise<Array<DenoisePartitionJobResult | ExactDedupJobResult>> {
  if (!tasks.length) return Promise.resolve([]);
  const workerCount = Math.max(1, Math.min(Math.floor(requestedWorkers) || 1, tasks.length));
  if (workerCount === 1) {
    return Promise.resolve(tasks.map((task) => {
      const result = task.kind === "denoise" ? runDenoisePartitionJob(task.job) : runExactDedupJob(task.job);
      onCompleted?.(result);
      return result;
    }));
  }
  return new Promise((resolve, reject) => {
    const children: Worker[] = [];
    const results: Array<DenoisePartitionJobResult | ExactDedupJobResult> = new Array(tasks.length);
    let nextTask = 0;
    let completed = 0;
    let settled = false;
    const close = () => children.forEach((child) => child.terminate());
    const fail = (error: Error) => {
      if (settled) return;
      settled = true;
      close();
      reject(error);
    };
    const dispatch = (child: Worker) => {
      if (nextTask >= tasks.length) return;
      const taskIndex = nextTask++;
      child.postMessage({ taskIndex, ...tasks[taskIndex] });
    };
    for (let index = 0; index < workerCount; index += 1) {
      const child = new Worker(new URL("./post-analysis-partition-worker.ts", import.meta.url), { type: "module" });
      child.onmessage = (event: MessageEvent<{ taskIndex: number; result?: DenoisePartitionJobResult | ExactDedupJobResult; error?: string }>) => {
        if (settled) return;
        if (event.data.error || !event.data.result) {
          fail(new Error(event.data.error || "A post-analysis partition worker returned no result."));
          return;
        }
        results[event.data.taskIndex] = event.data.result;
        onCompleted?.(event.data.result);
        completed += 1;
        if (completed === tasks.length) {
          settled = true;
          close();
          resolve(results);
        } else dispatch(child);
      };
      child.onerror = (event) => fail(new Error(event.message || "A post-analysis partition worker stopped unexpectedly."));
      children.push(child);
      dispatch(child);
    }
  });
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
    indelMergedVariants: result.indelMergedVariants,
    substitutionMergedVariants: result.substitutionMergedVariants,
    excludedAmbiguous: result.excludedAmbiguous,
    unresolvedRecords: result.unresolvedRecords,
    warnings: result.warnings,
  };
}

function refreshDoubleDLineageSummaries(result: LineageResult, dedup?: DedupResult, activeMask?: Uint8Array) {
  const summaries = new Map(result.summaries.map((summary) => {
    summary.doubleDPositiveMembers = 0;
    summary.doubleDPositiveAbundance = 0;
    return [summary.id, summary] as const;
  }));
  for (let index = 0; index < records.length; index += 1) {
    if (activeMask && !activeMask[index]) continue;
    const weight = dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
    if (!weight || !currentDoubleDMask?.[index]) continue;
    const summary = summaries.get(result.assignments[index]);
    if (!summary) continue;
    summary.doubleDPositiveMembers = (summary.doubleDPositiveMembers ?? 0) + 1;
    summary.doubleDPositiveAbundance = (summary.doubleDPositiveAbundance ?? 0) + weight;
  }
}

worker.onmessage = async (event: MessageEvent<Request>) => {
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
      currentDoubleDMask = request.doubleDMask;
      denoiseAccumulator = undefined;
      interned.clear();
      result = { expected };
    } else if (request.type === "ingest") {
      for (const row of request.rows) {
        const {ordinal,...values}=row;
        records.push(airrRowToPostAnalysisRecord({ordinal,values},intern));
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
      const shardCount = records.length >= 10_000 ? request.workers : 1;
      const plan = createExactDedupPlan(records, request.key, request.unresolvedPolicy, request.scope, request.respectConstantCall, shardCount);
      const taskResults = await runPartitionTasks(plan.jobs.map((job) => ({ kind: "exact" as const, job })), request.workers);
      currentDedup = finishExactDedupPlan(plan, taskResults as ExactDedupJobResult[]);
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
      const jobs = denoiseAccumulator.preparePartitionJobs();
      const variantWork = Math.max(1, jobs.reduce((total, job) => total + (job.options.mode === "fad"
        ? job.variants.length + job.variants.reduce((count, variant) => count + (variant.count >= job.options.minimumParentCount ? 1 : 0), 0)
        : job.variants.length), 0));
      let processedVariantWork = 0;
      worker.postMessage({ id: request.id, progress: { processed: 0, total: variantWork, phase: "variants" } });
      const parallelWorkers = jobs.reduce((total, job) => total + job.variants.length, 0) >= 500 ? request.workers : 1;
      const taskResults = await runPartitionTasks(jobs.map((job) => ({ kind: "denoise" as const, job })), parallelWorkers, (completed) => {
        if (!("variantWork" in completed)) return;
        processedVariantWork += completed.variantWork;
        worker.postMessage({ id: request.id, progress: { processed: Math.min(processedVariantWork, variantWork), total: variantWork, phase: "variants" } });
      });
      if (!jobs.length) worker.postMessage({ id: request.id, progress: { processed: variantWork, total: variantWork, phase: "variants" } });
      currentDedup = denoiseAccumulator.finishWithPartitionResults(taskResults as DenoisePartitionJobResult[], (processed, total, phase) => {
        worker.postMessage({ id: request.id, progress: { processed, total, phase } });
      }, false);
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
    } else if (request.type === "activeMask") {
      result = { mask: currentActiveMask?.slice() ?? null };
    } else if (request.type === "setCallOverrides") {
      const applied = applyCallOverrides(records, request.v, request.j, request.reassignmentPolicy, request.minimumPosterior, intern);
      currentLineages = undefined;
      result = applied;
    } else if (request.type === "lineages") {
      currentLineages = assignLineages(records, request.options, request.useDedup ? currentDedup : undefined, currentActiveMask, currentDoubleDMask);
      result = compactLineageResult(currentLineages);
    } else if (request.type === "query") {
      if (request.options.target === "trimmed" && !packedSketches) throw new Error("Build the VDJ sketch index before querying aligned sequences.");
      result = { hits: queryRecords(records, request.queries, request.options, packedSketches, currentActiveMask, currentLineages?.assignments, currentLineages?.lineageCount ?? 0) satisfies QueryHit[] };
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
    } else if (request.type === "lineageMembersMany") {
      if (!currentLineages) throw new Error("Run lineage assignment before opening lineage members.");
      const selected = new Set(request.lineageIds.filter((value) => value > 0));
      const byLineage = new Map<number, { lineageId: number; ordinals: number[]; total: number }>();
      selected.forEach((lineageId) => byLineage.set(lineageId, { lineageId, ordinals: [], total: 0 }));
      const limit = Math.max(1, request.limitPerLineage);
      const reservoirSlot = (ordinal: number, lineageId: number, seen: number) => {
        let value = (ordinal + 1) ^ Math.imul(lineageId + 1, 0x9e3779b1);
        value ^= value >>> 16;
        value = Math.imul(value, 0x7feb352d);
        value ^= value >>> 15;
        value = Math.imul(value, 0x846ca68b);
        value ^= value >>> 16;
        return (value >>> 0) % seen;
      };
      for (let ordinal = 0; ordinal < currentLineages.assignments.length; ordinal += 1) {
        if (currentActiveMask && !currentActiveMask[ordinal]) continue;
        const entry = byLineage.get(currentLineages.assignments[ordinal]);
        if (!entry) continue;
        entry.total += 1;
        if (entry.ordinals.length < limit) entry.ordinals.push(ordinal);
        else {
          const slot = reservoirSlot(ordinal, entry.lineageId, entry.total);
          if (slot < limit) entry.ordinals[slot] = ordinal;
        }
      }
      byLineage.forEach((entry) => entry.ordinals.sort((left, right) => left - right));
      result = { members: [...byLineage.values()] };
    } else if (request.type === "lineageNeighbours") {
      if (!currentLineages) throw new Error("Run lineage assignment before searching for lineage neighbours.");
      result = findLineageNeighbours(records, currentLineages.assignments, request.options, request.useDedup ? currentDedup : undefined, currentActiveMask);
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
    } else if (request.type === "dedupState") {
      if (!currentDedup) throw new Error("Run deduplication before saving its state.");
      result = { counts: currentDedup.counts, representatives: currentDedup.representatives };
    } else if (request.type === "restoreState") {
      if (request.activeMask && request.activeMask.length !== records.length) throw new Error("The saved working-set mask does not match the linked AIRR table.");
      currentActiveMask = request.activeMask ? request.activeMask.slice() : undefined;
      currentDedup = request.dedup ? { ...request.dedup.dashboard, counts: request.dedup.counts, representatives: request.dedup.representatives } : undefined;
      currentLineages = request.lineages ? { ...request.lineages.dashboard, assignments: request.lineages.assignments } : undefined;
      if (currentLineages) refreshDoubleDLineageSummaries(currentLineages, currentDedup, currentActiveMask);
      result = { restored: true, lineages: currentLineages ? compactLineageResult(currentLineages) : null };
    } else {
      records = [];
      currentDedup = undefined;
      currentLineages = undefined;
      packedSketches = undefined;
      currentActiveMask = undefined;
      currentDoubleDMask = undefined;
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
