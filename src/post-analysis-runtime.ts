import type { AirrResultStore, AirrScanRow } from "./result-store";
import type {
  CollapseMode,
  DedupKey,
  DenoiseOptions,
  ExpansionOptions,
  LineageOptions,
  LineageSummary,
  QueryHit,
  QueryOptions,
} from "./post-analysis-core";
import type { DatasetScope } from "./study-design";

export interface DedupDashboard {
  mode: CollapseMode;
  key: DedupKey;
  algorithm: string;
  inputRecords: number;
  inputAbundance: number;
  uniqueRecords: number;
  collapsedRecords: number;
  largestGroups: Array<{ ordinal: number; count: number }>;
  partitions: number;
  candidateComparisons: number;
  indelMergedVariants: number;
  substitutionMergedVariants: number;
  excludedAmbiguous: number;
  unresolvedRecords: number;
  warnings: string[];
}

export interface LineageDashboard {
  summaries: LineageSummary[];
  lineageCount: number;
  sizeHistogram: Array<{ label: string; count: number }>;
  vUsage: Array<{ call: string; lineages: number; abundance: number }>;
  jUsage: Array<{ call: string; lineages: number; abundance: number }>;
  assignedRecords: number;
  unassignedRecords: number;
  candidateComparisons: number;
  truncatedCandidates: number;
}

type WorkerResult = Record<string, unknown>;

export class PostAnalysisRuntime {
  private readonly worker = new Worker(new URL("./post-analysis-worker.ts", import.meta.url), { type: "module" });
  private readonly pending = new Map<number, { resolve: (value: unknown) => void; reject: (error: Error) => void }>();
  private nextId = 1;
  private indexed = false;
  private indexing: Promise<void> | null = null;
  private sketched = false;
  private sketching: Promise<void> | null = null;

  constructor(private readonly store: AirrResultStore) {
    this.worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string }>) => {
      const pending = this.pending.get(event.data.id);
      if (!pending) return;
      this.pending.delete(event.data.id);
      if (event.data.error) pending.reject(new Error(event.data.error));
      else pending.resolve(event.data.result);
    };
    this.worker.onerror = (event) => {
      const error = new Error(event.message || "The post-analysis worker stopped unexpectedly.");
      this.pending.forEach((pending) => pending.reject(error));
      this.pending.clear();
    };
  }

  private request<T>(message: Record<string, unknown>): Promise<T> {
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject });
      this.worker.postMessage({ id, ...message });
    });
  }

  ensureIndexed(onProgress?: (processed: number, total: number) => void, signal?: AbortSignal): Promise<void> {
    if (this.indexed) return Promise.resolve();
    if (this.indexing) return this.indexing;
    this.indexing = (async () => {
      await this.request({ type: "init", total: this.store.count });
      const fields = [
        "sequence_id", "sequence", "sequence_alignment", "locus", "v_call", "j_call",
        "cdr3", "cdr3_aa", "productive", "duplicate_count", "swig_dataset_id", "sample_id",
        "subject_id", "swig_cohort", "swig_timepoint",
      ];
      await this.store.scanAirrRows(fields, async (rows) => {
        if (signal?.aborted) throw new DOMException("Post-analysis was cancelled.", "AbortError");
        await this.request({ type: "ingest", rows: rows.map(toWorkerRow) });
      }, { batchSize: 2_000, onProgress, signal });
      this.indexed = true;
    })().finally(() => { this.indexing = null; });
    return this.indexing;
  }

  async deduplicate(key: DedupKey, unresolvedPolicy: "discard" | "retain" = "discard", scope: DatasetScope = "global"): Promise<DedupDashboard> {
    await this.ensureIndexed();
    return this.request<DedupDashboard>({ type: "dedup", key, unresolvedPolicy, scope });
  }

  async denoise(options: DenoiseOptions, onProgress?: (processed: number, total: number) => void, signal?: AbortSignal): Promise<DedupDashboard> {
    await this.ensureIndexed(onProgress, signal);
    await this.request({ type: "denoiseInit", options });
    await this.store.scanAirrRows(["sequence_alignment", "sequence", "v_sequence_start", "j_sequence_end"], async (rows) => {
      if (signal?.aborted) throw new DOMException("Denoising was cancelled.", "AbortError");
      await this.request({
        type: "denoiseIngest",
        rows: rows.map((row) => ({ ordinal: row.ordinal, sequence: denoiseVdjSequence(row) })),
      });
    }, { batchSize: 2_000, onProgress, signal });
    return this.request<DedupDashboard>({ type: "denoiseFinish" });
  }

  async applyDedupFilter(): Promise<{ mask: Uint8Array; retained: number }> {
    await this.ensureIndexed();
    return this.request({ type: "applyDedupFilter" });
  }

  async setActiveMask(mask: Uint8Array | null): Promise<{ retained: number }> {
    await this.ensureIndexed();
    return this.request({ type: "setActiveMask", mask });
  }

  async activeMask(): Promise<Uint8Array | null> {
    await this.ensureIndexed();
    const result = await this.request<{ mask: Uint8Array | null }>({ type: "activeMask" });
    return result.mask;
  }

  async assignLineages(options: LineageOptions, useDedup: boolean): Promise<LineageDashboard> {
    await this.ensureIndexed();
    return this.request<LineageDashboard>({ type: "lineages", options, useDedup });
  }

  async query(queries: string[], options: QueryOptions): Promise<QueryHit[]> {
    await this.ensureIndexed();
    if (options.target === "trimmed") await this.ensureSketches();
    const result = await this.request<{ hits: QueryHit[] }>({ type: "query", queries, options });
    return result.hits;
  }

  private ensureSketches(): Promise<void> {
    if (this.sketched) return Promise.resolve();
    if (this.sketching) return this.sketching;
    this.sketching = (async () => {
      await this.request({ type: "initSketches" });
      await this.store.scanAirrRows(["sequence_alignment", "sequence"], async (rows) => {
        await this.request({
          type: "ingestSketches",
          rows: rows.map((row) => ({ ordinal: row.ordinal, sequence: row.values.sequence_alignment || row.values.sequence })),
        });
      }, { batchSize: 2_000 });
      this.sketched = true;
    })().finally(() => { this.sketching = null; });
    return this.sketching;
  }

  async expand(seedOrdinals: number[], options: ExpansionOptions): Promise<{ ordinals: number[]; comparisons: number; capped: boolean }> {
    await this.ensureIndexed();
    return this.request({ type: "expand", seedOrdinals, options });
  }

  async lineageMembers(lineageId: number, offset = 0, limit = 500): Promise<{ ordinals: number[]; total: number }> {
    return this.request({ type: "lineageMembers", lineageId, offset, limit });
  }

  async lineageAssignments(): Promise<Int32Array> {
    const result = await this.request<{ assignments: Int32Array }>({ type: "lineageAssignments" });
    return result.assignments;
  }

  async dedupMembers(representative: number, offset = 0, limit = 500): Promise<{ ordinals: number[]; total: number }> {
    return this.request({ type: "dedupMembers", representative, offset, limit });
  }

  async dedupCounts(): Promise<Uint32Array> {
    const result = await this.request<{ counts: Uint32Array }>({ type: "dedupCounts" });
    return result.counts;
  }

  async dedupState(): Promise<{ counts: Uint32Array; representatives: Int32Array }> {
    return this.request({ type: "dedupState" });
  }

  async restoreState(state: {
    activeMask?: Uint8Array | null;
    dedup?: { dashboard: DedupDashboard; counts: Uint32Array; representatives: Int32Array };
    lineages?: { dashboard: LineageDashboard; assignments: Int32Array };
  }): Promise<void> {
    await this.ensureIndexed();
    await this.request({ type: "restoreState", ...state });
  }

  terminate() {
    this.worker.terminate();
    const error = new Error("Post-analysis was closed.");
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
  }
}

function toWorkerRow(row: AirrScanRow): WorkerResult {
  return {
    ordinal: row.ordinal,
    sequence_id: row.values.sequence_id,
    swig_dataset_id: row.values.swig_dataset_id,
    sample_id: row.values.sample_id,
    subject_id: row.values.subject_id,
    swig_cohort: row.values.swig_cohort,
    swig_timepoint: row.values.swig_timepoint,
    sequence: row.values.sequence,
    sequence_alignment: row.values.sequence_alignment,
    locus: row.values.locus,
    v_call: row.values.v_call,
    j_call: row.values.j_call,
    cdr3: row.values.cdr3,
    cdr3_aa: row.values.cdr3_aa,
    productive: row.values.productive,
    duplicate_count: row.values.duplicate_count,
  };
}

function denoiseVdjSequence(row: AirrScanRow): string {
  const raw = row.values.sequence ?? "";
  const start = Math.floor(Number(row.values.v_sequence_start));
  const end = Math.floor(Number(row.values.j_sequence_end));
  if (raw && Number.isFinite(start) && Number.isFinite(end) && start >= 1 && end >= start && end <= raw.length) {
    return raw.slice(start - 1, end);
  }
  return row.values.sequence_alignment || raw;
}
