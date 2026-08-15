import type { CompiledReferences } from "../reference-pack.ts";
import type { AirrResultStore } from "../result-store.ts";
import type {
  AlleleRefinementOptions,
  AlleleRefinementResult,
  RefinementSegment,
  SegmentRefinementResult,
} from "./types.ts";
import { refinementInputFields, toRefinementInputRow } from "./input.ts";

interface Progress {
  processed: number;
  total: number;
  phase: string;
  segment?: RefinementSegment;
}

interface Pending {
  resolve: (value: unknown) => void;
  reject: (error: Error) => void;
  onProgress?: (progress: Progress) => void;
}

export class AlleleRefinementRuntime {
  private readonly worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
  private readonly pending = new Map<number, Pending>();
  private nextId = 1;

  constructor() {
    this.worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string; progress?: Progress }>) => {
      const pending = this.pending.get(event.data.id);
      if (!pending) return;
      if (event.data.progress) {
        pending.onProgress?.(event.data.progress);
        return;
      }
      this.pending.delete(event.data.id);
      if (event.data.error) pending.reject(new Error(event.data.error));
      else pending.resolve(event.data.result);
    };
    this.worker.onerror = (event) => {
      const error = new Error(event.message || "The repertoire-level allele worker stopped unexpectedly.");
      this.pending.forEach((pending) => pending.reject(error));
      this.pending.clear();
    };
  }

  private request<T>(message: Record<string, unknown>, onProgress?: (progress: Progress) => void): Promise<T> {
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject, onProgress });
      this.worker.postMessage({ id, ...message });
    });
  }

  async run(
    store: AirrResultStore,
    references: CompiledReferences,
    options: AlleleRefinementOptions,
    includeMask: Uint8Array | null,
    onProgress?: (progress: Progress) => void,
    signal?: AbortSignal,
  ): Promise<AlleleRefinementResult> {
    if (includeMask && includeMask.length !== store.count) throw new Error("The allele-refinement mask does not match the AIRR record count.");
    await this.request({ type: "init", totalRecords: store.count, options });
    const segments: Partial<Record<RefinementSegment, SegmentRefinementResult>> = {};
    for (const segment of options.segments) {
      if (signal?.aborted) throw new DOMException("Allele refinement was cancelled.", "AbortError");
      const fasta = references[segment];
      if (!fasta.trim()) continue;
      await this.request({ type: "beginSegment", segment, fasta }, (progress) => onProgress?.({ ...progress, segment }));
      const fields = refinementInputFields(segment);
      await store.scanAirrRows(fields, async (rows) => {
        if (signal?.aborted) throw new DOMException("Allele refinement was cancelled.", "AbortError");
        await this.request({ type: "ingest", rows: rows.map((row) => toRefinementInputRow(row, segment)) });
      }, {
        batchSize: 2_000,
        includeMask: includeMask ?? undefined,
        signal,
        onProgress: (processed, total) => onProgress?.({ processed, total, phase: `Building sparse ${segment} evidence matrix`, segment }),
      });
      segments[segment] = await this.request<SegmentRefinementResult>({ type: "finishSegment" }, (progress) => onProgress?.({ ...progress, segment }));
    }
    let activeRecords = store.count;
    if (includeMask) {
      activeRecords = 0;
      for (const value of includeMask) activeRecords += value ? 1 : 0;
    }
    const warnings: string[] = [];
    if (options.model === "active-set") warnings.push("The fast hurdle model estimates repertoire-active usage, not literal genomic presence; a genomically present but silent allele is not identifiable from expressed reads alone.");
    if (options.weighting === "abundance") warnings.push("Read-abundance weighting lets clonal expansion influence the inferred mixture. Unique active records are the conservative default.");
    if (segments.D) warnings.push("D refinement is exploratory: short templated spans, exonuclease trimming, and N addition make D emissions less identifiable than V or J.");
    const truncated = Object.values(segments).reduce((sum, segment) => sum + (segment?.truncatedRows ?? 0), 0);
    if (truncated) warnings.push(`${truncated.toLocaleString()} sparse evidence rows reached the configured candidate cap. Explicit reported calls were retained; lower-weight neighbour candidates were truncated.`);
    const unconverged = Object.values(segments).flatMap((segment) => segment?.models ?? []).filter((model) => !model.converged).length;
    if (unconverged) warnings.push(`${unconverged.toLocaleString()} independent repertoire model${unconverged === 1 ? " did" : "s did"} not reach the configured convergence tolerance. Inspect the model export or increase the iteration limit.`);
    return { version: 1, options, totalRecords: store.count, activeRecords, segments, runAt: new Date().toISOString(), warnings };
  }

  terminate() {
    this.worker.terminate();
    const error = new Error("The repertoire-level allele worker was closed.");
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
  }
}
