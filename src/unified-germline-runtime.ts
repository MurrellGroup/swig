import type { AirrScanRow } from "./result-store.ts";
import type { UnifiedDashboard, UnifiedOptions } from "./unified-germline.ts";

interface Progress {
  processed: number;
  total: number;
  phase: string;
}

interface Pending {
  resolve: (value: unknown) => void;
  reject: (error: Error) => void;
  onProgress?: (progress: Progress) => void;
}

export class UnifiedGermlineRuntime {
  private worker: Worker | null = null;
  private readonly pending = new Map<number, Pending>();
  private nextId = 1;
  private generation = 0;
  private closed = false;

  constructor() {}

  private installWorker() {
    const generation = this.generation;
    const worker = new Worker(new URL("./unified-germline-worker.ts", import.meta.url), { type: "module" });
    this.worker = worker;
    worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string; progress?: Progress }>) => {
      if (generation !== this.generation) return;
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
    worker.onerror = (event) => {
      if (generation !== this.generation) return;
      const error = new Error(event.message || "The joint germline worker stopped unexpectedly.");
      this.pending.forEach((pending) => pending.reject(error));
      this.pending.clear();
      worker.terminate();
      if(this.worker===worker)this.worker=null;
    };
  }

  private request<T>(message: Record<string, unknown>, onProgress?: (progress: Progress) => void): Promise<T> {
    if (this.closed) return Promise.reject(new Error("The joint germline runtime is closed."));
    if (!this.worker) this.installWorker();
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject, onProgress });
      this.worker!.postMessage({ id, ...message });
    });
  }

  begin(referenceFasta: string, options: UnifiedOptions): Promise<unknown> {
    return this.request({ type: "begin", referenceFasta, options });
  }

  ingest(rows: readonly AirrScanRow[], assignments: Int32Array): Promise<unknown> {
    return this.request({
      type: "ingest",
      rows: rows.map((row) => ({ row: row.values, ordinal: row.ordinal, lineageId: assignments[row.ordinal] ?? 0 })),
    });
  }

  finish(onProgress?: (progress: Progress) => void): Promise<UnifiedDashboard> {
    return this.request<UnifiedDashboard>({ type: "finish" }, onProgress);
  }

  cancel() {
    if (this.closed) return;
    this.generation += 1;
    this.worker?.terminate();
    this.worker = null;
    const error = new DOMException("Personalized germline inference was cancelled.", "AbortError");
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
  }

  terminate() {
    this.closed = true;
    this.generation += 1;
    this.worker?.terminate();
    this.worker = null;
    const error = new Error("The joint germline runtime was closed.");
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
  }
}
