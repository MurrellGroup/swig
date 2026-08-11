import type { CompiledReferences } from "./reference-pack";

export interface ResultBatch {
  header: string;
  body: Uint8Array;
  count: number;
  processed: number;
  total: number | null;
}

export interface RunOptions {
  query: string | File;
  format: 1 | 2 | 3;
  references: CompiledReferences;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  workers: number;
  countHint?: number | null;
  subsample?: { size: number; seed: number };
  onProgress?: (stage: string, value: number) => void;
  onBatch?: (batch: ResultBatch) => void | Promise<void>;
  signal?: AbortSignal;
}

export interface RunResult {
  count: number;
  total: number;
  inputRecords: number;
  workers: number;
}

let requestId = 0;

export function runSwiftIg(options: RunOptions): Promise<RunResult> {
  const id = ++requestId;
  const worker = new Worker(new URL("./swiftig-worker.ts", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    let finished = false;
    const fail = (error: Error) => {
      if (finished) return;
      finished = true;
      options.signal?.removeEventListener("abort", abort);
      worker.terminate();
      reject(error);
    };
    const abort = () => fail(new DOMException("Analysis cancelled.", "AbortError"));
    if (options.signal?.aborted) {
      abort();
      return;
    }
    options.signal?.addEventListener("abort", abort, { once: true });

    worker.onmessage = (event) => {
      const message = event.data as {
        id: number;
        type: "progress" | "batch" | "result" | "error";
        stage?: string;
        value?: number;
        batch?: number;
        header?: string;
        body?: ArrayBuffer;
        count?: number;
        processed?: number;
        total?: number | null;
        workers?: number;
        inputRecords?: number;
        message?: string;
      };
      if (message.id !== id || finished) return;
      if (message.type === "progress") {
        options.onProgress?.(message.stage ?? "Working", message.value ?? 0);
        return;
      }
      if (message.type === "batch") {
        Promise.resolve(options.onBatch?.({
          header: message.header ?? "",
          body: new Uint8Array(message.body ?? new ArrayBuffer(0)),
          count: message.count ?? 0,
          processed: message.processed ?? 0,
          total: message.total ?? null,
        })).then(() => {
          worker.postMessage({ type: "ack", id, batch: message.batch ?? 0 });
        }).catch((error) => fail(error instanceof Error ? error : new Error(String(error))));
        return;
      }
      if (message.type === "error") {
        fail(new Error(message.message ?? "SwiftIG failed."));
        return;
      }
      finished = true;
      options.signal?.removeEventListener("abort", abort);
      worker.terminate();
      resolve({
        count: message.count ?? 0,
        total: message.total ?? 0,
        inputRecords: message.inputRecords ?? message.total ?? 0,
        workers: message.workers ?? 1,
      });
    };
    worker.onerror = (event) => fail(new Error(event.message || "The SwiftIG worker stopped unexpectedly."));
    worker.postMessage({
      type: "start",
      id,
      query: options.query,
      format: options.format,
      references: {
        V: options.references.V,
        D: options.references.D,
        J: options.references.J,
        C: options.references.C,
      },
      minimumIdentity: options.minimumIdentity,
      strand: options.strand,
      workers: options.workers,
      countHint: options.countHint ?? null,
      subsample: options.subsample,
    });
  });
}
