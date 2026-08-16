import type { FastTreeRun } from "./biowasm-runtime";

interface Pending {
  resolve: (value: unknown) => void;
  reject: (error: Error) => void;
  removeAbort?: () => void;
}

class BiowasmTaskRuntime {
  private worker!: Worker;
  private readonly pending = new Map<number, Pending>();
  private nextId = 1;
  private generation = 0;

  constructor() {
    this.installWorker();
  }

  private installWorker() {
    const generation = this.generation;
    this.worker = new Worker(new URL("./biowasm-task-worker.ts", import.meta.url), { type: "module" });
    this.worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string }>) => {
      if (generation !== this.generation) return;
      const pending = this.pending.get(event.data.id);
      if (!pending) return;
      this.pending.delete(event.data.id);
      pending.removeAbort?.();
      if (event.data.error) pending.reject(new Error(event.data.error));
      else pending.resolve(event.data.result);
    };
    this.worker.onerror = (event) => {
      if (generation !== this.generation) return;
      this.restart(new Error(event.message || "The Kalign/FastTree worker stopped unexpectedly."));
    };
  }

  private restart(error: Error) {
    this.generation += 1;
    this.worker.terminate();
    this.pending.forEach((pending) => {
      pending.removeAbort?.();
      pending.reject(error);
    });
    this.pending.clear();
    this.installWorker();
  }

  request<T>(message: Record<string, unknown>, signal?: AbortSignal): Promise<T> {
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      const abort = () => {
        if (this.pending.has(id)) this.restart(new DOMException("Kalign/FastTree was cancelled.", "AbortError"));
      };
      if (signal?.aborted) {
        reject(new DOMException("Kalign/FastTree was cancelled.", "AbortError"));
        return;
      }
      signal?.addEventListener("abort", abort, { once: true });
      this.pending.set(id, {
        resolve: resolve as (value: unknown) => void,
        reject,
        removeAbort: signal ? () => signal.removeEventListener("abort", abort) : undefined,
      });
      this.worker.postMessage({ id, ...message });
    });
  }
}

const runtime = new BiowasmTaskRuntime();

export function runKalignTask(fasta: string, signal?: AbortSignal): Promise<string> {
  return runtime.request({ type: "kalign", fasta }, signal);
}

export function runCodonAwareKalignTask(fasta: string, frames?: number[], signal?: AbortSignal): Promise<string> {
  return runtime.request({ type: "codonKalign", fasta, frames }, signal);
}

export function runFastTreeTask(alignedFasta: string, model: "gtr" | "jc" = "gtr", fast = false, signal?: AbortSignal): Promise<FastTreeRun> {
  return runtime.request({ type: "fastTree", fasta: alignedFasta, model, fast }, signal);
}
