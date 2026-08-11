import type { CompiledReferences } from "./reference-pack";

export interface RunOptions {
  query: string;
  format: 0 | 1 | 2 | 3;
  references: CompiledReferences;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  onProgress?: (stage: string, value: number) => void;
}

export interface RunResult {
  count: number;
  tsv: string;
}

let requestId = 0;

export function runSwiftIg(options: RunOptions): Promise<RunResult> {
  const id = ++requestId;
  const worker = new Worker(new URL("./swiftig-worker.ts", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    worker.onmessage = (event) => {
      const message = event.data as {
        id: number;
        type: "progress" | "result" | "error";
        stage?: string;
        value?: number;
        count?: number;
        tsv?: string;
        message?: string;
      };
      if (message.id !== id) return;
      if (message.type === "progress") {
        options.onProgress?.(message.stage ?? "Working", message.value ?? 0);
        return;
      }
      worker.terminate();
      if (message.type === "error") {
        reject(new Error(message.message ?? "SwiftIG failed."));
      } else {
        resolve({ count: message.count ?? 0, tsv: message.tsv ?? "" });
      }
    };
    worker.onerror = (event) => {
      worker.terminate();
      reject(new Error(event.message || "The SwiftIG worker stopped unexpectedly."));
    };
    worker.postMessage({
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
    });
  });
}
