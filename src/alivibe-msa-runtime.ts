import {
  assertAlivibeMsaResult,
  decodeAlivibeMsaSequences,
  encodeAlivibeMsaSequences,
} from "./alivibe-msa-codec";

export interface AlivibeMsaJob {
  result: Promise<string[]>;
  cancel: () => void;
}

export type AlivibeMsaScoringMode = "literal" | "nucleotide" | "amino-acid";

/**
 * Run one Alivibe-compatible MSA worker as an abortable Swig task. The worker
 * owns all mutable alignment state, so aborting before resolution cannot
 * install a partial result in the lineage workbench.
 */
export async function runAlivibeMsaTask(
  sequences: readonly string[],
  signal?: AbortSignal,
  iterations = 3,
  scoringMode: AlivibeMsaScoringMode = "nucleotide",
): Promise<string[]> {
  const job = createAlivibeMsaJob(sequences, iterations, scoringMode);
  const abort = () => job.cancel();
  if (signal?.aborted) abort();
  else signal?.addEventListener("abort", abort, { once: true });
  try {
    return await job.result;
  } finally {
    signal?.removeEventListener("abort", abort);
  }
}

export function createAlivibeMsaJob(
  sequences: readonly string[],
  iterations = 3,
  scoringMode: AlivibeMsaScoringMode = "nucleotide",
): AlivibeMsaJob {
  const inputSequences = sequences.map((sequence) => String(sequence));
  const input = encodeAlivibeMsaSequences(inputSequences);
  const worker = new Worker(new URL("./alivibe-msa-worker.ts", import.meta.url), { type: "module" });
  let settled = false;
  let rejectJob: ((reason: unknown) => void) | null = null;
  const result = new Promise<string[]>((resolve, reject) => {
    rejectJob = reject;
    worker.onmessage = (event) => {
      if (settled) return;
      const message = event.data as { type: "result" | "error"; result?: ArrayBuffer; message?: string };
      if (message.type === "error") {
        settled = true;
        worker.terminate();
        reject(new Error(message.message || "Alivibe MSA failed."));
        return;
      }
      try {
        const aligned = decodeAlivibeMsaSequences(message.result ?? new ArrayBuffer(0));
        assertAlivibeMsaResult(inputSequences, aligned);
        settled = true;
        worker.terminate();
        resolve(aligned);
      } catch (error) {
        settled = true;
        worker.terminate();
        reject(error);
      }
    };
    worker.onerror = (event) => {
      if (settled) return;
      settled = true;
      worker.terminate();
      reject(new Error(event.message || "The Alivibe MSA worker stopped unexpectedly."));
    };
    worker.postMessage({ type: "align", input, iterations, scoringMode }, [input]);
  });
  return {
    result,
    cancel: () => {
      if (settled) return;
      settled = true;
      worker.terminate();
      rejectJob?.(new DOMException("MSA alignment cancelled.", "AbortError"));
    },
  };
}
