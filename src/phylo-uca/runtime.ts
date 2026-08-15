import type { PhyloUcaInput, PhyloUcaProgress, PhyloUcaResult, PhyloUcaWorkerResponse } from "./types.ts";

let nextRequest = 0;

export function runPhyloUca(
  input: PhyloUcaInput,
  onProgress?: (progress: PhyloUcaProgress) => void,
  signal?: AbortSignal,
): Promise<PhyloUcaResult> {
  const id = ++nextRequest;
  const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    let finished = false;
    const stop = () => worker.terminate();
    const fail = (error: Error) => {
      if (finished) return;
      finished = true;
      signal?.removeEventListener("abort", abort);
      stop();
      reject(error);
    };
    const abort = () => fail(new DOMException("Phylogenetic UCA inference was cancelled.", "AbortError"));
    if (signal?.aborted) return abort();
    signal?.addEventListener("abort", abort, { once: true });
    worker.onerror = (event) => fail(new Error(event.message || "The phylogenetic UCA worker stopped."));
    worker.onmessage = (event: MessageEvent<PhyloUcaWorkerResponse>) => {
      const message = event.data;
      if (finished || message.id !== id) return;
      if (message.type === "progress") return onProgress?.(message.progress);
      if (message.type === "error") return fail(new Error(message.error));
      finished = true;
      signal?.removeEventListener("abort", abort);
      stop();
      resolve(message.result);
    };
    worker.postMessage({ id, input });
  });
}
