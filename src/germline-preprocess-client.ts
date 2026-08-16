import type {
  GermlineLocus,
  GermlinePreprocessReport,
  GermlineSegment,
  MetadataAllele,
} from "./germline-preprocess";

interface WorkerReply {
  report?: GermlinePreprocessReport;
  error?: string;
}

export function preprocessGermlinesInWorker(
  text: string,
  segment: GermlineSegment,
  templateTiers: MetadataAllele[][],
  allowedLoci: GermlineLocus[],
  signal?: AbortSignal,
): Promise<GermlinePreprocessReport> {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL("./germline-preprocess-worker.ts", import.meta.url), { type: "module" });
    let finished = false;
    const finish = (action: () => void) => {
      if (finished) return;
      finished = true;
      signal?.removeEventListener("abort", abort);
      worker.terminate();
      action();
    };
    const abort = () => finish(() => reject(new DOMException("Germline preprocessing was cancelled.", "AbortError")));
    if (signal?.aborted) {
      abort();
      return;
    }
    signal?.addEventListener("abort", abort, { once: true });
    worker.onmessage = (event: MessageEvent<WorkerReply>) => {
      finish(() => {
        if (event.data.report) resolve(event.data.report);
        else reject(new Error(event.data.error || "Germline preprocessing failed."));
      });
    };
    worker.onerror = (event) => {
      finish(() => reject(new Error(event.message || "Germline preprocessing worker failed.")));
    };
    worker.postMessage({ text, segment, templateTiers, allowedLoci });
  });
}
