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
): Promise<GermlinePreprocessReport> {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL("./germline-preprocess-worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (event: MessageEvent<WorkerReply>) => {
      worker.terminate();
      if (event.data.report) resolve(event.data.report);
      else reject(new Error(event.data.error || "Germline preprocessing failed."));
    };
    worker.onerror = (event) => {
      worker.terminate();
      reject(new Error(event.message || "Germline preprocessing worker failed."));
    };
    worker.postMessage({ text, segment, templateTiers, allowedLoci });
  });
}
