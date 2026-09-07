/// <reference lib="webworker" />

import {
  PersonalizedGermlineAccumulator,
} from "./personalized-germline.ts";

import {inferUnifiedGermline,DEFAULT_UNIFIED_OPTIONS,type UnifiedOptions} from "./unified-germline.ts";
let options:UnifiedOptions=DEFAULT_UNIFIED_OPTIONS;

interface InputRow {
  row: Record<string, string>;
  ordinal: number;
  lineageId: number;
}

type Request =
  | { id: number; type: "begin"; referenceFasta: string; options: UnifiedOptions }
  | { id: number; type: "ingest"; rows: InputRow[] }
  | { id: number; type: "finish" }
  | { id: number; type: "clear" };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let accumulator: PersonalizedGermlineAccumulator | null = null;

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.type === "begin") {
      options=request.options;
      accumulator = new PersonalizedGermlineAccumulator(request.referenceFasta);
      worker.postMessage({ id: request.id, result: { ready: true } });
      return;
    }
    if (!accumulator) throw new Error("Initialize personalized germline inference before adding AIRR evidence.");
    if (request.type === "ingest") {
      request.rows.forEach((item) => accumulator!.add(item.row, item.ordinal, item.lineageId));
      worker.postMessage({ id: request.id, result: { ingested: request.rows.length } });
      return;
    }
    if (request.type === "finish") {
      worker.postMessage({ id: request.id, progress: { processed: 0, total: 1, phase: "Testing inherited versus somatic sequence hypotheses" } });
      const result = inferUnifiedGermline(accumulator.researchSnapshot(),options,(processed, total) => {
        worker.postMessage({ id: request.id, progress: { processed, total: Math.max(1, total), phase: "Testing inherited versus somatic sequence hypotheses" } });
      });
      accumulator = null;
      worker.postMessage({ id: request.id, result });
      return;
    }
    accumulator = null;
    worker.postMessage({ id: request.id, result: { cleared: true } });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};
