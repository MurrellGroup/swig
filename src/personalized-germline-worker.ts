/// <reference lib="webworker" />

import {
  PersonalizedGermlineAccumulator,
  type PersonalizedGermlineDashboard,
  type PersonalizedGermlineOptions,
} from "./personalized-germline.ts";

interface InputRow {
  row: Record<string, string>;
  ordinal: number;
  lineageId: number;
}

type Request =
  | { id: number; type: "begin"; referenceFasta: string; options: PersonalizedGermlineOptions }
  | { id: number; type: "ingest"; rows: InputRow[] }
  | { id: number; type: "finish" }
  | { id: number; type: "clear" };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let accumulator: PersonalizedGermlineAccumulator | null = null;

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.type === "begin") {
      accumulator = new PersonalizedGermlineAccumulator(request.referenceFasta, request.options);
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
      worker.postMessage({ id: request.id, progress: { processed: 0, total: 1, phase: "Fitting lineage-weighted personalized V sets" } });
      const result: PersonalizedGermlineDashboard = accumulator.finish((processed, total) => {
        worker.postMessage({ id: request.id, progress: { processed, total: Math.max(1, total), phase: "Fitting lineage-weighted personalized V sets" } });
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
