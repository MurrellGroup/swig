/// <reference lib="webworker" />

import {
  prepareReferenceMsa,
  runChmm,
  threadSequenceToMsa,
  type ChmmOptions,
  type ReferenceMsa,
} from "./post-analysis-core";

interface ChmmRow {
  ordinal: number;
  call: string;
  sequenceAlignment: string;
  germlineAlignment: string;
}

interface ChmmWorkerResult {
  ordinal: number;
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  status: "evaluated" | "low_dfr" | "missing_alignment" | "missing_reference" | "error";
  error?: string;
}

type Request =
  | { id: number; type: "init"; msa: string; options: ChmmOptions; minDfr: number }
  | { id: number; type: "batch"; rows: ChmmRow[] };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let msa: ReferenceMsa | null = null;
let options: ChmmOptions | null = null;
let minDfr = 1;
const cache = new Map<string, Omit<ChmmWorkerResult, "ordinal">>();

function mismatchCount(query: string, germline: string): number {
  let count = 0;
  for (let index = 0; index < Math.min(query.length, germline.length); index += 1) {
    const left = query[index]?.toUpperCase();
    const right = germline[index]?.toUpperCase();
    if (!left || !right || left === "-" || right === "-" || left === "N" || right === "N") continue;
    if (left !== right) count += 1;
  }
  return count;
}

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.type === "init") {
      msa = prepareReferenceMsa(request.msa);
      options = request.options;
      minDfr = request.minDfr;
      cache.clear();
      worker.postMessage({ id: request.id, result: { references: msa.names.length, length: msa.length } });
      return;
    }
    if (!msa || !options) throw new Error("CHMMAIRRa worker is not initialized.");
    const results: ChmmWorkerResult[] = [];
    for (const row of request.rows) {
      if (!row.call || !row.sequenceAlignment || !row.germlineAlignment) {
        results.push({ ordinal: row.ordinal, probability: Number.NaN, dfr: 0, startingReference: "", recombinations: [], status: "missing_alignment" });
        continue;
      }
      const dfr = mismatchCount(row.sequenceAlignment, row.germlineAlignment);
      if (dfr < minDfr) {
        results.push({ ordinal: row.ordinal, probability: Number.NaN, dfr, startingReference: "", recombinations: [], status: "low_dfr" });
        continue;
      }
      const key = `${row.call}\u0000${row.sequenceAlignment}\u0000${row.germlineAlignment}`;
      const cached = cache.get(key);
      if (cached) {
        results.push({ ordinal: row.ordinal, ...cached });
        continue;
      }
      try {
        const threaded = threadSequenceToMsa(row.sequenceAlignment, row.germlineAlignment, row.call, msa);
        const result = runChmm(msa, threaded, row.sequenceAlignment, row.germlineAlignment, options);
        const value: Omit<ChmmWorkerResult, "ordinal"> = { ...result, status: "evaluated" };
        cache.set(key, value);
        results.push({ ordinal: row.ordinal, ...value });
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        const status = message.includes("missing") ? "missing_reference" : "error";
        const value: Omit<ChmmWorkerResult, "ordinal"> = { probability: Number.NaN, dfr, startingReference: "", recombinations: [], status, error: message };
        cache.set(key, value);
        results.push({ ordinal: row.ordinal, ...value });
      }
    }
    worker.postMessage({ id: request.id, result: { results } });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
