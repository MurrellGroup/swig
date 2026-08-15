/// <reference lib="webworker" />

import { SparseEvidenceAccumulator } from "./evidence.ts";
import { fitSparseAlleleModel } from "./model.ts";
import { buildReferenceAlleleGraph } from "./reference-graph.ts";
import type {
  AlleleRefinementOptions,
  RefinementInputRow,
  RefinementSegment,
  ReferenceAlleleGraph,
  SegmentRefinementResult,
} from "./types.ts";

type Request =
  | { id: number; type: "init"; totalRecords: number; options: AlleleRefinementOptions }
  | { id: number; type: "beginSegment"; segment: RefinementSegment; fasta: string }
  | { id: number; type: "ingest"; rows: RefinementInputRow[] }
  | { id: number; type: "finishSegment" }
  | { id: number; type: "clear" };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let totalRecords = 0;
let options: AlleleRefinementOptions | null = null;
let graph: ReferenceAlleleGraph | null = null;
let accumulator: SparseEvidenceAccumulator | null = null;

function transferResult(id: number, result: SegmentRefinementResult) {
  worker.postMessage({ id, result }, [
    result.mapNode.buffer,
    result.mapProbability.buffer,
    result.posteriorEntropy.buffer,
    result.localTopNode.buffer,
    result.localTopProbability.buffer,
  ]);
}

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.type === "init") {
      totalRecords = request.totalRecords;
      options = request.options;
      graph = null;
      accumulator = null;
      worker.postMessage({ id: request.id, result: { ready: true } });
      return;
    }
    if (!options || !totalRecords) throw new Error("Initialize repertoire-level allele refinement first.");
    if (request.type === "beginSegment") {
      worker.postMessage({ id: request.id, progress: { processed: 0, total: 1, phase: `Building ${request.segment} reference-neighbour graph` } });
      graph = buildReferenceAlleleGraph(request.fasta, request.segment, options.neighbourRadius);
      if (!graph.nodes.length) throw new Error(`The selected reference database has no usable ${request.segment} sequences.`);
      accumulator = new SparseEvidenceAccumulator(totalRecords, graph, options);
      worker.postMessage({ id: request.id, result: { nodes: graph.nodes.length, exactDuplicateLabels: graph.exactDuplicateLabels } });
      return;
    }
    if (request.type === "ingest") {
      if (!accumulator) throw new Error("Begin a reference segment before adding AIRR evidence.");
      let modeled = 0;
      request.rows.forEach((row) => { if (accumulator!.add(row)) modeled += 1; });
      worker.postMessage({ id: request.id, result: { ingested: request.rows.length, modeled } });
      return;
    }
    if (request.type === "finishSegment") {
      if (!accumulator || !graph) throw new Error("No allele-evidence matrix is available for fitting.");
      const matrix = accumulator.finish();
      accumulator = null;
      const modelLabel = options.model === "active-set" ? "hurdle active-set" : "Dirichlet";
      worker.postMessage({ id: request.id, progress: { processed: 0, total: Math.max(1, matrix.groupKeys.length), phase: `Fitting sparse ${graph.segment} ${modelLabel} models` } });
      const result = fitSparseAlleleModel(matrix, graph, options, totalRecords, (processed, total) => {
        worker.postMessage({ id: request.id, progress: { processed, total: Math.max(1, total), phase: `Fitting sparse ${graph!.segment} ${modelLabel} models` } });
      });
      graph = null;
      accumulator = null;
      transferResult(request.id, result);
      return;
    }
    graph = null;
    accumulator = null;
    options = null;
    totalRecords = 0;
    worker.postMessage({ id: request.id, result: { cleared: true } });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};
