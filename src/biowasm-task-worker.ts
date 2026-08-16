/// <reference lib="webworker" />

import { runCodonAwareKalign, runFastTree, runKalign } from "./biowasm-runtime";

type Request =
  | { id: number; type: "kalign"; fasta: string }
  | { id: number; type: "codonKalign"; fasta: string; frames?: number[] }
  | { id: number; type: "fastTree"; fasta: string; model: "gtr" | "jc"; fast: boolean };

const worker = self as unknown as DedicatedWorkerGlobalScope;

worker.onmessage = async (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    const result = request.type === "kalign"
      ? await runKalign(request.fasta)
      : request.type === "codonKalign"
        ? await runCodonAwareKalign(request.fasta, request.frames)
        : await runFastTree(request.fasta, request.model, request.fast);
    worker.postMessage({ id: request.id, result });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
