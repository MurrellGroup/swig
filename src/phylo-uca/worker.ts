/// <reference lib="webworker" />

import { inferPhyloUca } from "./inference.ts";
import type { PhyloUcaWorkerRequest, PhyloUcaWorkerResponse } from "./types.ts";

const worker = self as unknown as DedicatedWorkerGlobalScope;

worker.onmessage = (event: MessageEvent<PhyloUcaWorkerRequest>) => {
  const { id, input } = event.data;
  void inferPhyloUca(input, (progress) => {
    worker.postMessage({ id, type: "progress", progress } satisfies PhyloUcaWorkerResponse);
  }).then((result) => {
    worker.postMessage({ id, type: "result", result } satisfies PhyloUcaWorkerResponse);
  }).catch((error) => {
    worker.postMessage({ id, type: "error", error: error instanceof Error ? error.message : String(error) } satisfies PhyloUcaWorkerResponse);
  });
};

export {};
