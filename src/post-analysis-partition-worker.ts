/// <reference lib="webworker" />

import {
  runDenoisePartitionJob,
  runExactDedupJob,
  type DenoisePartitionJob,
  type ExactDedupJob,
} from "./post-analysis-core";

type Request =
  | { taskIndex: number; kind: "denoise"; job: DenoisePartitionJob }
  | { taskIndex: number; kind: "exact"; job: ExactDedupJob };

const worker = self as unknown as DedicatedWorkerGlobalScope;

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.kind === "denoise") {
      const result = runDenoisePartitionJob(request.job);
      worker.postMessage({ taskIndex: request.taskIndex, result }, [result.targets.buffer]);
    } else {
      const result = runExactDedupJob(request.job);
      worker.postMessage({ taskIndex: request.taskIndex, result }, [
        result.representatives.buffer,
        result.representativeOrdinals.buffer,
        result.representativeCounts.buffer,
      ]);
    }
  } catch (error) {
    worker.postMessage({ taskIndex: request.taskIndex, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
