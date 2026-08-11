/// <reference lib="webworker" />

import { streamSequenceBatches, type SequenceBatch, type SequenceFormat } from "./sequence-stream";

interface StartRequest {
  type: "start";
  id: number;
  query: string | File;
  format: SequenceFormat;
  references: { V: string; D: string; J: string; C: string };
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  workers: number;
  countHint: number | null;
}

interface AckRequest {
  type: "ack";
  id: number;
  batch: number;
}

interface ComputeBatch {
  batch: number;
  header: string;
  body: ArrayBuffer;
  count: number;
  milliseconds: number;
}

interface ComputeSlot {
  index: number;
  worker: Worker;
  busy: boolean;
  resolve?: (batch: ComputeBatch) => void;
  reject?: (error: Error) => void;
}

const acknowledgements = new Map<string, () => void>();
const encoder = new TextEncoder();

function postProgress(id: number, stage: string, value: number) {
  self.postMessage({ id, type: "progress", stage, value });
}

function waitForAcknowledgement(id: number, batch: number): Promise<void> {
  return new Promise((resolve) => acknowledgements.set(`${id}:${batch}`, resolve));
}

async function compileCore(): Promise<WebAssembly.Module> {
  const url = `${import.meta.env.BASE_URL}swiftig.wasm`;
  try {
    return await WebAssembly.compileStreaming(fetch(url));
  } catch {
    const response = await fetch(url);
    if (!response.ok) throw new Error("The SwiftIG WebAssembly core could not be loaded.");
    return WebAssembly.compile(await response.arrayBuffer());
  }
}

function effectiveWorkerCount(request: StartRequest): number {
  const requested = Math.max(1, Math.min(16, Math.floor(request.workers || 1)));
  if (request.countHint !== null && request.countHint <= 3) return 1;
  if (request.countHint !== null && request.countHint <= 500) return Math.min(2, requested);
  return requested;
}

function batchSize(countHint: number | null): number {
  if (countHint !== null && countHint <= 10) return 1;
  if (countHint !== null && countHint <= 1000) return 100;
  return 1000;
}

async function initializeSlot(
  index: number,
  module: WebAssembly.Module,
  references: StartRequest["references"],
): Promise<ComputeSlot> {
  const worker = new Worker(new URL("./swiftig-compute-worker.ts", import.meta.url), { type: "module" });
  const slot: ComputeSlot = { index, worker, busy: false };
  return new Promise((resolve, reject) => {
    const fail = (error: Error) => {
      worker.terminate();
      reject(error);
    };
    worker.onerror = (event) => fail(new Error(event.message || `Compute worker ${index + 1} stopped unexpectedly.`));
    worker.onmessage = (event) => {
      const message = event.data as { type: "ready" | "error"; genes?: number; message?: string };
      if (message.type === "error") {
        fail(new Error(message.message || `Compute worker ${index + 1} could not initialize.`));
        return;
      }
      resolve(slot);
    };
    worker.postMessage({ type: "initialize", worker: index, module, references });
  });
}

function installComputeHandler(slot: ComputeSlot, onFatal: (error: Error) => void) {
  slot.worker.onerror = (event) => {
    const error = new Error(event.message || `Compute worker ${slot.index + 1} stopped unexpectedly.`);
    slot.reject?.(error);
    onFatal(error);
  };
  slot.worker.onmessage = (event) => {
    const message = event.data as ComputeBatch & { type: "batch" | "error"; message?: string };
    if (message.type === "error") {
      const error = new Error(message.message || `Compute worker ${slot.index + 1} failed.`);
      slot.reject?.(error);
      slot.resolve = undefined;
      slot.reject = undefined;
      onFatal(error);
      return;
    }
    slot.resolve?.(message);
    slot.resolve = undefined;
    slot.reject = undefined;
  };
}

function annotate(slot: ComputeSlot, batch: SequenceBatch, request: StartRequest): Promise<ComputeBatch> {
  const query = encoder.encode(batch.text);
  return new Promise((resolve, reject) => {
    slot.resolve = resolve;
    slot.reject = reject;
    slot.worker.postMessage({
      type: "annotate",
      batch: batch.index,
      query: query.buffer,
      count: batch.count,
      format: batch.format,
      minimumIdentity: request.minimumIdentity,
      strand: request.strand,
    }, [query.buffer]);
  });
}

async function handleRequest(request: StartRequest) {
  const slots: ComputeSlot[] = [];
  try {
    const workerCount = effectiveWorkerCount(request);
    postProgress(request.id, `Loading WebAssembly for ${workerCount} parallel worker${workerCount === 1 ? "" : "s"}`, 0.03);
    const module = await compileCore();
    postProgress(request.id, `Indexing germlines in ${workerCount} worker${workerCount === 1 ? "" : "s"}`, 0.08);
    slots.push(...await Promise.all(Array.from(
      { length: workerCount },
      (_, index) => initializeSlot(index, module, request.references),
    )));

    let fatalError: Error | null = null;
    let wakeAvailability: (() => void) | null = null;
    let resolveFinished!: () => void;
    const finished = new Promise<void>((resolve) => {
      resolveFinished = resolve;
    });
    const pending = new Map<number, { result: ComputeBatch; slot: ComputeSlot }>();
    const maximumOutstanding = workerCount * 2;
    let nextCommit = 0;
    let parsed = 0;
    let committed = 0;
    let inputDone = false;
    let flushing = false;
    let lastProgress = 0.14;
    let bytesRead = 0;
    let totalBytes = typeof request.query === "string" ? request.query.length : request.query.size;

    const fail = (error: Error) => {
      if (fatalError) return;
      fatalError = error;
      wakeAvailability?.();
      wakeAvailability = null;
      resolveFinished();
    };
    slots.forEach((slot) => installComputeHandler(slot, fail));

    const report = (stage: string) => {
      const inputFraction = totalBytes ? Math.min(1, bytesRead / totalBytes) : 0;
      const completion = inputDone
        ? (parsed ? committed / parsed : 0)
        : Math.min(0.96, inputFraction * 0.92 + (parsed ? committed / parsed : 0) * 0.08);
      lastProgress = Math.max(lastProgress, 0.14 + completion * 0.82);
      postProgress(request.id, stage, Math.min(0.96, lastProgress));
    };

    const release = (slot: ComputeSlot) => {
      slot.busy = false;
      wakeAvailability?.();
      wakeAvailability = null;
    };

    const maybeFinish = () => {
      if (inputDone && committed === parsed && slots.every((slot) => !slot.busy) && !pending.size) {
        resolveFinished();
      }
    };

    const flush = async () => {
      if (flushing || fatalError) return;
      flushing = true;
      try {
        while (pending.has(nextCommit)) {
          const { result, slot } = pending.get(nextCommit)!;
          pending.delete(nextCommit);
          committed += result.count;
          const acknowledgement = waitForAcknowledgement(request.id, nextCommit);
          self.postMessage({
            id: request.id,
            type: "batch",
            batch: nextCommit,
            header: result.header,
            body: result.body,
            count: result.count,
            processed: committed,
            total: inputDone ? parsed : null,
            milliseconds: result.milliseconds,
          }, [result.body]);
          await acknowledgement;
          nextCommit += 1;
          wakeAvailability?.();
          wakeAvailability = null;
          report(inputDone
            ? `Annotated ${committed.toLocaleString()} of ${parsed.toLocaleString()} sequences`
            : `Streamed + annotated ${committed.toLocaleString()} sequences`);
        }
      } catch (error) {
        fail(error instanceof Error ? error : new Error(String(error)));
      } finally {
        flushing = false;
        maybeFinish();
      }
    };

    const acquire = async (batchIndex: number): Promise<ComputeSlot> => {
      while (true) {
        if (fatalError) throw fatalError;
        const available = slots.find((slot) => !slot.busy);
        if (available && batchIndex - nextCommit < maximumOutstanding) {
          available.busy = true;
          return available;
        }
        await new Promise<void>((resolve) => { wakeAvailability = resolve; });
      }
    };

    postProgress(request.id, `Ready · ${workerCount} parallel worker${workerCount === 1 ? "" : "s"}`, 0.14);
    for await (const batch of streamSequenceBatches({
      source: request.query,
      format: request.format,
      batchSize: batchSize(request.countHint),
      onProgress: (state) => {
        bytesRead = state.bytesRead;
        totalBytes = state.totalBytes;
        if (state.recordsRead && state.recordsRead % 1000 === 0) {
          report(`Streaming input · ${state.recordsRead.toLocaleString()} sequences parsed`);
        }
      },
    })) {
      if (fatalError) throw fatalError;
      parsed += batch.count;
      const slot = await acquire(batch.index);
      void annotate(slot, batch, request).then((result) => {
        pending.set(batch.index, { result, slot });
        release(slot);
        void flush();
      }).catch(fail);
    }
    inputDone = true;
    bytesRead = totalBytes;
    report(`Input complete · finishing ${parsed.toLocaleString()} sequences`);
    maybeFinish();
    await finished;
    if (fatalError) throw fatalError;
    postProgress(request.id, "Finalizing the local AIRR index", 0.98);
    self.postMessage({
      id: request.id,
      type: "result",
      count: committed,
      total: parsed,
      workers: workerCount,
    });
  } catch (error) {
    self.postMessage({
      id: request.id,
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  } finally {
    slots.forEach((slot) => slot.worker.terminate());
  }
}

self.addEventListener("message", (event: MessageEvent<StartRequest | AckRequest>) => {
  const message = event.data;
  if (message.type === "ack") {
    const key = `${message.id}:${message.batch}`;
    acknowledgements.get(key)?.();
    acknowledgements.delete(key);
    return;
  }
  void handleRequest(message);
});

export {};
