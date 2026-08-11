/// <reference lib="webworker" />

import { WASI } from "@bjorn3/browser_wasi_shim";

interface SwiftIgExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  swig_alloc: (size: number) => number;
  swig_free: (pointer: number) => void;
  swig_init_database: (
    vPointer: number, vSize: number, dPointer: number, dSize: number,
    jPointer: number, jSize: number, cPointer: number, cSize: number,
  ) => number;
  swig_annotate: (
    queryPointer: number, querySize: number, format: number,
    identityPerMille: number, strand: number,
  ) => number;
  swig_result_ptr: () => number;
  swig_result_len: () => number;
  swig_error_ptr: () => number;
  swig_error_len: () => number;
}

interface InitializeRequest {
  type: "initialize";
  worker: number;
  module: WebAssembly.Module;
  references: { V: string; D: string; J: string; C: string };
}

interface AnnotateRequest {
  type: "annotate";
  batch: number;
  query: ArrayBuffer;
  count: number;
  format: 1 | 2 | 3;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();
let runtime: SwiftIgExports | null = null;

function putBytes(exports: SwiftIgExports, bytes: Uint8Array): [number, number] {
  const pointer = exports.swig_alloc(bytes.byteLength);
  if (!pointer && bytes.byteLength) throw new Error("SwiftIG ran out of WebAssembly memory.");
  new Uint8Array(exports.memory.buffer, pointer, bytes.byteLength).set(bytes);
  return [pointer, bytes.byteLength];
}

function putText(exports: SwiftIgExports, value: string): [number, number] {
  return putBytes(exports, encoder.encode(value));
}

function readText(exports: SwiftIgExports, pointer: number, length: number): string {
  return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
}

function readError(exports: SwiftIgExports): string {
  return readText(exports, exports.swig_error_ptr(), exports.swig_error_len()) ||
    "SwiftIG could not complete the annotation.";
}

async function initialize(request: InitializeRequest) {
  const wasi = new WASI([], [], []);
  const instance = await WebAssembly.instantiate(request.module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  wasi.initialize(instance as WebAssembly.Instance & {
    exports: { memory: WebAssembly.Memory; _initialize?: () => unknown };
  });
  const exports = instance.exports as SwiftIgExports;
  const allocations = [
    putText(exports, request.references.V),
    putText(exports, request.references.D),
    putText(exports, request.references.J),
    putText(exports, request.references.C),
  ];
  try {
    const genes = exports.swig_init_database(
      allocations[0][0], allocations[0][1], allocations[1][0], allocations[1][1],
      allocations[2][0], allocations[2][1], allocations[3][0], allocations[3][1],
    );
    if (genes < 0) throw new Error(readError(exports));
    runtime = exports;
    self.postMessage({ type: "ready", worker: request.worker, genes });
  } finally {
    allocations.forEach(([pointer]) => exports.swig_free(pointer));
  }
}

function annotate(request: AnnotateRequest) {
  if (!runtime) throw new Error("The SwiftIG compute worker is not initialized.");
  const started = performance.now();
  const [queryPointer] = putBytes(runtime, new Uint8Array(request.query));
  let count: number;
  try {
    count = runtime.swig_annotate(
      queryPointer,
      request.query.byteLength,
      request.format,
      Math.round(request.minimumIdentity * 1000),
      request.strand,
    );
  } finally {
    runtime.swig_free(queryPointer);
  }
  if (count < 0) throw new Error(readError(runtime));
  if (count !== request.count) {
    throw new Error(`SwiftIG returned ${count} records for a ${request.count}-record input batch.`);
  }

  const pointer = runtime.swig_result_ptr();
  const length = runtime.swig_result_len();
  const wasmView = new Uint8Array(runtime.memory.buffer, pointer, length);
  const newline = wasmView.indexOf(10);
  if (newline < 0) throw new Error("SwiftIG returned an invalid AIRR table.");
  const header = decoder.decode(wasmView.subarray(0, newline)).replace(/\r$/, "");
  const body = wasmView.slice(newline + 1);
  self.postMessage({
    type: "batch",
    batch: request.batch,
    header,
    body: body.buffer,
    count,
    milliseconds: performance.now() - started,
  }, [body.buffer]);
}

self.addEventListener("message", (event: MessageEvent<InitializeRequest | AnnotateRequest>) => {
  const request = event.data;
  Promise.resolve().then(() => request.type === "initialize" ? initialize(request) : annotate(request)).catch((error) => {
    self.postMessage({
      type: "error",
      batch: request.type === "annotate" ? request.batch : -1,
      message: error instanceof Error ? error.message : String(error),
    });
  });
});

export {};
