/// <reference lib="webworker" />

import { WASI } from "@bjorn3/browser_wasi_shim";
import { applyBalancedDFilter, reconcileBalancedDoubleD } from "./balanced-calling-profile";
import type { AssignerStrategy, CallingProfile, DoubleDScreenOptions } from "./swiftig-runtime";

interface SwiftIgExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  swig_alloc: (size: number) => number;
  swig_free: (pointer: number) => void;
  swig_set_calling_profile: (profile: number) => number;
  swig_set_assigner_strategy: (strategy: number) => number;
  swig_init_database: (
    vPointer: number, vSize: number, dPointer: number, dSize: number,
    jPointer: number, jSize: number, cPointer: number, cSize: number,
  ) => number;
  swig_annotate: (
    queryPointer: number, querySize: number, format: number,
    identityPerMille: number, strand: number,
  ) => number;
  swig_annotate_double_d: (
    queryPointer: number, querySize: number, format: number,
    identityPerMille: number, strand: number, mode: number,
    minimumVjSpan: number, seedLength: number, pseudoTrim: number,
    maximumPseudoMismatches: number, minimumScoreGain: number,
  ) => number;
  swig_result_ptr: () => number;
  swig_result_len: () => number;
  swig_double_d_result_ptr: () => number;
  swig_double_d_result_len: () => number;
  swig_double_d_count: () => number;
  swig_error_ptr: () => number;
  swig_error_len: () => number;
}

interface InitializeRequest {
  type: "initialize";
  worker: number;
  module: WebAssembly.Module;
  references: { V: string; D: string; J: string; C: string };
  callingProfile: CallingProfile;
  assignerStrategy: AssignerStrategy;
}

interface AnnotateRequest {
  type: "annotate";
  batch: number;
  query: ArrayBuffer;
  count: number;
  format: 1 | 2 | 3;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  doubleD?: DoubleDScreenOptions;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();
let runtime: SwiftIgExports | null = null;
let activeCallingProfile: CallingProfile = "truth_optimized";

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
  const strategy = request.assignerStrategy === "riat_mp" ? 1 : request.assignerStrategy === "aer" ? 2 : 0;
  if (exports.swig_set_assigner_strategy(strategy) !== 0) {
    throw new Error("SwiftIG rejected the selected assignment strategy.");
  }
  if (exports.swig_set_calling_profile(
    request.callingProfile === "truth_optimized" ? 0 : 1,
  ) !== 0) {
    throw new Error("SwiftIG rejected the selected calling profile.");
  }
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
    activeCallingProfile = request.callingProfile;
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
    const doubleD = request.doubleD;
    count = doubleD && doubleD.mode !== "off"
      ? runtime.swig_annotate_double_d(
        queryPointer,
        request.query.byteLength,
        request.format,
        Math.round(request.minimumIdentity * 1000),
        request.strand,
        doubleD.mode === "all" ? 1 : 2,
        Math.round(doubleD.minimumVjSpan),
        Math.round(doubleD.seedLength),
        Math.round(doubleD.pseudoTrim),
        Math.round(doubleD.maximumPseudoMismatches),
        Math.round(doubleD.minimumScoreGain),
      )
      : runtime.swig_annotate(
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
  let body = wasmView.slice(newline + 1);
  const balanced = activeCallingProfile === "igblast_balanced"
    ? applyBalancedDFilter(header, body)
    : null;
  if (balanced) body = balanced.body;
  const message: Record<string, unknown> = {
    type: "batch",
    batch: request.batch,
    header,
    body: body.buffer,
    count,
    milliseconds: performance.now() - started,
  };
  const transfers: ArrayBuffer[] = [body.buffer];
  if (request.doubleD && request.doubleD.mode !== "off") {
    const doubleDView = new Uint8Array(
      runtime.memory.buffer,
      runtime.swig_double_d_result_ptr(),
      runtime.swig_double_d_result_len(),
    );
    const doubleDNewline = doubleDView.indexOf(10);
    if (doubleDNewline < 0) throw new Error("SwiftIG returned an invalid double-D evidence table.");
    const doubleDHeader = decoder.decode(doubleDView.subarray(0, doubleDNewline)).replace(/\r$/, "");
    let doubleDBody = doubleDView.slice(doubleDNewline + 1);
    if (balanced) {
      doubleDBody = reconcileBalancedDoubleD(
        doubleDHeader,
        doubleDBody,
        balanced.suppressedSequenceIds,
      );
    }
    message.doubleDHeader = doubleDHeader;
    message.doubleDBody = doubleDBody.buffer;
    message.doubleDCount = runtime.swig_double_d_count();
    transfers.push(doubleDBody.buffer);
  }
  self.postMessage(message, transfers);
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
