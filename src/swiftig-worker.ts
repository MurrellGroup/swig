/// <reference lib="webworker" />

import { WASI } from "@bjorn3/browser_wasi_shim";

interface WorkerRequest {
  id: number;
  query: string;
  format: 0 | 1 | 2 | 3;
  references: { V: string; D: string; J: string; C: string };
  minimumIdentity: number;
  strand: 0 | 1 | 2;
}

interface SwiftIgExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  _initialize: () => void;
  swig_alloc: (size: number) => number;
  swig_free: (pointer: number) => void;
  swig_init_database: (
    vPointer: number,
    vSize: number,
    dPointer: number,
    dSize: number,
    jPointer: number,
    jSize: number,
    cPointer: number,
    cSize: number,
  ) => number;
  swig_annotate: (
    queryPointer: number,
    querySize: number,
    format: number,
    identityPerMille: number,
    strand: number,
  ) => number;
  swig_result_ptr: () => number;
  swig_result_len: () => number;
  swig_error_ptr: () => number;
  swig_error_len: () => number;
}

let runtimePromise: Promise<SwiftIgExports> | null = null;
const encoder = new TextEncoder();
const decoder = new TextDecoder();

function progress(id: number, stage: string, value: number) {
  self.postMessage({ id, type: "progress", stage, value });
}

async function loadRuntime(): Promise<SwiftIgExports> {
  if (!runtimePromise) {
    runtimePromise = (async () => {
      const response = await fetch(`${import.meta.env.BASE_URL}swiftig.wasm`);
      if (!response.ok) throw new Error("The SwiftIG WebAssembly core could not be loaded.");
      const wasmModule = await WebAssembly.compile(await response.arrayBuffer());
      const wasi = new WASI([], [], []);
      const instance = await WebAssembly.instantiate(wasmModule, {
        wasi_snapshot_preview1: wasi.wasiImport,
      });
      wasi.initialize(instance as WebAssembly.Instance & {
        exports: { memory: WebAssembly.Memory; _initialize?: () => unknown };
      });
      return instance.exports as SwiftIgExports;
    })();
  }
  return runtimePromise;
}

function put(exports: SwiftIgExports, value: string): [number, number] {
  const bytes = encoder.encode(value);
  const pointer = exports.swig_alloc(bytes.length);
  if (!pointer && bytes.length) throw new Error("SwiftIG ran out of browser memory.");
  new Uint8Array(exports.memory.buffer, pointer, bytes.length).set(bytes);
  return [pointer, bytes.length];
}

function read(exports: SwiftIgExports, pointer: number, length: number): string {
  return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
}

function readError(exports: SwiftIgExports): string {
  return read(exports, exports.swig_error_ptr(), exports.swig_error_len()) ||
    "SwiftIG could not complete the annotation.";
}

self.onmessage = async (event: MessageEvent<WorkerRequest>) => {
  const request = event.data;
  try {
    progress(request.id, "Loading WebAssembly core", 0.08);
    const exports = await loadRuntime();
    progress(request.id, "Indexing V, D and J references", 0.22);
    const allocations = [
      put(exports, request.references.V),
      put(exports, request.references.D),
      put(exports, request.references.J),
      put(exports, request.references.C),
    ];
    const initialized = exports.swig_init_database(
      allocations[0][0], allocations[0][1],
      allocations[1][0], allocations[1][1],
      allocations[2][0], allocations[2][1],
      allocations[3][0], allocations[3][1],
    );
    allocations.forEach(([pointer]) => exports.swig_free(pointer));
    if (initialized < 0) throw new Error(readError(exports));
    progress(request.id, `Indexed ${initialized.toLocaleString()} germline alleles`, 0.46);
    const [queryPointer, querySize] = put(exports, request.query);
    const count = exports.swig_annotate(
      queryPointer,
      querySize,
      request.format,
      Math.round(request.minimumIdentity * 1000),
      request.strand,
    );
    exports.swig_free(queryPointer);
    if (count < 0) throw new Error(readError(exports));
    progress(request.id, `Annotated ${count.toLocaleString()} sequences`, 1);
    self.postMessage({
      id: request.id,
      type: "result",
      count,
      tsv: read(exports, exports.swig_result_ptr(), exports.swig_result_len()),
    });
  } catch (error) {
    self.postMessage({
      id: request.id,
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  }
};

export {};
