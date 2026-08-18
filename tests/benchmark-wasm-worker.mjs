import fs from "node:fs";
import { parentPort, workerData } from "node:worker_threads";

import { WASI } from "@bjorn3/browser_wasi_shim";

const encoder = new TextEncoder();
const decoder = new TextDecoder();
let runtime;

function put(bytes) {
  const pointer = runtime.swig_alloc(bytes.byteLength);
  if (!pointer && bytes.byteLength) throw new Error("WASM allocation failed");
  new Uint8Array(runtime.memory.buffer, pointer, bytes.byteLength).set(bytes);
  return [pointer, bytes.byteLength];
}

function errorText() {
  return decoder.decode(new Uint8Array(
    runtime.memory.buffer,
    runtime.swig_error_ptr(),
    runtime.swig_error_len(),
  ));
}

async function initialize(references) {
  const wasmBytes = fs.readFileSync(new URL("../public/swiftig.wasm", import.meta.url));
  const wasi = new WASI([], [], []);
  const module = await WebAssembly.compile(wasmBytes);
  const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
  wasi.initialize(instance);
  runtime = instance.exports;
  const strategy = workerData?.strategy === "riat_mp" ? 1 :
    workerData?.strategy === "aer" ? 2 : 0;
  if (runtime.swig_set_assigner_strategy(strategy) !== 0) {
    throw new Error("SwiftIG rejected the benchmark assignment strategy");
  }
  if (workerData?.reference && runtime.swig_set_optimized_kernels) {
    runtime.swig_set_optimized_kernels(0);
  }
  if (workerData?.referenceOutput && runtime.swig_set_optimized_output) {
    runtime.swig_set_optimized_output(0);
  }
  const allocations = [references.V, references.D, references.J, references.C].map((value) => put(encoder.encode(value)));
  const genes = runtime.swig_init_database(...allocations.flat());
  allocations.forEach(([pointer]) => runtime.swig_free(pointer));
  if (genes < 0) throw new Error(errorText());
  return genes;
}

function annotate(message) {
  const [pointer, length] = put(new Uint8Array(message.query));
  const count = runtime.swig_annotate(pointer, length, 1, 600, 0);
  runtime.swig_free(pointer);
  if (count < 0) throw new Error(errorText());
  const output = new Uint8Array(
    runtime.memory.buffer,
    runtime.swig_result_ptr(),
    runtime.swig_result_len(),
  ).slice();
  parentPort.postMessage({ type: "batch", id: message.id, count, output: output.buffer }, [output.buffer]);
}

parentPort.on("message", (message) => {
  Promise.resolve(message.type === "initialize" ? initialize(message.references) : annotate(message)).then((genes) => {
    if (message.type === "initialize") parentPort.postMessage({ type: "ready", genes });
  }).catch((error) => parentPort.postMessage({ type: "error", id: message.id, message: error.message }));
});
