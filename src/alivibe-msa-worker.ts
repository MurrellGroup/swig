/// <reference lib="webworker" />

import { WASI } from "@bjorn3/browser_wasi_shim";

interface AlivibeMsaExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  alivibe_msa_alloc: (size: number) => number;
  alivibe_msa_free: (pointer: number) => void;
  alivibe_msa_run: (pointer: number, length: number, iterations: number) => number;
  alivibe_msa_result_ptr: () => number;
  alivibe_msa_result_len: () => number;
  alivibe_msa_error_ptr: () => number;
  alivibe_msa_error_len: () => number;
}

interface AlignRequest {
  type: "align";
  input: ArrayBuffer;
  iterations: number;
}

const decoder = new TextDecoder();

function readText(runtime: AlivibeMsaExports, pointer: number, length: number): string {
  if (!pointer || !length) return "";
  return decoder.decode(new Uint8Array(runtime.memory.buffer, pointer, length));
}

async function align(request: AlignRequest): Promise<void> {
  const url = `${import.meta.env.BASE_URL}alivibe-msa.wasm`;
  let module: WebAssembly.Module;
  try {
    module = await WebAssembly.compileStreaming(fetch(url));
  } catch {
    const response = await fetch(url);
    if (!response.ok) throw new Error("The Alivibe MSA WebAssembly core could not be loaded.");
    module = await WebAssembly.compile(await response.arrayBuffer());
  }
  const wasi = new WASI([], [], []);
  const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
  wasi.initialize(instance as WebAssembly.Instance & {
    exports: { memory: WebAssembly.Memory; _initialize?: () => unknown };
  });
  const runtime = instance.exports as AlivibeMsaExports;
  const input = new Uint8Array(request.input);
  const pointer = runtime.alivibe_msa_alloc(input.byteLength);
  if (!pointer && input.byteLength) throw new Error("Alivibe MSA ran out of WebAssembly memory.");
  try {
    new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
    const count = runtime.alivibe_msa_run(pointer, input.byteLength, request.iterations);
    if (count < 0) {
      throw new Error(readText(
        runtime,
        runtime.alivibe_msa_error_ptr(),
        runtime.alivibe_msa_error_len(),
      ) || "Alivibe MSA failed.");
    }
    const result = new Uint8Array(
      runtime.memory.buffer,
      runtime.alivibe_msa_result_ptr(),
      runtime.alivibe_msa_result_len(),
    ).slice();
    self.postMessage({ type: "result", result: result.buffer, count }, [result.buffer]);
  } finally {
    runtime.alivibe_msa_free(pointer);
  }
}

self.addEventListener("message", (event: MessageEvent<AlignRequest>) => {
  void align(event.data).catch((error) => {
    self.postMessage({ type: "error", message: error instanceof Error ? error.message : String(error) });
  });
});

export {};
