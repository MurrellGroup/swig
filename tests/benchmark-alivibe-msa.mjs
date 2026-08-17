import assert from "node:assert/strict";
import fs from "node:fs";
import vm from "node:vm";
import { performance } from "node:perf_hooks";
import { WASI } from "@bjorn3/browser_wasi_shim";

import { decodeAlivibeMsaSequences, encodeAlivibeMsaSequences } from "../src/alivibe-msa-codec.ts";

function option(name, fallback) {
  const index = process.argv.indexOf(`--${name}`);
  return index >= 0 ? Number(process.argv[index + 1]) : fallback;
}

const rows = Math.max(2, Math.floor(option("rows", 80)));
const length = Math.max(40, Math.floor(option("length", 360)));
const repeats = Math.max(1, Math.floor(option("repeats", 3)));

const source = fs.readFileSync(new URL("../public/tools/nw.js", import.meta.url), "utf8");
const context = vm.createContext({});
vm.runInContext(`${source}\nthis.__refinedMSA = refinedMSA;`, context);
const refinedMsa = context.__refinedMSA;

const bytes = fs.readFileSync(new URL("../public/alivibe-msa.wasm", import.meta.url));
const wasi = new WASI([], [], []);
const module = await WebAssembly.compile(bytes);
const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
wasi.initialize(instance);
const runtime = instance.exports;

function wasmAlign(sequences) {
  const input = new Uint8Array(encodeAlivibeMsaSequences(sequences));
  const pointer = runtime.alivibe_msa_alloc(input.byteLength);
  try {
    new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
    if (runtime.alivibe_msa_run(pointer, input.byteLength, 3) < 0) throw new Error("WASM MSA failed.");
    return decodeAlivibeMsaSequences(new Uint8Array(
      runtime.memory.buffer,
      runtime.alivibe_msa_result_ptr(),
      runtime.alivibe_msa_result_len(),
    ).slice().buffer);
  } finally {
    runtime.alivibe_msa_free(pointer);
  }
}

let randomState = 0x6d2b79f5;
function random() {
  randomState = Math.imul(randomState ^ (randomState >>> 15), 1 | randomState);
  randomState ^= randomState + Math.imul(randomState ^ (randomState >>> 7), 61 | randomState);
  return ((randomState ^ (randomState >>> 14)) >>> 0) / 0x1_0000_0000;
}

const alphabet = "ACGT";
let reference = "";
for (let site = 0; site < length; site += 1) reference += alphabet[Math.floor(random() * alphabet.length)];
const sequences = [reference];
for (let row = 1; row < rows; row += 1) {
  const sequence = reference.split("");
  const substitutions = 3 + Math.floor(random() * 9);
  for (let change = 0; change < substitutions; change += 1) {
    sequence[Math.floor(random() * sequence.length)] = alphabet[Math.floor(random() * alphabet.length)];
  }
  if (row % 5 === 0) sequence.splice(80 + Math.floor(random() * Math.max(1, length - 160)), 3);
  if (row % 7 === 0) sequence.splice(80 + Math.floor(random() * Math.max(1, length - 160)), 0, "A", "C", "G");
  sequences.push(sequence.join(""));
}

// Warm both implementations outside the measured block.
const warm = sequences.slice(0, Math.min(8, sequences.length));
wasmAlign(warm);
Array.from(refinedMsa(Array.from(warm), 3));

const wasmTimes = [];
const javascriptTimes = [];
for (let repeat = 0; repeat < repeats; repeat += 1) {
  let started = performance.now();
  const wasm = wasmAlign(sequences);
  wasmTimes.push(performance.now() - started);
  started = performance.now();
  const javascript = Array.from(refinedMsa(Array.from(sequences), 3), String);
  javascriptTimes.push(performance.now() - started);
  assert.deepEqual(wasm, javascript);
}

const median = (values) => [...values].sort((left, right) => left - right)[Math.floor(values.length / 2)];
const wasmMedian = median(wasmTimes);
const javascriptMedian = median(javascriptTimes);
console.log(JSON.stringify({
  rows,
  inputLength: length,
  outputColumns: wasmAlign(sequences)[0].length,
  repeats,
  wasmMedianMs: Number(wasmMedian.toFixed(3)),
  javascriptMedianMs: Number(javascriptMedian.toFixed(3)),
  speedup: Number((javascriptMedian / wasmMedian).toFixed(2)),
}, null, 2));
