import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import vm from "node:vm";
import { WASI } from "@bjorn3/browser_wasi_shim";

import {
  assertAlivibeMsaResult,
  decodeAlivibeMsaSequences,
  encodeAlivibeMsaSequences,
} from "../src/alivibe-msa-codec.ts";

const oracleSource = fs.readFileSync(new URL("../public/tools/nw.js", import.meta.url), "utf8");
const oracleContext = vm.createContext({});
vm.runInContext(`${oracleSource}\nthis.__refinedMSA = refinedMSA;`, oracleContext);
const oracle = oracleContext.__refinedMSA;

const wasmBytes = fs.readFileSync(new URL("../public/alivibe-msa.wasm", import.meta.url));
const wasi = new WASI([], [], []);
const module = await WebAssembly.compile(wasmBytes);
const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
wasi.initialize(instance);
const runtime = instance.exports;

function oracleAlign(sequences, iterations = 3) {
  return Array.from(oracle(Array.from(sequences), iterations), (sequence) => String(sequence));
}

function wasmAlign(sequences, iterations = 3, validate = iterations === 3) {
  const input = new Uint8Array(encodeAlivibeMsaSequences(sequences));
  const pointer = runtime.alivibe_msa_alloc(input.byteLength);
  assert.ok(pointer || input.byteLength === 0);
  try {
    new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
    const count = runtime.alivibe_msa_run(pointer, input.byteLength, iterations);
    if (count < 0) {
      const message = new TextDecoder().decode(new Uint8Array(
        runtime.memory.buffer,
        runtime.alivibe_msa_error_ptr(),
        runtime.alivibe_msa_error_len(),
      ));
      throw new Error(message || "Alivibe MSA failed.");
    }
    const output = new Uint8Array(
      runtime.memory.buffer,
      runtime.alivibe_msa_result_ptr(),
      runtime.alivibe_msa_result_len(),
    ).slice();
    const aligned = decodeAlivibeMsaSequences(output.buffer);
    assert.equal(count, aligned.length);
    if (validate) assertAlivibeMsaResult(sequences, aligned);
    return aligned;
  } finally {
    runtime.alivibe_msa_free(pointer);
  }
}

function mutate(source, substitutions, deletionStart = -1, deletionLength = 0, insertionAt = -1, insertion = "") {
  const characters = source.split("");
  substitutions.forEach(([position, character]) => { characters[position] = character; });
  if (deletionStart >= 0) characters.splice(deletionStart, deletionLength);
  if (insertionAt >= 0) characters.splice(insertionAt, 0, ...insertion);
  return characters.join("");
}

const longNt = "ATGGACTGGACCTGGAGGATCCTCTTCTTGGTGGCAGCAGCTACAGGTGTCCACTCCCAGGTGCAGCTGCAGGAGTCGGGCCC";
const fixtures = [
  ["ACGTACGT"],
  ["ACGTACGT", "ACGTTACGT", "ACGTTCGT", "ACGACGT"],
  [
    longNt,
    mutate(longNt, [[7, "A"], [42, "T"]], 31, 3),
    mutate(longNt, [[15, "C"], [68, "A"]], -1, 0, 55, "GGA"),
    mutate(longNt, [[22, "T"], [50, "C"], [73, "G"]]),
  ],
  [
    "AAAAAAAAAAACCCCCCCCCCCGGGGGGGGGGGTTTTTTTTTTTAAAAAAAAAAA",
    "AAAAAAAAAAACCCCCCCCCCCGGGGGTTTTTTTTTTTAAAAAAAAAAA",
    "AAAAAAAAAAACCCCCGGGGGGGGGGGTTTTTTTTTTTAAAAAAAAAAA",
    "AAAAATAAAAACCCCCCCCCCCGGGGGGGGGGGTTTTTTTTTTTAAAAATAAAAA",
  ],
  [
    "MKTIIALSYIFCLVFADYKDDDDK",
    "MKTIIALSYIFCLVFAEYKDDDDK",
    "MKTILSYIFCLVFADYKDDD*K",
    "MKTIIALSYXXFCLVFADYKDDDDK",
  ],
  [
    longNt.slice(8),
    mutate(longNt, [[1, "C"]], 0, 8),
    `${longNt.slice(0, 45)}TTAC${longNt.slice(45)}`,
    longNt,
  ],
];

test("Alivibe WASM port matches the pinned refinedMSA output exactly", () => {
  fixtures.forEach((sequences, fixture) => {
    for (let iterations = 0; iterations <= 3; iterations += 1) {
      assert.deepEqual(wasmAlign(sequences, iterations), oracleAlign(sequences, iterations), `fixture ${fixture + 1}, pass ${iterations}`);
    }
  });
});

test("Alivibe WASM preserves the original order-sensitive behavior", () => {
  const sequences = fixtures[2];
  for (const order of [
    [0, 1, 2, 3],
    [3, 2, 1, 0],
    [1, 3, 0, 2],
  ]) {
    const reordered = order.map((index) => sequences[index]);
    assert.deepEqual(wasmAlign(reordered), oracleAlign(reordered));
  }
});

test("Alivibe WASM matches deterministic mutation/indel fixtures", () => {
  let seed = 0x5a17c9e3;
  const random = () => {
    seed ^= seed << 13;
    seed ^= seed >>> 17;
    seed ^= seed << 5;
    return (seed >>> 0) / 0x1_0000_0000;
  };
  const alphabet = "ACGT";
  for (let fixture = 0; fixture < 12; fixture += 1) {
    let reference = "";
    const length = fixture < 4 ? 18 + fixture * 3 : 70 + fixture * 2;
    for (let i = 0; i < length; i += 1) reference += alphabet[Math.floor(random() * alphabet.length)];
    const sequences = [reference];
    for (let row = 1; row < 5; row += 1) {
      const characters = reference.split("");
      for (let change = 0; change < 4; change += 1) {
        const position = Math.floor(random() * characters.length);
        characters[position] = alphabet[Math.floor(random() * alphabet.length)];
      }
      if (row % 2 === 0) characters.splice(Math.floor(random() * characters.length), 2);
      else characters.splice(Math.floor(random() * characters.length), 0, "A", "C");
      sequences.push(characters.join(""));
    }
    assert.deepEqual(wasmAlign(sequences), oracleAlign(sequences), `random fixture ${fixture + 1}`);
  }
});

test("Alivibe MSA codec rejects malformed and non-rectangular results", () => {
  const encoded = encodeAlivibeMsaSequences(["ACGT", "ACCT"]);
  assert.deepEqual(decodeAlivibeMsaSequences(encoded), ["ACGT", "ACCT"]);
  assert.throws(() => assertAlivibeMsaResult(["ACGT", "ACCT"], ["A-GT", "ACCT-"]), /non-rectangular/i);
  assert.throws(() => assertAlivibeMsaResult(["ACGT"], ["ACGT", "ACCT"]), /wrong number/i);
  assert.throws(() => decodeAlivibeMsaSequences(new ArrayBuffer(8)), /invalid result/i);
});
