import fs from "node:fs";
import { performance } from "node:perf_hooks";

import { WASI } from "@bjorn3/browser_wasi_shim";

function argumentsMap(argv) {
  const result = new Map();
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (!key.startsWith("--")) continue;
    const value = argv[index + 1]?.startsWith("--") ? "true" : argv[++index] ?? "true";
    result.set(key.slice(2), value);
  }
  return result;
}

const args = argumentsMap(process.argv.slice(2));
const required = (name) => {
  const value = args.get(name);
  if (!value) throw new Error(`Missing --${name}`);
  return value;
};

const paths = {
  queries: required("queries"),
  calls: required("calls"),
  V: required("v"),
  D: required("d"),
  J: required("j"),
  wasm: args.get("wasm") ?? new URL("../public/swiftig.wasm", import.meta.url),
};
const split = args.get("split") ?? "all";
if (!["all", "development", "holdout"].includes(split)) throw new Error(`Unknown split: ${split}`);
const batchSize = Number(args.get("batch") ?? 500);
const minimumIdentity = Number(args.get("minimum-identity") ?? 600);
const strand = Number(args.get("strand") ?? 0);
const outputPath = args.get("out") ?? "";
const sampleMod = Number(args.get("sample-mod") ?? 1);
const sampleResidue = Number(args.get("sample-residue") ?? 0);
const tuning = args.has("tuning") ? JSON.parse(args.get("tuning")) : null;
const callingProfile = args.get("calling-profile") ?? "truth_optimized";
if (!["truth_optimized", "igblast_compatible", "igblast_balanced", "r_optimized"].includes(callingProfile)) throw new Error(`Unknown calling profile: ${callingProfile}`);
const target = args.get("target") ?? "truth";
if (!["truth", "igblast", "both"].includes(target)) throw new Error(`Unknown target: ${target}`);

function splitAtoms(value) {
  return String(value ?? "").split(/[,/]/).map((item) => item.trim()).filter(Boolean);
}

function fnv1a(value) {
  let hash = 0x811c9dc5;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193) >>> 0;
  }
  return hash;
}

function selectedForSplit(id) {
  if (sampleMod > 1 && fnv1a(id) % sampleMod !== sampleResidue) return false;
  if (split === "all") return true;
  const development = (fnv1a(id) & 1) === 0;
  return split === "development" ? development : !development;
}

function parseFasta(text) {
  const records = [];
  let id = "";
  let sequence = "";
  const commit = () => {
    if (id && selectedForSplit(id)) records.push({ id, sequence });
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      commit();
      id = line.slice(1).trim().split(/\s+/, 1)[0];
      sequence = "";
    } else sequence += line.trim();
  }
  commit();
  return records;
}

function parseTable(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const header = lines.shift().split("\t");
  const index = Object.fromEntries(header.map((name, position) => [name, position]));
  const rows = new Map();
  for (const line of lines) {
    const values = line.split("\t");
    const id = values[index.sequence_id];
    if (!selectedForSplit(id)) continue;
    rows.set(id, Object.fromEntries(header.map((name, position) => [name, values[position] ?? ""])));
  }
  return rows;
}

function parseAirr(text) {
  const lines = text.trimEnd().split("\n");
  const header = lines.shift().replace(/\r$/, "").split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.replace(/\r$/, "").split("\t");
    return Object.fromEntries(header.map((name, index) => [name, values[index] ?? ""]));
  });
}

function metric() {
  return {
    records: 0,
    top: 0,
    fair: 0,
    exactSet: 0,
    called: 0,
    multiple: 0,
    callCount: 0,
    truthPositive: 0,
    truthPositiveTop: 0,
    truthPositiveFair: 0,
    truthNegative: 0,
    truthNegativeCorrect: 0,
  };
}

function updateMetric(state, truthValue, predictionValue) {
  const truth = new Set(splitAtoms(truthValue));
  const prediction = splitAtoms(predictionValue);
  state.records += 1;
  if (prediction.length) {
    state.called += 1;
    state.callCount += prediction.length;
    if (prediction.length > 1) state.multiple += 1;
  }
  if (!truth.size) {
    state.truthNegative += 1;
    if (!prediction.length) {
      state.top += 1;
      state.fair += 1;
      state.exactSet += 1;
      state.truthNegativeCorrect += 1;
    }
    return;
  }
  state.truthPositive += 1;
  if (prediction.length && truth.has(prediction[0])) {
    state.top += 1;
    state.truthPositiveTop += 1;
  }
  if (prediction.length) {
    const overlap = new Set(prediction.filter((call) => truth.has(call))).size;
    const credit = overlap / prediction.length;
    state.fair += credit;
    state.truthPositiveFair += credit;
  }
  if (prediction.length === truth.size && prediction.every((call) => truth.has(call))) state.exactSet += 1;
}

function summarize(state) {
  const divide = (numerator, denominator = state.records) => denominator ? numerator / denominator : 0;
  return {
    records: state.records,
    top: divide(state.top),
    fair: divide(state.fair),
    exactSet: divide(state.exactSet),
    callRate: divide(state.called),
    multiCallRate: divide(state.multiple),
    meanCallsAll: divide(state.callCount),
    meanCallsWhenCalled: divide(state.callCount, state.called),
    truthPositiveRecords: state.truthPositive,
    truthPositiveTop: divide(state.truthPositiveTop, state.truthPositive),
    truthPositiveFair: divide(state.truthPositiveFair, state.truthPositive),
    truthNegativeRecords: state.truthNegative,
    truthNegativeSpecificity: divide(state.truthNegativeCorrect, state.truthNegative),
  };
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();
const wasmBytes = fs.readFileSync(paths.wasm);
const wasi = new WASI([], [], []);
const module = await WebAssembly.compile(wasmBytes);
const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
wasi.initialize(instance);
const runtime = instance.exports;

if (typeof runtime.swig_set_calling_profile === "function") {
  const result = runtime.swig_set_calling_profile(
    callingProfile === "truth_optimized" ? 0 : callingProfile === "r_optimized" ? 2 : 1,
  );
  if (result !== 0) throw new Error("WASM rejected the requested calling profile");
} else if (callingProfile !== "truth_optimized") {
  throw new Error("This WASM binary does not expose swig_set_calling_profile");
}

if (tuning) {
  if (typeof runtime.swig_set_tuning_options !== "function") throw new Error("This WASM binary does not expose swig_set_tuning_options");
  const defaults = {
    dMatch: 2, dMismatch: -5, dGapOpen: -5, dGapExtend: -1, topD: 6, minDMatch: 7,
    jMatch: 2, jMismatch: -5, jGapOpen: -11, jGapExtend: -1, topJ: 3, minJLength: 10,
  };
  const options = { ...defaults, ...tuning };
  runtime.swig_set_tuning_options(
    options.dMatch, options.dMismatch, options.dGapOpen, options.dGapExtend, options.topD, options.minDMatch,
    options.jMatch, options.jMismatch, options.jGapOpen, options.jGapExtend, options.topJ, options.minJLength,
  );
}

function putString(value) {
  const bytes = encoder.encode(value);
  const pointer = runtime.swig_alloc(bytes.byteLength);
  if (!pointer && bytes.byteLength) throw new Error("WASM allocation failed");
  new Uint8Array(runtime.memory.buffer, pointer, bytes.byteLength).set(bytes);
  return [pointer, bytes.byteLength];
}

function errorText() {
  return decoder.decode(new Uint8Array(runtime.memory.buffer, runtime.swig_error_ptr(), runtime.swig_error_len()));
}

const referenceText = ["V", "D", "J"].map((segment) => fs.readFileSync(paths[segment], "utf8"));
const referenceAllocations = [...referenceText, ""].map(putString);
const geneCount = runtime.swig_init_database(...referenceAllocations.flat());
referenceAllocations.forEach(([pointer]) => runtime.swig_free(pointer));
if (geneCount < 0) throw new Error(errorText());

const truthById = parseTable(fs.readFileSync(paths.calls, "utf8"));
const records = parseFasta(fs.readFileSync(paths.queries, "utf8"));
if (truthById.size !== records.length) throw new Error(`Split row mismatch: ${truthById.size} truth rows and ${records.length} FASTA records`);

const targetNames = target === "both" ? ["truth", "igblast"] : [target];
const balancedDFilter = callingProfile === "igblast_balanced"
  ? (await import("../src/balanced-calling-profile.ts")).applyBalancedDFilter
  : null;
const metricsByTarget = Object.fromEntries(targetNames.map((name) => [name, { V: metric(), D: metric(), J: metric() }]));
const predictionFields = [
  "sequence_id", "v_call", "d_call", "j_call", "d_score", "d_identity", "d_cigar",
  "d_sequence_alignment", "d_germline_alignment", "d_alternatives",
  "v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start",
  "d_germline_start", "d_germline_end", "junction_length",
];
const predictionLines = [predictionFields.join("\t")];
const started = performance.now();
for (let offset = 0; offset < records.length; offset += batchSize) {
  const batch = records.slice(offset, offset + batchSize);
  const fasta = batch.map((record) => `>${record.id}\n${record.sequence}\n`).join("");
  const [pointer, length] = putString(fasta);
  const count = runtime.swig_annotate(pointer, length, 1, minimumIdentity, strand);
  runtime.swig_free(pointer);
  if (count < 0) throw new Error(errorText());
  const resultView = new Uint8Array(runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len());
  const newline = resultView.indexOf(10);
  if (newline < 0) throw new Error("SwiftIG returned an invalid AIRR table");
  const header = decoder.decode(resultView.subarray(0, newline)).replace(/\r$/, "");
  let body = resultView.slice(newline + 1);
  if (balancedDFilter) body = balancedDFilter(header, body).body;
  const airr = `${header}\n${decoder.decode(body)}`;
  const predictions = parseAirr(airr);
  if (predictions.length !== batch.length) throw new Error(`Batch output mismatch at ${offset}`);
  for (const prediction of predictions) {
    const truth = truthById.get(prediction.sequence_id);
    if (!truth) throw new Error(`No truth row for ${prediction.sequence_id}`);
    for (const targetName of targetNames) {
      for (const [segment, lower] of [["V", "v"], ["D", "d"], ["J", "j"]]) {
        const targetColumn = targetName === "truth" ? `true_${lower}_call` : `igblast_${lower}_call`;
        updateMetric(metricsByTarget[targetName][segment], truth[targetColumn], prediction[`${lower}_call`]);
      }
    }
    if (outputPath) predictionLines.push(predictionFields.map((field) => prediction[field] ?? "").join("\t"));
  }
  if ((offset + batch.length) % 5000 === 0 || offset + batch.length === records.length) {
    process.stderr.write(`\r${(offset + batch.length).toLocaleString()} / ${records.length.toLocaleString()}`);
  }
}
process.stderr.write("\n");
const elapsedSeconds = (performance.now() - started) / 1000;
if (outputPath) fs.writeFileSync(outputPath, `${predictionLines.join("\n")}\n`);
console.log(JSON.stringify({
  split,
  records: records.length,
  references: geneCount,
  wasm: String(paths.wasm),
  minimumIdentity,
  strand,
  sampleMod,
  sampleResidue,
  tuning,
  callingProfile,
  target,
  elapsedSeconds,
  readsPerSecond: records.length / elapsedSeconds,
  metrics: target === "both"
    ? Object.fromEntries(Object.entries(metricsByTarget).map(([targetName, targetMetrics]) => [
      targetName,
      Object.fromEntries(Object.entries(targetMetrics).map(([segment, value]) => [segment, summarize(value)])),
    ]))
    : Object.fromEntries(Object.entries(metricsByTarget[target]).map(([segment, value]) => [segment, summarize(value)])),
}, null, 2));
