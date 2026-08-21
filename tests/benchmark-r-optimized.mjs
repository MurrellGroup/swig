import fs from "node:fs";
import path from "node:path";
import zlib from "node:zlib";
import { performance } from "node:perf_hooks";

import { WASI } from "@bjorn3/browser_wasi_shim";

function parseArguments(argv) {
  const values = new Map();
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (!key.startsWith("--")) continue;
    const next = argv[index + 1];
    values.set(key.slice(2), next && !next.startsWith("--") ? argv[++index] : "true");
  }
  return values;
}

const args = parseArguments(process.argv.slice(2));
const required = (name) => {
  const value = args.get(name);
  if (!value) throw new Error(`Missing --${name}`);
  return path.resolve(value);
};
const paths = {
  truth: required("truth"),
  V: required("v"),
  D: required("d"),
  J: required("j"),
  wasm: path.resolve(args.get("wasm") ?? new URL("../public/swiftig.wasm", import.meta.url).pathname),
};
const split = args.get("split") ?? "train";
if (!["train", "validation", "test", "development", "all"].includes(split)) {
  throw new Error("--split must be train, validation, test, development, or all");
}
const profile = args.get("profile") ?? "truth_optimized";
const profileNumbers = { truth_optimized: 0, igblast_compatible: 1, igblast_balanced: 1, r_optimized: 2 };
if (!(profile in profileNumbers)) throw new Error(`Unknown profile: ${profile}`);
const batchSize = Number(args.get("batch") ?? 500);
const minimumIdentity = Number(args.get("minimum-identity") ?? 600);
const strand = Number(args.get("strand") ?? 1);
const outputPath = args.get("out") ? path.resolve(args.get("out")) : "";
const failuresPath = args.get("failures") ? path.resolve(args.get("failures")) : "";
const limit = Number(args.get("limit") ?? 0);
const sampleMod = Number(args.get("sample-mod") ?? 1);
const sampleResidue = Number(args.get("sample-residue") ?? 0);
const idsPath = args.get("ids") ? path.resolve(args.get("ids")) : "";
const selectedIds = idsPath
  ? new Set(readText(idsPath).split(/\s+/).map((value) => value.trim()).filter(Boolean))
  : null;
const excludedIds = new Set((args.get("exclude-id") ?? "").split(",").map((value) => value.trim()).filter(Boolean));
const tuning = args.has("tuning") ? JSON.parse(args.get("tuning")) : null;

function readText(file) {
  const bytes = fs.readFileSync(file);
  return file.endsWith(".gz") ? zlib.gunzipSync(bytes).toString("utf8") : bytes.toString("utf8");
}

function fnv1a(value) {
  let hash = 0x811c9dc5;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193) >>> 0;
  }
  return hash;
}

function selected(id) {
  if (selectedIds && !selectedIds.has(id)) return false;
  if (excludedIds.has(id)) return false;
  if (sampleMod > 1 && fnv1a(id) % sampleMod !== sampleResidue) return false;
  const residue = fnv1a(id) % 10;
  if (split === "all") return true;
  if (split === "train") return residue <= 5;
  if (split === "validation") return residue === 6 || residue === 7;
  if (split === "test") return residue >= 8;
  return residue <= 7;
}

function parseTable(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(header.map((name, index) => [name, values[index] ?? ""]));
  });
}

function parseFasta(text) {
  const records = [];
  let name = "";
  let sequence = "";
  const commit = () => {
    if (name) records.push({ name, sequence: sequence.toUpperCase().replace(/\s/g, "") });
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      commit();
      name = line.slice(1).trim().split(/\s+/, 1)[0];
      sequence = "";
    } else {
      sequence += line.trim();
    }
  }
  commit();
  return records;
}

function parseHeaderTruthFasta(text) {
  const rows = [];
  let header = "";
  let sequence = "";
  const commit = () => {
    if (!header) return;
    const tokens = header.trim().split(/\s+/);
    const fields = Object.fromEntries(tokens.slice(1).map((token) => {
      const separator = token.indexOf("=");
      return separator < 0 ? [token, ""] : [token.slice(0, separator), token.slice(separator + 1)];
    }));
    const row = {
      sequence_id: tokens[0],
      sequence: sequence.toUpperCase().replace(/\s/g, ""),
      v_call: fields.gold_v ?? "",
      d_call: fields.gold_d === "." ? "" : fields.gold_d ?? "",
      j_call: fields.gold_j ?? "",
      v_sequence_start: "", v_sequence_end: "",
      d_sequence_start: "", d_sequence_end: "",
      j_sequence_start: "", j_sequence_end: "",
    };
    const span = (fields.gold_span ?? "").match(/^(\d+)-(\d+)$/);
    const segment = fields.locus ?? fields.missed ?? "";
    if (span && ["V", "D", "J"].includes(segment)) {
      const lower = segment.toLowerCase();
      row[`${lower}_sequence_start`] = span[1];
      row[`${lower}_sequence_end`] = span[2];
    }
    rows.push(row);
  };
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.replace(/^\u200b/, "");
    if (line.startsWith(">")) {
      commit();
      header = line.slice(1);
      sequence = "";
    } else {
      sequence += line.trim();
    }
  }
  commit();
  return rows;
}

for (const file of (args.get("exclude-fasta") ?? "").split("|").map((value) => value.trim()).filter(Boolean)) {
  for (const row of parseHeaderTruthFasta(readText(path.resolve(file)))) {
    excludedIds.add(row.sequence_id.replace(/\|d=.*$/, ""));
  }
}

function callAtoms(value) {
  return String(value ?? "").split(/[,/]/).map((atom) => atom.trim()).filter(Boolean);
}

const referenceText = Object.fromEntries(["V", "D", "J"].map((segment) => [segment, readText(paths[segment])]));
const references = Object.fromEntries(Object.entries(referenceText).map(([segment, text]) => [segment, parseFasta(text)]));
const equivalence = Object.fromEntries(Object.entries(references).map(([segment, records]) => {
  const byName = new Map();
  const canonicalBySequence = new Map();
  for (const record of records) {
    if (!canonicalBySequence.has(record.sequence)) canonicalBySequence.set(record.sequence, record.name);
    byName.set(record.name, canonicalBySequence.get(record.sequence));
  }
  return [segment, byName];
}));

function canonicalAtom(segment, atom) {
  return equivalence[segment].get(atom) ?? atom;
}

function probabilityMap(segment, value) {
  const atoms = [...new Set(callAtoms(value).map((atom) => canonicalAtom(segment, atom)))];
  if (!atoms.length) return new Map([["__NONE__", 1]]);
  const probability = 1 / atoms.length;
  return new Map(atoms.map((atom) => [atom, probability]));
}

function scoreCall(segment, truthValue, predictionValue) {
  const truthAtoms = callAtoms(truthValue);
  if (truthAtoms.length > 1) return null;
  const truth = truthAtoms.length ? canonicalAtom(segment, truthAtoms[0]) : "__NONE__";
  const prediction = probabilityMap(segment, predictionValue);
  const truthProbability = prediction.get(truth) ?? 0;
  let squared = (truthProbability - 1) ** 2;
  for (const [atom, probability] of prediction) {
    if (atom !== truth) squared += probability ** 2;
  }
  return {
    brier: squared / 2,
    truthProbability,
    containsTruth: truthProbability > 0,
    exactSingleton: prediction.size === 1 && truthProbability === 1,
    predictedCount: prediction.size,
    truthAbsent: truth === "__NONE__",
    predictedAbsent: prediction.has("__NONE__"),
  };
}

function emptyCallMetric() {
  return {
    records: 0, brier: 0, truthProbability: 0, containsTruth: 0, exactSingleton: 0,
    predictedCount: 0, ambiguous: 0, truthAbsent: 0, trueNegative: 0,
  };
}

function addCallMetric(metric, score) {
  if (!score) return;
  metric.records += 1;
  metric.brier += score.brier;
  metric.truthProbability += score.truthProbability;
  metric.containsTruth += Number(score.containsTruth);
  metric.exactSingleton += Number(score.exactSingleton);
  metric.predictedCount += score.predictedCount;
  metric.ambiguous += Number(score.predictedCount > 1);
  metric.truthAbsent += Number(score.truthAbsent);
  metric.trueNegative += Number(score.truthAbsent && score.predictedAbsent);
}

function summarizeCall(metric) {
  const divide = (value, denominator = metric.records) => denominator ? value / denominator : null;
  return {
    records: metric.records,
    meanBrier: divide(metric.brier),
    meanTruthProbability: divide(metric.truthProbability),
    truthInReportedSet: divide(metric.containsTruth),
    exactSingleton: divide(metric.exactSingleton),
    callAmbiguityRate: divide(metric.ambiguous),
    meanReportedClasses: divide(metric.predictedCount),
    truthAbsentRecords: metric.truthAbsent,
    truthAbsentSpecificity: divide(metric.trueNegative, metric.truthAbsent),
  };
}

function emptyBoundaryMetric() {
  return { records: 0, sum: 0, absolute: 0, exact: 0, withinOne: 0, errors: [] };
}

function addBoundary(metric, truthValue, predictionValue) {
  if (truthValue === "" || predictionValue === "") return;
  const truth = Number(truthValue);
  const prediction = Number(predictionValue);
  if (!Number.isFinite(truth) || !Number.isFinite(prediction)) return;
  const error = prediction - truth;
  metric.records += 1;
  metric.sum += error;
  metric.absolute += Math.abs(error);
  metric.exact += Number(error === 0);
  metric.withinOne += Number(Math.abs(error) <= 1);
  metric.errors.push(error);
}

function quantile(values, probability) {
  if (!values.length) return null;
  const ordered = [...values].sort((a, b) => a - b);
  const index = (ordered.length - 1) * probability;
  const lower = Math.floor(index);
  const upper = Math.ceil(index);
  return ordered[lower] + (ordered[upper] - ordered[lower]) * (index - lower);
}

function summarizeBoundary(metric) {
  const divide = (value) => metric.records ? value / metric.records : null;
  return {
    records: metric.records,
    meanSignedError: divide(metric.sum),
    meanAbsoluteError: divide(metric.absolute),
    exact: divide(metric.exact),
    withinOne: divide(metric.withinOne),
    q025SignedError: quantile(metric.errors, 0.025),
    medianSignedError: quantile(metric.errors, 0.5),
    q975SignedError: quantile(metric.errors, 0.975),
    p95AbsoluteError: quantile(metric.errors.map(Math.abs), 0.95),
  };
}

const truthText = readText(paths.truth);
const allTruth = /^\s*\u200b?>/.test(truthText)
  ? parseHeaderTruthFasta(truthText)
  : parseTable(truthText);
let truth = allTruth.filter((row) => selected(row.sequence_id));
if (limit > 0) truth = truth.slice(0, limit);
const truthById = new Map(truth.map((row) => [row.sequence_id, row]));

const encoder = new TextEncoder();
const decoder = new TextDecoder();
const wasi = new WASI([], [], []);
const module = await WebAssembly.compile(fs.readFileSync(paths.wasm));
const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
wasi.initialize(instance);
const runtime = instance.exports;

const profileResult = runtime.swig_set_calling_profile(profileNumbers[profile]);
if (profileResult !== 0) throw new Error(`WASM rejected profile ${profile}`);
if (runtime.swig_set_assigner_strategy(3) !== 0) throw new Error("WASM rejected AER-R");
if (tuning) {
  const defaults = profile === "r_optimized" ? {
    vMatch: 2, vMismatch: -4, vGapOpen: -13, vGapExtend: -1,
    dMatch: 2, dMismatch: -4, dGapOpen: -13, dGapExtend: -1, topD: 2, minDMatch: 5,
    jMatch: 2, jMismatch: -4, jGapOpen: -17, jGapExtend: -2, topJ: 2, minJLength: 10,
  } : {
    vMatch: 2, vMismatch: -3, vGapOpen: -5, vGapExtend: -1,
    dMatch: 2, dMismatch: -3, dGapOpen: -13, dGapExtend: -1, topD: 2, minDMatch: 6,
    jMatch: 2, jMismatch: -3, jGapOpen: -17, jGapExtend: -2, topJ: 2, minJLength: 10,
  };
  const options = { ...defaults, ...tuning };
  const accepted = runtime.swig_set_tuning_options(
    options.dMatch, options.dMismatch, options.dGapOpen, options.dGapExtend,
    options.topD, options.minDMatch, options.jMatch, options.jMismatch,
    options.jGapOpen, options.jGapExtend, options.topJ, options.minJLength,
  );
  if (accepted !== 0) throw new Error("WASM rejected tuning options");
  if (profile === "r_optimized" || "vMatch" in tuning || "vMismatch" in tuning || "vGapOpen" in tuning || "vGapExtend" in tuning) {
    if (typeof runtime.swig_set_v_tuning_options !== "function") {
      throw new Error("WASM does not expose V tuning options");
    }
    const vAccepted = runtime.swig_set_v_tuning_options(
      options.vMatch,
      options.vMismatch,
      options.vGapOpen,
      options.vGapExtend,
      profile === "r_optimized" ? 1 : 0,
    );
    if (vAccepted !== 0) throw new Error("WASM rejected V tuning options");
  }
  if ("dPresencePenalty" in tuning) {
    if (typeof runtime.swig_set_aer_r_decision_tuning !== "function") {
      throw new Error("WASM does not expose AER-R decision tuning");
    }
    if (runtime.swig_set_aer_r_decision_tuning(tuning.dPresencePenalty) !== 0) {
      throw new Error("WASM rejected AER-R decision tuning");
    }
  }
}

function putText(value) {
  const bytes = encoder.encode(value);
  const pointer = runtime.swig_alloc(bytes.byteLength);
  if (!pointer && bytes.byteLength) throw new Error("WASM allocation failed");
  new Uint8Array(runtime.memory.buffer, pointer, bytes.byteLength).set(bytes);
  return [pointer, bytes.byteLength];
}

function errorText() {
  return decoder.decode(new Uint8Array(runtime.memory.buffer, runtime.swig_error_ptr(), runtime.swig_error_len()));
}

const referenceAllocations = [referenceText.V, referenceText.D, referenceText.J, ""].map(putText);
const genes = runtime.swig_init_database(...referenceAllocations.flat());
referenceAllocations.forEach(([pointer]) => runtime.swig_free(pointer));
if (genes < 0) throw new Error(errorText());

let balancedFilter = null;
if (profile === "igblast_balanced") {
  balancedFilter = (await import("../src/balanced-calling-profile.ts")).applyBalancedDFilter;
}

const callMetrics = { V: emptyCallMetric(), D: emptyCallMetric(), J: emptyCallMetric() };
const compatibleBoundaryMetrics = Object.fromEntries(
  ["v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end"]
    .map((field) => [field, emptyBoundaryMetric()]),
);
const allBoundaryMetrics = Object.fromEntries(
  Object.keys(compatibleBoundaryMetrics).map((field) => [field, emptyBoundaryMetric()]),
);
const predictionFields = [
  "sequence_id", "v_call", "d_call", "j_call",
  "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
  "j_sequence_start", "j_sequence_end", "v_score", "d_score", "j_score",
  "v_identity", "d_identity", "j_identity", "v_cigar", "d_cigar", "j_cigar",
  "v_sequence_alignment", "v_germline_alignment", "d_sequence_alignment",
  "d_germline_alignment", "j_sequence_alignment", "j_germline_alignment",
  "v_alternatives", "d_alternatives", "j_alternatives",
];
const predictionLines = [predictionFields.join("\t")];
const failures = [];
let tandemD = 0;
let processed = 0;
const started = performance.now();
for (let offset = 0; offset < truth.length; offset += batchSize) {
  const batch = truth.slice(offset, offset + batchSize);
  const fasta = batch.map((row) => `>${row.sequence_id}\n${row.sequence}\n`).join("");
  const [pointer, length] = putText(fasta);
  const count = runtime.swig_annotate(pointer, length, 1, minimumIdentity, strand);
  runtime.swig_free(pointer);
  if (count < 0) throw new Error(errorText());
  const result = new Uint8Array(runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len());
  const newline = result.indexOf(10);
  if (newline < 0) throw new Error("SwiftIG returned an invalid AIRR table");
  const header = decoder.decode(result.subarray(0, newline)).replace(/\r$/, "");
  let body = result.slice(newline + 1);
  if (balancedFilter) body = balancedFilter(header, body).body;
  const predictions = parseTable(`${header}\n${decoder.decode(body)}`);
  if (predictions.length !== batch.length) throw new Error(`Batch mismatch at ${offset}`);
  for (const prediction of predictions) {
    const actual = truthById.get(prediction.sequence_id);
    if (!actual) throw new Error(`Unexpected prediction ${prediction.sequence_id}`);
    const scores = {};
    for (const segment of ["V", "D", "J"]) {
      const lower = segment.toLowerCase();
      const score = scoreCall(segment, actual[`${lower}_call`], prediction[`${lower}_call`]);
      scores[segment] = score;
      if (segment === "D" && score === null) tandemD += 1;
      addCallMetric(callMetrics[segment], score);
      const compatible = score?.containsTruth ?? false;
      for (const suffix of ["sequence_start", "sequence_end"]) {
        const field = `${lower}_${suffix}`;
        addBoundary(allBoundaryMetrics[field], actual[field], prediction[field]);
        if (compatible) addBoundary(compatibleBoundaryMetrics[field], actual[field], prediction[field]);
      }
    }
    const boundaryErrors = Object.fromEntries(Object.keys(allBoundaryMetrics).map((field) => [
      field,
      actual[field] !== "" && prediction[field] !== "" ? Number(prediction[field]) - Number(actual[field]) : null,
    ]));
    if (["V", "D", "J"].some((segment) => scores[segment] && scores[segment].brier > 0) ||
        Object.values(boundaryErrors).some((error) => error !== null && Math.abs(error) > 1)) {
      failures.push({
        sequence_id: actual.sequence_id,
        truth: Object.fromEntries(["v_call", "d_call", "j_call", ...Object.keys(allBoundaryMetrics)].map((field) => [field, actual[field]])),
        prediction: Object.fromEntries(predictionFields.filter((field) => field !== "sequence_id").map((field) => [field, prediction[field] ?? ""])),
        scores,
        boundaryErrors,
      });
    }
    if (outputPath) predictionLines.push(predictionFields.map((field) => prediction[field] ?? "").join("\t"));
    processed += 1;
  }
  if (processed % 2_000 === 0 || processed === truth.length) {
    process.stderr.write(`\r${processed.toLocaleString()} / ${truth.length.toLocaleString()}`);
  }
}
process.stderr.write("\n");
const elapsedSeconds = (performance.now() - started) / 1_000;
if (outputPath) fs.writeFileSync(outputPath, `${predictionLines.join("\n")}\n`);
if (failuresPath) fs.writeFileSync(failuresPath, failures.map((record) => JSON.stringify(record)).join("\n") + "\n");

const report = {
  split,
  profile,
  records: truth.length,
  tandemDExcludedFromPrimaryScore: tandemD,
  references: genes,
  elapsedSeconds,
  readsPerSecond: truth.length / elapsedSeconds,
  callMetrics: Object.fromEntries(Object.entries(callMetrics).map(([segment, metric]) => [segment, summarizeCall(metric)])),
  boundaryMetrics: {
    allCalled: Object.fromEntries(Object.entries(allBoundaryMetrics).map(([field, metric]) => [field, summarizeBoundary(metric)])),
    truthCompatibleCall: Object.fromEntries(Object.entries(compatibleBoundaryMetrics).map(([field, metric]) => [field, summarizeBoundary(metric)])),
  },
  failureRecords: failures.length,
};
if (args.get("compact") === "true") {
  console.log(JSON.stringify({
    profile: report.profile,
    records: report.records,
    tuning,
    elapsedSeconds: report.elapsedSeconds,
    readsPerSecond: report.readsPerSecond,
    V: report.callMetrics.V,
    D: report.callMetrics.D,
    J: report.callMetrics.J,
    vEnd: report.boundaryMetrics.truthCompatibleCall.v_sequence_end,
    dStart: report.boundaryMetrics.truthCompatibleCall.d_sequence_start,
    dEnd: report.boundaryMetrics.truthCompatibleCall.d_sequence_end,
    jStart: report.boundaryMetrics.truthCompatibleCall.j_sequence_start,
  }));
} else {
  console.log(JSON.stringify(report, null, 2));
}
