import fs from "node:fs";
import path from "node:path";
import { performance } from "node:perf_hooks";

import { WASI } from "@bjorn3/browser_wasi_shim";

import {
  generateVdjDataset,
  parseFasta,
  recordsToFasta,
  referenceEquivalence,
} from "./vdj-simulator.mjs";

function argumentsMap(argv) {
  const values = new Map();
  for (let index = 0; index < argv.length; index += 1) {
    const name = argv[index];
    if (!name.startsWith("--")) continue;
    const next = argv[index + 1];
    values.set(name.slice(2), next && !next.startsWith("--") ? argv[++index] : "true");
  }
  return values;
}

const args = argumentsMap(process.argv.slice(2));
const root = new URL("../", import.meta.url);
const referenceRoot = args.get("references")
  ? path.resolve(args.get("references"))
  : new URL("../public/references/kimdb-1.1/Macaca_mulatta/IGH/", import.meta.url);
const resolveReference = (name) => fs.readFileSync(
  referenceRoot instanceof URL ? new URL(name, referenceRoot) : path.join(referenceRoot, name),
  "utf8",
);
const referenceText = {
  V: resolveReference("V.fasta"),
  D: resolveReference("D.fasta"),
  J: resolveReference("J.fasta"),
  C: "",
};
const references = {
  V: parseFasta(referenceText.V),
  D: parseFasta(referenceText.D),
  J: parseFasta(referenceText.J),
};
const count = Number(args.get("count") ?? 2_000);
const seed = Number(args.get("seed") ?? 913_771);
const batchSize = Number(args.get("batch") ?? 250);
const wasmPath = path.resolve(args.get("wasm") ?? new URL("../public/swiftig.wasm", import.meta.url).pathname);
const oraclePath = args.get("oracle-wasm") ? path.resolve(args.get("oracle-wasm")) : "";
const oracleCount = Number(args.get("oracle-count") ?? 40);
const outputPrefix = args.get("out") ? path.resolve(args.get("out")) : "";
if (!Number.isInteger(count) || count < 1) throw new Error("--count must be a positive integer");

const records = generateVdjDataset({
  ...references,
  count,
  seed,
  doubleDRate: Number(args.get("double-d-rate") ?? 0.025),
  indelRate: Number(args.get("indel-rate") ?? 0.0015),
  sequencingErrorRate: Number(args.get("error-rate") ?? 0.0008),
  ambiguousRate: Number(args.get("n-rate") ?? 0.00025),
  reverseRate: Number(args.get("reverse-rate") ?? 0.08),
});

function parseAirr(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const headers = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
  });
}

async function makeRuntime(binary, strategy) {
  const bytes = fs.readFileSync(binary);
  const wasi = new WASI([], [], []);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
  wasi.initialize(instance);
  const runtime = instance.exports;
  const encoder = new TextEncoder();
  const decoder = new TextDecoder();
  const put = (value) => {
    const bytesValue = encoder.encode(value);
    const pointer = runtime.swig_alloc(bytesValue.byteLength);
    if (!pointer && bytesValue.byteLength) throw new Error("WASM allocation failed");
    new Uint8Array(runtime.memory.buffer, pointer, bytesValue.byteLength).set(bytesValue);
    return [pointer, bytesValue.byteLength];
  };
  const error = () => decoder.decode(new Uint8Array(
    runtime.memory.buffer, runtime.swig_error_ptr(), runtime.swig_error_len(),
  ));
  if (runtime.swig_set_assigner_strategy(strategy) !== 0) {
    throw new Error(`WASM rejected strategy ${strategy}`);
  }
  if (runtime.swig_set_calling_profile(0) !== 0) throw new Error("WASM rejected the truth profile");
  const allocations = [referenceText.V, referenceText.D, referenceText.J, ""].map(put);
  const genes = runtime.swig_init_database(...allocations.flat());
  allocations.forEach(([pointer]) => runtime.swig_free(pointer));
  if (genes < 0) throw new Error(error());
  return {
    genes,
    annotate(fasta) {
      const [pointer, length] = put(fasta);
      const returned = runtime.swig_annotate(pointer, length, 1, 600, 0);
      runtime.swig_free(pointer);
      if (returned < 0) throw new Error(error());
      const result = decoder.decode(new Uint8Array(
        runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len(),
      ));
      return parseAirr(result);
    },
    annotateDoubleD(fasta) {
      const [pointer, length] = put(fasta);
      const returned = runtime.swig_annotate_double_d(
        pointer, length, 1, 600, 0,
        1, 40, 11, 5, 3, 8,
      );
      runtime.swig_free(pointer);
      if (returned < 0) throw new Error(error());
      const airr = decoder.decode(new Uint8Array(
        runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len(),
      ));
      const doubleD = decoder.decode(new Uint8Array(
        runtime.memory.buffer,
        runtime.swig_double_d_result_ptr(),
        runtime.swig_double_d_result_len(),
      ));
      return { airr: parseAirr(airr), doubleD: parseAirr(doubleD) };
    },
  };
}

async function annotateAll(binary, strategy, selected = records, progress = false) {
  const runtime = await makeRuntime(binary, strategy);
  const rows = [];
  const started = performance.now();
  for (let offset = 0; offset < selected.length; offset += batchSize) {
    const batch = selected.slice(offset, offset + batchSize);
    rows.push(...runtime.annotate(recordsToFasta(batch)));
    if (progress && (offset + batch.length === selected.length || (offset + batch.length) % 1_000 === 0)) {
      process.stderr.write(`\r${offset + batch.length} / ${selected.length}`);
    }
  }
  if (progress) process.stderr.write("\n");
  const seconds = (performance.now() - started) / 1_000;
  return { rows, seconds, readsPerSecond: selected.length / seconds, genes: runtime.genes };
}

async function screenDoubleD(binary, strategy, selected = records) {
  const runtime = await makeRuntime(binary, strategy);
  const rows = [];
  const started = performance.now();
  for (let offset = 0; offset < selected.length; offset += batchSize) {
    const batch = selected.slice(offset, offset + batchSize);
    const batchStarted = performance.now();
    rows.push(...runtime.annotateDoubleD(recordsToFasta(batch)).doubleD);
    const batchMilliseconds = performance.now() - batchStarted;
    if (args.get("profile-slow") === "true" && batchMilliseconds > 50) {
      process.stderr.write(`slow double-D strategy=${strategy} offset=${offset} count=${batch.length} ms=${batchMilliseconds.toFixed(1)} ids=${batch.map((record) => record.id).join(",")}\n`);
    }
  }
  const seconds = (performance.now() - started) / 1_000;
  return { rows, seconds, readsPerSecond: selected.length / seconds };
}

function callAtoms(value) {
  return String(value ?? "").split(/[,/]/).map((item) => item.trim()).filter(Boolean);
}

const equivalence = Object.fromEntries(Object.entries(references).map(([segment, values]) => [
  segment,
  referenceEquivalence(values),
]));

function equivalentCall(segment, truthName, predicted) {
  if (!truthName) return !callAtoms(predicted).length;
  const truthSequence = equivalence[segment].byName.get(truthName);
  return callAtoms(predicted).some((name) =>
    name === truthName || (truthSequence && equivalence[segment].byName.get(name) === truthSequence));
}

function orientedSpan(record, span) {
  if (!span || !record.truth.reversed) return span;
  return { start: record.sequence.length - span.end, end: record.sequence.length - span.start, bases: span.bases };
}

function predictedSpan(row, prefix) {
  const start = Number(row[`${prefix}_sequence_start`]);
  const end = Number(row[`${prefix}_sequence_end`]);
  return Number.isFinite(start) && Number.isFinite(end) && start > 0 && end >= start
    ? { start: start - 1, end }
    : null;
}

function overlap(left, right) {
  if (!left || !right) return 0;
  return Math.max(0, Math.min(left.end, right.end) - Math.max(left.start, right.start));
}

function evaluate(record, row) {
  const dTruth = record.truth.dCalls[0] ?? "";
  const dSpan = orientedSpan(record, record.truth.spans.D1);
  const dEligible = !record.truth.tandem && Boolean(dSpan && dSpan.bases >= 6);
  const dStrong = dEligible && dSpan.bases >= 10;
  const predictedD = predictedSpan(row, "d");
  const dCalled = callAtoms(row.d_call).length > 0;
  const compatibleD = new Set(record.truth.dCompatibleCalls?.[0] ?? [dTruth]);
  return {
    vExact: callAtoms(row.v_call).includes(record.truth.vCall),
    vEquivalent: equivalentCall("V", record.truth.vCall, row.v_call),
    jExact: callAtoms(row.j_call).includes(record.truth.jCall),
    jEquivalent: equivalentCall("J", record.truth.jCall, row.j_call),
    dEligible,
    dStrong,
    dCalled,
    dExact: dEligible && callAtoms(row.d_call).includes(dTruth),
    dEquivalent: dEligible && equivalentCall("D", dTruth, row.d_call),
    dCompatible: dEligible && callAtoms(row.d_call).some((call) => compatibleD.has(call)),
    dOverlap: dEligible && predictedD ? overlap(dSpan, predictedD) : 0,
    dTruthBases: dSpan?.bases ?? 0,
    dStartError: dEligible && predictedD ? Math.abs(predictedD.start - dSpan.start) : null,
    dEndError: dEligible && predictedD ? Math.abs(predictedD.end - dSpan.end) : null,
    falseD: !record.truth.tandem && (!dSpan || dSpan.bases < 6) && dCalled,
  };
}

function mean(values) {
  const finite = values.filter(Number.isFinite);
  return finite.length ? finite.reduce((sum, value) => sum + value, 0) / finite.length : null;
}

function quantile(values, probability) {
  const finite = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (!finite.length) return null;
  return finite[Math.min(finite.length - 1, Math.floor(probability * finite.length))];
}

function summarize(evaluations, timing) {
  const fraction = (predicate, selected = evaluations) => selected.length
    ? selected.filter(predicate).length / selected.length : null;
  const eligible = evaluations.filter((value) => value.dEligible);
  const strong = evaluations.filter((value) => value.dStrong);
  return {
    records: evaluations.length,
    seconds: timing.seconds,
    readsPerSecond: timing.readsPerSecond,
    vExact: fraction((value) => value.vExact),
    vSequenceEquivalent: fraction((value) => value.vEquivalent),
    jExact: fraction((value) => value.jExact),
    jSequenceEquivalent: fraction((value) => value.jEquivalent),
    dEligibleRecords: eligible.length,
    dDetection: fraction((value) => value.dCalled, eligible),
    dExact: fraction((value) => value.dExact, eligible),
    dSequenceEquivalent: fraction((value) => value.dEquivalent, eligible),
    dRetainedTractCompatible: fraction((value) => value.dCompatible, eligible),
    strongDRecords: strong.length,
    strongDDetection: fraction((value) => value.dCalled, strong),
    strongDSequenceEquivalent: fraction((value) => value.dEquivalent, strong),
    strongDRetainedTractCompatible: fraction((value) => value.dCompatible, strong),
    falseDOnSubthresholdTruth: fraction((value) => value.falseD),
    dBoundary: {
      meanStartAbsoluteError: mean(eligible.map((value) => value.dStartError)),
      p95StartAbsoluteError: quantile(eligible.map((value) => value.dStartError), 0.95),
      meanEndAbsoluteError: mean(eligible.map((value) => value.dEndError)),
      p95EndAbsoluteError: quantile(eligible.map((value) => value.dEndError), 0.95),
      meanTruthOverlap: mean(eligible.map((value) => value.dOverlap / value.dTruthBases)),
    },
  };
}

process.stderr.write("AER simulation benchmark\n");
const aer = await annotateAll(wasmPath, 2, records, true);
process.stderr.write("AER-R simulation benchmark\n");
const robust = await annotateAll(wasmPath, 3, records, true);
const aerRows = new Map(aer.rows.map((row) => [row.sequence_id, row]));
const robustRows = new Map(robust.rows.map((row) => [row.sequence_id, row]));
const aerEvaluation = records.map((record) => evaluate(record, aerRows.get(record.id)));
const robustEvaluation = records.map((record) => evaluate(record, robustRows.get(record.id)));

const cases = [];
for (let index = 0; index < records.length; index += 1) {
  const record = records[index];
  const a = aerEvaluation[index];
  const r = robustEvaluation[index];
  let category = "";
  if (a.dStrong && !a.dCompatible && r.dCompatible) category = "strong_D_recovered";
  else if (a.dEligible && !a.dCompatible && r.dCompatible) category = "D_recovered";
  else if (a.dCompatible && !r.dCompatible) category = "AER_R_D_regression";
  else if (r.dStrong && !r.dCompatible) category = "strong_D_unresolved";
  else if (a.vEquivalent && !r.vEquivalent) category = "AER_R_V_regression";
  else if (a.jEquivalent && !r.jEquivalent) category = "AER_R_J_regression";
  if (!category) continue;
  cases.push({
    category,
    id: record.id,
    truth: record.truth,
    sequence: record.sequence,
    aer: aerRows.get(record.id),
    aerRobust: robustRows.get(record.id),
  });
}

let oracle = null;
if (oraclePath) {
  const priority = [
    ...cases,
    ...records.filter((record, index) => robustEvaluation[index].dStrong && !robustEvaluation[index].dCompatible)
      .map((record) => ({ id: record.id })),
  ];
  const selectedIds = new Set(priority.slice(0, oracleCount).map((value) => value.id));
  for (const record of records) {
    if (selectedIds.size >= oracleCount) break;
    if (((Number(record.id.slice(4)) * 2654435761) >>> 0) % 97 < 4) selectedIds.add(record.id);
  }
  const selected = records.filter((record) => selectedIds.has(record.id)).slice(0, oracleCount);
  process.stderr.write(`Exhaustive oracle (${selected.length} selected records)\n`);
  const result = await annotateAll(oraclePath, 3, selected, true);
  const rows = new Map(result.rows.map((row) => [row.sequence_id, row]));
  let oracleRescues = 0;
  let productionMatches = 0;
  const classifications = [];
  for (const record of selected) {
    const production = evaluate(record, robustRows.get(record.id));
    const exhaustive = evaluate(record, rows.get(record.id));
    if (production.dCompatible === exhaustive.dCompatible &&
        production.vEquivalent === exhaustive.vEquivalent &&
        production.jEquivalent === exhaustive.jEquivalent) productionMatches += 1;
    if (!production.dCompatible && exhaustive.dCompatible) oracleRescues += 1;
    classifications.push({
      id: record.id,
      production,
      exhaustive,
      cause: !production.dCompatible && exhaustive.dCompatible
        ? "production_candidate_or_fast_path_pruning"
        : !exhaustive.dCompatible && exhaustive.dStrong
          ? "joint_scoring_or_identifiability"
          : "agreement",
      oracleRow: rows.get(record.id),
    });
  }
  oracle = {
    records: selected.length,
    seconds: result.seconds,
    productionExactAgreement: selected.length ? productionMatches / selected.length : null,
    DRescuesOverProduction: oracleRescues,
    classifications,
  };
}

let doubleD = null;
if (args.get("screen-double-d") === "true") {
  const baseline = await screenDoubleD(wasmPath, 2);
  const experimental = await screenDoubleD(wasmPath, 3);
  const summarizeDoubleD = (result) => {
    const byId = new Map(result.rows.map((row) => [row.sequence_id, row]));
    const tandem = records.filter((record) => record.truth.tandem);
    const single = records.filter((record) => !record.truth.tandem);
    let orderedCalls = 0;
    for (const record of tandem) {
      const row = byId.get(record.id);
      if (!row) continue;
      if (equivalentCall("D", record.truth.dCalls[0], row.d_call) &&
          equivalentCall("D", record.truth.dCalls[1], row.d2_call)) orderedCalls += 1;
    }
    return {
      seconds: result.seconds,
      readsPerSecond: result.readsPerSecond,
      tandemRecords: tandem.length,
      detection: tandem.length ? tandem.filter((record) => byId.has(record.id)).length / tandem.length : null,
      orderedSequenceEquivalentCalls: tandem.length ? orderedCalls / tandem.length : null,
      falsePositiveRate: single.length ? single.filter((record) => byId.has(record.id)).length / single.length : null,
    };
  };
  doubleD = { aer: summarizeDoubleD(baseline), aerRobust: summarizeDoubleD(experimental) };
}

const report = {
  generatedAt: new Date().toISOString(),
  simulator: {
    seed,
    records: count,
    references: { V: references.V.length, D: references.D.length, J: references.J.length },
    configuration: {
      doubleDRate: Number(args.get("double-d-rate") ?? 0.025),
      indelRate: Number(args.get("indel-rate") ?? 0.0015),
      sequencingErrorRate: Number(args.get("error-rate") ?? 0.0008),
      ambiguousRate: Number(args.get("n-rate") ?? 0.00025),
      reverseRate: Number(args.get("reverse-rate") ?? 0.08),
    },
    observed: {
      tandemD: records.filter((record) => record.truth.tandem).length,
      reversed: records.filter((record) => record.truth.reversed).length,
      indelBearing: records.filter((record) => record.truth.indelEvents > 0).length,
      meanShmMutations: mean(records.map((record) => record.truth.shmMutations)),
    },
  },
  aer: summarize(aerEvaluation, aer),
  aerRobust: summarize(robustEvaluation, robust),
  comparison: {
    throughputRatio: robust.readsPerSecond / aer.readsPerSecond,
    investigatedCases: cases.length,
    categories: Object.fromEntries([...new Set(cases.map((value) => value.category))]
      .map((category) => [category, cases.filter((value) => value.category === category).length])),
  },
  oracle,
  doubleD,
};

if (outputPrefix) {
  fs.mkdirSync(path.dirname(outputPrefix), { recursive: true });
  fs.writeFileSync(`${outputPrefix}.report.json`, `${JSON.stringify(report, null, 2)}\n`);
  fs.writeFileSync(`${outputPrefix}.cases.jsonl`, cases.map((value) => JSON.stringify(value)).join("\n") + (cases.length ? "\n" : ""));
  fs.writeFileSync(`${outputPrefix}.simulated.fasta`, recordsToFasta(records));
  fs.writeFileSync(`${outputPrefix}.truth.jsonl`, records.map((record) => JSON.stringify({ id: record.id, ...record.truth })).join("\n") + "\n");
}
process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
