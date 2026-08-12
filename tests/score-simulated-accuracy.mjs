import fs from "node:fs";

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
const callsPath = args.get("calls");
if (!callsPath) throw new Error("Missing --calls");
const method = String(args.get("method") ?? "predictions").toLowerCase();
const predictionPaths = String(args.get("predictions") ?? "").split(",").filter(Boolean);
if (method === "predictions" && !predictionPaths.length) {
  throw new Error("Supply --predictions file1.tsv[,file2.tsv,...] or --method igblast|swiftig|igformer");
}

function parseTable(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(header.map((name, index) => [name, values[index] ?? ""]));
  });
}

function splitAtoms(value) {
  return String(value ?? "").split(/[,/]/).map((item) => item.trim()).filter(Boolean);
}

function metric() {
  return {
    records: 0, top: 0, fair: 0, exactSet: 0, called: 0, multiple: 0, callCount: 0,
    truthPositive: 0, truthPositiveTop: 0, truthPositiveFair: 0,
    truthNegative: 0, truthNegativeCorrect: 0,
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

const truthRows = parseTable(fs.readFileSync(callsPath, "utf8"));
let predictions = null;
if (predictionPaths.length) {
  predictions = new Map();
  for (const path of predictionPaths) {
    for (const row of parseTable(fs.readFileSync(path, "utf8"))) {
      if (predictions.has(row.sequence_id)) throw new Error(`Duplicate prediction: ${row.sequence_id}`);
      predictions.set(row.sequence_id, row);
    }
  }
  if (predictions.size !== truthRows.length) {
    throw new Error(`Prediction coverage mismatch: ${predictions.size} predictions for ${truthRows.length} truths`);
  }
}

const prefix = method === "predictions" ? "" : `${method}_`;
const metrics = { V: metric(), D: metric(), J: metric() };
for (const truth of truthRows) {
  const prediction = predictions?.get(truth.sequence_id);
  if (predictions && !prediction) throw new Error(`Missing prediction: ${truth.sequence_id}`);
  for (const segment of ["v", "d", "j"]) {
    updateMetric(
      metrics[segment.toUpperCase()],
      truth[`true_${segment}_call`],
      prediction ? prediction[`${segment}_call`] : truth[`${prefix}${segment}_call`],
    );
  }
}

console.log(JSON.stringify({
  method,
  records: truthRows.length,
  predictions: predictionPaths,
  metrics: Object.fromEntries(Object.entries(metrics).map(([segment, value]) => [segment, summarize(value)])),
}, null, 2));
