import fs from "node:fs";

function table(path) {
  const lines = fs.readFileSync(path, "utf8").trimEnd().split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(header.map((name, index) => [name, values[index] ?? ""]));
  });
}

function atoms(value) {
  return [...new Set(String(value ?? "").split(/[,/]/).map((item) => item.trim()).filter(Boolean))];
}

function scores(referenceValue, predictionValue) {
  const reference = new Set(atoms(referenceValue));
  const prediction = atoms(predictionValue);
  if (!reference.size) return { top: prediction.length ? 0 : 1, fair: prediction.length ? 0 : 1 };
  return {
    top: prediction.length && reference.has(prediction[0]) ? 1 : 0,
    fair: prediction.length ? new Set(prediction.filter((call) => reference.has(call))).size / prediction.length : 0,
  };
}

function longestExact(query, germline) {
  let longest = 0;
  let current = 0;
  for (let index = 0; index < Math.min(query.length, germline.length); index += 1) {
    if (query[index] !== "-" && query[index] === germline[index]) longest = Math.max(longest, ++current);
    else current = 0;
  }
  return longest;
}

const calls = new Map(table(process.argv[2]).map((row) => [row.sequence_id, row]));
const predictions = table(process.argv[3]).map((prediction) => {
  const call = calls.get(prediction.sequence_id);
  if (!call) throw new Error(`Missing call row for ${prediction.sequence_id}`);
  const query = prediction.d_sequence_alignment ?? "";
  const germline = prediction.d_germline_alignment ?? "";
  const columns = Math.max(query.length, germline.length);
  const aligned = [...query].filter((base) => base !== "-").length;
  const gaps = [...query].filter((base) => base === "-").length + [...germline].filter((base) => base === "-").length;
  const feature = {
    exact: longestExact(query, germline),
    score: Number(prediction.d_score || 0),
    identity: Number(prediction.d_identity || 0),
    columns,
    aligned,
    gaps,
    calls: atoms(prediction.d_call).length,
    alternatives: String(prediction.d_alternatives ?? "").split(";").filter(Boolean).length,
    vjSpan: Math.max(0, Number(prediction.j_sequence_start || 0) - Number(prediction.v_sequence_end || 0)),
    leftFlank: Math.max(0, Number(prediction.d_sequence_start || 0) - Number(prediction.v_sequence_end || 0)),
    rightFlank: Math.max(0, Number(prediction.j_sequence_start || 0) - Number(prediction.d_sequence_end || 0)),
    junctionLength: Number(prediction.junction_length || 0),
  };
  return { prediction, call, feature };
});

const baseline = {
  truthTop: predictions.reduce((sum, row) => sum + scores(row.call.true_d_call, row.prediction.d_call).top, 0),
  truthFair: predictions.reduce((sum, row) => sum + scores(row.call.true_d_call, row.prediction.d_call).fair, 0),
  agreementTop: predictions.reduce((sum, row) => sum + scores(row.call.igblast_d_call, row.prediction.d_call).top, 0),
  agreementFair: predictions.reduce((sum, row) => sum + scores(row.call.igblast_d_call, row.prediction.d_call).fair, 0),
};

function evaluate(name, drop) {
  const totals = { ...baseline };
  let dropped = 0;
  let truthNegative = 0;
  let igblastNegative = 0;
  for (const row of predictions) {
    if (!row.prediction.d_call || !drop(row.feature)) continue;
    dropped += 1;
    if (!atoms(row.call.true_d_call).length) truthNegative += 1;
    if (!atoms(row.call.igblast_d_call).length) igblastNegative += 1;
    for (const [prefix, reference] of [["truth", row.call.true_d_call], ["agreement", row.call.igblast_d_call]]) {
      const before = scores(reference, row.prediction.d_call);
      const after = scores(reference, "");
      totals[`${prefix}Top`] += after.top - before.top;
      totals[`${prefix}Fair`] += after.fair - before.fair;
    }
  }
  const records = predictions.length;
  return {
    name,
    dropped,
    truthNegative,
    igblastNegative,
    truthTop: totals.truthTop / records,
    truthFair: totals.truthFair / records,
    agreementTop: totals.agreementTop / records,
    agreementFair: totals.agreementFair / records,
  };
}

const rules = [];
const add = (name, rule) => rules.push(evaluate(name, rule));
for (let score = 8; score <= 40; score += 1) add(`exact5 score<=${score}`, (f) => f.exact === 5 && f.score <= score);
for (let length = 5; length <= 24; length += 1) add(`exact5 columns<=${length}`, (f) => f.exact === 5 && f.columns <= length);
for (let identity = 0.60; identity <= 1.0001; identity += 0.025) add(`exact5 identity<=${identity.toFixed(3)}`, (f) => f.exact === 5 && f.identity <= identity + 1e-9);
for (let ratio = 0.5; ratio <= 2.5001; ratio += 0.1) add(`exact5 score/column<=${ratio.toFixed(1)}`, (f) => f.exact === 5 && f.columns && f.score / f.columns <= ratio + 1e-9);
for (let length = 5; length <= 18; length += 1) {
  for (let score = 8; score <= 32; score += 2) add(`exact5 columns<=${length} score<=${score}`, (f) => f.exact === 5 && f.columns <= length && f.score <= score);
}
for (let length = 5; length <= 18; length += 1) {
  for (let identity = 0.60; identity <= 1.0001; identity += 0.05) add(`exact5 columns<=${length} identity<=${identity.toFixed(2)}`, (f) => f.exact === 5 && f.columns <= length && f.identity <= identity + 1e-9);
}
for (const maximumCalls of [1, 2, 3, 4, 5]) add(`exact5 calls<=${maximumCalls}`, (f) => f.exact === 5 && f.calls <= maximumCalls);
for (const minimumCalls of [2, 3, 4, 5]) add(`exact5 calls>=${minimumCalls}`, (f) => f.exact === 5 && f.calls >= minimumCalls);
for (const gaps of [0, 1, 2, 3]) add(`exact5 gaps>=${gaps}`, (f) => f.exact === 5 && f.gaps >= gaps);
for (const feature of ["vjSpan", "leftFlank", "rightFlank", "junctionLength", "alternatives"]) {
  const values = [...new Set(predictions.map((row) => row.feature[feature]))].sort((a, b) => a - b);
  for (const threshold of values) {
    add(`exact5 ${feature}>=${threshold}`, (f) => f.exact === 5 && f[feature] >= threshold);
    add(`exact5 ${feature}<=${threshold}`, (f) => f.exact === 5 && f[feature] <= threshold);
    add(`exact5 calls>=3 ${feature}>=${threshold}`, (f) => f.exact === 5 && f.calls >= 3 && f[feature] >= threshold);
    add(`exact5 calls>=3 ${feature}<=${threshold}`, (f) => f.exact === 5 && f.calls >= 3 && f[feature] <= threshold);
  }
}
for (let callWeight = 1; callWeight <= 12; callWeight += 1) {
  const values = [...new Set(predictions.map((row) => callWeight * row.feature.calls - row.feature.vjSpan))].sort((a, b) => a - b);
  for (const threshold of values) add(
    `exact5 ${callWeight}*calls-vjSpan>=${threshold}`,
    (f) => f.exact === 5 && callWeight * f.calls - f.vjSpan >= threshold,
  );
}
for (let alternativeWeight = 1; alternativeWeight <= 8; alternativeWeight += 1) {
  const values = [...new Set(predictions.map((row) => alternativeWeight * row.feature.alternatives - row.feature.vjSpan))].sort((a, b) => a - b);
  for (const threshold of values) add(
    `exact5 ${alternativeWeight}*alternatives-vjSpan>=${threshold}`,
    (f) => f.exact === 5 && alternativeWeight * f.alternatives - f.vjSpan >= threshold,
  );
}

const thresholds = {
  truthTop: Number(process.argv[4] ?? 0),
  truthFair: Number(process.argv[5] ?? 0),
};
const feasible = rules.filter((rule) => rule.truthTop >= thresholds.truthTop && rule.truthFair >= thresholds.truthFair)
  .sort((a, b) => (b.agreementTop + b.agreementFair) - (a.agreementTop + a.agreementFair) || b.truthFair - a.truthFair);

const exactFive = predictions.filter((row) => row.feature.exact === 5 && row.prediction.d_call);
console.log(JSON.stringify({
  records: predictions.length,
  baseline: Object.fromEntries(Object.entries(baseline).map(([key, value]) => [key, value / predictions.length])),
  exactFive: {
    records: exactFive.length,
    truthNegative: exactFive.filter((row) => !atoms(row.call.true_d_call).length).length,
    igblastNegative: exactFive.filter((row) => !atoms(row.call.igblast_d_call).length).length,
  },
  exactFiveByCalls: Object.fromEntries([1, 2, 3, 4, 5, 6].map((count) => {
    const rows = exactFive.filter((row) => row.feature.calls === count || (count === 6 && row.feature.calls >= count));
    return [count === 6 ? "6+" : String(count), {
      records: rows.length,
      truthNegative: rows.filter((row) => !atoms(row.call.true_d_call).length).length,
      igblastNegative: rows.filter((row) => !atoms(row.call.igblast_d_call).length).length,
    }];
  })),
  thresholds,
  feasible: feasible.slice(0, 30),
}, null, 2));
