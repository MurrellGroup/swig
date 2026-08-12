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

function score(referenceValue, predictionValue) {
  const reference = new Set(atoms(referenceValue));
  const prediction = atoms(predictionValue);
  if (!reference.size) return { top: prediction.length ? 0 : 1, fair: prediction.length ? 0 : 1 };
  return {
    top: prediction.length && reference.has(prediction[0]) ? 1 : 0,
    fair: prediction.length
      ? new Set(prediction.filter((call) => reference.has(call))).size / prediction.length
      : 0,
  };
}

function agreementScore(referenceValue, predictionValue) {
  const reference = atoms(referenceValue);
  const prediction = atoms(predictionValue);
  const bothEmpty = !reference.length && !prediction.length;
  const overlap = prediction.filter((call) => reference.includes(call)).length;
  return {
    top: bothEmpty || (reference.length && prediction.length && reference[0] === prediction[0]) ? 1 : 0,
    fair: bothEmpty ? 1 : prediction.length ? overlap / prediction.length : 0,
  };
}

function longestExact(query, germline) {
  let current = 0;
  let longest = 0;
  for (let index = 0; index < Math.min(query.length, germline.length); index += 1) {
    if (query[index] !== "-" && query[index] === germline[index]) longest = Math.max(longest, ++current);
    else current = 0;
  }
  return longest;
}

const truthPath = process.argv[2];
const predictionPath = process.argv[3];
const requiredTruthTop = Number(process.argv[4] ?? 0);
const requiredTruthFair = Number(process.argv[5] ?? 0);
if (!truthPath || !predictionPath) {
  throw new Error("Usage: node select-balanced-d-rule.mjs truth.tsv predictions.tsv requiredTop requiredFair");
}

const truth = new Map(table(truthPath).map((row) => [row.sequence_id, row]));
const predictions = table(predictionPath);
if (truth.size !== predictions.length) throw new Error(`Coverage mismatch: ${truth.size} vs ${predictions.length}`);

const rows = predictions.map((prediction) => {
  const expected = truth.get(prediction.sequence_id);
  if (!expected) throw new Error(`Missing truth row: ${prediction.sequence_id}`);
  const query = prediction.d_sequence_alignment ?? "";
  const germline = prediction.d_germline_alignment ?? "";
  const originalTruth = score(expected.true_d_call, prediction.d_call);
  const originalAgreement = agreementScore(expected.igblast_d_call, prediction.d_call);
  const emptyTruth = score(expected.true_d_call, "");
  const emptyAgreement = agreementScore(expected.igblast_d_call, "");
  return {
    exact: longestExact(query, germline),
    calls: atoms(prediction.d_call).length,
    alternatives: String(prediction.d_alternatives ?? "").split(";").filter(Boolean).length,
    vjSpan: Math.max(0, Number(prediction.j_sequence_start || 0) - Number(prediction.v_sequence_end || 0)),
    dTruthTop: emptyTruth.top - originalTruth.top,
    dTruthFair: emptyTruth.fair - originalTruth.fair,
    dAgreementTop: emptyAgreement.top - originalAgreement.top,
    dAgreementFair: emptyAgreement.fair - originalAgreement.fair,
    originalTruth,
    originalAgreement,
    truthNegative: !atoms(expected.true_d_call).length,
    igblastNegative: !atoms(expected.igblast_d_call).length,
  };
});

const records = rows.length;
const baseline = rows.reduce((total, row) => {
  total.truthTop += row.originalTruth.top;
  total.truthFair += row.originalTruth.fair;
  total.agreementTop += row.originalAgreement.top;
  total.agreementFair += row.originalAgreement.fair;
  return total;
}, { truthTop: 0, truthFair: 0, agreementTop: 0, agreementFair: 0 });

const weak = rows.filter((row) => row.exact === 5 && row.calls);
const candidates = [];
const add = (name, predicate) => candidates.push({ name, predicate });
for (let span = 0; span <= 40; span += 1) add(`exact5 vjSpan<=${span}`, (row) => row.vjSpan <= span);
for (let weight = 1; weight <= 12; weight += 1) {
  for (let threshold = -30; threshold <= 50; threshold += 1) {
    add(`exact5 ${weight}*calls-vjSpan>=${threshold}`, (row) => weight * row.calls - row.vjSpan >= threshold);
  }
}
for (let weight = 1; weight <= 8; weight += 1) {
  for (let threshold = -30; threshold <= 50; threshold += 1) {
    add(`exact5 ${weight}*alternatives-vjSpan>=${threshold}`, (row) => weight * row.alternatives - row.vjSpan >= threshold);
  }
}

const evaluated = candidates.map(({ name, predicate }) => {
  const totals = { ...baseline };
  let dropped = 0;
  let truthNegative = 0;
  let igblastNegative = 0;
  for (const row of weak) {
    if (!predicate(row)) continue;
    dropped += 1;
    truthNegative += row.truthNegative ? 1 : 0;
    igblastNegative += row.igblastNegative ? 1 : 0;
    totals.truthTop += row.dTruthTop;
    totals.truthFair += row.dTruthFair;
    totals.agreementTop += row.dAgreementTop;
    totals.agreementFair += row.dAgreementFair;
  }
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
});

const feasible = evaluated
  .filter((result) => result.truthTop >= requiredTruthTop && result.truthFair >= requiredTruthFair)
  .sort((left, right) =>
    (right.agreementTop + right.agreementFair) - (left.agreementTop + left.agreementFair)
    || (right.truthTop + right.truthFair) - (left.truthTop + left.truthFair));

console.log(JSON.stringify({
  records,
  weakExactFive: weak.length,
  thresholds: { truthTop: requiredTruthTop, truthFair: requiredTruthFair },
  baseline: Object.fromEntries(Object.entries(baseline).map(([key, value]) => [key, value / records])),
  spanRules: evaluated.filter((result) => /^exact5 vjSpan<=/.test(result.name)),
  feasible: feasible.slice(0, 40),
}, null, 2));
