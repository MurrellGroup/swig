import fs from "node:fs";

function argsMap(argv) {
  const result = new Map();
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (!key.startsWith("--")) continue;
    result.set(key.slice(2), argv[++index] ?? "");
  }
  return result;
}

function table(path) {
  const lines = fs.readFileSync(path, "utf8").trimEnd().split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const fields = line.split("\t");
    return Object.fromEntries(header.map((name, index) => [name, fields[index] ?? ""]));
  });
}

function predictions(paths) {
  const rows = new Map();
  for (const path of paths.split(",").filter(Boolean)) {
    for (const row of table(path)) {
      if (rows.has(row.sequence_id)) throw new Error(`Duplicate prediction ${row.sequence_id}`);
      rows.set(row.sequence_id, row);
    }
  }
  return rows;
}

function atoms(value) {
  return [...new Set(String(value ?? "").split(/[,/]/).map((item) => item.trim()).filter(Boolean))];
}

function setEqual(a, b) {
  return a.length === b.length && a.every((value) => b.includes(value));
}

function comparison(call, reference) {
  const a = atoms(call);
  const b = atoms(reference);
  const overlap = a.filter((value) => b.includes(value)).length;
  const bothEmpty = !a.length && !b.length;
  const union = new Set([...a, ...b]).size;
  return {
    called: a.length > 0,
    referenceCalled: b.length > 0,
    firstExact: bothEmpty || (a.length > 0 && b.length > 0 && a[0] === b[0]),
    firstSupported: bothEmpty || (a.length > 0 && b.includes(a[0])),
    compatible: bothEmpty || overlap > 0,
    exactSet: setEqual(a, b),
    precision: bothEmpty ? 1 : a.length ? overlap / a.length : 0,
    recall: bothEmpty ? 1 : b.length ? overlap / b.length : 0,
    jaccard: bothEmpty ? 1 : union ? overlap / union : 0,
    callStatus: (a.length > 0) === (b.length > 0),
    atomCount: a.length,
  };
}

function emptyState() {
  return {
    records: 0, firstExact: 0, firstSupported: 0, compatible: 0, exactSet: 0,
    precision: 0, recall: 0, jaccard: 0, callStatus: 0, called: 0,
    referenceCalled: 0, atoms: 0,
  };
}

function add(state, result) {
  state.records += 1;
  for (const key of ["firstExact", "firstSupported", "compatible", "exactSet", "callStatus", "called", "referenceCalled"]) {
    state[key] += Number(result[key]);
  }
  for (const key of ["precision", "recall", "jaccard"]) state[key] += result[key];
  state.atoms += result.atomCount;
}

function summarize(state) {
  const fraction = (value) => value / state.records;
  return {
    records: state.records,
    firstExact: fraction(state.firstExact),
    firstSupported: fraction(state.firstSupported),
    compatible: fraction(state.compatible),
    exactSet: fraction(state.exactSet),
    directionalFair: fraction(state.precision),
    reciprocalCoverage: fraction(state.recall),
    jaccard: fraction(state.jaccard),
    callStatusAgreement: fraction(state.callStatus),
    callRate: fraction(state.called),
    referenceCallRate: fraction(state.referenceCalled),
    meanAtoms: fraction(state.atoms),
  };
}

function transitionState() {
  return {
    records: 0, changedSet: 0, changedFirst: 0, changedStatus: 0,
    firstToward: 0, firstAway: 0, compatibleToward: 0, compatibleAway: 0,
    exactToward: 0, exactAway: 0, fairHigher: 0, fairLower: 0,
    fairDelta: 0, jaccardDelta: 0,
    calledToEmptyIgEmpty: 0, calledToEmptyIgCalled: 0,
    emptyToCalledIgCalled: 0, emptyToCalledIgEmpty: 0,
  };
}

const args = argsMap(process.argv.slice(2));
for (const required of ["calls", "old", "new"]) {
  if (!args.get(required)) throw new Error(`Missing --${required}`);
}

const calls = table(args.get("calls"));
const oldRows = predictions(args.get("old"));
const newRows = predictions(args.get("new"));
if (oldRows.size !== calls.length || newRows.size !== calls.length) {
  throw new Error(`Coverage mismatch: calls=${calls.length}, old=${oldRows.size}, new=${newRows.size}`);
}

const output = {};
for (const segment of ["v", "d", "j"]) {
  const oldState = emptyState();
  const newState = emptyState();
  const oldReferenceCalled = emptyState();
  const newReferenceCalled = emptyState();
  const oldReferenceEmpty = emptyState();
  const newReferenceEmpty = emptyState();
  const oldBothCalled = emptyState();
  const newBothCalled = emptyState();
  const oldCommonCalled = emptyState();
  const newCommonCalled = emptyState();
  const transitions = transitionState();
  for (const row of calls) {
    const id = row.sequence_id;
    const oldCall = oldRows.get(id)?.[`${segment}_call`] ?? "";
    const newCall = newRows.get(id)?.[`${segment}_call`] ?? "";
    const reference = row[`igblast_${segment}_call`] ?? "";
    const oldResult = comparison(oldCall, reference);
    const newResult = comparison(newCall, reference);
    add(oldState, oldResult);
    add(newState, newResult);
    if (oldResult.referenceCalled) {
      add(oldReferenceCalled, oldResult);
      add(newReferenceCalled, newResult);
      if (oldResult.called) add(oldBothCalled, oldResult);
      if (newResult.called) add(newBothCalled, newResult);
      if (oldResult.called && newResult.called) {
        add(oldCommonCalled, oldResult);
        add(newCommonCalled, newResult);
      }
    } else {
      add(oldReferenceEmpty, oldResult);
      add(newReferenceEmpty, newResult);
    }
    transitions.records += 1;
    if (!setEqual(atoms(oldCall), atoms(newCall))) transitions.changedSet += 1;
    if ((atoms(oldCall)[0] ?? "") !== (atoms(newCall)[0] ?? "")) transitions.changedFirst += 1;
    if (oldResult.called !== newResult.called) transitions.changedStatus += 1;
    if (!oldResult.firstSupported && newResult.firstSupported) transitions.firstToward += 1;
    if (oldResult.firstSupported && !newResult.firstSupported) transitions.firstAway += 1;
    if (!oldResult.compatible && newResult.compatible) transitions.compatibleToward += 1;
    if (oldResult.compatible && !newResult.compatible) transitions.compatibleAway += 1;
    if (!oldResult.exactSet && newResult.exactSet) transitions.exactToward += 1;
    if (oldResult.exactSet && !newResult.exactSet) transitions.exactAway += 1;
    const fairDelta = newResult.precision - oldResult.precision;
    if (fairDelta > 1e-12) transitions.fairHigher += 1;
    if (fairDelta < -1e-12) transitions.fairLower += 1;
    transitions.fairDelta += fairDelta;
    transitions.jaccardDelta += newResult.jaccard - oldResult.jaccard;
    if (oldResult.called && !newResult.called) {
      if (atoms(reference).length) transitions.calledToEmptyIgCalled += 1;
      else transitions.calledToEmptyIgEmpty += 1;
    }
    if (!oldResult.called && newResult.called) {
      if (atoms(reference).length) transitions.emptyToCalledIgCalled += 1;
      else transitions.emptyToCalledIgEmpty += 1;
    }
  }
  output[segment.toUpperCase()] = {
    old: summarize(oldState),
    new: summarize(newState),
    strata: {
      referenceCalled: {
        old: summarize(oldReferenceCalled),
        new: summarize(newReferenceCalled),
      },
      referenceEmpty: oldReferenceEmpty.records ? {
        old: summarize(oldReferenceEmpty),
        new: summarize(newReferenceEmpty),
      } : null,
      bothCalled: {
        old: summarize(oldBothCalled),
        new: summarize(newBothCalled),
      },
      commonCalled: {
        old: summarize(oldCommonCalled),
        new: summarize(newCommonCalled),
      },
    },
    transitions: {
      ...transitions,
      changedSetRate: transitions.changedSet / transitions.records,
      changedFirstRate: transitions.changedFirst / transitions.records,
      changedStatusRate: transitions.changedStatus / transitions.records,
      meanFairDelta: transitions.fairDelta / transitions.records,
      meanJaccardDelta: transitions.jaccardDelta / transitions.records,
    },
  };
}

console.log(JSON.stringify(output, null, 2));
