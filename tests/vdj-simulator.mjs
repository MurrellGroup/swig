/**
 * Deterministic V(D)J stress simulator used only by development benchmarks.
 *
 * It is intentionally not a claim to reproduce a repertoire generator. The
 * aim is to exercise caller failure surfaces with plausible combinations of
 * skewed germline usage, exonuclease trimming, P/N junction sequence, SHM
 * hotspots, short indels, sequencing errors, partial reads, and occasional
 * tandem-D rearrangements while retaining exact per-base provenance.
 */

const DNA = "ACGT";
const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A", N: "N" };

export function parseFasta(text) {
  const records = [];
  let name = "";
  let sequence = "";
  const commit = () => {
    const normalized = sequence.toUpperCase().replace(/[^ACGTN]/g, "");
    if (name && normalized) records.push({ name, sequence: normalized });
  };
  for (const line of String(text).split(/\r?\n/)) {
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

export function seededRandom(seed = 1) {
  let state = Number(seed) >>> 0;
  return () => {
    state += 0x6d2b79f5;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 0x1_0000_0000;
  };
}

function integer(random, maximum) {
  return maximum > 0 ? Math.floor(random() * maximum) : 0;
}

function geometric(random, mean, cap) {
  if (mean <= 0 || cap <= 0) return 0;
  const p = 1 / (mean + 1);
  return Math.min(cap, Math.floor(Math.log1p(-random()) / Math.log1p(-p)));
}

function junctionLength(random, mean = 5, cap = 32) {
  if (random() < 0.18) return 0;
  // A two-component tail avoids making every challenge a short junction.
  return geometric(random, random() < 0.92 ? mean : mean * 2.8, cap);
}

function randomBase(random, frequencies = [0.27, 0.23, 0.25, 0.25]) {
  const value = random();
  let cumulative = 0;
  for (let index = 0; index < 4; index += 1) {
    cumulative += frequencies[index];
    if (value <= cumulative) return DNA[index];
  }
  return "T";
}

function randomSequence(random, length) {
  let sequence = "";
  for (let index = 0; index < length; index += 1) sequence += randomBase(random);
  return sequence;
}

function usageSampler(records, random, exponent = 0.72) {
  // Shuffle once, then apply a Zipf-like usage tail. This prevents FASTA order
  // from becoming a hidden truth prior while retaining realistic skew.
  const shuffled = records.map((record) => ({ record, key: random() }))
    .sort((left, right) => left.key - right.key)
    .map(({ record }) => record);
  const cumulative = [];
  let total = 0;
  for (let rank = 0; rank < shuffled.length; rank += 1) {
    total += 1 / Math.pow(rank + 3, exponent);
    cumulative.push(total);
  }
  return () => {
    const target = random() * total;
    let low = 0;
    let high = cumulative.length - 1;
    while (low < high) {
      const middle = (low + high) >>> 1;
      if (cumulative[middle] < target) low = middle + 1;
      else high = middle;
    }
    return shuffled[low];
  };
}

function transition(base, random) {
  if (base === "A") return random() < 0.68 ? "G" : (random() < 0.5 ? "C" : "T");
  if (base === "G") return random() < 0.68 ? "A" : (random() < 0.5 ? "C" : "T");
  if (base === "C") return random() < 0.68 ? "T" : (random() < 0.5 ? "A" : "G");
  if (base === "T") return random() < 0.68 ? "C" : (random() < 0.5 ? "A" : "G");
  return randomBase(random);
}

function isHotspot(sequence, index) {
  const base = sequence[index];
  const left = sequence[index - 1] ?? "N";
  const right = sequence[index + 1] ?? "N";
  // Compact AID-like WRC/GYW approximation; it is sufficient to make errors
  // clustered and context-dependent rather than iid.
  if (base === "C") return "AT".includes(left) && "CT".includes(right);
  if (base === "G") return "AG".includes(left) && "AT".includes(right);
  return false;
}

function cells(label, gene, sequence, referenceStart = 0) {
  return [...sequence].map((base, offset) => ({
    base,
    label,
    gene,
    reference: referenceStart + offset,
    mutated: false,
    inserted: false,
  }));
}

function appendN(target, random, label, length) {
  for (const base of randomSequence(random, length)) {
    target.push({ base, label, gene: "", reference: -1, mutated: false, inserted: false });
  }
}

function addPAndN(target, random, label, length, templateEdge = "") {
  const pLength = templateEdge && random() < 0.22 ? 1 + integer(random, Math.min(3, templateEdge.length)) : 0;
  if (pLength) {
    const p = [...templateEdge.slice(-pLength)].reverse().map((base) => COMPLEMENT[base] ?? "N").join("");
    for (const base of p) target.push({ base, label: `${label}P`, gene: "", reference: -1, mutated: false, inserted: false });
  }
  appendN(target, random, label, length);
  return pLength;
}

function applyShm(cellsValue, random, requestedRate) {
  const original = cellsValue.map((cell) => cell.base).join("");
  let mutations = 0;
  for (let index = 0; index < cellsValue.length; index += 1) {
    const cell = cellsValue[index];
    if (!/^[VDJ]/.test(cell.label) || !DNA.includes(cell.base)) continue;
    const edgeMultiplier = index < 18 || index + 18 >= cellsValue.length ? 0.72 : 1;
    const contextMultiplier = isHotspot(original, index) ? 3.2 : 0.72;
    const probability = Math.min(0.65, requestedRate * edgeMultiplier * contextMultiplier);
    if (random() >= probability) continue;
    cell.base = transition(cell.base, random);
    cell.mutated = true;
    mutations += 1;
  }
  return mutations;
}

function applyIndels(cellsValue, random, rate) {
  if (rate <= 0) return { cells: cellsValue, insertions: 0, deletions: 0 };
  const result = [];
  let insertions = 0;
  let deletions = 0;
  for (let index = 0; index < cellsValue.length; index += 1) {
    const cell = cellsValue[index];
    const templated = /^[VDJ]/.test(cell.label);
    if (templated && random() < rate) {
      const length = 1 + (random() < 0.82 ? 0 : integer(random, 3));
      if (random() < 0.52) {
        deletions += 1;
        index += Math.min(length - 1, cellsValue.length - index - 1);
        continue;
      }
      for (let offset = 0; offset < length; offset += 1) {
        result.push({
          base: randomBase(random), label: "ERRI", gene: "", reference: -1,
          mutated: false, inserted: true,
        });
      }
      insertions += 1;
    }
    result.push(cell);
  }
  return { cells: result, insertions, deletions };
}

function applySequencingNoise(cellsValue, random, substitutionRate, nRate) {
  let substitutions = 0;
  let ambiguous = 0;
  for (const cell of cellsValue) {
    if (random() < nRate) {
      cell.base = "N";
      ambiguous += 1;
    } else if (random() < substitutionRate && DNA.includes(cell.base)) {
      cell.base = transition(cell.base, random);
      substitutions += 1;
    }
  }
  return { substitutions, ambiguous };
}

function provenanceSpan(cellsValue, label) {
  const positions = [];
  for (let index = 0; index < cellsValue.length; index += 1) {
    if (cellsValue[index].label === label) positions.push(index);
  }
  return positions.length ? { start: positions[0], end: positions.at(-1) + 1, bases: positions.length } : null;
}

function reverseComplementCells(cellsValue) {
  return [...cellsValue].reverse().map((cell) => ({ ...cell, base: COMPLEMENT[cell.base] ?? "N" }));
}

export function generateVdjDataset({
  V,
  D,
  J,
  count = 5_000,
  seed = 1,
  doubleDRate = 0.025,
  indelRate = 0.0015,
  sequencingErrorRate = 0.0008,
  ambiguousRate = 0.00025,
  reverseRate = 0.08,
} = {}) {
  if (!V?.length || !J?.length) throw new Error("The simulator requires non-empty V and J references.");
  const random = seededRandom(seed);
  const takeV = usageSampler(V, random, 0.76);
  const takeD = D?.length ? usageSampler(D, random, 0.68) : () => null;
  const takeJ = usageSampler(J, random, 0.60);
  const records = [];

  for (let ordinal = 0; ordinal < count; ordinal += 1) {
    const v = takeV();
    const d1 = takeD();
    const tandem = Boolean(d1 && D.length > 1 && random() < doubleDRate);
    let d2 = tandem ? takeD() : null;
    if (d2?.name === d1?.name && D.length > 1) d2 = D[(D.indexOf(d1) + 1) % D.length];
    const j = takeJ();

    const vTrim3 = geometric(random, 3.2, Math.min(45, Math.max(0, v.sequence.length - 30)));
    const d1Trim5 = d1 ? geometric(random, 2.8, Math.max(0, d1.sequence.length - 1)) : 0;
    const d1Trim3 = d1 ? geometric(random, 3.0, Math.max(0, d1.sequence.length - d1Trim5 - 1)) : 0;
    const d2Trim5 = d2 ? geometric(random, 2.8, Math.max(0, d2.sequence.length - 1)) : 0;
    const d2Trim3 = d2 ? geometric(random, 3.0, Math.max(0, d2.sequence.length - d2Trim5 - 1)) : 0;
    const jTrim5 = geometric(random, 3.4, Math.min(35, Math.max(0, j.sequence.length - 18)));
    const n1Length = junctionLength(random, 5.2);
    const n2Length = junctionLength(random, 4.8);
    const nMiddleLength = tandem ? junctionLength(random, 3.5, 24) : 0;

    const retainedV = v.sequence.slice(0, v.sequence.length - vTrim3);
    const retainedD1 = d1?.sequence.slice(d1Trim5, d1.sequence.length - d1Trim3) ?? "";
    const retainedD2 = d2?.sequence.slice(d2Trim5, d2.sequence.length - d2Trim3) ?? "";
    const retainedJ = j.sequence.slice(jTrim5);
    let rearrangement = cells("V", v.name, retainedV, 0);
    addPAndN(rearrangement, random, "N1", n1Length, retainedV);
    if (d1) rearrangement.push(...cells("D1", d1.name, retainedD1, d1Trim5));
    if (d2) {
      addPAndN(rearrangement, random, "ND", nMiddleLength, retainedD1);
      rearrangement.push(...cells("D2", d2.name, retainedD2, d2Trim5));
    }
    addPAndN(rearrangement, random, "N2", n2Length, retainedD2 || retainedD1);
    rearrangement.push(...cells("J", j.name, retainedJ, jTrim5));

    const shmDraw = random();
    const shmRate = shmDraw < 0.20 ? 0.0025 : shmDraw < 0.50 ? 0.018 :
      shmDraw < 0.78 ? 0.055 : shmDraw < 0.94 ? 0.11 : 0.18;
    const shmMutations = applyShm(rearrangement, random, shmRate);
    const indels = applyIndels(rearrangement, random, indelRate);
    rearrangement = indels.cells;
    const noise = applySequencingNoise(
      rearrangement, random, sequencingErrorRate, ambiguousRate,
    );

    // Partial amplicons and non-reference primer/flank sequence are applied
    // after rearrangement generation so their provenance remains explicit.
    const trim5 = random() < 0.18 ? integer(random, Math.min(90, Math.max(1, retainedV.length - 30))) : 0;
    const trim3 = random() < 0.10 ? integer(random, Math.min(24, Math.max(1, retainedJ.length - 12))) : 0;
    rearrangement = rearrangement.slice(trim5, trim3 ? -trim3 : undefined);
    const flank5 = integer(random, 13);
    const flank3 = integer(random, 9);
    rearrangement = [
      ...cells("F5", "", randomSequence(random, flank5)),
      ...rearrangement,
      ...cells("F3", "", randomSequence(random, flank3)),
    ];
    const reversed = random() < reverseRate;
    if (reversed) rearrangement = reverseComplementCells(rearrangement);

    const spans = {
      V: provenanceSpan(rearrangement, "V"),
      D1: provenanceSpan(rearrangement, "D1"),
      D2: provenanceSpan(rearrangement, "D2"),
      J: provenanceSpan(rearrangement, "J"),
    };
    const sequence = rearrangement.map((cell) => cell.base).join("");
    const id = `sim_${String(ordinal + 1).padStart(7, "0")}`;
    records.push({
      id,
      sequence,
      truth: {
        vCall: v.name,
        dCalls: [d1?.name, d2?.name].filter(Boolean),
        dCompatibleCalls: [retainedD1, retainedD2].map((retained) =>
          retained.length >= 6
            ? D.filter((candidate) => candidate.sequence.includes(retained)).map((candidate) => candidate.name)
            : []),
        jCall: j.name,
        tandem,
        reversed,
        spans,
        retained: { V: retainedV.length, D1: retainedD1.length, D2: retainedD2.length, J: retainedJ.length },
        trims: { v3: vTrim3, d1_5: d1Trim5, d1_3: d1Trim3, d2_5: d2Trim5, d2_3: d2Trim3, j5: jTrim5, read5: trim5, read3: trim3 },
        junction: { n1: n1Length, nd: nMiddleLength, n2: n2Length },
        shmRate,
        shmMutations,
        indelEvents: indels.insertions + indels.deletions,
        sequencingSubstitutions: noise.substitutions,
        ambiguousBases: noise.ambiguous,
      },
    });
  }
  return records;
}

export function recordsToFasta(records) {
  return records.map((record) => `>${record.id}\n${record.sequence}\n`).join("");
}

export function referenceEquivalence(records) {
  const byName = new Map(records.map((record) => [record.name, record.sequence]));
  const namesBySequence = new Map();
  for (const record of records) {
    const names = namesBySequence.get(record.sequence) ?? [];
    names.push(record.name);
    namesBySequence.set(record.sequence, names);
  }
  return { byName, namesBySequence };
}
