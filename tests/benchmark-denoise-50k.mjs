import { performance } from "node:perf_hooks";
import { getHeapStatistics } from "node:v8";

import { DenoiseAccumulator, sequenceFingerprint } from "../src/post-analysis-core.ts";

const TOTAL = 50_000;
const PARTITIONS = 100;
const CENTROIDS_PER_PARTITION = 10;
const READS_PER_CENTROID = TOTAL / PARTITIONS / CENTROIDS_PER_PARTITION;

function sequence(seed, length = 120) {
  let state = seed >>> 0;
  let value = "";
  for (let index = 0; index < length; index += 1) {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    value += "ACGT"[state >>> 30];
  }
  return value;
}

const records = [];
const sequences = [];
let ordinal = 0;
for (let partition = 0; partition < PARTITIONS; partition += 1) {
  for (let centroid = 0; centroid < CENTROIDS_PER_PARTITION; centroid += 1) {
    const parent = sequence(partition * 10_000 + centroid * 101 + 17);
    const position = 30 + centroid;
    const replacement = parent[position] === "A" ? "C" : "A";
    const child = `${parent.slice(0, position)}${replacement}${parent.slice(position + 1)}`;
    for (let read = 0; read < READS_PER_CENTROID; read += 1) {
      const value = read === READS_PER_CENTROID - 1 ? child : parent;
      records.push({
        ordinal,
        sequenceId: `read_${ordinal}`,
        locus: "IGH",
        vCall: `IGHV${Math.floor(partition / 10) + 1}-1*01`,
        jCall: `IGHJ${partition % 10 + 1}*01`,
        cdr3Nt: value.slice(45, 75),
        cdr3Aa: "",
        productive: true,
        sequenceFingerprint: sequenceFingerprint(value),
        trimmedFingerprint: sequenceFingerprint(value),
      });
      sequences.push(value);
      ordinal += 1;
    }
  }
}

const before = getHeapStatistics().used_heap_size;
const started = performance.now();
const accumulator = new DenoiseAccumulator(records, {
  mode: "fad",
  errorRate: 0.00473,
  alpha: 0.01,
  callResolution: "allele",
  ambiguity: "strict",
  minimumParentCount: 2,
  ambiguousPolicy: "exclude",
  fadNeighborThreshold: 1,
  fadMethod: 2,
  expectedZeroErrorFraction: 1,
  maximumHammingDistance: 1,
  maxCandidatesPerVariant: 50_000,
});
for (let index = 0; index < sequences.length; index += 1) accumulator.add(index, sequences[index]);
const result = accumulator.finish();
const elapsedMs = performance.now() - started;
const after = getHeapStatistics().used_heap_size;

console.log(JSON.stringify({
  records: TOTAL,
  retained: result.uniqueRecords,
  partitions: result.partitions,
  verifiedCandidates: result.candidateComparisons,
  elapsedMs: Math.round(elapsedMs),
  recordsPerSecond: Math.round(TOTAL / (elapsedMs / 1_000)),
  heapDeltaMiB: Math.round((after - before) / 1024 / 1024),
}, null, 2));
