import assert from "node:assert/strict";
import test from "node:test";

import {
  assignLineages,
  deduplicate,
  expandSingleLinkage,
  minHashSketch,
  prepareReferenceMsa,
  queryRecords,
  runChmm,
  sequenceFingerprint,
  threadSequenceToMsa,
  type PostAnalysisRecord,
} from "../src/post-analysis-core.ts";

function record(
  ordinal: number,
  cdr3Nt: string,
  vCall = "IGHV1-2*01",
  jCall = "IGHJ4*02",
  sequence = `AAA${cdr3Nt}TTT`,
): PostAnalysisRecord {
  return {
    ordinal,
    sequenceId: `read_${ordinal}`,
    locus: "IGH",
    vCall,
    jCall,
    cdr3Nt,
    cdr3Aa: "CARDR",
    productive: true,
    sequenceFingerprint: sequenceFingerprint(sequence),
    trimmedFingerprint: sequenceFingerprint(sequence),
    trimmedSketch: minHashSketch(sequence),
  };
}

test("deduplication retains abundance and representative membership", () => {
  const records = [record(0, "TGTGCCAAA"), record(1, "TGTGCCAAA"), record(2, "TGTGCCTAA")];
  const result = deduplicate(records, "sequence");
  assert.equal(result.uniqueRecords, 2);
  assert.equal(result.collapsedRecords, 1);
  assert.deepEqual([...result.counts], [2, 0, 1]);
  assert.deepEqual([...result.representatives], [0, 0, 2]);
});

test("lineages use V/J ambiguity overlap, exact CDR3 length, and bounded Hamming distance", () => {
  const records = [
    record(0, "TGTGCCAAAA", "IGHV1-2*01,IGHV1-3*01"),
    record(1, "TGTGCCAAAT", "IGHV1-3*02"),
    record(2, "TGTGCCAAATG", "IGHV1-3*02"),
    record(3, "CCCCCCCCCC", "IGHV1-3*02"),
  ];
  const result = assignLineages(records, {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
  });
  assert.equal(result.lineageCount, 3);
  assert.equal(result.assignments[0], result.assignments[1]);
  assert.notEqual(result.assignments[1], result.assignments[2]);
  assert.notEqual(result.assignments[1], result.assignments[3]);
  assert.equal(result.summaries[0].abundance, 2);
});

test("sequence query and single-linkage expansion recover a transitive CDR3 neighborhood", () => {
  const records = [
    record(0, "AAAAAAAAAA"),
    record(1, "AAAAAAAATA"),
    record(2, "AAAAAATATA"),
    record(3, "CCCCCCCCCC"),
  ];
  const hits = queryRecords(records, ["AAAAAAAAAA"], {
    target: "cdr3_nt",
    metric: "hamming",
    identity: 0.9,
    maxResults: 100,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
  });
  assert.deepEqual(hits.map((hit) => hit.ordinal), [0, 1]);
  const expanded = expandSingleLinkage(records, [0], {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
    maxResults: 100,
  });
  assert.deepEqual(expanded.ordinals.sort((a, b) => a - b), [0, 1, 2]);
});

test("CHMMera BW port separates the published chimera and non-chimera example", () => {
  const msa = prepareReferenceMsa(">ref1\nACGTACGTACGT\n>ref2\nACCACCACCAAT\n");
  const options = {
    method: "BW" as const,
    priorProbability: 0.02,
    baseMutationProbability: 0.05,
    mutationRates: [0.005],
    mutationSwitchProbability: 0,
    detailed: true,
  };
  const chimera = runChmm(msa, "ACGTACACCAAT", "ACGTACACCAAT", "ACGTACGTACGT", options);
  const nonChimera = runChmm(msa, "ACCACCACCAGT", "ACCACCACCAGT", "ACCACCACCAAT", options);
  assert.ok(chimera.probability > 0.95, `${chimera.probability}`);
  assert.ok(nonChimera.probability < 0.1, `${nonChimera.probability}`);
  assert.equal(chimera.recombinations.length, 1);
});

test("CHMMera discretized Bayesian mode reports detailed parent switches", () => {
  const msa = prepareReferenceMsa(">left\nAAAAAAAAAAAA\n>right\nCCCCCCCCCCCC\n");
  const result = runChmm(msa, "CCCCCCAAAAAA", "CCCCCCAAAAAA", "AAAAAAAAAAAA", {
    method: "DB",
    priorProbability: 0.05,
    baseMutationProbability: 0.005,
    mutationRates: [0.005],
    mutationSwitchProbability: 0,
    detailed: true,
  });
  assert.ok(result.probability > 0.95, `${result.probability}`);
  assert.ok(result.recombinations.length >= 1);
  assert.notEqual(result.startReference, result.recombinations[0].right);
});

test("AIRR local alignments are threaded onto an uploaded reference MSA", () => {
  const msa = prepareReferenceMsa(">IGHV1*01\nAA--CCGG\n>IGHV2*01\nAATTCCGG\n");
  const threaded = threadSequenceToMsa("TCGG", "CCGG", "IGHV1*01", msa);
  assert.equal(threaded, "NN--TCGG");
});

test("50k post-analysis records retain bounded compact state and avoid all-pairs lineage comparisons", () => {
  let state = 0x12345678;
  const randomNt = () => {
    let value = "";
    for (let index = 0; index < 36; index += 1) {
      state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
      value += "ACGT"[state >>> 30];
    }
    return value;
  };
  const records: PostAnalysisRecord[] = [];
  for (let ordinal = 0; ordinal < 50_000; ordinal += 1) {
    const duplicate = ordinal % 5 === 4;
    const cdr3 = duplicate ? records[ordinal - 1].cdr3Nt : randomNt();
    records.push(record(ordinal, cdr3, `IGHV${ordinal % 20 + 1}*01`, `IGHJ${ordinal % 6 + 1}*01`, `AAA${cdr3}TTT`));
  }
  const collapsed = deduplicate(records, "sequence");
  assert.equal(collapsed.inputRecords, 50_000);
  assert.equal(collapsed.uniqueRecords, 40_000);
  const lineages = assignLineages(records, {
    identity: 0.85,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 50_000,
  }, collapsed);
  assert.equal(lineages.assignedRecords, 50_000);
  assert.equal(lineages.truncatedCandidates, 0);
  assert.ok(lineages.candidateComparisons < 50_000 * 100, `${lineages.candidateComparisons}`);
  assert.ok(lineages.lineageCount <= collapsed.uniqueRecords);
});
