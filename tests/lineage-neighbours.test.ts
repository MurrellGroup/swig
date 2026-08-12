import assert from "node:assert/strict";
import test from "node:test";

import { inferLineageGermline, quickAirrAlignment } from "../src/lineage-alignment.ts";
import { inferredGermlineIdentity } from "../src/lineage-neighbours.ts";
import { assignLineages, findLineageNeighbours, sequenceFingerprint, type PostAnalysisRecord } from "../src/post-analysis-core.ts";
import type { AirrDetailRow } from "../src/result-store.ts";
import { createSampleColorMap, sampleColor } from "../src/sample-colors.ts";

function detail(ordinal: number, sequence: string, germline: string, start = 1): AirrDetailRow {
  return {
    record: { ordinal, sequenceId: `read_${ordinal}`, sampleId: "sample_A", locus: "IGH" },
    values: {
      sequence_id: `read_${ordinal}`,
      sequence_alignment: sequence,
      germline_alignment: germline,
      v_germline_start: String(start),
      v_sequence_start: "1",
      sequence_frame: "1",
      locus: "IGH",
      sample_id: "sample_A",
    },
  } as AirrDetailRow;
}

test("lineage germline uses member votes instead of choosing the first row", () => {
  const rows = [
    detail(0, "AAATCCC", "CAANCCC"),
    detail(1, "AAATCCC", "AAANCCC"),
    detail(2, "AAATCCC", "AAANCCC"),
    detail(3, "AAATCCC", "AAANCCC"),
    detail(4, "AAACCCC", "AAANCCC"),
  ];
  const inferred = inferLineageGermline(rows);
  assert.equal(inferred.template, "AAANCCC");
  assert.equal(inferred.uca, "AAATCCC");
  assert.equal(inferred.rowsUsed, 5);
  assert.equal(inferred.knownColumns, 6);
  assert.equal(inferred.inferredColumns, 1);
});

test("reference quick view anchors rows on the V germline coordinate", () => {
  const fasta = quickAirrAlignment([
    detail(0, "ACGT", "ACGT", 1),
    detail(1, "CGT", "CGT", 2),
  ]);
  assert.match(fasta, />read_0__1\nACGT/);
  assert.match(fasta, />read_1__2\n-CGT/);
  assert.match(fasta, />__germline_N_masked__\nACGT/);
});

function record(ordinal: number, cdr3Nt: string, subjectId: string, vCall = "IGHV1-2*01"): PostAnalysisRecord {
  const sequence = `AAA${cdr3Nt}TTT`;
  return {
    ordinal,
    sequenceId: `read_${ordinal}`,
    sampleId: `${subjectId}_sample`,
    subjectId,
    locus: "IGH",
    vCall,
    jCall: "IGHJ4*02",
    cdr3Nt,
    cdr3Aa: "CARDR",
    productive: true,
    sequenceFingerprint: sequenceFingerprint(sequence),
    trimmedFingerprint: sequenceFingerprint(sequence),
  };
}

test("CDR3 neighbour search crosses assignment boundaries but preserves V/J and study constraints", () => {
  const records = [
    record(0, "AAAAAAAAAA", "donor_1"),
    record(1, "AAAAAAAACC", "donor_1"),
    record(2, "AAAAAAAACC", "donor_1", "IGHV3-7*01"),
    record(3, "AAAAAAAACC", "donor_2"),
  ];
  const options = { identity: 0.9, callResolution: "gene" as const, ambiguity: "overlap" as const, productiveOnly: true, requireSameLocus: true, maxCandidateComparisons: 10_000, scope: "subject" as const };
  const assigned = assignLineages(records, options);
  assert.notEqual(assigned.assignments[0], assigned.assignments[1]);
  const result = findLineageNeighbours(records, assigned.assignments, {
    ...options,
    sourceLineageIds: [assigned.assignments[0]],
    minimumIdentity: 0.8,
    maximumResults: 20,
  });
  assert.deepEqual(result.hits.map((hit) => hit.lineageId), [assigned.assignments[1]]);
  assert.equal(result.hits[0].cdr3Identity, 0.8);
});

test("inferred germline identity treats N as unknown while retaining indel cost", () => {
  assert.equal(inferredGermlineIdentity("ACNGT", "ACAGT", 0.9), 1);
  assert.equal(inferredGermlineIdentity("ACGT", "ACGGT", 0), 0.8);
});

test("sample colors are stable, editable associations", () => {
  const datasets = [
    { datasetId: "d1", inputName: "a.tsv", sampleId: "sample_A", subjectId: "donor", cohort: "case", timepoint: "d0" },
    { datasetId: "d2", inputName: "b.tsv", sampleId: "sample_B", subjectId: "donor", cohort: "case", timepoint: "d7" },
  ];
  const defaults = createSampleColorMap(datasets);
  const edited = createSampleColorMap(datasets, { ...defaults, sample_A: "#123456" });
  assert.equal(sampleColor("sample_A", edited), "#123456");
  assert.equal(sampleColor("sample_B", edited), defaults.sample_B);
});
