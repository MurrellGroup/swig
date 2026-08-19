import assert from "node:assert/strict";
import test from "node:test";

import { alignedSequenceFrameOffset, inferLineageGermline, lineageInputFasta, quickAirrAlignment } from "../src/lineage-alignment.ts";
import { inferredGermlineIdentity } from "../src/lineage-neighbours.ts";
import { assignLineages, findLineageNeighbours, sequenceFingerprint, type PostAnalysisRecord } from "../src/post-analysis-core.ts";
import type { AirrDetailRow } from "../src/result-store.ts";
import { categoricalValueColor, createSampleColorMap, sampleColor } from "../src/sample-colors.ts";

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

function doubleDDetail(ordinal: number, vIdentity = 0.98, jIdentity = 0.98): AirrDetailRow {
  const row = detail(ordinal, "AAAACCGGGTTCCCAATTTT", "AAAACCCCCCCCCCCCTTTT");
  Object.assign(row.values, {
    v_sequence_start: "1",
    v_sequence_end: "4",
    j_sequence_start: "17",
    j_sequence_end: "20",
    v_sequence_alignment: "AAAA",
    v_germline_alignment: "AAAA",
    j_sequence_alignment: "TTTT",
    j_germline_alignment: "TTTT",
    v_identity: String(vIdentity),
    j_identity: String(jIdentity),
    d_call: "IGHD1-1*01",
    d2_call: "IGHD2-2*01",
    d_sequence_start: "7",
    d_sequence_end: "9",
    d_germline_start: "1",
    d_germline_end: "3",
    d2_sequence_start: "12",
    d2_sequence_end: "14",
    d2_germline_start: "1",
    d2_germline_end: "3",
    d_sequence_alignment: "GGG",
    d_germline_alignment: "GGG",
    d2_sequence_alignment: "CCC",
    d2_germline_alignment: "CCC",
  });
  return row;
}

test("lineage germline defaults to the closest member and keeps consensus as an explicit mode", () => {
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
  assert.equal(inferred.method, "closest");
  assert.equal(inferred.selectedOrdinal, 1);
  assert.equal(inferred.rowsUsed, 1);
  const consensus = inferLineageGermline(rows, "consensus");
  assert.equal(consensus.template, "AAANCCC");
  assert.equal(consensus.uca, "AAATCCC");
  assert.equal(consensus.method, "consensus");
  assert.equal(consensus.rowsUsed, 5);
  assert.equal(consensus.knownColumns, 6);
  assert.equal(consensus.inferredColumns, 1);
});

test("closest-member ranking gives V and J identities equal segment weight", () => {
  const first=detail(0,"AAAANCCC","AAAANCCC");first.values.v_sequence_alignment="AAAAAAAAAA";first.values.v_germline_alignment="AAAAAAAAAA";first.values.j_sequence_alignment="AAAA";first.values.j_germline_alignment="AAAT";
  const second=detail(1,"CCCCNGGG","CCCCNGGG");second.values.v_sequence_alignment="AAAAAAAAAA";second.values.v_germline_alignment="AAAAAAAATA";second.values.j_sequence_alignment="AAAA";second.values.j_germline_alignment="AAAA";
  const inferred=inferLineageGermline([first,second]);
  assert.equal(inferred.selectedOrdinal,1);
  assert.equal(inferred.selectedVIdentity,0.9);
  assert.equal(inferred.selectedJIdentity,1);
  assert.equal(inferred.selectedVjIdentity,0.95);
  assert.equal(inferred.template,"CCCCNGGG");
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

test("reference quick view carries the AIRR biological phase after V-reference padding", () => {
  const row = detail(0, "A---ATGGCC", "A---ATGGCC", 3);
  // One query base precedes the first complete codon, while the AIRR reference
  // anchor contributes two left columns: (1 + 2) mod 3 = column 1/frame 1.
  row.values.v_sequence_start = "1";
  row.values.sequence_frame = "2";
  const input = lineageInputFasta([row]);
  assert.equal(input.frameAnchorUngappedOffset, 1);
  assert.equal(input.alignmentFrameOffset, 0);
  assert.equal(alignedSequenceFrameOffset("A---ATGGCC", 1, 2), 0);
});

test("Double-D lineage roots replace the baseline D with explicit D1 and D2", () => {
  const inferred = inferLineageGermline([doubleDDetail(0)]);
  assert.equal(inferred.template, "AAAANNGGGNNCCCNNTTTT");
  assert.equal(inferred.uca, "AAAACCGGGTTCCCAATTTT");
  assert.equal(inferred.doubleDTemplate, true);
  assert.equal(inferred.doubleDPositiveRows, 1);
  assert.equal(inferred.doubleDResolvedRows, 1);
  assert.equal(inferred.doubleDIncompleteRows, 0);
  assert.equal(inferred.selectedDCall, "IGHD1-1*01");
  assert.equal(inferred.selectedD2Call, "IGHD2-2*01");
  assert.match(quickAirrAlignment([doubleDDetail(0)]), />__germline_N_masked__\nAAAANNGGGNNCCCNNTTTT/);
});

test("Double-D coordinates project through gaps in the combined AIRR query", () => {
  const row = doubleDDetail(0);
  row.values.sequence_alignment = "AA-AACCGGGTTCCCAATTTT";
  row.values.germline_alignment = "AACAACCCCCCCCCCCCTTTT";
  const inferred = inferLineageGermline([row]);
  assert.equal(inferred.template, "AACAANNGGGNNCCCNNTTTT");
  assert.equal(inferred.uca, "AACAACCGGGTTCCCAATTTT");
});

test("Double-D lineages select a VDDJ-aware root member instead of a cleaner single-D composite", () => {
  const singleD = detail(0, "AAAACCGGGTTCCCAATTTT", "AAAANNNGGGNNNNNNTTTT");
  Object.assign(singleD.values, {
    v_sequence_start: "1", v_sequence_end: "4", j_sequence_start: "17", j_sequence_end: "20",
    v_sequence_alignment: "AAAA", v_germline_alignment: "AAAA", j_sequence_alignment: "TTTT", j_germline_alignment: "TTTT",
    v_identity: "1", j_identity: "1",
  });
  const inferred = inferLineageGermline([singleD, doubleDDetail(1, 0.98, 0.98)]);
  assert.equal(inferred.selectedOrdinal, 1);
  assert.equal(inferred.template, "AAAANNGGGNNCCCNNTTTT");
  assert.equal(inferred.doubleDTemplate, true);
  assert.equal(inferred.doubleDRowsUsed, 1);
});

test("Double-D consensus votes only VDDJ-aware templates and leaves ordinary VDJ behavior unchanged", () => {
  const ordinary = detail(0, "AAAACCGGGTTCCCAATTTT", "AAAANNNGGGNNNNNNTTTT");
  Object.assign(ordinary.values, { v_sequence_start: "1", v_sequence_end: "4", j_sequence_start: "17" });
  const doubleDConsensus = inferLineageGermline([ordinary, doubleDDetail(1)], "consensus");
  assert.equal(doubleDConsensus.template, "AAAANNGGGNNCCCNNTTTT");
  assert.equal(doubleDConsensus.rowsUsed, 1);
  assert.equal(doubleDConsensus.doubleDRowsUsed, 1);
  const ordinaryOnly = inferLineageGermline([ordinary], "consensus");
  assert.equal(ordinaryOnly.template, ordinary.values.germline_alignment);
  assert.equal(ordinaryOnly.doubleDTemplate, false);
  assert.equal(ordinaryOnly.rowsUsed, 1);
});

test("an incomplete imported D2 sidecar is reported and never partially projected", () => {
  const incomplete = doubleDDetail(0);
  incomplete.values.d2_germline_alignment = "";
  const inferred = inferLineageGermline([incomplete]);
  assert.equal(inferred.template, incomplete.values.germline_alignment);
  assert.equal(inferred.doubleDTemplate, false);
  assert.equal(inferred.doubleDPositiveRows, 1);
  assert.equal(inferred.doubleDResolvedRows, 0);
  assert.equal(inferred.doubleDIncompleteRows, 1);
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

test("indel-aware neighbours admit CDR3 length variation and can tolerate one noisy germline call", () => {
  const records = [
    record(0, "AAAAAAAAAA", "donor_1", "IGHV1-2*01"),
    record(1, "AAAAAAAAAAA", "donor_1", "IGHV1-2*01"),
    record(2, "AAAAAAAAAT", "donor_1", "IGHV3-7*01"),
  ];
  const assignments = Int32Array.from([1, 2, 3]);
  const base = {
    identity: 0.95,
    minimumIdentity: 0.9,
    callResolution: "gene" as const,
    ambiguity: "overlap" as const,
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
    maximumResults: 20,
    scope: "subject" as const,
    sourceLineageIds: [1],
  };
  const hamming = findLineageNeighbours(records, assignments, { ...base, metric: "hamming", callPolicy: "both" });
  assert.deepEqual(hamming.hits, []);
  const strictEdit = findLineageNeighbours(records, assignments, { ...base, metric: "edit", callPolicy: "both" });
  assert.deepEqual(strictEdit.hits.map((hit) => hit.lineageId), [2]);
  assert.equal(strictEdit.hits[0].cdr3Distance, 1);
  assert.equal(strictEdit.hits[0].cdr3LengthDelta, 1);
  assert.equal(strictEdit.hits[0].callAgreement, "both");
  const tolerant = findLineageNeighbours(records, assignments, { ...base, metric: "edit", callPolicy: "either" });
  assert.deepEqual(tolerant.hits.map((hit) => hit.lineageId).sort((left, right) => left - right), [2, 3]);
  assert.equal(tolerant.hits.find((hit) => hit.lineageId === 3)?.callAgreement, "j");
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
  assert.equal(categoricalValueColor("IGHG1"),categoricalValueColor("IGHG1"));
  assert.notEqual(categoricalValueColor("IGHG1"),categoricalValueColor("IGHM"));
  assert.equal(categoricalValueColor("Unassigned"),"#70817b");
});
