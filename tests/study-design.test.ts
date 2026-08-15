import test from "node:test";
import assert from "node:assert/strict";

import { assignLineages, deduplicate, sequenceFingerprint, type PostAnalysisRecord } from "../src/post-analysis-core.ts";
import { annotateAirrBatch, annotateDoubleDBatch, datasetScopeKey, DEFAULT_PIPELINE_PLAN, type DatasetManifestEntry } from "../src/study-design.ts";

const manifest: DatasetManifestEntry = {
  datasetId: "dataset_2",
  inputName: "day7.fastq.gz",
  sampleId: "sample_day7",
  subjectId: "donor_A",
  cohort: "vaccinated",
  timepoint: "day_7",
  compartment: "blood",
};

function record(ordinal: number, sampleId: string, subjectId: string, sequence = "ACGTACGT"): PostAnalysisRecord {
  return {
    ordinal,
    sequenceId: `read_${ordinal}`,
    datasetId: `dataset_${ordinal + 1}`,
    sampleId,
    subjectId,
    cohort: "cohort_1",
    timepoint: `t${ordinal}`,
    compartment: ordinal ? "lymph_node" : "blood",
    locus: "IGH",
    vCall: "IGHV1-1*01",
    jCall: "IGHJ4*01",
    cdr3Nt: "TGTGCCGAA",
    cdr3Aa: "CAE",
    productive: true,
    sequenceFingerprint: sequenceFingerprint(sequence),
    trimmedFingerprint: sequenceFingerprint(sequence),
  };
}

test("AIRR batches receive collision-safe IDs and explicit study metadata", () => {
  const result = annotateAirrBatch("sequence_id\tsequence\tsample_id", "read_1\tACGT\told_sample\n", manifest);
  const headers = result.header.split("\t");
  const values = result.body.trimEnd().split("\t");
  const row = Object.fromEntries(headers.map((field, index) => [field, values[index]]));
  assert.equal(row.sequence_id, "dataset_2::read_1");
  assert.equal(row.swig_source_sequence_id, "read_1");
  assert.equal(row.sample_id, "sample_day7");
  assert.equal(row.subject_id, "donor_A");
  assert.equal(row.swig_cohort, "vaccinated");
  assert.equal(row.swig_timepoint, "day_7");
  assert.equal(row.swig_compartment, "blood");
  assert.equal(headers.length, values.length);
});

test("double-D batches retain the relative record index while receiving study metadata", () => {
  const result = annotateDoubleDBatch("swig_batch_record_index\tsequence_id\td_call\td2_call", "0\tread_1\tD1\tD2\n", manifest);
  const headers = result.header.split("\t");
  const values = result.body.trimEnd().split("\t");
  const row = Object.fromEntries(headers.map((field, index) => [field, values[index]]));
  assert.equal(row.swig_batch_record_index, "0");
  assert.equal(row.sequence_id, "dataset_2::read_1");
  assert.equal(row.sample_id, "sample_day7");
  assert.equal(row.swig_compartment, "blood");
});

test("collapse scope treats technical libraries as one sample without crossing samples", () => {
  const records = [record(0, "sample_A", "donor_A"), record(1, "sample_A", "donor_A"), record(2, "sample_B", "donor_A")];
  const sampleScoped = deduplicate(records, "sequence", "discard", "sample");
  assert.equal(sampleScoped.uniqueRecords, 2);
  assert.equal(sampleScoped.counts[0], 2);
  assert.equal(sampleScoped.counts[2], 1);
  const datasetScoped = deduplicate(records, "sequence", "discard", "dataset");
  assert.equal(datasetScoped.uniqueRecords, 3);
});

test("lineages can span longitudinal samples within a donor but never cross donors", () => {
  const records = [record(0, "day_0", "donor_A"), record(1, "day_30", "donor_A"), record(2, "day_0", "donor_B")];
  const donorScoped = assignLineages(records, {
    identity: 1,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 1000,
    scope: "subject",
  });
  assert.equal(donorScoped.lineageCount, 2);
  assert.equal(donorScoped.assignments[0], donorScoped.assignments[1]);
  assert.notEqual(donorScoped.assignments[0], donorScoped.assignments[2]);
  assert.equal(datasetScopeKey(records[0], "subject"), datasetScopeKey(records[1], "subject"));
  const shared=donorScoped.summaries.find((summary)=>summary.id===donorScoped.assignments[0]);
  assert.deepEqual(shared?.sampleIds,["day_0","day_30"]);
  assert.deepEqual(shared?.compartments,["blood","lymph_node"]);
});

test("longitudinal/compartmental defaults keep collapse within sample and lineages within donor",()=>{
  assert.equal(DEFAULT_PIPELINE_PLAN.collapse.scope,"sample");
  assert.equal(DEFAULT_PIPELINE_PLAN.collapse.key,"trimmed");
  assert.equal(DEFAULT_PIPELINE_PLAN.collapse.respectConstantCall,true);
  assert.equal(DEFAULT_PIPELINE_PLAN.lineage.scope,"subject");
  assert.equal(DEFAULT_PIPELINE_PLAN.alleleRefinement.reassignmentPolicy,"confidence");
  assert.equal(DEFAULT_PIPELINE_PLAN.alleleRefinement.applyMinimumPosterior,0.8);
  assert.equal("excludedAlleles" in DEFAULT_PIPELINE_PLAN.chimera,false);
});
