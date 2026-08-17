import assert from "node:assert/strict";
import test from "node:test";

import { biologicalFrameOffset, biologicalSegmentAlignment, buildTrackFeatures, projectCodonAlignment, translateCodonSequence } from "../src/alignment-model.ts";
import { isProductiveLineageRow } from "../src/lineage-alignment.ts";
import { alignmentRetainedRows } from "../src/lineage-phylogeny.ts";
import type { AirrDetailRow } from "../src/result-store.ts";

test("the layered track model preserves AIRR region and germline coordinates", () => {
  const row = {
    sequence: "A".repeat(300),
    fwr1_start: "1", fwr1_end: "78",
    cdr1_start: "79", cdr1_end: "114",
    fwr2_start: "115", fwr2_end: "165",
    cdr2_start: "166", cdr2_end: "195",
    fwr3_start: "196", fwr3_end: "252",
    cdr3_start: "253", cdr3_end: "273",
    fwr4_start: "274", fwr4_end: "300",
    v_call: "IGHV1-1*01", v_sequence_start: "1", v_sequence_end: "252",
    d_call: "IGHD1-1*01", d_sequence_start: "255", d_sequence_end: "263",
    j_call: "IGHJ4*01", j_sequence_start: "267", j_sequence_end: "300",
  };
  const tracks = buildTrackFeatures(row);
  assert.deepEqual(tracks.regions.map((feature) => feature.key), ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]);
  assert.deepEqual(tracks.segments.map((feature) => feature.key), ["v", "d", "j"]);
  assert.equal(tracks.regions[0].left, 0);
  assert.equal(tracks.regions.at(-1)?.end, 300);
  assert.equal(tracks.segments[1].start, 255);
});

test("per-segment amino-acid alignment uses the rearrangement frame rather than optimizing a new frame", () => {
  const translated = biologicalSegmentAlignment("AAATGGGCT", "AAATGGGCT", 2, 1);
  assert.deepEqual(translated, { query: "MG", reference: "MG" });
  assert.equal(biologicalSegmentAlignment("AAATGGGCT", "AAATGGGCT", 2, 0), null);
  assert.equal(biologicalFrameOffset(5, 2), 0);
  assert.equal(biologicalFrameOffset(6, 2), 2);
});

test("codon projection preserves phase and emits only complete triplet columns", () => {
  const projected = projectCodonAlignment("AATGGGC", "-MG", 1);
  assert.equal(projected, "A-----ATGGGC");
  assert.equal(projected.length % 3, 0);
});

test("AA-guided MSA translation uses X for ambiguous codons and back-translates gaps as codon triplets", () => {
  assert.equal(translateCodonSequence("AATNNNTGC", 0), "NXC");
  assert.equal(translateCodonSequence("AATGGGC", 1), "MG");
  const dna = "ATGGGGCCC";
  const projected = projectCodonAlignment(dna, "M-GP", 0);
  assert.equal(projected, "ATG---GGGCCC");
  assert.equal(projected.replaceAll("-", ""), dna);
  assert.equal(projected.length % 3, 0);
  assert.equal(projected.slice(3, 6), "---");
  assert.doesNotMatch(projected, /^(?:---)|(?:---)$/);
});

test("the lineage productivity filter requires an explicit AIRR productive value", () => {
  const row = (value: string, fallback = "") => ({ values: { productive: value }, record: { productive: fallback } }) as unknown as AirrDetailRow;
  assert.equal(isProductiveLineageRow(row("T")), true);
  assert.equal(isProductiveLineageRow(row("true")), true);
  assert.equal(isProductiveLineageRow(row("", "T")), true);
  assert.equal(isProductiveLineageRow(row("F")), false);
  assert.equal(isProductiveLineageRow(row("")), false);
});

test("AIRR metadata for a row deleted from the MSA cannot leak into UCA preparation", () => {
  const rows = [4, 8, 12].map((ordinal) => ({ record: { ordinal }, values: {} })) as unknown as AirrDetailRow[];
  const retained = alignmentRetainedRows([
    { name: "member_a__5", sequence: "ACGT" },
    { name: "member_c__13", sequence: "ACGT" },
    { name: "__germline_N_masked__", sequence: "NNNN" },
  ], rows);
  assert.deepEqual(retained.map((row) => row.record.ordinal), [4, 12]);
});
