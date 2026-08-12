import assert from "node:assert/strict";
import test from "node:test";

import { applyBalancedDFilter, reconcileBalancedDoubleD } from "../src/balanced-calling-profile.ts";

const encoder = new TextEncoder();
const decoder = new TextDecoder();

test("the balanced profile clears only weak D calls in short V-J spans", () => {
  const header = [
    "sequence_id", "sequence", "rev_comp", "d_call", "d_score", "d_identity", "d_cigar",
    "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
    "d_germline_start", "d_germline_end", "j_sequence_start", "cdr3_start",
    "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment",
    "v_germline_alignment", "d_germline_alignment", "j_germline_alignment",
    "sequence_alignment", "sequence_alignment_aa", "germline_alignment", "germline_alignment_aa",
    "np1", "np2", "np1_length", "np2_length", "d_frame", "d_alternatives",
  ].join("\t");
  const row = [
    "q1", "AAACCCCCGGG", "F", "IGHD1*01", "10", "1", "5=", "1", "3", "4", "8",
    "1", "5", "9", "7", "AAA", "CCCCC", "GGG", "AAA", "CCCCC", "GGG",
    "AAACCCCCGGG", "", "AAACCCCCGGG", "", "", "", "0", "0", "1", "IGHD2*01|9|1|4|8|1|5",
  ].join("\t");
  const result = applyBalancedDFilter(header, encoder.encode(`${row}\n`));
  const values = decoder.decode(result.body).trimEnd().split("\t");
  const columns = new Map(header.split("\t").map((name, index) => [name, index]));
  assert.equal(values[columns.get("d_call")!], "");
  assert.equal(values[columns.get("np1")!], "CCCCC");
  assert.equal(values[columns.get("np1_length")!], "5");
  assert.equal(values[columns.get("sequence_alignment")!], "AAACCCCCGGG");
  assert.equal(values[columns.get("germline_alignment")!], "AAANNNNNGGG");
  assert.deepEqual([...result.suppressedSequenceIds], ["q1"]);
});

test("six-base support is retained and double-D standard calls are reconciled", () => {
  const header = "sequence_id\td_call\td_sequence_alignment\td_germline_alignment\tv_sequence_end\tj_sequence_start";
  const row = "q2\tIGHD1*01\tCCCCCC\tCCCCCC\t3\t9\n";
  const result = applyBalancedDFilter(header, encoder.encode(row));
  assert.equal(decoder.decode(result.body), row);
  assert.equal(result.suppressedSequenceIds.size, 0);

  const doubleDHeader = "sequence_id\tstandard_d_call\td_call\td2_call";
  const reconciled = reconcileBalancedDoubleD(
    doubleDHeader,
    encoder.encode("q1\tIGHD1*01\tIGHD2*01\tIGHD3*01\n"),
    new Set(["q1"]),
  );
  assert.equal(decoder.decode(reconciled), "q1\t\tIGHD2*01\tIGHD3*01\n");
});
