import assert from "node:assert/strict";
import test from "node:test";

import { inspectAlignment, validateCorrectedAlignment } from "../src/alignment-provenance.ts";
import { prepareFastTreeInput } from "../src/fasttree-input.ts";

const ORIGINAL = [
  ">read_a",
  "AC-GT",
  ">read_b",
  "A-CGT",
  ">__germline_N_masked__",
  "ACG-T",
  "",
].join("\n");

test("a complete gap-only alignment correction is accepted and fingerprinted", () => {
  const corrected = [
    ">read_a",
    "A-CGT",
    ">read_b",
    "AC-GT",
    ">__germline_N_masked__",
    "ACGT-",
    "",
  ].join("\n");
  const before = inspectAlignment(ORIGINAL);
  const after = validateCorrectedAlignment(ORIGINAL, corrected);
  assert.equal(after.rows, 3);
  assert.equal(after.columns, 5);
  assert.notEqual(after.fingerprint, before.fingerprint);
  assert.equal(after.fasta, corrected);
});

test("Alivibe correction may delete bad rows and nucleotide characters", () => {
  const selectedRows = [
    ">read_a",
    "AC-GT",
    ">__germline_N_masked__",
    "ACG-T",
    "",
  ].join("\n");
  const rowDeletion=validateCorrectedAlignment(ORIGINAL, selectedRows);
  assert.deepEqual(rowDeletion.removedRows,["read_b"]);

  const selectedColumns = [
    ">read_a",
    "AC-G",
    ">read_b",
    "A-CG",
    ">__germline_N_masked__",
    "ACG-",
    "",
  ].join("\n");
  const columnDeletion=validateCorrectedAlignment(ORIGINAL, selectedColumns);
  assert.equal(columnDeletion.removedNucleotides,3);
});

test("manual alignment import cannot mutate nucleotides or duplicate identifiers", () => {
  const mutated = ORIGINAL.replace("AC-GT", "AT-GT");
  assert.throws(() => validateCorrectedAlignment(ORIGINAL, mutated), /substitution/);
  const missingRoot=">read_a\nACGT\n>read_b\nACGT\n";
  assert.throws(()=>validateCorrectedAlignment(ORIGINAL,missingRoot),/must retain __germline_N_masked__/);
  const added=`${ORIGINAL}>new_row\nAC-GT\n`;
  assert.throws(()=>validateCorrectedAlignment(ORIGINAL,added),/unexpected or renamed rows/);
  assert.throws(() => inspectAlignment(">a\nACGT\n>a\nACGT\n"), /duplicate identifier/);
  assert.throws(() => inspectAlignment(">a\nAC?T\n>b\nACGT\n"), /unsupported nucleotide character/);
});

test("FastTree preparation preserves the exact corrected gaps and row order", () => {
  const prepared = prepareFastTreeInput(ORIGINAL, "gtr", true);
  assert.equal(prepared.inputFasta, ">0\nAC-GT\n>1\nA-CGT\n>2\nACG-T\n");
  assert.equal(prepared.alignmentFasta, ORIGINAL);
  assert.equal(prepared.command, "fasttree -nt -gtr -fastest input.fasta");
  assert.deepEqual(prepared.names, ["read_a", "read_b", "__germline_N_masked__"]);
});
