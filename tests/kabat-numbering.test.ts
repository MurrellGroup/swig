import assert from "node:assert/strict";
import test from "node:test";

import { inferKabatColumnsWithNumberer } from "../src/kabat-numbering.ts";

test("Kabat labels are projected back through codon-alignment gaps", () => {
  const sequence = `${"GCT".repeat(2)}---${"GCT".repeat(58)}`;
  const result = inferKabatColumnsWithNumberer([{ name: "member__1", sequence }], {
    number: (query) => {
      assert.equal(query.length, 60);
      return {
        chain: "H",
        confidence: 0.91,
        numbering: new Map([["1", "A"], ["2", "A"], ["3", "A"]]),
        query_start: 1,
        error: null,
      };
    },
  });
  assert.equal(result.chain, "H");
  assert.equal(result.labels[1], "1");
  assert.equal(result.labels[2], "");
  assert.equal(result.labels[3], "2");
  assert.equal(result.labels[4], "3");
});

test("Kabat projection rejects nucleotide alignments that are not codon-column aligned", () => {
  assert.throws(() => inferKabatColumnsWithNumberer([{ name: "member__1", sequence: "GCTA" }], { number: () => { throw new Error("unreachable"); } }), /codon-column/);
});
