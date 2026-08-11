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

test("Kabat projection warns but still numbers an alignment with a terminal partial codon", () => {
  const result = inferKabatColumnsWithNumberer([{ name: "member__1", sequence: `${"GCT".repeat(60)}G` }], {
    number: (query) => {
      assert.equal(query.length, 60);
      return {
        chain: "H",
        confidence: 0.88,
        numbering: new Map(Array.from({ length: 60 }, (_, index) => [String(index + 1), "A"])),
        query_start: 0,
        error: null,
      };
    },
  });
  assert.equal(result.labels.length, 61);
  assert.equal(result.partialCodonRecords, 1);
  assert.equal(result.stopCodons, 0);
  assert.match(result.warnings[0], /terminal partial codon/);
  assert.equal(result.numberedColumns, 60);
});

test("Kabat projection replaces stop codons for numbering and reports them", () => {
  const sequence = `${"GCT".repeat(30)}TAA${"GCT".repeat(29)}`;
  const result = inferKabatColumnsWithNumberer([{ name: "member__1", sequence }], {
    number: (query) => {
      assert.equal(query[30], "X");
      return {
        chain: "H",
        confidence: 0.77,
        numbering: new Map(Array.from({ length: 60 }, (_, index) => [String(index + 1), query[index]])),
        query_start: 0,
        error: null,
      };
    },
  });
  assert.equal(result.stopCodons, 1);
  assert.match(result.warnings[0], /stop codon/);
});
