import assert from "node:assert/strict";
import test from "node:test";

import { chimeraVisiblePositions, classifyChimeraQuerySite } from "../src/chimera-view-model.ts";

test("chimera visualization removes only columns that are gaps in query and both displayed parents", () => {
  assert.deepEqual(chimeraVisiblePositions("-A-C", "--GC", "---C"), [1, 2, 3]);
  assert.deepEqual(chimeraVisiblePositions("A--T", "A--T", "C--T"), [0, 3]);
});

test("parent-match highlighter categories depend only on literal nucleotide identity", () => {
  assert.equal(classifyChimeraQuerySite("A", "A", "C"), "parent_a");
  assert.equal(classifyChimeraQuerySite("C", "A", "C"), "parent_b");
  assert.equal(classifyChimeraQuerySite("A", "A", "A"), "neutral");
  assert.equal(classifyChimeraQuerySite("G", "A", "C"), "neutral");
  assert.equal(classifyChimeraQuerySite("-", "A", "C"), "neutral");
  assert.equal(classifyChimeraQuerySite("N", "A", "C"), "neutral");
});
