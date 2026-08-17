import assert from "node:assert/strict";
import test from "node:test";

import { prepareUploadedLineageTree } from "../src/uploaded-lineage-tree.ts";

test("uploaded lineage trees accept original sequence IDs and restore Swig alignment names", () => {
  const prepared = prepareUploadedLineageTree(
    "((alpha:0.1,beta:0.2):0.03,gamma:0.4);",
    ["alpha__1", "beta__2", "gamma__3"],
  );
  assert.match(prepared.canonicalNewick, /alpha__1/);
  assert.match(prepared.canonicalNewick, /beta__2/);
  assert.match(prepared.canonicalNewick, /gamma__3/);
});

test("uploaded lineage trees accept complete zero-based FastTree-style tips", () => {
  const prepared = prepareUploadedLineageTree(
    "((0:0.1,1:0.2):0.03,2:0.4);",
    ["alpha__1", "beta__2", "gamma__3"],
  );
  assert.match(prepared.rawNewick, /alpha__1/);
  assert.match(prepared.rawNewick, /gamma__3/);
});

test("ordinary display trees may retain the germline guide but are marked ineligible as observed-only UCA trees", () => {
  const prepared = prepareUploadedLineageTree(
    "((alpha:0.1,beta:0.2):0.03,(gamma:0.4,__germline_N_masked__:0.2):0.1);",
    ["alpha__1", "beta__2", "gamma__3"],
    true,
  );
  assert.equal(prepared.observedOnly, false);
  assert.match(prepared.rawNewick, /__germline_N_masked__/);
  assert.throws(
    () => prepareUploadedLineageTree(prepared.rawNewick, ["alpha__1", "beta__2", "gamma__3"]),
    /UCA inference accepts an observed-only tree/,
  );
});

test("uploaded lineage trees reject incomplete, duplicate, or length-free inputs", () => {
  assert.throws(
    () => prepareUploadedLineageTree("(alpha:0.1,beta:0.2);", ["alpha__1", "beta__2", "gamma__3"]),
    /has 2 tips.*has 3 biological rows/,
  );
  assert.throws(
    () => prepareUploadedLineageTree("((alpha:0.1,alpha:0.2):0.1,gamma:0.2);", ["alpha__1", "beta__2", "gamma__3"]),
    /more than once/,
  );
  assert.throws(
    () => prepareUploadedLineageTree("((alpha,beta),gamma);", ["alpha__1", "beta__2", "gamma__3"]),
    /needs branch lengths/,
  );
});
