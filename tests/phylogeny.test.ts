import assert from "node:assert/strict";
import test from "node:test";

import { extractNewick, parseNewick, rootOnOutgroup, serializeNewick } from "../src/phylogeny.ts";

function leafNames(node: ReturnType<typeof parseNewick>): string[] {
  return node.children.length ? node.children.flatMap(leafNames) : [node.name];
}

test("FastTree Newick extraction ignores diagnostics before and after the complete tree", () => {
  const output = [
    "FastTree Version 2.1.11 Double precision (No SSE3)",
    "Alignment: standard input",
    "Optimize all lengths: 1 rounds",
    "((0:0.010,1:0.020)0.91:0.030,GERMLINE_OUTGROUP:0.040);",
    "Total time: 0.12 seconds Unique: 3/3",
  ].join("\n");
  const newick = extractNewick(output);
  assert.equal(newick, "((0:0.010,1:0.020)0.91:0.030,GERMLINE_OUTGROUP:0.040);");
  assert.deepEqual(leafNames(parseNewick(output)).sort(), ["0", "1", "GERMLINE_OUTGROUP"]);
});

test("quoted Newick labels, comments, scientific lengths, and outgroup rooting round trip", () => {
  const parsed = parseNewick("(('read one':1.2e-3[&&NHX:S=1],read_two:0.2)95:0.3,GERMLINE_OUTGROUP:0.4)root;");
  assert.deepEqual(leafNames(parsed).sort(), ["GERMLINE_OUTGROUP", "read one", "read_two"]);
  const rooted = rootOnOutgroup(parsed, "GERMLINE_OUTGROUP");
  const roundTrip = parseNewick(`${serializeNewick(rooted)};`);
  assert.deepEqual(leafNames(roundTrip).sort(), ["GERMLINE_OUTGROUP", "read_one", "read_two"]);
});

test("incomplete FastTree output is rejected locally", () => {
  assert.throws(() => extractNewick("FastTree\n((a:0.1,b:0.2)\nTotal time"), /complete Newick tree/);
});
