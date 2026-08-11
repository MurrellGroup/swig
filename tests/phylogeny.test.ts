import assert from "node:assert/strict";
import test from "node:test";

import {
  collapseShortInternalBranches,
  extractNewick,
  layoutTree,
  parseNewick,
  rootOnOutgroup,
  serializeNewick,
  type TreeNode,
} from "../src/phylogeny.ts";

function leafNames(node: ReturnType<typeof parseNewick>): string[] {
  return node.children.length ? node.children.flatMap(leafNames) : [node.name];
}

function pairwiseLeafDistances(root: TreeNode): Map<string, number> {
  const adjacency = new Map<TreeNode, Array<{ node: TreeNode; length: number }>>();
  const leaves: TreeNode[] = [];
  const visit = (node: TreeNode) => {
    if (!adjacency.has(node)) adjacency.set(node, []);
    if (!node.children.length) leaves.push(node);
    for (const child of node.children) {
      adjacency.get(node)!.push({ node: child, length: child.length });
      adjacency.set(child, [{ node, length: child.length }]);
      visit(child);
    }
  };
  visit(root);
  const result = new Map<string, number>();
  for (let left = 0; left < leaves.length; left += 1) {
    const distances = new Map<TreeNode, number>([[leaves[left], 0]]);
    const stack: Array<{ node: TreeNode; parent: TreeNode | null }> = [{ node: leaves[left], parent: null }];
    while (stack.length) {
      const { node, parent } = stack.pop()!;
      for (const edge of adjacency.get(node) ?? []) {
        if (edge.node === parent) continue;
        distances.set(edge.node, (distances.get(node) ?? 0) + edge.length);
        stack.push({ node: edge.node, parent: node });
      }
    }
    for (let right = left + 1; right < leaves.length; right += 1) {
      const names = [leaves[left].name, leaves[right].name].sort();
      result.set(names.join("|"), distances.get(leaves[right]) ?? Number.NaN);
    }
  }
  return result;
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

test("phylogram layout preserves FastTree zero-length branches", () => {
  const parsed = parseNewick("((a:0,b:0.5):0,c:1);");
  const phylogram = layoutTree(parsed, 900, 24, "phylogram");
  const a = phylogram.nodes.find((node) => node.name === "a")!;
  const b = phylogram.nodes.find((node) => node.name === "b")!;
  const c = phylogram.nodes.find((node) => node.name === "c")!;
  assert.equal(phylogram.mode, "phylogram");
  assert.equal(a.x, 24);
  assert.ok(b.x > a.x);
  assert.ok(c.x > b.x);

  const cladogram = layoutTree(parsed, 900, 24, "cladogram");
  assert.equal(cladogram.mode, "cladogram");
  assert.ok(cladogram.nodes.find((node) => node.name === "a")!.x > 24);
});

test("germline rerooting preserves every pairwise patristic distance", () => {
  const parsed = parseNewick("(((a:0.1,b:0.2):0.3,c:0.4):0.2,GERMLINE_OUTGROUP:0.5);");
  const before = pairwiseLeafDistances(parsed);
  const after = pairwiseLeafDistances(rootOnOutgroup(parsed, "GERMLINE_OUTGROUP"));
  assert.deepEqual([...after.keys()].sort(), [...before.keys()].sort());
  for (const [pair, distance] of before) assert.ok(Math.abs((after.get(pair) ?? Number.NaN) - distance) < 1e-12, pair);
});

test("platform-dependent FastTreeDbl floor resolutions collapse to one stable polytomy", () => {
  const first = parseNewick("((((a:0.01,b:0.01):5e-9,c:0.02):6e-9,d:0.03):0.1,GERMLINE_OUTGROUP:0.2);");
  const second = parseNewick("(((a:0.01,(b:0.01,c:0.02):5e-9):6e-9,d:0.03):0.1,GERMLINE_OUTGROUP:0.2);");
  const firstStable = collapseShortInternalBranches(first);
  const secondStable = collapseShortInternalBranches(second);

  assert.equal(firstStable.collapsedEdges, 2);
  assert.equal(secondStable.collapsedEdges, 2);
  assert.equal(serializeNewick(firstStable.root), serializeNewick(secondStable.root));
  assert.deepEqual(leafNames(firstStable.root).sort(), ["GERMLINE_OUTGROUP", "a", "b", "c", "d"]);
});

test("short terminal branches are retained and supported internal lengths are not collapsed", () => {
  const parsed = parseNewick("(((a:0,b:5e-9):0.001,c:0.02):0.1,GERMLINE_OUTGROUP:0.2);");
  const stable = collapseShortInternalBranches(parsed);
  assert.equal(stable.collapsedEdges, 0);
  const reparsed = parseNewick(`${serializeNewick(stable.root)};`);
  assert.deepEqual(leafNames(reparsed).sort(), ["GERMLINE_OUTGROUP", "a", "b", "c"]);
  assert.equal(reparsed.children[0].children[0].length, 0.001);
});
