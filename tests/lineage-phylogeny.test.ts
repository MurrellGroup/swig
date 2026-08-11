import assert from "node:assert/strict";
import test from "node:test";

import {
  mapParsimonyMutations,
  motifCellMap,
  parseColumnSelection,
  translateAlignedNucleotides,
} from "../src/lineage-phylogeny.ts";
import { parseNewick } from "../src/phylogeny.ts";

const OUTGROUP = "__germline_N_masked__";

test("germline N is inferred at the UCA without creating an N mutation", () => {
  const tree = parseNewick(`((a:0.1,b:0.1):0.1,${OUTGROUP}:0.1);`);
  const result = mapParsimonyMutations(tree, new Map([
    ["a", "AAC"],
    ["b", "AGC"],
    [OUTGROUP, "ANC"],
  ]), OUTGROUP);
  assert.equal(result.ucaSequence, "AAC");
  assert.equal(result.score, 1);
  assert.deepEqual([...result.mutationsByClade.values()].flat().map((mutation) => `${mutation.column}:${mutation.from}>${mutation.to}`), ["1:A>G"]);
});

test("known germline bases constrain the UCA and gaps remain parsimony states", () => {
  const tree = parseNewick(`((a:0.1,b:0.1):0.1,${OUTGROUP}:0.1);`);
  const result = mapParsimonyMutations(tree, new Map([
    ["a", "A-A"],
    ["b", "AGA"],
    [OUTGROUP, "G-A"],
  ]), OUTGROUP);
  assert.equal(result.ucaSequence, "G-A");
  const changes = [...result.mutationsByClade.values()].flat().map((mutation) => `${mutation.column}:${mutation.from}>${mutation.to}`).sort();
  assert.deepEqual(changes, ["0:G>A", "1:->G"]);
});

test("codon translation preserves codon columns", () => {
  assert.equal(translateAlignedNucleotides("ATG---GCNTA-"), "M-XX");
});

test("custom selection accepts alignment ranges and Kabat-style labels", () => {
  assert.deepEqual(parseColumnSelection("1, 3-4, 31A", ["1", "2", "30", "31A", "31B"], 5), [0, 2, 3]);
});

test("motif mapping supports ungapped AA motifs and nucleotide IUPAC codes", () => {
  assert.deepEqual([...motifCellMap("CA-RDR", ["ARDR"], "aa")], [0, 1, 0, 1, 1, 1]);
  assert.deepEqual([...motifCellMap("AT-GC", ["RY"], "nt")], [1, 1, 0, 1, 1]);
});
