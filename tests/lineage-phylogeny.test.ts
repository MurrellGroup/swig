import assert from "node:assert/strict";
import test from "node:test";

import {
  aminoAcidBranchMutations,
  inferAlignedReadingFrame,
  mapParsimonyMutations,
  motifCellMap,
  parseColumnSelection,
  spacedColumnOffsets,
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

test("Alivibe frame offset preserves a complete gap codon without accepting split codons", () => {
  const alignment = ["A---ATGGCC", "A---ATGGCT", "A---ATGGCA"];
  const inferred = inferAlignedReadingFrame(alignment);
  assert.equal(inferred.offset, 1);
  assert.equal(inferred.completeGapCodons, 3);
  assert.equal(inferred.mixedGapCodons, 0);
  assert.equal(translateAlignedNucleotides(alignment[0], inferred.offset), "-MA");
  assert.equal(translateAlignedNucleotides("A--AATGGCC", 0)[0], "X");
});

test("amino-acid branch mapping omits synonymous nucleotide changes", () => {
  const mutations = aminoAcidBranchMutations("GCTAAATTT", "GCCAGATTC", "tip-a");
  assert.deepEqual(mutations, [{ column: 1, from: "K", to: "R", childClade: "tip-a" }]);
});

test("custom selection accepts alignment ranges and Kabat-style labels", () => {
  assert.deepEqual(parseColumnSelection("1, 3-4, 31A", ["1", "2", "30", "31A", "31B"], 5), [0, 2, 3]);
});

test("non-contiguous alignment runs receive a half-cell visual separator", () => {
  assert.deepEqual(spacedColumnOffsets([0, 1, 5, 6, 10], 10), [0, 10, 25, 35, 50]);
});

test("motif mapping supports ungapped AA motifs and nucleotide IUPAC codes", () => {
  assert.deepEqual([...motifCellMap("CA-RDR", ["ARDR"], "aa")], [0, 1, 0, 1, 1, 1]);
  assert.deepEqual([...motifCellMap("AT-GC", ["RY"], "nt")], [1, 1, 0, 1, 1]);
});
