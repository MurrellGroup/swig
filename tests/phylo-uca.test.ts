import assert from "node:assert/strict";
import test from "node:test";

import { defaultPhyloUcaOptions } from "../src/phylo-uca/defaults.ts";
import { compileGtr, HS5F_REVERSIBLE_GTR5 } from "../src/phylo-uca/gtr.ts";
import { phyloUcaHmmPosterior } from "../src/phylo-uca/hmm.ts";
import { inferPhyloUca } from "../src/phylo-uca/inference.ts";
import { prepareObservedOnlyAlignment, type PreparedPhyloUcaReferences } from "../src/phylo-uca/references.ts";
import { PhyloUcaTreeMessages, type ConditionalLikelihoodSurface } from "../src/phylo-uca/tree-messages.ts";

test("GTR4 and gap-aware GTR5 transition matrices are stochastic and reversible", () => {
  for (const includeGap of [false, true]) {
    const model = compileGtr(HS5F_REVERSIBLE_GTR5, includeGap);
    const identity = model.transition(0);
    const matrix = model.transition(0.173);
    for (let row = 0; row < model.dimension; row += 1) {
      let total = 0;
      for (let column = 0; column < model.dimension; column += 1) {
        assert.ok(Math.abs(identity[row * model.dimension + column] - Number(row === column)) < 1e-9);
        total += matrix[row * model.dimension + column];
        assert.ok(Math.abs(model.frequencies[row] * matrix[row * model.dimension + column] - model.frequencies[column] * matrix[column * model.dimension + row]) < 1e-8);
      }
      assert.ok(Math.abs(total - 1) < 1e-10);
    }
  }
});

test("automatic character selection uses GTR4 without gaps and GTR5 with an observed gap", () => {
  const tree = "((a:0.1,b:0.1):0.05,c:0.15);";
  const plain = new PhyloUcaTreeMessages(">a\nACGT\n>b\nACGT\n>c\nATGT\n", tree, HS5F_REVERSIBLE_GTR5, "auto");
  const gapped = new PhyloUcaTreeMessages(">a\nAC-T\n>b\nACGT\n>c\nATGT\n", tree, HS5F_REVERSIBLE_GTR5, "auto");
  assert.equal(plain.characterModel, "nucleotide-gtr4");
  assert.equal(plain.model.dimension, 4);
  assert.equal(gapped.characterModel, "gap-aware-gtr5");
  assert.equal(gapped.model.dimension, 5);
  const surface = gapped.conditionalLikelihoods(0, 0.04, 0.02);
  assert.equal(surface.stateCount, 5);
  assert.ok([...surface.logLikelihoods].every(Number.isFinite));
});

test("observed-only preparation removes the guide and only all-observed-gap columns", () => {
  const prepared = prepareObservedOnlyAlignment(">a\n-ACG-\n>b\n-ATG-\n>c\n-AC--\n>__germline_N_masked__\nAACGT\n", "__germline_N_masked__");
  assert.deepEqual(prepared.retainedColumns, [1, 2, 3]);
  assert.equal(prepared.rows, 3);
  assert.equal(prepared.guide, "ACG");
  assert.doesNotMatch(prepared.fasta, /germline/);
});

function exactSurface(sequence: string): ConditionalLikelihoodSurface {
  const values = new Float64Array(sequence.length * 4);
  for (let site = 0; site < sequence.length; site += 1) {
    const exact = "ACGT".indexOf(sequence[site]);
    for (let state = 0; state < 4; state += 1) values[site * 4 + state] = state === exact ? 0 : -12;
  }
  return { sites: sequence.length, stateCount: 4, logLikelihoods: values };
}

test("factorized HMM returns a complete V-to-J joint path", () => {
  const guide = "AAANNNNNNCCC";
  const references: PreparedPhyloUcaReferences = {
    v: [{ name: "IGHV1*01", sequence: "AAA", projection: "AAANNNNNNNNN", differences: 0, compared: 3, identity: 1, observedHypothesis: true }],
    d: [{ name: "IGHD1*01", sequence: "GG" }],
    j: [{ name: "IGHJ1*01", sequence: "CCC", projection: "NNNNNNNNNCCC", differences: 0, compared: 3, identity: 1, observedHypothesis: true }],
    vEndColumn: 2,
    jStartColumn: 9,
    guide,
    report: { locus: "IGH", v: ["IGHV1*01"], d: ["IGHD1*01"], j: ["IGHJ1*01"], totalVReferences: 1, totalDReferences: 1, totalJReferences: 1, observedVHypotheses: ["IGHV1*01"], observedJHypotheses: ["IGHJ1*01"], vCutoffDifferences: 0, jCutoffDifferences: 0, truncatedV: false, truncatedJ: false },
    warnings: [],
  };
  const options = defaultPhyloUcaOptions().hmm;
  const result = phyloUcaHmmPosterior(exactSurface("AAATGGTTTCCC"), references, { ...options, maximumDSegments: 1, minimumDMatch: 2 });
  assert.equal(result.mapVCall, "IGHV1*01");
  assert.equal(result.mapJCall, "IGHJ1*01");
  assert.equal(result.mapAlignedSequence.length, 12);
  assert.ok(Number.isFinite(result.logMarginalLikelihood));
  assert.equal(result.path[0].kind, "V");
  assert.equal(result.path.at(-1)?.kind, "J");
});

test("end-to-end inference roots at a zero-length UCA carrier and preserves the full pendant branch", async () => {
  const options = defaultPhyloUcaOptions();
  options.search = { ...options.search, fullHmmEdges: 3, edgeGridPoints: 2, branchGridPoints: 2, maximumUcaBranchLength: 0.05, localRefinementRounds: 0, localPosteriorPoints: 1 };
  options.hmm = { ...options.hmm, maximumDSegments: 1, minimumDMatch: 2 };
  const alignment = ">a__1\nAAATGGTTTCCC\n>b__2\nAAATGGCTTCCC\n>c__3\nAAATGGTTACCC\n>__germline_N_masked__\nAAANNNNNNCCC\n";
  const observed = prepareObservedOnlyAlignment(alignment, "__germline_N_masked__");
  const row = (ordinal: number, sequenceId: string) => ({ ordinal, sequenceId, locus: "IGH", values: { locus: "IGH", sequence_id: sequenceId, v_call: "IGHV1*01", j_call: "IGHJ1*01", v_sequence_alignment: "AAA", v_germline_alignment: "AAA", j_sequence_alignment: "CCC", j_germline_alignment: "CCC" } });
  const result = await inferPhyloUca({
    curatedAlignmentFasta: alignment,
    observedTreeNewick: "((a__1:0.01,b__2:0.01):0.01,c__3:0.02);",
    observedAlignmentFasta: observed.fasta,
    retainedColumns: observed.retainedColumns,
    germlineGuideName: "__germline_N_masked__",
    lineageRows: [row(0, "a"), row(1, "b"), row(2, "c")],
    references: { V: ">IGHV1*01\nAAA\n", D: ">IGHD1*01\nGG\n", J: ">IGHJ1*01\nCCC\n" },
    locus: "IGH",
    lineageLabel: "synthetic",
    alignmentFingerprint: "synthetic-fingerprint",
    options,
  });
  assert.equal(result.characterModel, "nucleotide-gtr4");
  assert.match(result.placedTreeNewick, /phylo_UCA:0(?:\.0+)?(?:[,)]|$)/);
  assert.equal(result.mapAlignedSequence.length, 12);
  assert.equal(result.posterior.length, 12);
  assert.equal(result.candidateReport.v[0], "IGHV1*01");
  assert.equal(result.candidateReport.j[0], "IGHJ1*01");
  assert.ok(result.bestPlacement.ucaBranchLength >= 0);
});
