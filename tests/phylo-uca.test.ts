import assert from "node:assert/strict";
import test from "node:test";

import { defaultPhyloUcaOptions } from "../src/phylo-uca/defaults.ts";
import { compileGtr, HS5F_REVERSIBLE_GTR5 } from "../src/phylo-uca/gtr.ts";
import { phyloUcaHmmPosterior } from "../src/phylo-uca/hmm.ts";
import { inferPhyloUca } from "../src/phylo-uca/inference.ts";
import { prepareObservedOnlyAlignment, type PreparedPhyloUcaReferences } from "../src/phylo-uca/references.ts";
import { alignmentGapSemantics, observedAlignedCharacterPartial, PhyloUcaTreeMessages, type ConditionalLikelihoodSurface } from "../src/phylo-uca/tree-messages.ts";
import { aminoAcidUcaLogoColumns, codonUcaLogoColumns, nucleotideUcaLogoColumns } from "../src/phylo-uca/logo.ts";
import { PHYLO_UCA_CODON_SYMBOLS, phyloUcaCodonStateIndex } from "../src/phylo-uca/codons.ts";
import { normalizeProbabilityVector } from "../src/probability-logo.ts";
import type { PhyloUcaCodonPosterior, PhyloUcaSitePosterior } from "../src/phylo-uca/types.ts";

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

test("terminal tip gaps are missing data while internal gaps are exact fifth-state observations", () => {
  const semantics = alignmentGapSemantics("--AC-GT--");
  assert.deepEqual(semantics, { firstObserved: 2, lastObserved: 6, internalGaps: 1, terminalMissingGaps: 4 });
  assert.deepEqual([...observedAlignedCharacterPartial("--AC-GT--", 0, 5, semantics)], [1, 1, 1, 1, 1]);
  assert.deepEqual([...observedAlignedCharacterPartial("--AC-GT--", 4, 5, semantics)], [0, 0, 0, 0, 1]);
  assert.deepEqual([...observedAlignedCharacterPartial("--AC-GT--", 8, 5, semantics)], [1, 1, 1, 1, 1]);
  assert.deepEqual([...observedAlignedCharacterPartial("--AC-GT--", 4, 4, semantics)], [1, 1, 1, 1]);
});

test("automatic character selection ignores leading and trailing gap runs", () => {
  const tree = "((a:0.1,b:0.1):0.05,c:0.15);";
  const terminalOnly = new PhyloUcaTreeMessages(">a\n--ACGT\n>b\nTTAC--\n>c\nGGACGT\n", tree, HS5F_REVERSIBLE_GTR5, "auto");
  assert.equal(terminalOnly.internalGapCount, 0);
  assert.equal(terminalOnly.terminalMissingGapCount, 4);
  assert.equal(terminalOnly.characterModel, "nucleotide-gtr4");
});

function posteriorSite(column: number, probabilities: [number, number, number, number, number]): PhyloUcaSitePosterior {
  return { alignmentColumn: column, probabilities, mapCharacter: "A", mapProbability: 1, entropyBits: 0, segment: "V" };
}

test("frequency-logo vectors are normalized and exact codon projection preserves complete versus split gaps", () => {
  const normalized = normalizeProbabilityVector([2, 3, 5, Number.NaN, -1]);
  assert.ok(Math.abs(normalized.reduce((sum, value) => sum + value, 0) - 1) < 1e-15);
  const nt = nucleotideUcaLogoColumns([posteriorSite(1, [2, 3, 5, 0, 0])]);
  assert.ok(Math.abs(nt[0].entries.reduce((sum, entry) => sum + entry.probability, 0) - 1) < 1e-15);
  const exactCodon = (codon: string): PhyloUcaCodonPosterior => {
    const probabilities = Array.from({ length: 125 }, () => 0);
    probabilities[PHYLO_UCA_CODON_SYMBOLS.indexOf(codon)] = 1;
    return { codonIndex: 1, alignmentColumns: [1, 2, 3], probabilities, mapCodon: codon, mapProbability: 1, entropyBits: 0 };
  };
  const amino = (codon: string) => aminoAcidUcaLogoColumns([exactCodon(codon)])[0];
  assert.equal(codonUcaLogoColumns([exactCodon("ATG")])[0].entries.find((entry) => entry.symbol === "ATG")?.probability, 1);
  assert.equal(amino("ATG").entries.find((entry) => entry.symbol === "M")?.probability, 1);
  assert.equal(amino("---").entries.find((entry) => entry.symbol === "-")?.probability, 1);
  assert.equal(amino("A-G").entries.find((entry) => entry.symbol === "X")?.probability, 1);
});

test("observed-only preparation removes the guide and columns missing at every observed tip", () => {
  const prepared = prepareObservedOnlyAlignment(">a\n-ACG-\n>b\n-ATG-\n>c\n-AC--\n>__germline_N_masked__\nAACGT\n", "__germline_N_masked__");
  assert.deepEqual(prepared.retainedColumns, [1, 2, 3]);
  assert.deepEqual(prepared.posteriorColumns, [0, 1, 2, 3, 4]);
  assert.equal(prepared.rows, 3);
  assert.equal(prepared.guide, "ACG");
  assert.doesNotMatch(prepared.fasta, /germline/);
  assert.match(prepared.posteriorFasta, />a\n-ACG-/);
});

test("observed-only preparation retains a column that is an internal gap in every tip", () => {
  const prepared = prepareObservedOnlyAlignment(">a\nA-C\n>b\nT-G\n>c\nG-T\n>__germline_N_masked__\nANC\n", "__germline_N_masked__");
  assert.deepEqual(prepared.retainedColumns, [0, 1, 2]);
  assert.match(prepared.fasta, /A-C/);
  const tree = new PhyloUcaTreeMessages(prepared.fasta, "((a:0.1,b:0.1):0.05,c:0.15);", HS5F_REVERSIBLE_GTR5, "auto");
  assert.equal(tree.characterModel, "gap-aware-gtr5");
  assert.equal(tree.internalGapCount, 3);
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
  assert.equal(result.codonPosterior.length, 4);
  for (const codon of result.codonPosterior) {
    assert.ok(Math.abs(codon.probabilities.reduce((sum, probability) => sum + probability, 0) - 1) < 1e-12);
    for (let position = 0; position < 3; position += 1) {
      const marginal = [0, 0, 0, 0, 0];
      for (let state = 0; state < 125; state += 1) {
        const character = position === 0 ? Math.floor(state / 25) : position === 1 ? Math.floor(state / 5) % 5 : state % 5;
        marginal[character] += codon.probabilities[state];
      }
      for (let character = 0; character < 5; character += 1) assert.ok(Math.abs(marginal[character] - result.probabilities[codon.startSite + position][character]) < 1e-10);
    }
  }
  const shifted = phyloUcaHmmPosterior(exactSurface("AAATGGTTTCCC"), references, { ...options, maximumDSegments: 1, minimumDMatch: 2 }, 1);
  assert.deepEqual(shifted.codonPosterior.map((codon) => codon.startSite), [1, 4, 7]);
});

test("exact codon posterior retains germline-candidate correlation that site products lose", () => {
  const surface: ConditionalLikelihoodSurface = { sites: 6, stateCount: 4, logLikelihoods: new Float64Array(24) };
  const candidate = (name: string, projection: string) => ({ name, sequence: projection.slice(0, 3), projection, differences: 0, compared: 3, identity: 1, observedHypothesis: true });
  const references: PreparedPhyloUcaReferences = {
    v: [candidate("V_A", "AAANNN"), candidate("V_G", "GGGNNN")],
    d: [],
    j: [candidate("J", "NNNCCC")],
    vEndColumn: 2,
    jStartColumn: 3,
    guide: "NNNNNN",
    report: { locus: "IGH", v: ["V_A", "V_G"], d: [], j: ["J"], totalVReferences: 2, totalDReferences: 0, totalJReferences: 1, observedVHypotheses: [], observedJHypotheses: [], vCutoffDifferences: 0, jCutoffDifferences: 0, truncatedV: false, truncatedJ: false },
    warnings: [],
  };
  const options = { ...defaultPhyloUcaOptions().hmm, maximumDSegments: 0, vTrimScale: 0.25, jTrimScale: 0.25, templateMismatchProbability: 0.001 };
  const result = phyloUcaHmmPosterior(surface, references, options, 0);
  const exact = result.codonPosterior[0].probabilities;
  const aaa = phyloUcaCodonStateIndex(0, 0, 0);
  const independentAaa = result.probabilities[0][0] * result.probabilities[1][0] * result.probabilities[2][0];
  assert.ok(exact[aaa] > 0.45);
  assert.ok(exact[aaa] > 3.5 * independentAaa);
  assert.ok(exact[phyloUcaCodonStateIndex(0, 2, 2)] < 0.001);
});

test("end-to-end inference roots at a zero-length UCA carrier and preserves the full pendant branch", async () => {
  const options = defaultPhyloUcaOptions();
  options.search = { ...options.search, fullHmmEdges: 3, edgeGridPoints: 2, branchGridPoints: 2, maximumUcaBranchLength: 0.05, localRefinementRounds: 0, localPosteriorPoints: 3 };
  options.hmm = { ...options.hmm, maximumDSegments: 1, minimumDMatch: 2 };
  const alignment = ">a__1\nAAATGGTTTCCC\n>b__2\nAAATGGCTTCCC\n>c__3\nAAATGGTTACCC\n>__germline_N_masked__\nAAANNNNNNCCC\n";
  const observed = prepareObservedOnlyAlignment(alignment, "__germline_N_masked__");
  const row = (ordinal: number, sequenceId: string) => ({ ordinal, sequenceId, locus: "IGH", values: { locus: "IGH", sequence_id: sequenceId, v_call: "IGHV1*01", j_call: "IGHJ1*01", v_sequence_alignment: "AAA", v_germline_alignment: "AAA", j_sequence_alignment: "CCC", j_germline_alignment: "CCC" } });
  const result = await inferPhyloUca({
    curatedAlignmentFasta: alignment,
    observedTreeNewick: "((a__1:0.01,b__2:0.01):0.01,c__3:0.02);",
    observedAlignmentFasta: observed.posteriorFasta,
    retainedColumns: observed.posteriorColumns,
    germlineGuideName: "__germline_N_masked__",
    lineageRows: [row(0, "a"), row(1, "b"), row(2, "c")],
    references: { V: ">IGHV1*01\nAAA\n", D: ">IGHD1*01\nGG\n", J: ">IGHJ1*01\nCCC\n" },
    locus: "IGH",
    lineageLabel: "synthetic",
    alignmentFingerprint: "synthetic-fingerprint",
    frameOffset: 0,
    options,
  });
  assert.equal(result.characterModel, "nucleotide-gtr4");
  assert.match(result.placedTreeNewick, /phylo_UCA:0(?:\.0+)?(?:[,)]|$)/);
  assert.equal(result.mapAlignedSequence.length, 12);
  assert.equal(result.posterior.length, 12);
  assert.equal(result.schema, 2);
  assert.equal(result.codonPosterior?.length, 4);
  for (const codon of result.codonPosterior ?? []) for (let position = 0; position < 3; position += 1) {
    const marginal = [0, 0, 0, 0, 0];
    for (let state = 0; state < 125; state += 1) {
      const character = position === 0 ? Math.floor(state / 25) : position === 1 ? Math.floor(state / 5) % 5 : state % 5;
      marginal[character] += codon.probabilities[state];
    }
    const site = result.posterior[codon.alignmentColumns[position] - 1];
    for (let character = 0; character < 5; character += 1) assert.ok(Math.abs(marginal[character] - site.probabilities[character]) < 1e-9);
  }
  assert.equal(result.candidateReport.v[0], "IGHV1*01");
  assert.equal(result.candidateReport.j[0], "IGHJ1*01");
  assert.ok(result.bestPlacement.ucaBranchLength >= 0);
});
