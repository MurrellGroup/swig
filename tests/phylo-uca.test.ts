import assert from "node:assert/strict";
import test from "node:test";

import { defaultPhyloUcaOptions } from "../src/phylo-uca/defaults.ts";
import { compileGtr, HS5F_REVERSIBLE_GTR5 } from "../src/phylo-uca/gtr.ts";
import { PhyloUcaHmmGibbsSampler, phyloUcaHmmPosterior } from "../src/phylo-uca/hmm.ts";
import { collapseAndOrderHmmAnnotationTracks } from "../src/phylo-uca/hmm-annotation-model.ts";
import { inferPhyloUca, vjNucleotideMixtureProfile } from "../src/phylo-uca/inference.ts";
import { prepareObservedOnlyAlignment, type PreparedPhyloUcaReferences } from "../src/phylo-uca/references.ts";
import { alignmentGapSemantics, observedAlignedCharacterPartial, PhyloUcaTreeMessages, type ConditionalLikelihoodSurface } from "../src/phylo-uca/tree-messages.ts";
import { aminoAcidUcaLogoColumns, codonUcaLogoColumns, nucleotideUcaLogoColumns } from "../src/phylo-uca/logo.ts";
import { PHYLO_UCA_CODON_SYMBOLS, phyloUcaCodonStateIndex } from "../src/phylo-uca/codons.ts";
import { phyloUcaBranchLengthGrid } from "../src/phylo-uca/search-grid.ts";
import { phyloUcaPriorPredictiveSummary } from "../src/phylo-uca/prior-predictive.ts";
import { normalizeProbabilityVector } from "../src/probability-logo.ts";
import type { PhyloUcaCodonPosterior, PhyloUcaHmmAnnotationTrack, PhyloUcaSitePosterior, PhyloUcaSegmentKind } from "../src/phylo-uca/types.ts";

function annotationTrack(
  id: string,
  kind: PhyloUcaSegmentKind,
  label: string,
  alignmentColumn: number,
  probabilities: [number, number, number, number, number],
  details: Partial<Pick<PhyloUcaHmmAnnotationTrack, "call" | "dOrdinal" | "registrationOffset">> = {},
): PhyloUcaHmmAnnotationTrack {
  return { id, kind, label, ...details, pure: probabilities.filter((value) => value > 0).length <= 1, points: [{ alignmentColumn, probabilities }], maximumWeight: probabilities.reduce((sum, value) => sum + value, 0) };
}

test("phylogenetic UCA defaults use zero leakage and continuous Gibbs/MH", () => {
  const options = defaultPhyloUcaOptions();
  assert.equal(options.hmm.templateMismatchProbability, 0);
  assert.equal(options.search.inferenceMode, "gibbs-mh");
  assert.equal(options.search.mcmcCollapsedRefreshInterval, 3);
  assert.equal(options.search.mcmcCollapsedInitializerMixture, 0.95);
  assert.equal(options.hmm.initialDProbability, 0.934);
  assert.equal(options.hmm.junctionNProbability, 0.973);
});

test("HMM annotation display collapses duplicate D registers into one mixed allele row", () => {
  const raw = [
    annotationTrack("V|IGHV1*01", "V", "V · IGHV1*01", 1, [1, 0, 0, 0, 0], { call: "IGHV1*01" }),
    annotationTrack("N0", "N", "N0", 2, [0.1, 0.2, 0.3, 0.4, 0]),
    annotationTrack("D|1|IGHD1*01|+2", "D", "D1 · IGHD1*01 · register +2", 3, [0.55, 0, 0, 0, 0], { call: "IGHD1*01", dOrdinal: 1, registrationOffset: 2 }),
    annotationTrack("D|1|IGHD1*01|+3", "D", "D1 · IGHD1*01 · register +3", 3, [0, 0.45, 0, 0, 0], { call: "IGHD1*01", dOrdinal: 1, registrationOffset: 3 }),
    annotationTrack("N1", "N", "N1", 4, [0, 0, 0, 1, 0]),
    annotationTrack("J|IGHJ1*01", "J", "J · IGHJ1*01", 5, [0, 1, 0, 0, 0], { call: "IGHJ1*01" }),
  ];
  const display = collapseAndOrderHmmAnnotationTracks(raw);
  assert.deepEqual(display.map((track) => track.label), ["V · IGHV1*01", "NT1", "D · IGHD1*01", "NT2", "J · IGHJ1*01"]);
  const d = display.find((track) => track.kind === "D");
  assert.ok(d);
  assert.equal(d.sourceTrackCount, 2);
  assert.deepEqual(d.sourceRegistrationOffsets, [2, 3]);
  assert.equal(d.pure, false);
  assert.deepEqual(d.points[0].probabilities, [0.55, 0.45, 0, 0, 0]);
});

test("HMM annotation display follows weighted V, NT, D-block, NT, J order", () => {
  const raw = [
    annotationTrack("J", "J", "J · J1", 10, [0, 1, 0, 0, 0], { call: "J1" }),
    annotationTrack("D2", "D", "D2 · D2", 7, [0, 0, 1, 0, 0], { call: "D2", dOrdinal: 2, registrationOffset: 5 }),
    annotationTrack("N2", "N", "N2", 8, [0, 0, 0, 1, 0]),
    annotationTrack("D1b", "D", "D1 · D1b", 5, [0, 1, 0, 0, 0], { call: "D1b", dOrdinal: 1, registrationOffset: 4 }),
    annotationTrack("N0", "N", "N0", 2, [1, 0, 0, 0, 0]),
    annotationTrack("D1a", "D", "D1 · D1a", 4, [1, 0, 0, 0, 0], { call: "D1a", dOrdinal: 1, registrationOffset: 3 }),
    annotationTrack("V", "V", "V · V1", 1, [1, 0, 0, 0, 0], { call: "V1" }),
    annotationTrack("N1", "N", "N1", 6, [0, 0, 0, 1, 0]),
  ];
  const display = collapseAndOrderHmmAnnotationTracks(raw);
  assert.deepEqual(display.map((track) => track.label), ["V · V1", "NT1", "D · D1a", "D · D1b", "NT2", "D · D2", "NT3", "J · J1"]);
  assert.deepEqual(display.map((track) => track.weightedCenter), [1, 2, 4, 5, 6, 7, 8, 10]);
});

test("marginalized annotation puts one summed NT mixture above V, D, and J", () => {
  const raw = [
    annotationTrack("J", "J", "J · J1", 10, [0, 0.8, 0, 0, 0], { call: "J1" }),
    annotationTrack("N0", "N", "N0", 2, [0.2, 0, 0, 0, 0]),
    annotationTrack("D", "D", "D · D1", 6, [0, 0, 0.7, 0, 0], { call: "D1", dOrdinal: 1 }),
    annotationTrack("N1", "N", "N1", 7, [0, 0.3, 0, 0, 0]),
    annotationTrack("V", "V", "V · V1", 1, [0, 0, 0, 1, 0], { call: "V1" }),
  ];
  const display = collapseAndOrderHmmAnnotationTracks(raw, "marginalized");
  assert.deepEqual(display.map((track) => track.kind), ["N", "V", "D", "J"]);
  assert.equal(display[0].sourceTrackCount, 2);
  assert.deepEqual(display[0].points.map((point) => point.alignmentColumn), [2, 7]);
  assert.equal(display[0].integratedWeight, 0.5);
});

test("V/J screen profile mixes allele nucleotides independently and omits the junction", () => {
  const candidate = (name: string, projection: string) => ({ name, sequence: projection, projection, differences: 0, compared: 1, identity: 1, observedHypothesis: true });
  const profile = vjNucleotideMixtureProfile({
    v: [candidate("V1", "AANNN"), candidate("V2", "ACNNN")],
    j: [candidate("J1", "NNNGG"), candidate("J2", "NNNGT")],
    vEndColumn: 1,
    jStartColumn: 3,
    guide: "NNNNN",
  });
  assert.deepEqual(profile, [[2, 0, 0, 0, 0], [1, 1, 0, 0, 0], null, [0, 0, 2, 0, 0], [0, 0, 1, 1, 0]]);
});

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
  const sampler = new PhyloUcaHmmGibbsSampler(references, { ...options, maximumDSegments: 1, minimumDMatch: 2 });
  const surface = exactSurface("AAATGGTTTCCC");
  const gibbs = sampler.draw(surface, () => 0.37);
  assert.ok(Math.abs(sampler.logMarginal(surface) - result.logMarginalLikelihood) < 1e-9, "reusable backward marginal must equal the ordinary forward likelihood");
  assert.ok(Math.abs(gibbs.logMarginalLikelihood - result.logMarginalLikelihood) < 1e-9, "D-state FFBS backward likelihood must equal the ordinary forward likelihood");
  assert.equal(result.mapVCall, "IGHV1*01");
  assert.equal(result.mapJCall, "IGHJ1*01");
  assert.equal(result.mapAlignedSequence.length, 12);
  assert.ok(Number.isFinite(result.logMarginalLikelihood));
  assert.equal(result.path[0].kind, "V");
  assert.equal(result.path.at(-1)?.kind, "J");
  const dStart = result.path.find((segment) => segment.kind === "D")?.startColumn;
  assert.equal(dStart, 4);
  assert.notEqual((dStart ?? 0) % 3, 0, "the nucleotide HMM must allow a V/D boundary inside a displayed codon");
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
  const options = { ...defaultPhyloUcaOptions().hmm, maximumDSegments: 0, vThreePrimeTrimContinuation: 0, jFivePrimeTrimContinuation: 0, templateMismatchProbability: 0.001 };
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
  assert.equal(result.schema, 6);
  assert.equal(result.bestPlacement.screenMode, "vj-mixture");
  assert.ok(Number.isFinite(result.bestPlacement.screenScore));
  assert.ok(result.placements.some((placement) => placement.edgeFraction > 0 && placement.edgeFraction < 1));
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
  assert.ok(result.hmmAnnotations);
  assert.ok(result.hmmAnnotations.viterbi.some((track) => track.kind === "V" && track.call === "IGHV1*01" && track.pure));
  assert.ok(result.hmmAnnotations.viterbi.some((track) => track.kind === "D" && track.call === "IGHD1*01" && track.pure));
  assert.ok(result.hmmAnnotations.viterbi.some((track) => track.kind === "J" && track.call === "IGHJ1*01" && track.pure));
  assert.ok(result.hmmAnnotations.viterbi.some((track) => track.kind === "N" && !track.pure));
  for (const track of [...result.hmmAnnotations.viterbi, ...result.hmmAnnotations.marginalized]) if (track.pure) {
    for (const point of track.points) assert.ok(point.probabilities.filter((probability) => probability > 1e-12).length <= 1, `${track.label} mixed template characters at column ${point.alignmentColumn}`);
  }
  for (let column = 1; column <= result.posterior.length; column += 1) {
    const total = result.hmmAnnotations.viterbi.flatMap((track) => track.points).filter((point) => point.alignmentColumn === column).reduce((sum, point) => sum + point.probabilities.reduce((inner, value) => inner + value, 0), 0);
    assert.ok(Math.abs(total - 1) < 1e-10, `Viterbi source tracks sum to ${total} at column ${column}`);
  }
  assert.ok(result.bestPlacement.ucaBranchLength >= 0);
});

test("exact HMM FFBS draws reproduce the forward-backward nucleotide posterior", () => {
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
  const options = { ...defaultPhyloUcaOptions().hmm, maximumDSegments: 0, vThreePrimeTrimContinuation: 0, jFivePrimeTrimContinuation: 0, templateMismatchProbability: 0.02 };
  const surface: ConditionalLikelihoodSurface = { sites: 6, stateCount: 4, logLikelihoods: Float64Array.from([
    0, -0.3, -0.6, -0.8, 0, -0.2, -0.7, -0.9, 0, -0.4, -0.1, -0.8,
    -0.5, -0.6, -0.2, 0, -0.4, -0.2, 0, -0.5, -0.3, -0.4, -0.2, 0,
  ]) };
  const exact = phyloUcaHmmPosterior(surface, references, options);
  const sampler = new PhyloUcaHmmGibbsSampler(references, options);
  assert.ok(Math.abs(sampler.logMarginal(surface) - exact.logMarginalLikelihood) < 1e-9);
  let state = 9137;
  const random = () => {
    state = Math.imul(state ^ state >>> 15, 1 | state);
    state ^= state + Math.imul(state ^ state >>> 7, 61 | state);
    return ((state ^ state >>> 14) >>> 0) / 4294967296;
  };
  const counts = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0]);
  for (let drawIndex = 0; drawIndex < 5000; drawIndex += 1) {
    const draw = sampler.draw(surface, random);
    assert.ok(Math.abs(draw.logMarginalLikelihood - exact.logMarginalLikelihood) < 1e-9);
    for (let site = 0; site < 6; site += 1) counts[site][draw.characterStates[site]] += 1;
  }
  for (let site = 0; site < 6; site += 1) for (let character = 0; character < 4; character += 1) {
    assert.ok(Math.abs(counts[site][character] / 5000 - exact.probabilities[site][character]) < 0.035, `FFBS mismatch at site ${site}, character ${character}`);
  }
});

test("grid inference reports its exact zero-plus-logarithmic pendant grid", async () => {
  const options = defaultPhyloUcaOptions();
  options.search = { ...options.search, inferenceMode: "grid-marginalization", fullHmmEdges: 1, edgeGridPoints: 3, branchGridPoints: 5, minimumPositiveUcaBranchLength: 0.00001, maximumUcaBranchLength: 0.05, localPosteriorPoints: 15 };
  options.hmm = { ...options.hmm, maximumDSegments: 1, minimumDMatch: 2 };
  const alignment = ">a__1\nAAATGGTTTCCC\n>b__2\nAAATGGCTTCCC\n>c__3\nAAATGGTTACCC\n>__germline_N_masked__\nAAANNNNNNCCC\n";
  const observed = prepareObservedOnlyAlignment(alignment, "__germline_N_masked__");
  const row = (ordinal: number, sequenceId: string) => ({ ordinal, sequenceId, locus: "IGH", values: { locus: "IGH", sequence_id: sequenceId, v_call: "IGHV1*01", j_call: "IGHJ1*01", v_sequence_alignment: "AAA", v_germline_alignment: "AAA", j_sequence_alignment: "CCC", j_germline_alignment: "CCC" } });
  const result = await inferPhyloUca({ curatedAlignmentFasta: alignment, observedTreeNewick: "((a__1:0.01,b__2:0.01):0.01,c__3:0.02);", observedAlignmentFasta: observed.posteriorFasta, retainedColumns: observed.posteriorColumns, germlineGuideName: "__germline_N_masked__", lineageRows: [row(0, "a"), row(1, "b"), row(2, "c")], references: { V: ">IGHV1*01\nAAA\n", D: ">IGHD1*01\nGG\n", J: ">IGHJ1*01\nCCC\n" }, locus: "IGH", lineageLabel: "grid", alignmentFingerprint: "grid", frameOffset: 0, options });
  assert.deepEqual(result.evaluatedUcaBranchLengths, phyloUcaBranchLengthGrid(options.search));
  assert.equal(result.placements.length, 15);
  assert.ok(result.placements.every((point) => result.evaluatedUcaBranchLengths?.includes(point.ucaBranchLength)));
  assert.equal(result.mcmcDiagnostics, undefined);
});

test("Gibbs/MH keeps pendant length continuous and emits mixing diagnostics", async () => {
  const options = defaultPhyloUcaOptions();
  options.search = { ...options.search, inferenceMode: "gibbs-mh", fullHmmEdges: 1, maximumUcaBranchLength: 0.05, mcmcIterations: 24, mcmcBurnIn: 4, mcmcThin: 2, mcmcMhStepsPerIteration: 3, mcmcBranchProposalScale: 0.009, mcmcPositionProposalScale: 0.23, mcmcGlobalJumpProbability: 0.2, mcmcSeed: 4141 };
  options.hmm = { ...options.hmm, maximumDSegments: 1, minimumDMatch: 2 };
  const alignment = ">a__1\nAAATGGTTTCCC\n>b__2\nAAATGGCTTCCC\n>c__3\nAAATGGTTACCC\n>__germline_N_masked__\nAAANNNNNNCCC\n";
  const observed = prepareObservedOnlyAlignment(alignment, "__germline_N_masked__");
  const row = (ordinal: number, sequenceId: string) => ({ ordinal, sequenceId, locus: "IGH", values: { locus: "IGH", sequence_id: sequenceId, v_call: "IGHV1*01", j_call: "IGHJ1*01", v_sequence_alignment: "AAA", v_germline_alignment: "AAA", j_sequence_alignment: "CCC", j_germline_alignment: "CCC" } });
  const result = await inferPhyloUca({ curatedAlignmentFasta: alignment, observedTreeNewick: "((a__1:0.01,b__2:0.01):0.01,c__3:0.02);", observedAlignmentFasta: observed.posteriorFasta, retainedColumns: observed.posteriorColumns, germlineGuideName: "__germline_N_masked__", lineageRows: [row(0, "a"), row(1, "b"), row(2, "c")], references: { V: ">IGHV1*01\nAAA\n", D: ">IGHD1*01\nGG\n", J: ">IGHJ1*01\nCCC\n" }, locus: "IGH", lineageLabel: "mcmc", alignmentFingerprint: "mcmc", frameOffset: 0, options });
  const diagnostics = result.mcmcDiagnostics;
  assert.ok(diagnostics);
  assert.equal(diagnostics.retainedSamples, 10);
  assert.equal(diagnostics.trace.length, 24);
  assert.equal(result.evaluatedUcaBranchLengths, undefined);
  const unrelatedGrid = phyloUcaBranchLengthGrid({ ...options.search, branchGridPoints: 7 });
  assert.ok(diagnostics.trace.some((point) => unrelatedGrid.every((gridValue) => Math.abs(gridValue - point.ucaBranchLength) > 1e-10)), "continuous MH values must not lie on a branch grid");
  assert.ok(diagnostics.branchProposals > 0);
  assert.ok(diagnostics.positionProposals > 0);
  assert.equal(diagnostics.gibbsDraws, diagnostics.iterations + diagnostics.collapsedAccepted, "a collapsed proposal may draw a new latent HMM/UCA only after acceptance");
  assert.ok((diagnostics.samplingMilliseconds ?? -1) >= 0);
  assert.ok((diagnostics.gibbsMilliseconds ?? -1) >= 0);
  assert.ok((diagnostics.collapsedMarginalMilliseconds ?? -1) >= 0);
  assert.ok((diagnostics.conditionalMhMilliseconds ?? -1) >= 0);
  assert.ok(result.placements.every((point) => point.localPosteriorWeight === 0.1));
  assert.equal(result.dCountPosterior?.reduce((sum, point) => sum + point.samples, 0), 10);
  assert.ok(Math.abs((result.dCountPosterior?.reduce((sum, point) => sum + point.probability, 0) ?? 0) - 1) < 1e-12);
});

test("J states can enter only on concrete projected J nucleotides", () => {
  const candidate = (name: string, sequence: string, projection: string) => ({ name, sequence, projection, differences: 0, compared: sequence.length, identity: 1, observedHypothesis: true });
  const references: PreparedPhyloUcaReferences = {
    v: [candidate("V", "AAA", "AAANNNNNN")],
    d: [],
    j: [candidate("J", "CCC", "NNNNNNCCC")],
    vEndColumn: 2,
    jStartColumn: 6,
    guide: "AAANNNCCC",
    report: { locus: "IGH", v: ["V"], d: [], j: ["J"], totalVReferences: 1, totalDReferences: 0, totalJReferences: 1, observedVHypotheses: [], observedJHypotheses: [], vCutoffDifferences: 0, jCutoffDifferences: 0, truncatedV: false, truncatedJ: false },
    warnings: [],
  };
  const result = phyloUcaHmmPosterior(exactSurface("AAATTTCCC"), references, { ...defaultPhyloUcaOptions().hmm, maximumDSegments: 0, vThreePrimeTrimContinuation: 0, jFivePrimeTrimContinuation: 0 });
  const jPoints = result.marginalTracks.filter((track) => track.kind === "J").flatMap((track) => track.points);
  assert.ok(jPoints.length > 0);
  assert.ok(jPoints.every((point) => point.alignmentColumn >= 7), "J occupancy must not leak into the pre-J junction");
});

test("one UCA-branch mismatch plus six D matches beats an all-NT explanation at zero leakage", () => {
  const candidate = (name: string, sequence: string, projection: string) => ({ name, sequence, projection, differences: 0, compared: sequence.length, identity: 1, observedHypothesis: true });
  const references: PreparedPhyloUcaReferences = {
    v: [candidate("V", "AAA", "AAANNNNNNNNNNNNNCCC")],
    d: [{ name: "D", sequence: "CCCCCGAAAAAA" }],
    j: [candidate("J", "CCC", "NNNNNNNNNNNNNNNCCC")],
    vEndColumn: 2,
    jStartColumn: 15,
    guide: "AAANNNNNNNNNNNNNCCC",
    report: { locus: "IGH", v: ["V"], d: ["D"], j: ["J"], totalVReferences: 1, totalDReferences: 1, totalJReferences: 1, observedVHypotheses: [], observedJHypotheses: [], vCutoffDifferences: 0, jCutoffDifferences: 0, truncatedV: false, truncatedJ: false },
    warnings: [],
  };
  const sequence = "AAACCCCCTAAAAAACCC";
  const logLikelihoods = new Float64Array(sequence.length * 4);
  for (let site = 0; site < sequence.length; site += 1) {
    const demanded = "ACGT".indexOf(sequence[site]);
    for (let base = 0; base < 4; base += 1) logLikelihoods[site * 4 + base] = base === demanded ? 0 : -10;
  }
  // At the first D base the tree favors T over the D-template G by ~212:1,
  // representing a substitution on the UCA-to-tree branch rather than leakage.
  logLikelihoods[8 * 4 + 2] = -Math.log(212);
  const result = phyloUcaHmmPosterior({ sites: sequence.length, stateCount: 4, logLikelihoods }, references, {
    ...defaultPhyloUcaOptions().hmm,
    maximumDSegments: 1,
    minimumDMatch: 5,
    templateMismatchProbability: 0,
    vThreePrimeTrimContinuation: 0,
    jFivePrimeTrimContinuation: 0,
  });
  const terminalDWeight = result.marginalTracks.filter((track) => track.kind === "D").flatMap((track) => track.points).filter((point) => point.alignmentColumn === 15).reduce((sum, point) => sum + point.probabilities.reduce((inner, value) => inner + value, 0), 0);
  assert.ok(terminalDWeight > 0.75, `expected D to beat the one-base N boundary alternative after six matches, observed ${terminalDWeight}`);
});

test("prior-predictive generator reproduces displayed means and D-count probabilities", () => {
  const options = defaultPhyloUcaOptions().hmm;
  const summary = phyloUcaPriorPredictiveSummary(options, [17, 20, 23, 31, 37], 100_000, 9191);
  const metric = (id: string) => summary.metrics.find((entry) => entry.id === id)!;
  assert.ok(Math.abs(metric("v3").mean - 3.04) < 0.08);
  assert.ok(Math.abs(metric("nRun").mean - 8.8) < 0.15);
  assert.ok(Math.abs(summary.dCountProbabilities[0] - (1 - options.initialDProbability)) < 0.005);
  assert.ok(summary.dCountProbabilities[2] > 0);
});
