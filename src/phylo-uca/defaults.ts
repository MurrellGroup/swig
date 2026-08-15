import { HS5F_REVERSIBLE_GTR5 } from "./gtr.ts";
import type { PhyloUcaOptions } from "./types.ts";

export function defaultPhyloUcaOptions(): PhyloUcaOptions {
  return {
    characterMode: "auto",
    model: { ...HS5F_REVERSIBLE_GTR5, frequencies: [...HS5F_REVERSIBLE_GTR5.frequencies], exchangeabilities: [...HS5F_REVERSIBLE_GTR5.exchangeabilities] },
    candidates: {
      vMaximumExtraDifferences: 6,
      jMaximumExtraDifferences: 4,
      vMinimumIdentity: 0.82,
      jMinimumIdentity: 0.75,
      maximumVCandidates: 48,
      maximumJCandidates: 24,
      retainObservedHypotheses: true,
    },
    hmm: {
      maximumDSegments: 3,
      minimumDMatch: 5,
      additionalDProbability: 0.015,
      dFivePrimeTrimContinuation: 0.72,
      dExitProbability: 0.28,
      meanNLength: 5,
      vTrimScale: 2.5,
      jTrimScale: 2.5,
      templateMismatchProbability: 0.003,
      junctionGapProbability: 0.015,
      unknownTemplateGapProbability: 0.01,
      nBaseFrequencies: [0.25, 0.25, 0.25, 0.25],
    },
    search: {
      screenMode: "vj-mixture",
      screenEdgeGridPoints: 5,
      fullHmmEdges: 6,
      edgeGridPoints: 3,
      branchGridPoints: 3,
      maximumUcaBranchLength: 0.3,
      branchPriorMean: 0.06,
      localRefinementRounds: 2,
      marginalizeLocally: true,
      localPosteriorPoints: 12,
      edgePrior: "uniform-length",
    },
  };
}
