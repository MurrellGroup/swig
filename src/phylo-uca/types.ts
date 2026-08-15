/**
 * Public contracts for the phylogenetic UCA module.
 *
 * Nothing in this directory imports React or Swig's result-store runtime. The
 * engine can therefore be moved to a separate worker/package without changing
 * its statistical interface.
 */

export const PHYLO_UCA_CHARACTERS = ["A", "C", "G", "T", "-"] as const;
export type PhyloUcaCharacter = typeof PHYLO_UCA_CHARACTERS[number];

export interface PhyloUcaGtrModel {
  /** Human-readable identifier stored in result provenance. */
  id: "hs5f-reversible" | "empirical-gtr" | "jc5" | "custom";
  label: string;
  /** A,C,G,T,gap stationary frequencies. Values are normalized by the engine. */
  frequencies: [number, number, number, number, number];
  /** AC,AG,AT,A-,CG,CT,C-,GT,G-,T- exchangeabilities. */
  exchangeabilities: [number, number, number, number, number, number, number, number, number, number];
  /** Free-form provenance suitable for JSON/session export. */
  provenance: string;
}

export interface PhyloUcaHmmOptions {
  /** Maximum number of D segments admitted by the recombination automaton. */
  maximumDSegments: number;
  minimumDMatch: number;
  /** P(add another D | a D has ended), before column-position normalization. */
  additionalDProbability: number;
  /** Geometric continuation probability for trimming bases from the 5' D end. */
  dFivePrimeTrimContinuation: number;
  /** Probability of exiting a D after the minimum emitted length. */
  dExitProbability: number;
  /** Expected number of N bases before a transition is attempted. */
  meanNLength: number;
  /** Width, in alignment columns, of the V trimming transition. */
  vTrimScale: number;
  /** Width, in alignment columns, of the J-entry transition. */
  jTrimScale: number;
  /** Robust leakage away from a deterministic germline nucleotide. */
  templateMismatchProbability: number;
  /** Prior probability that an N/junction alignment column is a gap. */
  junctionGapProbability: number;
  /** Prior probability of gap at a reference column whose projection is unknown. */
  unknownTemplateGapProbability: number;
  /** A,C,G,T conditional probabilities inside N regions. */
  nBaseFrequencies: [number, number, number, number];
}

export interface PhyloUcaCandidateOptions {
  vMaximumExtraDifferences: number;
  jMaximumExtraDifferences: number;
  vMinimumIdentity: number;
  jMinimumIdentity: number;
  maximumVCandidates: number;
  maximumJCandidates: number;
  /** Always retain observed selected/co-optimal/near-tied calls. */
  retainObservedHypotheses: boolean;
}

export interface PhyloUcaSearchOptions {
  /** Number of guide-ranked edges receiving the full recombination HMM. */
  fullHmmEdges: number;
  edgeGridPoints: number;
  branchGridPoints: number;
  maximumUcaBranchLength: number;
  branchPriorMean: number;
  localRefinementRounds: number;
  marginalizeLocally: boolean;
  localPosteriorPoints: number;
  edgePrior: "uniform-edge" | "uniform-length";
}

export interface PhyloUcaOptions {
  /** Auto selects GTR4 unless the observed curated alignment contains a gap. */
  characterMode: "auto" | "nucleotide-gtr4" | "gap-aware-gtr5";
  model: PhyloUcaGtrModel;
  hmm: PhyloUcaHmmOptions;
  candidates: PhyloUcaCandidateOptions;
  search: PhyloUcaSearchOptions;
}

export interface PhyloUcaReferenceRecord {
  name: string;
  sequence: string;
}

export interface PhyloUcaAirrRow {
  ordinal: number;
  sequenceId: string;
  locus: string;
  values: Record<string, string>;
}

export interface PhyloUcaInput {
  /** Exact user-curated nucleotide alignment, including the germline guide row. */
  curatedAlignmentFasta: string;
  /** FastTree inferred from observed rows only. */
  observedTreeNewick: string;
  /** Alignment actually supplied to observed-only FastTree after all-gap removal. */
  observedAlignmentFasta: string;
  /** Original alignment columns retained in observedAlignmentFasta. */
  retainedColumns: number[];
  germlineGuideName: string;
  lineageRows: PhyloUcaAirrRow[];
  references: { V: string; D: string; J: string };
  locus: string;
  lineageLabel: string;
  alignmentFingerprint: string;
  options: PhyloUcaOptions;
}

export type PhyloUcaProgressPhase =
  | "references"
  | "tree-messages"
  | "edge-screen"
  | "hmm-search"
  | "posterior"
  | "finalize";

export interface PhyloUcaProgress {
  phase: PhyloUcaProgressPhase;
  processed: number;
  total: number;
  detail: string;
}

export interface PhyloUcaPlacement {
  edgeId: string;
  endpointA: string;
  endpointB: string;
  distanceFromA: number;
  edgeLength: number;
  edgeFraction: number;
  ucaBranchLength: number;
  logMarginalLikelihood: number;
  logPosteriorScore: number;
  localPosteriorWeight: number;
  guideScore: number;
}

export type PhyloUcaSegmentKind = "V" | "N" | "D" | "J" | "unknown";

export interface PhyloUcaPathSegment {
  kind: PhyloUcaSegmentKind;
  call?: string;
  dOrdinal?: number;
  startColumn: number;
  endColumn: number;
  alignedSequence: string;
}

export interface PhyloUcaSitePosterior {
  /** One-based column in the original user-curated alignment. */
  alignmentColumn: number;
  probabilities: [number, number, number, number, number];
  mapCharacter: PhyloUcaCharacter;
  mapProbability: number;
  entropyBits: number;
  segment: PhyloUcaSegmentKind;
  call?: string;
}

export interface PhyloUcaCandidateReport {
  locus: string;
  v: string[];
  d: string[];
  j: string[];
  totalVReferences: number;
  totalDReferences: number;
  totalJReferences: number;
  observedVHypotheses: string[];
  observedJHypotheses: string[];
  vCutoffDifferences: number;
  jCutoffDifferences: number;
  truncatedV: boolean;
  truncatedJ: boolean;
}

export interface PhyloUcaResult {
  schema: 1;
  method: "fixed-tree-empirical-bayes-phylo-uca";
  lineageLabel: string;
  generatedAt: string;
  elapsedMs: number;
  characterModel: "nucleotide-gtr4" | "gap-aware-gtr5";
  model: PhyloUcaGtrModel;
  options: PhyloUcaOptions;
  observedTreeNewick: string;
  placedTreeNewick: string;
  observedAlignmentFasta: string;
  /** Zero-based original curated-alignment columns represented above. */
  retainedColumns: number[];
  alignmentFingerprint: string;
  bestPlacement: PhyloUcaPlacement;
  placements: PhyloUcaPlacement[];
  mapAlignedSequence: string;
  mapUngappedSequence: string;
  posteriorConsensusAligned: string;
  posterior: PhyloUcaSitePosterior[];
  path: PhyloUcaPathSegment[];
  candidateReport: PhyloUcaCandidateReport;
  mapVCall: string;
  mapDCalls: string[];
  mapJCall: string;
  logMarginalLikelihood: number;
  effectivePlacementCount: number;
  warnings: string[];
}

/** Portable/session-safe state for one exact lineage-alignment fingerprint. */
export interface PhyloUcaSavedState {
  lineageIds: number[];
  alignmentFingerprint: string;
  options: PhyloUcaOptions;
  result?: PhyloUcaResult;
}

export interface PhyloUcaWorkerRequest {
  id: number;
  input: PhyloUcaInput;
}

export type PhyloUcaWorkerResponse =
  | { id: number; type: "progress"; progress: PhyloUcaProgress }
  | { id: number; type: "result"; result: PhyloUcaResult }
  | { id: number; type: "error"; error: string };
