/**
 * Public contracts for the phylogenetic UCA module.
 *
 * Nothing in this directory imports React or Swig's result-store runtime. The
 * engine can therefore be moved to a separate worker/package without changing
 * its statistical interface.
 */

export const PHYLO_UCA_CHARACTERS = ["A", "C", "G", "T", "-"] as const;
export type PhyloUcaCharacter = typeof PHYLO_UCA_CHARACTERS[number];
export type PhyloUcaFrameOffset = 0 | 1 | 2;

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
  /** Minimum retained D-template run which is treated as identifiable. */
  minimumDMatch: number;
  /** P(use an identifiable first D segment). */
  initialDProbability: number;
  /** P(add another D | a D has ended). */
  additionalDProbability: number;
  /** P(a junction contains at least one non-templated nucleotide). */
  junctionNProbability: number;
  /** Geometric tail ratio P(k + 1 trims) / P(k trims) at the V 3' end. */
  vThreePrimeTrimContinuation: number;
  /** Geometric continuation probability for trimming bases from the 5' D end. */
  dFivePrimeTrimContinuation: number;
  /** Geometric tail ratio P(k + 1 trims) / P(k trims) at the D 3' end. */
  dThreePrimeTrimContinuation: number;
  /** Geometric tail ratio P(k + 1 trims) / P(k trims) at the J 5' end. */
  jFivePrimeTrimContinuation: number;
  /** Expected N-run length conditional on a non-empty N run. */
  meanNLength: number;
  /** P(N length = 1 | the N run is non-empty). */
  singleNProbability: number;
  /** Number of geometric phases in the positive N-length tail. */
  nLengthPhases: number;
  /** Extra alignment columns on each side of the observed V/J anchors in which D states are evaluated. */
  junctionSearchFlankColumns: number;
  /** Robust leakage away from a deterministic germline nucleotide. */
  templateMismatchProbability: number;
  /** Prior probability that an N/junction alignment column is a gap. */
  junctionGapProbability: number;
  /** Gap prior in leading/trailing alignment padding outside a projected V/J template. */
  terminalPaddingGapProbability: number;
  /** A,C,G,T conditional probabilities inside N regions. */
  nBaseFrequencies: [number, number, number, number];
  /** Legacy pre-audit alias retained only when reading older sessions. */
  dExitProbability?: number;
  /** Legacy pre-audit V-boundary sigmoid width. */
  vTrimScale?: number;
  /** Legacy pre-audit J-boundary sigmoid width. */
  jTrimScale?: number;
  /** Legacy name for terminalPaddingGapProbability. */
  unknownTemplateGapProbability?: number;
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
  /** How uncertainty in the tree attachment and pendant length is handled. */
  inferenceMode: "maximum-likelihood" | "grid-marginalization" | "gibbs-mh";
  /** Cheap model used only to rank starting edges; every retained point is re-evaluated by the full HMM. */
  screenMode: "vj-mixture" | "germline-guide";
  /** Interior attachment positions evaluated by the cheap screen on every edge. */
  screenEdgeGridPoints: number;
  /** Number of screen-ranked edges receiving the full recombination HMM; zero means every edge. */
  fullHmmEdges: number;
  edgeGridPoints: number;
  /** Number of points in the explicit, zero-plus-logarithmic pendant-length grid. */
  branchGridPoints: number;
  /** Smallest positive point in the explicit pendant-length grid. */
  minimumPositiveUcaBranchLength: number;
  maximumUcaBranchLength: number;
  branchPriorMean: number;
  /** Coordinate-ascent rounds for conditional-ML distance/length optimization. */
  mlOptimizationRounds: number;
  /** Unit-interval stopping tolerance for conditional-ML scalar searches. */
  mlOptimizationTolerance: number;
  localRefinementRounds: number;
  marginalizeLocally: boolean;
  localPosteriorPoints: number;
  edgePrior: "uniform-edge" | "uniform-length";
  /** Metropolis-within-Gibbs iterations, including burn-in. */
  mcmcIterations: number;
  mcmcBurnIn: number;
  mcmcThin: number;
  /** Cheap placement/length MH updates performed for each exact HMM Gibbs draw. */
  mcmcMhStepsPerIteration: number;
  /** Reflected random-walk scale for pendant length, in substitutions/site. */
  mcmcBranchProposalScale: number;
  /** Reflected within-edge random-walk scale as a fraction of that edge. */
  mcmcPositionProposalScale: number;
  /** Probability that a placement MH update proposes a global edge jump. */
  mcmcGlobalJumpProbability: number;
  /** Fraction of global jumps centered on each edge's V/J-screen optimum rather than uniform position. */
  mcmcGlobalPositionMixture: number;
  /** Circular half-width of the screen-centered component of a global within-edge proposal. */
  mcmcGlobalPositionScale: number;
  /** Fraction of collapsed branch proposals drawn from a short-branch interval. */
  mcmcGlobalBranchMixture: number;
  /** Upper bound of the short-branch component, in substitutions/site. */
  mcmcGlobalBranchMaximum: number;
  /** Fraction of collapsed edge proposals informed by the already-computed full-HMM initializer scores. */
  mcmcCollapsedInitializerMixture: number;
  /** Iterations between exact collapsed placement refreshes; zero disables them. */
  mcmcCollapsedRefreshInterval: number;
  /** Reproducible 32-bit sampler seed. */
  mcmcSeed: number;
}

export interface PhyloUcaOptions {
  /** Auto selects GTR4 unless an observed tip contains an internal gap. */
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
  /** Full-width observed alignment used for tree messages and the HMM posterior. */
  observedAlignmentFasta: string;
  /** Original alignment columns represented in observedAlignmentFasta. */
  retainedColumns: number[];
  germlineGuideName: string;
  lineageRows: PhyloUcaAirrRow[];
  references: { V: string; D: string; J: string };
  locus: string;
  lineageLabel: string;
  alignmentFingerprint: string;
  /** Zero-based first nucleotide column of the selected codon frame. */
  frameOffset: PhyloUcaFrameOffset;
  options: PhyloUcaOptions;
}

export type PhyloUcaProgressPhase =
  | "references"
  | "tree-messages"
  | "edge-screen"
  | "hmm-search"
  | "mcmc"
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
  /** Log quadrature mass used only by explicit grid marginalization. */
  integrationLogWeight?: number;
  /** Cheap edge-screen score. This never substitutes for logMarginalLikelihood. */
  screenScore?: number;
  screenMode?: "vj-mixture" | "germline-guide";
  /** Legacy name retained so schema-1–3 sessions remain readable. */
  guideScore: number;
}

export interface PhyloUcaMcmcTracePoint {
  iteration: number;
  edgeId: string;
  edgeFraction: number;
  ucaBranchLength: number;
  /** Conditional log density for placement/length given the sampled UCA. */
  conditionalLogTarget: number;
  /** Full-HMM marginal likelihood at the pre-MH Gibbs state. */
  logMarginalLikelihood: number;
  retained: boolean;
}

export interface PhyloUcaMcmcDiagnostics {
  iterations: number;
  burnIn: number;
  thin: number;
  retainedSamples: number;
  mhStepsPerIteration: number;
  seed: number;
  branchProposals: number;
  branchAccepted: number;
  positionProposals: number;
  positionAccepted: number;
  globalProposals: number;
  globalAccepted: number;
  collapsedProposals: number;
  collapsedAccepted: number;
  edgeSwitches: number;
  branchEffectiveSampleSize: number;
  logTargetEffectiveSampleSize: number;
  /** Wall time spent in the sampling loop, excluding screening and initialization. */
  samplingMilliseconds?: number;
  /** Exact HMM conditional draws; each also evaluates the full marginal likelihood. */
  gibbsDraws?: number;
  gibbsMilliseconds?: number;
  /** Proposed collapsed full-HMM marginal evaluations, excluding accepted conditional redraws. */
  collapsedMarginalMilliseconds?: number;
  /** Cheap placement proposals evaluated conditional on one sampled UCA/path. */
  conditionalMhMilliseconds?: number;
  trace: PhyloUcaMcmcTracePoint[];
}

export interface PhyloUcaDCountPosteriorPoint {
  dCount: number;
  probability: number;
  /** Retained joint Gibbs draws having this D count. */
  samples: number;
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

export interface PhyloUcaHmmAnnotationPoint {
  /** One-based column in the original user-curated alignment. */
  alignmentColumn: number;
  /**
   * Unnormalized A/C/G/T/gap masses. Their sum is this track's posterior
   * occupancy at the column, rather than one.
   */
  probabilities: [number, number, number, number, number];
}

export interface PhyloUcaHmmAnnotationTrack {
  /** Stable across locally marginalized placement hypotheses. */
  id: string;
  kind: PhyloUcaSegmentKind;
  label: string;
  call?: string;
  dOrdinal?: number;
  /** D register = alignment-site index minus reference-D position. */
  registrationOffset?: number;
  /** True only when every occupied column has one fixed template character. */
  pure: boolean;
  points: PhyloUcaHmmAnnotationPoint[];
  maximumWeight: number;
}

export interface PhyloUcaHmmAnnotations {
  /** Allele groups must reach this occupancy in some column to be displayed. */
  minimumDisplayedWeight: number;
  /** Viterbi state path at the single best placement. */
  viterbi: PhyloUcaHmmAnnotationTrack[];
  /** Forward-backward state occupancy mixed across retained placements. */
  marginalized: PhyloUcaHmmAnnotationTrack[];
  omittedMarginalTrackCount: number;
  omittedMarginalMaximumWeight: number;
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

export interface PhyloUcaCodonPosterior {
  /** One-based codon ordinal in the selected alignment frame. */
  codonIndex: number;
  /** One-based original alignment columns forming this codon. */
  alignmentColumns: [number, number, number];
  /**
   * Joint probabilities in base-5 lexicographic A/C/G/T/gap order:
   * index = 25 * first + 5 * second + third.
   */
  probabilities: number[];
  mapCodon: string;
  mapProbability: number;
  entropyBits: number;
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
  schema: 1 | 2 | 3 | 4 | 5 | 6;
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
  /** Present in schema 2; schema-1 sessions implicitly used frame 0. */
  frameOffset?: PhyloUcaFrameOffset;
  bestPlacement: PhyloUcaPlacement;
  placements: PhyloUcaPlacement[];
  /** Exact full-HMM pendant-length grid in grid-marginalization mode. */
  evaluatedUcaBranchLengths?: number[];
  /** Present only for Metropolis-within-Gibbs inference. */
  mcmcDiagnostics?: PhyloUcaMcmcDiagnostics;
  /** Retained-draw posterior over the number of D segments; schema 6 Gibbs/MH results. */
  dCountPosterior?: PhyloUcaDCountPosteriorPoint[];
  mapAlignedSequence: string;
  mapUngappedSequence: string;
  posteriorConsensusAligned: string;
  posterior: PhyloUcaSitePosterior[];
  /** Exact three-column HMM/placement posterior; present in schema 2. */
  codonPosterior?: PhyloUcaCodonPosterior[];
  /** HMM-derived V/D/J/N annotation tracks; present in schema 3. */
  hmmAnnotations?: PhyloUcaHmmAnnotations;
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
  /** Codon posterior is valid only for this selected alignment frame. */
  frameOffset?: PhyloUcaFrameOffset;
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
