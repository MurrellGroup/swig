import type { DatasetScope } from "../study-design.ts";

export type RefinementSegment = "V" | "D" | "J";
export type AlleleReassignmentPolicy = "best" | "confidence";
export type RefinementWeighting = "unique" | "abundance";

export interface AlleleRefinementOptions {
  /** Genotype/usage pools never cross this study boundary. Subject is the scientific default. */
  scope: DatasetScope;
  segments: RefinementSegment[];
  /** Symmetric Dirichlet pseudo-count assigned to every modelled reference element. */
  alphaPerAllele: number;
  /** Irreducible per-SNP neighbour evidence odds at zero estimated SHM. */
  baselineNeighbourOdds: number;
  /** Multiplier on the mechanistic mu/[3(1-mu)] SHM contribution. */
  shmLeakageSensitivity: number;
  /** Conservative cap on the per-SNP neighbour evidence odds. */
  maximumNeighbourOdds: number;
  /** Read-level SHM estimates are clamped here before adapting neighbour evidence. */
  maximumShm: number;
  /** Maximum reference-graph edit radius retained in each sparse row. */
  neighbourRadius: number;
  /** Converts SwiftIG score differences into relative evidence for retained near-tied hits. */
  alternativeScoreTemperature: number;
  /** Fallback evidence weight when a retained alternative lacks a usable score. */
  unscoredAlternativeWeight: number;
  /** Protects memory when an unusually dense reference neighbourhood is encountered. */
  maxCandidatesPerRow: number;
  /** Active representatives normally count once; abundance weighting is an explicit alternative. */
  weighting: RefinementWeighting;
  maxIterations: number;
  convergenceTolerance: number;
}

export const DEFAULT_ALLELE_REFINEMENT_OPTIONS: AlleleRefinementOptions = {
  scope: "subject",
  segments: ["V", "J"],
  alphaPerAllele: 0.1,
  baselineNeighbourOdds: 0.01,
  shmLeakageSensitivity: 1,
  maximumNeighbourOdds: 0.25,
  maximumShm: 0.3,
  neighbourRadius: 2,
  alternativeScoreTemperature: 2,
  unscoredAlternativeWeight: 0.25,
  maxCandidatesPerRow: 32,
  weighting: "unique",
  maxIterations: 100,
  convergenceTolerance: 1e-6,
};

export interface ReferenceAlleleNode {
  /** Stable position in this segment's reference graph. */
  index: number;
  segment: RefinementSegment;
  locus: string;
  /** All database identifiers with an identical nucleotide sequence. */
  names: string[];
  sequence: string;
}

export interface ReferenceNeighbour {
  index: number;
  distance: number;
  /** True when distance is a same-length nucleotide substitution count. */
  substitutionOnly: boolean;
}

export interface ReferenceAlleleGraph {
  segment: RefinementSegment;
  nodes: ReferenceAlleleNode[];
  callToNode: Map<string, number>;
  neighbours: ReferenceNeighbour[][];
  exactDuplicateLabels: number;
}

export interface RefinementInputRow {
  ordinal: number;
  sequenceId: string;
  datasetId: string;
  sampleId: string;
  subjectId: string;
  locus: string;
  call: string;
  score: number | null;
  identity: number | null;
  /** Best-reference read SHM estimate, normally 1 - best V-region identity. */
  shm: number | null;
  alternatives: string;
  abundance: number;
}

export interface SparseEvidenceMatrix {
  /** One more than the number of rows. */
  rowOffsets: Uint32Array;
  /** Reference-node index for every non-zero matrix entry. */
  columns: Uint32Array;
  /** Natural log of the local, pre-repertoire assignment evidence. */
  logEvidence: Float32Array;
  ordinals: Uint32Array;
  weights: Float32Array;
  groupKeys: string[];
  rowGroups: Uint32Array;
  localTop: Int32Array;
  localTopProbability: Float32Array;
  skippedRows: number;
  truncatedRows: number;
}

export interface AllelePosteriorSummary {
  nodeIndex: number;
  names: string[];
  sequenceLength: number;
  posteriorMean: number;
  expectedAssignments: number;
  localEvidenceAssignments: number;
  posteriorSd: number;
}

export interface RefinementModelSummary {
  key: string;
  scopeValue: string;
  locus: string;
  segment: RefinementSegment;
  rows: number;
  effectiveRows: number;
  nonZeros: number;
  /** Complete known-locus reference-node count under the Dirichlet prior. */
  databaseNodes: number;
  /** Prior-only nodes absent from every sparse candidate row in this pool. */
  inactivePriorNodes: number;
  alleles: AllelePosteriorSummary[];
  iterations: number;
  converged: boolean;
  finalMaximumChange: number;
}

export interface SegmentRefinementResult {
  segment: RefinementSegment;
  nodes: ReferenceAlleleNode[];
  /** MAP reference node for every AIRR ordinal; -1 means the row was not modelled. */
  mapNode: Int32Array;
  mapProbability: Float32Array;
  posteriorEntropy: Float32Array;
  localTopNode: Int32Array;
  localTopProbability: Float32Array;
  /** Fitted-model index for every AIRR ordinal; -1 means the row was not modelled. */
  modelIndex?: Int32Array;
  /** Unique-record or duplicate_count weight used by the fitted model. */
  assignmentWeight?: Float32Array;
  models: RefinementModelSummary[];
  modeledRows: number;
  changedMapRows: number;
  skippedRows: number;
  matrixNonZeros: number;
  truncatedRows: number;
  exactDuplicateLabels: number;
}

export interface AlleleRefinementResult {
  version: 1;
  options: AlleleRefinementOptions;
  totalRecords: number;
  activeRecords: number;
  segments: Partial<Record<RefinementSegment, SegmentRefinementResult>>;
  runAt: string;
  warnings: string[];
}

export interface SavedSegmentRefinement {
  segment: RefinementSegment;
  nodes: ReferenceAlleleNode[];
  mapNode: { type: "i32"; length: number; base64: string };
  mapProbability: { type: "f32"; length: number; base64: string };
  posteriorEntropy: { type: "f32"; length: number; base64: string };
  localTopNode: { type: "i32"; length: number; base64: string };
  localTopProbability: { type: "f32"; length: number; base64: string };
  modelIndex?: { type: "i32"; length: number; base64: string };
  assignmentWeight?: { type: "f32"; length: number; base64: string };
  models: RefinementModelSummary[];
  modeledRows: number;
  changedMapRows: number;
  skippedRows: number;
  matrixNonZeros: number;
  truncatedRows: number;
  exactDuplicateLabels: number;
}

export interface SavedAlleleRefinement {
  version: 1;
  options: AlleleRefinementOptions;
  totalRecords: number;
  activeRecords: number;
  segments: Partial<Record<RefinementSegment, SavedSegmentRefinement>>;
  runAt: string;
  warnings: string[];
  applied: boolean;
  reassignmentPolicy?: AlleleReassignmentPolicy;
  applyMinimumPosterior: number;
}
