import type { ConditionalLikelihoodSurface } from "./tree-messages.ts";
import type {
  PhyloUcaCharacter,
  PhyloUcaFrameOffset,
  PhyloUcaHmmAnnotationTrack,
  PhyloUcaHmmOptions,
  PhyloUcaPathSegment,
  PhyloUcaSegmentKind,
} from "./types.ts";
import type { PreparedPhyloUcaReferences } from "./references.ts";
import { normalizeProbabilityVector } from "../probability-logo.ts";
import { PHYLO_UCA_CODON_STATE_COUNT, phyloUcaCodonStateIndex } from "./codons.ts";

const NEGATIVE_INFINITY = Number.NEGATIVE_INFINITY;
const CHARACTERS = ["A", "C", "G", "T", "-"] as const;
export const PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT = 0.01;
export const PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT = 0.001;

interface HmmState {
  kind: Exclude<PhyloUcaSegmentKind, "unknown">;
  call?: string;
  dOrdinal?: number;
  projection?: string;
  character?: string;
  /** Emission category for a fixed-character state such as D. */
  fixedCategory?: number;
  dUsed: number;
  dRun: number;
  dPosition: number;
  nMode?: "single" | "tail";
  nPhase?: number;
  templateIndex?: number;
  templateFirst: number;
  templateLast: number;
}

interface StateCatalog {
  states: HmmState[];
  v: number[];
  nSingle: number[];
  nTail: number[][];
  nonDStates: number[];
  dStates: number[];
  dByCount: number[][];
  dEntries: number[][];
  dEntryLogPrior: Float64Array;
  dContinue: Int32Array;
  dExitProbability: Float64Array;
  dStayLog: Float64Array;
  dExitLog: Float64Array;
  j: number[];
  /** Conditional V-exit hazard, indexed by site × V candidate. */
  vExitProbability: Float64Array;
  /** J allele/5'-trim log prior, indexed by site × J candidate. */
  jEntryLogPrior: Float64Array;
  sites: number;
  minimumDMatch: number;
}

export interface PhyloUcaHmmPosterior {
  logMarginalLikelihood: number;
  probabilities: Array<[number, number, number, number, number]>;
  /** Exact three-column probabilities conditional on the HMM and placement. */
  codonPosterior: Array<{ startSite: number; probabilities: number[] }>;
  mapAlignedSequence: string;
  posteriorConsensusAligned: string;
  stateKinds: PhyloUcaSegmentKind[];
  stateCalls: Array<string | undefined>;
  stateDOrdinals: Array<number | undefined>;
  path: PhyloUcaPathSegment[];
  /** Best-placement Viterbi state path, expressed as compact source tracks. */
  viterbiTracks: PhyloUcaHmmAnnotationTrack[];
  /** Best-placement forward-backward source occupancy. */
  marginalTracks: PhyloUcaHmmAnnotationTrack[];
  /** Track identities screened below visualization-only occupancy thresholds. */
  omittedMarginalTrackIds: string[];
  omittedMarginalMaximumWeight: number;
  mapVCall: string;
  mapDCalls: string[];
  mapJCall: string;
}

export interface PhyloUcaHmmGibbsDraw {
  /** Exact forward marginal at the placement used for this Gibbs update. */
  logMarginalLikelihood: number;
  /** One joint draw of UCA characters and the recombination-HMM path. */
  alignedSequence: string;
  characterStates: Int8Array;
  stateKinds: PhyloUcaSegmentKind[];
  stateCalls: Array<string | undefined>;
  stateDOrdinals: Array<number | undefined>;
  path: PhyloUcaPathSegment[];
  /** One-hot source occupancy for this draw, in the ordinary track schema. */
  tracks: PhyloUcaHmmAnnotationTrack[];
  mapVCall: string;
  mapDCalls: string[];
  mapJCall: string;
}

function boundedProbability(value: number, fallback: number): number {
  if (!Number.isFinite(value)) return fallback;
  return Math.max(0, Math.min(1, value));
}

function logAdd(left: number, right: number): number {
  if (left === NEGATIVE_INFINITY) return right;
  if (right === NEGATIVE_INFINITY) return left;
  const maximum = Math.max(left, right);
  return maximum + Math.log(Math.exp(left - maximum) + Math.exp(right - maximum));
}

function logSum(values: Iterable<number>): number {
  let total = NEGATIVE_INFINITY;
  for (const value of values) total = logAdd(total, value);
  return total;
}

function normalizedLogs(weights: readonly number[]): number[] {
  const safe = weights.map((weight) => Number.isFinite(weight) && weight > 0 ? weight : 0);
  const total = safe.reduce((sum, weight) => sum + weight, 0);
  if (!(total > 0)) return safe.map(() => NEGATIVE_INFINITY);
  return safe.map((weight) => weight > 0 ? Math.log(weight / total) : NEGATIVE_INFINITY);
}

function normalizedProjectionCharacter(projection: string | undefined, site: number): string {
  return (projection?.[site] ?? "N").toUpperCase().replace("U", "T").replace(".", "-");
}

function projectedNucleotideSites(projection: string): number[] {
  const sites: number[] = [];
  for (let site = 0; site < projection.length; site += 1) if (/[ACGT]/.test(normalizedProjectionCharacter(projection, site))) sites.push(site);
  return sites;
}

function fallbackReferencePositions(projection: string): number[] {
  const result = Array.from({ length: projection.length }, () => -1);
  let position = 0;
  for (let site = 0; site < projection.length; site += 1) {
    if (!/[ACGT]/.test(normalizedProjectionCharacter(projection, site))) continue;
    result[site] = position;
    position += 1;
  }
  return result;
}

function buildStateCatalog(references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): StateCatalog {
  const states: HmmState[] = [];
  const sites = references.guide.length || references.v[0]?.projection.length || references.j[0]?.projection.length || 0;
  const v = references.v.map((candidate, templateIndex) => {
    const index = states.length;
    const nucleotideSites = projectedNucleotideSites(candidate.projection);
    states.push({ kind: "V", call: candidate.name, projection: candidate.projection, dUsed: 0, dRun: 0, dPosition: -1, templateIndex, templateFirst: nucleotideSites[0] ?? 0, templateLast: nucleotideSites.at(-1) ?? -1 });
    return index;
  });
  const maximumD = Math.max(0, Math.min(5, Math.floor(options.maximumDSegments)));
  const nSingle = Array.from({ length: maximumD + 1 }, (_, dUsed) => {
    const index = states.length;
    states.push({ kind: "N", dUsed, dRun: 0, dPosition: -1, nMode: "single", nPhase: 0, templateFirst: 0, templateLast: sites - 1 });
    return index;
  });
  const nPhaseCount = Math.max(1, Math.min(4, Math.floor(options.nLengthPhases ?? 2)));
  const nTail = Array.from({ length: maximumD + 1 }, (_, dUsed) => Array.from({ length: nPhaseCount }, (_, nPhase) => {
    const index = states.length;
    states.push({ kind: "N", dUsed, dRun: 0, dPosition: -1, nMode: "tail", nPhase, templateFirst: 0, templateLast: sites - 1 });
    return index;
  }));
  const minimumDMatch = Math.max(1, Math.min(12, Math.floor(options.minimumDMatch)));
  const dByCount = Array.from({ length: maximumD + 1 }, () => [] as number[]);
  const dEntries = Array.from({ length: maximumD + 1 }, () => [] as number[]);
  const dKeys = new Map<string, number>();
  for (let dUsed = 1; dUsed <= maximumD; dUsed += 1) {
    for (const record of references.d) {
      if (record.sequence.length < minimumDMatch) continue;
      for (let position = 0; position < record.sequence.length; position += 1) {
        for (let run = 1; run <= minimumDMatch; run += 1) {
          const start = position - run + 1;
          if (start < 0 || record.sequence.length - start < minimumDMatch) continue;
          // Once the minimum is reached, paths with earlier starts deliberately
          // coalesce into one state: their probabilities have already been summed.
          if (run === minimumDMatch && position < minimumDMatch - 1) continue;
          const index = states.length;
          const dCharacter = record.sequence[position].toUpperCase().replace("U", "T");
          const dCategory = CHARACTERS.indexOf(dCharacter as PhyloUcaCharacter);
          states.push({
            kind: "D",
            call: record.name,
            character: record.sequence[position],
            fixedCategory: dCategory >= 0 ? dCategory : 5,
            dOrdinal: dUsed,
            dUsed,
            dRun: run,
            dPosition: position,
            templateFirst: 0,
            templateLast: sites - 1,
          });
          dKeys.set(`${dUsed}|${record.name}|${position}|${run}`, index);
          dByCount[dUsed].push(index);
          if (run === 1 && record.sequence.length - position >= minimumDMatch) dEntries[dUsed].push(index);
        }
      }
    }
  }
  const j = references.j.map((candidate, templateIndex) => {
    const index = states.length;
    const nucleotideSites = projectedNucleotideSites(candidate.projection);
    states.push({ kind: "J", call: candidate.name, projection: candidate.projection, dUsed: 0, dRun: 0, dPosition: -1, templateIndex, templateFirst: nucleotideSites[0] ?? sites, templateLast: nucleotideSites.at(-1) ?? sites - 1 });
    return index;
  });
  const dContinue = new Int32Array(states.length);
  dContinue.fill(-1);
  for (const stateIndex of dByCount.flat()) {
    const state = states[stateIndex];
    const nextRun = Math.min(minimumDMatch, state.dRun + 1);
    dContinue[stateIndex] = dKeys.get(`${state.dUsed}|${state.call}|${state.dPosition + 1}|${nextRun}`) ?? -1;
  }
  const dEntryLogPrior = new Float64Array(states.length);
  dEntryLogPrior.fill(NEGATIVE_INFINITY);
  const continuation = boundedProbability(options.dFivePrimeTrimContinuation, 0.8574);
  const eligibleDRecords = references.d.filter((record) => record.sequence.length >= minimumDMatch);
  for (let dUsed = 1; dUsed <= maximumD; dUsed += 1) {
    for (const record of eligibleDRecords) {
      const entries = dEntries[dUsed].filter((stateIndex) => states[stateIndex].call === record.name);
      const raw = entries.map((stateIndex) => Math.pow(continuation, states[stateIndex].dPosition));
      const total = raw.reduce((sum, value) => sum + value, 0);
      entries.forEach((stateIndex, entry) => {
        if (raw[entry] > 0 && total > 0) dEntryLogPrior[stateIndex] = -Math.log(eligibleDRecords.length) + Math.log(raw[entry] / total);
      });
    }
  }

  const dExitProbability = new Float64Array(states.length);
  const dThreePrimeContinuation = boundedProbability(options.dThreePrimeTrimContinuation ?? (options.dExitProbability === undefined ? 0.8471 : 1 - options.dExitProbability), 0.8471);
  const dLengths = new Map(references.d.map((record) => [record.name, record.sequence.length]));
  for (const stateIndex of dByCount.flat()) {
    const state = states[stateIndex];
    if (state.dRun < minimumDMatch) continue;
    const remainingAfterCurrent = Math.max(0, (dLengths.get(state.call ?? "") ?? state.dPosition + 1) - state.dPosition - 1);
    const currentWeight = Math.pow(dThreePrimeContinuation, remainingAfterCurrent);
    let remainingWeight = 0;
    for (let trim = 0; trim <= remainingAfterCurrent; trim += 1) remainingWeight += Math.pow(dThreePrimeContinuation, trim);
    dExitProbability[stateIndex] = remainingWeight > 0 ? currentWeight / remainingWeight : Number(state.dPosition >= (dLengths.get(state.call ?? "") ?? 1) - 1);
  }
  const dStayLog = new Float64Array(states.length);
  const dExitLog = new Float64Array(states.length);
  dStayLog.fill(NEGATIVE_INFINITY);
  dExitLog.fill(NEGATIVE_INFINITY);
  for (const stateIndex of dByCount.flat()) {
    const state = states[stateIndex];
    if (state.dRun < minimumDMatch) continue;
    const exitProbability = dExitProbability[stateIndex];
    if (dContinue[stateIndex] >= 0 && exitProbability < 1) dStayLog[stateIndex] = Math.log1p(-exitProbability);
    if (exitProbability > 0) dExitLog[stateIndex] = Math.log(exitProbability);
  }

  const vExitProbability = new Float64Array(Math.max(1, sites * Math.max(1, v.length)));
  const vContinuation = boundedProbability(options.vThreePrimeTrimContinuation, 0.7527);
  references.v.forEach((candidate, candidateIndex) => {
    const possibleSites = projectedNucleotideSites(candidate.projection);
    const referencePositions = candidate.referencePositions?.length === candidate.projection.length ? candidate.referencePositions : fallbackReferencePositions(candidate.projection);
    const weights = possibleSites.map((site) => Math.pow(vContinuation, Math.max(0, candidate.sequence.length - 1 - referencePositions[site])));
    let tail = 0;
    for (let entry = possibleSites.length - 1; entry >= 0; entry -= 1) {
      tail += weights[entry];
      if (tail > 0) vExitProbability[possibleSites[entry] * Math.max(1, v.length) + candidateIndex] = weights[entry] / tail;
    }
  });

  const jEntryLogPrior = new Float64Array(Math.max(1, sites * Math.max(1, j.length)));
  jEntryLogPrior.fill(NEGATIVE_INFINITY);
  const jContinuation = boundedProbability(options.jFivePrimeTrimContinuation, 0.8708);
  references.j.forEach((candidate, candidateIndex) => {
    const possibleSites = projectedNucleotideSites(candidate.projection);
    const referencePositions = candidate.referencePositions?.length === candidate.projection.length ? candidate.referencePositions : fallbackReferencePositions(candidate.projection);
    const weights = possibleSites.map((site) => Math.pow(jContinuation, Math.max(0, referencePositions[site])));
    const total = weights.reduce((sum, value) => sum + value, 0);
    possibleSites.forEach((site, entry) => {
      if (weights[entry] > 0 && total > 0) jEntryLogPrior[site * Math.max(1, j.length) + candidateIndex] = Math.log(weights[entry] / total);
    });
  });
  const dStates = dByCount.flat();
  const dStateSet = new Set(dStates);
  const nonDStates = states.map((_, index) => index).filter((index) => !dStateSet.has(index));
  return { states, v, nSingle, nTail, nonDStates, dStates, dByCount, dEntries, dEntryLogPrior, dContinue, dExitProbability, dStayLog, dExitLog, j, vExitProbability, jEntryLogPrior, sites, minimumDMatch };
}

function statePrior(state: HmmState, site: number, stateCount: 4 | 5, options: PhyloUcaHmmOptions): Float64Array {
  const result = new Float64Array(stateCount);
  const baseRaw = options.nBaseFrequencies.map((value) => Number.isFinite(value) && value > 0 ? value : 1e-9);
  const baseTotal = baseRaw.reduce((sum, value) => sum + value, 0);
  const bases = baseRaw.map((value) => value / baseTotal);
  const exact = state.kind === "N" ? "N" : state.kind === "D" ? state.character ?? "N" : state.projection?.[site] ?? "N";
  const normalized = exact.toUpperCase().replace("U", "T").replace(".", "-");
  const exactIndex = CHARACTERS.indexOf(normalized as PhyloUcaCharacter);
  const mismatch = Math.max(0, Math.min(0.25, Number.isFinite(options.templateMismatchProbability) ? options.templateMismatchProbability : 0));
  if (exactIndex >= 0 && exactIndex < 4) {
    const gapLeak = stateCount === 5 ? mismatch * 0.02 : 0;
    const nucleotideLeak = mismatch - gapLeak;
    for (let index = 0; index < 4; index += 1) result[index] = index === exactIndex ? 1 - mismatch : nucleotideLeak / 3;
    if (stateCount === 5) result[4] = gapLeak;
    return result;
  }
  if (normalized === "-" && stateCount === 5) {
    result[4] = 1 - mismatch;
    for (let index = 0; index < 4; index += 1) result[index] = mismatch * bases[index];
    return result;
  }
  // In GTR4 a projected reference gap is deliberately treated as unknown: the
  // observed alignment contains no gap character to identify an indel state.
  const gapProbability = stateCount === 5
    ? Math.max(0, Math.min(0.5, state.kind === "N" ? options.junctionGapProbability : options.terminalPaddingGapProbability ?? options.unknownTemplateGapProbability ?? 0.01))
    : 0;
  for (let index = 0; index < 4; index += 1) result[index] = (1 - gapProbability) * bases[index];
  if (stateCount === 5) result[4] = gapProbability;
  return result;
}

function stateAllowed(state: HmmState, site: number): boolean {
  if (state.kind === "N" || state.kind === "D") return true;
  const character = normalizedProjectionCharacter(state.projection, site);
  if (/[ACGT-]/.test(character)) return true;
  // Unknown projection is admitted only as terminal alignment padding. It is
  // never admitted before J entry or after V exit, so it cannot masquerade as
  // a uniform-emission junction state.
  return state.kind === "V" ? site < state.templateFirst : site > state.templateLast;
}

const TEMPLATE_UNKNOWN_CATEGORY = 5;
const JUNCTION_CATEGORY = 6;

function stateCategory(state: HmmState, site: number): number {
  if (state.kind === "N") return JUNCTION_CATEGORY;
  if (state.fixedCategory !== undefined) return state.fixedCategory;
  const character = state.kind === "D" ? (state.character ?? "N").toUpperCase().replace("U", "T").replace(".", "-") : normalizedProjectionCharacter(state.projection, site);
  const exact = CHARACTERS.indexOf(character as PhyloUcaCharacter);
  return exact >= 0 ? exact : TEMPLATE_UNKNOWN_CATEGORY;
}

function categoryPriors(stateCount: 4 | 5, options: PhyloUcaHmmOptions): Float64Array[] {
  const template = (character: string): HmmState => ({ kind: "D", character, dUsed: 1, dRun: 1, dPosition: 0, templateFirst: 0, templateLast: 0 });
  const junction: HmmState = { kind: "N", dUsed: 0, dRun: 0, dPosition: -1, templateFirst: 0, templateLast: 0 };
  return ["A", "C", "G", "T", "-", "N"].map((character) => statePrior(template(character), 0, stateCount, options)).concat([statePrior(junction, 0, stateCount, options)]);
}

interface EmissionCache {
  /** Site-major emission log likelihoods for A,C,G,T,gap,unknown-template,N. */
  categoryValues: Float64Array;
  categoryCount: number;
  siteOffsets: Float64Array;
  catalog: StateCatalog;
  window: JunctionWindow;
  prior: (state: number, site: number) => Float64Array;
}

interface JunctionWindow { start: number; end: number }

function junctionWindow(sites: number, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): JunctionWindow {
  // D/N uncertainty is local to the V–J junction. Keeping the all-D automaton
  // inactive elsewhere avoids scanning tens of thousands of D states across
  // framework columns. The exposed flank is a computational search bound, not
  // a hidden biological prior; increase it for unusually deep V/J trimming.
  const legacyFlank = Math.ceil(4 * Math.max(0.25, Math.max(options.vTrimScale ?? 2.5, options.jTrimScale ?? 2.5))) + 4;
  const flank = Math.max(0, Math.min(sites, Math.floor(options.junctionSearchFlankColumns ?? legacyFlank)));
  return {
    start: Math.max(0, references.vEndColumn - flank),
    end: Math.min(sites - 1, references.jStartColumn + flank),
  };
}

function buildEmissionCache(surface: ConditionalLikelihoodSurface, catalog: StateCatalog, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): EmissionCache {
  const stateCount = surface.stateCount;
  const siteOffsets = new Float64Array(surface.sites);
  const priors = categoryPriors(stateCount, options);
  const categoryValues = new Float64Array(surface.sites * priors.length);
  categoryValues.fill(NEGATIVE_INFINITY);
  const prior = (state: number, site: number) => priors[stateCategory(catalog.states[state], site)];
  const window = junctionWindow(surface.sites, references, options);
  for (let site = 0; site < surface.sites; site += 1) {
    const likelihoodOffset = site * stateCount;
    let siteMaximum = NEGATIVE_INFINITY;
    for (let character = 0; character < stateCount; character += 1) siteMaximum = Math.max(siteMaximum, surface.logLikelihoods[likelihoodOffset + character]);
    siteOffsets[site] = siteMaximum;
    for (let category = 0; category < priors.length; category += 1) {
      const probabilities = priors[category];
      let emission = NEGATIVE_INFINITY;
      for (let character = 0; character < stateCount; character += 1) {
        if (probabilities[character] > 0) emission = logAdd(emission, Math.log(probabilities[character]) + surface.logLikelihoods[likelihoodOffset + character] - siteMaximum);
      }
      categoryValues[site * priors.length + category] = emission;
    }
  }
  return { categoryValues, categoryCount: priors.length, siteOffsets, catalog, window, prior };
}

function emissionValue(emissions: EmissionCache, stateIndex: number, site: number): number {
  const state = emissions.catalog.states[stateIndex];
  if (state.kind === "D") {
    if (site < emissions.window.start || site > emissions.window.end) return NEGATIVE_INFINITY;
  } else if (!stateAllowed(state, site)) return NEGATIVE_INFINITY;
  return emissions.categoryValues[site * emissions.categoryCount + stateCategory(state, site)];
}

interface TransitionContext {
  catalog: StateCatalog;
  references: PreparedPhyloUcaReferences;
  options: PhyloUcaHmmOptions;
  window: JunctionWindow;
}

function routingLogs(context: TransitionContext, source: "V" | "N" | "D", dUsed: number): [number, number, number] {
  const { options, catalog } = context;
  const nProbability = boundedProbability(options.junctionNProbability, 0.973);
  const dProbability = dUsed === 0
    ? boundedProbability(options.initialDProbability, 0.934)
    : boundedProbability(options.additionalDProbability, 0.00125);
  const canUseAnotherD = dUsed < catalog.dEntries.length - 1 && (catalog.dEntries[dUsed + 1]?.length ?? 0) > 0;
  const nextDProbability = canUseAnotherD ? dProbability : 0;
  if (source === "N") return normalizedLogs([0, nextDProbability, 1 - nextDProbability]) as [number, number, number];
  return normalizedLogs([
    nProbability,
    (1 - nProbability) * nextDProbability,
    (1 - nProbability) * (1 - nextDProbability),
  ]) as [number, number, number];
}

function nDurationParameters(options: PhyloUcaHmmOptions): { single: number; tailStay: number } {
  const single = boundedProbability(options.singleNProbability, 0.027);
  const phases = Math.max(1, Math.min(4, Math.floor(options.nLengthPhases ?? 2)));
  const conditionalMean = Math.max(1, Number.isFinite(options.meanNLength) ? options.meanNLength : 8.8);
  const tailMean = single < 1 ? Math.max(phases, (conditionalMean - single) / (1 - single)) : phases;
  return { single, tailStay: Math.max(0, Math.min(1 - 1e-9, 1 - phases / tailMean)) };
}

function vExitProbability(catalog: StateCatalog, state: HmmState, site: number): number {
  if (state.templateIndex === undefined || !/[ACGT]/.test(normalizedProjectionCharacter(state.projection, site))) return 0;
  return catalog.vExitProbability[site * Math.max(1, catalog.v.length) + state.templateIndex] ?? 0;
}

function jEntryLogPrior(catalog: StateCatalog, state: HmmState, site: number): number {
  if (state.templateIndex === undefined || site < 0 || site >= catalog.sites) return NEGATIVE_INFINITY;
  return catalog.jEntryLogPrior[site * Math.max(1, catalog.j.length) + state.templateIndex] ?? NEGATIVE_INFINITY;
}

type CombineMode = "sum" | "max";

function combineDestination(destination: Float64Array, index: number, value: number, mode: CombineMode, sources?: Int32Array, source = -1): void {
  if (mode === "sum") destination[index] = logAdd(destination[index], value);
  else if (value > destination[index]) {
    destination[index] = value;
    if (sources) sources[index] = source;
  }
}

function aggregate(values: Float64Array, indexes: readonly number[], mode: CombineMode): { score: number; source: number } {
  let score = NEGATIVE_INFINITY;
  let source = -1;
  for (const index of indexes) {
    if (mode === "sum") score = logAdd(score, values[index]);
    else if (values[index] > score) {
      score = values[index];
      source = index;
    }
  }
  return { score, source };
}

function distributeD(destination: Float64Array, context: TransitionContext, dUsed: number, score: number, routeLog: number, mode: CombineMode, sources?: Int32Array, source = -1): void {
  if (score === NEGATIVE_INFINITY || routeLog === NEGATIVE_INFINITY) return;
  for (const target of context.catalog.dEntries[dUsed] ?? []) {
    combineDestination(destination, target, score + routeLog + context.catalog.dEntryLogPrior[target], mode, sources, source);
  }
}

function distributeN(destination: Float64Array, context: TransitionContext, dUsed: number, score: number, routeLog: number, mode: CombineMode, sources?: Int32Array, source = -1): void {
  if (score === NEGATIVE_INFINITY || routeLog === NEGATIVE_INFINITY) return;
  const duration = nDurationParameters(context.options);
  if (duration.single > 0) combineDestination(destination, context.catalog.nSingle[dUsed], score + routeLog + Math.log(duration.single), mode, sources, source);
  if (duration.single < 1) combineDestination(destination, context.catalog.nTail[dUsed][0], score + routeLog + Math.log1p(-duration.single), mode, sources, source);
}

function nEntryFuture(future: Float64Array, context: TransitionContext, dUsed: number): number {
  const duration = nDurationParameters(context.options);
  return logSum([
    duration.single > 0 ? Math.log(duration.single) + future[context.catalog.nSingle[dUsed]] : NEGATIVE_INFINITY,
    duration.single < 1 ? Math.log1p(-duration.single) + future[context.catalog.nTail[dUsed][0]] : NEGATIVE_INFINITY,
  ]);
}

function distributeJ(destination: Float64Array, context: TransitionContext, targetSite: number, score: number, routeLog: number, mode: CombineMode, sources?: Int32Array, source = -1): void {
  if (score === NEGATIVE_INFINITY || routeLog === NEGATIVE_INFINITY || !context.catalog.j.length) return;
  const candidateLog = -Math.log(context.catalog.j.length);
  for (const target of context.catalog.j) {
    const entryLog = jEntryLogPrior(context.catalog, context.catalog.states[target], targetSite);
    if (entryLog !== NEGATIVE_INFINITY) combineDestination(destination, target, score + routeLog + candidateLog + entryLog, mode, sources, source);
  }
}

function transitionForward(sourceValues: Float64Array, site: number, context: TransitionContext, mode: CombineMode, sources?: Int32Array): Float64Array {
  const { catalog, options, window } = context;
  const destination = new Float64Array(catalog.states.length);
  destination.fill(NEGATIVE_INFINITY);
  sources?.fill(-1);
  const vHubValues = new Float64Array(catalog.v.length);
  for (let entry = 0; entry < catalog.v.length; entry += 1) {
    const stateIndex = catalog.v[entry];
    const state = catalog.states[stateIndex];
    const exitProbability = vExitProbability(catalog, state, site);
    if (stateAllowed(state, site + 1) && exitProbability < 1) combineDestination(destination, stateIndex, sourceValues[stateIndex] + Math.log1p(-exitProbability), mode, sources, stateIndex);
    vHubValues[entry] = exitProbability > 0 ? sourceValues[stateIndex] + Math.log(exitProbability) : NEGATIVE_INFINITY;
  }
  const vHub = aggregate(vHubValues, catalog.v.map((_, index) => index), mode);
  const vSource = vHub.source >= 0 ? catalog.v[vHub.source] : -1;
  const vRoutes = routingLogs(context, "V", 0);
  distributeN(destination, context, 0, vHub.score, vRoutes[0], mode, sources, vSource);
  distributeD(destination, context, 1, vHub.score, vRoutes[1], mode, sources, vSource);
  distributeJ(destination, context, site + 1, vHub.score, vRoutes[2], mode, sources, vSource);

  const nDuration = nDurationParameters(options);
  for (let dUsed = 0; dUsed < catalog.nTail.length; dUsed += 1) {
    const routes = routingLogs(context, "N", dUsed);
    const singleState = catalog.nSingle[dUsed];
    distributeD(destination, context, dUsed + 1, sourceValues[singleState], routes[1], mode, sources, singleState);
    distributeJ(destination, context, site + 1, sourceValues[singleState], routes[2], mode, sources, singleState);
    const phases = catalog.nTail[dUsed];
    for (let phase = 0; phase < phases.length; phase += 1) {
      const state = phases[phase];
      if (nDuration.tailStay > 0) combineDestination(destination, state, sourceValues[state] + Math.log(nDuration.tailStay), mode, sources, state);
      if (phase + 1 < phases.length) combineDestination(destination, phases[phase + 1], sourceValues[state] + Math.log1p(-nDuration.tailStay), mode, sources, state);
      else {
        const hub = sourceValues[state] + Math.log1p(-nDuration.tailStay);
        distributeD(destination, context, dUsed + 1, hub, routes[1], mode, sources, state);
        distributeJ(destination, context, site + 1, hub, routes[2], mode, sources, state);
      }
    }
  }

  if (site >= window.start && site <= window.end) for (let dUsed = 1; dUsed < catalog.dByCount.length; dUsed += 1) {
    let exitScore = NEGATIVE_INFINITY;
    let exitSource = -1;
    for (const stateIndex of catalog.dByCount[dUsed]) {
      const state = catalog.states[stateIndex];
      const next = catalog.dContinue[stateIndex];
      if (state.dRun < catalog.minimumDMatch) {
        if (next >= 0) combineDestination(destination, next, sourceValues[stateIndex], mode, sources, stateIndex);
        continue;
      }
      const forcedExit = next < 0 || site >= window.end;
      if (next >= 0) combineDestination(destination, next, sourceValues[stateIndex] + (forcedExit ? NEGATIVE_INFINITY : catalog.dStayLog[stateIndex]), mode, sources, stateIndex);
      const candidate = sourceValues[stateIndex] + (forcedExit ? 0 : catalog.dExitLog[stateIndex]);
      if (mode === "sum") exitScore = logAdd(exitScore, candidate);
      else if (candidate > exitScore) {
        exitScore = candidate;
        exitSource = stateIndex;
      }
    }
    const routes = routingLogs(context, "D", dUsed);
    distributeN(destination, context, dUsed, exitScore, routes[0], mode, sources, exitSource);
    distributeD(destination, context, dUsed + 1, exitScore, routes[1], mode, sources, exitSource);
    distributeJ(destination, context, site + 1, exitScore, routes[2], mode, sources, exitSource);
  }
  for (const state of catalog.j) if (stateAllowed(catalog.states[state], site + 1)) combineDestination(destination, state, sourceValues[state], mode, sources, state);
  return destination;
}

function weightedDestination(values: Float64Array, indexes: readonly number[], logWeights: (index: number) => number): number {
  let result = NEGATIVE_INFINITY;
  for (const index of indexes) result = logAdd(result, values[index] + logWeights(index));
  return result;
}

function transitionBackward(future: Float64Array, site: number, context: TransitionContext, destination?: Float64Array): Float64Array {
  const { catalog, options, window } = context;
  const source = destination ?? new Float64Array(catalog.states.length);
  source.fill(NEGATIVE_INFINITY);
  const jFuture = weightedDestination(future, catalog.j, (index) => -Math.log(catalog.j.length) + jEntryLogPrior(catalog, catalog.states[index], site + 1));
  const dFuture = Array.from({ length: catalog.dEntries.length }, () => NEGATIVE_INFINITY);
  // D entry is impossible when the destination column is outside the exposed
  // junction window; skip those otherwise all-negative-infinity register scans.
  if (site + 1 >= window.start && site + 1 <= window.end) {
    for (let dUsed = 0; dUsed < catalog.dEntries.length; dUsed += 1) {
      dFuture[dUsed] = weightedDestination(future, catalog.dEntries[dUsed], (index) => catalog.dEntryLogPrior[index]);
    }
  }
  const vRoutes = routingLogs(context, "V", 0);
  const vHubFuture = logSum([
    vRoutes[0] + nEntryFuture(future, context, 0),
    vRoutes[1] + (dFuture[1] ?? NEGATIVE_INFINITY),
    vRoutes[2] + jFuture,
  ]);
  for (const stateIndex of catalog.v) {
    const state = catalog.states[stateIndex];
    const exitProbability = vExitProbability(catalog, state, site);
    const stay = stateAllowed(state, site + 1) && exitProbability < 1 ? Math.log1p(-exitProbability) + future[stateIndex] : NEGATIVE_INFINITY;
    const exit = exitProbability > 0 ? Math.log(exitProbability) + vHubFuture : NEGATIVE_INFINITY;
    source[stateIndex] = logAdd(stay, exit);
  }

  const nDuration = nDurationParameters(options);
  for (let dUsed = 0; dUsed < catalog.nTail.length; dUsed += 1) {
    const routes = routingLogs(context, "N", dUsed);
    const exitFuture = logSum([routes[1] + (dFuture[dUsed + 1] ?? NEGATIVE_INFINITY), routes[2] + jFuture]);
    source[catalog.nSingle[dUsed]] = exitFuture;
    const phases = catalog.nTail[dUsed];
    for (let phase = 0; phase < phases.length; phase += 1) {
      const state = phases[phase];
      const advance = phase + 1 < phases.length ? future[phases[phase + 1]] : exitFuture;
      source[state] = logAdd(Math.log(nDuration.tailStay) + future[state], Math.log1p(-nDuration.tailStay) + advance);
    }
  }

  if (site >= window.start && site <= window.end) for (let dUsed = 1; dUsed < catalog.dByCount.length; dUsed += 1) {
    const routes = routingLogs(context, "D", dUsed);
    const exitFuture = logSum([
      routes[0] + nEntryFuture(future, context, dUsed),
      routes[1] + (dFuture[dUsed + 1] ?? NEGATIVE_INFINITY),
      routes[2] + jFuture,
    ]);
    for (const stateIndex of catalog.dByCount[dUsed]) {
      const state = catalog.states[stateIndex];
      const next = catalog.dContinue[stateIndex];
      if (state.dRun < catalog.minimumDMatch) source[stateIndex] = next >= 0 ? future[next] : NEGATIVE_INFINITY;
      else {
        source[stateIndex] = next < 0 || site >= window.end
          ? exitFuture
          : logAdd(catalog.dStayLog[stateIndex] + future[next], catalog.dExitLog[stateIndex] + exitFuture);
      }
    }
  }
  for (const state of catalog.j) source[state] = stateAllowed(catalog.states[state], site + 1) ? future[state] : NEGATIVE_INFINITY;
  return source;
}

function initialize(catalog: StateCatalog, emissions: EmissionCache): Float64Array {
  const values = new Float64Array(catalog.states.length);
  values.fill(NEGATIVE_INFINITY);
  const prior = -Math.log(catalog.v.length);
  for (const state of catalog.v) values[state] = prior + emissionValue(emissions, state, 0);
  return values;
}

function terminalLogLikelihood(values: Float64Array, catalog: StateCatalog): number {
  return logSum(catalog.j.map((state) => values[state]));
}

function sampleLogCategorical(logWeights: ArrayLike<number>, random: () => number): number {
  let maximum = NEGATIVE_INFINITY;
  for (let index = 0; index < logWeights.length; index += 1) maximum = Math.max(maximum, logWeights[index]);
  if (!Number.isFinite(maximum)) throw new Error("The UCA Gibbs sampler encountered an empty conditional distribution.");
  let total = 0;
  for (let index = 0; index < logWeights.length; index += 1) if (Number.isFinite(logWeights[index])) total += Math.exp(logWeights[index] - maximum);
  let threshold = Math.max(0, Math.min(1 - Number.EPSILON, random())) * total;
  for (let index = 0; index < logWeights.length; index += 1) {
    if (!Number.isFinite(logWeights[index])) continue;
    threshold -= Math.exp(logWeights[index] - maximum);
    if (threshold <= 0) return index;
  }
  for (let index = logWeights.length - 1; index >= 0; index -= 1) if (Number.isFinite(logWeights[index])) return index;
  throw new Error("The UCA Gibbs sampler could not draw a finite state.");
}

/**
 * FFBS needs one backward row per site, but D states are possible only in the
 * bounded junction window. Persisting a dense sites × all-states matrix made
 * the ordinary framework columns dominate both allocation traffic and memory.
 * This store is exactly equivalent while retaining D values only where they
 * can be reached.
 */
class SparseBackwardRows {
  private readonly nonD: Float64Array;
  private readonly d: Float64Array;
  private readonly nonDPosition: Int32Array;
  private readonly dPosition: Int32Array;
  private readonly nonDStates: number[];
  private readonly dStates: number[];
  private readonly window: JunctionWindow;
  private readonly firstDState: number;

  constructor(catalog: StateCatalog, window: JunctionWindow, sites: number) {
    this.nonDStates = catalog.nonDStates;
    this.dStates = catalog.dStates;
    this.window = window;
    this.firstDState = this.dStates[0] ?? -1;
    this.nonD = new Float64Array(sites * this.nonDStates.length);
    this.d = new Float64Array(Math.max(0, window.end - window.start + 1) * this.dStates.length);
    this.nonDPosition = new Int32Array(catalog.states.length);
    this.dPosition = new Int32Array(catalog.states.length);
    this.nonDPosition.fill(-1);
    this.dPosition.fill(-1);
    this.nonDStates.forEach((state, position) => { this.nonDPosition[state] = position; });
    this.dStates.forEach((state, position) => { this.dPosition[state] = position; });
    for (let position = 0; position < this.dStates.length; position += 1) {
      if (this.dStates[position] !== this.firstDState + position) throw new Error("The UCA HMM D-state block is not contiguous.");
    }
  }

  set(site: number, values: Float64Array): void {
    const nonDOffset = site * this.nonDStates.length;
    for (let position = 0; position < this.nonDStates.length; position += 1) this.nonD[nonDOffset + position] = values[this.nonDStates[position]];
    if (site < this.window.start || site > this.window.end) return;
    const dOffset = (site - this.window.start) * this.dStates.length;
    if (this.dStates.length) this.d.set(values.subarray(this.firstDState, this.firstDState + this.dStates.length), dOffset);
  }

  get(site: number, state: number): number {
    const nonDPosition = this.nonDPosition[state];
    if (nonDPosition >= 0) return this.nonD[site * this.nonDStates.length + nonDPosition];
    const dPosition = this.dPosition[state];
    if (dPosition < 0 || site < this.window.start || site > this.window.end) return NEGATIVE_INFINITY;
    return this.d[(site - this.window.start) * this.dStates.length + dPosition];
  }
}

/**
 * Reusable exact FFBS sampler for p(HMM path, UCA characters | placement, tree).
 * The large D-state catalog is constructed once. Each draw performs one
 * backward HMM pass, samples one coherent path, and samples UCA characters
 * from the tree/template conditional at the visited states.
 */
export class PhyloUcaHmmGibbsSampler {
  private readonly catalog: StateCatalog;
  private readonly references: PreparedPhyloUcaReferences;
  private readonly options: PhyloUcaHmmOptions;
  private backwardRows?: SparseBackwardRows;

  constructor(
    references: PreparedPhyloUcaReferences,
    options: PhyloUcaHmmOptions,
  ) {
    this.references = references;
    this.options = options;
    this.catalog = buildStateCatalog(references, options);
    if (!this.catalog.v.length || !this.catalog.j.length) throw new Error("The UCA HMM Gibbs sampler requires at least one V and one J candidate.");
  }

  /** Exact full-HMM marginal using the reusable catalog and bounded D window. */
  logMarginal(surface: ConditionalLikelihoodSurface): number {
    const { catalog, references, options } = this;
    const emissions = buildEmissionCache(surface, catalog, references, options);
    const stateTotal = catalog.states.length;
    const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
    let beta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    beta.fill(NEGATIVE_INFINITY);
    for (const state of catalog.j) beta[state] = 0;
    const future: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    let nextBeta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    for (let site = surface.sites - 2; site >= 0; site -= 1) {
      future.fill(NEGATIVE_INFINITY);
      for (const state of catalog.nonDStates) future[state] = emissionValue(emissions, state, site + 1) + beta[state];
      if (site + 1 >= context.window.start && site + 1 <= context.window.end) {
        for (const state of catalog.dStates) future[state] = emissionValue(emissions, state, site + 1) + beta[state];
      }
      const reusable = beta;
      beta = transitionBackward(future, site, context, nextBeta);
      nextBeta = reusable;
    }
    const initial = initialize(catalog, emissions);
    for (let state = 0; state < stateTotal; state += 1) initial[state] += beta[state];
    return logSum(initial) + emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
  }

  draw(surface: ConditionalLikelihoodSurface, random: () => number = Math.random): PhyloUcaHmmGibbsDraw {
    const { catalog, references, options } = this;
    const emissions = buildEmissionCache(surface, catalog, references, options);
    const stateTotal = catalog.states.length;
    const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
    if (surface.sites !== catalog.sites) throw new Error("The UCA Gibbs surface length changed after sampler construction.");
    const betaRows = this.backwardRows ??= new SparseBackwardRows(catalog, context.window, surface.sites);
    let beta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    beta.fill(NEGATIVE_INFINITY);
    for (const state of catalog.j) beta[state] = 0;
    betaRows.set(surface.sites - 1, beta);
    const future: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    let nextBeta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
    for (let site = surface.sites - 2; site >= 0; site -= 1) {
      future.fill(NEGATIVE_INFINITY);
      for (const state of catalog.nonDStates) future[state] = emissionValue(emissions, state, site + 1) + beta[state];
      if (site + 1 >= context.window.start && site + 1 <= context.window.end) {
        for (const state of catalog.dStates) future[state] = emissionValue(emissions, state, site + 1) + beta[state];
      }
      const reusable = beta;
      beta = transitionBackward(future, site, context, nextBeta);
      nextBeta = reusable;
      betaRows.set(site, beta);
    }

    const initial = initialize(catalog, emissions);
    const initialConditional = new Float64Array(stateTotal);
    for (let state = 0; state < stateTotal; state += 1) initialConditional[state] = initial[state] + betaRows.get(0, state);
    const logMarginalLikelihood = logSum(initialConditional) + emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
    const statePath = new Int32Array(surface.sites);
    statePath[0] = sampleLogCategorical(initialConditional, random);

    const addDestination = (destinations: number[], weights: number[], target: number, transitionLog: number, site: number) => {
      if (target < 0 || !Number.isFinite(transitionLog)) return;
      const future = emissionValue(emissions, target, site + 1) + betaRows.get(site + 1, target);
      if (!Number.isFinite(future)) return;
      destinations.push(target);
      weights.push(transitionLog + future);
    };
    const addNEntry = (destinations: number[], weights: number[], dUsed: number, transitionLog: number, site: number) => {
      const duration = nDurationParameters(options);
      if (duration.single > 0) addDestination(destinations, weights, catalog.nSingle[dUsed], transitionLog + Math.log(duration.single), site);
      if (duration.single < 1) addDestination(destinations, weights, catalog.nTail[dUsed][0], transitionLog + Math.log1p(-duration.single), site);
    };

    for (let site = 0; site < surface.sites - 1; site += 1) {
      const sourceIndex = statePath[site];
      const source = catalog.states[sourceIndex];
      const destinations: number[] = [];
      const weights: number[] = [];
      if (source.kind === "V") {
        const exitProbability = vExitProbability(catalog, source, site);
        if (stateAllowed(source, site + 1) && exitProbability < 1) addDestination(destinations, weights, sourceIndex, Math.log1p(-exitProbability), site);
        const routes = routingLogs(context, "V", 0);
        const hub = exitProbability > 0 ? Math.log(exitProbability) : NEGATIVE_INFINITY;
        addNEntry(destinations, weights, 0, hub + routes[0], site);
        for (const target of catalog.dEntries[1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
        const candidateLog = -Math.log(catalog.j.length);
        for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog + jEntryLogPrior(catalog, catalog.states[target], site + 1), site);
      } else if (source.kind === "N") {
        const routes = routingLogs(context, "N", source.dUsed);
        let hub = 0;
        if (source.nMode === "tail") {
          const duration = nDurationParameters(options);
          addDestination(destinations, weights, sourceIndex, Math.log(duration.tailStay), site);
          const phases = catalog.nTail[source.dUsed];
          const phase = source.nPhase ?? 0;
          if (phase + 1 < phases.length) addDestination(destinations, weights, phases[phase + 1], Math.log1p(-duration.tailStay), site);
          else hub = Math.log1p(-duration.tailStay);
          if (phase + 1 < phases.length) hub = NEGATIVE_INFINITY;
        }
        if (hub !== NEGATIVE_INFINITY) {
          for (const target of catalog.dEntries[source.dUsed + 1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
          const candidateLog = -Math.log(catalog.j.length);
          for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog + jEntryLogPrior(catalog, catalog.states[target], site + 1), site);
        }
      } else if (source.kind === "D") {
        const next = catalog.dContinue[sourceIndex];
        if (source.dRun < catalog.minimumDMatch) addDestination(destinations, weights, next, 0, site);
        else {
          const forcedExit = next < 0 || site >= context.window.end;
          if (next >= 0) addDestination(destinations, weights, next, forcedExit ? NEGATIVE_INFINITY : catalog.dStayLog[sourceIndex], site);
          const routes = routingLogs(context, "D", source.dUsed);
          const hub = forcedExit ? 0 : catalog.dExitLog[sourceIndex];
          addNEntry(destinations, weights, source.dUsed, hub + routes[0], site);
          for (const target of catalog.dEntries[source.dUsed + 1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
          const candidateLog = -Math.log(catalog.j.length);
          for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog + jEntryLogPrior(catalog, catalog.states[target], site + 1), site);
        }
      } else if (stateAllowed(source, site + 1)) addDestination(destinations, weights, sourceIndex, 0, site);
      const selected = sampleLogCategorical(weights, random);
      statePath[site + 1] = destinations[selected];
    }

    const characterStates = new Int8Array(surface.sites);
    let alignedSequence = "";
    const stateKinds: PhyloUcaSegmentKind[] = [];
    const stateCalls: Array<string | undefined> = [];
    const stateDOrdinals: Array<number | undefined> = [];
    const trackAccumulator = new Map<string, TrackAccumulator>();
    const siteTotals = new Float64Array(surface.sites);
    siteTotals.fill(1);
    for (let site = 0; site < surface.sites; site += 1) {
      const stateIndex = statePath[site];
      const conditional = conditionalCharacterProbabilities(surface, emissions, stateIndex, site);
      const character = sampleLogCategorical(Float64Array.from(conditional, (probability) => probability > 0 ? Math.log(probability) : NEGATIVE_INFINITY), random);
      characterStates[site] = character;
      alignedSequence += CHARACTERS[character];
      const oneHot = new Float64Array(5);
      oneHot[character] = 1;
      addTrackMass(trackAccumulator, catalog.states[stateIndex], site, 1, oneHot);
      const state = catalog.states[stateIndex];
      stateKinds.push(state.kind);
      stateCalls.push(state.call);
      stateDOrdinals.push(state.dOrdinal);
    }
    const path = pathSegments(catalog.states, statePath, alignedSequence);
    return {
      logMarginalLikelihood,
      alignedSequence,
      characterStates,
      stateKinds,
      stateCalls,
      stateDOrdinals,
      path,
      tracks: finalizeTracks(trackAccumulator, siteTotals),
      mapVCall: path.find((segment) => segment.kind === "V")?.call ?? "",
      mapDCalls: path.filter((segment) => segment.kind === "D").map((segment) => segment.call ?? ""),
      mapJCall: path.find((segment) => segment.kind === "J")?.call ?? "",
    };
  }
}

function forward(
  surface: ConditionalLikelihoodSurface,
  references: PreparedPhyloUcaReferences,
  options: PhyloUcaHmmOptions,
  keepRows: boolean,
  mode: CombineMode = "sum",
): { logLikelihood: number; catalog: StateCatalog; emissions: EmissionCache; rows?: Float64Array; backpointers?: Int32Array; endState?: number } {
  const catalog = buildStateCatalog(references, options);
  if (!catalog.v.length || !catalog.j.length) throw new Error("The UCA HMM requires at least one V and one J candidate.");
  const emissions = buildEmissionCache(surface, catalog, references, options);
  const stateTotal = catalog.states.length;
  let current = initialize(catalog, emissions);
  const rows = keepRows ? new Float64Array(surface.sites * stateTotal) : undefined;
  rows?.set(current, 0);
  const backpointers = mode === "max" ? new Int32Array(surface.sites * stateTotal) : undefined;
  backpointers?.fill(-1);
  const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
  for (let site = 0; site < surface.sites - 1; site += 1) {
    const sources = mode === "max" ? new Int32Array(stateTotal) : undefined;
    const next = transitionForward(current, site, context, mode, sources);
    for (let state = 0; state < stateTotal; state += 1) if (next[state] !== NEGATIVE_INFINITY) next[state] += emissionValue(emissions, state, site + 1);
    if (backpointers && sources) backpointers.set(sources, (site + 1) * stateTotal);
    current = next;
    rows?.set(current, (site + 1) * stateTotal);
  }
  let endState = -1;
  let terminal = NEGATIVE_INFINITY;
  if (mode === "sum") terminal = terminalLogLikelihood(current, catalog);
  else for (const state of catalog.j) if (current[state] > terminal) {
    terminal = current[state];
    endState = state;
  }
  const siteOffset = emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
  return { logLikelihood: terminal + siteOffset, catalog, emissions, rows, backpointers, endState };
}

export function phyloUcaHmmLogMarginal(surface: ConditionalLikelihoodSurface, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): number {
  return forward(surface, references, options, false).logLikelihood;
}

function conditionalCharacterProbabilities(surface: ConditionalLikelihoodSurface, emissions: EmissionCache, state: number, site: number): Float64Array {
  const prior = emissions.prior(state, site);
  const result = new Float64Array(5);
  const offset = site * surface.stateCount;
  let maximum = NEGATIVE_INFINITY;
  for (let character = 0; character < surface.stateCount; character += 1) {
    if (prior[character] > 0) maximum = Math.max(maximum, Math.log(prior[character]) + surface.logLikelihoods[offset + character]);
  }
  let total = 0;
  for (let character = 0; character < surface.stateCount; character += 1) {
    if (!(prior[character] > 0)) continue;
    result[character] = Math.exp(Math.log(prior[character]) + surface.logLikelihoods[offset + character] - maximum);
    total += result[character];
  }
  if (total > 0) for (let character = 0; character < 5; character += 1) result[character] /= total;
  return result;
}

interface TrackDescriptor {
  id: string;
  kind: PhyloUcaSegmentKind;
  label: string;
  call?: string;
  dOrdinal?: number;
  registrationOffset?: number;
  pure: boolean;
  characterIndex?: number;
}

interface TrackAccumulator extends TrackDescriptor {
  points: Map<number, Float64Array>;
}

function normalizedTemplateCharacter(state: HmmState, site: number): string {
  if (state.kind === "N") return "N";
  return (state.kind === "D" ? state.character : state.projection?.[site])?.toUpperCase().replace("U", "T").replace(".", "-") ?? "N";
}

/**
 * V and J projections have one fixed register. D paths retain a register key
 * (alignment site minus D-reference position), so each allele row remains a
 * pure template even while start/trim uncertainty is marginalized. Unknown
 * projected template columns are diverted to a mixed unresolved row rather
 * than contaminating an allele row with a nucleotide mixture.
 */
function trackDescriptor(state: HmmState, site: number): TrackDescriptor {
  if (state.kind === "N") {
    return {
      id: `N|${state.dUsed}`,
      kind: "N",
      label: `N${state.dUsed}`,
      dOrdinal: state.dUsed || undefined,
      pure: false,
    };
  }
  const characterIndex = CHARACTERS.indexOf(normalizedTemplateCharacter(state, site) as PhyloUcaCharacter);
  if (characterIndex < 0) {
    if (state.kind === "V" || state.kind === "J") return {
      id: `${state.kind}|${state.call ?? "?"}|terminal-padding`,
      kind: state.kind,
      label: `${state.kind} · ${state.call ?? "?"} · terminal alignment padding`,
      call: state.call,
      pure: false,
    };
    const boundaryLabel = `N · unresolved ${state.kind}${state.dOrdinal ?? ""}`;
    return {
      id: `N|boundary|${state.kind}|${state.dUsed}`,
      kind: "N",
      label: boundaryLabel,
      dOrdinal: state.dOrdinal,
      pure: false,
    };
  }
  if (state.kind === "D") {
    const registrationOffset = site - state.dPosition;
    const sign = registrationOffset >= 0 ? "+" : "";
    return {
      id: `D|${state.dOrdinal ?? state.dUsed}|${state.call ?? "?"}|${registrationOffset}`,
      kind: "D",
      label: `D${state.dOrdinal ?? state.dUsed} · ${state.call ?? "?"} · register ${sign}${registrationOffset}`,
      call: state.call,
      dOrdinal: state.dOrdinal,
      registrationOffset,
      pure: true,
      characterIndex,
    };
  }
  return {
    id: `${state.kind}|${state.call ?? "?"}`,
    kind: state.kind,
    label: `${state.kind} · ${state.call ?? "?"}`,
    call: state.call,
    pure: true,
    characterIndex,
  };
}

function trackAlleleGroup(descriptor: TrackDescriptor): string {
  return descriptor.call ? `${descriptor.kind}|${descriptor.dOrdinal ?? 0}|${descriptor.call}` : descriptor.id;
}

function selectedMarginalTrackIds(
  metadata: ReadonlyMap<string, TrackDescriptor>,
  trackMaximums: ReadonlyMap<string, number>,
  alleleGroupMaximums: ReadonlyMap<string, number>,
): Set<string> {
  const selected = new Set<string>();
  const selectedGroups = new Set<string>();
  const strongestByGroup = new Map<string, { id: string; maximum: number }>();
  for (const descriptor of metadata.values()) {
    const maximum = trackMaximums.get(descriptor.id) ?? 0;
    const group = trackAlleleGroup(descriptor);
    const previous = strongestByGroup.get(group);
    if (!previous || maximum > previous.maximum) strongestByGroup.set(group, { id: descriptor.id, maximum });
    if (!descriptor.call) {
      if (maximum >= PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT) {
        selected.add(descriptor.id);
        selectedGroups.add(group);
      }
      continue;
    }
    if ((alleleGroupMaximums.get(group) ?? 0) < PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT) continue;
    if (descriptor.kind !== "D" || maximum >= PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT) {
      selected.add(descriptor.id);
      selectedGroups.add(group);
    }
  }
  // A D allele can exceed the group threshold while its mass is split across
  // many individually tiny registers. Retain the strongest pure register so a
  // non-negligible allele is never absent from the annotation.
  for (const [group, maximum] of alleleGroupMaximums) {
    if (maximum < PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT) continue;
    if (selectedGroups.has(group)) continue;
    const strongest = strongestByGroup.get(group);
    if (strongest) selected.add(strongest.id);
  }
  return selected;
}

function addTrackMass(
  tracks: Map<string, TrackAccumulator>,
  state: HmmState,
  site: number,
  weight: number,
  conditional: Float64Array,
): void {
  if (!(weight > 0) || !Number.isFinite(weight)) return;
  const descriptor = trackDescriptor(state, site);
  let track = tracks.get(descriptor.id);
  if (!track) {
    track = { ...descriptor, points: new Map() };
    tracks.set(descriptor.id, track);
  }
  let masses = track.points.get(site);
  if (!masses) {
    masses = new Float64Array(5);
    track.points.set(site, masses);
  }
  if (descriptor.pure && descriptor.characterIndex !== undefined) masses[descriptor.characterIndex] += weight;
  else for (let character = 0; character < 5; character += 1) masses[character] += weight * conditional[character];
}

function finalizeTracks(tracks: Map<string, TrackAccumulator>, siteTotals: Float64Array): PhyloUcaHmmAnnotationTrack[] {
  const kindOrder: Record<PhyloUcaSegmentKind, number> = { V: 0, N: 1, D: 2, J: 3, unknown: 4 };
  return [...tracks.values()].map((track) => {
    let maximumWeight = 0;
    const points = [...track.points.entries()].sort(([left], [right]) => left - right).flatMap(([site, raw]) => {
      const denominator = siteTotals[site];
      if (!(denominator > 0)) return [];
      const probabilities = Array.from(raw, (value) => value / denominator) as [number, number, number, number, number];
      const total = probabilities.reduce((sum, value) => sum + value, 0);
      if (!(total > 1e-14)) return [];
      maximumWeight = Math.max(maximumWeight, total);
      return [{ alignmentColumn: site + 1, probabilities }];
    });
    return {
      id: track.id,
      kind: track.kind,
      label: track.label,
      call: track.call,
      dOrdinal: track.dOrdinal,
      registrationOffset: track.registrationOffset,
      pure: track.pure,
      points,
      maximumWeight,
    };
  }).filter((track) => track.points.length > 0).sort((left, right) => {
    const kind = kindOrder[left.kind] - kindOrder[right.kind];
    if (kind) return kind;
    const ordinal = (left.dOrdinal ?? 0) - (right.dOrdinal ?? 0);
    if (ordinal) return ordinal;
    const start = (left.points[0]?.alignmentColumn ?? 0) - (right.points[0]?.alignmentColumn ?? 0);
    return start || left.label.localeCompare(right.label);
  });
}

/** Log q(character | HMM state, phylogenetic column data), character-major. */
function conditionalCharacterLogMatrix(surface: ConditionalLikelihoodSurface, emissions: EmissionCache, site: number, stateTotal: number): Float64Array {
  const result = new Float64Array(5 * stateTotal);
  result.fill(NEGATIVE_INFINITY);
  for (let state = 0; state < stateTotal; state += 1) {
    const conditional = conditionalCharacterProbabilities(surface, emissions, state, site);
    for (let character = 0; character < 5; character += 1) {
      if (conditional[character] > 0) result[character * stateTotal + state] = Math.log(conditional[character]);
    }
  }
  return result;
}

function addSiteEmission(values: Float64Array, emissions: EmissionCache, site: number, stateTotal: number): void {
  for (let state = 0; state < stateTotal; state += 1) {
    if (values[state] !== NEGATIVE_INFINITY) values[state] += emissionValue(emissions, state, site);
  }
}

/**
 * Exact P(x_i,x_{i+1},x_{i+2} | tree data, HMM, placement).
 *
 * The alpha row already contains the integrated emission at i. Multiplying it
 * by q_i(x_i | state_i,data_i), advancing the ordinary HMM, and finally
 * contracting against beta at i+2 fixes the three emitted characters while
 * summing over every germline candidate and recombination path. This preserves
 * the within-codon dependence that a product of site marginals discards.
 */
function exactCodonPosterior(
  surface: ConditionalLikelihoodSurface,
  emissions: EmissionCache,
  rows: Float64Array,
  betaAtEnd: Float64Array,
  startSite: number,
  context: TransitionContext,
): number[] {
  const stateTotal = context.catalog.states.length;
  const firstLogs = conditionalCharacterLogMatrix(surface, emissions, startSite, stateTotal);
  const secondLogs = conditionalCharacterLogMatrix(surface, emissions, startSite + 1, stateTotal);
  const thirdLogs = conditionalCharacterLogMatrix(surface, emissions, startSite + 2, stateTotal);
  const codonLogs = new Float64Array(PHYLO_UCA_CODON_STATE_COUNT);
  codonLogs.fill(NEGATIVE_INFINITY);
  const firstRow = startSite * stateTotal;

  const activeCharacters = surface.stateCount;
  for (let first = 0; first < activeCharacters; first += 1) {
    const fixedFirst = new Float64Array(stateTotal);
    for (let state = 0; state < stateTotal; state += 1) fixedFirst[state] = rows[firstRow + state] + firstLogs[first * stateTotal + state];
    const atSecond = transitionForward(fixedFirst, startSite, context, "sum");
    addSiteEmission(atSecond, emissions, startSite + 1, stateTotal);

    for (let second = 0; second < activeCharacters; second += 1) {
      const fixedSecond = new Float64Array(stateTotal);
      for (let state = 0; state < stateTotal; state += 1) fixedSecond[state] = atSecond[state] + secondLogs[second * stateTotal + state];
      const atThird = transitionForward(fixedSecond, startSite + 1, context, "sum");
      addSiteEmission(atThird, emissions, startSite + 2, stateTotal);

      for (let third = 0; third < activeCharacters; third += 1) {
        let score = NEGATIVE_INFINITY;
        const logOffset = third * stateTotal;
        for (let state = 0; state < stateTotal; state += 1) {
          score = logAdd(score, atThird[state] + thirdLogs[logOffset + state] + betaAtEnd[state]);
        }
        codonLogs[phyloUcaCodonStateIndex(first, second, third)] = score;
      }
    }
  }

  let maximum = NEGATIVE_INFINITY;
  for (const score of codonLogs) maximum = Math.max(maximum, score);
  if (!Number.isFinite(maximum)) throw new Error(`The exact UCA codon posterior at alignment columns ${startSite + 1}-${startSite + 3} has no finite probability mass.`);
  const masses = Array.from(codonLogs, (score) => Number.isFinite(score) ? Math.exp(score - maximum) : 0);
  return normalizeProbabilityVector(masses);
}

function pathSegments(states: readonly HmmState[], path: Int32Array, sequence: string): PhyloUcaPathSegment[] {
  const segments: PhyloUcaPathSegment[] = [];
  for (let column = 0; column < path.length; column += 1) {
    const state = states[path[column]];
    const previous = segments.at(-1);
    if (previous && previous.kind === state.kind && previous.call === state.call && previous.dOrdinal === state.dOrdinal) {
      previous.endColumn = column;
      previous.alignedSequence += sequence[column];
    } else segments.push({
      kind: state.kind,
      call: state.call,
      dOrdinal: state.dOrdinal,
      startColumn: column,
      endColumn: column,
      alignedSequence: sequence[column],
    });
  }
  return segments;
}

export function phyloUcaHmmPosterior(
  surface: ConditionalLikelihoodSurface,
  references: PreparedPhyloUcaReferences,
  options: PhyloUcaHmmOptions,
  frameOffset: PhyloUcaFrameOffset = 0,
): PhyloUcaHmmPosterior {
  const summed = forward(surface, references, options, true, "sum");
  const viterbi = forward(surface, references, options, false, "max");
  const { catalog, emissions, rows } = summed;
  if (!rows || !viterbi.backpointers || viterbi.endState === undefined || viterbi.endState < 0) throw new Error("The UCA recombination HMM found no complete V-to-J path.");
  const stateTotal = catalog.states.length;
  const logZWithoutSiteOffsets = summed.logLikelihood - emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
  const posterior = Array.from({ length: surface.sites }, () => [0, 0, 0, 0, 0] as [number, number, number, number, number]);
  const codonPosterior: Array<{ startSite: number; probabilities: number[] }> = [];
  const marginalTrackMetadata = new Map<string, TrackDescriptor>();
  const marginalTrackMaximums = new Map<string, number>();
  const marginalAlleleGroupMaximums = new Map<string, number>();
  let beta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
  beta.fill(NEGATIVE_INFINITY);
  for (const state of catalog.j) beta[state] = 0;
  const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
  for (let site = surface.sites - 1; site >= 0; site -= 1) {
    const rowOffset = site * stateTotal;
    const siteTrackWeights = new Map<string, number>();
    const siteAlleleGroupWeights = new Map<string, number>();
    for (let state = 0; state < stateTotal; state += 1) {
      const gammaLog = rows[rowOffset + state] + beta[state] - logZWithoutSiteOffsets;
      if (!Number.isFinite(gammaLog) || gammaLog < -745) continue;
      const weight = Math.exp(gammaLog);
      const conditional = conditionalCharacterProbabilities(surface, emissions, state, site);
      for (let character = 0; character < 5; character += 1) posterior[site][character] += weight * conditional[character];
      if (weight >= 1e-8) {
        const descriptor = trackDescriptor(catalog.states[state], site);
        marginalTrackMetadata.set(descriptor.id, descriptor);
        siteTrackWeights.set(descriptor.id, (siteTrackWeights.get(descriptor.id) ?? 0) + weight);
        const group = trackAlleleGroup(descriptor);
        siteAlleleGroupWeights.set(group, (siteAlleleGroupWeights.get(group) ?? 0) + weight);
      }
    }
    for (const [id, weight] of siteTrackWeights) marginalTrackMaximums.set(id, Math.max(marginalTrackMaximums.get(id) ?? 0, weight));
    for (const [group, weight] of siteAlleleGroupWeights) marginalAlleleGroupMaximums.set(group, Math.max(marginalAlleleGroupMaximums.get(group) ?? 0, weight));
    const total = posterior[site].reduce((sum, value) => sum + value, 0);
    if (!(total > 0) || !Number.isFinite(total)) throw new Error(`The UCA HMM posterior at alignment column ${site + 1} has no finite probability mass.`);
    posterior[site] = normalizeProbabilityVector(posterior[site]) as [number, number, number, number, number];
    if (site >= frameOffset + 2 && (site - frameOffset) % 3 === 2) {
      codonPosterior.push({
        startSite: site - 2,
        probabilities: exactCodonPosterior(surface, emissions, rows, beta, site - 2, context),
      });
    }
    if (site > 0) {
      const future = new Float64Array(stateTotal);
      for (let state = 0; state < stateTotal; state += 1) future[state] = emissionValue(emissions, state, site) + beta[state];
      beta = transitionBackward(future, site - 1, context);
    }
  }
  codonPosterior.reverse();
  const selectedTrackIds = selectedMarginalTrackIds(marginalTrackMetadata, marginalTrackMaximums, marginalAlleleGroupMaximums);
  const omittedMarginalTrackIds = [...marginalTrackMetadata.keys()].filter((id) => !selectedTrackIds.has(id));
  let omittedMarginalMaximumWeight = 0;
  for (const id of omittedMarginalTrackIds) omittedMarginalMaximumWeight = Math.max(omittedMarginalMaximumWeight, marginalTrackMaximums.get(id) ?? 0);
  const marginalTrackAccumulator = new Map<string, TrackAccumulator>();
  const marginalTrackSiteTotals = new Float64Array(surface.sites);
  beta = new Float64Array(stateTotal);
  beta.fill(NEGATIVE_INFINITY);
  for (const state of catalog.j) beta[state] = 0;
  for (let site = surface.sites - 1; site >= 0; site -= 1) {
    const rowOffset = site * stateTotal;
    let stateWeightTotal = 0;
    for (let state = 0; state < stateTotal; state += 1) {
      const gammaLog = rows[rowOffset + state] + beta[state] - logZWithoutSiteOffsets;
      if (!Number.isFinite(gammaLog) || gammaLog < -745) continue;
      const weight = Math.exp(gammaLog);
      stateWeightTotal += weight;
      if (weight < 1e-8) continue;
      const descriptor = trackDescriptor(catalog.states[state], site);
      if (!selectedTrackIds.has(descriptor.id)) continue;
      const conditional = conditionalCharacterProbabilities(surface, emissions, state, site);
      addTrackMass(marginalTrackAccumulator, catalog.states[state], site, weight, conditional);
    }
    marginalTrackSiteTotals[site] = stateWeightTotal;
    if (site > 0) {
      const future = new Float64Array(stateTotal);
      for (let state = 0; state < stateTotal; state += 1) future[state] = emissionValue(emissions, state, site) + beta[state];
      beta = transitionBackward(future, site - 1, context);
    }
  }
  const statePath = new Int32Array(surface.sites);
  statePath[surface.sites - 1] = viterbi.endState;
  for (let site = surface.sites - 1; site > 0; site -= 1) {
    const prior = viterbi.backpointers[site * stateTotal + statePath[site]];
    if (prior < 0) throw new Error(`The UCA Viterbi traceback stopped at alignment column ${site + 1}.`);
    statePath[site - 1] = prior;
  }
  let mapAlignedSequence = "";
  let posteriorConsensusAligned = "";
  const stateKinds: PhyloUcaSegmentKind[] = [];
  const stateCalls: Array<string | undefined> = [];
  const stateDOrdinals: Array<number | undefined> = [];
  const viterbiTrackAccumulator = new Map<string, TrackAccumulator>();
  const viterbiTrackSiteTotals = new Float64Array(surface.sites);
  for (let site = 0; site < surface.sites; site += 1) {
    const stateIndex = statePath[site];
    const conditional = conditionalCharacterProbabilities(surface, emissions, stateIndex, site);
    viterbiTrackSiteTotals[site] = 1;
    addTrackMass(viterbiTrackAccumulator, catalog.states[stateIndex], site, 1, conditional);
    let jointCharacter = 0;
    for (let character = 1; character < 5; character += 1) if (conditional[character] > conditional[jointCharacter]) jointCharacter = character;
    mapAlignedSequence += CHARACTERS[jointCharacter];
    let marginalCharacter = 0;
    for (let character = 1; character < 5; character += 1) if (posterior[site][character] > posterior[site][marginalCharacter]) marginalCharacter = character;
    posteriorConsensusAligned += CHARACTERS[marginalCharacter];
    const state = catalog.states[stateIndex];
    stateKinds.push(state.kind);
    stateCalls.push(state.call);
    stateDOrdinals.push(state.dOrdinal);
  }
  const segments = pathSegments(catalog.states, statePath, mapAlignedSequence);
  const viterbiTracks = finalizeTracks(viterbiTrackAccumulator, viterbiTrackSiteTotals);
  const marginalTracks = finalizeTracks(marginalTrackAccumulator, marginalTrackSiteTotals);
  return {
    logMarginalLikelihood: summed.logLikelihood,
    probabilities: posterior,
    codonPosterior,
    mapAlignedSequence,
    posteriorConsensusAligned,
    stateKinds,
    stateCalls,
    stateDOrdinals,
    path: segments,
    viterbiTracks,
    marginalTracks,
    omittedMarginalTrackIds,
    omittedMarginalMaximumWeight,
    mapVCall: segments.find((segment) => segment.kind === "V")?.call ?? "",
    mapDCalls: segments.filter((segment) => segment.kind === "D").map((segment) => segment.call ?? ""),
    mapJCall: segments.find((segment) => segment.kind === "J")?.call ?? "",
  };
}
