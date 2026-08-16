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
const PROBABILITY_FLOOR = 1e-300;
const CHARACTERS = ["A", "C", "G", "T", "-"] as const;
export const PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT = 0.01;
export const PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT = 0.001;

interface HmmState {
  kind: Exclude<PhyloUcaSegmentKind, "unknown">;
  call?: string;
  dOrdinal?: number;
  projection?: string;
  character?: string;
  dUsed: number;
  dRun: number;
  dPosition: number;
}

interface StateCatalog {
  states: HmmState[];
  v: number[];
  n: number[];
  dByCount: number[][];
  dEntries: number[][];
  dEntryLogPrior: Float64Array;
  dContinue: Int32Array;
  j: number[];
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

function clampProbability(value: number, fallback: number): number {
  if (!Number.isFinite(value)) return fallback;
  return Math.max(1e-9, Math.min(1 - 1e-9, value));
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

function sigmoid(value: number): number {
  if (value >= 0) return 1 / (1 + Math.exp(-value));
  const exponential = Math.exp(value);
  return exponential / (1 + exponential);
}

function normalizedLogs(weights: readonly number[]): number[] {
  const safe = weights.map((weight) => Number.isFinite(weight) && weight > 0 ? weight : 0);
  const total = safe.reduce((sum, weight) => sum + weight, 0);
  if (!(total > 0)) return safe.map(() => NEGATIVE_INFINITY);
  return safe.map((weight) => weight > 0 ? Math.log(weight / total) : NEGATIVE_INFINITY);
}

function buildStateCatalog(references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): StateCatalog {
  const states: HmmState[] = [];
  const v = references.v.map((candidate) => {
    const index = states.length;
    states.push({ kind: "V", call: candidate.name, projection: candidate.projection, dUsed: 0, dRun: 0, dPosition: -1 });
    return index;
  });
  const maximumD = Math.max(0, Math.min(5, Math.floor(options.maximumDSegments)));
  const n = Array.from({ length: maximumD + 1 }, (_, dUsed) => {
    const index = states.length;
    states.push({ kind: "N", dUsed, dRun: 0, dPosition: -1 });
    return index;
  });
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
          states.push({
            kind: "D",
            call: record.name,
            character: record.sequence[position],
            dOrdinal: dUsed,
            dUsed,
            dRun: run,
            dPosition: position,
          });
          dKeys.set(`${dUsed}|${record.name}|${position}|${run}`, index);
          dByCount[dUsed].push(index);
          if (run === 1 && record.sequence.length - position >= minimumDMatch) dEntries[dUsed].push(index);
        }
      }
    }
  }
  const j = references.j.map((candidate) => {
    const index = states.length;
    states.push({ kind: "J", call: candidate.name, projection: candidate.projection, dUsed: 0, dRun: 0, dPosition: -1 });
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
  const continuation = clampProbability(options.dFivePrimeTrimContinuation, 0.72);
  for (let dUsed = 1; dUsed <= maximumD; dUsed += 1) {
    const raw = dEntries[dUsed].map((stateIndex) => {
      const state = states[stateIndex];
      return Math.pow(continuation, state.dPosition);
    });
    const total = raw.reduce((sum, value) => sum + value, 0);
    dEntries[dUsed].forEach((stateIndex, entry) => {
      dEntryLogPrior[stateIndex] = Math.log(Math.max(PROBABILITY_FLOOR, raw[entry] / Math.max(PROBABILITY_FLOOR, total)));
    });
  }
  return { states, v, n, dByCount, dEntries, dEntryLogPrior, dContinue, j, minimumDMatch };
}

function statePrior(state: HmmState, site: number, stateCount: 4 | 5, options: PhyloUcaHmmOptions): Float64Array {
  const result = new Float64Array(stateCount);
  const baseRaw = options.nBaseFrequencies.map((value) => Number.isFinite(value) && value > 0 ? value : 1e-9);
  const baseTotal = baseRaw.reduce((sum, value) => sum + value, 0);
  const bases = baseRaw.map((value) => value / baseTotal);
  const exact = state.kind === "N" ? "N" : state.kind === "D" ? state.character ?? "N" : state.projection?.[site] ?? "N";
  const normalized = exact.toUpperCase().replace("U", "T").replace(".", "-");
  const exactIndex = CHARACTERS.indexOf(normalized as PhyloUcaCharacter);
  const mismatch = Math.max(1e-9, Math.min(0.25, options.templateMismatchProbability));
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
    ? Math.max(1e-9, Math.min(0.5, state.kind === "N" ? options.junctionGapProbability : options.unknownTemplateGapProbability))
    : 0;
  for (let index = 0; index < 4; index += 1) result[index] = (1 - gapProbability) * bases[index];
  if (stateCount === 5) result[4] = gapProbability;
  return result;
}

const TEMPLATE_UNKNOWN_CATEGORY = 5;
const JUNCTION_CATEGORY = 6;

function stateCategory(state: HmmState, site: number): number {
  if (state.kind === "N") return JUNCTION_CATEGORY;
  const character = (state.kind === "D" ? state.character : state.projection?.[site])?.toUpperCase().replace("U", "T").replace(".", "-") ?? "N";
  const exact = CHARACTERS.indexOf(character as PhyloUcaCharacter);
  return exact >= 0 ? exact : TEMPLATE_UNKNOWN_CATEGORY;
}

function categoryPriors(stateCount: 4 | 5, options: PhyloUcaHmmOptions): Float64Array[] {
  const template = (character: string): HmmState => ({ kind: "D", character, dUsed: 1, dRun: 1, dPosition: 0 });
  const junction: HmmState = { kind: "N", dUsed: 0, dRun: 0, dPosition: -1 };
  return ["A", "C", "G", "T", "-", "N"].map((character) => statePrior(template(character), 0, stateCount, options)).concat([statePrior(junction, 0, stateCount, options)]);
}

interface EmissionCache {
  values: Float64Array;
  siteOffsets: Float64Array;
  prior: (state: number, site: number) => Float64Array;
}

interface JunctionWindow { start: number; end: number }

function junctionWindow(sites: number, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): JunctionWindow {
  // D/N uncertainty is local to the V–J junction. Keeping the all-D automaton
  // inactive elsewhere is exact under the fixed boundary prior up to this
  // deliberately generous (>4 logistic scale) truncation and avoids scanning
  // tens of thousands of D states across framework columns.
  const left = Math.ceil(4 * Math.max(0.25, options.vTrimScale)) + 4;
  const right = Math.ceil(4 * Math.max(0.25, options.jTrimScale)) + 4;
  return {
    start: Math.max(0, references.vEndColumn - left),
    end: Math.min(sites - 1, references.jStartColumn + right),
  };
}

function buildEmissionCache(surface: ConditionalLikelihoodSurface, catalog: StateCatalog, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): EmissionCache {
  const stateCount = surface.stateCount;
  const stateTotal = catalog.states.length;
  const values = new Float64Array(surface.sites * stateTotal);
  values.fill(NEGATIVE_INFINITY);
  const siteOffsets = new Float64Array(surface.sites);
  const priors = categoryPriors(stateCount, options);
  const prior = (state: number, site: number) => priors[stateCategory(catalog.states[state], site)];
  const window = junctionWindow(surface.sites, references, options);
  const dStates = catalog.dByCount.flat();
  const nonDStates = catalog.states.map((state, index) => state.kind === "D" ? -1 : index).filter((index) => index >= 0);
  for (let site = 0; site < surface.sites; site += 1) {
    const likelihoodOffset = site * stateCount;
    let siteMaximum = NEGATIVE_INFINITY;
    for (let character = 0; character < stateCount; character += 1) siteMaximum = Math.max(siteMaximum, surface.logLikelihoods[likelihoodOffset + character]);
    siteOffsets[site] = siteMaximum;
    const categoryEmissions = new Float64Array(priors.length);
    categoryEmissions.fill(NEGATIVE_INFINITY);
    for (let category = 0; category < priors.length; category += 1) {
      const probabilities = priors[category];
      let emission = NEGATIVE_INFINITY;
      for (let character = 0; character < stateCount; character += 1) {
        if (probabilities[character] > 0) emission = logAdd(emission, Math.log(probabilities[character]) + surface.logLikelihoods[likelihoodOffset + character] - siteMaximum);
      }
      categoryEmissions[category] = emission;
    }
    for (const state of nonDStates) values[site * stateTotal + state] = categoryEmissions[stateCategory(catalog.states[state], site)];
    if (site >= window.start && site <= window.end) for (const state of dStates) values[site * stateTotal + state] = categoryEmissions[stateCategory(catalog.states[state], site)];
  }
  return { values, siteOffsets, prior };
}

interface TransitionContext {
  catalog: StateCatalog;
  references: PreparedPhyloUcaReferences;
  options: PhyloUcaHmmOptions;
  window: JunctionWindow;
}

function routingLogs(context: TransitionContext, site: number, source: "V" | "N" | "D", dUsed: number): [number, number, number] {
  const { references, options, catalog } = context;
  const readiness = sigmoid((site - references.jStartColumn) / Math.max(0.25, options.jTrimScale));
  const canD = site + 1 >= context.window.start && site + 1 <= context.window.end && dUsed < catalog.dEntries.length - 1 && catalog.dEntries[dUsed + 1].length > 0;
  if (source === "V") return normalizedLogs([0.24, canD ? 0.72 * (1 - 0.65 * readiness) : 0, 0.015 + 1.8 * readiness]) as [number, number, number];
  if (source === "N") {
    const dWeight = !canD ? 0 : dUsed === 0 ? 0.8 * (1 - 0.7 * readiness) : options.additionalDProbability * (1.2 - readiness);
    return normalizedLogs([0, dWeight, 0.08 + 2.2 * readiness]) as [number, number, number];
  }
  const dWeight = canD ? options.additionalDProbability * (1.2 - readiness) : 0;
  return normalizedLogs([0.65, dWeight, 0.1 + 2.1 * readiness]) as [number, number, number];
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

function distributeJ(destination: Float64Array, context: TransitionContext, score: number, routeLog: number, mode: CombineMode, sources?: Int32Array, source = -1): void {
  if (score === NEGATIVE_INFINITY || routeLog === NEGATIVE_INFINITY || !context.catalog.j.length) return;
  const candidateLog = -Math.log(context.catalog.j.length);
  for (const target of context.catalog.j) combineDestination(destination, target, score + routeLog + candidateLog, mode, sources, source);
}

function transitionForward(sourceValues: Float64Array, site: number, context: TransitionContext, mode: CombineMode, sources?: Int32Array): Float64Array {
  const { catalog, references, options, window } = context;
  const destination = new Float64Array(catalog.states.length);
  destination.fill(NEGATIVE_INFINITY);
  sources?.fill(-1);
  const vExit = Math.max(1e-7, Math.min(0.995, sigmoid((site - references.vEndColumn) / Math.max(0.25, options.vTrimScale))));
  const vHubValues = new Float64Array(catalog.v.length);
  for (let entry = 0; entry < catalog.v.length; entry += 1) {
    const state = catalog.v[entry];
    combineDestination(destination, state, sourceValues[state] + Math.log1p(-vExit), mode, sources, state);
    vHubValues[entry] = sourceValues[state] + Math.log(vExit);
  }
  const vHub = aggregate(vHubValues, catalog.v.map((_, index) => index), mode);
  const vSource = vHub.source >= 0 ? catalog.v[vHub.source] : -1;
  const vRoutes = routingLogs(context, site, "V", 0);
  combineDestination(destination, catalog.n[0], vHub.score + vRoutes[0], mode, sources, vSource);
  distributeD(destination, context, 1, vHub.score, vRoutes[1], mode, sources, vSource);
  distributeJ(destination, context, vHub.score, vRoutes[2], mode, sources, vSource);

  const nStay = Math.max(1e-7, Math.min(0.999, Math.max(0, options.meanNLength) / (Math.max(0, options.meanNLength) + 1)));
  for (let dUsed = 0; dUsed < catalog.n.length; dUsed += 1) {
    const state = catalog.n[dUsed];
    combineDestination(destination, state, sourceValues[state] + Math.log(nStay), mode, sources, state);
    const routes = routingLogs(context, site, "N", dUsed);
    const hub = sourceValues[state] + Math.log1p(-nStay);
    distributeD(destination, context, dUsed + 1, hub, routes[1], mode, sources, state);
    distributeJ(destination, context, hub, routes[2], mode, sources, state);
  }

  const configuredDExit = clampProbability(options.dExitProbability, 0.28);
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
      const exitProbability = next < 0 || site >= window.end ? 1 : configuredDExit;
      if (next >= 0) combineDestination(destination, next, sourceValues[stateIndex] + Math.log1p(-exitProbability), mode, sources, stateIndex);
      const candidate = sourceValues[stateIndex] + Math.log(exitProbability);
      if (mode === "sum") exitScore = logAdd(exitScore, candidate);
      else if (candidate > exitScore) {
        exitScore = candidate;
        exitSource = stateIndex;
      }
    }
    const routes = routingLogs(context, site, "D", dUsed);
    combineDestination(destination, catalog.n[dUsed], exitScore + routes[0], mode, sources, exitSource);
    distributeD(destination, context, dUsed + 1, exitScore, routes[1], mode, sources, exitSource);
    distributeJ(destination, context, exitScore, routes[2], mode, sources, exitSource);
  }
  for (const state of catalog.j) combineDestination(destination, state, sourceValues[state], mode, sources, state);
  return destination;
}

function weightedDestination(values: Float64Array, indexes: readonly number[], logWeights: (index: number) => number): number {
  let result = NEGATIVE_INFINITY;
  for (const index of indexes) result = logAdd(result, values[index] + logWeights(index));
  return result;
}

function transitionBackward(future: Float64Array, site: number, context: TransitionContext): Float64Array {
  const { catalog, references, options, window } = context;
  const source = new Float64Array(catalog.states.length);
  source.fill(NEGATIVE_INFINITY);
  const jFuture = weightedDestination(future, catalog.j, () => -Math.log(catalog.j.length));
  const dFuture = catalog.dEntries.map((entries) => weightedDestination(future, entries, (index) => catalog.dEntryLogPrior[index]));
  const vExit = Math.max(1e-7, Math.min(0.995, sigmoid((site - references.vEndColumn) / Math.max(0.25, options.vTrimScale))));
  const vRoutes = routingLogs(context, site, "V", 0);
  const vHubFuture = logSum([
    vRoutes[0] + future[catalog.n[0]],
    vRoutes[1] + (dFuture[1] ?? NEGATIVE_INFINITY),
    vRoutes[2] + jFuture,
  ]);
  for (const state of catalog.v) source[state] = logAdd(Math.log1p(-vExit) + future[state], Math.log(vExit) + vHubFuture);

  const nStay = Math.max(1e-7, Math.min(0.999, Math.max(0, options.meanNLength) / (Math.max(0, options.meanNLength) + 1)));
  for (let dUsed = 0; dUsed < catalog.n.length; dUsed += 1) {
    const state = catalog.n[dUsed];
    const routes = routingLogs(context, site, "N", dUsed);
    const exitFuture = logSum([routes[1] + (dFuture[dUsed + 1] ?? NEGATIVE_INFINITY), routes[2] + jFuture]);
    source[state] = logAdd(Math.log(nStay) + future[state], Math.log1p(-nStay) + exitFuture);
  }

  const configuredDExit = clampProbability(options.dExitProbability, 0.28);
  if (site >= window.start && site <= window.end) for (let dUsed = 1; dUsed < catalog.dByCount.length; dUsed += 1) {
    const routes = routingLogs(context, site, "D", dUsed);
    const exitFuture = logSum([
      routes[0] + future[catalog.n[dUsed]],
      routes[1] + (dFuture[dUsed + 1] ?? NEGATIVE_INFINITY),
      routes[2] + jFuture,
    ]);
    for (const stateIndex of catalog.dByCount[dUsed]) {
      const state = catalog.states[stateIndex];
      const next = catalog.dContinue[stateIndex];
      if (state.dRun < catalog.minimumDMatch) source[stateIndex] = next >= 0 ? future[next] : NEGATIVE_INFINITY;
      else {
        const exitProbability = next < 0 || site >= window.end ? 1 : configuredDExit;
        source[stateIndex] = next < 0
          ? exitFuture
          : logAdd(Math.log1p(-exitProbability) + future[next], Math.log(exitProbability) + exitFuture);
      }
    }
  }
  for (const state of catalog.j) source[state] = future[state];
  return source;
}

function initialize(catalog: StateCatalog, emissions: EmissionCache): Float64Array {
  const values = new Float64Array(catalog.states.length);
  values.fill(NEGATIVE_INFINITY);
  const prior = -Math.log(catalog.v.length);
  for (const state of catalog.v) values[state] = prior + emissions.values[state];
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
 * Reusable exact FFBS sampler for p(HMM path, UCA characters | placement, tree).
 * The large D-state catalog is constructed once. Each draw performs one
 * backward HMM pass, samples one coherent path, and samples UCA characters
 * from the tree/template conditional at the visited states.
 */
export class PhyloUcaHmmGibbsSampler {
  private readonly catalog: StateCatalog;
  private readonly references: PreparedPhyloUcaReferences;
  private readonly options: PhyloUcaHmmOptions;

  constructor(
    references: PreparedPhyloUcaReferences,
    options: PhyloUcaHmmOptions,
  ) {
    this.references = references;
    this.options = options;
    this.catalog = buildStateCatalog(references, options);
    if (!this.catalog.v.length || !this.catalog.j.length) throw new Error("The UCA HMM Gibbs sampler requires at least one V and one J candidate.");
  }

  draw(surface: ConditionalLikelihoodSurface, random: () => number = Math.random): PhyloUcaHmmGibbsDraw {
    const { catalog, references, options } = this;
    const emissions = buildEmissionCache(surface, catalog, references, options);
    const stateTotal = catalog.states.length;
    const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
    const betaRows = new Float64Array(surface.sites * stateTotal);
    betaRows.fill(NEGATIVE_INFINITY);
    const lastOffset = (surface.sites - 1) * stateTotal;
    for (const state of catalog.j) betaRows[lastOffset + state] = 0;
    let beta: Float64Array = betaRows.slice(lastOffset, lastOffset + stateTotal);
    for (let site = surface.sites - 2; site >= 0; site -= 1) {
      const future = new Float64Array(stateTotal);
      const emissionOffset = (site + 1) * stateTotal;
      for (let state = 0; state < stateTotal; state += 1) future[state] = emissions.values[emissionOffset + state] + beta[state];
      beta = transitionBackward(future, site, context);
      betaRows.set(beta, site * stateTotal);
    }

    const initial = initialize(catalog, emissions);
    const initialConditional = new Float64Array(stateTotal);
    for (let state = 0; state < stateTotal; state += 1) initialConditional[state] = initial[state] + betaRows[state];
    const logMarginalLikelihood = logSum(initialConditional) + emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
    const statePath = new Int32Array(surface.sites);
    statePath[0] = sampleLogCategorical(initialConditional, random);

    const addDestination = (destinations: number[], weights: number[], target: number, transitionLog: number, site: number) => {
      if (target < 0 || !Number.isFinite(transitionLog)) return;
      const future = emissions.values[(site + 1) * stateTotal + target] + betaRows[(site + 1) * stateTotal + target];
      if (!Number.isFinite(future)) return;
      destinations.push(target);
      weights.push(transitionLog + future);
    };

    for (let site = 0; site < surface.sites - 1; site += 1) {
      const sourceIndex = statePath[site];
      const source = catalog.states[sourceIndex];
      const destinations: number[] = [];
      const weights: number[] = [];
      if (source.kind === "V") {
        const vExit = Math.max(1e-7, Math.min(0.995, sigmoid((site - references.vEndColumn) / Math.max(0.25, options.vTrimScale))));
        addDestination(destinations, weights, sourceIndex, Math.log1p(-vExit), site);
        const routes = routingLogs(context, site, "V", 0);
        const hub = Math.log(vExit);
        addDestination(destinations, weights, catalog.n[0], hub + routes[0], site);
        for (const target of catalog.dEntries[1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
        const candidateLog = -Math.log(catalog.j.length);
        for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog, site);
      } else if (source.kind === "N") {
        const nStay = Math.max(1e-7, Math.min(0.999, Math.max(0, options.meanNLength) / (Math.max(0, options.meanNLength) + 1)));
        addDestination(destinations, weights, sourceIndex, Math.log(nStay), site);
        const routes = routingLogs(context, site, "N", source.dUsed);
        const hub = Math.log1p(-nStay);
        for (const target of catalog.dEntries[source.dUsed + 1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
        const candidateLog = -Math.log(catalog.j.length);
        for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog, site);
      } else if (source.kind === "D") {
        const next = catalog.dContinue[sourceIndex];
        if (source.dRun < catalog.minimumDMatch) addDestination(destinations, weights, next, 0, site);
        else {
          const configuredDExit = clampProbability(options.dExitProbability, 0.28);
          const exitProbability = next < 0 || site >= context.window.end ? 1 : configuredDExit;
          if (next >= 0) addDestination(destinations, weights, next, Math.log1p(-exitProbability), site);
          const routes = routingLogs(context, site, "D", source.dUsed);
          const hub = Math.log(exitProbability);
          addDestination(destinations, weights, catalog.n[source.dUsed], hub + routes[0], site);
          for (const target of catalog.dEntries[source.dUsed + 1] ?? []) addDestination(destinations, weights, target, hub + routes[1] + catalog.dEntryLogPrior[target], site);
          const candidateLog = -Math.log(catalog.j.length);
          for (const target of catalog.j) addDestination(destinations, weights, target, hub + routes[2] + candidateLog, site);
        }
      } else addDestination(destinations, weights, sourceIndex, 0, site);
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
    const emissionOffset = (site + 1) * stateTotal;
    for (let state = 0; state < stateTotal; state += 1) if (next[state] !== NEGATIVE_INFINITY) next[state] += emissions.values[emissionOffset + state];
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
    const boundaryLabel = state.kind === "V" ? "N · V-trim boundary" : state.kind === "J" ? "N · J-entry boundary" : `N · unresolved ${state.kind}${state.dOrdinal ?? ""}`;
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
  const offset = site * stateTotal;
  for (let state = 0; state < stateTotal; state += 1) {
    if (values[state] !== NEGATIVE_INFINITY) values[state] += emissions.values[offset + state];
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
      const emissionOffset = site * stateTotal;
      for (let state = 0; state < stateTotal; state += 1) future[state] = emissions.values[emissionOffset + state] + beta[state];
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
      const emissionOffset = site * stateTotal;
      for (let state = 0; state < stateTotal; state += 1) future[state] = emissions.values[emissionOffset + state] + beta[state];
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
