import type { ConditionalLikelihoodSurface } from "./tree-messages.ts";
import type {
  PhyloUcaCharacter,
  PhyloUcaHmmOptions,
  PhyloUcaPathSegment,
  PhyloUcaSegmentKind,
} from "./types.ts";
import type { PreparedPhyloUcaReferences } from "./references.ts";

const NEGATIVE_INFINITY = Number.NEGATIVE_INFINITY;
const PROBABILITY_FLOOR = 1e-300;
const CHARACTERS = ["A", "C", "G", "T", "-"] as const;

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
  mapAlignedSequence: string;
  posteriorConsensusAligned: string;
  stateKinds: PhyloUcaSegmentKind[];
  stateCalls: Array<string | undefined>;
  stateDOrdinals: Array<number | undefined>;
  path: PhyloUcaPathSegment[];
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

export function phyloUcaHmmPosterior(surface: ConditionalLikelihoodSurface, references: PreparedPhyloUcaReferences, options: PhyloUcaHmmOptions): PhyloUcaHmmPosterior {
  const summed = forward(surface, references, options, true, "sum");
  const viterbi = forward(surface, references, options, false, "max");
  const { catalog, emissions, rows } = summed;
  if (!rows || !viterbi.backpointers || viterbi.endState === undefined || viterbi.endState < 0) throw new Error("The UCA recombination HMM found no complete V-to-J path.");
  const stateTotal = catalog.states.length;
  const logZWithoutSiteOffsets = summed.logLikelihood - emissions.siteOffsets.reduce((sum, value) => sum + value, 0);
  const posterior = Array.from({ length: surface.sites }, () => [0, 0, 0, 0, 0] as [number, number, number, number, number]);
  let beta: Float64Array<ArrayBufferLike> = new Float64Array(stateTotal);
  beta.fill(NEGATIVE_INFINITY);
  for (const state of catalog.j) beta[state] = 0;
  const context = { catalog, references, options, window: junctionWindow(surface.sites, references, options) };
  for (let site = surface.sites - 1; site >= 0; site -= 1) {
    const rowOffset = site * stateTotal;
    for (let state = 0; state < stateTotal; state += 1) {
      const gammaLog = rows[rowOffset + state] + beta[state] - logZWithoutSiteOffsets;
      if (gammaLog < -40 || !Number.isFinite(gammaLog)) continue;
      const weight = Math.exp(gammaLog);
      const conditional = conditionalCharacterProbabilities(surface, emissions, state, site);
      for (let character = 0; character < 5; character += 1) posterior[site][character] += weight * conditional[character];
    }
    const total = posterior[site].reduce((sum, value) => sum + value, 0);
    if (total > 0) for (let character = 0; character < 5; character += 1) posterior[site][character] /= total;
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
  for (let site = 0; site < surface.sites; site += 1) {
    const stateIndex = statePath[site];
    const conditional = conditionalCharacterProbabilities(surface, emissions, stateIndex, site);
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
  return {
    logMarginalLikelihood: summed.logLikelihood,
    probabilities: posterior,
    mapAlignedSequence,
    posteriorConsensusAligned,
    stateKinds,
    stateCalls,
    stateDOrdinals,
    path: segments,
    mapVCall: segments.find((segment) => segment.kind === "V")?.call ?? "",
    mapDCalls: segments.filter((segment) => segment.kind === "D").map((segment) => segment.call ?? ""),
    mapJCall: segments.find((segment) => segment.kind === "J")?.call ?? "",
  };
}
