import type { PhyloUcaHmmOptions } from "./types.ts";

export interface PhyloUcaPriorPredictiveMetric {
  id: string;
  label: string;
  unit: "nt" | "segments";
  mean: number;
  median: number;
  p90: number;
  p95: number;
  observations: number;
}

export interface PhyloUcaPriorPredictiveSummary {
  draws: number;
  metrics: PhyloUcaPriorPredictiveMetric[];
  dCountProbabilities: number[];
}

function probability(value: number, fallback: number): number {
  return Number.isFinite(value) ? Math.max(0, Math.min(1, value)) : fallback;
}

function seededRandom(seed: number): () => number {
  let state = seed >>> 0;
  return () => {
    state += 0x6d2b79f5;
    let value = state;
    value = Math.imul(value ^ value >>> 15, value | 1);
    value ^= value + Math.imul(value ^ value >>> 7, value | 61);
    return ((value ^ value >>> 14) >>> 0) / 4294967296;
  };
}

/** Draw k >= 0 from (1-c)c^k. The HMM itself uses exact finite support. */
function geometricTail(continuationRaw: number, random: () => number): number {
  const continuation = Math.min(0.999999, probability(continuationRaw, 0));
  if (continuation <= 0) return 0;
  return Math.floor(Math.log1p(-random()) / Math.log(continuation));
}

/** Draw k in 0..maximum with probability proportional to c^k. */
function finiteGeometricTail(continuationRaw: number, maximum: number, random: () => number): number {
  const limit = Math.max(0, Math.floor(maximum));
  const continuation = Math.min(0.999999, probability(continuationRaw, 0));
  if (limit === 0 || continuation <= 0) return 0;
  let total = 0;
  let weight = 1;
  for (let value = 0; value <= limit; value += 1) {
    total += weight;
    weight *= continuation;
  }
  let threshold = random() * total;
  weight = 1;
  for (let value = 0; value <= limit; value += 1) {
    threshold -= weight;
    if (threshold <= 0) return value;
    weight *= continuation;
  }
  return limit;
}

function summarize(id: string, label: string, unit: "nt" | "segments", values: number[]): PhyloUcaPriorPredictiveMetric {
  if (!values.length) return { id, label, unit, mean: 0, median: 0, p90: 0, p95: 0, observations: 0 };
  const sorted = [...values].sort((left, right) => left - right);
  const quantile = (q: number) => sorted[Math.min(sorted.length - 1, Math.max(0, Math.ceil(q * sorted.length) - 1))];
  return {
    id,
    label,
    unit,
    mean: values.reduce((sum, value) => sum + value, 0) / values.length,
    median: quantile(0.5),
    p90: quantile(0.9),
    p95: quantile(0.95),
    observations: values.length,
  };
}

/**
 * Fast deterministic prior-predictive audit of the HMM route, finite D-trim,
 * and N-run distributions. V/J deletion summaries use their exposed
 * geometric tail moments because the panel passes active D lengths only. It
 * generates recombination scenarios, not tree sequences, and never enters
 * inference runtime.
 */
export function phyloUcaPriorPredictiveSummary(
  options: PhyloUcaHmmOptions,
  dReferenceLengths: readonly number[],
  draws = 50_000,
  seed = 0x51a6c9,
): PhyloUcaPriorPredictiveSummary {
  const random = seededRandom(seed);
  const maximumD = Math.max(0, Math.min(5, Math.floor(options.maximumDSegments)));
  const minimumD = Math.max(1, Math.min(12, Math.floor(options.minimumDMatch)));
  const dLengths = dReferenceLengths.filter((length) => Number.isFinite(length) && length >= minimumD).map(Math.floor);
  const initialD = dLengths.length && maximumD ? probability(options.initialDProbability, 0.934) : 0;
  const additionalD = probability(options.additionalDProbability, 0.00125);
  const nProbability = probability(options.junctionNProbability, 0.973);
  const nMean = Math.max(1, Number.isFinite(options.meanNLength) ? options.meanNLength : 8.8);
  const singleN = probability(options.singleNProbability, 0.027);
  const nPhases = Math.max(1, Math.min(4, Math.floor(options.nLengthPhases ?? 2)));
  const nTailMean = singleN < 1 ? Math.max(nPhases, (nMean - singleN) / (1 - singleN)) : nPhases;
  const nTailContinuation = 1 - nPhases / nTailMean;
  const samples: Record<string, number[]> = {
    v3: [], d5: [], d3: [], j5: [], retainedD: [], nRun: [], totalN: [], dCount: [], junction: [],
  };
  const dCountProbabilities = Array.from({ length: maximumD + 1 }, () => 0);

  for (let draw = 0; draw < Math.max(1, Math.floor(draws)); draw += 1) {
    const v3 = geometricTail(options.vThreePrimeTrimContinuation, random);
    const j5 = geometricTail(options.jFivePrimeTrimContinuation, random);
    let dCount = random() < initialD ? 1 : 0;
    while (dCount > 0 && dCount < maximumD && random() < additionalD) dCount += 1;
    dCountProbabilities[dCount] += 1;
    let retainedDTotal = 0;
    for (let d = 0; d < dCount; d += 1) {
      const length = dLengths[Math.floor(random() * dLengths.length)];
      // This is exactly the HMM's sequential finite-support prior: choose a
      // valid 5' start within the allele, then choose a valid 3' exit given
      // that start. No rejection fallback can distort either tail.
      const d5 = finiteGeometricTail(options.dFivePrimeTrimContinuation, length - minimumD, random);
      const d3 = finiteGeometricTail(options.dThreePrimeTrimContinuation, length - minimumD - d5, random);
      const retained = length - d5 - d3;
      samples.d5.push(d5);
      samples.d3.push(d3);
      samples.retainedD.push(retained);
      retainedDTotal += retained;
    }
    let totalN = 0;
    for (let junction = 0; junction < dCount + 1; junction += 1) if (random() < nProbability) {
      let length = 1;
      if (random() >= singleN) {
        length = 0;
        for (let phase = 0; phase < nPhases; phase += 1) length += 1 + geometricTail(nTailContinuation, random);
      }
      samples.nRun.push(length);
      totalN += length;
    }
    samples.v3.push(v3);
    samples.j5.push(j5);
    samples.totalN.push(totalN);
    samples.dCount.push(dCount);
    samples.junction.push(totalN + retainedDTotal);
  }
  const actualDraws = samples.dCount.length;
  for (let count = 0; count < dCountProbabilities.length; count += 1) dCountProbabilities[count] /= actualDraws;
  return {
    draws: actualDraws,
    dCountProbabilities,
    metrics: [
      summarize("v3", "V 3′ deletion", "nt", samples.v3),
      summarize("d5", "D 5′ deletion", "nt", samples.d5),
      summarize("retainedD", "Retained D template", "nt", samples.retainedD),
      summarize("d3", "D 3′ deletion", "nt", samples.d3),
      summarize("j5", "J 5′ deletion", "nt", samples.j5),
      summarize("nRun", "N per occupied junction", "nt", samples.nRun),
      summarize("totalN", "Total N across junctions", "nt", samples.totalN),
      summarize("dCount", "D segments", "segments", samples.dCount),
      summarize("junction", "Modeled V–J junction span", "nt", samples.junction),
    ],
  };
}
