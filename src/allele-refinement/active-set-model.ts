import type {
  AllelePosteriorSummary,
  AlleleRefinementOptions,
  ReferenceAlleleGraph,
  RefinementModelSummary,
  SegmentRefinementResult,
  SparseEvidenceMatrix,
} from "./types.ts";

interface EvidencePattern {
  columns: number[];
  evidence: number[];
  weight: number;
}

interface IncidentPattern {
  pattern: number;
  evidence: number;
}

interface MixtureFit {
  theta: Float64Array;
  counts: Float64Array;
  iterations: number;
  converged: boolean;
  finalMaximumChange: number;
}

interface GroupFit {
  summary: RefinementModelSummary;
  rowIndices: number[];
  localByNode: Map<number, number>;
  active: Uint8Array;
  theta: Float64Array;
}

function parseGroupKey(key: string): { scopeValue: string; locus: string; segment: "V" | "D" | "J" } {
  const [scopeValue = "", locus = "", segment = "V"] = key.split("\u0000");
  return { scopeValue, locus, segment: segment as "V" | "D" | "J" };
}

function logistic(value: number): number {
  if (value >= 0) {
    const inverse = Math.exp(-value);
    return 1 / (1 + inverse);
  }
  const exponential = Math.exp(value);
  return exponential / (1 + exponential);
}

function logSumExp(values: readonly number[]): number {
  const maximum = Math.max(...values);
  if (!Number.isFinite(maximum)) return maximum;
  let sum = 0;
  for (const value of values) sum += Math.exp(value - maximum);
  return maximum + Math.log(sum);
}

/**
 * Small deterministic quadrature for the positive-frequency slab. A gamma
 * shape below one gives the empirically required long left tail. The prior
 * mean is tied to the expected active fraction of the complete locus-matched
 * database, while the integration floor is numerical rather than a count or
 * retention threshold.
 */
function positiveFrequencyQuadrature(
  databaseNodes: number,
  options: AlleleRefinementOptions,
): Array<{ frequency: number; logWeight: number }> {
  const points = Math.max(4, Math.min(64, Math.floor(options.activeSetQuadraturePoints)));
  const floor = Math.max(1e-12, Math.min(0.05, options.activeSetFrequencyFloor));
  const ceiling = 0.95;
  const shape = Math.max(0.05, Math.min(10, options.activeSetTailShape));
  const activeFraction = Math.max(1e-6, Math.min(1 - 1e-6, options.activeSetPriorActiveFraction));
  const expectedActive = Math.max(2, databaseNodes * activeFraction);
  const mean = Math.min(0.5, 1 / expectedActive);
  const scale = mean / shape;
  const logFloor = Math.log(floor);
  const logCeiling = Math.log(ceiling);
  const raw: Array<{ frequency: number; logWeight: number }> = [];
  for (let index = 0; index < points; index += 1) {
    const lower = Math.exp(logFloor + (logCeiling - logFloor) * index / points);
    const upper = Math.exp(logFloor + (logCeiling - logFloor) * (index + 1) / points);
    const frequency = Math.sqrt(lower * upper);
    const logDensity = (shape - 1) * Math.log(frequency) - frequency / scale;
    raw.push({ frequency, logWeight: logDensity + Math.log(Math.max(Number.MIN_VALUE, upper - lower)) });
  }
  const normalization = logSumExp(raw.map((entry) => entry.logWeight));
  return raw.map((entry) => ({ ...entry, logWeight: entry.logWeight - normalization }));
}

function patternMixtures(
  patterns: readonly EvidencePattern[],
  active: Uint8Array,
  theta: Float64Array,
): Float64Array {
  const mixtures = new Float64Array(patterns.length);
  patterns.forEach((pattern, patternIndex) => {
    let mixture = 0;
    for (let offset = 0; offset < pattern.columns.length; offset += 1) {
      const column = pattern.columns[offset];
      if (active[column]) mixture += theta[column] * pattern.evidence[offset];
    }
    mixtures[patternIndex] = mixture;
  });
  return mixtures;
}

/** Ordinary sparse EM restricted to the exact current active set. */
function fitActiveFrequencies(
  patterns: readonly EvidencePattern[],
  active: Uint8Array,
  initialTheta: Float64Array,
  options: AlleleRefinementOptions,
): MixtureFit {
  const theta = Float64Array.from(initialTheta);
  let thetaSum = 0;
  for (let index = 0; index < theta.length; index += 1) {
    if (!active[index]) theta[index] = 0;
    thetaSum += theta[index];
  }
  if (!(thetaSum > 0)) {
    for (let index = 0; index < theta.length; index += 1) if (active[index]) theta[index] = 1;
    thetaSum = theta.reduce((sum, value) => sum + value, 0);
  }
  for (let index = 0; index < theta.length; index += 1) theta[index] /= thetaSum;

  const maximumIterations = Math.max(1, Math.floor(options.maxIterations));
  const tolerance = Math.max(1e-12, options.convergenceTolerance);
  let counts = new Float64Array(theta.length);
  let converged = false;
  let finalMaximumChange = Number.POSITIVE_INFINITY;
  let iterations = 0;
  for (; iterations < maximumIterations; iterations += 1) {
    counts = new Float64Array(theta.length);
    for (const pattern of patterns) {
      let normalizer = 0;
      for (let offset = 0; offset < pattern.columns.length; offset += 1) {
        const column = pattern.columns[offset];
        if (active[column]) normalizer += theta[column] * pattern.evidence[offset];
      }
      if (!(normalizer > 0)) {
        let fallback = -1;
        let fallbackEvidence = -1;
        for (let offset = 0; offset < pattern.columns.length; offset += 1) {
          const column = pattern.columns[offset];
          if (active[column] && pattern.evidence[offset] > fallbackEvidence) {
            fallback = column;
            fallbackEvidence = pattern.evidence[offset];
          }
        }
        if (fallback >= 0) counts[fallback] += pattern.weight;
        continue;
      }
      for (let offset = 0; offset < pattern.columns.length; offset += 1) {
        const column = pattern.columns[offset];
        if (!active[column]) continue;
        counts[column] += pattern.weight * theta[column] * pattern.evidence[offset] / normalizer;
      }
    }
    let total = 0;
    for (let index = 0; index < counts.length; index += 1) if (active[index]) total += counts[index];
    if (!(total > 0)) break;
    finalMaximumChange = 0;
    for (let index = 0; index < theta.length; index += 1) {
      const next = active[index] ? counts[index] / total : 0;
      finalMaximumChange = Math.max(finalMaximumChange, Math.abs(next - theta[index]));
      theta[index] = next;
    }
    if (finalMaximumChange <= tolerance) {
      converged = true;
      iterations += 1;
      break;
    }
  }
  return { theta, counts, iterations, converged, finalMaximumChange };
}

/**
 * Approximate posterior inclusion probability from an exact one-dimensional
 * conditional likelihood and the heavy-tailed positive-frequency quadrature.
 * Only patterns containing the candidate are visited; all other records enter
 * through one aggregated log(1-f) term.
 */
function inclusionProbability(
  candidate: number,
  active: Uint8Array,
  theta: Float64Array,
  mixtures: Float64Array,
  patterns: readonly EvidencePattern[],
  incident: readonly IncidentPattern[],
  totalWeight: number,
  quadrature: readonly { frequency: number; logWeight: number }[],
  priorLogOdds: number,
): number {
  const isActive = Boolean(active[candidate]);
  const candidateWeight = isActive ? theta[candidate] : 0;
  const remaining = 1 - candidateWeight;
  if (isActive && remaining <= 1e-12) return 1;
  let incidentWeight = 0;
  const ratios: Array<{ weight: number; ratio: number }> = [];
  for (const entry of incident) {
    const pattern = patterns[entry.pattern];
    incidentWeight += pattern.weight;
    const baseline = isActive
      ? (mixtures[entry.pattern] - candidateWeight * entry.evidence) / remaining
      : mixtures[entry.pattern];
    if (!(baseline > 1e-300)) return 1;
    ratios.push({ weight: pattern.weight, ratio: entry.evidence / baseline });
  }
  const outsideWeight = Math.max(0, totalWeight - incidentWeight);
  const components = quadrature.map(({ frequency, logWeight }) => {
    let delta = outsideWeight * Math.log1p(-frequency);
    for (const entry of ratios) {
      const relative = (1 - frequency) + frequency * entry.ratio;
      if (!(relative > 0)) return Number.NEGATIVE_INFINITY;
      delta += entry.weight * Math.log(relative);
    }
    return logWeight + delta;
  });
  return logistic(priorLogOdds + logSumExp(components));
}

function fitGroup(
  matrix: SparseEvidenceMatrix,
  graph: ReferenceAlleleGraph,
  options: AlleleRefinementOptions,
  key: string,
  rowIndices: number[],
): GroupFit {
  const parsed = parseGroupKey(key);
  const nodeSet = new Set<number>();
  let nonZeros = 0;
  let effectiveRows = 0;
  for (const row of rowIndices) {
    effectiveRows += matrix.weights[row];
    const begin = matrix.rowOffsets[row];
    const end = matrix.rowOffsets[row + 1];
    nonZeros += end - begin;
    for (let offset = begin; offset < end; offset += 1) nodeSet.add(matrix.columns[offset]);
  }
  const nodes = [...nodeSet].sort((left, right) => left - right);
  const localByNode = new Map(nodes.map((node, index) => [node, index] as const));
  const databaseNodes = graph.nodes.reduce((count, node) => count + Number(node.locus === parsed.locus || (!node.locus && nodeSet.has(node.index))), 0);
  const inactivePriorNodes = Math.max(0, databaseNodes - nodes.length);

  const localEvidenceCounts = new Float64Array(nodes.length);
  const patternsByKey = new Map<string, EvidencePattern>();
  for (const row of rowIndices) {
    const columns: number[] = [];
    const evidence: number[] = [];
    let patternKey = "";
    for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
      const local = localByNode.get(matrix.columns[offset]);
      if (local === undefined) continue;
      const value = Math.exp(matrix.logEvidence[offset]);
      columns.push(local);
      evidence.push(value);
      localEvidenceCounts[local] += value * matrix.weights[row];
      patternKey += `${local}:${matrix.logEvidence[offset].toFixed(6)},`;
    }
    const previous = patternsByKey.get(patternKey);
    if (previous) previous.weight += matrix.weights[row];
    else patternsByKey.set(patternKey, { columns, evidence, weight: matrix.weights[row] });
  }
  const patterns = [...patternsByKey.values()];
  const incident = Array.from({ length: nodes.length }, () => [] as IncidentPattern[]);
  patterns.forEach((pattern, patternIndex) => pattern.columns.forEach((column, offset) => {
    incident[column].push({ pattern: patternIndex, evidence: pattern.evidence[offset] });
  }));

  const active = new Uint8Array(nodes.length);
  active.fill(1);
  let theta: Float64Array = Float64Array.from(localEvidenceCounts);
  let totalInitial = theta.reduce((sum, value) => sum + value, 0);
  if (!(totalInitial > 0)) totalInitial = Math.max(1, theta.length);
  for (let index = 0; index < theta.length; index += 1) theta[index] = theta[index] > 0 ? theta[index] / totalInitial : 1 / Math.max(1, theta.length);

  const activeFraction = Math.max(1e-6, Math.min(1 - 1e-6, options.activeSetPriorActiveFraction));
  const priorLogOdds = Math.log(activeFraction) - Math.log1p(-activeFraction);
  const threshold = Math.max(0, Math.min(1, options.activeSetInclusionThreshold));
  const quadrature = positiveFrequencyQuadrature(databaseNodes, options);
  const totalWeight = patterns.reduce((sum, pattern) => sum + pattern.weight, 0);
  const inclusion = new Float64Array(nodes.length);
  inclusion.fill(1);
  let totalIterations = 0;
  let selectionConverged = false;
  let lastFit = fitActiveFrequencies(patterns, active, theta, options);
  totalIterations += lastFit.iterations;
  theta = lastFit.theta;

  const maximumSweeps = Math.max(1, Math.min(50, Math.floor(options.activeSetMaxSweeps)));
  for (let sweep = 0; sweep < maximumSweeps; sweep += 1) {
    const mixtures = patternMixtures(patterns, active, theta);
    const candidates: Array<{ index: number; probability: number }> = [];
    for (let index = 0; index < nodes.length; index += 1) {
      if (!active[index]) continue;
      const probability = inclusionProbability(index, active, theta, mixtures, patterns, incident[index], totalWeight, quadrature, priorLogOdds);
      inclusion[index] = probability;
      if (probability < threshold) candidates.push({ index, probability });
    }
    candidates.sort((left, right) => left.probability - right.probability || theta[left.index] - theta[right.index] || left.index - right.index);
    const activePerPattern = patterns.map((pattern) => pattern.columns.reduce((count, column) => count + Number(active[column]), 0));
    let removed = 0;
    for (const candidate of candidates) {
      if (!active[candidate.index]) continue;
      if (incident[candidate.index].some((entry) => activePerPattern[entry.pattern] <= 1)) continue;
      active[candidate.index] = 0;
      theta[candidate.index] = 0;
      for (const entry of incident[candidate.index]) activePerPattern[entry.pattern] -= 1;
      removed += 1;
    }
    if (!removed) {
      selectionConverged = true;
      break;
    }
    lastFit = fitActiveFrequencies(patterns, active, theta, options);
    totalIterations += lastFit.iterations;
    theta = lastFit.theta;
  }

  // Ensure the reported frequencies and active-component inclusion scores are
  // based on the final selected set, including when the sweep guard was hit.
  lastFit = fitActiveFrequencies(patterns, active, theta, options);
  totalIterations += lastFit.iterations;
  theta = lastFit.theta;
  const finalMixtures = patternMixtures(patterns, active, theta);
  for (let index = 0; index < nodes.length; index += 1) {
    if (active[index]) inclusion[index] = inclusionProbability(index, active, theta, finalMixtures, patterns, incident[index], totalWeight, quadrature, priorLogOdds);
  }

  const alleles: AllelePosteriorSummary[] = nodes.map((node, index) => {
    const mean = active[index] ? theta[index] : 0;
    return {
      nodeIndex: node,
      names: [...graph.nodes[node].names],
      sequenceLength: graph.nodes[node].sequence.length,
      posteriorMean: mean,
      expectedAssignments: active[index] ? lastFit.counts[index] : 0,
      localEvidenceAssignments: localEvidenceCounts[index],
      posteriorSd: Math.sqrt(Math.max(0, mean * (1 - mean) / Math.max(1, effectiveRows + 1))),
      inclusionProbability: inclusion[index],
      active: Boolean(active[index]),
    };
  }).sort((left, right) => Number(right.active) - Number(left.active)
    || right.posteriorMean - left.posteriorMean
    || (right.inclusionProbability ?? 0) - (left.inclusionProbability ?? 0)
    || left.names.join(",").localeCompare(right.names.join(",")));

  return {
    rowIndices,
    localByNode,
    active,
    theta,
    summary: {
      key,
      scopeValue: parsed.scopeValue,
      locus: parsed.locus,
      segment: parsed.segment,
      rows: rowIndices.length,
      effectiveRows,
      nonZeros,
      databaseNodes,
      inactivePriorNodes,
      inferenceModel: "active-set",
      activeAlleles: active.reduce((sum, value) => sum + Number(value), 0),
      alleles,
      iterations: totalIterations,
      converged: lastFit.converged && selectionConverged,
      finalMaximumChange: lastFit.finalMaximumChange,
    },
  };
}

/** Fast exact-zero hurdle alternative to the continuous Dirichlet model. */
export function fitSparseActiveSetAlleleModel(
  matrix: SparseEvidenceMatrix,
  graph: ReferenceAlleleGraph,
  options: AlleleRefinementOptions,
  totalRecords: number,
  onProgress?: (completedGroups: number, totalGroups: number) => void,
): SegmentRefinementResult {
  const mapNode = new Int32Array(totalRecords);
  const localTopNode = new Int32Array(totalRecords);
  mapNode.fill(-1);
  localTopNode.fill(-1);
  const mapProbability = new Float32Array(totalRecords);
  const localTopProbability = new Float32Array(totalRecords);
  const posteriorEntropy = new Float32Array(totalRecords);
  const modelIndex = new Int32Array(totalRecords);
  modelIndex.fill(-1);
  const assignmentWeight = new Float32Array(totalRecords);
  const rowsByGroup = matrix.groupKeys.map(() => [] as number[]);
  for (let row = 0; row < matrix.ordinals.length; row += 1) rowsByGroup[matrix.rowGroups[row]].push(row);
  const models: RefinementModelSummary[] = [];
  let changedMapRows = 0;

  for (let group = 0; group < rowsByGroup.length; group += 1) {
    const fit = fitGroup(matrix, graph, options, matrix.groupKeys[group], rowsByGroup[group]);
    models.push(fit.summary);
    for (const row of fit.rowIndices) {
      const ordinal = matrix.ordinals[row];
      modelIndex[ordinal] = group;
      assignmentWeight[ordinal] = matrix.weights[row];
      let maximum = Number.NEGATIVE_INFINITY;
      for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
        const local = fit.localByNode.get(matrix.columns[offset]);
        if (local === undefined || !fit.active[local] || !(fit.theta[local] > 0)) continue;
        maximum = Math.max(maximum, matrix.logEvidence[offset] + Math.log(fit.theta[local]));
      }
      let normalizer = 0;
      if (Number.isFinite(maximum)) {
        for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
          const local = fit.localByNode.get(matrix.columns[offset]);
          if (local === undefined || !fit.active[local] || !(fit.theta[local] > 0)) continue;
          normalizer += Math.exp(matrix.logEvidence[offset] + Math.log(fit.theta[local]) - maximum);
        }
      }
      let topNode = matrix.localTop[row];
      let topProbability = 1;
      let entropy = 0;
      if (normalizer > 0) {
        topNode = -1;
        topProbability = -1;
        for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
          const node = matrix.columns[offset];
          const local = fit.localByNode.get(node);
          if (local === undefined || !fit.active[local] || !(fit.theta[local] > 0)) continue;
          const probability = Math.exp(matrix.logEvidence[offset] + Math.log(fit.theta[local]) - maximum) / normalizer;
          if (probability > topProbability || (probability === topProbability && node < topNode)) {
            topNode = node;
            topProbability = probability;
          }
          if (probability > 0) entropy -= probability * Math.log(probability);
        }
      }
      mapNode[ordinal] = topNode;
      mapProbability[ordinal] = Math.max(0, topProbability);
      posteriorEntropy[ordinal] = entropy;
      localTopNode[ordinal] = matrix.localTop[row];
      localTopProbability[ordinal] = matrix.localTopProbability[row];
      if (topNode !== matrix.localTop[row]) changedMapRows += 1;
    }
    onProgress?.(group + 1, rowsByGroup.length);
  }

  return {
    segment: graph.segment,
    nodes: graph.nodes,
    mapNode,
    mapProbability,
    posteriorEntropy,
    localTopNode,
    localTopProbability,
    modelIndex,
    assignmentWeight,
    models,
    modeledRows: matrix.ordinals.length,
    changedMapRows,
    skippedRows: matrix.skippedRows,
    matrixNonZeros: matrix.columns.length,
    truncatedRows: matrix.truncatedRows,
    exactDuplicateLabels: graph.exactDuplicateLabels,
  };
}
