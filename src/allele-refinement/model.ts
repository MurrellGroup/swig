import type {
  AllelePosteriorSummary,
  AlleleRefinementOptions,
  ReferenceAlleleGraph,
  RefinementModelSummary,
  SegmentRefinementResult,
  SparseEvidenceMatrix,
} from "./types.ts";

/** Accurate enough for variational mixture updates over positive Dirichlet parameters. */
export function digamma(value: number): number {
  let x = value;
  let result = 0;
  while (x < 8) {
    result -= 1 / x;
    x += 1;
  }
  const inverse = 1 / x;
  const square = inverse * inverse;
  return result + Math.log(x) - 0.5 * inverse - square * (1 / 12 - square * (1 / 120 - square / 252));
}

interface GroupFit {
  summary: RefinementModelSummary;
  rowIndices: number[];
  nodes: number[];
  expectedLog: Float64Array;
}

function parseGroupKey(key: string): { scopeValue: string; locus: string; segment: "V" | "D" | "J" } {
  const [scopeValue = "", locus = "", segment = "V"] = key.split("\u0000");
  return { scopeValue, locus, segment: segment as "V" | "D" | "J" };
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
  // Prior-only database nodes remain implicit so a large reference set is not
  // materialized once per donor. Their alpha mass still enters the complete
  // Dirichlet normalization. Unlabelled custom nodes enter only through an
  // observed sparse candidate row because their receptor locus is unknown.
  const databaseNodes = graph.nodes.reduce((count, node) => count + Number(node.locus === parsed.locus || (!node.locus && nodeSet.has(node.index))), 0);
  const inactivePriorNodes = Math.max(0, databaseNodes - nodes.length);
  const local = new Map(nodes.map((node, index) => [node, index] as const));
  const alpha = Math.max(1e-9, options.alphaPerAllele);
  const gamma = new Float64Array(nodes.length);
  const localEvidenceCounts = new Float64Array(nodes.length);
  gamma.fill(alpha);
  for (const row of rowIndices) {
    const begin = matrix.rowOffsets[row];
    const end = matrix.rowOffsets[row + 1];
    const weight = matrix.weights[row];
    for (let offset = begin; offset < end; offset += 1) {
      const index = local.get(matrix.columns[offset]);
      if (index === undefined) continue;
      const value = Math.exp(matrix.logEvidence[offset]) * weight;
      gamma[index] += value;
      localEvidenceCounts[index] += value;
    }
  }
  const patterns = new Map<string, { columns: number[]; logEvidence: number[]; weight: number }>();
  for (const row of rowIndices) {
    const begin = matrix.rowOffsets[row];
    const end = matrix.rowOffsets[row + 1];
    const columns: number[] = [];
    const evidence: number[] = [];
    let key = "";
    for (let offset = begin; offset < end; offset += 1) {
      const column = matrix.columns[offset];
      const logValue = matrix.logEvidence[offset];
      columns.push(column);
      evidence.push(logValue);
      key += `${column}:${logValue.toFixed(6)},`;
    }
    const previous = patterns.get(key);
    if (previous) previous.weight += matrix.weights[row];
    else patterns.set(key, { columns, logEvidence: evidence, weight: matrix.weights[row] });
  }
  let converged = false;
  let finalMaximumChange = Number.POSITIVE_INFINITY;
  let iterations = 0;
  for (; iterations < options.maxIterations; iterations += 1) {
    let gammaSum = inactivePriorNodes * alpha;
    for (const value of gamma) gammaSum += value;
    const expectedLog = Float64Array.from(gamma, (value) => digamma(value) - digamma(gammaSum));
    const counts = new Float64Array(nodes.length);
    for (const pattern of patterns.values()) {
      let maximum = Number.NEGATIVE_INFINITY;
      for (let offset = 0; offset < pattern.columns.length; offset += 1) {
        const index = local.get(pattern.columns[offset])!;
        maximum = Math.max(maximum, pattern.logEvidence[offset] + expectedLog[index]);
      }
      let normalizer = 0;
      for (let offset = 0; offset < pattern.columns.length; offset += 1) normalizer += Math.exp(pattern.logEvidence[offset] + expectedLog[local.get(pattern.columns[offset])!] - maximum);
      for (let offset = 0; offset < pattern.columns.length; offset += 1) {
        const probability = Math.exp(pattern.logEvidence[offset] + expectedLog[local.get(pattern.columns[offset])!] - maximum) / normalizer;
        counts[local.get(pattern.columns[offset])!] += pattern.weight * probability;
      }
    }
    finalMaximumChange = 0;
    for (let index = 0; index < gamma.length; index += 1) {
      const next = alpha + counts[index];
      finalMaximumChange = Math.max(finalMaximumChange, Math.abs(next - gamma[index]) / Math.max(1, gamma[index]));
      gamma[index] = next;
    }
    if (finalMaximumChange <= options.convergenceTolerance) {
      converged = true;
      iterations += 1;
      break;
    }
  }
  let gammaSum = inactivePriorNodes * alpha;
  for (const value of gamma) gammaSum += value;
  const expectedLog = Float64Array.from(gamma, (value) => digamma(value) - digamma(gammaSum));
  const alleles: AllelePosteriorSummary[] = nodes.map((node, index) => {
    const mean = gamma[index] / gammaSum;
    const variance = gamma[index] * (gammaSum - gamma[index]) / (gammaSum * gammaSum * (gammaSum + 1));
    return {
      nodeIndex: node,
      names: [...graph.nodes[node].names],
      sequenceLength: graph.nodes[node].sequence.length,
      posteriorMean: mean,
      expectedAssignments: Math.max(0, gamma[index] - alpha),
      localEvidenceAssignments: localEvidenceCounts[index],
      posteriorSd: Math.sqrt(Math.max(0, variance)),
    };
  }).sort((left, right) => right.posteriorMean - left.posteriorMean || left.names.join(",").localeCompare(right.names.join(",")));
  return {
    rowIndices,
    nodes,
    expectedLog,
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
      alleles,
      iterations,
      converged,
      finalMaximumChange,
    },
  };
}

export function fitSparseAlleleModel(
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
    const local = new Map(fit.nodes.map((node, index) => [node, index] as const));
    fit.rowIndices.forEach((row) => {
      const ordinal = matrix.ordinals[row];
      modelIndex[ordinal] = group;
      assignmentWeight[ordinal] = matrix.weights[row];
      const begin = matrix.rowOffsets[row];
      let topNode = -1;
      let topProbability = -1;
      let entropy = 0;
      let maximum = Number.NEGATIVE_INFINITY;
      for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) maximum = Math.max(maximum, matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])!]);
      let normalizer = 0;
      for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) normalizer += Math.exp(matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])!] - maximum);
      for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) {
        const probability = Math.exp(matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])!] - maximum) / normalizer;
        const node = matrix.columns[offset];
        if (probability > topProbability || (probability === topProbability && node < topNode)) {
          topNode = node;
          topProbability = probability;
        }
        if (probability > 0) entropy -= probability * Math.log(probability);
      }
      mapNode[ordinal] = topNode;
      mapProbability[ordinal] = Math.max(0, topProbability);
      posteriorEntropy[ordinal] = entropy;
      localTopNode[ordinal] = matrix.localTop[row];
      localTopProbability[ordinal] = matrix.localTopProbability[row];
      if (topNode !== matrix.localTop[row]) changedMapRows += 1;
    });
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
