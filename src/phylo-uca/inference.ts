import { parseFasta } from "../post-analysis-core.ts";
import { phyloUcaHmmLogMarginal, phyloUcaHmmPosterior } from "./hmm.ts";
import { preparePhyloUcaReferences } from "./references.ts";
import { PhyloUcaTreeMessages } from "./tree-messages.ts";
import type {
  PhyloUcaInput,
  PhyloUcaPlacement,
  PhyloUcaProgress,
  PhyloUcaResult,
  PhyloUcaSitePosterior,
} from "./types.ts";

const NEGATIVE_INFINITY = Number.NEGATIVE_INFINITY;
const CHARACTERS = ["A", "C", "G", "T", "-"] as const;

function entropyBits(probabilities: readonly number[]): number {
  let entropy = 0;
  for (const probability of probabilities) if (probability > 0) entropy -= probability * Math.log2(probability);
  return entropy;
}

function normalizedWeights(scores: readonly number[]): number[] {
  if (!scores.length) return [];
  const maximum = Math.max(...scores);
  const weights = scores.map((score) => Math.exp(score - maximum));
  const total = weights.reduce((sum, weight) => sum + weight, 0);
  return weights.map((weight) => weight / total);
}

function uniqueNumbers(values: readonly number[], tolerance = 1e-10): number[] {
  const sorted = [...values].filter(Number.isFinite).sort((left, right) => left - right);
  return sorted.filter((value, index) => index === 0 || Math.abs(value - sorted[index - 1]) > tolerance);
}

function linearGrid(points: number, maximum: number): number[] {
  if (points <= 1 || !(maximum > 0)) return [0];
  return Array.from({ length: points }, (_, index) => maximum * index / (points - 1));
}

function branchGrid(points: number, maximum: number): number[] {
  if (points <= 1 || !(maximum > 0)) return [0];
  return Array.from({ length: points }, (_, index) => maximum * Math.pow(index / (points - 1), 2));
}

function placementKey(edge: number, distance: number, branch: number): string {
  return `${edge}|${distance.toPrecision(11)}|${branch.toPrecision(11)}`;
}

export async function inferPhyloUca(input: PhyloUcaInput, onProgress?: (progress: PhyloUcaProgress) => void): Promise<PhyloUcaResult> {
  const started = performance.now();
  const progress = (phase: PhyloUcaProgress["phase"], processed: number, total: number, detail: string) => onProgress?.({ phase, processed, total, detail });
  progress("references", 0, 1, "Preparing broad V, D, and J hypothesis sets");
  const observedRecords = parseFasta(input.observedAlignmentFasta, true);
  if (!observedRecords.length) throw new Error("The observed-only alignment is empty.");
  const observedColumns = observedRecords[0].sequence.length;
  if (!observedColumns || observedRecords.some((record) => record.sequence.length !== observedColumns)) throw new Error("The observed-only alignment is not rectangular.");
  if (input.retainedColumns.length !== observedColumns) throw new Error("The retained-column map does not match the observed-only alignment.");
  const curated = parseFasta(input.curatedAlignmentFasta, true);
  const originalColumns = curated[0]?.sequence.length ?? observedColumns;
  const guideRecord = curated.find((record) => record.name === input.germlineGuideName);
  const projectedGuide = input.retainedColumns.map((column) => guideRecord?.sequence[column]?.toUpperCase().replace("U", "T").replace(".", "-") ?? "N").join("");
  const references = preparePhyloUcaReferences(projectedGuide, input.lineageRows, input.references, input.locus, input.options.candidates);
  progress("references", 1, 1, `${references.v.length} V, ${references.d.length} D, and ${references.j.length} J candidates retained`);

  let messageCount = 0;
  progress("tree-messages", 0, 1, "Computing directed Felsenstein messages on the fixed observed tree");
  const tree = new PhyloUcaTreeMessages(
    input.observedAlignmentFasta,
    input.observedTreeNewick,
    input.options.model,
    input.options.characterMode,
    (complete, total) => {
      messageCount = complete;
      progress("tree-messages", complete, total, "Computing all directed half-edge messages");
    },
  );
  if (!messageCount) progress("tree-messages", 1, 1, "Directed messages ready");
  const warnings = [...references.warnings];
  const observedHasGap = observedRecords.some((record) => record.sequence.includes("-") || record.sequence.includes("."));
  if (input.options.characterMode === "nucleotide-gtr4" && observedHasGap) warnings.push("The alignment contains gaps but the advanced GTR4 override was selected; observed gaps were treated as missing characters. Auto mode would use the gap-aware GTR5 approximation.");
  warnings.push(tree.characterModel === "gap-aware-gtr5"
    ? "The observed alignment contains gaps, so the phylogenetic likelihood used the explicit A/C/G/T/gap fixed-alignment model. This is not a continuous-time insertion/deletion process."
    : "The observed alignment contains no gaps after removal of all-gap columns, so the phylogenetic likelihood used ordinary four-state nucleotide GTR.");

  const guideBranchSamples = [0, Math.min(0.02, input.options.search.maximumUcaBranchLength), Math.min(0.08, input.options.search.maximumUcaBranchLength)];
  const edgeScores: Array<{ edge: number; guideScore: number }> = [];
  progress("edge-screen", 0, tree.edges.length, "Screening every edge with the germline guide");
  for (let edgeIndex = 0; edgeIndex < tree.edges.length; edgeIndex += 1) {
    const edge = tree.edges[edgeIndex];
    let guideScore = NEGATIVE_INFINITY;
    for (const branch of guideBranchSamples) {
      const surface = tree.conditionalLikelihoods(edgeIndex, edge.length / 2, branch);
      guideScore = Math.max(guideScore, tree.guideScore(surface, references.guide));
    }
    edgeScores.push({ edge: edgeIndex, guideScore });
    progress("edge-screen", edgeIndex + 1, tree.edges.length, `Screened ${edge.id}`);
    if ((edgeIndex & 3) === 3) await Promise.resolve();
  }
  edgeScores.sort((left, right) => right.guideScore - left.guideScore);
  const selectedEdges = edgeScores.slice(0, Math.max(1, Math.min(tree.edges.length, Math.floor(input.options.search.fullHmmEdges))));
  const evaluated = new Map<string, PhyloUcaPlacement>();
  const totalEdgeLength = tree.edges.reduce((sum, edge) => sum + Math.max(1e-12, edge.length), 0);
  const branchMean = Math.max(1e-5, input.options.search.branchPriorMean);
  const evaluate = (edgeIndex: number, distance: number, branch: number, guideScore: number): PhyloUcaPlacement => {
    const edge = tree.edges[edgeIndex];
    const safeDistance = Math.max(0, Math.min(edge.length, distance));
    const safeBranch = Math.max(0, Math.min(input.options.search.maximumUcaBranchLength, branch));
    const key = placementKey(edgeIndex, safeDistance, safeBranch);
    const previous = evaluated.get(key);
    if (previous) return previous;
    const surface = tree.conditionalLikelihoods(edgeIndex, safeDistance, safeBranch);
    const likelihood = phyloUcaHmmLogMarginal(surface, references, input.options.hmm);
    const edgePrior = input.options.search.edgePrior === "uniform-length"
      ? Math.log(Math.max(1e-12, edge.length) / Math.max(1e-12, totalEdgeLength))
      : -Math.log(tree.edges.length);
    const branchPrior = -safeBranch / branchMean - Math.log(branchMean);
    const placement: PhyloUcaPlacement = {
      edgeId: edge.id,
      endpointA: edge.endpointA,
      endpointB: edge.endpointB,
      distanceFromA: safeDistance,
      edgeLength: edge.length,
      edgeFraction: edge.length > 0 ? safeDistance / edge.length : 0,
      ucaBranchLength: safeBranch,
      logMarginalLikelihood: likelihood,
      logPosteriorScore: likelihood + edgePrior + branchPrior,
      localPosteriorWeight: 0,
      guideScore,
    };
    evaluated.set(key, placement);
    return placement;
  };

  const coarseTasks = selectedEdges.flatMap(({ edge, guideScore }) => {
    const treeEdge = tree.edges[edge];
    return linearGrid(Math.max(2, input.options.search.edgeGridPoints), treeEdge.length).flatMap((distance) =>
      branchGrid(Math.max(2, input.options.search.branchGridPoints), input.options.search.maximumUcaBranchLength).map((branch) => ({ edge, distance, branch, guideScore }))
    );
  });
  progress("hmm-search", 0, coarseTasks.length, "Evaluating the recombination HMM on guide-ranked edges");
  for (let index = 0; index < coarseTasks.length; index += 1) {
    const task = coarseTasks[index];
    evaluate(task.edge, task.distance, task.branch, task.guideScore);
    progress("hmm-search", index + 1, coarseTasks.length, `Full HMM placement ${index + 1} of ${coarseTasks.length}`);
    if ((index & 1) === 1) await Promise.resolve();
  }
  for (let round = 0; round < Math.max(0, input.options.search.localRefinementRounds); round += 1) {
    const leaders = [...evaluated.values()].sort((left, right) => right.logPosteriorScore - left.logPosteriorScore).slice(0, 2);
    const tasks: Array<{ edge: number; distance: number; branch: number; guideScore: number }> = [];
    for (const leader of leaders) {
      const edge = tree.edges.find((candidate) => candidate.id === leader.edgeId)!;
      const edgeStep = edge.length / Math.max(2, input.options.search.edgeGridPoints - 1) / Math.pow(2, round + 1);
      const branchStep = input.options.search.maximumUcaBranchLength / Math.max(2, input.options.search.branchGridPoints - 1) / Math.pow(2, round + 1);
      const edgeIndex = edge.index;
      for (const distance of uniqueNumbers([leader.distanceFromA - edgeStep, leader.distanceFromA, leader.distanceFromA + edgeStep])) {
        for (const branch of uniqueNumbers([leader.ucaBranchLength - branchStep, leader.ucaBranchLength, leader.ucaBranchLength + branchStep])) {
          if (distance >= 0 && distance <= edge.length && branch >= 0 && branch <= input.options.search.maximumUcaBranchLength) tasks.push({ edge: edgeIndex, distance, branch, guideScore: leader.guideScore });
        }
      }
    }
    for (let index = 0; index < tasks.length; index += 1) {
      const task = tasks[index];
      evaluate(task.edge, task.distance, task.branch, task.guideScore);
      progress("hmm-search", index + 1, tasks.length, `Local placement refinement ${round + 1}`);
      if ((index & 1) === 1) await Promise.resolve();
    }
  }
  const placements = [...evaluated.values()].sort((left, right) => right.logPosteriorScore - left.logPosteriorScore);
  if (!placements.length) throw new Error("No UCA attachment placement could be evaluated.");
  const localCount = input.options.search.marginalizeLocally ? Math.max(1, Math.min(placements.length, input.options.search.localPosteriorPoints)) : 1;
  const local = placements.slice(0, localCount);
  const weights = normalizedWeights(local.map((placement) => placement.logPosteriorScore));
  local.forEach((placement, index) => { placement.localPosteriorWeight = weights[index]; });
  const mixture = Array.from({ length: observedColumns }, () => [0, 0, 0, 0, 0] as [number, number, number, number, number]);
  let bestPosterior: ReturnType<typeof phyloUcaHmmPosterior> | null = null;
  progress("posterior", 0, local.length, "Marginalizing UCA nucleotides around the best placement");
  for (let index = 0; index < local.length; index += 1) {
    const placement = local[index];
    const edge = tree.edges.find((candidate) => candidate.id === placement.edgeId)!;
    const surface = tree.conditionalLikelihoods(edge.index, placement.distanceFromA, placement.ucaBranchLength);
    const posterior = phyloUcaHmmPosterior(surface, references, input.options.hmm);
    if (index === 0) bestPosterior = posterior;
    for (let site = 0; site < observedColumns; site += 1) for (let character = 0; character < 5; character += 1) mixture[site][character] += weights[index] * posterior.probabilities[site][character];
    progress("posterior", index + 1, local.length, `Integrated placement ${index + 1} of ${local.length}`);
    await Promise.resolve();
  }
  if (!bestPosterior) throw new Error("The best UCA placement did not produce a posterior sequence.");
  let mapAlignedSequence = "-".repeat(originalColumns);
  let posteriorConsensusAligned = "-".repeat(originalColumns);
  const mapCharacters = [...mapAlignedSequence];
  const consensusCharacters = [...posteriorConsensusAligned];
  const posteriorByOriginal = new Map<number, PhyloUcaSitePosterior>();
  for (let site = 0; site < observedColumns; site += 1) {
    const original = input.retainedColumns[site];
    mapCharacters[original] = bestPosterior.mapAlignedSequence[site];
    let mapCharacter = 0;
    for (let character = 1; character < 5; character += 1) if (mixture[site][character] > mixture[site][mapCharacter]) mapCharacter = character;
    consensusCharacters[original] = CHARACTERS[mapCharacter];
    const probabilities = mixture[site];
    posteriorByOriginal.set(original, {
      alignmentColumn: original + 1,
      probabilities,
      mapCharacter: CHARACTERS[mapCharacter],
      mapProbability: probabilities[mapCharacter],
      entropyBits: entropyBits(probabilities),
      segment: bestPosterior.stateKinds[site],
      call: bestPosterior.stateCalls[site],
    });
  }
  mapAlignedSequence = mapCharacters.join("");
  posteriorConsensusAligned = consensusCharacters.join("");
  const posterior: PhyloUcaSitePosterior[] = [];
  for (let column = 0; column < originalColumns; column += 1) posterior.push(posteriorByOriginal.get(column) ?? {
    alignmentColumn: column + 1,
    probabilities: [0, 0, 0, 0, 1],
    mapCharacter: "-",
    mapProbability: 1,
    entropyBits: 0,
    segment: "unknown",
  });
  const path = bestPosterior.path.map((segment) => ({
    ...segment,
    startColumn: (input.retainedColumns[segment.startColumn] ?? segment.startColumn) + 1,
    endColumn: (input.retainedColumns[segment.endColumn] ?? segment.endColumn) + 1,
  }));
  const bestPlacement = placements[0];
  const bestEdge = tree.edges.find((edge) => edge.id === bestPlacement.edgeId)!;
  const effectivePlacementCount = Math.exp(-weights.reduce((sum, weight) => weight > 0 ? sum + weight * Math.log(weight) : sum, 0));
  progress("finalize", 1, 1, "Preparing UCA sequence, placement tree, and provenance");
  return {
    schema: 1,
    method: "fixed-tree-empirical-bayes-phylo-uca",
    lineageLabel: input.lineageLabel,
    generatedAt: new Date().toISOString(),
    elapsedMs: performance.now() - started,
    characterModel: tree.characterModel,
    model: input.options.model,
    options: input.options,
    observedTreeNewick: input.observedTreeNewick,
    placedTreeNewick: tree.placedTreeNewick(bestEdge.index, bestPlacement.distanceFromA, bestPlacement.ucaBranchLength),
    observedAlignmentFasta: input.observedAlignmentFasta,
    retainedColumns: [...input.retainedColumns],
    alignmentFingerprint: input.alignmentFingerprint,
    bestPlacement,
    placements,
    mapAlignedSequence,
    mapUngappedSequence: mapAlignedSequence.replaceAll("-", ""),
    posteriorConsensusAligned,
    posterior,
    path,
    candidateReport: references.report,
    mapVCall: bestPosterior.mapVCall,
    mapDCalls: bestPosterior.mapDCalls,
    mapJCall: bestPosterior.mapJCall,
    logMarginalLikelihood: bestPosterior.logMarginalLikelihood,
    effectivePlacementCount,
    warnings,
  };
}
