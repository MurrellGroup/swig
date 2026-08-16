import { parseFasta } from "../post-analysis-core.ts";
import {
  PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT,
  PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT,
  PhyloUcaHmmGibbsSampler,
  phyloUcaHmmPosterior,
} from "./hmm.ts";
import { preparePhyloUcaReferences } from "./references.ts";
import { PhyloUcaTreeMessages } from "./tree-messages.ts";
import { normalizeProbabilityVector } from "../probability-logo.ts";
import { PHYLO_UCA_CODON_STATE_COUNT, PHYLO_UCA_CODON_SYMBOLS } from "./codons.ts";
import { gridCellWidths, phyloUcaBranchLengthGrid } from "./search-grid.ts";
import type {
  PhyloUcaCodonPosterior,
  PhyloUcaDCountPosteriorPoint,
  PhyloUcaHmmAnnotationTrack,
  PhyloUcaInput,
  PhyloUcaPlacement,
  PhyloUcaProgress,
  PhyloUcaResult,
  PhyloUcaMcmcDiagnostics,
  PhyloUcaSitePosterior,
} from "./types.ts";

const NEGATIVE_INFINITY = Number.NEGATIVE_INFINITY;
const CHARACTERS = ["A", "C", "G", "T", "-"] as const;

interface AnnotationTrackMixture {
  metadata: Omit<PhyloUcaHmmAnnotationTrack, "points" | "maximumWeight">;
  points: Map<number, Float64Array>;
}

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

function addAnnotationTrackMixture(
  mixture: Map<string, AnnotationTrackMixture>,
  tracks: readonly PhyloUcaHmmAnnotationTrack[],
  placementWeight: number,
): void {
  for (const track of tracks) {
    let destination = mixture.get(track.id);
    if (!destination) {
      const { points: _points, maximumWeight: _maximumWeight, ...metadata } = track;
      destination = { metadata, points: new Map() };
      mixture.set(track.id, destination);
    }
    for (const point of track.points) {
      let masses = destination.points.get(point.alignmentColumn);
      if (!masses) {
        masses = new Float64Array(5);
        destination.points.set(point.alignmentColumn, masses);
      }
      for (let character = 0; character < 5; character += 1) masses[character] += placementWeight * point.probabilities[character];
    }
  }
}

function finalizeAnnotationTrackMixture(mixture: Map<string, AnnotationTrackMixture>): PhyloUcaHmmAnnotationTrack[] {
  return [...mixture.values()].map(({ metadata, points: rawPoints }) => {
    let maximumWeight = 0;
    const points = [...rawPoints.entries()].sort(([left], [right]) => left - right).map(([alignmentColumn, raw]) => {
      const probabilities = Array.from(raw) as [number, number, number, number, number];
      maximumWeight = Math.max(maximumWeight, probabilities.reduce((sum, value) => sum + value, 0));
      return { alignmentColumn, probabilities };
    });
    return { ...metadata, points, maximumWeight };
  });
}

function projectAnnotationTracks(tracks: readonly PhyloUcaHmmAnnotationTrack[], retainedColumns: readonly number[]): PhyloUcaHmmAnnotationTrack[] {
  return tracks.map((track) => ({
    ...track,
    points: track.points.map((point) => ({
      ...point,
      alignmentColumn: (retainedColumns[point.alignmentColumn - 1] ?? point.alignmentColumn - 1) + 1,
    })),
  }));
}

function annotationAlleleGroup(track: PhyloUcaHmmAnnotationTrack): string {
  return track.call ? `${track.kind}|${track.dOrdinal ?? 0}|${track.call}` : track.id;
}

function pruneMarginalAnnotationTracks(tracks: readonly PhyloUcaHmmAnnotationTrack[]): {
  retained: PhyloUcaHmmAnnotationTrack[];
  omittedCount: number;
  omittedMaximumWeight: number;
} {
  const groups = new Map<string, PhyloUcaHmmAnnotationTrack[]>();
  for (const track of tracks) {
    const key = annotationAlleleGroup(track);
    const group = groups.get(key) ?? [];
    group.push(track);
    groups.set(key, group);
  }
  const retainedIds = new Set<string>();
  for (const group of groups.values()) {
    if (!group[0].call) {
      for (const track of group) if (track.maximumWeight >= PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT) retainedIds.add(track.id);
      continue;
    }
    const columnTotals = new Map<number, number>();
    for (const track of group) for (const point of track.points) {
      const total = point.probabilities.reduce((sum, value) => sum + value, 0);
      columnTotals.set(point.alignmentColumn, (columnTotals.get(point.alignmentColumn) ?? 0) + total);
    }
    let groupMaximum = 0;
    for (const total of columnTotals.values()) groupMaximum = Math.max(groupMaximum, total);
    if (groupMaximum < PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT) continue;
    const eligible = group.filter((track) => track.maximumWeight >= PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT);
    for (const track of eligible.length ? eligible : [...group].sort((left, right) => right.maximumWeight - left.maximumWeight).slice(0, 1)) retainedIds.add(track.id);
  }
  const retained = tracks.filter((track) => retainedIds.has(track.id));
  const omitted = tracks.filter((track) => !retainedIds.has(track.id));
  let omittedMaximumWeight = 0;
  for (const track of omitted) omittedMaximumWeight = Math.max(omittedMaximumWeight, track.maximumWeight);
  return {
    retained,
    omittedCount: omitted.length,
    omittedMaximumWeight,
  };
}

function uniqueNumbers(values: readonly number[], tolerance = 1e-10): number[] {
  const sorted = [...values].filter(Number.isFinite).sort((left, right) => left - right);
  return sorted.filter((value, index) => index === 0 || Math.abs(value - sorted[index - 1]) > tolerance);
}

function linearGrid(points: number, maximum: number): number[] {
  if (points <= 1 || !(maximum > 0)) return [0];
  return Array.from({ length: points }, (_, index) => maximum * index / (points - 1));
}

function placementKey(edge: number, distance: number, branch: number): string {
  return `${edge}|${distance.toPrecision(11)}|${branch.toPrecision(11)}`;
}

function maximizeUnitInterval(
  objective: (value: number) => number,
  tolerance: number,
  maximumIterations = 24,
): { value: number; score: number } {
  const cache = new Map<string, number>();
  const evaluate = (value: number) => {
    const bounded = Math.max(0, Math.min(1, value));
    const key = bounded.toPrecision(14);
    const previous = cache.get(key);
    if (previous !== undefined) return previous;
    const score = objective(bounded);
    cache.set(key, score);
    return score;
  };
  const inversePhi = (Math.sqrt(5) - 1) / 2;
  let left = 0;
  let right = 1;
  let middleLeft = right - inversePhi * (right - left);
  let middleRight = left + inversePhi * (right - left);
  let leftScore = evaluate(middleLeft);
  let rightScore = evaluate(middleRight);
  evaluate(0);
  evaluate(1);
  for (let iteration = 0; iteration < maximumIterations && right - left > Math.max(1e-5, tolerance); iteration += 1) {
    if (leftScore >= rightScore) {
      right = middleRight;
      middleRight = middleLeft;
      rightScore = leftScore;
      middleLeft = right - inversePhi * (right - left);
      leftScore = evaluate(middleLeft);
    } else {
      left = middleLeft;
      middleLeft = middleRight;
      leftScore = rightScore;
      middleRight = left + inversePhi * (right - left);
      rightScore = evaluate(middleRight);
    }
  }
  let bestValue = 0;
  let bestScore = Number.NEGATIVE_INFINITY;
  for (const [raw, score] of cache) if (score > bestScore) {
    bestValue = Number(raw);
    bestScore = score;
  }
  return { value: bestValue, score: bestScore };
}

function seededRandom(rawSeed: number): () => number {
  let state = (Number.isFinite(rawSeed) ? Math.floor(rawSeed) : 1729) >>> 0;
  return () => {
    state += 0x6d2b79f5;
    let value = state;
    value = Math.imul(value ^ value >>> 15, value | 1);
    value ^= value + Math.imul(value ^ value >>> 7, value | 61);
    return ((value ^ value >>> 14) >>> 0) / 4294967296;
  };
}

function reflected(value: number, maximum: number): number {
  if (!(maximum > 0)) return 0;
  const period = 2 * maximum;
  const wrapped = ((value % period) + period) % period;
  return wrapped <= maximum ? wrapped : period - wrapped;
}

function effectiveSampleSize(values: readonly number[]): number {
  const count = values.length;
  if (count < 2) return count;
  const mean = values.reduce((sum, value) => sum + value, 0) / count;
  const centered = values.map((value) => value - mean);
  const variance = centered.reduce((sum, value) => sum + value * value, 0) / count;
  if (!(variance > 0)) return count;
  let correlationSum = 0;
  for (let lag = 1; lag < Math.min(count, 200); lag += 1) {
    let covariance = 0;
    for (let index = 0; index + lag < count; index += 1) covariance += centered[index] * centered[index + lag];
    const correlation = covariance / ((count - lag) * variance);
    if (!(correlation > 0)) break;
    correlationSum += correlation;
  }
  return Math.max(1, Math.min(count, count / (1 + 2 * correlationSum)));
}

function fixedUcaLogLikelihood(surface: ReturnType<PhyloUcaTreeMessages["conditionalLikelihoods"]>, characters: Int8Array): number {
  let total = 0;
  for (let site = 0; site < surface.sites; site += 1) {
    const character = characters[site];
    if (character < 0 || character >= surface.stateCount) continue;
    total += surface.logLikelihoods[site * surface.stateCount + character];
  }
  return total;
}

/** Equal-weight, independent nucleotide mixtures over retained V and J alleles. */
export function vjNucleotideMixtureProfile(references: Pick<ReturnType<typeof preparePhyloUcaReferences>, "v" | "j" | "vEndColumn" | "jStartColumn" | "guide">): Array<[number, number, number, number, number] | null> {
  const profile: Array<[number, number, number, number, number] | null> = Array.from({ length: references.guide.length }, () => null);
  const add = (candidates: typeof references.v, start: number, end: number) => {
    for (let column = Math.max(0, start); column <= Math.min(profile.length - 1, end); column += 1) {
      const masses: [number, number, number, number, number] = [0, 0, 0, 0, 0];
      for (const candidate of candidates) {
        const state = ["A", "C", "G", "T", "-"].indexOf(candidate.projection[column] ?? "N");
        if (state >= 0) masses[state] += 1;
      }
      if (masses.some((value) => value > 0)) profile[column] = masses;
    }
  };
  add(references.v, 0, references.vEndColumn);
  add(references.j, references.jStartColumn, profile.length - 1);
  return profile;
}

export async function inferPhyloUca(input: PhyloUcaInput, onProgress?: (progress: PhyloUcaProgress) => void): Promise<PhyloUcaResult> {
  const started = performance.now();
  const progress = (phase: PhyloUcaProgress["phase"], processed: number, total: number, detail: string) => onProgress?.({ phase, processed, total, detail });
  const frameOffset = input.frameOffset === 1 || input.frameOffset === 2 ? input.frameOffset : 0;
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
  if (tree.terminalMissingGapCount > 0) warnings.push(`${tree.terminalMissingGapCount.toLocaleString()} leading/trailing tip-gap characters were treated as missing sequence coverage, not as a fifth character.`);
  if (input.options.characterMode === "nucleotide-gtr4" && tree.internalGapCount > 0) warnings.push(`The alignment contains ${tree.internalGapCount.toLocaleString()} internal tip-gap characters but the advanced GTR4 override was selected; those internal gaps were treated as missing. Auto mode would use gap-aware GTR5.`);
  warnings.push(tree.characterModel === "gap-aware-gtr5"
    ? `The observed alignment contains ${tree.internalGapCount.toLocaleString()} internal tip-gap characters, so those internal positions used the explicit A/C/G/T/gap fixed-alignment model. This is not a continuous-time insertion/deletion process.`
    : tree.internalGapCount > 0
      ? "The forced four-state model treated internal as well as terminal gaps as missing data."
      : "The observed alignment contains no internal gaps, so the phylogenetic likelihood used ordinary four-state nucleotide GTR; terminal tip gaps, including columns missing at every tip, were missing data.");

  const search = input.options.search;
  const inferenceMode = search.inferenceMode ?? (search.marginalizeLocally ? "grid-marginalization" : "maximum-likelihood");
  const screenMode = search.screenMode ?? "vj-mixture";
  const maximumUcaBranchLength = Math.max(0, Number.isFinite(search.maximumUcaBranchLength) ? search.maximumUcaBranchLength : 0.3);
  const branchFromUnit = (value: number) => maximumUcaBranchLength * Math.pow(Math.max(0, Math.min(1, value)), 4);
  const screenEdgePoints = Math.max(3, Math.floor(search.screenEdgeGridPoints ?? 5));
  const guideBranchSamples = uniqueNumbers([
    0,
    Math.min(maximumUcaBranchLength, Math.max(1e-5, search.minimumPositiveUcaBranchLength ?? 1e-5)),
    Math.min(maximumUcaBranchLength, 0.001),
    Math.min(maximumUcaBranchLength, 0.004),
    Math.min(maximumUcaBranchLength, 0.01),
    Math.min(maximumUcaBranchLength, 0.02),
    Math.min(maximumUcaBranchLength, 0.08),
  ]);
  const vjProfile = vjNucleotideMixtureProfile(references);
  const screenScore = (edge: number, distance: number, branch: number) => {
    const surface = tree.conditionalLikelihoods(edge, distance, branch);
    return screenMode === "vj-mixture" ? tree.nucleotideMixtureScore(surface, vjProfile) : tree.guideScore(surface, references.guide);
  };
  const edgeScores: Array<{ edge: number; guideScore: number; distance: number; branch: number }> = [];
  const requestedFullHmmEdges = Math.floor(search.fullHmmEdges);
  const screenRefinementCount = requestedFullHmmEdges <= 0
    ? tree.edges.length
    : Math.min(tree.edges.length, Math.max(requestedFullHmmEdges, 8));
  const screenTotal = tree.edges.length + screenRefinementCount;
  progress("edge-screen", 0, screenTotal, screenMode === "vj-mixture" ? "Screening every edge with the independent-site V/J nucleotide mixture" : "Screening every edge with the single germline guide");
  for (let edgeIndex = 0; edgeIndex < tree.edges.length; edgeIndex += 1) {
    const edge = tree.edges[edgeIndex];
    let guideScore = NEGATIVE_INFINITY;
    let bestDistance = edge.length / 2;
    let bestBranch = 0;
    for (const distance of linearGrid(screenEdgePoints, edge.length)) for (const branch of guideBranchSamples) {
      const score = screenScore(edgeIndex, distance, branch);
      if (score > guideScore) {
        guideScore = score;
        bestDistance = distance;
        bestBranch = branch;
      }
    }
    edgeScores.push({ edge: edgeIndex, guideScore, distance: bestDistance, branch: bestBranch });
    progress("edge-screen", edgeIndex + 1, screenTotal, `Screened ${edge.id}`);
    if ((edgeIndex & 3) === 3) await Promise.resolve();
  }
  edgeScores.sort((left, right) => right.guideScore - left.guideScore);
  // Refine the cheap V/J-only surface continuously. This is only an initializer:
  // every retained ML/grid point is still evaluated with the complete HMM.
  for (let index = 0; index < screenRefinementCount; index += 1) {
    const candidate = edgeScores[index];
    const edge = tree.edges[candidate.edge];
    let distance = candidate.distance;
    let branch = candidate.branch;
    for (let round = 0; round < 2; round += 1) {
      const branchOptimum = maximizeUnitInterval((coordinate) => screenScore(candidate.edge, distance, branchFromUnit(coordinate)), 0.004, 16);
      branch = branchFromUnit(branchOptimum.value);
      const distanceOptimum = maximizeUnitInterval((fraction) => screenScore(candidate.edge, edge.length * fraction, branch), 0.004, 16);
      distance = edge.length * distanceOptimum.value;
    }
    candidate.distance = distance;
    candidate.branch = branch;
    candidate.guideScore = screenScore(candidate.edge, distance, branch);
    progress("edge-screen", tree.edges.length + index + 1, screenTotal, `Refined ${edge.id} with the V/J-only surface`);
    if ((index & 1) === 1) await Promise.resolve();
  }
  edgeScores.sort((left, right) => right.guideScore - left.guideScore);
  const selectedEdges = edgeScores.slice(0, requestedFullHmmEdges <= 0 ? tree.edges.length : Math.max(1, Math.min(tree.edges.length, requestedFullHmmEdges)));
  const screenByEdge = new Map(edgeScores.map((entry) => [entry.edge, entry]));
  const evaluated = new Map<string, PhyloUcaPlacement>();
  // The germline/D automaton is invariant across tree placements. Reuse it for
  // ML, grid, and Gibbs evaluations instead of rebuilding tens of thousands of
  // states at every proposed point.
  const reusableHmm = new PhyloUcaHmmGibbsSampler(references, input.options.hmm);
  const totalEdgeLength = tree.edges.reduce((sum, edge) => sum + Math.max(1e-12, edge.length), 0);
  const branchMean = Math.max(1e-5, Number.isFinite(search.branchPriorMean) ? search.branchPriorMean : 0.06);
  const edgePriorLog = (edgeIndex: number) => (search.edgePrior ?? "uniform-length") === "uniform-length"
    ? Math.log(Math.max(1e-12, tree.edges[edgeIndex].length) / Math.max(1e-12, totalEdgeLength))
    : -Math.log(tree.edges.length);
  const branchPriorLog = (branch: number) => -branch / branchMean - Math.log(branchMean);
  const placement = (edgeIndex: number, distance: number, branch: number, guideScore: number, likelihood: number): PhyloUcaPlacement => {
    const edge = tree.edges[edgeIndex];
    const safeDistance = Math.max(0, Math.min(edge.length, distance));
    const safeBranch = Math.max(0, Math.min(maximumUcaBranchLength, branch));
    return {
      edgeId: edge.id,
      endpointA: edge.endpointA,
      endpointB: edge.endpointB,
      distanceFromA: safeDistance,
      edgeLength: edge.length,
      edgeFraction: edge.length > 0 ? safeDistance / edge.length : 0,
      ucaBranchLength: safeBranch,
      logMarginalLikelihood: likelihood,
      logPosteriorScore: likelihood + edgePriorLog(edgeIndex) + branchPriorLog(safeBranch),
      localPosteriorWeight: 0,
      screenScore: guideScore,
      screenMode,
      guideScore,
    };
  };
  const evaluate = (edgeIndex: number, distance: number, branch: number, guideScore = screenByEdge.get(edgeIndex)?.guideScore ?? NEGATIVE_INFINITY): PhyloUcaPlacement => {
    const edge = tree.edges[edgeIndex];
    const safeDistance = Math.max(0, Math.min(edge.length, distance));
    const safeBranch = Math.max(0, Math.min(maximumUcaBranchLength, branch));
    const key = placementKey(edgeIndex, safeDistance, safeBranch);
    const previous = evaluated.get(key);
    if (previous) return previous;
    const surface = tree.conditionalLikelihoods(edgeIndex, safeDistance, safeBranch);
    const result = placement(edgeIndex, safeDistance, safeBranch, guideScore, reusableHmm.logMarginal(surface));
    evaluated.set(key, result);
    return result;
  };

  let placements: PhyloUcaPlacement[] = [];
  let local: PhyloUcaPlacement[] = [];
  let weights: number[] = [];
  let evaluatedUcaBranchLengths: number[] | undefined;
  let mcmcDiagnostics: PhyloUcaMcmcDiagnostics | undefined;
  let dCountPosterior: PhyloUcaDCountPosteriorPoint[] | undefined;
  const mixture = Array.from({ length: observedColumns }, () => [0, 0, 0, 0, 0] as [number, number, number, number, number]);
  let codonMixture: Array<{ startSite: number; probabilities: number[] }> = [];
  let bestPosterior: ReturnType<typeof phyloUcaHmmPosterior> | null = null;
  const marginalTrackMixture = new Map<string, AnnotationTrackMixture>();
  const prefilteredMarginalTrackIds = new Set<string>();
  let prefilteredMarginalMaximumWeight = 0;

  if (inferenceMode === "maximum-likelihood") {
    const initialTasks = selectedEdges.flatMap((entry) => uniqueNumbers([0, entry.branch]).map((branch) => ({ ...entry, branch })));
    progress("hmm-search", 0, initialTasks.length, "Evaluating full-HMM starting points for continuous conditional-ML optimization");
    for (let index = 0; index < initialTasks.length; index += 1) {
      const task = initialTasks[index];
      evaluate(task.edge, task.distance, task.branch, task.guideScore);
      progress("hmm-search", index + 1, initialTasks.length, `Full-HMM ML initializer ${index + 1} of ${initialTasks.length}`);
      if ((index & 1) === 1) await Promise.resolve();
    }
    // Once an edge is admitted by the user-controlled screen breadth, optimize
    // every admitted edge. Dropping edges after one full-HMM initializer can
    // miss a narrow near-zero pendant optimum—the failure this route exists to
    // avoid.
    const refineEdgeIds = new Set(selectedEdges.map((entry) => tree.edges[entry.edge].id));
    const rounds = Math.max(1, Math.floor(search.mlOptimizationRounds ?? 2));
    const tolerance = Math.max(1e-5, search.mlOptimizationTolerance ?? 0.002);
    let refined = 0;
    for (const edgeId of refineEdgeIds) {
      const edge = tree.edges.find((candidate) => candidate.id === edgeId)!;
      let current = [...evaluated.values()].filter((point) => point.edgeId === edgeId).sort((left, right) => right.logMarginalLikelihood - left.logMarginalLikelihood)[0];
      for (let round = 0; round < rounds; round += 1) {
        const branchOptimum = maximizeUnitInterval((coordinate) => evaluate(edge.index, current.distanceFromA, branchFromUnit(coordinate)).logMarginalLikelihood, tolerance);
        current = evaluate(edge.index, current.distanceFromA, branchFromUnit(branchOptimum.value));
        const distanceOptimum = maximizeUnitInterval((fraction) => evaluate(edge.index, edge.length * fraction, current.ucaBranchLength).logMarginalLikelihood, tolerance);
        current = evaluate(edge.index, edge.length * distanceOptimum.value, current.ucaBranchLength);
      }
      refined += 1;
      progress("hmm-search", refined, refineEdgeIds.size, `Continuously optimized ${edge.id} under the full HMM`);
      await Promise.resolve();
    }
    placements = [...evaluated.values()].sort((left, right) => right.logMarginalLikelihood - left.logMarginalLikelihood);
    if (!placements.length) throw new Error("No UCA attachment placement could be evaluated.");
    local = placements.slice(0, 1);
    weights = [1];
    local[0].localPosteriorWeight = 1;
    warnings.push("Conditional-ML mode reports the single full-HMM likelihood optimum; it does not average over tree attachment or UCA branch length, and placement/branch priors do not affect the optimum.");
  } else if (inferenceMode === "grid-marginalization") {
    const branchLengths = phyloUcaBranchLengthGrid(search);
    evaluatedUcaBranchLengths = branchLengths;
    const branchWidths = gridCellWidths(branchLengths, 0, maximumUcaBranchLength);
    const edgeFractions = linearGrid(Math.max(3, Math.floor(search.edgeGridPoints)), 1);
    const fractionWidths = gridCellWidths(edgeFractions, 0, 1);
    const tasks = selectedEdges.flatMap((entry) => edgeFractions.flatMap((fraction, fractionIndex) => branchLengths.map((branch, branchIndex) => ({ ...entry, fraction, fractionIndex, branch, branchIndex }))));
    progress("hmm-search", 0, tasks.length, "Evaluating the explicit full-HMM attachment × pendant-length quadrature grid");
    for (let index = 0; index < tasks.length; index += 1) {
      const task = tasks[index];
      const point = evaluate(task.edge, tree.edges[task.edge].length * task.fraction, task.branch, task.guideScore);
      point.integrationLogWeight = Math.log(fractionWidths[task.fractionIndex]) + Math.log(branchWidths[task.branchIndex]);
      progress("hmm-search", index + 1, tasks.length, `Full-HMM grid point ${index + 1} of ${tasks.length}`);
      if ((index & 1) === 1) await Promise.resolve();
    }
    placements = [...evaluated.values()].sort((left, right) =>
      right.logPosteriorScore + (right.integrationLogWeight ?? 0) - left.logPosteriorScore - (left.integrationLogWeight ?? 0));
    if (!placements.length) throw new Error("No UCA attachment placement could be evaluated.");
    const localCount = Math.max(1, Math.min(placements.length, Math.floor(search.localPosteriorPoints ?? 12)));
    local = placements.slice(0, localCount);
    weights = normalizedWeights(local.map((point) => point.logPosteriorScore + (point.integrationLogWeight ?? 0)));
    local.forEach((point, index) => { point.localPosteriorWeight = weights[index]; });
    if (localCount < placements.length) {
      const allWeights = normalizedWeights(placements.map((point) => point.logPosteriorScore + (point.integrationLogWeight ?? 0)));
      const retainedMass = allWeights.slice(0, localCount).reduce((sum, value) => sum + value, 0);
      warnings.push(`Grid posterior output retained the leading ${localCount.toLocaleString()} of ${placements.length.toLocaleString()} evaluated points (${(100 * retainedMass).toFixed(3)}% of quadrature mass) for character/HMM marginalization.`);
    }
  } else {
    const random = seededRandom(search.mcmcSeed ?? 1729);
    const iterations = Math.max(2, Math.floor(search.mcmcIterations ?? 160));
    const burnIn = Math.max(0, Math.min(iterations - 1, Math.floor(search.mcmcBurnIn ?? 40)));
    const thin = Math.max(1, Math.floor(search.mcmcThin ?? 2));
    const mhSteps = Math.max(1, Math.floor(search.mcmcMhStepsPerIteration ?? 4));
    const collapsedRefreshInterval = Math.max(0, Math.floor(search.mcmcCollapsedRefreshInterval ?? 3));
    const sampler = reusableHmm;
    let currentEdge = edgeScores[0].edge;
    let currentFraction = tree.edges[currentEdge].length > 0 ? edgeScores[0].distance / tree.edges[currentEdge].length : 0;
    let currentBranch = edgeScores[0].branch;
    const maximumScreen = Math.max(...edgeScores.map((entry) => entry.guideScore));
    const edgeProposalWeights = edgeScores.map((entry) => Math.exp(Math.max(-30, (entry.guideScore - maximumScreen) / 8)) + 1e-12);
    const edgeProposalTotal = edgeProposalWeights.reduce((sum, value) => sum + value, 0);
    const proposalProbability = new Map(edgeScores.map((entry, index) => [entry.edge, edgeProposalWeights[index] / edgeProposalTotal]));
    const sampleGlobalEdge = () => {
      let threshold = random() * edgeProposalTotal;
      for (let index = 0; index < edgeScores.length; index += 1) {
        threshold -= edgeProposalWeights[index];
        if (threshold <= 0) return edgeScores[index].edge;
      }
      return edgeScores.at(-1)!.edge;
    };
    const globalPositionMixture = Math.max(0, Math.min(1, search.mcmcGlobalPositionMixture ?? 0.85));
    const globalPositionScale = Math.max(1e-6, Math.min(0.5, search.mcmcGlobalPositionScale ?? 0.18));
    const wrappedUnit = (value: number) => ((value % 1) + 1) % 1;
    const screenFraction = (edgeIndex: number) => {
      const edge = tree.edges[edgeIndex];
      return edge.length > 0 ? wrappedUnit((screenByEdge.get(edgeIndex)?.distance ?? edge.length / 2) / edge.length) : 0;
    };
    const circularDistance = (left: number, right: number) => {
      const distance = Math.abs(wrappedUnit(left) - wrappedUnit(right));
      return Math.min(distance, 1 - distance);
    };
    const globalPositionDensity = (edgeIndex: number, fraction: number) => {
      const focused = circularDistance(fraction, screenFraction(edgeIndex)) <= globalPositionScale ? 1 / (2 * globalPositionScale) : 0;
      return (1 - globalPositionMixture) + globalPositionMixture * focused;
    };
    const sampleGlobalFraction = (edgeIndex: number) => random() < globalPositionMixture
      ? wrappedUnit(screenFraction(edgeIndex) + (2 * random() - 1) * globalPositionScale)
      : random();
    const logProposalDensity = (density: number) => density > 0 && Number.isFinite(density) ? Math.log(density) : NEGATIVE_INFINITY;
    const globalBranchMixture = Math.max(0, Math.min(1, search.mcmcGlobalBranchMixture ?? 0.9));
    const globalBranchMaximum = Math.max(1e-9, Math.min(maximumUcaBranchLength, search.mcmcGlobalBranchMaximum ?? 0.03));
    const sampleGlobalBranch = () => {
      if (!(maximumUcaBranchLength > 0)) return 0;
      return random() < globalBranchMixture ? random() * globalBranchMaximum : random() * maximumUcaBranchLength;
    };
    const globalBranchDensity = (branch: number) => {
      if (!(maximumUcaBranchLength > 0) || branch < 0 || branch > maximumUcaBranchLength) return 0;
      return (1 - globalBranchMixture) / maximumUcaBranchLength
        + (branch <= globalBranchMaximum ? globalBranchMixture / globalBranchMaximum : 0);
    };
    // The V/J screen determines the bounded initializer set, but the actual
    // starting state is selected under the complete HMM. This avoids spending
    // burn-in escaping a screen optimum whose junction explanation is poor.
    progress("hmm-search", 0, selectedEdges.length, "Selecting a Gibbs/MH initializer under the complete recombination HMM");
    let initializerScore = NEGATIVE_INFINITY;
    const initializerScores = new Map<number, number>();
    for (let index = 0; index < selectedEdges.length; index += 1) {
      const candidate = selectedEdges[index];
      const edge = tree.edges[candidate.edge];
      const fraction = edge.length > 0 ? candidate.distance / edge.length : 0;
      const surface = tree.conditionalLikelihoods(candidate.edge, candidate.distance, candidate.branch);
      const score = sampler.logMarginal(surface) + edgePriorLog(candidate.edge) + branchPriorLog(candidate.branch);
      initializerScores.set(candidate.edge, score);
      if (score > initializerScore) {
        initializerScore = score;
        currentEdge = candidate.edge;
        currentFraction = fraction;
        currentBranch = candidate.branch;
      }
      progress("hmm-search", index + 1, selectedEdges.length, `Full-HMM initializer ${index + 1} of ${selectedEdges.length}`);
      await Promise.resolve();
    }
    // Reuse those full-HMM evaluations to make the exact collapsed refresh a
    // materially better independence proposal. A residual guide-screen
    // component keeps every edge reachable, and its exact density enters the
    // Hastings ratio below; this affects efficiency, never the target model.
    const initializerMixture = Math.max(0, Math.min(1, search.mcmcCollapsedInitializerMixture ?? 0.95));
    const finiteInitializerScores = [...initializerScores.values()].filter(Number.isFinite);
    const maximumInitializerScore = finiteInitializerScores.length ? Math.max(...finiteInitializerScores) : 0;
    const rawInitializerWeights = new Map<number, number>();
    let initializerWeightTotal = 0;
    for (const [edgeIndex, score] of initializerScores) {
      const weight = Number.isFinite(score) ? Math.exp(Math.max(-60, score - maximumInitializerScore)) : 0;
      rawInitializerWeights.set(edgeIndex, weight);
      initializerWeightTotal += weight;
    }
    const collapsedProposalProbability = new Map<number, number>();
    for (const edge of tree.edges) {
      const broad = proposalProbability.get(edge.index) ?? 0;
      const informed = initializerWeightTotal > 0 ? (rawInitializerWeights.get(edge.index) ?? 0) / initializerWeightTotal : broad;
      collapsedProposalProbability.set(edge.index, (1 - initializerMixture) * broad + initializerMixture * informed);
    }
    const sampleCollapsedEdge = () => {
      let threshold = random();
      for (const edge of tree.edges) {
        threshold -= collapsedProposalProbability.get(edge.index) ?? 0;
        if (threshold <= 0) return edge.index;
      }
      return tree.edges.at(-1)!.index;
    };
    const collapsedProposalLogDensity = (edgeIndex: number, fraction: number, branch: number) =>
      logProposalDensity(collapsedProposalProbability.get(edgeIndex) ?? 0)
      + logProposalDensity(globalPositionDensity(edgeIndex, fraction))
      + logProposalDensity(globalBranchDensity(branch));
    const retained: Array<{ point: PhyloUcaPlacement; draw: ReturnType<PhyloUcaHmmGibbsSampler["draw"]> }> = [];
    const trace: PhyloUcaMcmcDiagnostics["trace"] = [];
    let branchProposals = 0;
    let branchAccepted = 0;
    let positionProposals = 0;
    let positionAccepted = 0;
    let globalProposals = 0;
    let globalAccepted = 0;
    let collapsedProposals = 0;
    let collapsedAccepted = 0;
    let edgeSwitches = 0;
    let gibbsDraws = 0;
    let gibbsMilliseconds = 0;
    let collapsedMarginalMilliseconds = 0;
    let conditionalMhMilliseconds = 0;
    let currentEdgeObject = tree.edges[currentEdge];
    let currentSurface = tree.conditionalLikelihoods(currentEdge, currentEdgeObject.length * currentFraction, currentBranch);
    const samplingStarted = performance.now();
    progress("mcmc", 0, iterations, "Running exact HMM Gibbs draws with continuous tree-position and branch-length MH updates");
    for (let iteration = 0; iteration < iterations; iteration += 1) {
      const gibbsStarted = performance.now();
      let draw = sampler.draw(currentSurface, random);
      gibbsMilliseconds += performance.now() - gibbsStarted;
      gibbsDraws += 1;
      if (collapsedRefreshInterval > 0 && iteration > 0 && iteration % collapsedRefreshInterval === 0) {
        collapsedProposals += 1;
        const proposedEdge = sampleCollapsedEdge();
        const proposedFraction = sampleGlobalFraction(proposedEdge);
        const proposedBranch = sampleGlobalBranch();
        const proposedEdgeObject = tree.edges[proposedEdge];
        const collapsedStarted = performance.now();
        const proposedSurface = tree.conditionalLikelihoods(proposedEdge, proposedEdgeObject.length * proposedFraction, proposedBranch);
        const proposedMarginalLikelihood = sampler.logMarginal(proposedSurface);
        collapsedMarginalMilliseconds += performance.now() - collapsedStarted;
        const currentTarget = draw.logMarginalLikelihood + edgePriorLog(currentEdge) + branchPriorLog(currentBranch);
        const proposedTarget = proposedMarginalLikelihood + edgePriorLog(proposedEdge) + branchPriorLog(proposedBranch);
        const logHastings = collapsedProposalLogDensity(currentEdge, currentFraction, currentBranch)
          - collapsedProposalLogDensity(proposedEdge, proposedFraction, proposedBranch);
        if (Math.log(Math.max(Number.MIN_VALUE, random())) < Math.min(0, proposedTarget - currentTarget + logHastings)) {
          collapsedAccepted += 1;
          if (proposedEdge !== currentEdge) edgeSwitches += 1;
          currentEdge = proposedEdge;
          currentFraction = proposedFraction;
          currentBranch = proposedBranch;
          currentEdgeObject = proposedEdgeObject;
          currentSurface = proposedSurface;
          const acceptedGibbsStarted = performance.now();
          draw = sampler.draw(proposedSurface, random);
          gibbsMilliseconds += performance.now() - acceptedGibbsStarted;
          gibbsDraws += 1;
        }
      }
      let conditionalTarget = fixedUcaLogLikelihood(currentSurface, draw.characterStates) + edgePriorLog(currentEdge) + branchPriorLog(currentBranch);
      const keep = iteration >= burnIn && (iteration - burnIn) % thin === 0;
      trace.push({
        iteration: iteration + 1,
        edgeId: currentEdgeObject.id,
        edgeFraction: currentFraction,
        ucaBranchLength: currentBranch,
        conditionalLogTarget: conditionalTarget,
        logMarginalLikelihood: draw.logMarginalLikelihood,
        retained: keep,
      });
      if (keep) retained.push({
        point: placement(currentEdge, currentEdgeObject.length * currentFraction, currentBranch, screenByEdge.get(currentEdge)?.guideScore ?? NEGATIVE_INFINITY, draw.logMarginalLikelihood),
        draw,
      });

      for (let step = 0; step < mhSteps; step += 1) {
        let proposedEdge = currentEdge;
        let proposedFraction = currentFraction;
        let proposedBranch = currentBranch;
        let logHastings = 0;
        let proposalKind: "branch" | "position" | "global";
        if (random() < Math.max(0, Math.min(1, search.mcmcGlobalJumpProbability ?? 0.12))) {
          proposalKind = "global";
          globalProposals += 1;
          proposedEdge = sampleGlobalEdge();
          proposedFraction = sampleGlobalFraction(proposedEdge);
          logHastings = logProposalDensity(proposalProbability.get(currentEdge) ?? 0)
            + logProposalDensity(globalPositionDensity(currentEdge, currentFraction))
            - logProposalDensity(proposalProbability.get(proposedEdge) ?? 0)
            - logProposalDensity(globalPositionDensity(proposedEdge, proposedFraction));
        } else if (random() < 0.5) {
          proposalKind = "branch";
          branchProposals += 1;
          const scale = Math.max(1e-9, search.mcmcBranchProposalScale ?? 0.004);
          proposedBranch = reflected(currentBranch + (2 * random() - 1) * scale, maximumUcaBranchLength);
        } else {
          proposalKind = "position";
          positionProposals += 1;
          const scale = Math.max(1e-9, search.mcmcPositionProposalScale ?? 0.18);
          proposedFraction = reflected(currentFraction + (2 * random() - 1) * scale, 1);
        }
        const proposedEdgeObject = tree.edges[proposedEdge];
        const conditionalStarted = performance.now();
        const proposedSurface = tree.conditionalLikelihoods(proposedEdge, proposedEdgeObject.length * proposedFraction, proposedBranch);
        const proposedTarget = fixedUcaLogLikelihood(proposedSurface, draw.characterStates) + edgePriorLog(proposedEdge) + branchPriorLog(proposedBranch);
        conditionalMhMilliseconds += performance.now() - conditionalStarted;
        if (Math.log(Math.max(Number.MIN_VALUE, random())) < Math.min(0, proposedTarget - conditionalTarget + logHastings)) {
          if (proposalKind === "branch") branchAccepted += 1;
          else if (proposalKind === "position") positionAccepted += 1;
          else {
            globalAccepted += 1;
            if (proposedEdge !== currentEdge) edgeSwitches += 1;
          }
          currentEdge = proposedEdge;
          currentFraction = proposedFraction;
          currentBranch = proposedBranch;
          currentEdgeObject = proposedEdgeObject;
          currentSurface = proposedSurface;
          conditionalTarget = proposedTarget;
        }
      }
      progress("mcmc", iteration + 1, iterations, `MCMC iteration ${iteration + 1} of ${iterations}; ${retained.length} retained draws`);
      await Promise.resolve();
    }
    const samplingMilliseconds = performance.now() - samplingStarted;
    if (!retained.length) throw new Error("The Gibbs/MH settings retained no posterior draws; reduce burn-in or thinning.");
    const sampleWeight = 1 / retained.length;
    const dCountSamples = Array.from({ length: Math.max(0, Math.floor(input.options.hmm.maximumDSegments)) + 1 }, () => 0);
    for (const sample of retained) {
      const dCount = sample.draw.path.filter((segment) => segment.kind === "D").length;
      while (dCountSamples.length <= dCount) dCountSamples.push(0);
      dCountSamples[dCount] += 1;
      sample.point.localPosteriorWeight = sampleWeight;
      addAnnotationTrackMixture(marginalTrackMixture, sample.draw.tracks, sampleWeight);
      for (let site = 0; site < observedColumns; site += 1) mixture[site][sample.draw.characterStates[site]] += sampleWeight;
      if (!codonMixture.length) for (let startSite = frameOffset; startSite + 2 < observedColumns; startSite += 3) {
        codonMixture.push({ startSite, probabilities: Array.from({ length: PHYLO_UCA_CODON_STATE_COUNT }, () => 0) });
      }
      for (const codon of codonMixture) {
        const first = sample.draw.characterStates[codon.startSite];
        const second = sample.draw.characterStates[codon.startSite + 1];
        const third = sample.draw.characterStates[codon.startSite + 2];
        codon.probabilities[25 * first + 5 * second + third] += sampleWeight;
      }
    }
    dCountPosterior = dCountSamples.map((samples, dCount) => ({ dCount, samples, probability: samples / retained.length }));
    placements = retained.map((sample) => sample.point).sort((left, right) => right.logMarginalLikelihood - left.logMarginalLikelihood);
    local = [...placements];
    weights = local.map(() => sampleWeight);
    const bestPlacement = placements[0];
    const bestEdge = tree.edges.find((edge) => edge.id === bestPlacement.edgeId)!;
    bestPosterior = phyloUcaHmmPosterior(tree.conditionalLikelihoods(bestEdge.index, bestPlacement.distanceFromA, bestPlacement.ucaBranchLength), references, input.options.hmm, frameOffset);
    const retainedTrace = trace.filter((point) => point.retained);
    const marginalTargets = retainedTrace.map((point) => {
      const edge = tree.edges.find((candidate) => candidate.id === point.edgeId)!;
      return point.logMarginalLikelihood + edgePriorLog(edge.index) + branchPriorLog(point.ucaBranchLength);
    });
    const branchEffectiveSampleSize = effectiveSampleSize(retainedTrace.map((point) => point.ucaBranchLength));
    const logTargetEffectiveSampleSize = effectiveSampleSize(marginalTargets);
    mcmcDiagnostics = {
      iterations,
      burnIn,
      thin,
      retainedSamples: retained.length,
      mhStepsPerIteration: mhSteps,
      seed: search.mcmcSeed ?? 1729,
      branchProposals,
      branchAccepted,
      positionProposals,
      positionAccepted,
      globalProposals,
      globalAccepted,
      collapsedProposals,
      collapsedAccepted,
      edgeSwitches,
      branchEffectiveSampleSize,
      logTargetEffectiveSampleSize,
      samplingMilliseconds,
      gibbsDraws,
      gibbsMilliseconds,
      collapsedMarginalMilliseconds,
      conditionalMhMilliseconds,
      trace,
    };
    if (Math.min(branchEffectiveSampleSize, logTargetEffectiveSampleSize) < 20) warnings.push(`This Gibbs/MH run retained ${retained.length} draws but at least one placement diagnostic had ESS below 20 (branch ${branchEffectiveSampleSize.toFixed(1)}, log target ${logTargetEffectiveSampleSize.toFixed(1)}). Treat placement-marginal summaries as provisional or increase iterations/collapsed-refresh frequency.`);
    warnings.push(`Gibbs/MH used continuous pendant lengths and continuous within-edge attachment fractions. No branch-length or attachment grid was used; each MH likelihood used the exact GTR transition at the proposed value and the tree's cached directed half-edge messages. Global jumps used an explicitly Hastings-corrected mixture of V/J-screen-centered and uniform within-edge positions${collapsedRefreshInterval > 0 ? `, with an exact collapsed refresh every ${collapsedRefreshInterval} iterations` : ""}.`);
  }

  if (inferenceMode !== "gibbs-mh") {
    progress("posterior", 0, local.length, inferenceMode === "maximum-likelihood" ? "Computing the exact posterior at the conditional-ML placement" : "Marginalizing UCA nucleotides and exact codons over retained grid points");
    for (let index = 0; index < local.length; index += 1) {
      const point = local[index];
      const edge = tree.edges.find((candidate) => candidate.id === point.edgeId)!;
      const surface = tree.conditionalLikelihoods(edge.index, point.distanceFromA, point.ucaBranchLength);
      const posterior = phyloUcaHmmPosterior(surface, references, input.options.hmm, frameOffset);
      if (index === 0) bestPosterior = posterior;
      addAnnotationTrackMixture(marginalTrackMixture, posterior.marginalTracks, weights[index]);
      for (const id of posterior.omittedMarginalTrackIds) prefilteredMarginalTrackIds.add(id);
      prefilteredMarginalMaximumWeight = Math.max(prefilteredMarginalMaximumWeight, posterior.omittedMarginalMaximumWeight);
      for (let site = 0; site < observedColumns; site += 1) for (let character = 0; character < 5; character += 1) mixture[site][character] += weights[index] * posterior.probabilities[site][character];
      if (!codonMixture.length) codonMixture = posterior.codonPosterior.map((codon) => ({ startSite: codon.startSite, probabilities: Array.from({ length: PHYLO_UCA_CODON_STATE_COUNT }, () => 0) }));
      if (posterior.codonPosterior.length !== codonMixture.length) throw new Error("UCA placements produced incompatible codon posterior dimensions.");
      for (let codon = 0; codon < codonMixture.length; codon += 1) {
        const source = posterior.codonPosterior[codon];
        const destination = codonMixture[codon];
        if (source.startSite !== destination.startSite || source.probabilities.length !== PHYLO_UCA_CODON_STATE_COUNT) throw new Error("UCA placements produced incompatible codon state orderings.");
        for (let state = 0; state < PHYLO_UCA_CODON_STATE_COUNT; state += 1) destination.probabilities[state] += weights[index] * source.probabilities[state];
      }
      progress("posterior", index + 1, local.length, `Integrated nucleotide and codon posterior for placement ${index + 1} of ${local.length}`);
      await Promise.resolve();
    }
  }
  if (!bestPosterior || !placements.length) throw new Error("The best UCA placement did not produce a posterior sequence.");
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
    const mixtureTotal = mixture[site].reduce((sum, value) => sum + value, 0);
    if (!(mixtureTotal > 0) || !Number.isFinite(mixtureTotal)) throw new Error(`The UCA posterior at alignment column ${original + 1} has no finite probability mass.`);
    const probabilities = normalizeProbabilityVector(mixture[site]) as [number, number, number, number, number];
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
  const codonPosterior: PhyloUcaCodonPosterior[] = codonMixture.map((codon, index) => {
    const probabilities = normalizeProbabilityVector(codon.probabilities);
    let mapState = 0;
    for (let state = 1; state < probabilities.length; state += 1) if (probabilities[state] > probabilities[mapState]) mapState = state;
    const columns = [0, 1, 2].map((offset) => (input.retainedColumns[codon.startSite + offset] ?? codon.startSite + offset) + 1) as [number, number, number];
    return {
      codonIndex: index + 1,
      alignmentColumns: columns,
      probabilities,
      mapCodon: PHYLO_UCA_CODON_SYMBOLS[mapState],
      mapProbability: probabilities[mapState],
      entropyBits: entropyBits(probabilities),
    };
  });
  const path = bestPosterior.path.map((segment) => ({
    ...segment,
    startColumn: (input.retainedColumns[segment.startColumn] ?? segment.startColumn) + 1,
    endColumn: (input.retainedColumns[segment.endColumn] ?? segment.endColumn) + 1,
  }));
  const viterbiAnnotationTracks = projectAnnotationTracks(bestPosterior.viterbiTracks, input.retainedColumns);
  const allMarginalAnnotationTracks = projectAnnotationTracks(finalizeAnnotationTrackMixture(marginalTrackMixture), input.retainedColumns);
  const prunedMarginalAnnotation = pruneMarginalAnnotationTracks(allMarginalAnnotationTracks);
  const retainedMarginalTrackIds = new Set(prunedMarginalAnnotation.retained.map((track) => track.id));
  for (const id of retainedMarginalTrackIds) prefilteredMarginalTrackIds.delete(id);
  for (const track of allMarginalAnnotationTracks) if (!retainedMarginalTrackIds.has(track.id)) prefilteredMarginalTrackIds.add(track.id);
  const bestPlacement = placements[0];
  const bestEdge = tree.edges.find((edge) => edge.id === bestPlacement.edgeId)!;
  const effectivePlacementCount = Math.exp(-weights.reduce((sum, weight) => weight > 0 ? sum + weight * Math.log(weight) : sum, 0));
  progress("finalize", 1, 1, "Preparing UCA sequence, placement tree, and provenance");
  return {
    schema: 6,
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
    frameOffset,
    bestPlacement,
    placements,
    evaluatedUcaBranchLengths,
    mcmcDiagnostics,
    dCountPosterior,
    mapAlignedSequence,
    mapUngappedSequence: mapAlignedSequence.replaceAll("-", ""),
    posteriorConsensusAligned,
    posterior,
    codonPosterior,
    hmmAnnotations: {
      minimumDisplayedWeight: PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT,
      viterbi: viterbiAnnotationTracks,
      marginalized: prunedMarginalAnnotation.retained,
      omittedMarginalTrackCount: prefilteredMarginalTrackIds.size,
      omittedMarginalMaximumWeight: Math.max(prefilteredMarginalMaximumWeight, prunedMarginalAnnotation.omittedMaximumWeight),
    },
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
