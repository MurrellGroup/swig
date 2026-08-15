import { parseFasta } from "../post-analysis-core.ts";
import {
  PHYLO_UCA_ANNOTATION_ALLELE_MINIMUM_WEIGHT,
  PHYLO_UCA_ANNOTATION_REGISTER_MINIMUM_WEIGHT,
  phyloUcaHmmLogMarginal,
  phyloUcaHmmPosterior,
} from "./hmm.ts";
import { preparePhyloUcaReferences } from "./references.ts";
import { PhyloUcaTreeMessages } from "./tree-messages.ts";
import { normalizeProbabilityVector } from "../probability-logo.ts";
import { PHYLO_UCA_CODON_STATE_COUNT, PHYLO_UCA_CODON_SYMBOLS } from "./codons.ts";
import type {
  PhyloUcaCodonPosterior,
  PhyloUcaHmmAnnotationTrack,
  PhyloUcaInput,
  PhyloUcaPlacement,
  PhyloUcaProgress,
  PhyloUcaResult,
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

function branchGrid(points: number, maximum: number): number[] {
  if (points <= 1 || !(maximum > 0)) return [0];
  return Array.from({ length: points }, (_, index) => maximum * Math.pow(index / (points - 1), 2));
}

function placementKey(edge: number, distance: number, branch: number): string {
  return `${edge}|${distance.toPrecision(11)}|${branch.toPrecision(11)}`;
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

  const screenMode = input.options.search.screenMode ?? "vj-mixture";
  // Endpoints plus at least one branch-interior attachment are always tested.
  const screenEdgePoints = Math.max(3, Math.floor(input.options.search.screenEdgeGridPoints ?? 5));
  const guideBranchSamples = uniqueNumbers([0, Math.min(0.02, input.options.search.maximumUcaBranchLength), Math.min(0.08, input.options.search.maximumUcaBranchLength)]);
  const vjProfile = vjNucleotideMixtureProfile(references);
  const edgeScores: Array<{ edge: number; guideScore: number; distance: number; branch: number }> = [];
  progress("edge-screen", 0, tree.edges.length, screenMode === "vj-mixture" ? "Screening every edge with the independent-site V/J nucleotide mixture" : "Screening every edge with the single germline guide");
  for (let edgeIndex = 0; edgeIndex < tree.edges.length; edgeIndex += 1) {
    const edge = tree.edges[edgeIndex];
    let guideScore = NEGATIVE_INFINITY;
    let bestDistance = edge.length / 2;
    let bestBranch = 0;
    for (const distance of linearGrid(screenEdgePoints, edge.length)) for (const branch of guideBranchSamples) {
      const surface = tree.conditionalLikelihoods(edgeIndex, distance, branch);
      const score = screenMode === "vj-mixture" ? tree.nucleotideMixtureScore(surface, vjProfile) : tree.guideScore(surface, references.guide);
      if (score > guideScore) {
        guideScore = score;
        bestDistance = distance;
        bestBranch = branch;
      }
    }
    edgeScores.push({ edge: edgeIndex, guideScore, distance: bestDistance, branch: bestBranch });
    progress("edge-screen", edgeIndex + 1, tree.edges.length, `Screened ${edge.id}`);
    if ((edgeIndex & 3) === 3) await Promise.resolve();
  }
  edgeScores.sort((left, right) => right.guideScore - left.guideScore);
  const requestedFullHmmEdges = Math.floor(input.options.search.fullHmmEdges);
  const selectedEdges = edgeScores.slice(0, requestedFullHmmEdges <= 0 ? tree.edges.length : Math.max(1, Math.min(tree.edges.length, requestedFullHmmEdges)));
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
      screenScore: guideScore,
      screenMode,
      guideScore,
    };
    evaluated.set(key, placement);
    return placement;
  };

  const coarseTasks = selectedEdges.flatMap(({ edge, guideScore, distance: screenDistance, branch: screenBranch }) => {
    const treeEdge = tree.edges[edge];
    const uniform = linearGrid(Math.max(3, input.options.search.edgeGridPoints), treeEdge.length).flatMap((distance) =>
      branchGrid(Math.max(2, input.options.search.branchGridPoints), input.options.search.maximumUcaBranchLength).map((branch) => ({ edge, distance, branch, guideScore }))
    );
    return [...uniform, { edge, distance: screenDistance, branch: screenBranch, guideScore }];
  });
  progress("hmm-search", 0, coarseTasks.length, "Evaluating the full recombination HMM on screen-ranked branch points");
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
  let codonMixture: Array<{ startSite: number; probabilities: number[] }> = [];
  let bestPosterior: ReturnType<typeof phyloUcaHmmPosterior> | null = null;
  const marginalTrackMixture = new Map<string, AnnotationTrackMixture>();
  const prefilteredMarginalTrackIds = new Set<string>();
  let prefilteredMarginalMaximumWeight = 0;
  progress("posterior", 0, local.length, "Marginalizing UCA nucleotides and exact codons around the best placement");
  for (let index = 0; index < local.length; index += 1) {
    const placement = local[index];
    const edge = tree.edges.find((candidate) => candidate.id === placement.edgeId)!;
    const surface = tree.conditionalLikelihoods(edge.index, placement.distanceFromA, placement.ucaBranchLength);
    const posterior = phyloUcaHmmPosterior(surface, references, input.options.hmm, frameOffset);
    if (index === 0) bestPosterior = posterior;
    addAnnotationTrackMixture(marginalTrackMixture, posterior.marginalTracks, weights[index]);
    for (const id of posterior.omittedMarginalTrackIds) prefilteredMarginalTrackIds.add(id);
    prefilteredMarginalMaximumWeight = Math.max(prefilteredMarginalMaximumWeight, posterior.omittedMarginalMaximumWeight);
    for (let site = 0; site < observedColumns; site += 1) for (let character = 0; character < 5; character += 1) mixture[site][character] += weights[index] * posterior.probabilities[site][character];
    if (!codonMixture.length) codonMixture = posterior.codonPosterior.map((codon) => ({ startSite: codon.startSite, probabilities: Array.from({ length: PHYLO_UCA_CODON_STATE_COUNT }, () => 0) }));
    if (posterior.codonPosterior.length !== codonMixture.length) throw new Error("Local UCA placements produced incompatible codon posterior dimensions.");
    for (let codon = 0; codon < codonMixture.length; codon += 1) {
      const source = posterior.codonPosterior[codon];
      const destination = codonMixture[codon];
      if (source.startSite !== destination.startSite || source.probabilities.length !== PHYLO_UCA_CODON_STATE_COUNT) throw new Error("Local UCA placements produced incompatible codon state orderings.");
      for (let state = 0; state < PHYLO_UCA_CODON_STATE_COUNT; state += 1) destination.probabilities[state] += weights[index] * source.probabilities[state];
    }
    progress("posterior", index + 1, local.length, `Integrated nucleotide and codon posterior for placement ${index + 1} of ${local.length}`);
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
    schema: 4,
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
