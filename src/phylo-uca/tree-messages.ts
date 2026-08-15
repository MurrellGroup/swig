import { parseNewick, serializeNewick, type TreeNode } from "../phylogeny.ts";
import { parseFasta } from "../post-analysis-core.ts";
import {
  compileGtr,
  observedCharacterPartial,
  transportLikelihood,
  type ReversibleCharacterModel,
} from "./gtr.ts";
import type { PhyloUcaGtrModel, PhyloUcaOptions } from "./types.ts";

const LIKELIHOOD_FLOOR = 1e-300;

interface GraphNode {
  id: number;
  name: string;
  edges: number[];
  leaf: boolean;
}

export interface PhyloUcaTreeEdge {
  id: string;
  index: number;
  a: number;
  b: number;
  length: number;
  endpointA: string;
  endpointB: string;
}

interface CavityMessage {
  values: Float64Array;
  scales: Float64Array;
}

function directedKey(node: number, blockedEdge: number): string {
  return `${node}:${blockedEdge}`;
}

function maximum(values: ArrayLike<number>, offset: number, length: number): number {
  let result = 0;
  for (let index = 0; index < length; index += 1) result = Math.max(result, values[offset + index]);
  return result;
}

function normalizeSite(values: Float64Array, offset: number, dimension: number): number {
  const normalizer = maximum(values, offset, dimension);
  if (!(normalizer > 0)) {
    for (let state = 0; state < dimension; state += 1) values[offset + state] = LIKELIHOOD_FLOOR;
    return Math.log(LIKELIHOOD_FLOOR);
  }
  for (let state = 0; state < dimension; state += 1) values[offset + state] /= normalizer;
  return Math.log(normalizer);
}

function safeLength(value: number): number {
  return Number.isFinite(value) && value > 0 ? value : 0;
}

export interface ConditionalLikelihoodSurface {
  /** Site-major log likelihoods; state order A,C,G,T[,gap]. */
  logLikelihoods: Float64Array;
  sites: number;
  stateCount: 4 | 5;
}

export class PhyloUcaTreeMessages {
  readonly sites: number;
  readonly model: ReversibleCharacterModel;
  readonly nodes: GraphNode[];
  readonly edges: PhyloUcaTreeEdge[];
  readonly alignmentNames: string[];
  readonly alignmentSequences: string[];
  readonly characterModel: "nucleotide-gtr4" | "gap-aware-gtr5";
  private readonly cavities = new Map<string, CavityMessage>();

  constructor(
    alignmentFasta: string,
    newick: string,
    model: PhyloUcaGtrModel,
    characterMode: PhyloUcaOptions["characterMode"] = "auto",
    onMessage?: (complete: number, total: number) => void,
  ) {
    const records = parseFasta(alignmentFasta, true);
    if (records.length < 2) throw new Error("Phylogenetic UCA inference needs at least two observed aligned sequences.");
    this.sites = records[0].sequence.length;
    if (!this.sites || records.some((record) => record.sequence.length !== this.sites)) {
      throw new Error("The observed-only phylogenetic alignment is not rectangular.");
    }
    this.alignmentNames = records.map((record) => record.name);
    this.alignmentSequences = records.map((record) => record.sequence.toUpperCase().replaceAll(".", "-").replaceAll("U", "T"));
    const observedGap = this.alignmentSequences.some((sequence) => sequence.includes("-"));
    this.characterModel = characterMode === "auto"
      ? observedGap ? "gap-aware-gtr5" : "nucleotide-gtr4"
      : characterMode;
    this.model = compileGtr(model, this.characterModel === "gap-aware-gtr5");
    this.nodes = [];
    this.edges = [];
    const root = parseNewick(newick);
    const add = (tree: TreeNode, parent: number | null): number => {
      const id = this.nodes.length;
      this.nodes.push({ id, name: tree.name, edges: [], leaf: tree.children.length === 0 });
      if (parent !== null) {
        const index = this.edges.length;
        const edge: PhyloUcaTreeEdge = {
          id: `edge_${index + 1}`,
          index,
          a: parent,
          b: id,
          length: safeLength(tree.length),
          endpointA: "",
          endpointB: "",
        };
        this.edges.push(edge);
        this.nodes[parent].edges.push(index);
        this.nodes[id].edges.push(index);
      }
      for (const child of tree.children) add(child, id);
      return id;
    };
    add(root, null);
    const alignmentByName = new Map(records.map((record, index) => [record.name, index]));
    const treeLeaves = this.nodes.filter((node) => node.leaf);
    const missing = treeLeaves.filter((node) => !alignmentByName.has(node.name)).map((node) => node.name);
    const absent = records.filter((record) => !treeLeaves.some((node) => node.name === record.name)).map((record) => record.name);
    if (missing.length || absent.length) {
      throw new Error(`Observed tree/alignment tip names differ${missing.length ? `; absent from alignment: ${missing.slice(0, 4).join(", ")}` : ""}${absent.length ? `; absent from tree: ${absent.slice(0, 4).join(", ")}` : ""}.`);
    }
    const nodeLabel = (node: GraphNode): string => node.leaf ? node.name : `internal_${node.id + 1}`;
    for (const edge of this.edges) {
      edge.endpointA = nodeLabel(this.nodes[edge.a]);
      edge.endpointB = nodeLabel(this.nodes[edge.b]);
    }
    const total = this.edges.length * 2;
    let completed = 0;
    for (const edge of this.edges) {
      this.cavity(edge.a, edge.index, alignmentByName);
      completed += 1;
      onMessage?.(completed, total);
      this.cavity(edge.b, edge.index, alignmentByName);
      completed += 1;
      onMessage?.(completed, total);
    }
  }

  private other(edge: PhyloUcaTreeEdge, node: number): number {
    return edge.a === node ? edge.b : edge.a;
  }

  private cavity(nodeId: number, blockedEdge: number, alignmentByName: ReadonlyMap<string, number>): CavityMessage {
    const key = directedKey(nodeId, blockedEdge);
    const cached = this.cavities.get(key);
    if (cached) return cached;
    const node = this.nodes[nodeId];
    const dimension = this.model.dimension;
    const values = new Float64Array(this.sites * dimension);
    const scales = new Float64Array(this.sites);
    if (node.leaf) {
      const row = alignmentByName.get(node.name);
      if (row === undefined) throw new Error(`Tree tip ${node.name} has no alignment row.`);
      const sequence = this.alignmentSequences[row];
      for (let site = 0; site < this.sites; site += 1) {
        const partial = observedCharacterPartial(sequence[site], dimension);
        values.set(partial, site * dimension);
      }
    } else {
      values.fill(1);
    }
    // Publish a placeholder before recursion. A proper tree cannot revisit this
    // directed half-edge; the placeholder turns malformed cycles into a finite
    // error rather than unbounded recursion.
    const result = { values, scales };
    this.cavities.set(key, result);
    for (const edgeIndex of node.edges) {
      if (edgeIndex === blockedEdge) continue;
      const edge = this.edges[edgeIndex];
      const neighbor = this.other(edge, nodeId);
      const incoming = this.cavity(neighbor, edgeIndex, alignmentByName);
      const transition = this.model.transition(edge.length);
      const source = new Float64Array(dimension);
      const transported = new Float64Array(dimension);
      for (let site = 0; site < this.sites; site += 1) {
        const offset = site * dimension;
        for (let state = 0; state < dimension; state += 1) source[state] = incoming.values[offset + state];
        transportLikelihood(transition, source, transported);
        let scale = incoming.scales[site];
        const transportMaximum = maximum(transported, 0, dimension);
        if (transportMaximum > 0) {
          scale += Math.log(transportMaximum);
          for (let state = 0; state < dimension; state += 1) transported[state] /= transportMaximum;
        }
        for (let state = 0; state < dimension; state += 1) values[offset + state] *= transported[state];
        scales[site] += scale;
      }
    }
    for (let site = 0; site < this.sites; site += 1) scales[site] += normalizeSite(values, site * dimension, dimension);
    return result;
  }

  conditionalLikelihoods(edgeIndex: number, rawDistanceFromA: number, rawUcaBranchLength: number): ConditionalLikelihoodSurface {
    const edge = this.edges[edgeIndex];
    if (!edge) throw new Error(`Unknown placement edge ${edgeIndex}.`);
    const distanceFromA = Math.max(0, Math.min(edge.length, Number.isFinite(rawDistanceFromA) ? rawDistanceFromA : 0));
    const distanceFromB = edge.length - distanceFromA;
    const branchLength = Math.max(0, Number.isFinite(rawUcaBranchLength) ? rawUcaBranchLength : 0);
    const cavityA = this.cavities.get(directedKey(edge.a, edgeIndex));
    const cavityB = this.cavities.get(directedKey(edge.b, edgeIndex));
    if (!cavityA || !cavityB) throw new Error(`Edge ${edge.id} is missing a directed likelihood message.`);
    const transitionA = this.model.transition(distanceFromA);
    const transitionB = this.model.transition(distanceFromB);
    const pendantTransition = this.model.transition(branchLength);
    const dimension = this.model.dimension;
    const logLikelihoods = new Float64Array(this.sites * dimension);
    const sourceA = new Float64Array(dimension);
    const sourceB = new Float64Array(dimension);
    const sideA = new Float64Array(dimension);
    const sideB = new Float64Array(dimension);
    const point = new Float64Array(dimension);
    const uca = new Float64Array(dimension);
    for (let site = 0; site < this.sites; site += 1) {
      const offset = site * dimension;
      for (let state = 0; state < dimension; state += 1) {
        sourceA[state] = cavityA.values[offset + state];
        sourceB[state] = cavityB.values[offset + state];
      }
      transportLikelihood(transitionA, sourceA, sideA);
      transportLikelihood(transitionB, sourceB, sideB);
      let pointMaximum = 0;
      for (let state = 0; state < dimension; state += 1) {
        point[state] = sideA[state] * sideB[state];
        pointMaximum = Math.max(pointMaximum, point[state]);
      }
      const commonScale = cavityA.scales[site] + cavityB.scales[site] + Math.log(Math.max(LIKELIHOOD_FLOOR, pointMaximum));
      for (let state = 0; state < dimension; state += 1) point[state] /= Math.max(LIKELIHOOD_FLOOR, pointMaximum);
      transportLikelihood(pendantTransition, point, uca);
      for (let state = 0; state < dimension; state += 1) {
        logLikelihoods[offset + state] = Math.log(Math.max(LIKELIHOOD_FLOOR, uca[state])) + commonScale;
      }
    }
    return { logLikelihoods, sites: this.sites, stateCount: dimension };
  }

  guideScore(surface: ConditionalLikelihoodSurface, alignedGuide: string): number {
    if (alignedGuide.length !== this.sites) throw new Error("The germline placement guide has different columns from the observed-only alignment.");
    let score = 0;
    const nucleotideTotal = this.model.frequencies[0] + this.model.frequencies[1] + this.model.frequencies[2] + this.model.frequencies[3];
    for (let site = 0; site < this.sites; site += 1) {
      const character = alignedGuide[site]?.toUpperCase().replace("U", "T").replace(".", "-") ?? "N";
      const exact = ["A", "C", "G", "T", "-"].indexOf(character);
      const offset = site * surface.stateCount;
      if (exact >= 0 && exact < surface.stateCount) {
        score += surface.logLikelihoods[offset + exact];
        continue;
      }
      let maximumLog = Number.NEGATIVE_INFINITY;
      for (let state = 0; state < 4; state += 1) maximumLog = Math.max(maximumLog, surface.logLikelihoods[offset + state]);
      let sum = 0;
      for (let state = 0; state < 4; state += 1) {
        sum += this.model.frequencies[state] / nucleotideTotal * Math.exp(surface.logLikelihoods[offset + state] - maximumLog);
      }
      score += maximumLog + Math.log(Math.max(LIKELIHOOD_FLOOR, sum));
    }
    return score;
  }

  placedTreeNewick(edgeIndex: number, rawDistanceFromA: number, rawUcaBranchLength: number, ucaName = "phylo_UCA"): string {
    const selected = this.edges[edgeIndex];
    if (!selected) throw new Error(`Unknown placement edge ${edgeIndex}.`);
    const distanceFromA = Math.max(0, Math.min(selected.length, rawDistanceFromA));
    const distanceFromB = selected.length - distanceFromA;
    const orient = (nodeId: number, parentId: number, length: number): TreeNode => {
      const node = this.nodes[nodeId];
      return {
        name: node.leaf ? node.name : "",
        length,
        children: node.edges.flatMap((edgeIndexValue) => {
          const edge = this.edges[edgeIndexValue];
          const next = this.other(edge, nodeId);
          return next === parentId ? [] : [orient(next, nodeId, edge.length)];
        }),
      };
    };
    // Root at the inferred UCA. As elsewhere in Swig, the named UCA leaf is a
    // zero-length label/sequence carrier; the complete pendant length belongs
    // to the branch from the UCA root to the observed-tree attachment point.
    const attachment: TreeNode = {
      name: "",
      length: Math.max(0, rawUcaBranchLength),
      children: [
        orient(selected.a, selected.b, distanceFromA),
        orient(selected.b, selected.a, distanceFromB),
      ],
    };
    const root: TreeNode = {
      name: "",
      length: 0,
      children: [
        { name: ucaName, length: 0, children: [] },
        attachment,
      ],
    };
    return `${serializeNewick(root)};`;
  }
}
