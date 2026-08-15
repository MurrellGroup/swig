import { parseFasta } from "../post-analysis-core.ts";
import type {
  ReferenceAlleleGraph,
  ReferenceAlleleNode,
  ReferenceNeighbour,
  RefinementSegment,
} from "./types.ts";

function normalizeReferenceSequence(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}

export function referenceLocus(name: string): string {
  return name.toUpperCase().match(/^(IGH|IGK|IGL|TRA|TRB|TRD|TRG)/)?.[1] ?? "";
}

/** Banded Levenshtein distance with an early cutoff; only tiny radii are used here. */
export function boundedReferenceDistance(left: string, right: string, maximum: number): number | null {
  if (Math.abs(left.length - right.length) > maximum) return null;
  if (left === right) return 0;
  if (left.length === right.length) {
    let distance = 0;
    for (let index = 0; index < left.length; index += 1) {
      if (left[index] === right[index]) continue;
      distance += 1;
      if (distance > maximum) return null;
    }
    return distance;
  }
  const columns = right.length + 1;
  let previous = new Uint16Array(columns);
  let current = new Uint16Array(columns);
  for (let column = 0; column < columns; column += 1) previous[column] = column;
  for (let row = 1; row <= left.length; row += 1) {
    current.fill(maximum + 1);
    current[0] = row;
    const begin = Math.max(1, row - maximum);
    const end = Math.min(right.length, row + maximum);
    let rowMinimum = maximum + 1;
    for (let column = begin; column <= end; column += 1) {
      const substitution = previous[column - 1] + (left[row - 1] === right[column - 1] ? 0 : 1);
      const deletion = previous[column] + 1;
      const insertion = current[column - 1] + 1;
      const value = Math.min(substitution, deletion, insertion);
      current[column] = value;
      rowMinimum = Math.min(rowMinimum, value);
    }
    if (rowMinimum > maximum) return null;
    [previous, current] = [current, previous];
  }
  const distance = previous[right.length];
  return distance <= maximum ? distance : null;
}

export function buildReferenceAlleleGraph(
  fasta: string,
  segment: RefinementSegment,
  maximumDistance: number,
): ReferenceAlleleGraph {
  const bySequence = new Map<string, ReferenceAlleleNode>();
  const callToNode = new Map<string, number>();
  let exactDuplicateLabels = 0;
  for (const record of parseFasta(fasta, true)) {
    const name = record.name.split(/\s+/, 1)[0]?.trim();
    const sequence = normalizeReferenceSequence(record.sequence);
    if (!name || !sequence) continue;
    const locus = referenceLocus(name);
    const key = `${locus}\u0000${sequence}`;
    let node = bySequence.get(key);
    if (!node) {
      node = { index: bySequence.size, segment, locus, names: [], sequence };
      bySequence.set(key, node);
    } else exactDuplicateLabels += 1;
    if (!node.names.includes(name)) node.names.push(name);
  }
  const nodes = [...bySequence.values()];
  nodes.forEach((node, index) => {
    node.index = index;
    node.names.sort();
    node.names.forEach((name) => callToNode.set(name, index));
  });
  const neighbours: ReferenceNeighbour[][] = nodes.map(() => []);
  const radius = Math.max(0, Math.min(4, Math.floor(maximumDistance)));
  if (radius > 0) {
    const byLocusLength = new Map<string, number[]>();
    nodes.forEach((node, index) => {
      const key = `${node.locus}\u0000${node.sequence.length}`;
      const values = byLocusLength.get(key);
      if (values) values.push(index); else byLocusLength.set(key, [index]);
    });
    const candidates = new Set<number>();
    const addCandidate = (left: number, right: number) => {
      const low = Math.min(left, right); const high = Math.max(left, right);
      candidates.add(low * nodes.length + high);
    };
    // A Hamming-radius-d pair must share at least one of d+1 exact blocks.
    // This replaces the dominant all-pairs scan for ordinary equal-length
    // allele sets while preserving exact verification below.
    for (const indices of byLocusLength.values()) {
      if (indices.length < 2) continue;
      const length = nodes[indices[0]].sequence.length;
      const blocks = Math.min(radius + 1, Math.max(1, length));
      for (let block = 0; block < blocks; block += 1) {
        const begin = Math.floor(block * length / blocks);
        const end = Math.floor((block + 1) * length / blocks);
        const buckets = new Map<string, number[]>();
        for (const index of indices) {
          const key = nodes[index].sequence.slice(begin, end);
          const previous = buckets.get(key);
          if (previous) {
            for (const other of previous) addCandidate(other, index);
            previous.push(index);
          } else buckets.set(key, [index]);
        }
      }
    }
    // Length-changing neighbours are uncommon and cannot use the exact block
    // guarantee without generating deletion variants. Limit their exact
    // pairwise scan to same-locus length groups inside the tiny edit radius.
    const groups = [...byLocusLength.entries()].map(([key, indices]) => {
      const split = key.lastIndexOf("\u0000");
      return { locus: key.slice(0, split), length: Number(key.slice(split + 1)), indices };
    });
    for (let leftGroup = 0; leftGroup < groups.length; leftGroup += 1) {
      for (let rightGroup = leftGroup + 1; rightGroup < groups.length; rightGroup += 1) {
        const left = groups[leftGroup]; const right = groups[rightGroup];
        if (left.locus !== right.locus || left.length === right.length || Math.abs(left.length - right.length) > radius) continue;
        for (const leftIndex of left.indices) for (const rightIndex of right.indices) addCandidate(leftIndex, rightIndex);
      }
    }
    for (const encoded of candidates) {
      const left = Math.floor(encoded / nodes.length);
      const right = encoded % nodes.length;
      const distance = boundedReferenceDistance(nodes[left].sequence, nodes[right].sequence, radius);
      if (distance === null || distance === 0) continue;
      const substitutionOnly = nodes[left].sequence.length === nodes[right].sequence.length;
      neighbours[left].push({ index: right, distance, substitutionOnly });
      neighbours[right].push({ index: left, distance, substitutionOnly });
    }
  }
  neighbours.forEach((values) => values.sort((left, right) => left.distance - right.distance || left.index - right.index));
  return { segment, nodes, callToNode, neighbours, exactDuplicateLabels };
}

export function referenceNodeLabel(node: ReferenceAlleleNode): string {
  return node.names.join(",");
}
