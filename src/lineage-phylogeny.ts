import type { AirrDetailRow } from "./result-store";
import type { FastaRecord } from "./post-analysis-core";
import type { TreeNode } from "./phylogeny";

export const VARIABLE_REGIONS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"] as const;
export type VariableRegion = typeof VARIABLE_REGIONS[number];

export interface BranchMutation {
  column: number;
  from: string;
  to: string;
  childClade: string;
}

export interface ParsimonyMap {
  ucaSequence: string;
  score: number;
  sequencesByClade: Map<string, string>;
  mutationsByClade: Map<string, BranchMutation[]>;
}

const STATES = ["A", "C", "G", "T", "-"] as const;
const NUCLEOTIDE_MASKS: Record<string, number> = {
  A: 1, C: 2, G: 4, T: 8, U: 8,
  R: 1 | 4, Y: 2 | 8, S: 2 | 4, W: 1 | 8,
  K: 4 | 8, M: 1 | 2, B: 2 | 4 | 8, D: 1 | 4 | 8,
  H: 1 | 2 | 8, V: 1 | 2 | 4, N: 1 | 2 | 4 | 8,
  X: 1 | 2 | 4 | 8, "?": 1 | 2 | 4 | 8, "-": 16, ".": 16,
};
const INF = 1_000_000;

function stateMask(value: string): number {
  return NUCLEOTIDE_MASKS[value.toUpperCase()] ?? 15;
}

function leafNames(node: TreeNode, target: string[] = []): string[] {
  if (!node.children.length) target.push(node.name);
  else node.children.forEach((child) => leafNames(child, target));
  return target;
}

/** A topology-stable node identifier, independent of Newick child order. */
export function cladeSignature(node: TreeNode): string {
  return leafNames(node).sort().join("\u0000");
}

function selectOutgroup(root: TreeNode, outgroupName: string): { ingroup: TreeNode; outgroup?: TreeNode } {
  const direct = root.children.find((child) => !child.children.length && child.name === outgroupName);
  if (direct && root.children.length === 2) return { ingroup: root.children.find((child) => child !== direct)!, outgroup: direct };
  return { ingroup: root };
}

function chooseState(costs: number[], preferred = -1): number {
  let best = Math.min(...costs);
  if (preferred >= 0 && costs[preferred] === best) return preferred;
  for (let state = 0; state < costs.length; state += 1) if (costs[state] === best) return state;
  return 0;
}

/**
 * Equal-cost Sankoff reconstruction rooted on the N-masked germline.
 *
 * Known germline bases constrain the inferred UCA. An N is an unknown member
 * of {A,C,G,T}: its UCA base is inferred from the descendants and no synthetic
 * N→base event is emitted. A gap is a fifth state, so indels remain visible.
 */
export function mapParsimonyMutations(
  root: TreeNode,
  alignedByName: Map<string, string>,
  outgroupName: string,
): ParsimonyMap {
  if (!alignedByName.size) throw new Error("Parsimony mapping requires an alignment.");
  const widths = new Set([...alignedByName.values()].map((sequence) => sequence.length));
  if (widths.size !== 1) throw new Error("Parsimony mapping requires equal-length aligned sequences.");
  const columns = widths.values().next().value as number;
  const { ingroup, outgroup } = selectOutgroup(root, outgroupName);
  const germline = alignedByName.get(outgroupName)?.toUpperCase() ?? "N".repeat(columns);
  const nodes: TreeNode[] = [];
  const collect = (node: TreeNode) => { nodes.push(node); node.children.forEach(collect); };
  collect(ingroup);
  const sequenceBuffers = new Map(nodes.map((node) => [node, [] as string[]]));
  const mutationBuffers = new Map<string, BranchMutation[]>();
  const ucaBuffer: string[] = [];
  let score = 0;

  for (let column = 0; column < columns; column += 1) {
    const costs = new Map<TreeNode, number[]>();
    const postorder = (node: TreeNode): number[] => {
      if (!node.children.length) {
        const mask = stateMask(alignedByName.get(node.name)?.[column] ?? "N");
        const leafCost = STATES.map((_, state) => mask & (1 << state) ? 0 : INF);
        costs.set(node, leafCost);
        return leafCost;
      }
      const childCosts = node.children.map(postorder);
      const current = STATES.map((_, parentState) => childCosts.reduce((sum, child) => {
        let best = INF;
        for (let childState = 0; childState < STATES.length; childState += 1) {
          best = Math.min(best, child[childState] + (childState === parentState ? 0 : 1));
        }
        return sum + best;
      }, 0));
      costs.set(node, current);
      return current;
    };
    const rootCosts = postorder(ingroup);
    const germlineMask = stateMask(germline[column] ?? "N");
    let ucaState = -1;
    if (germlineMask && (germlineMask & (germlineMask - 1)) === 0) {
      ucaState = Math.log2(germlineMask);
    } else {
      const constrained = STATES.map((_, candidate) => {
        if (!(germlineMask & (1 << candidate))) return INF;
        return Math.min(...rootCosts.map((value, ingroupState) => value + (ingroupState === candidate ? 0 : 1)));
      });
      ucaState = chooseState(constrained);
    }
    ucaBuffer.push(STATES[ucaState]);
    const ingroupTransitionCosts = rootCosts.map((value, state) => value + (state === ucaState ? 0 : 1));
    const ingroupState = chooseState(ingroupTransitionCosts, ucaState);
    score += ingroupTransitionCosts[ingroupState];
    if (ingroupState !== ucaState) {
      const signature = cladeSignature(ingroup);
      const mutations = mutationBuffers.get(signature) ?? [];
      mutations.push({ column, from: STATES[ucaState], to: STATES[ingroupState], childClade: signature });
      mutationBuffers.set(signature, mutations);
    }

    const preorder = (node: TreeNode, state: number) => {
      sequenceBuffers.get(node)!.push(STATES[state]);
      for (const child of node.children) {
        const childCost = costs.get(child)!;
        const transitionCost = childCost.map((value, childState) => value + (childState === state ? 0 : 1));
        const childState = chooseState(transitionCost, state);
        if (childState !== state) {
          const signature = cladeSignature(child);
          const mutations = mutationBuffers.get(signature) ?? [];
          mutations.push({ column, from: STATES[state], to: STATES[childState], childClade: signature });
          mutationBuffers.set(signature, mutations);
        }
        preorder(child, childState);
      }
    };
    preorder(ingroup, ingroupState);
  }

  const ucaSequence = ucaBuffer.join("");
  const sequencesByClade = new Map<string, string>();
  for (const [node, buffer] of sequenceBuffers) sequencesByClade.set(cladeSignature(node), buffer.join(""));
  if (outgroup) sequencesByClade.set(cladeSignature(root), ucaSequence);
  return { ucaSequence, score, sequencesByClade, mutationsByClade: mutationBuffers };
}

const CODONS: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S", TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
  CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P", CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T", AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A", GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

export type AlignmentFrameOffset = 0 | 1 | 2;

export interface AlignmentFrameEvidence {
  offset: AlignmentFrameOffset;
  completeGapCodons: number;
  mixedGapCodons: number;
  stopCodons: number;
  ambiguousCodons: number;
}

function normalizedFrameOffset(value: number): AlignmentFrameOffset {
  return value === 1 || value === 2 ? value : 0;
}

/**
 * Infer the shared codon phase of a nucleotide MSA without redefining a mixed
 * base/gap triplet as a gap. Codon-preserving editors emit complete `---`
 * triplets in the selected phase, so mixed gap codons are the primary signal;
 * stop codons and unresolved codons are conservative tie-breakers.
 */
export function inferAlignedReadingFrame(
  sequences: Iterable<string>,
  preferred: AlignmentFrameOffset = 0,
): AlignmentFrameEvidence {
  const normalized = [...sequences].map((sequence) => sequence.toUpperCase().replaceAll(".", "-"));
  const evidence = ([0, 1, 2] as AlignmentFrameOffset[]).map((offset): AlignmentFrameEvidence => {
    let completeGapCodons = 0;
    let mixedGapCodons = 0;
    let stopCodons = 0;
    let ambiguousCodons = 0;
    for (const sequence of normalized) {
      for (let index = offset; index + 2 < sequence.length; index += 3) {
        const codon = sequence.slice(index, index + 3);
        if (codon === "---") completeGapCodons += 1;
        else if (codon.includes("-")) mixedGapCodons += 1;
        else if (/[^ACGT]/.test(codon)) ambiguousCodons += 1;
        else if (CODONS[codon] === "*") stopCodons += 1;
      }
    }
    return { offset, completeGapCodons, mixedGapCodons, stopCodons, ambiguousCodons };
  });
  return evidence.sort((left, right) =>
    left.mixedGapCodons - right.mixedGapCodons ||
    left.stopCodons - right.stopCodons ||
    right.completeGapCodons - left.completeGapCodons ||
    left.ambiguousCodons - right.ambiguousCodons ||
    Number(right.offset === preferred) - Number(left.offset === preferred) ||
    left.offset - right.offset
  )[0];
}

export function translateAlignedNucleotides(sequence: string, frameOffset: number = 0): string {
  let result = "";
  for (let index = normalizedFrameOffset(frameOffset); index < sequence.length; index += 3) {
    const codon = sequence.slice(index, index + 3).toUpperCase();
    result += codon === "---" ? "-" : codon.includes("-") || /[^ACGT]/.test(codon) || codon.length < 3 ? "X" : CODONS[codon] ?? "X";
  }
  return result;
}

/**
 * Project a reconstructed nucleotide branch onto amino-acid replacements.
 * Multiple nucleotide events in one codon become one amino-acid event;
 * synonymous events and codons with an unresolved translation are omitted.
 */
export function aminoAcidBranchMutations(parentSequence: string, childSequence: string, childClade = "", frameOffset: number = 0): BranchMutation[] {
  if (parentSequence.length !== childSequence.length) throw new Error("Amino-acid branch mapping requires equal-length nucleotide states.");
  const parent = translateAlignedNucleotides(parentSequence, frameOffset);
  const child = translateAlignedNucleotides(childSequence, frameOffset);
  const mutations: BranchMutation[] = [];
  for (let column = 0; column < Math.min(parent.length, child.length); column += 1) {
    const from = parent[column];
    const to = child[column];
    if (from === to || from === "X" || to === "X") continue;
    mutations.push({ column, from, to, childClade });
  }
  return mutations;
}

export function ordinalFromAlignmentName(name: string): number | null {
  const match = name.match(/__(\d+)$/);
  return match ? Number(match[1]) - 1 : null;
}

function numeric(row: AirrDetailRow, key: string): number | null {
  const value = Number(row.values[key]);
  return Number.isFinite(value) && value > 0 ? Math.floor(value) : null;
}

function rowRegion(row: AirrDetailRow, position: number): VariableRegion | null {
  for (const region of VARIABLE_REGIONS) {
    const start = numeric(row, `${region}_start`);
    const end = numeric(row, `${region}_end`);
    if (start !== null && end !== null && position >= start && position <= end) return region;
  }
  return null;
}

function coordinateOffset(row: AirrDetailRow, alignedSequence: string): number {
  const ungapped = alignedSequence.replace(/-/g, "");
  const raw = row.values.sequence?.toUpperCase().replace(/[^ACGTN]/g, "") ?? "";
  const found = raw.indexOf(ungapped);
  if (found >= 0) return found;
  return Math.max(0, (numeric(row, "v_sequence_start") ?? 1) - 1);
}

/** Majority AIRR region label at every aligned nucleotide column. */
export function alignmentRegionMap(records: FastaRecord[], rows: AirrDetailRow[]): Array<VariableRegion | null> {
  const rowByOrdinal = new Map(rows.map((row) => [row.record.ordinal, row]));
  const columns = Math.max(0, ...records.map((record) => record.sequence.length));
  const votes = Array.from({ length: columns }, () => new Map<VariableRegion, number>());
  for (const record of records) {
    const ordinal = ordinalFromAlignmentName(record.name);
    const row = ordinal === null ? undefined : rowByOrdinal.get(ordinal);
    if (!row) continue;
    const offset = coordinateOffset(row, record.sequence);
    let residue = 0;
    for (let column = 0; column < record.sequence.length; column += 1) {
      if (record.sequence[column] === "-") continue;
      const region = rowRegion(row, offset + residue + 1);
      residue += 1;
      if (!region) continue;
      votes[column].set(region, (votes[column].get(region) ?? 0) + 1);
    }
  }
  return votes.map((columnVotes) => {
    let best: VariableRegion | null = null;
    let count = 0;
    for (const region of VARIABLE_REGIONS) {
      const next = columnVotes.get(region) ?? 0;
      if (next > count) { best = region; count = next; }
    }
    return best;
  });
}

export function aminoAcidRegionMap(nucleotideRegions: Array<VariableRegion | null>, frameOffset: number = 0): Array<VariableRegion | null> {
  const result: Array<VariableRegion | null> = [];
  for (let index = normalizedFrameOffset(frameOffset); index < nucleotideRegions.length; index += 3) {
    const votes = new Map<VariableRegion, number>();
    nucleotideRegions.slice(index, index + 3).forEach((region) => {
      if (region) votes.set(region, (votes.get(region) ?? 0) + 1);
    });
    result.push(VARIABLE_REGIONS.reduce<VariableRegion | null>((best, region) => (votes.get(region) ?? 0) > (best ? votes.get(best) ?? 0 : 0) ? region : best, null));
  }
  return result;
}

export type RegionPreset = "full" | "variable" | "cdrs" | VariableRegion | "custom";

export function columnsForRegionPreset(regions: Array<VariableRegion | null>, preset: RegionPreset): number[] {
  if (preset === "full") return regions.map((_, index) => index);
  if (preset === "variable") return regions.flatMap((region, index) => region ? [index] : []);
  if (preset === "cdrs") return regions.flatMap((region, index) => region?.startsWith("cdr") ? [index] : []);
  if (preset === "custom") return [];
  return regions.flatMap((region, index) => region === preset ? [index] : []);
}

/** Display offsets with a half-cell separator at every non-contiguous run. */
export function spacedColumnOffsets(columns: number[], cellWidth: number): number[] {
  let discontinuityOffset = 0;
  return columns.map((column, index) => {
    if (index > 0 && column !== columns[index - 1] + 1) discontinuityOffset += cellWidth * 0.5;
    return index * cellWidth + discontinuityOffset;
  });
}

/** Parse 1-based coordinates, ranges, or numbering labels such as 31, 31A, 32-35B. */
export function parseColumnSelection(text: string, labels: string[], maximum: number): number[] {
  const labelIndex = new Map(labels.map((label, index) => [label.toUpperCase(), index]));
  const resolve = (token: string): number | null => {
    const normalized = token.trim().toUpperCase();
    if (!normalized) return null;
    if (labelIndex.has(normalized)) return labelIndex.get(normalized)!;
    const numericValue = Number(normalized);
    return Number.isInteger(numericValue) && numericValue >= 1 && numericValue <= maximum ? numericValue - 1 : null;
  };
  const selected = new Set<number>();
  for (const token of text.split(/[\s,;]+/).filter(Boolean)) {
    const range = token.match(/^([0-9]+[A-Z]?)-([0-9]+[A-Z]?)$/i);
    if (range) {
      const start = resolve(range[1]);
      const end = resolve(range[2]);
      if (start !== null && end !== null) for (let index = Math.min(start, end); index <= Math.max(start, end); index += 1) selected.add(index);
      continue;
    }
    const value = resolve(token);
    if (value !== null) selected.add(value);
  }
  return [...selected].sort((left, right) => left - right);
}

const IUPAC_MATCH: Record<string, string> = {
  A: "A", C: "C", G: "G", T: "T", U: "T", R: "AG", Y: "CT", S: "CG", W: "AT",
  K: "GT", M: "AC", B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT", X: "ACGT", ".": "ACGT",
};

/** Return 1-based motif indices for aligned cells; zero means no motif. */
export function motifCellMap(sequence: string, motifs: string[], alphabet: "nt" | "aa"): Int16Array {
  const result = new Int16Array(sequence.length);
  const ungapped: string[] = [];
  const alignedColumns: number[] = [];
  [...sequence.toUpperCase()].forEach((value, column) => {
    if (value === "-") return;
    ungapped.push(value);
    alignedColumns.push(column);
  });
  const query = ungapped.join("");
  motifs.forEach((rawMotif, motifIndex) => {
    const motif = rawMotif.toUpperCase().replace(/\s|-/g, "");
    if (!motif) return;
    for (let start = 0; start + motif.length <= query.length; start += 1) {
      let matches = true;
      for (let offset = 0; offset < motif.length; offset += 1) {
        const expected = motif[offset];
        const observed = query[start + offset];
        if (alphabet === "nt" ? !(IUPAC_MATCH[expected] ?? expected).includes(observed) : expected !== "X" && expected !== "." && expected !== observed) {
          matches = false;
          break;
        }
      }
      if (!matches) continue;
      for (let offset = 0; offset < motif.length; offset += 1) {
        const column = alignedColumns[start + offset];
        if (!result[column]) result[column] = motifIndex + 1;
      }
    }
  });
  return result;
}
