import { buildSparseEvidenceRow } from "./evidence.ts";
import type {
  AlleleRefinementOptions,
  ReferenceAlleleGraph,
  ReferenceAlleleNode,
  RefinementInputRow,
  RefinementModelSummary,
} from "./types.ts";

export interface ReferenceKernelAlternative {
  nodeIndex: number;
  names: string[];
  sequence: string;
  probability: number;
  distance: number;
  substitutionOnly: boolean;
}

export interface ReferenceKernelInspection {
  primary: ReferenceAlleleNode;
  primaryProbability: number;
  alternativeProbability: number;
  alternatives: ReferenceKernelAlternative[];
  truncated: boolean;
}

export interface AlignedReferenceRow {
  nodeIndex: number;
  names: string[];
  sequence: string;
  probability: number;
  primary: boolean;
  distance: number;
  substitutionOnly: boolean;
}

export interface ReferenceKernelAlignment {
  rows: AlignedReferenceRow[];
  columns: number;
  strippedAllGapColumns: number;
}

export interface AssignmentShiftDatum {
  nodeIndex: number;
  label: string;
  names: string[];
  before: number;
  after: number;
  delta: number;
  beforeAssignments: number;
  afterAssignments: number;
}

function inspectionRow(node: ReferenceAlleleNode, assumedShm: number): RefinementInputRow {
  return {
    ordinal: 0,
    sequenceId: "reference-kernel-inspection",
    datasetId: "reference-kernel-inspection",
    sampleId: "reference-kernel-inspection",
    subjectId: "reference-kernel-inspection",
    locus: node.locus,
    call: node.names[0] ?? "",
    score: 0,
    identity: 1 - assumedShm,
    shm: assumedShm,
    alternatives: "",
    abundance: 1,
  };
}

/**
 * Evaluate the exact local evidence kernel used by the repertoire model for a
 * selected database sequence. The selected node is the only literal call;
 * every other non-zero entry is therefore induced by reference-neighbour
 * leakage and is directly interpretable as a parameter consequence.
 */
export function inspectReferenceEvidenceKernel(
  graph: ReferenceAlleleGraph,
  selectedNode: number,
  assumedShm: number,
  options: AlleleRefinementOptions,
): ReferenceKernelInspection | null {
  const primary = graph.nodes[selectedNode];
  if (!primary) return null;
  const evidence = buildSparseEvidenceRow(
    inspectionRow(primary, Math.max(0, Math.min(0.95, assumedShm))),
    graph,
    options,
  );
  if (!evidence) return null;
  const neighbourByNode = new Map(graph.neighbours[selectedNode].map((neighbour) => [neighbour.index, neighbour] as const));
  const primaryProbability = evidence.entries.find((entry) => entry.node === selectedNode)?.weight ?? 0;
  const alternatives = evidence.entries
    .filter((entry) => entry.node !== selectedNode)
    .map((entry) => {
      const node = graph.nodes[entry.node];
      const neighbour = neighbourByNode.get(entry.node);
      return {
        nodeIndex: entry.node,
        names: [...node.names],
        sequence: node.sequence,
        probability: entry.weight,
        distance: neighbour?.distance ?? 0,
        substitutionOnly: neighbour?.substitutionOnly ?? false,
      };
    })
    .sort((left, right) => right.probability - left.probability
      || left.distance - right.distance
      || left.names.join(",").localeCompare(right.names.join(",")));
  return {
    primary,
    primaryProbability,
    alternativeProbability: alternatives.reduce((sum, alternative) => sum + alternative.probability, 0),
    alternatives,
    truncated: evidence.truncated,
  };
}

interface PairwiseAlignment {
  anchor: string;
  query: string;
}

/** Deterministic global alignment used only for the small diagnostic view. */
function alignPairToAnchor(anchor: string, query: string): PairwiseAlignment {
  const rows = anchor.length + 1;
  const columns = query.length + 1;
  const scores = new Int32Array(rows * columns);
  const trace = new Uint8Array(rows * columns); // 1 diagonal, 2 up, 3 left
  const at = (row: number, column: number) => row * columns + column;
  const gap = -2;
  for (let row = 1; row < rows; row += 1) {
    scores[at(row, 0)] = row * gap;
    trace[at(row, 0)] = 2;
  }
  for (let column = 1; column < columns; column += 1) {
    scores[at(0, column)] = column * gap;
    trace[at(0, column)] = 3;
  }
  for (let row = 1; row < rows; row += 1) {
    for (let column = 1; column < columns; column += 1) {
      const diagonal = scores[at(row - 1, column - 1)] + (anchor[row - 1] === query[column - 1] ? 2 : -1);
      const up = scores[at(row - 1, column)] + gap;
      const left = scores[at(row, column - 1)] + gap;
      // Diagonal preference keeps substitutions localized. The other tie is
      // deterministic, which prevents the display from jumping across renders.
      if (diagonal >= up && diagonal >= left) {
        scores[at(row, column)] = diagonal;
        trace[at(row, column)] = 1;
      } else if (up >= left) {
        scores[at(row, column)] = up;
        trace[at(row, column)] = 2;
      } else {
        scores[at(row, column)] = left;
        trace[at(row, column)] = 3;
      }
    }
  }
  const alignedAnchor: string[] = [];
  const alignedQuery: string[] = [];
  let row = anchor.length;
  let column = query.length;
  while (row > 0 || column > 0) {
    const direction = trace[at(row, column)];
    if (direction === 1) {
      alignedAnchor.push(anchor[row - 1]);
      alignedQuery.push(query[column - 1]);
      row -= 1;
      column -= 1;
    } else if (direction === 2 || column === 0) {
      alignedAnchor.push(anchor[row - 1]);
      alignedQuery.push("-");
      row -= 1;
    } else {
      alignedAnchor.push("-");
      alignedQuery.push(query[column - 1]);
      column -= 1;
    }
  }
  return { anchor: alignedAnchor.reverse().join(""), query: alignedQuery.reverse().join("") };
}

interface AnchorProjection {
  insertions: string[];
  bases: string[];
}

function projectOntoAnchor(anchorLength: number, alignment: PairwiseAlignment): AnchorProjection {
  const insertions = Array.from({ length: anchorLength + 1 }, () => "");
  const bases = Array.from({ length: anchorLength }, () => "-");
  let anchorPosition = 0;
  for (let column = 0; column < alignment.anchor.length; column += 1) {
    const anchorBase = alignment.anchor[column];
    const queryBase = alignment.query[column];
    if (anchorBase === "-") insertions[anchorPosition] += queryBase;
    else {
      bases[anchorPosition] = queryBase;
      anchorPosition += 1;
    }
  }
  return { insertions, bases };
}

function materializeProjection(anchor: string, projection: AnchorProjection, insertionWidths: number[]): string {
  let sequence = "";
  for (let position = 0; position <= anchor.length; position += 1) {
    const insertion = projection.insertions[position] ?? "";
    sequence += insertion + "-".repeat(Math.max(0, insertionWidths[position] - insertion.length));
    if (position < anchor.length) sequence += projection.bases[position] ?? "-";
  }
  return sequence;
}

/**
 * Construct one anchor-coherent reference alignment. Pairwise alignments are
 * merged through insertion slots around the selected reference rather than
 * progressively realigning previously placed sequences.
 */
export function alignReferenceKernelInspection(inspection: ReferenceKernelInspection): ReferenceKernelAlignment {
  const anchor = inspection.primary.sequence;
  const projections = inspection.alternatives.map((alternative) => projectOntoAnchor(
    anchor.length,
    alignPairToAnchor(anchor, alternative.sequence),
  ));
  const insertionWidths = Array.from({ length: anchor.length + 1 }, (_, position) => Math.max(
    0,
    ...projections.map((projection) => projection.insertions[position]?.length ?? 0),
  ));
  const anchorProjection: AnchorProjection = {
    insertions: Array.from({ length: anchor.length + 1 }, () => ""),
    bases: [...anchor],
  };
  const unstripped = [
    {
      nodeIndex: inspection.primary.index,
      names: [...inspection.primary.names],
      sequence: materializeProjection(anchor, anchorProjection, insertionWidths),
      probability: inspection.primaryProbability,
      primary: true,
      distance: 0,
      substitutionOnly: true,
    },
    ...inspection.alternatives.map((alternative, index) => ({
      nodeIndex: alternative.nodeIndex,
      names: [...alternative.names],
      sequence: materializeProjection(anchor, projections[index], insertionWidths),
      probability: alternative.probability,
      primary: false,
      distance: alternative.distance,
      substitutionOnly: alternative.substitutionOnly,
    })),
  ];
  const width = unstripped[0]?.sequence.length ?? 0;
  const retainedColumns: number[] = [];
  for (let column = 0; column < width; column += 1) {
    if (unstripped.some((row) => row.sequence[column] !== "-")) retainedColumns.push(column);
  }
  const rows = unstripped.map((row) => ({
    ...row,
    sequence: retainedColumns.map((column) => row.sequence[column]).join(""),
  }));
  return {
    rows,
    columns: retainedColumns.length,
    strippedAllGapColumns: width - retainedColumns.length,
  };
}

/** Frequencies before and after repertoire pooling, excluding Dirichlet prior mass. */
export function assignmentShiftData(model: RefinementModelSummary): AssignmentShiftDatum[] {
  const beforeTotal = model.alleles.reduce((sum, allele) => sum + Math.max(0, allele.localEvidenceAssignments), 0);
  const afterTotal = model.alleles.reduce((sum, allele) => sum + Math.max(0, allele.expectedAssignments), 0);
  return model.alleles.map((allele) => {
    const before = beforeTotal > 0 ? Math.max(0, allele.localEvidenceAssignments) / beforeTotal : 0;
    const after = afterTotal > 0 ? Math.max(0, allele.expectedAssignments) / afterTotal : 0;
    return {
      nodeIndex: allele.nodeIndex,
      label: allele.names.join(", "),
      names: [...allele.names],
      before,
      after,
      delta: after - before,
      beforeAssignments: Math.max(0, allele.localEvidenceAssignments),
      afterAssignments: Math.max(0, allele.expectedAssignments),
    };
  }).sort((left, right) => right.after - left.after
    || right.before - left.before
    || left.label.localeCompare(right.label, undefined, { numeric: true }));
}
