import assert from "node:assert/strict";
import test from "node:test";

import { adaptiveNeighbourOdds, buildSparseEvidenceRow, sparseEvidenceMatrix } from "../src/allele-refinement/evidence.ts";
import { alignReferenceKernelInspection, hardAssignmentShiftData, inspectReferenceEvidenceKernel, survivingAlleleReference } from "../src/allele-refinement/diagnostics.ts";
import { fitSparseAlleleModel } from "../src/allele-refinement/model.ts";
import { buildReferenceAlleleGraph, boundedReferenceDistance } from "../src/allele-refinement/reference-graph.ts";
import { DEFAULT_ALLELE_REFINEMENT_OPTIONS, type RefinementInputRow, type SegmentRefinementResult } from "../src/allele-refinement/types.ts";
import { applyCallOverrides } from "../src/allele-refinement/apply.ts";
import { restoreAlleleRefinement, saveAlleleRefinement } from "../src/allele-refinement/serialization.ts";
import { modelSummaryTable, writeRefinementSidecar } from "../src/allele-refinement/export.ts";
import type { AirrResultStore } from "../src/result-store.ts";

const fasta = [
  ">IGHV1*01", "AAAAAAAAAA",
  ">IGHV1*01_duplicate_label", "AAAAAAAAAA",
  ">IGHV1*02", "AAAAAAAACA",
  ">IGHV1*03", "AAAAAAACCA",
  ">IGKV1*01", "AAAAAAAAAA",
].join("\n") + "\n";

const options = { ...DEFAULT_ALLELE_REFINEMENT_OPTIONS, segments: ["V" as const] };

function row(ordinal: number, call: string, alternatives = "", subjectId = "donor_1"): RefinementInputRow {
  return {
    ordinal,
    sequenceId: `r${ordinal}`,
    datasetId: "dataset_1",
    sampleId: "sample_1",
    subjectId,
    locus: "IGH",
    call,
    score: 100,
    identity: 0.98,
    shm: 0.02,
    alternatives,
    abundance: 1,
  };
}

test("reference graph collapses sequence-identical labels and retains bounded nucleotide neighbours", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  assert.equal(graph.nodes.length, 4);
  assert.equal(graph.exactDuplicateLabels, 1);
  assert.equal(graph.callToNode.get("IGHV1*01"), graph.callToNode.get("IGHV1*01_duplicate_label"));
  const first = graph.callToNode.get("IGHV1*01")!;
  const second = graph.callToNode.get("IGHV1*02")!;
  const third = graph.callToNode.get("IGHV1*03")!;
  assert.equal(graph.neighbours[first].find((value) => value.index === second)?.distance, 1);
  assert.equal(graph.neighbours[first].find((value) => value.index === third)?.distance, 2);
  assert.equal(graph.neighbours[first].some((value) => value.index === graph.callToNode.get("IGKV1*01")), false);
  assert.equal(boundedReferenceDistance("AAAA", "AATA", 1), 1);
  assert.equal(boundedReferenceDistance("AAAA", "AATT", 1), null);
  assert.equal(boundedReferenceDistance("AAAA", "AAA", 1), 1);
});

test("one-SNP neighbour leakage increases with SHM but retains an explicit floor and cap", () => {
  const baseline = adaptiveNeighbourOdds(0, options);
  const lowShm = adaptiveNeighbourOdds(0.03, options);
  const highShm = adaptiveNeighbourOdds(0.12, options);
  const capped = adaptiveNeighbourOdds(0.9, options);
  assert.equal(baseline, options.baselineNeighbourOdds);
  assert.ok(lowShm > baseline);
  assert.ok(highShm > lowShm);
  assert.ok(capped <= options.maximumNeighbourOdds);
  assert.equal(adaptiveNeighbourOdds(0.9, { ...options, shmLeakageSensitivity: 10 }), options.maximumNeighbourOdds);
  assert.equal(adaptiveNeighbourOdds(0.12, { ...options, shmLeakageSensitivity: 0 }), baseline);
});

test("the interactive error-model diagnostic uses the fitted sparse kernel and excludes the primary from alternative bars", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const selected = graph.callToNode.get("IGHV1*01")!;
  const zeroShm = inspectReferenceEvidenceKernel(graph, selected, 0, options)!;
  const highShm = inspectReferenceEvidenceKernel(graph, selected, 0.12, options)!;
  const second = graph.callToNode.get("IGHV1*02")!;
  const zeroAlternative = zeroShm.alternatives.find((value) => value.nodeIndex === second)!;
  const highAlternative = highShm.alternatives.find((value) => value.nodeIndex === second)!;
  assert.ok(highAlternative.probability > zeroAlternative.probability);
  assert.equal(highShm.alternatives.some((value) => value.nodeIndex === selected), false);
  assert.ok(Math.abs(highShm.primaryProbability + highShm.alternativeProbability - 1) < 1e-12);

  const alignment = alignReferenceKernelInspection(highShm);
  assert.equal(new Set(alignment.rows.map((value) => value.sequence.length)).size, 1);
  for (let column = 0; column < alignment.columns; column += 1) {
    assert.equal(alignment.rows.every((value) => value.sequence[column] === "-"), false);
  }
});

test("hard assignment projection counts MAP reassignments and exposes vanished alleles", () => {
  const result: SegmentRefinementResult = {
    segment: "V",
    nodes: [
      { index: 0, segment: "V", locus: "IGH", names: ["IGHV1*01"], sequence: "AAAA" },
      { index: 1, segment: "V", locus: "IGH", names: ["IGHV1*02"], sequence: "AAAT" },
      { index: 2, segment: "V", locus: "IGH", names: ["IGHV1*03"], sequence: "AATT" },
    ],
    mapNode: Int32Array.of(1, 1, 1),
    mapProbability: Float32Array.of(0.95, 0.7, 0.9),
    posteriorEntropy: Float32Array.of(0.1, 0.5, 0.2),
    localTopNode: Int32Array.of(0, 0, 1),
    localTopProbability: Float32Array.of(0.8, 0.8, 0.9),
    modelIndex: Int32Array.of(0, 0, 0),
    assignmentWeight: Float32Array.of(1, 1, 1),
    models: [{
      key: "donor_1\u0000IGH\u0000V", scopeValue: "donor_1", locus: "IGH", segment: "V",
      rows: 3, effectiveRows: 3, nonZeros: 6, databaseNodes: 3, inactivePriorNodes: 1,
      iterations: 5, converged: true, finalMaximumChange: 1e-7,
      alleles: [
        { nodeIndex: 0, names: ["IGHV1*01"], sequenceLength: 4, posteriorMean: 0.25, posteriorSd: 0.1, localEvidenceAssignments: 2, expectedAssignments: 0.5 },
        { nodeIndex: 1, names: ["IGHV1*02"], sequenceLength: 4, posteriorMean: 0.75, posteriorSd: 0.1, localEvidenceAssignments: 1, expectedAssignments: 2.5 },
      ],
    }],
    modeledRows: 3, changedMapRows: 2, skippedRows: 0, matrixNonZeros: 6,
    truncatedRows: 0, exactDuplicateLabels: 0,
  };
  const best = hardAssignmentShiftData(result, 0, "best", 0.8)!;
  const bestFirst = best.rows.find((value) => value.label === "IGHV1*01")!;
  const bestSecond = best.rows.find((value) => value.label === "IGHV1*02")!;
  assert.deepEqual([bestFirst.before, bestFirst.after, bestFirst.vanishes], [2, 0, true]);
  assert.deepEqual([bestSecond.before, bestSecond.after], [1, 3]);
  assert.equal(best.changedAssignments, 2);
  assert.equal(best.vanishedAlleles, 1);
  const allReferences = survivingAlleleReference(result, 0, "best", 0.8, 0)!;
  assert.equal(allReferences.retainedNodes, 3);
  assert.match(allReferences.fasta, />IGHV1\*01 post_reassignment_reads=0\nAAAA/);
  assert.match(allReferences.fasta, />IGHV1\*03 post_reassignment_reads=0\nAATT/);
  const surviving = survivingAlleleReference(result, 0, "best", 0.8, 1)!;
  assert.equal(surviving.retainedNodes, 1);
  assert.equal(surviving.excludedNodes, 2);
  assert.doesNotMatch(surviving.fasta, /IGHV1\*01/);
  assert.match(surviving.fasta, />IGHV1\*02 post_reassignment_reads=3\nAAAT/);

  const confidence = hardAssignmentShiftData(result, 0, "confidence", 0.8)!;
  assert.deepEqual(confidence.rows.map((value) => [value.label, value.before, value.after]), [
    ["IGHV1*02", 1, 2],
    ["IGHV1*01", 2, 1],
  ]);
  assert.equal(confidence.changedAssignments, 1);
  assert.equal(confidence.heldBelowConfidence, 1);
  assert.equal(confidence.vanishedAlleles, 0);
});

test("literal co-optimal calls begin equal and graph neighbours receive geometric non-zero evidence", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const evidence = buildSparseEvidenceRow(row(0, "IGHV1*01,IGHV1*02"), graph, options)!;
  const byNode = new Map(evidence.entries.map((entry) => [entry.node, entry.weight]));
  const first = graph.callToNode.get("IGHV1*01")!;
  const second = graph.callToNode.get("IGHV1*02")!;
  const third = graph.callToNode.get("IGHV1*03")!;
  assert.equal(byNode.get(first), byNode.get(second));
  assert.ok((byNode.get(third) ?? 0) > 0);
  assert.ok((byNode.get(third) ?? 0) < (byNode.get(first) ?? 0));
});

test("sparse Dirichlet updates resolve ambiguous rows using repertoire-wide support", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const rows = [
    ...Array.from({ length: 80 }, (_, ordinal) => row(ordinal, "IGHV1*01")),
    ...Array.from({ length: 20 }, (_, index) => row(80 + index, "IGHV1*01,IGHV1*02")),
  ].map((value) => buildSparseEvidenceRow(value, graph, options)!).filter(Boolean);
  const result = fitSparseAlleleModel(sparseEvidenceMatrix(rows), graph, options, rows.length);
  const first = graph.callToNode.get("IGHV1*01")!;
  for (let ordinal = 80; ordinal < 100; ordinal += 1) {
    assert.equal(result.mapNode[ordinal], first);
    assert.ok(result.mapProbability[ordinal] > 0.98);
  }
  assert.equal(result.models.length, 1);
  assert.equal(result.models[0].converged, true);
  assert.equal(result.matrixNonZeros, rows.reduce((sum, value) => sum + value.entries.length, 0));
});

test("Dirichlet pools are isolated by donor", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const source = [
    ...Array.from({ length: 30 }, (_, ordinal) => row(ordinal, "IGHV1*01", "", "donor_A")),
    row(30, "IGHV1*01,IGHV1*02", "", "donor_A"),
    ...Array.from({ length: 30 }, (_, index) => row(31 + index, "IGHV1*02", "", "donor_B")),
    row(61, "IGHV1*01,IGHV1*02", "", "donor_B"),
  ];
  const rows = source.map((value) => buildSparseEvidenceRow(value, graph, options)!).filter(Boolean);
  const result = fitSparseAlleleModel(sparseEvidenceMatrix(rows), graph, options, source.length);
  assert.equal(result.models.length, 2);
  assert.equal(graph.nodes[result.mapNode[30]].names.includes("IGHV1*01"), true);
  assert.equal(graph.nodes[result.mapNode[61]].names.includes("IGHV1*02"), true);
});

test("Dirichlet normalization retains implicit prior mass for every locus-matched database node", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const noNeighbours = { ...options, neighbourRadius: 0 };
  const evidence = buildSparseEvidenceRow(row(0, "IGHV1*01"), graph, noNeighbours)!;
  const result = fitSparseAlleleModel(sparseEvidenceMatrix([evidence]), graph, noNeighbours, 1);
  const model = result.models[0];
  assert.equal(model.databaseNodes, 3);
  assert.equal(model.inactivePriorNodes, 2);
  assert.equal(model.alleles.length, 1);
  assert.ok(Math.abs(model.alleles[0].posteriorMean - 1.1 / 1.3) < 1e-9);
  const exported = modelSummaryTable({
    version: 1, options: noNeighbours, totalRecords: 1, activeRecords: 1,
    segments: { V: result }, runAt: "2026-08-15T00:00:00.000Z", warnings: [],
  }, "tsv");
  assert.match(exported, /database_nodes\tinactive_prior_nodes/);
  assert.match(exported, /\t3\t2\t/);
});

test("posterior call policies are explicit and resetting restores the immutable AIRR calls", () => {
  const records = [
    { ordinal: 0, vCall: "IGHV1*01,IGHV1*02", jCall: "IGHJ4*01", originalVCall: "IGHV1*01,IGHV1*02", originalJCall: "IGHJ4*01" },
    { ordinal: 1, vCall: "IGHV1*01", jCall: "IGHJ4*01", originalVCall: "IGHV1*01", originalJCall: "IGHJ4*01" },
  ];
  const v = { labels: ["IGHV1*01", "IGHV1*02"], mapNode: Int32Array.of(1, 1), probability: Float32Array.of(0.95, 0.7) };
  const applied = applyCallOverrides(records, v, undefined, "confidence", 0.8);
  assert.equal(applied.changedV, 1);
  assert.equal(applied.policy, "confidence");
  assert.equal(records[0].vCall, "IGHV1*02");
  assert.equal(records[1].vCall, "IGHV1*01");
  const everyMap = applyCallOverrides(records, v, undefined, "best", 0.8);
  assert.equal(everyMap.changedV, 2);
  assert.equal(records[1].vCall, "IGHV1*02");
  applyCallOverrides(records, undefined, undefined, "confidence", 0.8);
  assert.equal(records[0].vCall, "IGHV1*01,IGHV1*02");
});

test("saved sessions round-trip sparse posterior vectors and apply state", () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const evidenceRows = [row(0, "IGHV1*01"), row(1, "IGHV1*01,IGHV1*02")]
    .map((value) => buildSparseEvidenceRow(value, graph, options)!).filter(Boolean);
  const segment = fitSparseAlleleModel(sparseEvidenceMatrix(evidenceRows), graph, options, 2);
  const original = { version: 1 as const, options, totalRecords: 2, activeRecords: 2, segments: { V: segment }, runAt: "2026-08-15T00:00:00.000Z", warnings: [] };
  const saved = saveAlleleRefinement(original, true, "confidence", 0.83);
  const restored = restoreAlleleRefinement(saved);
  assert.deepEqual([...restored.segments.V!.mapNode], [...segment.mapNode]);
  assert.deepEqual([...restored.segments.V!.mapProbability], [...segment.mapProbability]);
  assert.deepEqual([...restored.segments.V!.modelIndex!], [...segment.modelIndex!]);
  assert.deepEqual([...restored.segments.V!.assignmentWeight!], [...segment.assignmentWeight!]);
  assert.equal(saved.applied, true);
  assert.equal(saved.reassignmentPolicy, "confidence");
  assert.equal(saved.applyMinimumPosterior, 0.83);
});

test("posterior sidecar streams the complete sparse responsibility vector", async () => {
  const graph = buildReferenceAlleleGraph(fasta, "V", 2);
  const inputs = [row(0, "IGHV1*01"), row(1, "IGHV1*01,IGHV1*02")];
  const evidenceRows = inputs.map((value) => buildSparseEvidenceRow(value, graph, options)!).filter(Boolean);
  const segment = fitSparseAlleleModel(sparseEvidenceMatrix(evidenceRows), graph, options, 2);
  const result = { version: 1 as const, options, totalRecords: 2, activeRecords: 2, segments: { V: segment }, runAt: "2026-08-15T00:00:00.000Z", warnings: [] };
  const store = {
    scanAirrRows: async (_fields: readonly string[], onBatch: (rows: Array<{ ordinal: number; values: Record<string,string> }>) => void | Promise<void>) => onBatch(inputs.map((input) => ({
      ordinal: input.ordinal,
      values: {
        sequence_id: input.sequenceId, swig_dataset_id: input.datasetId, sample_id: input.sampleId,
        subject_id: input.subjectId, locus: input.locus, duplicate_count: "1", v_call: input.call,
        v_score: String(input.score), v_identity: String(input.identity), v_alternatives: input.alternatives,
      },
    }))),
  } as unknown as AirrResultStore;
  let output = "";
  await writeRefinementSidecar(store, result, "confidence", 0.8, "tsv", async (part) => { output += typeof part === "string" ? part : ""; });
  const lines = output.trim().split("\n");
  const headers = lines[0].split("\t");
  const sequence = headers.indexOf("sequence_id");
  const posterior = headers.indexOf("repertoire_posterior");
  const candidate = headers.indexOf("candidate_call");
  assert.ok(headers.includes("reassignment_policy"));
  assert.ok(headers.includes("map_selected_by_policy"));
  const rows = lines.slice(1).map((line) => line.split("\t"));
  for (const id of ["r0", "r1"]) {
    const selected = rows.filter((values) => values[sequence] === id);
    const total = selected.reduce((sum, values) => sum + Number(values[posterior]), 0);
    assert.ok(Math.abs(total - 1) < 1e-6);
  }
  assert.ok(new Set(rows.filter((values) => values[sequence] === "r1").map((values) => values[candidate])).size >= 2);
});
