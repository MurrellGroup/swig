import { tableHeader, tableRow, type TableExportFormat } from "../export-formats.ts";
import type { AirrDetailRow, AirrResultStore } from "../result-store.ts";
import { posteriorMapPassesPolicy } from "./apply.ts";
import type { AlleleReassignmentPolicy, AlleleRefinementResult, RefinementSegment, SegmentRefinementResult } from "./types.ts";
import { buildReferenceAlleleGraph } from "./reference-graph.ts";
import { buildSparseEvidenceRow, parseAlternativeEvidence } from "./evidence.ts";
import { refinementInputFields, toRefinementInputRow } from "./input.ts";
import { digamma } from "./model.ts";

const SEGMENTS: RefinementSegment[] = ["V", "D", "J"];

function assignment(result: SegmentRefinementResult | undefined, ordinal: number) {
  if (!result) return null;
  const node = result.mapNode[ordinal] ?? -1;
  if (node < 0 || !result.nodes[node]) return null;
  const local = result.localTopNode[ordinal] ?? -1;
  return {
    call: result.nodes[node].names.join(","),
    probability: result.mapProbability[ordinal] ?? 0,
    entropy: result.posteriorEntropy[ordinal] ?? 0,
    localCall: local >= 0 && result.nodes[local] ? result.nodes[local].names.join(",") : "",
    localProbability: result.localTopProbability[ordinal] ?? 0,
  };
}

export function refinedCall(
  result: AlleleRefinementResult | null,
  segment: RefinementSegment,
  ordinal: number,
  policy: AlleleReassignmentPolicy,
  minimumPosterior: number,
): string | null {
  const value = assignment(result?.segments[segment], ordinal);
  return value && posteriorMapPassesPolicy(policy, value.probability, minimumPosterior) ? value.call : null;
}

export function refineDetailRows(
  rows: AirrDetailRow[],
  result: AlleleRefinementResult | null,
  policy: AlleleReassignmentPolicy,
  minimumPosterior: number,
  applied: boolean,
): AirrDetailRow[] {
  if (!result || !applied) return rows;
  return rows.map((row) => {
    const v = refinedCall(result, "V", row.record.ordinal, policy, minimumPosterior);
    const d = refinedCall(result, "D", row.record.ordinal, policy, minimumPosterior);
    const j = refinedCall(result, "J", row.record.ordinal, policy, minimumPosterior);
    if (!v && !d && !j) return row;
    const values = { ...row.values };
    if (v) { values.swig_original_v_call = values.v_call ?? ""; values.swig_repertoire_v_call = v; values.v_call = v; }
    if (d) { values.swig_original_d_call = values.d_call ?? ""; values.swig_repertoire_d_call = d; values.d_call = d; }
    if (j) { values.swig_original_j_call = values.j_call ?? ""; values.swig_repertoire_j_call = j; values.j_call = j; }
    return {
      values,
      record: {
        ...row.record,
        vCall: v ?? row.record.vCall,
        dCall: d ?? row.record.dCall,
        jCall: j ?? row.record.jCall,
      },
    };
  });
}

export const REFINEMENT_SIDECAR_FIELDS = [
  "sequence_id", "swig_dataset_id", "sample_id", "subject_id", "locus", "segment",
  "original_call", "candidate_call", "candidate_source", "local_evidence",
  "repertoire_posterior", "pool_posterior_mean", "is_map", "map_probability",
  "posterior_entropy", "reassignment_policy", "minimum_posterior", "map_selected_by_policy",
];

function resultGraph(result: SegmentRefinementResult, radius: number) {
  const fasta = result.nodes.flatMap((node) => node.names.map((name) => `>${name}\n${node.sequence}\n`)).join("");
  const graph = buildReferenceAlleleGraph(fasta, result.segment, radius);
  if (graph.nodes.length !== result.nodes.length || graph.nodes.some((node, index) => node.sequence !== result.nodes[index].sequence)) {
    throw new Error(`Could not reconstruct the saved ${result.segment} reference-neighbour graph for posterior export.`);
  }
  return graph;
}

function expectedLogByModel(result: SegmentRefinementResult, alpha: number) {
  return new Map(result.models.map((model) => {
    const gamma = model.alleles.map((allele) => Math.max(1e-12, alpha + allele.expectedAssignments));
    const total = (model.inactivePriorNodes ?? 0) * alpha + gamma.reduce((sum, value) => sum + value, 0);
    return [model.key, {
      expectedLog: new Map(model.alleles.map((allele, index) => [allele.nodeIndex, digamma(gamma[index]) - digamma(total)] as const)),
      posteriorMean: new Map(model.alleles.map((allele) => [allele.nodeIndex, allele.posteriorMean] as const)),
    }] as const;
  }));
}

export async function writeRefinementSidecar(
  store: AirrResultStore,
  result: AlleleRefinementResult,
  policy: AlleleReassignmentPolicy,
  minimumPosterior: number,
  format: TableExportFormat,
  write: (part: string | Blob | Uint8Array) => Promise<void>,
  includeMask?: Uint8Array,
) {
  const header = tableHeader(REFINEMENT_SIDECAR_FIELDS, format);
  if (header) await write(header);
  for (const segment of SEGMENTS) {
    const segmentResult = result.segments[segment];
    if (!segmentResult) continue;
    const graph = resultGraph(segmentResult, result.options.neighbourRadius);
    const models = expectedLogByModel(segmentResult, Math.max(1e-9, result.options.alphaPerAllele));
    await store.scanAirrRows(refinementInputFields(segment), async (rows) => {
      let body = "";
      for (const row of rows) {
        const input = toRefinementInputRow(row, segment);
        const sparse = buildSparseEvidenceRow(input, graph, result.options);
        if (!sparse) continue;
        const model = models.get(sparse.groupKey);
        if (!model) continue;
        let maximum = Number.NEGATIVE_INFINITY;
        for (const entry of sparse.entries) maximum = Math.max(maximum, Math.log(entry.weight) + (model.expectedLog.get(entry.node) ?? Number.NEGATIVE_INFINITY));
        let normalizer = 0;
        for (const entry of sparse.entries) normalizer += Math.exp(Math.log(entry.weight) + (model.expectedLog.get(entry.node) ?? Number.NEGATIVE_INFINITY) - maximum);
        const reported = new Set(input.call.split(",").map((value) => value.trim()).filter(Boolean).map((call) => graph.callToNode.get(call)).filter((value): value is number => value !== undefined));
        const alternatives = new Set(parseAlternativeEvidence(input.alternatives).map((value) => graph.callToNode.get(value.call)).filter((value): value is number => value !== undefined));
        for (const entry of sparse.entries) {
          const posterior = Math.exp(Math.log(entry.weight) + (model.expectedLog.get(entry.node) ?? Number.NEGATIVE_INFINITY) - maximum) / normalizer;
          body += tableRow(REFINEMENT_SIDECAR_FIELDS, {
            sequence_id: input.sequenceId, swig_dataset_id: input.datasetId, sample_id: input.sampleId,
            subject_id: input.subjectId, locus: input.locus, segment, original_call: input.call,
            candidate_call: graph.nodes[entry.node].names.join(","),
            candidate_source: reported.has(entry.node) ? "reported" : alternatives.has(entry.node) ? "retained_alternative" : "reference_neighbour",
            local_evidence: entry.weight, repertoire_posterior: posterior,
            pool_posterior_mean: model.posteriorMean.get(entry.node) ?? 0,
            is_map: segmentResult.mapNode[row.ordinal] === entry.node,
            map_probability: segmentResult.mapProbability[row.ordinal] ?? 0,
            posterior_entropy: segmentResult.posteriorEntropy[row.ordinal] ?? 0,
            reassignment_policy: policy,
            minimum_posterior: policy === "confidence" ? Math.max(0, Math.min(1, minimumPosterior)) : "",
            map_selected_by_policy: segmentResult.mapNode[row.ordinal] === entry.node
              && posteriorMapPassesPolicy(policy, segmentResult.mapProbability[row.ordinal] ?? 0, minimumPosterior),
          }, format);
        }
      }
      if (body) await write(body);
    }, { batchSize: 250, includeMask });
  }
}

export async function writeRefinedAirr(
  store: AirrResultStore,
  result: AlleleRefinementResult,
  policy: AlleleReassignmentPolicy,
  minimumPosterior: number,
  format: TableExportFormat,
  write: (part: string | Blob | Uint8Array) => Promise<void>,
  includeMask?: Uint8Array,
) {
  const custom = ["swig_reassignment_policy", "swig_reassignment_minimum_posterior", ...SEGMENTS.flatMap((segment) => {
    const lower = segment.toLowerCase();
    return [`swig_original_${lower}_call`, `swig_repertoire_${lower}_call`, `swig_repertoire_${lower}_probability`, `swig_repertoire_${lower}_entropy`];
  })];
  const fields = [...store.airrHeaders, ...custom.filter((field) => !store.airrHeaders.includes(field))];
  const header = tableHeader(fields, format);
  if (header) await write(header);
  await store.scanAirrRows(store.airrHeaders, async (rows) => {
    let body = "";
    for (const row of rows) {
      const values: Record<string, string | number> = { ...row.values };
      values.swig_reassignment_policy = policy;
      values.swig_reassignment_minimum_posterior = policy === "confidence" ? Math.max(0, Math.min(1, minimumPosterior)) : "";
      for (const segment of SEGMENTS) {
        const lower = segment.toLowerCase();
        const value = assignment(result.segments[segment], row.ordinal);
        if (!value) continue;
        values[`swig_original_${lower}_call`] = row.values[`${lower}_call`] ?? "";
        values[`swig_repertoire_${lower}_call`] = value.call;
        values[`swig_repertoire_${lower}_probability`] = value.probability;
        values[`swig_repertoire_${lower}_entropy`] = value.entropy;
        if (posteriorMapPassesPolicy(policy, value.probability, minimumPosterior)) values[`${lower}_call`] = value.call;
      }
      body += tableRow(fields, values, format);
    }
    if (body) await write(body);
  }, { batchSize: 2_000, includeMask });
}

export function modelSummaryTable(result: AlleleRefinementResult, format: TableExportFormat): string {
  const fields = ["pool", "locus", "segment", "reference_labels", "posterior_mean", "posterior_sd", "expected_assignments", "local_evidence_assignments", "rows", "effective_rows", "database_nodes", "inactive_prior_nodes", "iterations", "converged"];
  let output = tableHeader(fields, format);
  for (const segment of SEGMENTS) for (const model of result.segments[segment]?.models ?? []) for (const allele of model.alleles) {
    output += tableRow(fields, {
      pool: model.scopeValue, locus: model.locus, segment: model.segment,
      reference_labels: allele.names.join(","), posterior_mean: allele.posteriorMean,
      posterior_sd: allele.posteriorSd, expected_assignments: allele.expectedAssignments,
      local_evidence_assignments: allele.localEvidenceAssignments, rows: model.rows,
      effective_rows: model.effectiveRows, database_nodes: model.databaseNodes ?? model.alleles.length,
      inactive_prior_nodes: model.inactivePriorNodes ?? 0,
      iterations: model.iterations, converged: model.converged,
    }, format);
  }
  return output;
}
