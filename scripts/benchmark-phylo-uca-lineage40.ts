import { readFileSync } from "node:fs";
import { spawnSync } from "node:child_process";
import { performance } from "node:perf_hooks";

import { defaultPhyloUcaOptions } from "../src/phylo-uca/defaults.ts";
import { inferPhyloUca } from "../src/phylo-uca/inference.ts";
import { parseFasta } from "../src/post-analysis-core.ts";
import type { PhyloUcaOptions, PhyloUcaResult } from "../src/phylo-uca/types.ts";

const alignmentPath = process.argv[2];
const fastTreePath = process.argv[3];
if (!alignmentPath || !fastTreePath) {
  throw new Error("Usage: node --experimental-strip-types scripts/benchmark-phylo-uca-lineage40.ts ALIGNMENT_FASTA FASTTREEDBL");
}

const V_GUIDE = "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAAAGA";
const J_GUIDE = "TGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG";
const GUIDE_NAME = "__germline_N_masked__";

async function kiarva(segment: "IGHV" | "IGHD" | "IGHJ"): Promise<string> {
  const response = await fetch(`https://kiarva.scilifelab.se/api/fasta/genomic?file_name=${segment}`, { headers: { "X-api-key": "kiarvafrontend" } });
  if (!response.ok) throw new Error(`KIARVA ${segment} download failed: HTTP ${response.status}`);
  return response.text();
}

function summarize(result: PhyloUcaResult, wallMs: number) {
  const dOccupancy = new Map<number, number>();
  for (const track of result.hmmAnnotations?.marginalized ?? []) {
    if (track.kind !== "D" || track.call !== "IGHD3-3*01") continue;
    for (const point of track.points) dOccupancy.set(point.alignmentColumn, (dOccupancy.get(point.alignmentColumn) ?? 0) + point.probabilities.reduce((sum, value) => sum + value, 0));
  }
  const supportedColumns = [...dOccupancy].filter(([, mass]) => mass >= 0.1).map(([column]) => column).sort((left, right) => left - right);
  return {
    mode: result.options.search.inferenceMode,
    wallSeconds: wallMs / 1000,
    engineSeconds: result.elapsedMs / 1000,
    bestEdge: result.bestPlacement.edgeId,
    edgeFraction: result.bestPlacement.edgeFraction,
    ucaBranchLength: result.bestPlacement.ucaBranchLength,
    logMarginalLikelihood: result.logMarginalLikelihood,
    placementsUsed: result.placements.filter((point) => point.localPosteriorWeight > 0).length,
    d3_3_at_least_10pct: supportedColumns.length ? [supportedColumns[0], supportedColumns.at(-1)] : [],
    d3_3_columns_346_356: Array.from({ length: 11 }, (_, offset) => ({ column: 346 + offset, occupancy: dOccupancy.get(346 + offset) ?? 0 })),
    mcmc: result.mcmcDiagnostics ? {
      retained: result.mcmcDiagnostics.retainedSamples,
      branchAcceptance: result.mcmcDiagnostics.branchAccepted / Math.max(1, result.mcmcDiagnostics.branchProposals),
      positionAcceptance: result.mcmcDiagnostics.positionAccepted / Math.max(1, result.mcmcDiagnostics.positionProposals),
      globalAcceptance: result.mcmcDiagnostics.globalAccepted / Math.max(1, result.mcmcDiagnostics.globalProposals),
      edgeSwitches: result.mcmcDiagnostics.edgeSwitches,
      branchEss: result.mcmcDiagnostics.branchEffectiveSampleSize,
      logTargetEss: result.mcmcDiagnostics.logTargetEffectiveSampleSize,
      distinctPendantLengths: new Set(result.mcmcDiagnostics.trace.map((point) => point.ucaBranchLength.toPrecision(12))).size,
    } : undefined,
    grid: result.evaluatedUcaBranchLengths,
  };
}

const observedAlignmentFasta = readFileSync(alignmentPath, "utf8");
const observedRecords = parseFasta(observedAlignmentFasta, true);
const columns = observedRecords[0]?.sequence.length ?? 0;
if (!columns || observedRecords.some((record) => record.sequence.length !== columns)) throw new Error("Benchmark alignment is not rectangular.");
const guide = V_GUIDE + "N".repeat(columns - V_GUIDE.length - J_GUIDE.length) + J_GUIDE;
if (guide.length !== columns) throw new Error(`Benchmark guide has ${guide.length} columns, expected ${columns}.`);
const curatedAlignmentFasta = `${observedAlignmentFasta.trim()}\n>${GUIDE_NAME}\n${guide}\n`;
const treeRun = spawnSync(fastTreePath, ["-quiet", "-nt", "-gtr"], { input: observedAlignmentFasta, encoding: "utf8", maxBuffer: 16 * 1024 * 1024 });
if (treeRun.status !== 0 || !treeRun.stdout.trim().endsWith(";")) throw new Error(`FastTreeDbl failed (${treeRun.status}): ${treeRun.stderr}`);
const [V, D, J] = await Promise.all([kiarva("IGHV"), kiarva("IGHD"), kiarva("IGHJ")]);
const rows = observedRecords.map((record, ordinal) => ({
  ordinal,
  sequenceId: record.name,
  locus: "IGH",
  values: {
    locus: "IGH",
    sequence_id: record.name,
    v_call: "IGHV3-30*18,IGHV3-30-5*01",
    j_call: "IGHJ3*02",
    v_sequence_alignment: record.sequence.slice(0, V_GUIDE.length),
    v_germline_alignment: V_GUIDE,
    j_sequence_alignment: record.sequence.slice(-J_GUIDE.length),
    j_germline_alignment: J_GUIDE,
  },
}));

const common = defaultPhyloUcaOptions();
common.hmm = { ...common.hmm, templateMismatchProbability: 0 };
common.search = { ...common.search, fullHmmEdges: 6, screenEdgeGridPoints: 5, maximumUcaBranchLength: 0.3 };
const requestedMode = process.argv[4] as PhyloUcaOptions["search"]["inferenceMode"] | undefined;
const modes: PhyloUcaOptions["search"]["inferenceMode"][] = requestedMode ? [requestedMode] : ["maximum-likelihood", "grid-marginalization", "gibbs-mh"];
const summaries = [];
for (const mode of modes) {
  const options: PhyloUcaOptions = {
    ...common,
    model: { ...common.model, frequencies: [...common.model.frequencies], exchangeabilities: [...common.model.exchangeabilities] },
    candidates: { ...common.candidates },
    hmm: { ...common.hmm, nBaseFrequencies: [...common.hmm.nBaseFrequencies] },
    search: { ...common.search, inferenceMode: mode },
  };
  const started = performance.now();
  const result = await inferPhyloUca({
    curatedAlignmentFasta,
    observedTreeNewick: treeRun.stdout.trim(),
    observedAlignmentFasta,
    retainedColumns: Array.from({ length: columns }, (_, column) => column),
    germlineGuideName: GUIDE_NAME,
    lineageRows: rows,
    references: { V, D, J },
    locus: "IGH",
    lineageLabel: "swig-study-lineage-40",
    alignmentFingerprint: "provided-lineage-40",
    frameOffset: 0,
    options,
  });
  summaries.push(summarize(result, performance.now() - started));
  console.error(`${mode}: ${summaries.at(-1)!.wallSeconds.toFixed(3)} s`);
}
console.log(JSON.stringify({ records: observedRecords.length, columns, fastTreeStderr: treeRun.stderr.trim(), summaries }, null, 2));
