import { parseFasta } from "./post-analysis-core.ts";
import { biologicalFrameOffset } from "./alignment-model.ts";
import type { AirrDetailRow } from "./result-store.ts";

export const GERMLINE_OUTGROUP = "__germline_N_masked__";

function safeName(value: string, fallback: string): string {
  return (value || fallback).replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
}

function normalizedAlignment(value: string, fallback = ""): string {
  return (value || fallback).toUpperCase().replaceAll(".", "-").replaceAll("U", "T").replace(/[^ACGTN-]/g, "N");
}

export interface AnchoredLineageRow {
  row: AirrDetailRow;
  name: string;
  sequence: string;
  germline: string;
  frame: number;
}

export interface InferredLineageGermline {
  /** V/D/J-aware root template; unresolved recombination columns remain N. */
  template: string;
  /** Template with unresolved N columns filled only by an unweighted member consensus. */
  uca: string;
  /** Supported endpoint-trimmed forms used for lineage-neighbour comparison. */
  trimmedTemplate: string;
  trimmedUca: string;
  rowsUsed: number;
  columns: number;
  startColumn: number;
  endColumn: number;
  knownColumns: number;
  inferredColumns: number;
  conflictingColumns: number;
  minimumEndpointSupport: number;
}

/**
 * Put AIRR pairwise alignments onto a common V-reference origin. The AIRR
 * `v_germline_start` coordinate supplies left padding; all remaining columns
 * retain the caller's insertion/gap structure. This is deliberately a quick
 * reference-anchored view rather than a de-novo MSA.
 */
export function referenceAnchoredLineageRows(rows: AirrDetailRow[]): AnchoredLineageRow[] {
  if (!rows.length) throw new Error("The selected lineage has no retrievable AIRR records.");
  const provisional = rows.map((row, index) => {
    const rawSequence = normalizedAlignment(row.values.sequence_alignment, row.values.sequence);
    const rawGermline = normalizedAlignment(row.values.germline_alignment, rawSequence.replace(/[ACGT]/g, "N"));
    const localColumns = Math.max(rawSequence.length, rawGermline.length);
    const reportedStart = Math.floor(Number(row.values.v_germline_start));
    const left = Number.isFinite(reportedStart) && reportedStart >= 1 && reportedStart <= 10_000 ? reportedStart - 1 : 0;
    return {
      row,
      name: `${safeName(row.values.sequence_id, `sequence_${index + 1}`)}__${row.record.ordinal + 1}`,
      sequence: "-".repeat(left) + rawSequence.padEnd(localColumns, "-"),
      germline: "-".repeat(left) + rawGermline.padEnd(localColumns, "-"),
      frame: biologicalFrameOffset(Number(row.values.v_sequence_start) || 1, Number(row.values.sequence_frame) || 1),
    };
  });
  const columns = Math.max(...provisional.map((record) => Math.max(record.sequence.length, record.germline.length)));
  return provisional.map((record) => ({
    ...record,
    sequence: record.sequence.padEnd(columns, "-"),
    germline: record.germline.padEnd(columns, "-"),
  }));
}

function plurality(counts: readonly number[]): { base: string; count: number; tied: boolean } {
  const bases = ["A", "C", "G", "T"];
  let maximum = 0;
  let maximumIndex = 0;
  let tied = false;
  for (let index = 0; index < counts.length; index += 1) {
    if (counts[index] > maximum) {
      maximum = counts[index];
      maximumIndex = index;
      tied = false;
    } else if (counts[index] === maximum && maximum > 0) tied = true;
  }
  return { base: bases[maximumIndex], count: maximum, tied };
}

const BASE_INDEX: Record<string, number | undefined> = { A: 0, C: 1, G: 2, T: 3 };

/**
 * Infer one lineage germline without selecting a closest read. Known V/D/J
 * bases are voted from each member's AIRR germline alignment. A known base is
 * retained only with >=80% agreement among informative members. At columns
 * that all germline alignments mark N, an unweighted >=60% query consensus is
 * exposed separately as `uca`; the tree template deliberately remains N so
 * parsimony can treat it as unknown. Endpoint trimming requires coverage from
 * at least 20% of the loaded unique representatives.
 */
export function inferLineageGermline(rows: AirrDetailRow[]): InferredLineageGermline {
  const anchored = referenceAnchoredLineageRows(rows);
  const minimumEndpointSupport = Math.max(1, Math.ceil(anchored.length * 0.2));
  const template: string[] = [];
  const uca: string[] = [];
  const supported: boolean[] = [];
  let knownColumns = 0;
  let inferredColumns = 0;
  let conflictingColumns = 0;

  for (let column = 0; column < anchored[0].sequence.length; column += 1) {
    const germlineCounts = [0, 0, 0, 0];
    const queryCounts = [0, 0, 0, 0];
    let coverage = 0;
    for (const record of anchored) {
      const germlineBase = record.germline[column] ?? "-";
      const queryBase = record.sequence[column] ?? "-";
      const germlineIndex = BASE_INDEX[germlineBase];
      const queryIndex = BASE_INDEX[queryBase];
      if (germlineIndex !== undefined) germlineCounts[germlineIndex] += 1;
      if (queryIndex !== undefined) queryCounts[queryIndex] += 1;
      if (germlineBase !== "-" || queryBase !== "-") coverage += 1;
    }
    supported.push(coverage >= minimumEndpointSupport);
    const germlineVotes = germlineCounts.reduce((sum, value) => sum + value, 0);
    const queryVotes = queryCounts.reduce((sum, value) => sum + value, 0);
    const known = plurality(germlineCounts);
    const observed = plurality(queryCounts);
    if (germlineVotes && !known.tied && known.count / germlineVotes >= 0.8) {
      template.push(known.base);
      uca.push(known.base);
      knownColumns += 1;
    } else {
      if (germlineVotes) conflictingColumns += 1;
      template.push(coverage ? "N" : "-");
      if (queryVotes && !observed.tied && observed.count / queryVotes >= 0.6) {
        uca.push(observed.base);
        inferredColumns += 1;
      } else uca.push(coverage ? "N" : "-");
    }
  }
  let start = supported.findIndex(Boolean);
  if (start < 0) start = 0;
  let end = supported.length;
  while (end > start && !supported[end - 1]) end -= 1;
  return {
    template: template.join(""),
    uca: uca.join(""),
    trimmedTemplate: template.slice(start, end).join(""),
    trimmedUca: uca.slice(start, end).join(""),
    rowsUsed: anchored.length,
    columns: template.length,
    startColumn: start,
    endColumn: end,
    knownColumns,
    inferredColumns,
    conflictingColumns,
    minimumEndpointSupport,
  };
}

export function lineageInputFasta(rows: AirrDetailRow[]): { fasta: string; frames: number[]; germline: string } {
  const anchored = referenceAnchoredLineageRows(rows);
  const inferred = inferLineageGermline(rows);
  const germline = inferred.template.padEnd(anchored[0].sequence.length, "-");
  const records = anchored.map((record) => ({
    name: record.name,
    sequence: record.sequence,
    frame: record.frame,
  }));
  records.push({ name: GERMLINE_OUTGROUP, sequence: germline, frame: records[0]?.frame ?? 0 });
  return {
    fasta: records.map((record) => `>${record.name}\n${record.sequence}`).join("\n") + "\n",
    frames: records.map((record) => record.frame),
    germline,
  };
}

export function quickAirrAlignment(rows: AirrDetailRow[]): string {
  const input = lineageInputFasta(rows);
  const records = parseFasta(input.fasta, true);
  const maximum = Math.max(...records.map((record) => record.sequence.length));
  return records.map((record) => `>${record.name}\n${record.sequence.padEnd(maximum, "-")}`).join("\n") + "\n";
}
