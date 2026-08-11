import { parseFasta } from "./post-analysis-core";
import { biologicalFrameOffset } from "./alignment-model";
import type { AirrDetailRow } from "./result-store";

export const GERMLINE_OUTGROUP = "__germline_N_masked__";

function safeName(value: string, fallback: string): string {
  return (value || fallback).replace(/[\s():;,]/g, "_");
}

function nMaskedGermline(row: AirrDetailRow): string {
  const existing = row.values.germline_alignment.toUpperCase().replaceAll(".", "-").replace(/[^ACGTN-]/g, "N");
  if (existing) return existing;
  const sequence = row.values.sequence_alignment || row.values.sequence;
  return sequence.toUpperCase().replace(/[ACGT]/g, "N").replace(/[^N-]/g, "N");
}

export function lineageInputFasta(rows: AirrDetailRow[]): { fasta: string; frames: number[]; germline: string } {
  if (!rows.length) throw new Error("The selected lineage has no retrievable AIRR records.");
  const germline = nMaskedGermline(rows[0]);
  const records = rows.map((row, index) => ({
    name: `${safeName(row.values.sequence_id, `sequence_${index + 1}`)}__${row.record.ordinal + 1}`,
    sequence: (row.values.sequence_alignment || row.values.sequence).toUpperCase().replace(/[^ACGTN-]/g, "N"),
    frame: biologicalFrameOffset(Number(row.values.v_sequence_start) || 1, Number(row.values.sequence_frame) || 1),
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
