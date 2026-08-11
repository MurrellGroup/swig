import { inspectAlignment } from "./alignment-provenance.ts";

export interface PreparedFastTreeInput {
  inputFasta: string;
  alignmentFasta: string;
  names: string[];
  flags: string;
  command: string;
  rows: number;
  columns: number;
  fingerprint: string;
}

export function prepareFastTreeInput(alignedFasta: string, model: "gtr" | "jc" = "gtr", fast = false): PreparedFastTreeInput {
  const inspected = inspectAlignment(alignedFasta, 3);
  const inputFasta = inspected.records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
  const flags = ["-nt", model === "gtr" ? "-gtr" : "", fast ? "-fastest" : ""].filter(Boolean).join(" ");
  return {
    inputFasta,
    alignmentFasta: inspected.fasta,
    names: inspected.records.map((record) => record.name),
    flags,
    command: `fasttree ${flags} input.fasta`,
    rows: inspected.rows,
    columns: inspected.columns,
    fingerprint: inspected.fingerprint,
  };
}
