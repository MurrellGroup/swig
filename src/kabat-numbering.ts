import { translateAlignedNucleotides } from "./lineage-phylogeny.ts";
import type { FastaRecord } from "./post-analysis-core.ts";

const GERMLINE_OUTGROUP = "__germline_N_masked__";

export interface KabatColumnMap {
  labels: string[];
  chain: "H" | "K" | "L";
  confidence: number;
  numberedColumns: number;
  contributingSequences: number;
  warnings: string[];
  partialCodonRecords: number;
  stopCodons: number;
}

interface KabatNumberingResult {
  chain: string | null;
  confidence: number | null;
  numbering: Map<string, string> | null;
  query_start: number | null;
  error: string | null;
}

export interface KabatNumberer {
  number: (sequence: string) => KabatNumberingResult;
}

/**
 * Number several aligned members and vote per alignment column. This avoids
 * attaching the entire alignment's coordinates to one sequence-specific gap
 * pattern while preserving Kabat insertion labels (for example 35A).
 */
export function inferKabatColumnsWithNumberer(records: FastaRecord[], annotator: KabatNumberer, maximumContributors = 24, frameOffset = 0): KabatColumnMap {
  if (!records.length) throw new Error("Kabat numbering requires a lineage alignment.");
  const offset = frameOffset === 1 || frameOffset === 2 ? frameOffset : 0;
  const partialCodonRecords = records.filter((record) => Math.max(0, record.sequence.length - offset) % 3 !== 0).length;
  const translatedRecords = records.map((record) => translateAlignedNucleotides(record.sequence, offset));
  const stopCodons = translatedRecords.reduce((total, sequence) => total + [...sequence].filter((residue) => residue === "*").length, 0);
  const warnings: string[] = [];
  if (partialCodonRecords || stopCodons) {
    const details = [
      partialCodonRecords ? `${partialCodonRecords} sequence${partialCodonRecords === 1 ? " has" : "s have"} a terminal partial codon` : "",
      stopCodons ? `${stopCodons} stop codon${stopCodons === 1 ? "" : "s"}` : "",
    ].filter(Boolean).join("; ");
    warnings.push(`Partial or stop codons detected (${details}). Kabat numbering may be unreliable near those columns.`);
  }
  const width = Math.max(...records.map((record) => Math.ceil(Math.max(0, record.sequence.length - offset) / 3)));
  const votes = Array.from({ length: width }, () => new Map<string, number>());
  const chainVotes = new Map<"H" | "K" | "L", number>();
  let confidence = 0;
  let contributingSequences = 0;
  for (let recordIndex = 0; recordIndex < records.length; recordIndex += 1) {
    const record = records[recordIndex];
    if (record.name === GERMLINE_OUTGROUP || contributingSequences >= maximumContributors) continue;
    const aligned = translatedRecords[recordIndex];
    const translatedLength = Math.max(0, record.sequence.length - offset);
    const partialColumn = translatedLength % 3 ? Math.floor(translatedLength / 3) : -1;
    const columns: number[] = [];
    const sequence = [...aligned].filter((residue, column) => {
      // A terminal one- or two-base overhang is displayed as X but must not be
      // presented to the domain numberer as a complete amino acid.
      if (residue === "-" || column === partialColumn) return false;
      columns.push(column);
      return true;
    }).join("").replace(/\*/g, "X");
    if (sequence.length < 55) continue;
    const result = annotator.number(sequence);
    if (result.error || !result.numbering || result.query_start === null || !result.chain || !["H", "K", "L"].includes(result.chain)) continue;
    const entries = [...result.numbering.keys()];
    for (let index = 0; index < entries.length; index += 1) {
      const column = columns[result.query_start + index];
      if (column === undefined) continue;
      const label = entries[index];
      votes[column].set(label, (votes[column].get(label) ?? 0) + 1);
    }
    chainVotes.set(result.chain as "H" | "K" | "L", (chainVotes.get(result.chain as "H" | "K" | "L") ?? 0) + 1);
    confidence += result.confidence ?? 0;
    contributingSequences += 1;
  }
  if (!contributingSequences) throw new Error("No aligned member could be assigned a confident IGH, IGK, or IGL Kabat numbering.");
  const chain = [...chainVotes].sort((left, right) => right[1] - left[1])[0][0];
  const labels = votes.map((columnVotes) => [...columnVotes].sort((left, right) => right[1] - left[1] || left[0].localeCompare(right[0]))[0]?.[0] ?? "");
  return {
    labels,
    chain,
    confidence: confidence / contributingSequences,
    numberedColumns: labels.filter(Boolean).length,
    contributingSequences,
    warnings,
    partialCodonRecords,
    stopCodons,
  };
}

export async function inferKabatColumns(records: FastaRecord[], maximumContributors = 24, frameOffset = 0): Promise<KabatColumnMap> {
  const { Annotator, initializeImmunum } = await import("virtual:immunum-browser");
  await initializeImmunum();
  const annotator = new Annotator(["H", "K", "L"], "kabat", 0.25);
  try {
    return inferKabatColumnsWithNumberer(records, annotator, maximumContributors, frameOffset);
  } finally {
    annotator.free();
  }
}
