import { translateAlignedNucleotides } from "./lineage-phylogeny.ts";
import type { FastaRecord } from "./post-analysis-core.ts";

const GERMLINE_OUTGROUP = "__germline_N_masked__";

export interface KabatColumnMap {
  labels: string[];
  chain: "H" | "K" | "L";
  confidence: number;
  numberedColumns: number;
  contributingSequences: number;
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
export function inferKabatColumnsWithNumberer(records: FastaRecord[], annotator: KabatNumberer, maximumContributors = 24): KabatColumnMap {
  if (!records.length) throw new Error("Kabat numbering requires a lineage alignment.");
  if (records.some((record) => record.sequence.length % 3 !== 0)) {
    throw new Error("Kabat numbering requires a codon-column nucleotide alignment. Use the codon-aware aligner or import a codon-preserving correction.");
  }
  const width = Math.max(...records.map((record) => record.sequence.length / 3));
  const votes = Array.from({ length: width }, () => new Map<string, number>());
  const chainVotes = new Map<"H" | "K" | "L", number>();
  let confidence = 0;
  let contributingSequences = 0;
  for (const record of records) {
    if (record.name === GERMLINE_OUTGROUP || contributingSequences >= maximumContributors) continue;
    const aligned = translateAlignedNucleotides(record.sequence);
    const columns: number[] = [];
    const sequence = [...aligned].filter((residue, column) => {
      if (residue === "-") return false;
      columns.push(column);
      return true;
    }).join("");
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
  };
}

export async function inferKabatColumns(records: FastaRecord[], maximumContributors = 24): Promise<KabatColumnMap> {
  const { Annotator, initializeImmunum } = await import("virtual:immunum-browser");
  await initializeImmunum();
  const annotator = new Annotator(["H", "K", "L"], "kabat", 0.25);
  try {
    return inferKabatColumnsWithNumberer(records, annotator, maximumContributors);
  } finally {
    annotator.free();
  }
}
