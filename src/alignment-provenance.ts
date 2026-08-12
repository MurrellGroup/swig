import { parseFasta, type FastaRecord } from "./post-analysis-core.ts";

export interface AlignmentInspection {
  fasta: string;
  records: FastaRecord[];
  rows: number;
  columns: number;
  fingerprint: string;
}

export interface AlignmentCorrectionInspection extends AlignmentInspection {
  removedRows: string[];
  removedNucleotides: number;
}

function assertFastaAlphabet(text: string): void {
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith(">")) continue;
    const unsupported = line.replace(/\s/g, "").match(/[^ACGTUNRYKMSWBDHV.\-]/i)?.[0];
    if (unsupported) throw new Error(`The alignment contains the unsupported nucleotide character ${JSON.stringify(unsupported)}.`);
  }
}

function fnv1a64(value: string): string {
  let hash = 0xcbf29ce484222325n;
  const mask = 0xffffffffffffffffn;
  for (const byte of new TextEncoder().encode(value)) {
    hash ^= BigInt(byte);
    hash = (hash * 0x100000001b3n) & mask;
  }
  return hash.toString(16).padStart(16, "0");
}

export function inspectAlignment(text: string, minimumRows = 2): AlignmentInspection {
  assertFastaAlphabet(text);
  const normalizedText = text.split(/\r?\n/).map((line) => line.startsWith(">") ? line : line.replace(/U/gi, "T")).join("\n");
  const records = parseFasta(normalizedText, true);
  if (records.length < minimumRows) throw new Error(`The alignment must contain at least ${minimumRows} sequences.`);
  const columns = records[0]?.sequence.length ?? 0;
  if (!columns || records.some((record) => record.sequence.length !== columns)) {
    throw new Error("Every alignment record must have the same non-zero aligned length.");
  }
  const names = new Set<string>();
  for (const record of records) {
    if (!record.name) throw new Error("Every alignment record must have a FASTA identifier.");
    if (names.has(record.name)) throw new Error(`The alignment contains the duplicate identifier ${record.name}.`);
    names.add(record.name);
  }
  const fasta = records.map((record) => `>${record.name}\n${record.sequence}`).join("\n") + "\n";
  return {
    fasta,
    records,
    rows: records.length,
    columns,
    fingerprint: fnv1a64(fasta),
  };
}

function ungapped(sequence: string): string {
  return sequence.replaceAll("-", "").replaceAll("U", "T");
}

function isSubsequence(candidate: string, original: string): boolean {
  let offset = 0;
  for (const base of candidate) {
    offset = original.indexOf(base, offset);
    if (offset < 0) return false;
    offset += 1;
  }
  return true;
}

/**
 * Manual correction may move gaps, delete biological rows, and remove bases or
 * alignment columns. It may not add/rename rows or introduce/substitute bases.
 * The N-masked germline is retained because downstream rooting depends on it.
 */
export function validateCorrectedAlignment(currentText: string, correctedText: string): AlignmentCorrectionInspection {
  const current = inspectAlignment(currentText);
  const corrected = inspectAlignment(correctedText);
  const currentByName = new Map(current.records.map((record) => [record.name, record.sequence]));
  const correctedNames = new Set(corrected.records.map((record) => record.name));
  const missing = current.records.filter((record) => !correctedNames.has(record.name)).map((record) => record.name);
  const added = corrected.records.filter((record) => !currentByName.has(record.name)).map((record) => record.name);
  if (added.length) throw new Error(`The corrected alignment contains unexpected or renamed rows: ${added.slice(0, 4).join(", ")}${added.length > 4 ? "…" : ""}. New biological sequences cannot be introduced during alignment correction.`);
  if (missing.includes("__germline_N_masked__")) throw new Error("The corrected alignment must retain __germline_N_masked__ so the lineage tree can be rooted reproducibly.");
  let removedNucleotides = 0;
  for (const record of corrected.records) {
    const original = currentByName.get(record.name)!;
    const originalUngapped = ungapped(original);
    const correctedUngapped = ungapped(record.sequence);
    if (!isSubsequence(correctedUngapped, originalUngapped)) throw new Error(`The ungapped sequence for ${record.name} contains a substitution, inserted base, or changed order. Row/base deletion and gap editing are allowed; new nucleotide content is not.`);
    removedNucleotides += originalUngapped.length - correctedUngapped.length;
  }
  return { ...corrected, removedRows: missing, removedNucleotides };
}
