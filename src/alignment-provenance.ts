import { parseFasta, type FastaRecord } from "./post-analysis-core.ts";

export interface AlignmentInspection {
  fasta: string;
  records: FastaRecord[];
  rows: number;
  columns: number;
  fingerprint: string;
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

/**
 * Manual correction may move or add gap columns, but it must not silently add,
 * remove, truncate, rename, or mutate biological sequences.
 */
export function validateCorrectedAlignment(currentText: string, correctedText: string): AlignmentInspection {
  const current = inspectAlignment(currentText);
  const corrected = inspectAlignment(correctedText);
  const currentByName = new Map(current.records.map((record) => [record.name, record.sequence]));
  const correctedNames = new Set(corrected.records.map((record) => record.name));
  const missing = current.records.filter((record) => !correctedNames.has(record.name)).map((record) => record.name);
  const added = corrected.records.filter((record) => !currentByName.has(record.name)).map((record) => record.name);
  if (missing.length || added.length) {
    const details = [
      missing.length ? `missing: ${missing.slice(0, 4).join(", ")}${missing.length > 4 ? "…" : ""}` : "",
      added.length ? `unexpected: ${added.slice(0, 4).join(", ")}${added.length > 4 ? "…" : ""}` : "",
    ].filter(Boolean).join("; ");
    throw new Error(`The corrected alignment must contain exactly the original rows (${details}). A selected Alivibe fragment is not a complete alignment.`);
  }
  for (const record of corrected.records) {
    const original = currentByName.get(record.name)!;
    if (ungapped(record.sequence) !== ungapped(original)) {
      throw new Error(`The ungapped sequence for ${record.name} changed or was truncated. Alignment correction may change gaps, but not nucleotides or row contents.`);
    }
  }
  return corrected;
}
