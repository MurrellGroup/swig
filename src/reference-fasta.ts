export interface ReferenceFastaRecord {
  /** Complete FASTA header without the leading `>`. */
  header: string;
  /** Exact identifier used by AIRR calls: the first whitespace-delimited token. */
  name: string;
  /** Sequence with line wrapping and whitespace removed; gaps and IUPAC symbols are retained. */
  sequence: string;
}

export interface FilteredReferenceFasta {
  fasta: string;
  total: number;
  retained: number;
  excluded: number;
  excludedNames: string[];
  unmatchedExclusions: string[];
}

export function parseReferenceFasta(text: string): ReferenceFastaRecord[] {
  const records: ReferenceFastaRecord[] = [];
  let header = "";
  let sequence: string[] = [];
  const commit = () => {
    if (!header) return;
    const name = header.trim().split(/\s+/, 1)[0] ?? "";
    const joined = sequence.join("").replace(/\s/g, "").toUpperCase();
    if (name && joined) records.push({ header: header.trim(), name, sequence: joined });
  };
  for (const rawLine of text.split(/\r?\n/)) {
    if (rawLine.startsWith(">")) {
      commit();
      header = rawLine.slice(1);
      sequence = [];
    } else if (header) sequence.push(rawLine);
  }
  commit();
  return records;
}

export function serializeReferenceFasta(records: readonly ReferenceFastaRecord[]): string {
  return records.map((record) => `>${record.header}\n${record.sequence}\n`).join("");
}

export function referenceFastaNames(text: string): string[] {
  const seen = new Set<string>();
  const names: string[] = [];
  for (const record of parseReferenceFasta(text)) {
    if (seen.has(record.name)) continue;
    seen.add(record.name);
    names.push(record.name);
  }
  return names;
}

/** Remove exact FASTA identifiers while retaining headers, metadata, and alignment gaps. */
export function filterReferenceFasta(text: string, excluded: Iterable<string>): FilteredReferenceFasta {
  const requested = new Set([...excluded].map((name) => name.trim()).filter(Boolean));
  const records = parseReferenceFasta(text);
  const present = new Set(records.map((record) => record.name));
  const retained = records.filter((record) => !requested.has(record.name));
  const excludedNames = [...requested].filter((name) => present.has(name)).sort((left, right) => left.localeCompare(right, undefined, { numeric: true }));
  return {
    fasta: serializeReferenceFasta(retained),
    total: records.length,
    retained: retained.length,
    excluded: records.length - retained.length,
    excludedNames,
    unmatchedExclusions: [...requested].filter((name) => !present.has(name)).sort((left, right) => left.localeCompare(right, undefined, { numeric: true })),
  };
}
