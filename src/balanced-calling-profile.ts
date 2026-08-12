const decoder = new TextDecoder();
const encoder = new TextEncoder();

export interface BalancedDFilterResult {
  body: Uint8Array<ArrayBuffer>;
  suppressedSequenceIds: Set<string>;
}

function longestExactRun(query: string, germline: string): number {
  let current = 0;
  let longest = 0;
  for (let index = 0; index < Math.min(query.length, germline.length); index += 1) {
    if (query[index] !== "-" && query[index] === germline[index]) {
      current += 1;
      longest = Math.max(longest, current);
    } else current = 0;
  }
  return longest;
}

const complement: Record<string, string> = {
  A: "T", C: "G", G: "C", T: "A", U: "A",
  R: "Y", Y: "R", M: "K", K: "M", S: "S", W: "W",
  B: "V", V: "B", D: "H", H: "D", N: "N", "-": "-",
};

function reverseComplement(sequence: string): string {
  let result = "";
  for (let index = sequence.length - 1; index >= 0; index -= 1) {
    result += complement[sequence[index].toUpperCase()] ?? "N";
  }
  return result;
}

const aminoAcids =
  "FFLLSSSSYY**CC*W" +
  "LLLLPPPPHHQQRRRR" +
  "IIIMTTTTNNKKSSRR" +
  "VVVVAAAADDEEGGGG";
const baseIndex: Record<string, number> = { T: 0, C: 1, A: 2, G: 3 };

function translateAlignment(sequence: string, frame: number): string {
  let protein = "";
  for (let index = frame; index + 2 < sequence.length; index += 3) {
    const codon = sequence.slice(index, index + 3).toUpperCase();
    if (codon.includes("-")) {
      protein += "-";
      continue;
    }
    const a = baseIndex[codon[0]];
    const b = baseIndex[codon[1]];
    const c = baseIndex[codon[2]];
    protein += a === undefined || b === undefined || c === undefined
      ? "X"
      : aminoAcids[a * 16 + b * 4 + c];
  }
  return protein;
}

function indexColumns(header: string): Map<string, number> {
  return new Map(header.split("\t").map((name, index) => [name, index]));
}

function field(fields: string[], columns: Map<string, number>, name: string): string {
  const index = columns.get(name);
  return index === undefined ? "" : fields[index] ?? "";
}

function setField(fields: string[], columns: Map<string, number>, name: string, value = "") {
  const index = columns.get(name);
  if (index !== undefined) fields[index] = value;
}

function rebuildCompositeAlignment(fields: string[], columns: Map<string, number>) {
  const sequence = field(fields, columns, "sequence");
  const oriented = field(fields, columns, "rev_comp") === "T" ? reverseComplement(sequence) : sequence;
  const vEnd = Number(field(fields, columns, "v_sequence_end"));
  const jStart = Number(field(fields, columns, "j_sequence_start"));
  const gap = Number.isFinite(vEnd) && Number.isFinite(jStart) && jStart > vEnd
    ? oriented.slice(vEnd, jStart - 1)
    : "";
  const query = field(fields, columns, "v_sequence_alignment") + gap +
    field(fields, columns, "j_sequence_alignment");
  const germline = field(fields, columns, "v_germline_alignment") + "N".repeat(gap.length) +
    field(fields, columns, "j_germline_alignment");
  setField(fields, columns, "sequence_alignment", query);
  setField(fields, columns, "germline_alignment", germline);
  setField(fields, columns, "np1", gap);
  setField(fields, columns, "np1_length", String(Math.max(0, jStart - 1 - vEnd)));
  setField(fields, columns, "np2");
  setField(fields, columns, "np2_length");

  const vStart = Number(field(fields, columns, "v_sequence_start"));
  const cdr3Start = Number(field(fields, columns, "cdr3_start"));
  let cysColumn = -1;
  if (Number.isFinite(vStart) && Number.isFinite(cdr3Start)) {
    const cysPosition = cdr3Start - 4;
    let queryPosition = vStart - 1;
    for (let column = 0; column < query.length; column += 1) {
      if (query[column] === "-") continue;
      if (queryPosition === cysPosition) {
        cysColumn = column;
        break;
      }
      queryPosition += 1;
    }
  }
  if (cysColumn >= 0) {
    const frame = cysColumn % 3;
    setField(fields, columns, "sequence_alignment_aa", translateAlignment(query, frame));
    setField(fields, columns, "germline_alignment_aa", translateAlignment(germline, frame));
  } else {
    setField(fields, columns, "sequence_alignment_aa");
    setField(fields, columns, "germline_alignment_aa");
  }
}

/**
 * Apply the calibrated IgBLAST-balanced evidence rule to an AIRR batch.
 *
 * Only D hits whose strongest alignment has exactly a five-nucleotide exact
 * run are candidates. Such a hit is removed when
 * `j_sequence_start - v_sequence_end <= 11`; all longer exact runs and all
 * five-base hits in longer V-J spans are retained. D-dependent AIRR fields
 * and the stitched V-J alignment are rebuilt so the row stays internally
 * consistent.
 */
export function applyBalancedDFilter(header: string, body: Uint8Array<ArrayBuffer>): BalancedDFilterResult {
  const columns = indexColumns(header);
  const required = [
    "sequence_id", "d_call", "d_sequence_alignment", "d_germline_alignment",
    "v_sequence_end", "j_sequence_start",
  ];
  if (required.some((name) => !columns.has(name))) {
    throw new Error("The AIRR table is missing fields required by the IgBLAST-balanced D rule.");
  }
  const text = decoder.decode(body);
  const hadTrailingNewline = text.endsWith("\n");
  const lines = text.split("\n");
  if (hadTrailingNewline) lines.pop();
  const suppressedSequenceIds = new Set<string>();
  const dFields = [
    "d_call", "d_score", "d_identity", "d_cigar",
    "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
    "d_sequence_alignment", "d_germline_alignment", "d_frame", "d_alternatives",
  ];
  const transformed = lines.map((line) => {
    if (!line) return line;
    const fields = line.replace(/\r$/, "").split("\t");
    if (!field(fields, columns, "d_call")) return line.replace(/\r$/, "");
    const exact = longestExactRun(
      field(fields, columns, "d_sequence_alignment"),
      field(fields, columns, "d_germline_alignment"),
    );
    const vjSpan = Math.max(0,
      Number(field(fields, columns, "j_sequence_start")) -
      Number(field(fields, columns, "v_sequence_end")),
    );
    if (exact !== 5 || vjSpan > 11) return line.replace(/\r$/, "");
    suppressedSequenceIds.add(field(fields, columns, "sequence_id"));
    dFields.forEach((name) => setField(fields, columns, name));
    rebuildCompositeAlignment(fields, columns);
    return fields.join("\t");
  });
  if (!suppressedSequenceIds.size) return { body, suppressedSequenceIds };
  return {
    body: encoder.encode(transformed.join("\n") + (hadTrailingNewline ? "\n" : "")),
    suppressedSequenceIds,
  };
}

export function reconcileBalancedDoubleD(
  header: string,
  body: Uint8Array<ArrayBuffer>,
  suppressedSequenceIds: ReadonlySet<string>,
): Uint8Array<ArrayBuffer> {
  if (!suppressedSequenceIds.size) return body;
  const columns = indexColumns(header);
  const idIndex = columns.get("sequence_id");
  const standardDIndex = columns.get("standard_d_call");
  if (idIndex === undefined || standardDIndex === undefined) return body;
  const text = decoder.decode(body);
  const hadTrailingNewline = text.endsWith("\n");
  const lines = text.split("\n");
  if (hadTrailingNewline) lines.pop();
  let changed = false;
  const transformed = lines.map((line) => {
    const fields = line.replace(/\r$/, "").split("\t");
    if (!suppressedSequenceIds.has(fields[idIndex] ?? "")) return line.replace(/\r$/, "");
    fields[standardDIndex] = "";
    changed = true;
    return fields.join("\t");
  });
  return changed
    ? encoder.encode(transformed.join("\n") + (hadTrailingNewline ? "\n" : ""))
    : body;
}
