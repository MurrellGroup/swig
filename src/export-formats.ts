export type TableExportFormat = "tsv" | "csv" | "jsonl";
export type AlignmentExportFormat = "fasta" | "clustal" | "phylip" | "stockholm" | "nexus";

export interface NamedSequence {
  name: string;
  sequence: string;
}

function csvCell(value: string | number | boolean | null | undefined): string {
  const text = value == null ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

export function tableExtension(format: TableExportFormat): string {
  return format === "jsonl" ? ".jsonl" : format === "csv" ? ".csv" : ".tsv";
}

export function tableMime(format: TableExportFormat): string {
  return format === "jsonl" ? "application/x-ndjson" : format === "csv" ? "text/csv" : "text/tab-separated-values";
}

export function tableHeader(fields: string[], format: TableExportFormat): string {
  if (format === "jsonl") return "";
  return `${format === "csv" ? fields.map(csvCell).join(",") : fields.join("\t")}\n`;
}

export function tableRow(
  fields: string[],
  values: Record<string, string | number | boolean | null | undefined>,
  format: TableExportFormat,
): string {
  if (format === "jsonl") {
    return `${JSON.stringify(Object.fromEntries(fields.map((field) => [field, values[field] ?? ""]))) }\n`;
  }
  const separator = format === "csv" ? "," : "\t";
  return `${fields.map((field) => format === "csv" ? csvCell(values[field]) : String(values[field] ?? "").replace(/[\t\r\n]/g, " ")).join(separator)}\n`;
}

function normalize(records: NamedSequence[]): NamedSequence[] {
  const result = records.map((record, index) => ({
    name: record.name.trim() || `sequence_${index + 1}`,
    sequence: record.sequence.replace(/\s/g, "").toUpperCase(),
  }));
  if (!result.length) throw new Error("The alignment contains no sequences.");
  const length = result[0].sequence.length;
  if (result.some((record) => record.sequence.length !== length)) {
    throw new Error("Alignment exports require sequences with equal aligned lengths.");
  }
  return result;
}

function safeName(name: string, index: number): string {
  return name.replace(/\s+/g, "_").replace(/[^A-Za-z0-9_.|:+-]/g, "_") || `sequence_${index + 1}`;
}

function nexusName(name: string): string {
  return `'${name.replace(/'/g, "''")}'`;
}

export function alignmentExtension(format: AlignmentExportFormat): string {
  return ({ fasta: ".fasta", clustal: ".aln", phylip: ".phy", stockholm: ".sto", nexus: ".nex" })[format];
}

export function alignmentText(input: NamedSequence[], format: AlignmentExportFormat): string {
  const records = normalize(input);
  const length = records[0].sequence.length;
  if (format === "fasta") return records.map((record) => `>${record.name}\n${record.sequence}\n`).join("");
  if (format === "phylip") {
    return `${records.length} ${length}\n${records.map((record, index) => `${safeName(record.name, index)} ${record.sequence}`).join("\n")}\n`;
  }
  if (format === "stockholm") {
    return `# STOCKHOLM 1.0\n${records.map((record, index) => `${safeName(record.name, index).padEnd(Math.max(...records.map((item, itemIndex) => safeName(item.name, itemIndex).length)) + 2)}${record.sequence}`).join("\n")}\n//\n`;
  }
  if (format === "nexus") {
    const datatype = records.every((record) => /^[ACGTUNRYKMSWBDHV?.-]*$/.test(record.sequence)) ? "DNA" : "PROTEIN";
    return `#NEXUS\n\nbegin data;\n  dimensions ntax=${records.length} nchar=${length};\n  format datatype=${datatype} gap=- missing=?;\n  matrix\n${records.map((record) => `    ${nexusName(record.name)} ${record.sequence}`).join("\n")}\n  ;\nend;\n`;
  }
  const width = 60;
  const nameWidth = Math.min(40, Math.max(12, ...records.map((record, index) => safeName(record.name, index).length)));
  let output = "CLUSTAL W multiple sequence alignment\n\n";
  for (let start = 0; start < length; start += width) {
    for (let index = 0; index < records.length; index += 1) {
      output += `${safeName(records[index].name, index).slice(0, nameWidth).padEnd(nameWidth + 2)}${records[index].sequence.slice(start, start + width)}\n`;
    }
    output += "\n";
  }
  return output;
}

export function treeNexus(newick: string, treeName = "swig_tree"): string {
  const clean = newick.trim().replace(/;+$/, "");
  if (!clean.startsWith("(")) throw new Error("A complete Newick tree is required.");
  return `#NEXUS\n\nbegin trees;\n  tree ${safeName(treeName, 0)} = [&R] ${clean};\nend;\n`;
}
