import type {
  GzipMemberRange,
  GzipMemberSource,
  SequenceFormat,
} from "./sequence-stream";

const SCAN_CHUNK_BYTES = 4 * 1024 * 1024;
const MAX_HEADER_BYTES = 1024 * 1024;
const PREVIEW_BYTES = 128 * 1024;

interface ParsedGzipHeader {
  offset: number;
  name: string;
}

export interface GzipMemberInspection extends GzipMemberRange {
  compressedBytes: number;
  uncompressedBytes: number | null;
  headerName: string;
  firstLine: string;
  firstRecordId: string;
  startsAtRecordBoundary: boolean;
}

function gzipDecoder(): TransformStream<Uint8Array, Uint8Array> {
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot inspect gzip members. Decompress the file first.");
  }
  return new DecompressionStream("gzip") as TransformStream<Uint8Array, Uint8Array>;
}

async function fixedHeaderCandidates(file: File): Promise<number[]> {
  const candidates = new Set<number>();
  for (let offset = 0; offset < file.size; offset += SCAN_CHUNK_BYTES) {
    const start = Math.max(0, offset - 3);
    const end = Math.min(file.size, offset + SCAN_CHUNK_BYTES);
    const bytes = new Uint8Array(await file.slice(start, end).arrayBuffer());
    for (let index = 0; index + 3 < bytes.length; index += 1) {
      if (
        bytes[index] === 0x1f
        && bytes[index + 1] === 0x8b
        && bytes[index + 2] === 0x08
        && (bytes[index + 3] & 0xe0) === 0
      ) candidates.add(start + index);
    }
  }
  return [...candidates].sort((a, b) => a - b);
}

function zeroTerminatedEnd(bytes: Uint8Array, start: number): number {
  for (let index = start; index < bytes.length; index += 1) {
    if (bytes[index] === 0) return index;
  }
  return -1;
}

async function parseGzipHeader(file: File, offset: number): Promise<ParsedGzipHeader | null> {
  const available = Math.min(MAX_HEADER_BYTES, file.size - offset);
  if (available < 10) return null;
  const bytes = new Uint8Array(await file.slice(offset, offset + available).arrayBuffer());
  if (
    bytes[0] !== 0x1f
    || bytes[1] !== 0x8b
    || bytes[2] !== 0x08
    || (bytes[3] & 0xe0) !== 0
  ) return null;

  const flags = bytes[3];
  let position = 10;
  if (flags & 0x04) {
    if (position + 2 > bytes.length) return null;
    const extraLength = bytes[position] | (bytes[position + 1] << 8);
    position += 2 + extraLength;
    if (position > bytes.length) return null;
  }

  let name = "";
  if (flags & 0x08) {
    const end = zeroTerminatedEnd(bytes, position);
    if (end < 0) return null;
    name = new TextDecoder().decode(bytes.slice(position, end));
    position = end + 1;
  }
  if (flags & 0x10) {
    const end = zeroTerminatedEnd(bytes, position);
    if (end < 0) return null;
    position = end + 1;
  }
  if (flags & 0x02) position += 2;
  if (position > bytes.length) return null;
  return { offset, name };
}

function firstNonemptyLine(text: string): string {
  return text.split(/\r?\n/g).find((line) => line.trim().length > 0) ?? "";
}

function startsAtRecordBoundary(line: string, format: SequenceFormat): boolean {
  if (format === 1) return line.startsWith(">");
  if (format === 2) return line.startsWith("@");
  const delimiter = line.includes("\t") ? "\t" : ",";
  return line.split(delimiter).includes("sequence");
}

function recordId(line: string, format: SequenceFormat): string {
  if ((format === 1 && line.startsWith(">")) || (format === 2 && line.startsWith("@"))) {
    return line.slice(1).trim().split(/\s+/, 1)[0] ?? "";
  }
  return "";
}

async function inspectMember(
  file: File,
  header: ParsedGzipHeader,
  end: number,
  format: SequenceFormat,
): Promise<GzipMemberInspection> {
  const reader = file.slice(header.offset, end).stream().pipeThrough(gzipDecoder()).getReader();
  const preview: Uint8Array[] = [];
  let previewLength = 0;
  let uncompressedBytes = 0;
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      uncompressedBytes += value.byteLength;
      if (previewLength < PREVIEW_BYTES) {
        const take = Math.min(value.byteLength, PREVIEW_BYTES - previewLength);
        preview.push(value.slice(0, take));
        previewLength += take;
      }
    }
  } finally {
    reader.releaseLock();
  }

  const prefix = new Uint8Array(previewLength);
  let position = 0;
  for (const chunk of preview) {
    prefix.set(chunk, position);
    position += chunk.byteLength;
  }
  const firstLine = firstNonemptyLine(new TextDecoder().decode(prefix));
  return {
    start: header.offset,
    end,
    compressedBytes: end - header.offset,
    uncompressedBytes,
    headerName: header.name,
    firstLine,
    firstRecordId: recordId(firstLine, format),
    startsAtRecordBoundary: startsAtRecordBoundary(firstLine, format),
  };
}

/**
 * Locate and validate every gzip member. Single-member inputs take the cheap
 * path (header scan only); multi-member candidates are each decompressed once
 * so an accidental gzip signature inside DEFLATE data cannot split a sample.
 */
export async function inspectGzipMembers(
  file: File,
  format: SequenceFormat,
): Promise<GzipMemberInspection[]> {
  const fixedCandidates = await fixedHeaderCandidates(file);
  const parsed = (await Promise.all(fixedCandidates.map((offset) => parseGzipHeader(file, offset))))
    .filter((header): header is ParsedGzipHeader => header !== null);
  if (!parsed.length || parsed[0].offset !== 0) {
    throw new Error("The selected .gz file does not begin with a valid gzip header.");
  }
  if (parsed.length === 1) {
    return [{
      start: 0,
      end: file.size,
      compressedBytes: file.size,
      uncompressedBytes: null,
      headerName: parsed[0].name,
      firstLine: "",
      firstRecordId: "",
      startsAtRecordBoundary: true,
    }];
  }

  const headerByOffset = new Map(parsed.map((header) => [header.offset, header]));
  const candidateEnds = [...parsed.slice(1).map((header) => header.offset), file.size];
  const members: GzipMemberInspection[] = [];
  let start = 0;
  let lastError: unknown;
  while (start < file.size) {
    const header = headerByOffset.get(start);
    if (!header) throw new Error(`Could not validate the gzip member beginning at byte ${start.toLocaleString()}.`);
    let accepted: GzipMemberInspection | null = null;
    for (const end of candidateEnds) {
      if (end <= start) continue;
      try {
        accepted = await inspectMember(file, header, end, format);
        break;
      } catch (error) {
        lastError = error;
      }
    }
    if (!accepted) {
      const detail = lastError instanceof Error ? ` ${lastError.message}` : "";
      throw new Error(`The gzip member beginning at byte ${start.toLocaleString()} is invalid.${detail}`);
    }
    members.push(accepted);
    start = accepted.end;
  }
  return members;
}

export function gzipMemberSource(file: File, members: GzipMemberRange[]): GzipMemberSource {
  return {
    kind: "gzip-members",
    file,
    members: members.map(({ start, end }) => ({ start, end })),
  };
}

export function gzipMemberFile(file: File, member: GzipMemberRange, name: string): File {
  return new File([file.slice(member.start, member.end)], name, {
    type: file.type || "application/gzip",
    lastModified: file.lastModified,
  });
}

function basename(name: string): string {
  return name.split(/[\\/]/g).pop() ?? name;
}

function gzipSuffix(name: string): string {
  return name.match(/\.(?:fa|fasta|fna|fas|fq|fastq|tsv|csv|txt)\.gz$/i)?.[0]
    ?? name.match(/\.gz$/i)?.[0]
    ?? ".gz";
}

function withoutGzipSuffix(name: string): string {
  const suffix = gzipSuffix(name);
  return name.slice(0, Math.max(0, name.length - suffix.length));
}

export function suggestedGzipMemberName(
  originalName: string,
  member: GzipMemberInspection,
  index: number,
): string {
  const suffix = gzipSuffix(originalName);
  if (member.headerName) {
    const named = basename(member.headerName).replace(/[^A-Za-z0-9_.-]+/g, "_");
    if (/\.gz$/i.test(named)) return named;
    return `${named}${suffix}`;
  }
  const accession = member.firstRecordId.match(/^([A-Za-z]{3}\d+)\.\d+$/)?.[1];
  if (accession) return `${accession}${suffix}`;
  const stem = withoutGzipSuffix(basename(originalName)).replace(/[^A-Za-z0-9_.-]+/g, "_") || "dataset";
  return `${stem}.member-${index + 1}${suffix}`;
}
