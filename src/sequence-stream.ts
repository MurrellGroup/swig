export type SequenceFormat = 1 | 2 | 3;

export interface SequenceBatch {
  index: number;
  text: string;
  count: number;
  format: SequenceFormat;
}

export interface SequenceStreamProgress {
  bytesRead: number;
  totalBytes: number;
  recordsRead: number;
  recordsSelected: number;
  maxBatchCharacters: number;
  maxCarryCharacters: number;
}

export interface SequenceSubsample {
  size: number;
  seed: number;
}

export interface SequenceStreamOptions {
  source: string | File;
  format: SequenceFormat;
  batchSize?: number;
  subsample?: SequenceSubsample;
  signal?: AbortSignal;
  onProgress?: (progress: SequenceStreamProgress) => void;
}

const STRING_CHUNK_SIZE = 256 * 1024;

function abortIfNeeded(signal?: AbortSignal) {
  if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
}

async function* decodedChunks(
  source: string | File,
  signal: AbortSignal | undefined,
  onBytes: (bytesRead: number, totalBytes: number) => void,
): AsyncGenerator<string> {
  if (typeof source === "string") {
    const total = source.length;
    for (let offset = 0; offset < total; offset += STRING_CHUNK_SIZE) {
      abortIfNeeded(signal);
      const next = Math.min(offset + STRING_CHUNK_SIZE, total);
      yield source.slice(offset, next);
      onBytes(next, total);
      await Promise.resolve();
    }
    return;
  }

  const total = source.size;
  let bytesRead = 0;
  const measured = source.stream().pipeThrough(new TransformStream<Uint8Array, Uint8Array>({
    transform(chunk, controller) {
      abortIfNeeded(signal);
      bytesRead += chunk.byteLength;
      onBytes(bytesRead, total);
      controller.enqueue(chunk);
    },
  }));
  const bytes = source.name.toLowerCase().endsWith(".gz")
    ? measured.pipeThrough(makeGzipDecoder())
    : measured;
  const reader = bytes.getReader();
  const decoder = new TextDecoder();
  try {
    while (true) {
      abortIfNeeded(signal);
      const { done, value } = await reader.read();
      if (done) break;
      const text = decoder.decode(value, { stream: true });
      if (text) yield text;
    }
    const tail = decoder.decode();
    if (tail) yield tail;
  } finally {
    reader.releaseLock();
  }
}

function makeGzipDecoder(): TransformStream<Uint8Array, Uint8Array> {
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot stream gzip input. Decompress the file first.");
  }
  return new DecompressionStream("gzip") as TransformStream<Uint8Array, Uint8Array>;
}

async function* lines(
  source: string | File,
  signal: AbortSignal | undefined,
  onBytes: (bytesRead: number, totalBytes: number) => void,
  onCarry: (characters: number) => void,
): AsyncGenerator<string> {
  let carry = "";
  for await (const chunk of decodedChunks(source, signal, onBytes)) {
    carry += chunk;
    let start = 0;
    while (true) {
      const end = carry.indexOf("\n", start);
      if (end < 0) break;
      yield carry.slice(start, end).replace(/\r$/, "");
      start = end + 1;
    }
    carry = carry.slice(start);
    onCarry(carry.length);
  }
  if (carry) yield carry.replace(/\r$/, "");
}

async function* fastaRecords(source: AsyncIterable<string>): AsyncGenerator<string> {
  let header = "";
  let sequence: string[] = [];
  for await (const line of source) {
    if (line.startsWith(">")) {
      if (header) {
        const joined = sequence.join("");
        if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
        yield `${header}\n${joined}\n`;
      }
      header = line;
      sequence = [];
      continue;
    }
    if (!line.trim()) continue;
    if (!header) throw new Error("Expected a FASTA header beginning with '>'.");
    sequence.push(line.replace(/\s/g, ""));
  }
  if (header) {
    const joined = sequence.join("");
    if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
    yield `${header}\n${joined}\n`;
  }
}

async function* fastqRecords(source: AsyncIterable<string>): AsyncGenerator<string> {
  let header = "";
  let plus = "";
  let sequence: string[] = [];
  let quality: string[] = [];
  let sequenceLength = 0;
  let qualityLength = 0;
  let state: "header" | "sequence" | "quality" = "header";

  for await (const line of source) {
    if (state === "header") {
      if (!line.trim()) continue;
      if (!line.startsWith("@")) throw new Error("Expected a FASTQ header beginning with '@'.");
      header = line;
      plus = "";
      sequence = [];
      quality = [];
      sequenceLength = 0;
      qualityLength = 0;
      state = "sequence";
      continue;
    }
    if (state === "sequence") {
      if (line.startsWith("+")) {
        if (!sequenceLength) throw new Error(`Empty FASTQ sequence: ${header.slice(1).trim() || "unnamed sequence"}.`);
        plus = line;
        state = "quality";
      } else {
        const normalized = line.replace(/\s/g, "");
        sequence.push(normalized);
        sequenceLength += normalized.length;
      }
      continue;
    }

    quality.push(line);
    qualityLength += line.length;
    if (qualityLength > sequenceLength) {
      throw new Error(`FASTQ sequence/quality length mismatch: ${header.slice(1).trim() || "unnamed sequence"}.`);
    }
    if (qualityLength === sequenceLength) {
      yield `${header}\n${sequence.join("")}\n${plus}\n${quality.join("")}\n`;
      state = "header";
    }
  }

  if (state === "sequence") throw new Error(`The FASTQ input ended before a '+' line for ${header.slice(1).trim() || "a sequence"}.`);
  if (state === "quality") throw new Error(`The FASTQ quality string is truncated for ${header.slice(1).trim() || "a sequence"}.`);
}

async function* airrRecords(source: AsyncIterable<string>): AsyncGenerator<{ header: string; row: string }> {
  let header = "";
  let delimiter = "\t";
  let sequenceColumn = -1;
  for await (const line of source) {
    if (!header) {
      if (!line.trim()) continue;
      header = line;
      delimiter = header.includes("\t") ? "\t" : ",";
      sequenceColumn = header.split(delimiter).indexOf("sequence");
      if (sequenceColumn < 0) throw new Error("AIRR input requires a 'sequence' column.");
      continue;
    }
    if (!line.trim()) continue;
    const values = line.split(delimiter);
    if (values[sequenceColumn]) yield { header, row: line };
  }
}

interface StreamRecord {
  ordinal: number;
  text: string;
  header?: string;
}

function seededRandom(seed: number): () => number {
  let value = Math.trunc(seed) >>> 0;
  return () => {
    value = (value + 0x6d2b79f5) >>> 0;
    let mixed = value;
    mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
    mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
    return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
  };
}

export async function* streamSequenceBatches(options: SequenceStreamOptions): AsyncGenerator<SequenceBatch> {
  const batchSize = Math.max(1, Math.floor(options.batchSize ?? 1000));
  const subsampleSize = options.subsample
    ? Math.max(1, Math.floor(options.subsample.size))
    : 0;
  let bytesRead = 0;
  let totalBytes = typeof options.source === "string" ? options.source.length : options.source.size;
  let recordsRead = 0;
  let recordsSelected = 0;
  let maxBatchCharacters = 0;
  let maxCarryCharacters = 0;
  const report = () => options.onProgress?.({
    bytesRead,
    totalBytes,
    recordsRead,
    recordsSelected,
    maxBatchCharacters,
    maxCarryCharacters,
  });
  const inputLines = lines(
    options.source,
    options.signal,
    (read, total) => {
      bytesRead = read;
      totalBytes = total;
      report();
    },
    (characters) => {
      maxCarryCharacters = Math.max(maxCarryCharacters, characters);
    },
  );

  const parsedRecords = async function* (): AsyncGenerator<StreamRecord> {
    if (options.format === 3) {
      for await (const record of airrRecords(inputLines)) {
        abortIfNeeded(options.signal);
        recordsRead += 1;
        if (recordsRead % 1000 === 0) report();
        yield { ordinal: recordsRead - 1, text: `${record.row}\n`, header: record.header };
      }
      return;
    }
    const recordSource = options.format === 1 ? fastaRecords(inputLines) : fastqRecords(inputLines);
    for await (const record of recordSource) {
      abortIfNeeded(options.signal);
      recordsRead += 1;
      if (recordsRead % 1000 === 0) report();
      yield { ordinal: recordsRead - 1, text: record };
    }
  };

  let selected: AsyncIterable<StreamRecord> = parsedRecords();
  if (subsampleSize) {
    const random = seededRandom(options.subsample?.seed ?? 1);
    const reservoir: StreamRecord[] = [];
    for await (const record of selected) {
      if (reservoir.length < subsampleSize) {
        reservoir.push(record);
      } else {
        const replacement = Math.floor(random() * recordsRead);
        if (replacement < subsampleSize) reservoir[replacement] = record;
      }
      recordsSelected = Math.min(recordsRead, subsampleSize);
    }
    reservoir.sort((a, b) => a.ordinal - b.ordinal);
    selected = (async function* () { yield* reservoir; })();
    report();
  }

  let batchIndex = 0;
  let records: StreamRecord[] = [];
  const emit = (): SequenceBatch => {
    const header = options.format === 3 ? records[0]?.header : undefined;
    const body = records.map((record) => record.text).join("");
    const text = header ? `${header}\n${body}` : body;
    maxBatchCharacters = Math.max(maxBatchCharacters, text.length);
    const batch = { index: batchIndex++, text, count: records.length, format: options.format };
    records = [];
    report();
    return batch;
  };

  for await (const record of selected) {
    abortIfNeeded(options.signal);
    records.push(record);
    recordsSelected = subsampleSize ? Math.min(recordsRead, subsampleSize) : recordsRead;
    if (records.length === batchSize) yield emit();
  }
  if (records.length) yield emit();

  bytesRead = totalBytes;
  recordsSelected = subsampleSize ? Math.min(recordsRead, subsampleSize) : recordsRead;
  report();
  if (!recordsRead) throw new Error("No sequence records were found in the input.");
}
