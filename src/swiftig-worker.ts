/// <reference lib="webworker" />

import { WASI } from "@bjorn3/browser_wasi_shim";

interface StartRequest {
  type: "start";
  id: number;
  query: string | File;
  format: 0 | 1 | 2 | 3;
  references: { V: string; D: string; J: string; C: string };
  minimumIdentity: number;
  strand: 0 | 1 | 2;
}

interface AckRequest {
  type: "ack";
  id: number;
  batch: number;
}

interface SwiftIgExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  _initialize: () => void;
  swig_alloc: (size: number) => number;
  swig_free: (pointer: number) => void;
  swig_init_database: (
    vPointer: number, vSize: number, dPointer: number, dSize: number,
    jPointer: number, jSize: number, cPointer: number, cSize: number,
  ) => number;
  swig_annotate: (
    queryPointer: number, querySize: number, format: number,
    identityPerMille: number, strand: number,
  ) => number;
  swig_result_ptr: () => number;
  swig_result_len: () => number;
  swig_error_ptr: () => number;
  swig_error_len: () => number;
}

interface QueryBatch {
  text: string;
  count: number;
}

let runtimePromise: Promise<SwiftIgExports> | null = null;
const acknowledgements = new Map<string, () => void>();
const encoder = new TextEncoder();
const decoder = new TextDecoder();

function progress(id: number, stage: string, value: number) {
  self.postMessage({ id, type: "progress", stage, value });
}

async function loadRuntime(): Promise<SwiftIgExports> {
  if (!runtimePromise) {
    runtimePromise = (async () => {
      const response = await fetch(`${import.meta.env.BASE_URL}swiftig.wasm`);
      if (!response.ok) throw new Error("The SwiftIG WebAssembly core could not be loaded.");
      const wasmModule = await WebAssembly.compile(await response.arrayBuffer());
      const wasi = new WASI([], [], []);
      const instance = await WebAssembly.instantiate(wasmModule, {
        wasi_snapshot_preview1: wasi.wasiImport,
      });
      wasi.initialize(instance as WebAssembly.Instance & {
        exports: { memory: WebAssembly.Memory; _initialize?: () => unknown };
      });
      return instance.exports as SwiftIgExports;
    })();
  }
  return runtimePromise;
}

function put(exports: SwiftIgExports, value: string): [number, number] {
  const bytes = encoder.encode(value);
  const pointer = exports.swig_alloc(bytes.length);
  if (!pointer && bytes.length) throw new Error("SwiftIG ran out of browser memory.");
  new Uint8Array(exports.memory.buffer, pointer, bytes.length).set(bytes);
  return [pointer, bytes.length];
}

function read(exports: SwiftIgExports, pointer: number, length: number): string {
  return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
}

function readError(exports: SwiftIgExports): string {
  return read(exports, exports.swig_error_ptr(), exports.swig_error_len()) ||
    "SwiftIG could not complete the annotation.";
}

function lineAt(input: string, start: number): { value: string; next: number } {
  const end = input.indexOf("\n", start);
  const stop = end < 0 ? input.length : end;
  return {
    value: input.slice(start, stop).replace(/\r$/, ""),
    next: end < 0 ? input.length : end + 1,
  };
}

function* fastaRecords(input: string): Generator<string> {
  const starts = /^>/gm;
  let previous: number | null = null;
  for (let match = starts.exec(input); match; match = starts.exec(input)) {
    if (previous !== null) yield input.slice(previous, match.index);
    previous = match.index;
  }
  if (previous !== null) yield input.slice(previous);
}

function* fastqRecords(input: string): Generator<string> {
  let position = 0;
  while (position < input.length) {
    while (position < input.length) {
      const peek = lineAt(input, position);
      if (peek.value.trim()) break;
      position = peek.next;
    }
    if (position >= input.length) return;
    const start = position;
    const header = lineAt(input, position);
    if (!header.value.startsWith("@")) throw new Error("Expected a FASTQ header beginning with '@'.");
    position = header.next;
    let sequenceLength = 0;
    let foundPlus = false;
    while (position < input.length) {
      const line = lineAt(input, position);
      position = line.next;
      if (line.value.startsWith("+")) {
        foundPlus = true;
        break;
      }
      sequenceLength += line.value.replace(/\s/g, "").length;
    }
    if (!foundPlus) throw new Error("The FASTQ input ended before a '+' line.");
    let qualityLength = 0;
    while (position < input.length && qualityLength < sequenceLength) {
      const line = lineAt(input, position);
      position = line.next;
      qualityLength += line.value.length;
    }
    if (qualityLength !== sequenceLength) throw new Error("A FASTQ sequence and quality string differ in length.");
    yield input.slice(start, position);
  }
}

function* airrRecords(input: string): Generator<{ header: string; row: string }> {
  let position = 0;
  let header = "";
  while (position < input.length && !header) {
    const line = lineAt(input, position);
    position = line.next;
    header = line.value;
  }
  if (!header) return;
  const delimiter = header.includes("\t") ? "\t" : ",";
  const columns = header.split(delimiter);
  const sequenceColumn = columns.indexOf("sequence");
  if (sequenceColumn < 0) throw new Error("AIRR input requires a 'sequence' column.");
  while (position < input.length) {
    const line = lineAt(input, position);
    position = line.next;
    if (!line.value.trim()) continue;
    const values = line.value.split(delimiter);
    if (!values[sequenceColumn]) continue;
    yield { header, row: line.value };
  }
}

function resolvedFormat(input: string, format: number): 1 | 2 | 3 {
  if (format === 1 || format === 2 || format === 3) return format;
  const first = input.trimStart()[0];
  return first === ">" ? 1 : first === "@" ? 2 : 3;
}

function* queryBatches(input: string, format: number, batchSize: number): Generator<QueryBatch> {
  const actual = resolvedFormat(input, format);
  if (actual === 3) {
    let header = "";
    let rows: string[] = [];
    for (const record of airrRecords(input)) {
      header = record.header;
      rows.push(record.row);
      if (rows.length === batchSize) {
        yield { text: `${header}\n${rows.join("\n")}\n`, count: rows.length };
        rows = [];
      }
    }
    if (rows.length) yield { text: `${header}\n${rows.join("\n")}\n`, count: rows.length };
    return;
  }
  const source = actual === 1 ? fastaRecords(input) : fastqRecords(input);
  let records: string[] = [];
  for (const record of source) {
    records.push(record.endsWith("\n") ? record : `${record}\n`);
    if (records.length === batchSize) {
      yield { text: records.join(""), count: records.length };
      records = [];
    }
  }
  if (records.length) yield { text: records.join(""), count: records.length };
}

function countQueries(input: string, format: number): number {
  let count = 0;
  for (const batch of queryBatches(input, format, 2000)) count += batch.count;
  return count;
}

async function readQuery(query: string | File, id: number): Promise<string> {
  if (typeof query === "string") return query;
  progress(id, query.name.toLowerCase().endsWith(".gz") ? "Decompressing input" : "Reading input", 0.03);
  if (query.name.toLowerCase().endsWith(".gz")) {
    if (!("DecompressionStream" in globalThis)) {
      throw new Error("This browser cannot decompress gzip input.");
    }
    const stream = query.stream().pipeThrough(new DecompressionStream("gzip"));
    return new Response(stream).text();
  }
  const reader = query.stream().getReader();
  const textDecoder = new TextDecoder();
  let consumed = 0;
  let output = "";
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    consumed += value.byteLength;
    output += textDecoder.decode(value, { stream: true });
    progress(id, "Reading input", 0.03 + 0.05 * (consumed / Math.max(query.size, 1)));
  }
  return output + textDecoder.decode();
}

function waitForAcknowledgement(id: number, batch: number): Promise<void> {
  return new Promise((resolve) => acknowledgements.set(`${id}:${batch}`, resolve));
}

async function handleRequest(request: StartRequest) {
  try {
    const queryText = await readQuery(request.query, request.id);
    progress(request.id, "Counting sequences", 0.09);
    const total = countQueries(queryText, request.format);
    if (!total) throw new Error("No sequence records were found in the input.");
    progress(request.id, `Found ${total.toLocaleString()} sequences`, 0.12);

    const exports = await loadRuntime();
    progress(request.id, "Indexing V, D and J references", 0.16);
    const allocations = [
      put(exports, request.references.V), put(exports, request.references.D),
      put(exports, request.references.J), put(exports, request.references.C),
    ];
    const initialized = exports.swig_init_database(
      allocations[0][0], allocations[0][1], allocations[1][0], allocations[1][1],
      allocations[2][0], allocations[2][1], allocations[3][0], allocations[3][1],
    );
    allocations.forEach(([pointer]) => exports.swig_free(pointer));
    if (initialized < 0) throw new Error(readError(exports));
    progress(request.id, `Indexed ${initialized.toLocaleString()} germline alleles`, 0.24);

    const batchSize = total <= 10 ? 1 : total <= 1000 ? 50 : 500;
    let processed = 0;
    let batchNumber = 0;
    for (const batch of queryBatches(queryText, request.format, batchSize)) {
      const [queryPointer, querySize] = put(exports, batch.text);
      const count = exports.swig_annotate(
        queryPointer, querySize, request.format,
        Math.round(request.minimumIdentity * 1000), request.strand,
      );
      exports.swig_free(queryPointer);
      if (count < 0) throw new Error(readError(exports));
      const tsv = read(exports, exports.swig_result_ptr(), exports.swig_result_len());
      const breakAt = tsv.indexOf("\n");
      if (breakAt < 0) throw new Error("SwiftIG returned an invalid AIRR table.");
      processed += count;
      const acknowledgement = waitForAcknowledgement(request.id, batchNumber);
      self.postMessage({
        id: request.id,
        type: "batch",
        batch: batchNumber,
        header: tsv.slice(0, breakAt).replace(/\r$/, ""),
        body: tsv.slice(breakAt + 1),
        count,
        processed,
        total,
      });
      await acknowledgement;
      batchNumber += 1;
      progress(
        request.id,
        `Annotated ${processed.toLocaleString()} of ${total.toLocaleString()}`,
        0.24 + 0.72 * (processed / total),
      );
    }
    progress(request.id, "Finalizing local result index", 0.98);
    self.postMessage({ id: request.id, type: "result", count: processed, total });
  } catch (error) {
    self.postMessage({
      id: request.id,
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  }
}

self.addEventListener("message", (event: MessageEvent<StartRequest | AckRequest>) => {
  const message = event.data;
  if (message.type === "ack") {
    const key = `${message.id}:${message.batch}`;
    acknowledgements.get(key)?.();
    acknowledgements.delete(key);
    return;
  }
  void handleRequest(message);
});

export {};

