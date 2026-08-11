export interface AirrIndexRecord {
  ordinal: number;
  chunk: number;
  line: number;
  sequenceId: string;
  locus: string;
  vCall: string;
  dCall: string;
  jCall: string;
  productive: string;
  cdr3: string;
  cdr3Aa: string;
  junctionAa: string;
  vIdentity: number | null;
  dIdentity: number | null;
  jIdentity: number | null;
  cdr3AaLength: number | null;
  vjInFrame: string;
  stopCodon: string;
  completeVdj: string;
  revComp: string;
}

export interface ResultFilters {
  sequenceId: string;
  cdr3: string;
  locus: string;
  productive: string;
  vCall: string;
  dCall: string;
  jCall: string;
  minVIdentity: number;
  minDIdentity: number;
  minJIdentity: number;
  minCdr3AaLength: number;
  maxCdr3AaLength: number;
  vjInFrame: string;
  stopCodon: string;
  completeVdj: string;
  revComp: string;
  hasD: boolean;
  hasCdr3: boolean;
}

export interface FacetValue {
  value: string;
  count: number;
}

export interface ResultFacets {
  loci: FacetValue[];
  productive: FacetValue[];
  vCalls: FacetValue[];
  dCalls: FacetValue[];
  jCalls: FacetValue[];
}

export interface ResultPage {
  rows: AirrIndexRecord[];
  hasMore: boolean;
  totalMatches: number | null;
  scanned: number;
}

type AirrRow = Record<string, string>;

interface PackedIndexRecord {
  o: number;
  c: number;
  n: number;
  i: string;
  l: string;
  v: string;
  d: string;
  j: string;
  p: string;
  r: string;
  a: string;
  u: string;
  vi: number | null;
  di: number | null;
  ji: number | null;
  z: number | null;
  f: string;
  s: string;
  q: string;
  x: string;
}

interface ChunkRecord {
  index: number;
  storage: "indexed" | "external";
  data?: Blob;
  compressed?: boolean;
  start?: number;
  length?: number;
}

interface ManifestRecord {
  key: "manifest";
  headerLine: string;
  chunks: number;
  records: number;
  outputBytes: number;
}

export interface AirrOutputWritable {
  write: (data: string | Blob | Uint8Array) => Promise<void>;
  close: () => Promise<void>;
  abort?: () => Promise<void>;
}

export interface AirrOutputHandle {
  getFile: () => Promise<File>;
}

export interface DirectAirrOutput {
  handle: AirrOutputHandle;
  writable: AirrOutputWritable;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();
const MAX_FALLBACK_BLOB_BYTES = 256 * 1024 * 1024;

function requestResult<T>(request: IDBRequest<T>): Promise<T> {
  return new Promise((resolve, reject) => {
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error ?? new Error("The local result index failed."));
  });
}

function transactionDone(transaction: IDBTransaction): Promise<void> {
  return new Promise((resolve, reject) => {
    transaction.oncomplete = () => resolve();
    transaction.onerror = () => reject(transaction.error ?? new Error("The local result index failed."));
    transaction.onabort = () => reject(transaction.error ?? new Error("The local result index was interrupted."));
  });
}

function bump(map: Map<string, number>, value: string) {
  if (!value) return;
  map.set(value, (map.get(value) ?? 0) + 1);
}

function facet(map: Map<string, number>): FacetValue[] {
  return [...map.entries()]
    .map(([value, count]) => ({ value, count }))
    .sort((a, b) => b.count - a.count || a.value.localeCompare(b.value));
}

function numeric(value: string): number | null {
  const parsed = Number(value);
  return value !== "" && Number.isFinite(parsed) ? parsed : null;
}

function packRecord(record: AirrIndexRecord): PackedIndexRecord {
  return {
    o: record.ordinal, c: record.chunk, n: record.line, i: record.sequenceId,
    l: record.locus, v: record.vCall, d: record.dCall, j: record.jCall,
    p: record.productive, r: record.cdr3, a: record.cdr3Aa, u: record.junctionAa,
    vi: record.vIdentity, di: record.dIdentity, ji: record.jIdentity,
    z: record.cdr3AaLength, f: record.vjInFrame, s: record.stopCodon,
    q: record.completeVdj, x: record.revComp,
  };
}

function unpackRecord(record: PackedIndexRecord): AirrIndexRecord {
  return {
    ordinal: record.o, chunk: record.c, line: record.n, sequenceId: record.i,
    locus: record.l, vCall: record.v, dCall: record.d, jCall: record.j,
    productive: record.p, cdr3: record.r, cdr3Aa: record.a, junctionAa: record.u,
    vIdentity: record.vi, dIdentity: record.di, jIdentity: record.ji,
    cdr3AaLength: record.z, vjInFrame: record.f, stopCodon: record.s,
    completeVdj: record.q, revComp: record.x,
  };
}

export const EMPTY_FILTERS: ResultFilters = {
  sequenceId: "",
  cdr3: "",
  locus: "",
  productive: "",
  vCall: "",
  dCall: "",
  jCall: "",
  minVIdentity: 0,
  minDIdentity: 0,
  minJIdentity: 0,
  minCdr3AaLength: 0,
  maxCdr3AaLength: 0,
  vjInFrame: "",
  stopCodon: "",
  completeVdj: "",
  revComp: "",
  hasD: false,
  hasCdr3: false,
};

export class AirrResultStore {
  readonly databaseName = `swig-results-${Date.now()}-${crypto.randomUUID()}`;
  private readonly database: Promise<IDBDatabase>;
  private readonly facetMaps = {
    loci: new Map<string, number>(),
    productive: new Map<string, number>(),
    vCalls: new Map<string, number>(),
    dCalls: new Map<string, number>(),
    jCalls: new Map<string, number>(),
  };
  private headerLine = "";
  private headers: string[] = [];
  private nextChunk = 0;
  private nextOrdinal = 0;
  private assigned = 0;
  private productive = 0;
  private withCdr3 = 0;
  private readonly directOutput?: DirectAirrOutput;
  private outputByteCount = 0;
  private finalized = false;

  constructor(directOutput?: DirectAirrOutput) {
    this.directOutput = directOutput;
    this.database = new Promise((resolve, reject) => {
      const request = indexedDB.open(this.databaseName, 1);
      request.onupgradeneeded = () => {
        const database = request.result;
        database.createObjectStore("chunks", { keyPath: "index" });
        database.createObjectStore("meta", { keyPath: "key" });
        const records = database.createObjectStore("records", { keyPath: "o" });
        records.createIndex("locus", "l");
        records.createIndex("productive", "p");
        records.createIndex("vCall", "v");
        records.createIndex("dCall", "d");
        records.createIndex("jCall", "j");
      };
      request.onsuccess = () => resolve(request.result);
      request.onerror = () => reject(request.error ?? new Error("Could not create the local result index."));
    });
  }

  get count(): number {
    return this.nextOrdinal;
  }

  get airrHeaders(): string[] {
    return this.headers;
  }

  get outputBytes(): number {
    return this.outputByteCount;
  }

  get streamedDirectly(): boolean {
    return Boolean(this.directOutput);
  }

  get summary() {
    return { assigned: this.assigned, productive: this.productive, withCdr3: this.withCdr3 };
  }

  facets(): ResultFacets {
    return {
      loci: facet(this.facetMaps.loci),
      productive: facet(this.facetMaps.productive),
      vCalls: facet(this.facetMaps.vCalls),
      dCalls: facet(this.facetMaps.dCalls),
      jCalls: facet(this.facetMaps.jCalls),
    };
  }

  async appendBatch(headerLine: string, body: string | Uint8Array): Promise<void> {
    if (this.finalized) throw new Error("Cannot append to a finalized AIRR result store.");
    const normalizedHeader = headerLine.replace(/\r$/, "");
    if (!this.headerLine) {
      this.headerLine = normalizedHeader;
      this.headers = this.headerLine.split("\t");
      const headerBytes = encoder.encode(`${this.headerLine}\n`);
      this.outputByteCount += headerBytes.byteLength;
      if (this.directOutput) await this.directOutput.writable.write(headerBytes);
    } else if (this.headerLine !== normalizedHeader) {
      throw new Error("SwiftIG returned inconsistent AIRR columns between batches.");
    }

    const bodyBytes = typeof body === "string" ? encoder.encode(body) : body;
    const bodyText = typeof body === "string" ? body : decoder.decode(body);
    const lines = bodyText.split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
    if (!lines.length) return;
    const chunkIndex = this.nextChunk;
    let chunk: ChunkRecord;
    if (this.directOutput) {
      const start = this.outputByteCount;
      await this.directOutput.writable.write(bodyBytes);
      chunk = { index: chunkIndex, storage: "external", start, length: bodyBytes.byteLength };
    } else {
      let data: Blob;
      let compressed = false;
      if ("CompressionStream" in globalThis) {
        const raw = new Blob([bodyBytes.slice().buffer]);
        data = await new Response(raw.stream().pipeThrough(new CompressionStream("gzip"))).blob();
        compressed = true;
      } else {
        data = new Blob([bodyBytes.slice().buffer]);
      }
      chunk = { index: chunkIndex, storage: "indexed", data, compressed };
    }
    this.outputByteCount += bodyBytes.byteLength;

    const database = await this.database;
    const transaction = database.transaction(["chunks", "records"], "readwrite");
    const chunks = transaction.objectStore("chunks");
    const records = transaction.objectStore("records");
    chunks.put(chunk);

    const positions = Object.fromEntries(this.headers.map((name, index) => [name, index]));
    const at = (values: string[], name: string) => values[positions[name]] ?? "";
    for (let line = 0; line < lines.length; line += 1) {
      const values = lines[line].split("\t");
      const record: AirrIndexRecord = {
        ordinal: this.nextOrdinal++,
        chunk: chunkIndex,
        line,
        sequenceId: at(values, "sequence_id"),
        locus: at(values, "locus"),
        vCall: at(values, "v_call"),
        dCall: at(values, "d_call"),
        jCall: at(values, "j_call"),
        productive: at(values, "productive"),
        cdr3: at(values, "cdr3"),
        cdr3Aa: at(values, "cdr3_aa"),
        junctionAa: at(values, "junction_aa"),
        vIdentity: numeric(at(values, "v_identity")),
        dIdentity: numeric(at(values, "d_identity")),
        jIdentity: numeric(at(values, "j_identity")),
        cdr3AaLength: at(values, "cdr3_aa") ? at(values, "cdr3_aa").length : null,
        vjInFrame: at(values, "vj_in_frame"),
        stopCodon: at(values, "stop_codon"),
        completeVdj: at(values, "complete_vdj"),
        revComp: at(values, "rev_comp"),
      };
      records.put(packRecord(record));
      if (record.vCall && record.jCall) this.assigned += 1;
      if (record.productive === "T") this.productive += 1;
      if (record.cdr3 || record.cdr3Aa) this.withCdr3 += 1;
      bump(this.facetMaps.loci, record.locus);
      bump(this.facetMaps.productive, record.productive);
      bump(this.facetMaps.vCalls, record.vCall);
      bump(this.facetMaps.dCalls, record.dCall);
      bump(this.facetMaps.jCalls, record.jCall);
    }
    try {
      await transactionDone(transaction);
      this.nextChunk += 1;
    } catch (error) {
      if (error instanceof DOMException && error.name === "QuotaExceededError") {
        throw new Error("Browser storage is full. Re-run with ‘stream AIRR to disk’ so output is never retained in browser storage.");
      }
      throw error;
    }
  }

  async page(
    filters: ResultFilters,
    offset: number,
    limit: number,
    onScan?: (scanned: number) => void,
    signal?: AbortSignal,
  ): Promise<ResultPage> {
    const database = await this.database;
    const transaction = database.transaction("records", "readonly");
    const store = transaction.objectStore("records");
    const normalizedId = filters.sequenceId.trim().toLowerCase();
    const normalizedCdr3 = filters.cdr3.trim().toUpperCase();
    const filtered = Boolean(
      normalizedId || normalizedCdr3 || filters.locus || filters.productive ||
      filters.vCall || filters.dCall || filters.jCall || filters.minVIdentity ||
      filters.minDIdentity || filters.minJIdentity || filters.minCdr3AaLength ||
      filters.maxCdr3AaLength || filters.vjInFrame || filters.stopCodon ||
      filters.completeVdj || filters.revComp || filters.hasD || filters.hasCdr3,
    );

    if (!filtered) {
      if (offset >= this.count) return { rows: [], hasMore: false, totalMatches: this.count, scanned: 0 };
      const upper = Math.min(this.count - 1, offset + limit - 1);
      const packed = await requestResult(store.getAll(IDBKeyRange.bound(offset, upper), limit)) as PackedIndexRecord[];
      const rows = packed.map(unpackRecord);
      return { rows, hasMore: offset + rows.length < this.count, totalMatches: this.count, scanned: rows.length };
    }

    let source: IDBObjectStore | IDBIndex = store;
    let range: IDBKeyRange | undefined;
    const exactCandidates: Array<[keyof ResultFilters, string]> = [
      ["vCall", "vCall"], ["jCall", "jCall"], ["dCall", "dCall"],
      ["locus", "locus"], ["productive", "productive"],
    ];
    for (const [filterName, indexName] of exactCandidates) {
      const value = filters[filterName];
      if (typeof value === "string" && value) {
        source = store.index(indexName);
        range = IDBKeyRange.only(value);
        break;
      }
    }

    const matches = (record: AirrIndexRecord) => {
      if (normalizedId && !record.sequenceId.toLowerCase().includes(normalizedId)) return false;
      if (normalizedCdr3 && !record.cdr3.toUpperCase().includes(normalizedCdr3) &&
        !record.cdr3Aa.toUpperCase().includes(normalizedCdr3)) return false;
      if (filters.locus && record.locus !== filters.locus) return false;
      if (filters.productive && record.productive !== filters.productive) return false;
      if (filters.vCall && record.vCall !== filters.vCall) return false;
      if (filters.dCall && record.dCall !== filters.dCall) return false;
      if (filters.jCall && record.jCall !== filters.jCall) return false;
      if (filters.minVIdentity && (record.vIdentity ?? 0) < filters.minVIdentity) return false;
      if (filters.minDIdentity && (record.dIdentity ?? 0) < filters.minDIdentity) return false;
      if (filters.minJIdentity && (record.jIdentity ?? 0) < filters.minJIdentity) return false;
      if (filters.minCdr3AaLength && (record.cdr3AaLength ?? 0) < filters.minCdr3AaLength) return false;
      if (filters.maxCdr3AaLength && (record.cdr3AaLength ?? Number.POSITIVE_INFINITY) > filters.maxCdr3AaLength) return false;
      if (filters.vjInFrame && record.vjInFrame !== filters.vjInFrame) return false;
      if (filters.stopCodon && record.stopCodon !== filters.stopCodon) return false;
      if (filters.completeVdj && record.completeVdj !== filters.completeVdj) return false;
      if (filters.revComp && record.revComp !== filters.revComp) return false;
      if (filters.hasD && !record.dCall) return false;
      if (filters.hasCdr3 && !record.cdr3 && !record.cdr3Aa) return false;
      return true;
    };

    return new Promise((resolve, reject) => {
      const rows: AirrIndexRecord[] = [];
      let scanned = 0;
      let matched = 0;
      let resolved = false;
      const request = source.openCursor(range);
      request.onerror = () => reject(request.error ?? new Error("Could not search the local result index."));
      request.onsuccess = () => {
        if (signal?.aborted) {
          if (!resolved) {
            resolved = true;
            resolve({ rows: [], hasMore: false, totalMatches: null, scanned });
          }
          return;
        }
        const cursor = request.result;
        if (!cursor) {
          if (!resolved) resolve({ rows, hasMore: false, totalMatches: matched, scanned });
          return;
        }
        scanned += 1;
        if (scanned % 2500 === 0) onScan?.(scanned);
        const record = unpackRecord(cursor.value as PackedIndexRecord);
        if (matches(record)) {
          if (matched >= offset && rows.length < limit) rows.push(record);
          matched += 1;
          if (rows.length === limit && matched > offset + limit) {
            resolved = true;
            resolve({ rows, hasMore: true, totalMatches: null, scanned });
            return;
          }
        }
        cursor.continue();
      };
    });
  }

  async detail(record: AirrIndexRecord): Promise<AirrRow> {
    const database = await this.database;
    const transaction = database.transaction("chunks", "readonly");
    const chunk = await requestResult(transaction.objectStore("chunks").get(record.chunk)) as ChunkRecord | undefined;
    if (!chunk) throw new Error("That AIRR result is no longer present in the local index.");
    const text = await this.chunkText(chunk);
    const line = text.split("\n")[record.line]?.replace(/\r$/, "") ?? "";
    const values = line.split("\t");
    return Object.fromEntries(this.headers.map((header, index) => [header, values[index] ?? ""]));
  }

  private async chunkBlob(chunk: ChunkRecord): Promise<Blob> {
    if (chunk.storage === "external") {
      if (!this.directOutput || chunk.start === undefined || chunk.length === undefined) {
        throw new Error("The streamed AIRR output file is no longer available.");
      }
      const file = await this.directOutput.handle.getFile();
      return file.slice(chunk.start, chunk.start + chunk.length);
    }
    if (!chunk.data) throw new Error("An indexed AIRR output batch is missing.");
    if (!chunk.compressed) return chunk.data;
    if (!("DecompressionStream" in globalThis)) {
      throw new Error("This browser cannot decompress the local AIRR result batch.");
    }
    return new Response(chunk.data.stream().pipeThrough(new DecompressionStream("gzip"))).blob();
  }

  private async chunkText(chunk: ChunkRecord): Promise<string> {
    return (await this.chunkBlob(chunk)).text();
  }

  async finalize(): Promise<void> {
    if (this.finalized) return;
    const database = await this.database;
    const transaction = database.transaction("meta", "readwrite");
    transaction.objectStore("meta").put({
      key: "manifest",
      headerLine: this.headerLine,
      chunks: this.nextChunk,
      records: this.nextOrdinal,
      outputBytes: this.outputByteCount,
    } satisfies ManifestRecord);
    await transactionDone(transaction);
    await this.directOutput?.writable.close();
    this.finalized = true;
  }

  async abort(): Promise<void> {
    if (!this.finalized) {
      try {
        await this.directOutput?.writable.abort?.();
      } catch {
        // The browser may already have discarded an interrupted output stream.
      }
    }
  }

  streamingDownloadUrl(baseUrl: string, name: string): string {
    const url = new URL(`${baseUrl}__swig_download__`, globalThis.location?.href);
    url.searchParams.set("database", this.databaseName);
    url.searchParams.set("name", name);
    return url.href;
  }

  async airrBlob(): Promise<Blob> {
    if (this.directOutput) return this.directOutput.handle.getFile();
    if (this.outputByteCount > MAX_FALLBACK_BLOB_BYTES) {
      throw new Error("This AIRR table is too large for a memory-backed download. Enable the streaming download worker or re-run with direct-to-disk output.");
    }
    const database = await this.database;
    const parts: BlobPart[] = [`${this.headerLine}\n`];
    for (let index = 0; index < this.nextChunk; index += 1) {
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(index)) as ChunkRecord | undefined;
      if (chunk) parts.push(await this.chunkBlob(chunk));
    }
    return new Blob(parts, { type: "text/tab-separated-values;charset=utf-8" });
  }

  async writeAirr(write: (part: string | Blob | Uint8Array) => Promise<void>): Promise<void> {
    if (this.directOutput) {
      await write(await this.directOutput.handle.getFile());
      return;
    }
    await write(`${this.headerLine}\n`);
    const database = await this.database;
    for (let index = 0; index < this.nextChunk; index += 1) {
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(index)) as ChunkRecord | undefined;
      if (chunk) await write(await this.chunkBlob(chunk));
    }
  }

  async clear(): Promise<void> {
    await this.abort();
    const database = await this.database;
    database.close();
    await new Promise<void>((resolve) => {
      const request = indexedDB.deleteDatabase(this.databaseName);
      request.onsuccess = () => resolve();
      request.onerror = () => resolve();
      request.onblocked = () => resolve();
    });
  }
}
