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

interface ChunkRecord {
  index: number;
  text: string;
}

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

  constructor() {
    this.database = new Promise((resolve, reject) => {
      const request = indexedDB.open(this.databaseName, 1);
      request.onupgradeneeded = () => {
        const database = request.result;
        database.createObjectStore("chunks", { keyPath: "index" });
        const records = database.createObjectStore("records", { keyPath: "ordinal" });
        records.createIndex("sequenceId", "sequenceId");
        records.createIndex("locus", "locus");
        records.createIndex("productive", "productive");
        records.createIndex("vCall", "vCall");
        records.createIndex("dCall", "dCall");
        records.createIndex("jCall", "jCall");
        records.createIndex("vjInFrame", "vjInFrame");
        records.createIndex("stopCodon", "stopCodon");
        records.createIndex("completeVdj", "completeVdj");
        records.createIndex("revComp", "revComp");
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

  async appendBatch(headerLine: string, body: string): Promise<void> {
    if (!this.headerLine) {
      this.headerLine = headerLine.replace(/\r$/, "");
      this.headers = this.headerLine.split("\t");
    } else if (this.headerLine !== headerLine.replace(/\r$/, "")) {
      throw new Error("SwiftIG returned inconsistent AIRR columns between batches.");
    }

    const lines = body.split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
    if (!lines.length) return;
    const database = await this.database;
    const transaction = database.transaction(["chunks", "records"], "readwrite");
    const chunks = transaction.objectStore("chunks");
    const records = transaction.objectStore("records");
    const chunkIndex = this.nextChunk++;
    chunks.put({ index: chunkIndex, text: body } satisfies ChunkRecord);

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
      records.put(record);
      if (record.vCall && record.jCall) this.assigned += 1;
      if (record.productive === "T") this.productive += 1;
      if (record.cdr3 || record.cdr3Aa) this.withCdr3 += 1;
      bump(this.facetMaps.loci, record.locus);
      bump(this.facetMaps.productive, record.productive);
      bump(this.facetMaps.vCalls, record.vCall);
      bump(this.facetMaps.dCalls, record.dCall);
      bump(this.facetMaps.jCalls, record.jCall);
    }
    await transactionDone(transaction);
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
      const rows = await requestResult(store.getAll(IDBKeyRange.bound(offset, upper), limit)) as AirrIndexRecord[];
      return { rows, hasMore: offset + rows.length < this.count, totalMatches: this.count, scanned: rows.length };
    }

    let source: IDBObjectStore | IDBIndex = store;
    let range: IDBKeyRange | undefined;
    const exactCandidates: Array<[keyof ResultFilters, string]> = [
      ["vCall", "vCall"], ["jCall", "jCall"], ["dCall", "dCall"],
      ["locus", "locus"], ["productive", "productive"],
      ["vjInFrame", "vjInFrame"], ["stopCodon", "stopCodon"],
      ["completeVdj", "completeVdj"], ["revComp", "revComp"],
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
        const record = cursor.value as AirrIndexRecord;
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
    const line = chunk.text.split("\n")[record.line]?.replace(/\r$/, "") ?? "";
    const values = line.split("\t");
    return Object.fromEntries(this.headers.map((header, index) => [header, values[index] ?? ""]));
  }

  async airrBlob(): Promise<Blob> {
    const database = await this.database;
    const transaction = database.transaction("chunks", "readonly");
    const chunks = await requestResult(transaction.objectStore("chunks").getAll()) as ChunkRecord[];
    const parts: BlobPart[] = [`${this.headerLine}\n`, ...chunks.map((chunk) => chunk.text)];
    return new Blob(parts, { type: "text/tab-separated-values;charset=utf-8" });
  }

  async writeAirr(write: (part: string) => Promise<void>): Promise<void> {
    await write(`${this.headerLine}\n`);
    const database = await this.database;
    for (let index = 0; index < this.nextChunk; index += 1) {
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(index)) as ChunkRecord | undefined;
      if (chunk) await write(chunk.text);
    }
  }

  async clear(): Promise<void> {
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
