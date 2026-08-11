export type CallResolution = "gene" | "allele";
export type AmbiguityPolicy = "overlap" | "top" | "strict";
export type DedupKey = "sequence" | "trimmed" | "cdr3" | "rearrangement";

export interface PostAnalysisRecord {
  ordinal: number;
  sequenceId: string;
  locus: string;
  vCall: string;
  jCall: string;
  cdr3Nt: string;
  cdr3Aa: string;
  productive: boolean;
  sequenceFingerprint: string;
  trimmedFingerprint: string;
  trimmedSketch?: Uint32Array;
}

export interface DedupResult {
  key: DedupKey;
  inputRecords: number;
  uniqueRecords: number;
  collapsedRecords: number;
  representatives: Int32Array;
  counts: Uint32Array;
  largestGroups: Array<{ ordinal: number; count: number }>;
}

export interface LineageOptions {
  identity: number;
  callResolution: CallResolution;
  ambiguity: AmbiguityPolicy;
  productiveOnly: boolean;
  requireSameLocus: boolean;
  maxCandidateComparisons: number;
}

export interface LineageSummary {
  id: number;
  representativeOrdinal: number;
  uniqueMembers: number;
  abundance: number;
  locus: string;
  vCalls: string[];
  jCalls: string[];
  cdr3Length: number;
}

export interface LineageResult {
  assignments: Int32Array;
  summaries: LineageSummary[];
  lineageCount: number;
  sizeHistogram: Array<{ label: string; count: number }>;
  vUsage: Array<{ call: string; lineages: number; abundance: number }>;
  jUsage: Array<{ call: string; lineages: number; abundance: number }>;
  assignedRecords: number;
  unassignedRecords: number;
  candidateComparisons: number;
  truncatedCandidates: number;
}

export type QueryTarget = "cdr3_nt" | "cdr3_aa" | "trimmed";
export type QueryMetric = "exact" | "substring" | "hamming" | "edit" | "sketch";

export interface QueryConstraint {
  locus?: string;
  vCall?: string;
  jCall?: string;
}

export interface QueryOptions {
  target: QueryTarget;
  metric: QueryMetric;
  identity: number;
  maxResults: number;
  locus?: string;
  vCall?: string;
  jCall?: string;
  callResolution: CallResolution;
  ambiguity: AmbiguityPolicy;
  productiveOnly: boolean;
  queryConstraints?: QueryConstraint[];
}

export interface QueryHit {
  ordinal: number;
  queryIndex: number;
  score: number;
  distance: number;
  matched: string;
}

const HASH_SEEDS = [
  0x811c9dc5, 0x9e3779b9, 0x85ebca6b, 0xc2b2ae35,
  0x27d4eb2f, 0x165667b1, 0xd3a2646c, 0xfd7046c5,
] as const;

function mix32(value: number): number {
  value ^= value >>> 16;
  value = Math.imul(value, 0x7feb352d);
  value ^= value >>> 15;
  value = Math.imul(value, 0x846ca68b);
  value ^= value >>> 16;
  return value >>> 0;
}

export function hashSequence(value: string, seed: number = HASH_SEEDS[0]): number {
  let hash = seed >>> 0;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193);
  }
  return mix32(hash);
}

export function sequenceFingerprint(value: string): string {
  const normalized = normalizeNt(value);
  return `${normalized.length}:${[0, 2, 4, 6].map((index) => hashSequence(normalized, HASH_SEEDS[index]).toString(36)).join(":")}`;
}

export function minHashSketch(value: string, k = 7): Uint32Array {
  const normalized = normalizeNt(value).replaceAll("N", "");
  const sketch = new Uint32Array(HASH_SEEDS.length);
  sketch.fill(0xffffffff);
  const width = Math.max(1, Math.min(k, normalized.length));
  if (!normalized.length) return sketch;
  for (let index = 0; index <= normalized.length - width; index += 1) {
    const kmer = normalized.slice(index, index + width);
    for (let seed = 0; seed < HASH_SEEDS.length; seed += 1) {
      const hash = hashSequence(kmer, HASH_SEEDS[seed]);
      if (hash < sketch[seed]) sketch[seed] = hash;
    }
  }
  return sketch;
}

export function sketchSimilarity(left: Uint32Array, right: Uint32Array): number {
  const length = Math.min(left.length, right.length);
  if (!length) return 0;
  let equal = 0;
  let informative = 0;
  for (let index = 0; index < length; index += 1) {
    if (left[index] === 0xffffffff && right[index] === 0xffffffff) continue;
    informative += 1;
    if (left[index] === right[index]) equal += 1;
  }
  return informative ? equal / informative : 0;
}

export function normalizeNt(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}

export function normalizeAa(value: string): string {
  return value.toUpperCase().replace(/[^A-Z*]/g, "");
}

export function normalizeCall(call: string, resolution: CallResolution): string {
  const clean = call.trim().toUpperCase();
  return resolution === "gene" ? clean.replace(/\*.*$/, "") : clean;
}

export function callSet(value: string, resolution: CallResolution, ambiguity: AmbiguityPolicy): string[] {
  const values = value.split(",").map((call) => normalizeCall(call, resolution)).filter(Boolean);
  const selected = ambiguity === "top" ? values.slice(0, 1) : [...new Set(values)].sort();
  return selected;
}

function callsCompatible(
  left: string,
  right: string,
  resolution: CallResolution,
  ambiguity: AmbiguityPolicy,
): boolean {
  const a = callSet(left, resolution, ambiguity);
  const b = callSet(right, resolution, ambiguity);
  if (!a.length || !b.length) return false;
  if (ambiguity === "strict") return a.length === b.length && a.every((value, index) => value === b[index]);
  const values = new Set(a);
  return b.some((value) => values.has(value));
}

export function hammingDistanceWithin(left: string, right: string, maximum: number, ambiguousN = true): number {
  if (left.length !== right.length) return maximum + 1;
  let distance = 0;
  for (let index = 0; index < left.length; index += 1) {
    if (left[index] !== right[index] && (!ambiguousN || (left[index] !== "N" && right[index] !== "N"))) {
      distance += 1;
      if (distance > maximum) return distance;
    }
  }
  return distance;
}

export function bandedEditDistance(left: string, right: string, maximum: number): number {
  if (Math.abs(left.length - right.length) > maximum) return maximum + 1;
  let previous = new Int32Array(right.length + 1);
  let current = new Int32Array(right.length + 1);
  for (let column = 0; column <= right.length; column += 1) previous[column] = column;
  for (let row = 1; row <= left.length; row += 1) {
    current.fill(maximum + 1);
    current[0] = row;
    const start = Math.max(1, row - maximum);
    const end = Math.min(right.length, row + maximum);
    let rowMinimum = maximum + 1;
    for (let column = start; column <= end; column += 1) {
      const substitution = left[row - 1] === right[column - 1] ? 0 : 1;
      current[column] = Math.min(
        previous[column] + 1,
        current[column - 1] + 1,
        previous[column - 1] + substitution,
      );
      rowMinimum = Math.min(rowMinimum, current[column]);
    }
    if (rowMinimum > maximum) return maximum + 1;
    [previous, current] = [current, previous];
  }
  return previous[right.length];
}

function dedupKey(record: PostAnalysisRecord, key: DedupKey): string {
  if (key === "sequence") return record.sequenceFingerprint;
  if (key === "trimmed") return record.trimmedFingerprint;
  if (key === "cdr3") return `${record.locus}\u0000${record.cdr3Nt}`;
  return `${record.locus}\u0000${record.vCall}\u0000${record.jCall}\u0000${record.cdr3Nt}`;
}

export function deduplicate(records: PostAnalysisRecord[], key: DedupKey): DedupResult {
  const representatives = new Int32Array(records.length);
  const counts = new Uint32Array(records.length);
  const seen = new Map<string, number>();
  for (let index = 0; index < records.length; index += 1) {
    const value = dedupKey(records[index], key);
    const previous = seen.get(value);
    if (previous === undefined) {
      seen.set(value, index);
      representatives[index] = index;
      counts[index] = 1;
    } else {
      representatives[index] = previous;
      counts[previous] += 1;
    }
  }
  const largestGroups = [...seen.values()]
    .filter((ordinal) => counts[ordinal] > 1)
    .sort((a, b) => counts[b] - counts[a] || a - b)
    .slice(0, 100)
    .map((ordinal) => ({ ordinal, count: counts[ordinal] }));
  return {
    key,
    inputRecords: records.length,
    uniqueRecords: seen.size,
    collapsedRecords: records.length - seen.size,
    representatives,
    counts,
    largestGroups,
  };
}

class UnionFind {
  readonly parent: Int32Array;
  private readonly rank: Uint8Array;

  constructor(size: number) {
    this.parent = new Int32Array(size);
    this.rank = new Uint8Array(size);
    for (let index = 0; index < size; index += 1) this.parent[index] = index;
  }

  find(value: number): number {
    let root = value;
    while (this.parent[root] !== root) root = this.parent[root];
    while (this.parent[value] !== value) {
      const next = this.parent[value];
      this.parent[value] = root;
      value = next;
    }
    return root;
  }

  union(left: number, right: number) {
    let a = this.find(left);
    let b = this.find(right);
    if (a === b) return;
    if (this.rank[a] < this.rank[b]) [a, b] = [b, a];
    this.parent[b] = a;
    if (this.rank[a] === this.rank[b]) this.rank[a] += 1;
  }
}

function blocks(value: string, count: number): Array<{ index: number; value: string }> {
  const result: Array<{ index: number; value: string }> = [];
  for (let index = 0; index < count; index += 1) {
    const start = Math.floor(index * value.length / count);
    const end = Math.floor((index + 1) * value.length / count);
    result.push({ index, value: value.slice(start, end) });
  }
  return result;
}

function compatibleRecords(left: PostAnalysisRecord, right: PostAnalysisRecord, options: LineageOptions): boolean {
  if (options.requireSameLocus && left.locus !== right.locus) return false;
  return callsCompatible(left.vCall, right.vCall, options.callResolution, options.ambiguity) &&
    callsCompatible(left.jCall, right.jCall, options.callResolution, options.ambiguity);
}

function recordIndexTokens(record: PostAnalysisRecord, options: LineageOptions): string[] {
  const v = callSet(record.vCall, options.callResolution, options.ambiguity);
  const j = callSet(record.jCall, options.callResolution, options.ambiguity);
  if (!v.length || !j.length) return [];
  if (options.ambiguity === "strict") return [`${v.join("+")}\u0001${j.join("+")}`];
  return v.flatMap((vCall) => j.map((jCall) => `${vCall}\u0001${jCall}`));
}

export function assignLineages(
  records: PostAnalysisRecord[],
  options: LineageOptions,
  dedup?: DedupResult,
  activeMask?: Uint8Array,
): LineageResult {
  const union = new UnionFind(records.length);
  const bucket = new Map<string, number[]>();
  const exact = new Map<string, number>();
  let candidateComparisons = 0;
  let truncatedCandidates = 0;
  let assignedRecords = 0;
  let activeAbundance = 0;

  for (let index = 0; index < records.length; index += 1) {
    const record = records[index];
    const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : 1;
    if (!weight) continue;
    activeAbundance += weight;
    const cdr3 = normalizeNt(record.cdr3Nt);
    const tokens = recordIndexTokens(record, options);
    if (!cdr3 || !tokens.length || (options.productiveOnly && !record.productive)) continue;
    assignedRecords += weight;
    const distanceLimit = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
    const prefix = `${options.requireSameLocus ? record.locus : "*"}\u0000${cdr3.length}\u0000`;
    const exactSignature = `${prefix}${tokens.join("|")}\u0000${cdr3}`;
    const exactRepresentative = exact.get(exactSignature);
    if (exactRepresentative !== undefined) {
      union.union(index, exactRepresentative);
      continue;
    }
    exact.set(exactSignature, index);

    const candidates = new Set<number>();
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
        for (const candidate of bucket.get(key) ?? []) {
          candidates.add(candidate);
          if (candidates.size >= options.maxCandidateComparisons) break;
        }
        if (candidates.size >= options.maxCandidateComparisons) break;
      }
      if (candidates.size >= options.maxCandidateComparisons) break;
    }
    if (candidates.size >= options.maxCandidateComparisons) truncatedCandidates += 1;
    for (const candidate of candidates) {
      candidateComparisons += 1;
      const other = records[candidate];
      if (!compatibleRecords(record, other, options)) continue;
      if (hammingDistanceWithin(cdr3, normalizeNt(other.cdr3Nt), distanceLimit, false) <= distanceLimit) {
        union.union(index, candidate);
      }
    }
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
        const values = bucket.get(key);
        if (values) values.push(index);
        else bucket.set(key, [index]);
      }
    }
  }

  const eligible = new Uint8Array(records.length);
  const uniqueByRoot = new Uint32Array(records.length);
  const abundanceByRoot = new Float64Array(records.length);
  const representativeByRoot = new Int32Array(records.length);
  representativeByRoot.fill(-1);
  for (let index = 0; index < records.length; index += 1) {
    const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : 1;
    if (!weight) continue;
    const record = records[index];
    if (!record.cdr3Nt || !recordIndexTokens(record, options).length || (options.productiveOnly && !record.productive)) continue;
    const root = union.find(index);
    eligible[index] = 1;
    uniqueByRoot[root] += 1;
    abundanceByRoot[root] += weight;
    if (representativeByRoot[root] < 0) representativeByRoot[root] = index;
  }

  const assignments = new Int32Array(records.length);
  const rootIds = new Int32Array(records.length);
  const top: LineageSummary[] = [];
  const sizeBins = new Uint32Array(7);
  const vUsageMap = new Map<string, { lineages: number; abundance: number }>();
  const jUsageMap = new Map<string, { lineages: number; abundance: number }>();
  let lineageCount = 0;
  const smaller = (left: LineageSummary, right: LineageSummary) =>
    left.abundance < right.abundance ||
    (left.abundance === right.abundance && left.uniqueMembers < right.uniqueMembers) ||
    (left.abundance === right.abundance && left.uniqueMembers === right.uniqueMembers && left.id > right.id);
  const retainTop = (summary: LineageSummary) => {
    if (top.length < 1_000) {
      top.push(summary);
      let index = top.length - 1;
      while (index > 0) {
        const parent = Math.floor((index - 1) / 2);
        if (!smaller(top[index], top[parent])) break;
        [top[index], top[parent]] = [top[parent], top[index]];
        index = parent;
      }
      return;
    }
    if (!smaller(top[0], summary)) return;
    top[0] = summary;
    let index = 0;
    while (true) {
      const left = index * 2 + 1;
      const right = left + 1;
      let smallest = index;
      if (left < top.length && smaller(top[left], top[smallest])) smallest = left;
      if (right < top.length && smaller(top[right], top[smallest])) smallest = right;
      if (smallest === index) break;
      [top[index], top[smallest]] = [top[smallest], top[index]];
      index = smallest;
    }
  };
  const addUsage = (map: Map<string, { lineages: number; abundance: number }>, call: string, abundance: number) => {
    const value = map.get(call);
    if (value) {
      value.lineages += 1;
      value.abundance += abundance;
    } else map.set(call, { lineages: 1, abundance });
  };
  for (let root = 0; root < records.length; root += 1) {
    if (!uniqueByRoot[root]) continue;
    const id = ++lineageCount;
    rootIds[root] = id;
    const representativeOrdinal = representativeByRoot[root];
    const representative = records[representativeOrdinal];
    const abundance = abundanceByRoot[root];
    const bin = abundance === 1 ? 0 : abundance <= 3 ? 1 : abundance <= 9 ? 2 : abundance <= 24 ? 3 : abundance <= 99 ? 4 : abundance <= 499 ? 5 : 6;
    sizeBins[bin] += 1;
    const vCalls = callSet(representative.vCall, options.callResolution, options.ambiguity);
    const jCalls = callSet(representative.jCall, options.callResolution, options.ambiguity);
    vCalls.forEach((call) => addUsage(vUsageMap, call, abundance));
    jCalls.forEach((call) => addUsage(jUsageMap, call, abundance));
    retainTop({
      id,
      representativeOrdinal,
      uniqueMembers: uniqueByRoot[root],
      abundance,
      locus: representative.locus,
      vCalls,
      jCalls,
      cdr3Length: representative.cdr3Nt.length,
    });
  }
  for (let index = 0; index < records.length; index += 1) {
    if (eligible[index]) assignments[index] = rootIds[union.find(index)];
  }
  if (dedup) {
    for (let index = 0; index < dedup.representatives.length; index += 1) {
      assignments[index] = assignments[dedup.representatives[index]];
    }
  }
  const usage = (map: Map<string, { lineages: number; abundance: number }>) => [...map.entries()]
    .map(([call, value]) => ({ call, ...value }))
    .sort((a, b) => b.abundance - a.abundance || b.lineages - a.lineages || a.call.localeCompare(b.call));

  return {
    assignments,
    summaries: top.sort((a, b) => b.abundance - a.abundance || b.uniqueMembers - a.uniqueMembers || a.id - b.id),
    lineageCount,
    sizeHistogram: ["1", "2–3", "4–9", "10–24", "25–99", "100–499", "500+"].map((label, index) => ({ label, count: sizeBins[index] })),
    vUsage: usage(vUsageMap),
    jUsage: usage(jUsageMap),
    assignedRecords,
    unassignedRecords: activeAbundance - assignedRecords,
    candidateComparisons,
    truncatedCandidates,
  };
}

function queryValue(record: PostAnalysisRecord, target: QueryTarget): string {
  if (target === "cdr3_nt") return normalizeNt(record.cdr3Nt);
  if (target === "cdr3_aa") return normalizeAa(record.cdr3Aa);
  return "";
}

function matchesOneQueryConstraint(record: PostAnalysisRecord, constraint: QueryConstraint | undefined, options: QueryOptions): boolean {
  if (constraint?.locus && record.locus !== constraint.locus) return false;
  if (constraint?.vCall && !callsCompatible(record.vCall, constraint.vCall, options.callResolution, options.ambiguity)) return false;
  if (constraint?.jCall && !callsCompatible(record.jCall, constraint.jCall, options.callResolution, options.ambiguity)) return false;
  return true;
}

function matchesQueryConstraints(record: PostAnalysisRecord, options: QueryOptions): boolean {
  if (options.productiveOnly && !record.productive) return false;
  return matchesOneQueryConstraint(record, options, options);
}

export function queryRecords(
  records: PostAnalysisRecord[],
  queries: string[],
  options: QueryOptions,
  packedSketches?: Uint32Array,
  activeMask?: Uint8Array,
): QueryHit[] {
  const normalizedQueries = queries.map((query) => options.target === "cdr3_aa" ? normalizeAa(query) : normalizeNt(query));
  const querySketches = options.target === "trimmed" ? normalizedQueries.map((query) => minHashSketch(query)) : [];
  const hits: QueryHit[] = [];
  for (let recordIndex = 0; recordIndex < records.length; recordIndex += 1) {
    if (activeMask && !activeMask[recordIndex]) continue;
    const record = records[recordIndex];
    if (!matchesQueryConstraints(record, options)) continue;
    for (let queryIndex = 0; queryIndex < normalizedQueries.length; queryIndex += 1) {
      if (!matchesOneQueryConstraint(record, options.queryConstraints?.[queryIndex], options)) continue;
      const query = normalizedQueries[queryIndex];
      if (!query) continue;
      if (options.target === "trimmed") {
        const sketch = record.trimmedSketch ?? packedSketches?.subarray(record.ordinal * HASH_SEEDS.length, (record.ordinal + 1) * HASH_SEEDS.length);
        if (!sketch) continue;
        const score = sketchSimilarity(sketch, querySketches[queryIndex]);
        if (score >= options.identity) hits.push({ ordinal: record.ordinal, queryIndex, score, distance: Math.round((1 - score) * 1000), matched: "VDJ k-mer sketch" });
        continue;
      }
      const target = queryValue(record, options.target);
      if (!target) continue;
      let distance = Number.POSITIVE_INFINITY;
      let score = 0;
      if (options.metric === "exact") {
        if (target !== query) continue;
        distance = 0;
        score = 1;
      } else if (options.metric === "substring") {
        if (!target.includes(query) && !query.includes(target)) continue;
        distance = Math.abs(target.length - query.length);
        score = Math.min(target.length, query.length) / Math.max(target.length, query.length);
      } else if (options.metric === "hamming") {
        if (target.length !== query.length) continue;
        const maximum = Math.floor((1 - options.identity) * target.length + 1e-9);
        distance = hammingDistanceWithin(target, query, maximum, options.target !== "cdr3_aa");
        if (distance > maximum) continue;
        score = 1 - distance / target.length;
      } else {
        const maximum = Math.floor((1 - options.identity) * Math.max(target.length, query.length) + 1e-9);
        distance = bandedEditDistance(target, query, maximum);
        if (distance > maximum) continue;
        score = 1 - distance / Math.max(target.length, query.length);
      }
      hits.push({ ordinal: record.ordinal, queryIndex, score, distance, matched: target });
    }
  }
  return hits.sort((a, b) => b.score - a.score || a.distance - b.distance || a.ordinal - b.ordinal).slice(0, options.maxResults);
}

export interface ExpansionOptions extends LineageOptions {
  maxResults: number;
}

export function expandSingleLinkage(
  records: PostAnalysisRecord[],
  seedOrdinals: number[],
  options: ExpansionOptions,
  activeMask?: Uint8Array,
): { ordinals: number[]; comparisons: number; capped: boolean } {
  const byLength = new Map<number, Map<string, number[]>>();
  for (let index = 0; index < records.length; index += 1) {
    if (activeMask && !activeMask[index]) continue;
    const record = records[index];
    const cdr3 = normalizeNt(record.cdr3Nt);
    if (!cdr3 || (options.productiveOnly && !record.productive)) continue;
    const maximum = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, maximum + 1));
    let indexMap = byLength.get(cdr3.length);
    if (!indexMap) {
      indexMap = new Map();
      byLength.set(cdr3.length, indexMap);
    }
    for (const block of blocks(cdr3, blockCount)) {
      const key = `${block.index}\u0000${block.value}`;
      const values = indexMap.get(key);
      if (values) values.push(index);
      else indexMap.set(key, [index]);
    }
  }

  const visited = new Set(seedOrdinals.filter((ordinal) => ordinal >= 0 && ordinal < records.length && (!activeMask || Boolean(activeMask[ordinal]))));
  const queue = [...visited];
  let comparisons = 0;
  for (let head = 0; head < queue.length && visited.size < options.maxResults; head += 1) {
    const current = queue[head];
    const source = records[current];
    const cdr3 = normalizeNt(source.cdr3Nt);
    const maximum = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, maximum + 1));
    const candidates = new Set<number>();
    const indexMap = byLength.get(cdr3.length);
    for (const block of blocks(cdr3, blockCount)) {
      for (const candidate of indexMap?.get(`${block.index}\u0000${block.value}`) ?? []) candidates.add(candidate);
    }
    for (const candidate of candidates) {
      if (visited.has(candidate)) continue;
      comparisons += 1;
      if (!compatibleRecords(source, records[candidate], options)) continue;
      if (hammingDistanceWithin(cdr3, normalizeNt(records[candidate].cdr3Nt), maximum, false) <= maximum) {
        visited.add(candidate);
        queue.push(candidate);
        if (visited.size >= options.maxResults) break;
      }
    }
  }
  return { ordinals: [...visited], comparisons, capped: visited.size >= options.maxResults };
}

export interface FastaRecord {
  name: string;
  sequence: string;
}

export function parseFasta(text: string, aligned = false): FastaRecord[] {
  const records: FastaRecord[] = [];
  let name = "";
  let sequence: string[] = [];
  const commit = () => {
    if (!name) return;
    const joined = sequence.join("").toUpperCase().replaceAll(".", "-").replace(/\s/g, "");
    const normalized = aligned ? joined.replace(/[^ACGTNRYKMSWBDHV-]/g, "N") : normalizeNt(joined);
    records.push({ name, sequence: normalized });
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      commit();
      name = line.slice(1).trim().split(/\s+/, 1)[0] ?? "";
      sequence = [];
    } else if (name) {
      sequence.push(line);
    }
  }
  commit();
  return records;
}

export interface ReferenceMsa {
  names: string[];
  sequences: string[];
  lookup: Map<string, number>;
  length: number;
}

function nameAliases(name: string): string[] {
  const fields = name.split("|").filter(Boolean);
  const gene = fields.find((field) => /(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)[VJ]/i.test(field));
  return [...new Set([name, gene ?? "", name.split(/\s+/, 1)[0] ?? ""].filter(Boolean))];
}

export function prepareReferenceMsa(text: string): ReferenceMsa {
  const records = parseFasta(text, true);
  if (records.length < 2) throw new Error("CHMMAIRRa requires at least two aligned reference sequences.");
  const length = records[0].sequence.length;
  if (!length || records.some((record) => record.sequence.length !== length)) {
    throw new Error("Reference MSA records must all have the same aligned length.");
  }
  const lookup = new Map<string, number>();
  records.forEach((record, index) => nameAliases(record.name).forEach((alias) => lookup.set(alias, index)));
  return { names: records.map((record) => record.name), sequences: records.map((record) => record.sequence), lookup, length };
}

function topCallName(value: string): string {
  return value.split(",", 1)[0]?.trim() ?? "";
}

export function threadSequenceToMsa(
  sequenceAlignment: string,
  germlineAlignment: string,
  call: string,
  msa: ReferenceMsa,
): string {
  const callName = topCallName(call);
  let referenceIndex = msa.lookup.get(callName);
  if (referenceIndex === undefined) {
    const gene = callName.replace(/\*.*$/, "");
    const match = [...msa.lookup.entries()].find(([name]) => name.replace(/\*.*$/, "") === gene);
    referenceIndex = match?.[1];
  }
  if (referenceIndex === undefined) throw new Error(`Reference MSA is missing ${callName || "the assigned allele"}.`);
  const gappedFullReference = msa.sequences[referenceIndex];
  const degappedReference = gappedFullReference.replaceAll("-", "");
  const localGermline = germlineAlignment.toUpperCase().replaceAll(".", "-");
  // CHMMera maps every non-ACGT observation other than a gap to its
  // non-informative N state. Apply that conversion before threading so IUPAC
  // query symbols cannot become accidental informative matches.
  const localQuery = sequenceAlignment.toUpperCase().replaceAll(".", "-").replace(/[^ACGT-]/g, "N");
  const degappedLocal = localGermline.replaceAll("-", "");
  const matchStart = degappedReference.indexOf(degappedLocal);
  if (matchStart < 0) throw new Error(`The ${callName} local germline alignment cannot be located in its MSA reference.`);
  const withoutInsertions: string[] = [];
  for (let index = 0; index < localGermline.length; index += 1) {
    if (localGermline[index] !== "-") withoutInsertions.push(localQuery[index] ?? "N");
  }
  const query = [..."N".repeat(matchStart), ...withoutInsertions];
  const threaded = new Array<string>(msa.length).fill("-");
  let position = 0;
  for (let index = 0; index < msa.length; index += 1) {
    if (position >= query.length) break;
    if (gappedFullReference[index] !== "-") threaded[index] = query[position++] ?? "N";
  }
  return threaded.join("");
}

export interface ChmmOptions {
  method: "BW" | "DB";
  priorProbability: number;
  baseMutationProbability: number;
  mutationRates: number[];
  mutationSwitchProbability: number;
  detailed: boolean;
  tracePath?: boolean;
}

export interface ChmmResult {
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  referencePath?: Int32Array;
}

export function chmmairraDistanceFromReference(query: string, germline: string): number {
  const normalizedQuery = query.toUpperCase().replaceAll(".", "-");
  const normalizedGermline = germline.toUpperCase().replaceAll(".", "-");
  let mismatches = 0;
  // This intentionally follows CHMMAIRRa.jl's add_DFR_column!: zip the local
  // alignments and count every unequal character, including gap and N sites.
  for (let index = 0; index < Math.min(normalizedQuery.length, normalizedGermline.length); index += 1) {
    if (normalizedQuery[index] !== normalizedGermline[index]) mismatches += 1;
  }
  return mismatches;
}

function emission(reference: string, observation: string, mutation: number): number {
  if (observation === "-" || observation === "N") return 1;
  return reference === observation ? 1 - mutation : mutation / 3;
}

function approximateEmissions(msa: ReferenceMsa, observation: string, mutations: Float64Array): Float64Array {
  const values = new Float64Array(msa.sequences.length * msa.length);
  for (let state = 0; state < msa.sequences.length; state += 1) {
    for (let site = 0; site < msa.length; site += 1) {
      values[state * msa.length + site] = emission(msa.sequences[state][site], observation[site], mutations[state]);
    }
  }
  return values;
}

export function estimateMutationProbabilities(msa: ReferenceMsa, observation: string, prior: number, baseMutation: number): Float64Array {
  const states = msa.sequences.length;
  const length = msa.length;
  const mutations = new Float64Array(states);
  mutations.fill(baseMutation);
  if (states < 2 || length < 1) return mutations;
  const emissions = approximateEmissions(msa, observation, mutations);
  const alpha = new Float64Array(states * length);
  const beta = new Float64Array(states * length);
  const scale = new Float64Array(length);
  scale[0] = 1;
  const switchProbability = prior / length;
  const same = 1 - switchProbability;
  const different = switchProbability / (states - 1);
  for (let state = 0; state < states; state += 1) alpha[state * length] = emissions[state * length] / states;
  for (let site = 0; site < length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += alpha[state * length + site];
    let nextTotal = 0;
    for (let state = 0; state < states; state += 1) {
      const value = ((total - alpha[state * length + site]) * different + alpha[state * length + site] * same) * emissions[state * length + site + 1];
      alpha[state * length + site + 1] = value;
      nextTotal += value;
    }
    scale[site + 1] = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) alpha[state * length + site + 1] *= scale[site + 1];
  }
  for (let state = 0; state < states; state += 1) beta[state * length + length - 1] = scale[length - 1];
  for (let site = length - 2; site >= 0; site -= 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += beta[state * length + site + 1] * emissions[state * length + site + 1];
    for (let state = 0; state < states; state += 1) {
      const own = beta[state * length + site + 1] * emissions[state * length + site + 1];
      beta[state * length + site] = (different * (total - own) + same * own) * scale[site];
    }
  }
  const mutated = new Float64Array(states);
  const sameCount = new Float64Array(states);
  mutated.fill(2);
  sameCount.fill(10);
  for (let site = 0; site < length; site += 1) {
    if (observation[site] === "-" || observation[site] === "N") continue;
    let normalization = 0;
    for (let state = 0; state < states; state += 1) normalization += alpha[state * length + site] * beta[state * length + site];
    if (!normalization) continue;
    for (let state = 0; state < states; state += 1) {
      const posterior = alpha[state * length + site] * beta[state * length + site] / normalization;
      if (observation[site] === msa.sequences[state][site]) sameCount[state] += posterior;
      else mutated[state] += posterior;
    }
  }
  for (let state = 0; state < states; state += 1) mutations[state] = mutated[state] / (mutated[state] + sameCount[state]);
  return mutations;
}

function approximateProbability(msa: ReferenceMsa, observation: string, prior: number, mutations: Float64Array): number {
  const states = msa.sequences.length;
  if (states < 2) return 0;
  const emissions = approximateEmissions(msa, observation, mutations);
  const nonChimeric = new Float64Array(states);
  const chimeric = new Float64Array(states);
  const switchProbability = prior / msa.length;
  const same = 1 - switchProbability;
  const different = switchProbability / (states - 1);
  for (let state = 0; state < states; state += 1) nonChimeric[state] = emissions[state * msa.length] / states;
  for (let site = 0; site < msa.length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
    let nextTotal = 0;
    for (let state = 0; state < states; state += 1) {
      const emissionValue = emissions[state * msa.length + site + 1];
      const own = nonChimeric[state] + chimeric[state];
      chimeric[state] = ((total - own) * different + chimeric[state] * same) * emissionValue;
      nonChimeric[state] = nonChimeric[state] * same * emissionValue;
      nextTotal += chimeric[state] + nonChimeric[state];
    }
    const scaling = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) {
      chimeric[state] *= scaling;
      nonChimeric[state] *= scaling;
    }
  }
  let chimericTotal = 0;
  let total = 0;
  for (let state = 0; state < states; state += 1) {
    chimericTotal += chimeric[state];
    total += chimeric[state] + nonChimeric[state];
  }
  return total ? chimericTotal / total : 0;
}

function fullProbability(msa: ReferenceMsa, observation: string, options: ChmmOptions): number {
  const references = msa.sequences.length;
  const rates = options.mutationRates.length ? options.mutationRates : [0.005];
  const rateCount = rates.length;
  const states = references * rateCount;
  if (references < 2) return 0;
  const nonChimeric = new Float64Array(states);
  const chimeric = new Float64Array(states);
  const nextNonChimeric = new Float64Array(states);
  const nextChimeric = new Float64Array(states);
  const switchProbability = options.priorProbability / msa.length;
  const self = 1 - switchProbability - options.mutationSwitchProbability;
  const differentReference = switchProbability / ((references - 1) * rateCount);
  const differentMutation = rateCount === 1 ? 0 : options.mutationSwitchProbability / (rateCount - 1);
  const emissionAt = (state: number, site: number) => {
    const reference = Math.floor(state / rateCount);
    const rate = state % rateCount;
    return emission(msa.sequences[reference][site], observation[site], rates[rate]);
  };
  for (let state = 0; state < states; state += 1) nonChimeric[state] = emissionAt(state, 0) / states;
  for (let site = 0; site < msa.length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
    let nextTotal = 0;
    for (let reference = 0; reference < references; reference += 1) {
      let nonReference = 0;
      let chimReference = 0;
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        nonReference += nonChimeric[state];
        chimReference += chimeric[state];
      }
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        const emit = emissionAt(state, site + 1);
        nextNonChimeric[state] = (nonChimeric[state] * self + (nonReference - nonChimeric[state]) * differentMutation) * emit;
        nextChimeric[state] = (chimeric[state] * self + (chimReference - chimeric[state]) * differentMutation + (total - chimReference - nonReference) * differentReference) * emit;
        nextTotal += nextNonChimeric[state] + nextChimeric[state];
      }
    }
    const scaling = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) {
      nonChimeric[state] = nextNonChimeric[state] * scaling;
      chimeric[state] = nextChimeric[state] * scaling;
    }
  }
  let chimericTotal = 0;
  let total = 0;
  for (let state = 0; state < states; state += 1) {
    chimericTotal += chimeric[state];
    total += chimeric[state] + nonChimeric[state];
  }
  return total ? chimericTotal / total : 0;
}

function fullViterbi(msa: ReferenceMsa, observation: string, options: ChmmOptions): Pick<ChmmResult, "startingReference" | "recombinations" | "referencePath"> {
  const references = msa.sequences.length;
  const rates = options.mutationRates.length ? options.mutationRates : [0.005];
  const rateCount = rates.length;
  const states = references * rateCount;
  const length = msa.length;
  if (references < 2) return {
    startingReference: msa.names[0] ?? "",
    recombinations: [],
    referencePath: options.tracePath ? new Int32Array(length) : undefined,
  };
  const previous = new Float64Array(states);
  const next = new Float64Array(states);
  const from = new Int32Array(states * length);
  const switchProbability = options.priorProbability / length;
  const logSelf = Math.log(1 - switchProbability - options.mutationSwitchProbability);
  const logDifferentReference = Math.log(switchProbability / ((references - 1) * rateCount));
  const logDifferentMutation = rateCount === 1 ? Number.NEGATIVE_INFINITY : Math.log(options.mutationSwitchProbability / (rateCount - 1));
  const emissionAt = (state: number, site: number) => emission(msa.sequences[Math.floor(state / rateCount)][site], observation[site], rates[state % rateCount]);
  for (let state = 0; state < states; state += 1) {
    previous[state] = Math.log(emissionAt(state, 0) / states);
    from[state * length] = state;
  }
  for (let site = 0; site < length - 1; site += 1) {
    const bestInReference = new Int32Array(references);
    for (let reference = 0; reference < references; reference += 1) {
      let best = reference * rateCount;
      for (let rate = 1; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        if (previous[state] > previous[best]) best = state;
      }
      bestInReference[reference] = best;
    }
    let bestReference = 0;
    let secondReference = references > 1 ? 1 : 0;
    if (previous[bestInReference[secondReference]] > previous[bestInReference[bestReference]]) [bestReference, secondReference] = [secondReference, bestReference];
    for (let reference = 2; reference < references; reference += 1) {
      if (previous[bestInReference[reference]] > previous[bestInReference[bestReference]]) {
        secondReference = bestReference;
        bestReference = reference;
      } else if (previous[bestInReference[reference]] > previous[bestInReference[secondReference]]) secondReference = reference;
    }
    for (let reference = 0; reference < references; reference += 1) {
      const outside = bestInReference[reference === bestReference ? secondReference : bestReference];
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        let selected = state;
        let score = previous[state] + logSelf;
        const referenceScore = previous[outside] + logDifferentReference;
        if (referenceScore > score) {
          selected = outside;
          score = referenceScore;
        }
        if (rateCount > 1) {
          let otherRate = reference * rateCount + (rate === 0 ? 1 : 0);
          for (let candidateRate = 0; candidateRate < rateCount; candidateRate += 1) {
            const candidate = reference * rateCount + candidateRate;
            if (candidate !== state && previous[candidate] > previous[otherRate]) otherRate = candidate;
          }
          const mutationScore = previous[otherRate] + logDifferentMutation;
          if (mutationScore > score) {
            selected = otherRate;
            score = mutationScore;
          }
        }
        from[state * length + site + 1] = selected;
        next[state] = score + Math.log(emissionAt(state, site + 1));
      }
    }
    previous.set(next);
  }
  let current = 0;
  for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
  const recombinations: Array<{ position: number; left: string; right: string }> = [];
  const referencePath = options.tracePath ? new Int32Array(length) : undefined;
  if (referencePath) referencePath[length - 1] = Math.floor(current / rateCount);
  for (let site = length - 1; site >= 1; site -= 1) {
    const parent = from[current * length + site];
    const leftReference = Math.floor(parent / rateCount);
    const rightReference = Math.floor(current / rateCount);
    if (leftReference !== rightReference) recombinations.push({ position: site + 1, left: msa.names[leftReference], right: msa.names[rightReference] });
    current = parent;
    if (referencePath) referencePath[site - 1] = Math.floor(current / rateCount);
  }
  return { startingReference: msa.names[Math.floor(current / rateCount)], recombinations: recombinations.reverse(), referencePath };
}

function approximateViterbi(msa: ReferenceMsa, observation: string, prior: number, mutations: Float64Array, tracePath = false): Pick<ChmmResult, "startingReference" | "recombinations" | "referencePath"> {
  const states = msa.sequences.length;
  const length = msa.length;
  const emissions = approximateEmissions(msa, observation, mutations);
  const previous = new Float64Array(states);
  const next = new Float64Array(states);
  const from = new Int32Array(states * length);
  const same = Math.log(1 - prior / length);
  const different = Math.log((prior / length) / (states - 1));
  for (let state = 0; state < states; state += 1) {
    previous[state] = Math.log(emissions[state * length] / states);
    from[state * length] = state;
  }
  for (let site = 0; site < length - 1; site += 1) {
    let best = 0;
    for (let state = 1; state < states; state += 1) if (previous[state] > previous[best]) best = state;
    for (let state = 0; state < states; state += 1) {
      const stay = previous[state] + same;
      const change = previous[best] + different;
      const remain = state === best || stay > change;
      from[state * length + site + 1] = remain ? state : best;
      next[state] = (remain ? stay : change) + Math.log(emissions[state * length + site + 1]);
    }
    previous.set(next);
  }
  let current = 0;
  for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
  const reverse: Array<{ position: number; left: string; right: string }> = [];
  const referencePath = tracePath ? new Int32Array(length) : undefined;
  if (referencePath) referencePath[length - 1] = current;
  for (let site = length - 1; site >= 1; site -= 1) {
    const parent = from[current * length + site];
    if (parent !== current) reverse.push({ position: site + 1, left: msa.names[parent], right: msa.names[current] });
    current = parent;
    if (referencePath) referencePath[site - 1] = current;
  }
  return { startingReference: msa.names[current], recombinations: reverse.reverse(), referencePath };
}

export function runChmm(
  msa: ReferenceMsa,
  threadedObservation: string,
  localSequenceAlignment: string,
  localGermlineAlignment: string,
  options: ChmmOptions,
): ChmmResult {
  if (threadedObservation.length !== msa.length) throw new Error("Threaded query and reference MSA lengths differ.");
  const dfr = chmmairraDistanceFromReference(localSequenceAlignment, localGermlineAlignment);
  if (options.method === "DB") {
    const path = options.detailed ? fullViterbi(msa, threadedObservation, options) : { startingReference: "", recombinations: [] };
    return {
      probability: fullProbability(msa, threadedObservation, options),
      dfr,
      ...path,
    };
  }
  const mutations = estimateMutationProbabilities(msa, threadedObservation, options.priorProbability, options.baseMutationProbability);
  const probability = approximateProbability(msa, threadedObservation, options.priorProbability, mutations);
  const path = options.detailed ? approximateViterbi(msa, threadedObservation, options.priorProbability, mutations, options.tracePath) : { startingReference: "", recombinations: [] };
  return { probability, dfr, ...path };
}
