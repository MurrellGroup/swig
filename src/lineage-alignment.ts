import { parseFasta } from "./post-analysis-core.ts";
import { biologicalFrameOffset } from "./alignment-model.ts";
import type { AirrDetailRow } from "./result-store.ts";

export const GERMLINE_OUTGROUP = "__germline_N_masked__";

export type LineageGermlineMethod = "closest" | "consensus";

function safeName(value: string, fallback: string): string {
  return (value || fallback).replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
}

function normalizedAlignment(value: string, fallback = ""): string {
  return (value || fallback).toUpperCase().replaceAll(".", "-").replaceAll("U", "T").replace(/[^ACGTN-]/g, "N");
}

export interface AnchoredLineageRow {
  row: AirrDetailRow;
  name: string;
  sequence: string;
  germline: string;
  frame: number;
  /** A sparse Double-D call was present on the AIRR detail row. */
  doubleDPositive: boolean;
  /** Both D alignments were safely projected into the combined AIRR columns. */
  doubleDGermlineApplied: boolean;
  dCall?: string;
  d2Call?: string;
}

export interface InferredLineageGermline {
  method: LineageGermlineMethod;
  /** V/D/J-aware root template; unresolved recombination columns remain N. */
  template: string;
  /** Template with unresolved N columns filled by the selected reconstruction method. */
  uca: string;
  /** Supported endpoint-trimmed forms used for lineage-neighbour comparison. */
  trimmedTemplate: string;
  trimmedUca: string;
  rowsUsed: number;
  columns: number;
  startColumn: number;
  endColumn: number;
  knownColumns: number;
  inferredColumns: number;
  conflictingColumns: number;
  minimumEndpointSupport: number;
  /** Rows carrying a sparse supported D2 call, whether or not projection succeeded. */
  doubleDPositiveRows: number;
  /** Double-D rows whose D1 and D2 germlines were projected into the root template. */
  doubleDResolvedRows: number;
  /** Positive rows omitted from VDDJ reconstruction because their coordinates were incomplete. */
  doubleDIncompleteRows: number;
  /** True when the returned root was reconstructed from V-D1-D2-J-aware rows. */
  doubleDTemplate: boolean;
  /** Number of VDDJ-aware rows used by the selected reconstruction method. */
  doubleDRowsUsed: number;
  /** AIRR member used by closest-member reconstruction. */
  selectedOrdinal?: number;
  selectedSequenceId?: string;
  selectedDCall?: string;
  selectedD2Call?: string;
  /** Equal-weight V/J identity across informative matched segment columns. */
  selectedVjIdentity?: number;
  selectedVIdentity?: number;
  selectedJIdentity?: number;
  selectedComparedColumns?: number;
}

interface DoubleDGermlineProjection {
  germline: string;
  positive: boolean;
  applied: boolean;
  dCall?: string;
  d2Call?: string;
}

function positiveInteger(value: unknown): number | null {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 && parsed <= 10_000_000 ? parsed : null;
}

/**
 * Replace the baseline single-D portion of a combined AIRR germline alignment
 * with the independently screened D1 and D2 alignments. Sequence coordinates
 * are projected through the gapped combined query, so V/J indels do not shift
 * either D. The current Double-D screener emits ungapped D alignments; the gap
 * branch below also accepts a future gapped sidecar when the corresponding
 * combined AIRR gap column exists.
 */
function projectDoubleDGermline(
  values: Record<string, string>,
  rawSequence: string,
  rawGermline: string,
): DoubleDGermlineProjection {
  const dCall = String(values.d_call || "").trim();
  const d2Call = String(values.d2_call || "").trim();
  if (!d2Call) return { germline: rawGermline, positive: false, applied: false };

  const unresolved = (): DoubleDGermlineProjection => ({
    germline: rawGermline,
    positive: true,
    applied: false,
    dCall: dCall || undefined,
    d2Call,
  });
  // Absolute AIRR sequence coordinates can be mapped only against the actual
  // combined alignment. Falling back to an unaligned full-length `sequence`
  // would silently shift every D column when V begins internally.
  if (!String(values.sequence_alignment || "").trim() || !String(values.germline_alignment || "").trim()) return unresolved();
  const alignmentStart = positiveInteger(values.v_sequence_start);
  const vEnd = positiveInteger(values.v_sequence_end);
  const jStart = positiveInteger(values.j_sequence_start);
  const dStart = positiveInteger(values.d_sequence_start);
  const dEnd = positiveInteger(values.d_sequence_end);
  const d2Start = positiveInteger(values.d2_sequence_start);
  const d2End = positiveInteger(values.d2_sequence_end);
  if (alignmentStart === null || vEnd === null || jStart === null || dStart === null || dEnd === null ||
      d2Start === null || d2End === null || alignmentStart > vEnd || vEnd >= dStart || dStart > dEnd ||
      dEnd >= d2Start || d2Start > d2End || d2End >= jStart) return unresolved();

  const sequence = rawSequence.padEnd(Math.max(rawSequence.length, rawGermline.length), "-");
  const corrected = [...rawGermline.padEnd(sequence.length, "-")];
  const coordinateToColumn = new Map<number, number>();
  let queryCoordinate = alignmentStart;
  for (let column = 0; column < sequence.length; column += 1) {
    if (sequence[column] === "-") continue;
    coordinateToColumn.set(queryCoordinate, column);
    queryCoordinate += 1;
  }
  const vEndColumn = coordinateToColumn.get(vEnd);
  const jStartColumn = coordinateToColumn.get(jStart);
  if (vEndColumn === undefined || jStartColumn === undefined || vEndColumn >= jStartColumn) return unresolved();

  // Remove every baseline single-D reference state before adding the supported
  // D1/D2 pair. Query-bearing insertion sites remain N; gap columns remain gaps.
  for (let column = vEndColumn + 1; column < jStartColumn; column += 1) {
    corrected[column] = sequence[column] === "-" ? "-" : "N";
  }

  const occupiedGapColumns = new Set<number>();
  const projectSegment = (start: number, end: number, queryValue: string, germlineValue: string): boolean => {
    const alignedQuery = normalizedAlignment(queryValue);
    const alignedGermline = normalizedAlignment(germlineValue);
    if (!alignedQuery || !alignedGermline || alignedQuery.length !== alignedGermline.length) return false;
    let coordinate = start;
    let previousColumn = (coordinateToColumn.get(start) ?? vEndColumn + 1) - 1;
    for (let index = 0; index < alignedQuery.length; index += 1) {
      const queryBase = alignedQuery[index];
      const germlineBase = alignedGermline[index];
      if (queryBase !== "-") {
        const column = coordinateToColumn.get(coordinate);
        if (column === undefined) return false;
        const combinedBase = sequence[column];
        if (BASE_INDEX[queryBase] !== undefined && BASE_INDEX[combinedBase] !== undefined && queryBase !== combinedBase) return false;
        corrected[column] = germlineBase;
        previousColumn = column;
        coordinate += 1;
        continue;
      }

      // A reference insertion can be represented only if that gap column is
      // already present in the combined AIRR alignment shared by every row.
      const nextColumn = coordinateToColumn.get(coordinate) ?? jStartColumn;
      let gapColumn = -1;
      for (let column = previousColumn + 1; column < nextColumn; column += 1) {
        if (sequence[column] === "-" && !occupiedGapColumns.has(column)) {
          gapColumn = column;
          break;
        }
      }
      if (gapColumn < 0) return false;
      corrected[gapColumn] = germlineBase;
      occupiedGapColumns.add(gapColumn);
      previousColumn = gapColumn;
    }
    return coordinate === end + 1;
  };

  if (!projectSegment(dStart, dEnd, values.d_sequence_alignment, values.d_germline_alignment) ||
      !projectSegment(d2Start, d2End, values.d2_sequence_alignment, values.d2_germline_alignment)) return unresolved();
  return {
    germline: corrected.join(""),
    positive: true,
    applied: true,
    dCall: dCall || undefined,
    d2Call,
  };
}

/**
 * Put AIRR pairwise alignments onto a common V-reference origin. The AIRR
 * `v_germline_start` coordinate supplies left padding; all remaining columns
 * retain the caller's insertion/gap structure. This is deliberately a quick
 * reference-anchored view rather than a de-novo MSA.
 */
export function referenceAnchoredLineageRows(rows: AirrDetailRow[]): AnchoredLineageRow[] {
  if (!rows.length) throw new Error("The selected lineage has no retrievable AIRR records.");
  const provisional = rows.map((row, index) => {
    const rawSequence = normalizedAlignment(row.values.sequence_alignment, row.values.sequence);
    const rawGermline = normalizedAlignment(row.values.germline_alignment, rawSequence.replace(/[ACGT]/g, "N"));
    const doubleD = projectDoubleDGermline(row.values, rawSequence, rawGermline);
    const localColumns = Math.max(rawSequence.length, doubleD.germline.length);
    const reportedStart = Math.floor(Number(row.values.v_germline_start));
    const left = Number.isFinite(reportedStart) && reportedStart >= 1 && reportedStart <= 10_000 ? reportedStart - 1 : 0;
    return {
      row,
      name: `${safeName(row.values.sequence_id, `sequence_${index + 1}`)}__${row.record.ordinal + 1}`,
      sequence: "-".repeat(left) + rawSequence.padEnd(localColumns, "-"),
      germline: "-".repeat(left) + doubleD.germline.padEnd(localColumns, "-"),
      frame: biologicalFrameOffset(Number(row.values.v_sequence_start) || 1, Number(row.values.sequence_frame) || 1),
      doubleDPositive: doubleD.positive,
      doubleDGermlineApplied: doubleD.applied,
      dCall: doubleD.dCall,
      d2Call: doubleD.d2Call,
    };
  });
  const columns = Math.max(...provisional.map((record) => Math.max(record.sequence.length, record.germline.length)));
  return provisional.map((record) => ({
    ...record,
    sequence: record.sequence.padEnd(columns, "-"),
    germline: record.germline.padEnd(columns, "-"),
  }));
}

function plurality(counts: readonly number[]): { base: string; count: number; tied: boolean } {
  const bases = ["A", "C", "G", "T"];
  let maximum = 0;
  let maximumIndex = 0;
  let tied = false;
  for (let index = 0; index < counts.length; index += 1) {
    if (counts[index] > maximum) {
      maximum = counts[index];
      maximumIndex = index;
      tied = false;
    } else if (counts[index] === maximum && maximum > 0) tied = true;
  }
  return { base: bases[maximumIndex], count: maximum, tied };
}

const BASE_INDEX: Record<string, number | undefined> = { A: 0, C: 1, G: 2, T: 3 };

interface SegmentIdentity {
  identity: number | null;
  matches: number;
  compared: number;
}

function segmentIdentity(queryValue: string, germlineValue: string, reportedValue: string | number | null | undefined): SegmentIdentity {
  const query = normalizedAlignment(queryValue);
  const germline = normalizedAlignment(germlineValue);
  let matches = 0;
  let compared = 0;
  const columns = Math.max(query.length, germline.length);
  for (let column = 0; column < columns; column += 1) {
    const queryBase = query[column];
    const germlineBase = germline[column];
    if (BASE_INDEX[queryBase] === undefined || BASE_INDEX[germlineBase] === undefined) continue;
    compared += 1;
    if (queryBase === germlineBase) matches += 1;
  }
  if (compared) return { identity: matches / compared, matches, compared };
  const reported = Number(reportedValue);
  return Number.isFinite(reported) && reported >= 0 && reported <= 1
    ? { identity: reported, matches: 0, compared: 0 }
    : { identity: null, matches: 0, compared: 0 };
}

interface ClosestMemberScore {
  row: AirrDetailRow;
  anchored: AnchoredLineageRow;
  v: SegmentIdentity;
  j: SegmentIdentity;
  equalWeightIdentity: number;
  combinedIdentity: number;
  compared: number;
  informativeSegments: number;
}

function closestMember(anchored: AnchoredLineageRow[]): ClosestMemberScore {
  const scores = anchored.map((record) => {
    const values = record.row.values;
    const v = segmentIdentity(values.v_sequence_alignment, values.v_germline_alignment, values.v_identity ?? record.row.record.vIdentity);
    const j = segmentIdentity(values.j_sequence_alignment, values.j_germline_alignment, values.j_identity ?? record.row.record.jIdentity);
    const segmentIdentities = [v.identity, j.identity].filter((value): value is number => value !== null);
    let compared = v.compared + j.compared;
    let matches = v.matches + j.matches;
    let equalWeightIdentity = segmentIdentities.length ? segmentIdentities.reduce((sum, value) => sum + value, 0) / segmentIdentities.length : -1;
    let combinedIdentity = compared ? matches / compared : equalWeightIdentity;
    if (!segmentIdentities.length) {
      // Older imported AIRR tables may omit segment-specific alignment fields.
      // Fall back to known V/D/J columns in the combined AIRR alignment.
      const fallback = segmentIdentity(values.sequence_alignment || values.sequence, values.germline_alignment, null);
      compared = fallback.compared;
      matches = fallback.matches;
      equalWeightIdentity = fallback.identity ?? -1;
      combinedIdentity = fallback.identity ?? -1;
    }
    return { row: record.row, anchored: record, v, j, equalWeightIdentity, combinedIdentity, compared, informativeSegments: segmentIdentities.length };
  });
  scores.sort((left, right) =>
    right.informativeSegments - left.informativeSegments ||
    right.equalWeightIdentity - left.equalWeightIdentity ||
    right.combinedIdentity - left.combinedIdentity ||
    right.compared - left.compared ||
    left.row.record.ordinal - right.row.record.ordinal,
  );
  return scores[0];
}

function inferClosestLineageGermline(rows: AirrDetailRow[]): InferredLineageGermline {
  const anchored = referenceAnchoredLineageRows(rows);
  const doubleDResolved = anchored.filter((record) => record.doubleDGermlineApplied);
  // If this lineage contains a supported VDDJ architecture, do not let a
  // baseline single-D member silently replace it merely because its V/J ends
  // are marginally less mutated. Rank the VDDJ-aware members by the same rule.
  const selected = closestMember(doubleDResolved.length ? doubleDResolved : anchored);
  const template = [...selected.anchored.germline];
  const uca = [...selected.anchored.germline];
  const supported = template.map((base, column) => base !== "-" || selected.anchored.sequence[column] !== "-");
  let knownColumns = 0;
  let inferredColumns = 0;
  for (let column = 0; column < template.length; column += 1) {
    if (BASE_INDEX[template[column]] !== undefined) {
      knownColumns += 1;
      continue;
    }
    if (template[column] === "N" && BASE_INDEX[selected.anchored.sequence[column]] !== undefined) {
      uca[column] = selected.anchored.sequence[column];
      inferredColumns += 1;
    }
  }
  let start = supported.findIndex(Boolean);
  if (start < 0) start = 0;
  let end = supported.length;
  while (end > start && !supported[end - 1]) end -= 1;
  const selectedSequenceId = selected.row.values.sequence_id || selected.row.record.sequenceId;
  return {
    method: "closest",
    template: template.join(""),
    uca: uca.join(""),
    trimmedTemplate: template.slice(start, end).join(""),
    trimmedUca: uca.slice(start, end).join(""),
    rowsUsed: 1,
    columns: template.length,
    startColumn: start,
    endColumn: end,
    knownColumns,
    inferredColumns,
    conflictingColumns: 0,
    minimumEndpointSupport: 1,
    doubleDPositiveRows: anchored.filter((record) => record.doubleDPositive).length,
    doubleDResolvedRows: doubleDResolved.length,
    doubleDIncompleteRows: anchored.filter((record) => record.doubleDPositive && !record.doubleDGermlineApplied).length,
    doubleDTemplate: selected.anchored.doubleDGermlineApplied,
    doubleDRowsUsed: selected.anchored.doubleDGermlineApplied ? 1 : 0,
    selectedOrdinal: selected.row.record.ordinal,
    selectedSequenceId,
    selectedDCall: selected.anchored.dCall,
    selectedD2Call: selected.anchored.d2Call,
    selectedVjIdentity: selected.equalWeightIdentity >= 0 ? selected.equalWeightIdentity : undefined,
    selectedVIdentity: selected.v.identity ?? undefined,
    selectedJIdentity: selected.j.identity ?? undefined,
    selectedComparedColumns: selected.compared,
  };
}

/**
 * Infer one lineage germline without selecting a closest read. Known V/D/J
 * bases are voted from each member's AIRR germline alignment. A known base is
 * retained only with >=80% agreement among informative members. At columns
 * that all germline alignments mark N, an unweighted >=60% query consensus is
 * exposed separately as `uca`; the tree template deliberately remains N so
 * parsimony can treat it as unknown. Endpoint trimming requires coverage from
 * at least 20% of the loaded unique representatives.
 */
export function inferConsensusLineageGermline(rows: AirrDetailRow[]): InferredLineageGermline {
  const anchored = referenceAnchoredLineageRows(rows);
  const doubleDResolved = anchored.filter((record) => record.doubleDGermlineApplied);
  // Mixing the unchanged baseline single-D composite into a VDDJ junction vote
  // can erase D2. Once a supported D2 architecture exists, use those members
  // for the complete root vote; all ordinary lineages retain the prior path.
  const votingRows = doubleDResolved.length ? doubleDResolved : anchored;
  const minimumEndpointSupport = Math.max(1, Math.ceil(votingRows.length * 0.2));
  const template: string[] = [];
  const uca: string[] = [];
  const supported: boolean[] = [];
  let knownColumns = 0;
  let inferredColumns = 0;
  let conflictingColumns = 0;

  for (let column = 0; column < anchored[0].sequence.length; column += 1) {
    const germlineCounts = [0, 0, 0, 0];
    const queryCounts = [0, 0, 0, 0];
    let coverage = 0;
    for (const record of votingRows) {
      const germlineBase = record.germline[column] ?? "-";
      const queryBase = record.sequence[column] ?? "-";
      const germlineIndex = BASE_INDEX[germlineBase];
      const queryIndex = BASE_INDEX[queryBase];
      if (germlineIndex !== undefined) germlineCounts[germlineIndex] += 1;
      if (queryIndex !== undefined) queryCounts[queryIndex] += 1;
      if (germlineBase !== "-" || queryBase !== "-") coverage += 1;
    }
    supported.push(coverage >= minimumEndpointSupport);
    const germlineVotes = germlineCounts.reduce((sum, value) => sum + value, 0);
    const queryVotes = queryCounts.reduce((sum, value) => sum + value, 0);
    const known = plurality(germlineCounts);
    const observed = plurality(queryCounts);
    if (germlineVotes && !known.tied && known.count / germlineVotes >= 0.8) {
      template.push(known.base);
      uca.push(known.base);
      knownColumns += 1;
    } else {
      if (germlineVotes) conflictingColumns += 1;
      template.push(coverage ? "N" : "-");
      if (queryVotes && !observed.tied && observed.count / queryVotes >= 0.6) {
        uca.push(observed.base);
        inferredColumns += 1;
      } else uca.push(coverage ? "N" : "-");
    }
  }
  let start = supported.findIndex(Boolean);
  if (start < 0) start = 0;
  let end = supported.length;
  while (end > start && !supported[end - 1]) end -= 1;
  return {
    method: "consensus",
    template: template.join(""),
    uca: uca.join(""),
    trimmedTemplate: template.slice(start, end).join(""),
    trimmedUca: uca.slice(start, end).join(""),
    rowsUsed: votingRows.length,
    columns: template.length,
    startColumn: start,
    endColumn: end,
    knownColumns,
    inferredColumns,
    conflictingColumns,
    minimumEndpointSupport,
    doubleDPositiveRows: anchored.filter((record) => record.doubleDPositive).length,
    doubleDResolvedRows: doubleDResolved.length,
    doubleDIncompleteRows: anchored.filter((record) => record.doubleDPositive && !record.doubleDGermlineApplied).length,
    doubleDTemplate: doubleDResolved.length > 0,
    doubleDRowsUsed: doubleDResolved.length,
  };
}

/**
 * Reconstruct a lineage root using the least-mutated member by default.
 * Consensus voting remains available as an explicit alternative.
 */
export function inferLineageGermline(
  rows: AirrDetailRow[],
  method: LineageGermlineMethod = "closest",
): InferredLineageGermline {
  return method === "consensus" ? inferConsensusLineageGermline(rows) : inferClosestLineageGermline(rows);
}

export function lineageInputFasta(
  rows: AirrDetailRow[],
  method: LineageGermlineMethod = "closest",
): { fasta: string; frames: number[]; germline: string; inferred: InferredLineageGermline } {
  const anchored = referenceAnchoredLineageRows(rows);
  const inferred = inferLineageGermline(rows, method);
  const germline = inferred.template.padEnd(anchored[0].sequence.length, "-");
  const records = anchored.map((record) => ({
    name: record.name,
    sequence: record.sequence,
    frame: record.frame,
  }));
  records.push({ name: GERMLINE_OUTGROUP, sequence: germline, frame: records[0]?.frame ?? 0 });
  return {
    fasta: records.map((record) => `>${record.name}\n${record.sequence}`).join("\n") + "\n",
    frames: records.map((record) => record.frame),
    germline,
    inferred,
  };
}

export function quickAirrAlignment(rows: AirrDetailRow[], method: LineageGermlineMethod = "closest"): string {
  const input = lineageInputFasta(rows, method);
  const records = parseFasta(input.fasta, true);
  const maximum = Math.max(...records.map((record) => record.sequence.length));
  return records.map((record) => `>${record.name}\n${record.sequence.padEnd(maximum, "-")}`).join("\n") + "\n";
}
