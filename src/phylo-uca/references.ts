import { parseFasta } from "../post-analysis-core.ts";
import type {
  PhyloUcaAirrRow,
  PhyloUcaCandidateOptions,
  PhyloUcaCandidateReport,
  PhyloUcaReferenceRecord,
} from "./types.ts";

export interface ProjectedGermlineCandidate {
  name: string;
  sequence: string;
  projection: string;
  differences: number;
  compared: number;
  identity: number;
  observedHypothesis: boolean;
  plausibleForMember?: boolean;
  bestMemberDifferences?: number;
}

export interface PreparedPhyloUcaReferences {
  v: ProjectedGermlineCandidate[];
  d: PhyloUcaReferenceRecord[];
  j: ProjectedGermlineCandidate[];
  vEndColumn: number;
  jStartColumn: number;
  guide: string;
  report: PhyloUcaCandidateReport;
  warnings: string[];
}

export interface ObservedOnlyAlignment {
  fasta: string;
  retainedColumns: number[];
  guide: string;
  rows: number;
  columns: number;
  names: string[];
  sequences: string[];
}

function normalizeSequence(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replaceAll(".", "-").replace(/[^ACGTNRYKMSWBDHV?\-]/g, "N");
}

export function prepareObservedOnlyAlignment(curatedFasta: string, guideName: string): ObservedOnlyAlignment {
  const records = parseFasta(curatedFasta, true);
  if (!records.length) throw new Error("The lineage alignment is empty.");
  const columns = records[0].sequence.length;
  if (!columns || records.some((record) => record.sequence.length !== columns)) throw new Error("The curated lineage alignment is not rectangular.");
  const guideRecord = records.find((record) => record.name === guideName);
  const observed = records.filter((record) => record.name !== guideName);
  if (observed.length < 3) throw new Error("Phylogenetic UCA placement needs at least three observed sequences after the germline guide is removed.");
  const retainedColumns: number[] = [];
  for (let column = 0; column < columns; column += 1) {
    if (observed.some((record) => {
      const character = normalizeSequence(record.sequence)[column];
      return character !== "-" && character !== "?";
    })) retainedColumns.push(column);
  }
  if (!retainedColumns.length) throw new Error("The observed sequences contain no informative alignment columns.");
  const project = (sequence: string) => {
    const normalized = normalizeSequence(sequence);
    return retainedColumns.map((column) => normalized[column] ?? "-").join("");
  };
  const fasta = observed.map((record) => `>${record.name}\n${project(record.sequence)}`).join("\n") + "\n";
  return {
    fasta,
    retainedColumns,
    guide: project(guideRecord?.sequence ?? "-".repeat(columns)),
    rows: observed.length,
    columns: retainedColumns.length,
    names: observed.map((record) => record.name),
    sequences: observed.map((record) => project(record.sequence)),
  };
}

export function parseReferenceFasta(fasta: string, locus = ""): PhyloUcaReferenceRecord[] {
  if (!fasta.trim()) return [];
  const records = parseFasta(fasta, true).map((record) => ({
    name: record.name.split(/\s+/, 1)[0],
    sequence: normalizeSequence(record.sequence).replace(/[^ACGTN]/g, ""),
  })).filter((record) => record.name && record.sequence);
  if (!locus) return records;
  const normalizedLocus = locus.toUpperCase();
  const matching = records.filter((record) => record.name.toUpperCase().startsWith(normalizedLocus));
  return matching.length ? matching : records;
}

interface AlignmentMap {
  score: number;
  queryToTarget: Int32Array;
  matches: number;
  compared: number;
}

/** Align the complete query to the best substring of target, with free target ends. */
function semiGlobalMap(queryRaw: string, targetRaw: string): AlignmentMap {
  const query = queryRaw.toUpperCase();
  const target = targetRaw.toUpperCase();
  const rows = query.length + 1;
  const columns = target.length + 1;
  const scores = new Int32Array(rows * columns);
  const trace = new Uint8Array(rows * columns);
  const gap = -4;
  for (let row = 1; row < rows; row += 1) {
    scores[row * columns] = row * gap;
    trace[row * columns] = 1;
  }
  const substitution = (left: string, right: string) => left === "N" || right === "N" ? 0 : left === right ? 3 : -2;
  for (let row = 1; row < rows; row += 1) {
    for (let column = 1; column < columns; column += 1) {
      const diagonal = scores[(row - 1) * columns + column - 1] + substitution(query[row - 1], target[column - 1]);
      const up = scores[(row - 1) * columns + column] + gap;
      const left = scores[row * columns + column - 1] + gap;
      const index = row * columns + column;
      if (diagonal >= up && diagonal >= left) {
        scores[index] = diagonal;
        trace[index] = 0;
      } else if (up >= left) {
        scores[index] = up;
        trace[index] = 1;
      } else {
        scores[index] = left;
        trace[index] = 2;
      }
    }
  }
  let row = query.length;
  let column = 0;
  let score = Number.NEGATIVE_INFINITY;
  for (let candidate = 0; candidate < columns; candidate += 1) {
    const value = scores[row * columns + candidate];
    if (value > score) {
      score = value;
      column = candidate;
    }
  }
  const queryToTarget = new Int32Array(query.length);
  queryToTarget.fill(-1);
  let matches = 0;
  let compared = 0;
  while (row > 0) {
    const direction = trace[row * columns + column];
    if (direction === 0 && column > 0) {
      queryToTarget[row - 1] = column - 1;
      if (query[row - 1] !== "N" && target[column - 1] !== "N") {
        compared += 1;
        if (query[row - 1] === target[column - 1]) matches += 1;
      }
      row -= 1;
      column -= 1;
    } else if (direction === 1 || column === 0) {
      row -= 1;
    } else {
      column -= 1;
    }
  }
  return { score, queryToTarget, matches, compared };
}

function splitCalls(value: string): string[] {
  return value.split(",").map((call) => call.trim()).filter(Boolean);
}

function parseAlternativeCalls(value: string): string[] {
  return value.split(";").map((entry) => entry.split("|", 1)[0]?.trim()).filter(Boolean);
}

function observedHypotheses(rows: readonly PhyloUcaAirrRow[], segment: "v" | "j"): string[] {
  const calls = new Set<string>();
  for (const row of rows) {
    splitCalls(row.values[`${segment}_call`] ?? "").forEach((call) => calls.add(call));
    parseAlternativeCalls(row.values[`${segment}_alternatives`] ?? "").forEach((call) => calls.add(call));
  }
  return [...calls].sort();
}

function guideCoordinates(guide: string): { ungapped: string; columns: number[] } {
  const characters: string[] = [];
  const columns: number[] = [];
  for (let column = 0; column < guide.length; column += 1) {
    const character = guide[column];
    if (character === "-") continue;
    characters.push(/[ACGT]/.test(character) ? character : "N");
    columns.push(column);
  }
  return { ungapped: characters.join(""), columns };
}

function mappedSegmentColumns(segmentRaw: string, guide: string): number[] {
  const segment = normalizeSequence(segmentRaw).replace(/[^ACGTN]/g, "");
  if (!segment) return [];
  const coordinates = guideCoordinates(guide);
  if (!coordinates.ungapped) return [];
  const aligned = semiGlobalMap(segment, coordinates.ungapped);
  const columns: number[] = [];
  for (let query = 0; query < aligned.queryToTarget.length; query += 1) {
    const target = aligned.queryToTarget[query];
    if (target >= 0 && coordinates.columns[target] !== undefined) columns.push(coordinates.columns[target]);
  }
  return columns;
}

function projectObservedSegment(segmentRaw: string, guide: string): string {
  const segment = normalizeSequence(segmentRaw).replace(/[^ACGTN]/g, "");
  const projection = Array.from({ length: guide.length }, () => "N");
  if (!segment) return projection.join("");
  const coordinates = guideCoordinates(guide);
  if (!coordinates.ungapped) return projection.join("");
  const aligned = semiGlobalMap(segment, coordinates.ungapped);
  for (let query = 0; query < aligned.queryToTarget.length; query += 1) {
    const target = aligned.queryToTarget[query];
    if (target >= 0 && coordinates.columns[target] !== undefined) projection[coordinates.columns[target]] = segment[query];
  }
  return projection.join("");
}

function addMemberPlausibility(
  candidates: ProjectedGermlineCandidate[],
  rows: readonly PhyloUcaAirrRow[],
  segment: "v" | "j",
  guide: string,
  startColumn: number,
  endColumn: number,
  extraDifferences: number,
  minimumIdentity: number,
): void {
  const members = rows.map((row) => projectObservedSegment(row.values[`${segment}_sequence_alignment`] ?? "", guide)).filter((projection) => /[ACGT]/.test(projection));
  if (!members.length) return;
  for (const member of members) {
    const scores = candidates.map((candidate) => {
      let differences = 0;
      let compared = 0;
      for (let column = startColumn; column <= endColumn; column += 1) {
        const observed = member[column];
        const germline = candidate.projection[column];
        if (!/[ACGT]/.test(observed) || !/[ACGT]/.test(germline)) continue;
        compared += 1;
        if (observed !== germline) differences += 1;
      }
      return { candidate, differences, compared, identity: compared ? (compared - differences) / compared : 0 };
    });
    const best = scores.reduce((minimum, score) => score.compared ? Math.min(minimum, score.differences) : minimum, Number.POSITIVE_INFINITY);
    for (const score of scores) {
      if (!score.compared || score.identity < minimumIdentity || score.differences > best + Math.max(0, extraDifferences)) continue;
      score.candidate.plausibleForMember = true;
      score.candidate.bestMemberDifferences = Math.min(score.candidate.bestMemberDifferences ?? Number.POSITIVE_INFINITY, score.differences);
    }
  }
}

function inferSegmentBoundaries(guide: string, rows: readonly PhyloUcaAirrRow[]): { vEnd: number; jStart: number; warnings: string[] } {
  const vEnds: number[] = [];
  const jStarts: number[] = [];
  for (const row of rows) {
    const vColumns = mappedSegmentColumns(row.values.v_germline_alignment ?? "", guide);
    const jColumns = mappedSegmentColumns(row.values.j_germline_alignment ?? "", guide);
    if (vColumns.length) vEnds.push(Math.max(...vColumns));
    if (jColumns.length) jStarts.push(Math.min(...jColumns));
  }
  const median = (values: number[]): number => {
    const sorted = [...values].sort((left, right) => left - right);
    return sorted[Math.floor(sorted.length / 2)];
  };
  const warnings: string[] = [];
  let vEnd = vEnds.length ? median(vEnds) : -1;
  let jStart = jStarts.length ? median(jStarts) : -1;
  if (vEnd < 0 || jStart < 0 || vEnd >= jStart) {
    const unknown = [...guide].map((character, column) => character === "N" ? column : -1).filter((column) => column >= 0);
    if (unknown.length) {
      vEnd = Math.max(0, unknown[0] - 1);
      jStart = Math.min(guide.length - 1, unknown[unknown.length - 1] + 1);
    } else {
      vEnd = Math.max(0, Math.floor(guide.length * 0.72));
      jStart = Math.max(vEnd + 1, Math.floor(guide.length * 0.84));
    }
    warnings.push("V/J boundary anchors were incomplete; Swig used the N-junction span (or a conservative alignment-fraction fallback). Inspect the reported boundaries before interpreting a partial alignment.");
  }
  return { vEnd, jStart, warnings };
}

function projectCandidate(
  record: PhyloUcaReferenceRecord,
  guide: string,
  startColumn: number,
  endColumn: number,
  observed: ReadonlySet<string>,
): ProjectedGermlineCandidate {
  const anchorCharacters: string[] = [];
  const anchorColumns: number[] = [];
  for (let column = startColumn; column <= endColumn; column += 1) {
    const character = guide[column];
    if (!/[ACGT]/.test(character)) continue;
    anchorCharacters.push(character);
    anchorColumns.push(column);
  }
  const aligned = semiGlobalMap(anchorCharacters.join(""), record.sequence);
  const projection = Array.from({ length: guide.length }, () => "N");
  let differences = 0;
  let compared = 0;
  for (let anchor = 0; anchor < anchorColumns.length; anchor += 1) {
    const target = aligned.queryToTarget[anchor];
    const column = anchorColumns[anchor];
    if (target < 0) {
      projection[column] = "-";
      differences += 1;
      compared += 1;
      continue;
    }
    const character = record.sequence[target] ?? "N";
    projection[column] = /[ACGT]/.test(character) ? character : "N";
    if (/[ACGT]/.test(character)) {
      compared += 1;
      if (character !== guide[column]) differences += 1;
    }
  }
  return {
    ...record,
    projection: projection.join(""),
    differences,
    compared,
    identity: compared ? (compared - differences) / compared : 0,
    observedHypothesis: observed.has(record.name),
  };
}

function selectCandidates(
  candidates: ProjectedGermlineCandidate[],
  observed: ReadonlySet<string>,
  extraDifferences: number,
  minimumIdentity: number,
  cap: number,
): { selected: ProjectedGermlineCandidate[]; cutoff: number; truncated: boolean } {
  const best = candidates.reduce((minimum, candidate) => Math.min(minimum, candidate.differences), Number.POSITIVE_INFINITY);
  const cutoff = Number.isFinite(best) ? best + Math.max(0, extraDifferences) : Number.POSITIVE_INFINITY;
  const eligible = candidates.filter((candidate) => candidate.plausibleForMember || candidate.differences <= cutoff && candidate.identity >= minimumIdentity || observed.has(candidate.name));
  eligible.sort((left, right) => Number(right.observedHypothesis) - Number(left.observedHypothesis) || Number(Boolean(right.plausibleForMember)) - Number(Boolean(left.plausibleForMember)) || (left.bestMemberDifferences ?? left.differences) - (right.bestMemberDifferences ?? right.differences) || left.differences - right.differences || right.identity - left.identity || left.name.localeCompare(right.name));
  if (eligible.length <= cap) return { selected: eligible, cutoff, truncated: false };
  const forced = eligible.filter((candidate) => candidate.observedHypothesis);
  const retained = [...forced, ...eligible.filter((candidate) => !candidate.observedHypothesis).slice(0, Math.max(0, cap - forced.length))];
  return { selected: retained, cutoff, truncated: retained.length < eligible.length };
}

export function preparePhyloUcaReferences(
  guideRaw: string,
  rows: readonly PhyloUcaAirrRow[],
  references: { V: string; D: string; J: string },
  locus: string,
  options: PhyloUcaCandidateOptions,
): PreparedPhyloUcaReferences {
  const guide = normalizeSequence(guideRaw);
  const boundaries = inferSegmentBoundaries(guide, rows);
  const allV = parseReferenceFasta(references.V, locus);
  const allD = parseReferenceFasta(references.D, locus);
  const allJ = parseReferenceFasta(references.J, locus);
  if (!allV.length || !allJ.length) throw new Error(`The active ${locus || "lineage"} reference composition has no usable V or J records.`);
  const observedV = new Set(observedHypotheses(rows, "v"));
  const observedJ = new Set(observedHypotheses(rows, "j"));
  const projectedV = allV.map((record) => projectCandidate(record, guide, 0, boundaries.vEnd, observedV));
  const projectedJ = allJ.map((record) => projectCandidate(record, guide, boundaries.jStart, guide.length - 1, observedJ));
  // Re-run a deliberately broad, fixed-alignment candidate score for every
  // lineage member, then take the union. This avoids allowing the least-mutated
  // guide or the initially reported top call to prune a plausible V/J allele.
  addMemberPlausibility(projectedV, rows, "v", guide, 0, boundaries.vEnd, options.vMaximumExtraDifferences, options.vMinimumIdentity);
  addMemberPlausibility(projectedJ, rows, "j", guide, boundaries.jStart, guide.length - 1, options.jMaximumExtraDifferences, options.jMinimumIdentity);
  const selectedV = selectCandidates(projectedV, observedV, options.vMaximumExtraDifferences, options.vMinimumIdentity, options.maximumVCandidates);
  const selectedJ = selectCandidates(projectedJ, observedJ, options.jMaximumExtraDifferences, options.jMinimumIdentity, options.maximumJCandidates);
  if (!selectedV.selected.length || !selectedJ.selected.length) throw new Error("Broad germline screening retained no V or J candidate. Increase the candidate difference window or lower the minimum identity.");
  const warnings = [...boundaries.warnings];
  if (selectedV.truncated) warnings.push(`The V candidate cap retained ${selectedV.selected.length} of ${projectedV.filter((candidate) => candidate.plausibleForMember || candidate.differences <= selectedV.cutoff && candidate.identity >= options.vMinimumIdentity || candidate.observedHypothesis).length} plausible candidates. Observed hypotheses were never removed.`);
  if (selectedJ.truncated) warnings.push(`The J candidate cap retained ${selectedJ.selected.length} plausible candidates. Observed hypotheses were never removed.`);
  if (!allD.length && locus.endsWith("H")) warnings.push(`The active ${locus} reference composition contains no D records; the HMM can use a direct V–J path but cannot evaluate a templated D.`);
  const report: PhyloUcaCandidateReport = {
    locus,
    v: selectedV.selected.map((candidate) => candidate.name),
    d: allD.map((candidate) => candidate.name),
    j: selectedJ.selected.map((candidate) => candidate.name),
    totalVReferences: allV.length,
    totalDReferences: allD.length,
    totalJReferences: allJ.length,
    observedVHypotheses: [...observedV].sort(),
    observedJHypotheses: [...observedJ].sort(),
    vCutoffDifferences: selectedV.cutoff,
    jCutoffDifferences: selectedJ.cutoff,
    truncatedV: selectedV.truncated,
    truncatedJ: selectedJ.truncated,
  };
  return { v: selectedV.selected, d: allD, j: selectedJ.selected, vEndColumn: boundaries.vEnd, jStartColumn: boundaries.jStart, guide, report, warnings };
}
