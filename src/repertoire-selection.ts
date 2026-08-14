import type { AirrResultStore, DoubleDEvidenceRecord } from "./result-store";

export type SelectionTriState = "any" | "yes" | "no";
export type MotifTarget = "sequence" | "cdr3_nt" | "cdr3_aa" | "junction_aa";
export type MotifSyntax = "substring" | "iupac" | "regex";

export interface RepertoireSelectionOptions {
  sequenceId: string;
  datasetId: string;
  sampleId: string;
  subjectId: string;
  cohort: string;
  timepoint: string;
  compartment: string;
  locus: string;
  vCall: string;
  vCallIncludeAmbiguous: boolean;
  d1Call: string;
  d1CallIncludeAmbiguous: boolean;
  d2Call: string;
  jCall: string;
  jCallIncludeAmbiguous: boolean;
  cCall: string;
  cCallIncludeAmbiguous: boolean;
  isotype: string;
  cdr3Nt: string;
  cdr3Aa: string;
  motif: string;
  motifTarget: MotifTarget;
  motifSyntax: MotifSyntax;
  motifMode: "any" | "all";
  productive: SelectionTriState;
  completeVdj: SelectionTriState;
  vjInFrame: SelectionTriState;
  stopCodon: SelectionTriState;
  hasD: SelectionTriState;
  hasCdr3: SelectionTriState;
  doubleD: "any" | "positive" | "negative";
  minCdr3NtLength: number;
  maxCdr3NtLength: number;
  minCdr3AaLength: number;
  maxCdr3AaLength: number;
  minVIdentity: number;
  minJIdentity: number;
  minVMutation: number;
  maxVMutation: number;
}

export const DEFAULT_REPERTOIRE_SELECTION: RepertoireSelectionOptions = {
  sequenceId: "", datasetId: "", sampleId: "", subjectId: "", cohort: "", timepoint: "", compartment: "", locus: "", vCall: "", vCallIncludeAmbiguous: false, d1Call: "", d1CallIncludeAmbiguous: false, d2Call: "", jCall: "", jCallIncludeAmbiguous: false, cCall: "", cCallIncludeAmbiguous: false, isotype: "",
  cdr3Nt: "", cdr3Aa: "", motif: "", motifTarget: "cdr3_aa", motifSyntax: "substring", motifMode: "any",
  productive: "any", completeVdj: "any", vjInFrame: "any", stopCodon: "any", hasD: "any", hasCdr3: "any",
  doubleD: "any", minCdr3NtLength: 0, maxCdr3NtLength: 0, minCdr3AaLength: 0, maxCdr3AaLength: 0,
  minVIdentity: 0, minJIdentity: 0, minVMutation: 0, maxVMutation: 0,
};

export interface RepertoireSelectionResult {
  mask: Uint8Array;
  inputRecords: number;
  retainedRecords: number;
  discardedRecords: number;
  summary: string;
}

const IUPAC: Record<string, string> = {
  A: "A", C: "C", G: "G", T: "[TU]", U: "[TU]", R: "[AG]", Y: "[CTU]", S: "[GC]", W: "[ATU]",
  K: "[GTU]", M: "[AC]", B: "[CGTU]", D: "[AGTU]", H: "[ACTU]", V: "[ACG]", N: "[ACGTU]", X: ".",
};

function tokens(value: string): string[] {
  return value.split(/[\n,;]+/).map((item) => item.trim()).filter(Boolean);
}

function callMatches(value: string, query: string, includeAmbiguous: boolean): boolean {
  const wanted = tokens(query).map((item) => item.toUpperCase());
  if (!wanted.length) return true;
  const assignments = value.split(",").map((item) => item.trim().toUpperCase()).filter(Boolean);
  if (!includeAmbiguous && assignments.length > 1) return false;
  return wanted.some((item) => assignments.some((assignment) => assignment === item || assignment.replace(/\*.*$/, "") === item || assignment.includes(item)));
}

function facetMatches(value: string, query: string): boolean {
  const wanted = tokens(query).map((item) => item.toUpperCase());
  if (!wanted.length) return true;
  const observed = value.toUpperCase();
  return wanted.some((item) => observed === item || observed.includes(item));
}

function boolValue(value: string): boolean {
  return /^(T|TRUE|YES|1)$/i.test(value.trim());
}

function triMatches(value: string, wanted: SelectionTriState): boolean {
  return wanted === "any" || boolValue(value) === (wanted === "yes");
}

function motifRegexes(options: RepertoireSelectionOptions): RegExp[] {
  return tokens(options.motif).map((motif) => {
    if (options.motifSyntax === "regex") return new RegExp(motif, "i");
    if (options.motifSyntax === "iupac") {
      const source = [...motif.toUpperCase()].map((symbol) => IUPAC[symbol] ?? symbol.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("");
      return new RegExp(source, "i");
    }
    return new RegExp(motif.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), "i");
  });
}

export function validateRepertoireSelection(options: RepertoireSelectionOptions): string[] {
  const errors: string[] = [];
  for (const [minimum, maximum, label] of [
    [options.minCdr3NtLength, options.maxCdr3NtLength, "CDR3 nucleotide length"],
    [options.minCdr3AaLength, options.maxCdr3AaLength, "CDR3 amino-acid length"],
    [options.minVMutation, options.maxVMutation, "V mutation fraction"],
  ] as const) if (maximum > 0 && minimum > maximum) errors.push(`${label}: the minimum exceeds the maximum.`);
  try { motifRegexes(options); } catch (error) { errors.push(`Motif expression: ${error instanceof Error ? error.message : String(error)}`); }
  return errors;
}

function number(value: string): number {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : 0;
}

export function repertoireRowMatches(
  row: Record<string, string>,
  options: RepertoireSelectionOptions,
  doubleD?: DoubleDEvidenceRecord,
): boolean {
  const text = (name: string) => row[name] ?? "";
  if (options.sequenceId && !text("sequence_id").toLowerCase().includes(options.sequenceId.trim().toLowerCase())) return false;
  if (!facetMatches(text("swig_dataset_id"), options.datasetId) || !facetMatches(text("sample_id"), options.sampleId) ||
    !facetMatches(text("subject_id"), options.subjectId) || !facetMatches(text("swig_cohort"), options.cohort) ||
    !facetMatches(text("swig_timepoint"), options.timepoint) || !facetMatches(text("swig_compartment"), options.compartment)) return false;
  if (options.locus && !facetMatches(text("locus"), options.locus)) return false;
  if (!callMatches(text("v_call"), options.vCall, options.vCallIncludeAmbiguous) || !callMatches(doubleD?.values.d_call || text("d_call"), options.d1Call, options.d1CallIncludeAmbiguous) ||
    !facetMatches(doubleD?.values.d2_call || text("d2_call"), options.d2Call) || !callMatches(text("j_call"), options.jCall, options.jCallIncludeAmbiguous) ||
    !callMatches(text("c_call"), options.cCall, options.cCallIncludeAmbiguous) || !facetMatches(text("isotype"), options.isotype)) return false;
  const cdr3Nt = (text("cdr3") || text("junction")).toUpperCase();
  const cdr3Aa = (text("cdr3_aa") || text("junction_aa")).toUpperCase();
  if (options.cdr3Nt && !cdr3Nt.includes(options.cdr3Nt.replace(/\s/g, "").toUpperCase())) return false;
  if (options.cdr3Aa && !cdr3Aa.includes(options.cdr3Aa.replace(/\s/g, "").toUpperCase())) return false;
  if (!triMatches(text("productive"), options.productive) || !triMatches(text("complete_vdj"), options.completeVdj) ||
    !triMatches(text("vj_in_frame"), options.vjInFrame) || !triMatches(text("stop_codon"), options.stopCodon)) return false;
  const dPresent = Boolean(doubleD?.values.d_call || text("d_call"));
  if (options.hasD !== "any" && dPresent !== (options.hasD === "yes")) return false;
  if (options.hasCdr3 !== "any" && Boolean(cdr3Nt || cdr3Aa) !== (options.hasCdr3 === "yes")) return false;
  const isDoubleD = Boolean(doubleD?.values.d2_call || text("d2_call"));
  if (options.doubleD === "positive" && !isDoubleD) return false;
  if (options.doubleD === "negative" && isDoubleD) return false;
  if (options.minCdr3NtLength && cdr3Nt.length < options.minCdr3NtLength) return false;
  if (options.maxCdr3NtLength && cdr3Nt.length > options.maxCdr3NtLength) return false;
  if (options.minCdr3AaLength && cdr3Aa.length < options.minCdr3AaLength) return false;
  if (options.maxCdr3AaLength && cdr3Aa.length > options.maxCdr3AaLength) return false;
  const vIdentity = number(text("v_identity"));
  const jIdentity = number(text("j_identity"));
  if (options.minVIdentity && vIdentity < options.minVIdentity) return false;
  if (options.minJIdentity && jIdentity < options.minJIdentity) return false;
  const vMutation = vIdentity ? 1 - vIdentity : 0;
  if (options.minVMutation && vMutation < options.minVMutation) return false;
  if (options.maxVMutation && (!vIdentity || vMutation > options.maxVMutation)) return false;
  const expressions = motifRegexes(options);
  if (expressions.length) {
    const target = options.motifTarget === "sequence" ? (text("sequence_alignment") || text("sequence")) :
      options.motifTarget === "cdr3_nt" ? cdr3Nt : options.motifTarget === "cdr3_aa" ? cdr3Aa : text("junction_aa");
    const matched = expressions.map((expression) => expression.test(target));
    if (options.motifMode === "all" ? matched.some((value) => !value) : matched.every((value) => !value)) return false;
  }
  return true;
}

function countMask(mask: Uint8Array | undefined, count: number): number {
  if (!mask) return count;
  let total = 0;
  for (let ordinal = 0; ordinal < Math.min(mask.length, count); ordinal += 1) total += mask[ordinal] ? 1 : 0;
  return total;
}

export function selectionSummary(options: RepertoireSelectionOptions): string {
  const labels: string[] = [];
  if (options.datasetId) labels.push(`dataset ${options.datasetId}`);
  if (options.sampleId) labels.push(`sample ${options.sampleId}`);
  if (options.subjectId) labels.push(`donor ${options.subjectId}`);
  if (options.cohort) labels.push(`cohort ${options.cohort}`);
  if (options.timepoint) labels.push(`timepoint ${options.timepoint}`);
  if (options.compartment) labels.push(`compartment ${options.compartment}`);
  if (options.doubleD !== "any") labels.push(options.doubleD === "positive" ? "double-D positive" : "double-D negative");
  if (options.locus) labels.push(`locus ${options.locus}`);
  if (options.vCall) labels.push(`V ${options.vCall}`);
  if (options.d1Call) labels.push(`D1 ${options.d1Call}`);
  if (options.d2Call) labels.push(`D2 ${options.d2Call}`);
  if (options.jCall) labels.push(`J ${options.jCall}`);
  if (options.cdr3Nt || options.cdr3Aa) labels.push("CDR3 substring");
  if (options.motif) labels.push(`${options.motifTarget.replace("_", " ")} motif`);
  return labels.length ? labels.join(" · ") : "all records in the upstream working set";
}

export async function selectRepertoire(
  store: AirrResultStore,
  options: RepertoireSelectionOptions,
  baseMask?: Uint8Array,
  onProgress?: (processed: number, retained: number) => void,
): Promise<RepertoireSelectionResult> {
  const errors = validateRepertoireSelection(options);
  if (errors.length) throw new Error(errors.join(" "));
  const mask = new Uint8Array(store.count);
  const evidence = await store.doubleDRecords();
  const byOrdinal = new Map(evidence.map((record) => [record.ordinal, record]));
  const inputRecords = countMask(baseMask, store.count);

  // Positive double-D selection is sparse: resolve only evidence-positive rows instead of scanning a million-row table.
  if (options.doubleD === "positive") {
    const eligible = evidence.filter((record) => !baseMask || baseMask[record.ordinal]);
    for (let offset = 0; offset < eligible.length; offset += 1000) {
      const batch = eligible.slice(offset, offset + 1000);
      const details = await store.detailMany(batch.map((record) => record.ordinal));
      for (const detail of details) if (repertoireRowMatches(detail.values, options, byOrdinal.get(detail.record.ordinal))) mask[detail.record.ordinal] = 1;
      onProgress?.(Math.min(eligible.length, offset + batch.length), mask.reduce((sum, value) => sum + value, 0));
    }
  } else {
    let retained = 0;
    await store.scanAirrRows(
      ["sequence_id", "swig_dataset_id", "sample_id", "subject_id", "swig_cohort", "swig_timepoint", "swig_compartment", "sequence", "sequence_alignment", "locus", "v_call", "d_call", "j_call", "c_call", "isotype", "cdr3", "cdr3_aa", "junction", "junction_aa", "productive", "complete_vdj", "vj_in_frame", "stop_codon", "v_identity", "j_identity"],
      async (rows) => {
        for (const row of rows) {
          if (baseMask && !baseMask[row.ordinal]) continue;
          if (repertoireRowMatches(row.values, options, byOrdinal.get(row.ordinal))) { mask[row.ordinal] = 1; retained += 1; }
        }
      },
      { batchSize: 2500, includeMask: baseMask, onProgress: (processed) => onProgress?.(processed, retained) },
    );
  }
  const retainedRecords = countMask(mask, store.count);
  return { mask, inputRecords, retainedRecords, discardedRecords: inputRecords - retainedRecords, summary: selectionSummary(options) };
}
