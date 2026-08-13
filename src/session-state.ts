import type { CompiledReferences, ScopeKey } from "./reference-pack";
import type { RepertoireSelectionOptions } from "./repertoire-selection";
import type { MissingAlleleDashboard, MissingAlleleOptions } from "./germline-evidence";
import type { ShmDashboard, ShmMetricKey } from "./shm-analysis";
import type { CallingProfile } from "./swiftig-runtime";
import type { DatasetManifestEntry, PipelinePlan, StudyDesign } from "./study-design";
import type { SampleColorMap } from "./sample-colors";
import type { LineageGermlineMethod } from "./lineage-alignment";
import type { AlignmentFrameOffset } from "./lineage-phylogeny";

export const SWIG_SESSION_SCHEMA = 1 as const;

export interface LinkedAirrBinding {
  name: string;
  size: number;
  lastModified: number;
  records: number;
  headers: string[];
  fingerprint: string;
}

export interface SessionVector {
  type: "u8" | "u16" | "u32" | "i32" | "f32";
  length: number;
  base64: string;
}

export interface CollapseReplayState {
  mode: "exact" | "fad" | "conservative" | "indel";
  options: Record<string, string | number | boolean>;
  counts?: SessionVector;
  representatives?: SessionVector;
  representativeMask?: SessionVector;
  dashboard?: Record<string, unknown>;
}

export interface ChimeraReplayState {
  options: Record<string, string | number | boolean>;
  msa?: string;
  dashboard?: Record<string, unknown>;
  filterThreshold: number;
  probabilities?: SessionVector;
  dfr?: SessionVector;
  retainedMask?: SessionVector;
}

export interface LineageReplayState {
  options: Record<string, string | number | boolean>;
  assignments?: SessionVector;
  dashboard?: Record<string, unknown>;
}

export interface PostAnalysisSessionSnapshot {
  workingStages: Array<{ id: "dedup" | "chimera" | "selection"; label: string; input: number; retained: number; discarded: number; detail: string }>;
  activeMask?: SessionVector;
  collapse?: CollapseReplayState;
  chimera?: ChimeraReplayState;
  selection?: { options: RepertoireSelectionOptions; mask?: SessionVector; baseMask?: SessionVector };
  lineage?: LineageReplayState;
  /** Original lineage IDs currently opened together in the workbench. */
  selectedLineageIds?: number[];
  lineageGermlineMethod?: LineageGermlineMethod;
  query?: Record<string, unknown>;
  alignment?: { fasta: string; source: string; selectedLineageId?: number; frameOffset?: AlignmentFrameOffset };
  /** Shared nucleotide-column offset used for codon translation in the current lineage MSA. */
  alignmentFrameOffset?: AlignmentFrameOffset;
  /** Only manual/corrected alignments are retained; generated alignments are reproducible from AIRR rows. */
  editedAlignments?: Array<{ key: string; lineageIds: number[]; fasta: string; source: string; frameOffset?: AlignmentFrameOffset; savedAt: string }>;
  lineageMerges?: Array<{ id: string; label: string; originalLineageIds: number[]; createdAt: string }>;
  tree?: { rawNewick: string; rootedNewick: string; stableNewick: string; source: string; lineageIds?: number[]; run?: Record<string, unknown> };
  shm?: { metric: ShmMetricKey; dashboard: ShmDashboard; sampleOrder?: string[] };
  missingAlleles?: { options: MissingAlleleOptions; dashboard: MissingAlleleDashboard; selectedCandidateIds?: string[] };
}

export interface SwigSession {
  schema: typeof SWIG_SESSION_SCHEMA;
  application: "Swig";
  applicationVersion: string;
  savedAt: string;
  linkedAirr: LinkedAirrBinding;
  analysis: {
    inputName: string;
    species: string;
    scope: ScopeKey;
    workers: number;
    callingProfile?: CallingProfile;
    minimumIdentity: number;
    strand: 0 | 1 | 2;
    references: CompiledReferences;
    doubleD?: Record<string, string | number | boolean>;
    datasets?: DatasetManifestEntry[];
    studyDesign?: StudyDesign;
    pipeline?: PipelinePlan;
    sampleColors?: SampleColorMap;
  };
  doubleD: Array<{ ordinal: number; values: Record<string, string> }>;
  postAnalysis: PostAnalysisSessionSnapshot;
}

function bytesToBase64(bytes: Uint8Array): string {
  let output = "";
  for (let offset = 0; offset < bytes.length; offset += 0x8000) output += String.fromCharCode(...bytes.subarray(offset, offset + 0x8000));
  return btoa(output);
}

function base64ToBytes(value: string): Uint8Array {
  const binary = atob(value); const result = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) result[index] = binary.charCodeAt(index);
  return result;
}

export function packSessionVector(value: Uint8Array | Uint16Array | Uint32Array | Int32Array | Float32Array): SessionVector {
  const type = value instanceof Uint8Array ? "u8" : value instanceof Uint16Array ? "u16" : value instanceof Uint32Array ? "u32" : value instanceof Int32Array ? "i32" : "f32";
  return { type, length: value.length, base64: bytesToBase64(new Uint8Array(value.buffer, value.byteOffset, value.byteLength)) };
}

export function unpackSessionVector(value: SessionVector): Uint8Array | Uint16Array | Uint32Array | Int32Array | Float32Array {
  const bytes = base64ToBytes(value.base64); const copy = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
  const result = value.type === "u8" ? new Uint8Array(copy) : value.type === "u16" ? new Uint16Array(copy) : value.type === "u32" ? new Uint32Array(copy) : value.type === "i32" ? new Int32Array(copy) : new Float32Array(copy);
  if (result.length !== value.length) throw new Error("The saved session contains a malformed analysis vector.");
  return result;
}

export function validateSession(value: unknown): SwigSession {
  if (!value || typeof value !== "object") throw new Error("This file is not a Swig session.");
  const session = value as Partial<SwigSession>;
  if (session.schema !== SWIG_SESSION_SCHEMA || session.application !== "Swig" || !session.linkedAirr || !session.analysis || !session.postAnalysis) {
    throw new Error("This session uses an unsupported or incomplete Swig schema.");
  }
  return session as SwigSession;
}

export async function encodeSession(session: SwigSession): Promise<Blob> {
  const json = new TextEncoder().encode(JSON.stringify(session));
  if ("CompressionStream" in globalThis) {
    const compressed = await new Response(new Blob([json]).stream().pipeThrough(new CompressionStream("gzip"))).blob();
    return new Blob([compressed], { type: "application/gzip" });
  }
  return new Blob([json], { type: "application/json" });
}

export async function decodeSession(file: Blob): Promise<SwigSession> {
  const bytes = new Uint8Array(await file.arrayBuffer());
  let text: string;
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    if (!("DecompressionStream" in globalThis)) throw new Error("This browser cannot read gzip-compressed Swig sessions.");
    text = await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text();
  } else text = new TextDecoder().decode(bytes);
  try { return validateSession(JSON.parse(text)); } catch (error) {
    if (error instanceof SyntaxError) throw new Error("The selected file is not valid JSON or gzip-compressed JSON.");
    throw error;
  }
}

export function linkedAirrMatches(session: SwigSession, binding: LinkedAirrBinding): string[] {
  const expected=session.linkedAirr;const errors:string[]=[];
  if(binding.records!==expected.records)errors.push(`record count ${binding.records.toLocaleString()} does not match ${expected.records.toLocaleString()}`);
  if(binding.fingerprint!==expected.fingerprint)errors.push("content fingerprint differs");
  if(binding.headers.join("\t")!==expected.headers.join("\t"))errors.push("AIRR columns differ");
  return errors;
}

export function sessionBaseName(name:string):string{return (name.replace(/\.swig-session(?:\.json)?(?:\.gz)?$/i,"").replace(/(\.airr)?\.(tsv|csv|txt)(\.gz)?$/i,"")||"swig")+".swig-session.json.gz";}
