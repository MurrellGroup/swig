import {
  ChangeEvent,
  DragEvent,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import { AlignmentViewer } from "./alignment-view";
import { CommitNumberInput } from "./commit-number-input";
import { CommitTextInput } from "./commit-text-input";
import type { GermlinePreprocessReport, MetadataAllele } from "./germline-preprocess";
import { preprocessGermlinesInWorker } from "./germline-preprocess-client";
import { RepertoireDashboard } from "./repertoire-charts";
import { PostAnalysisWorkbench, type PostAnalysisSessionHandle } from "./post-analysis";
import { tableExtension, type TableExportFormat } from "./export-formats";
import {
  allelesForScope,
  allelesToFasta,
  availableScopes,
  compileReferences,
  loadReferencePack,
  lociForScope,
  makeDemoFasta,
  type ReferencePack,
  type ReferenceSpecies,
  type CompiledReferences,
  type LocusKey,
  type ScopeKey,
  type SegmentKey,
} from "./reference-pack";
import {
  composeReferenceOverrides,
  referenceCellKey,
  segmentAppliesToLocus,
} from "./reference-composition";
import {
  collectionsForDatabase,
  databasesForCell,
  databaseOptionsFor,
  DEFAULT_DATABASE_ID,
  loadCollectionSegment,
  type ReferenceDatabase,
} from "./reference-catalog";
import {
  AirrResultStore,
  EMPTY_FILTERS,
  inferIsotype,
  type AirrOutputHandle,
  type AirrOutputWritable,
  type AirrIndexRecord,
  type DoubleDEvidenceRecord,
  type DirectAirrOutput,
  type ResultFacets,
  type ResultFilters,
  type ResultPage,
} from "./result-store";
import {
  runSwiftIg,
  withAnalysisWebLock,
  type DoubleDScreenMode,
  type DoubleDScreenOptions,
  type CallingProfile,
} from "./swiftig-runtime";
import { streamSequenceBatches } from "./sequence-stream";
import { prepareReferenceMsa } from "./post-analysis-core";
import { createSampleColorMap, sampleColor, sampleIds, type SampleColorMap } from "./sample-colors";
import { decodeSession, encodeSession, linkedAirrMatches, sessionBaseName, SWIG_SESSION_SCHEMA, type PostAnalysisSessionSnapshot, type SwigSession } from "./session-state";
import {
  annotateAirrBatch,
  annotateDoubleDBatch,
  DATASET_SCOPE_LABELS,
  DEFAULT_PIPELINE_PLAN,
  stableDatasetSeed,
  studyScopeDefaults,
  type DatasetManifestEntry,
  type DatasetScope,
  type PipelinePlan,
  type StudyDesign,
} from "./study-design";

type AppPage = "home" | "analyze" | "results";
type InputFormat = "FASTA" | "FASTQ" | "AIRR TSV";
type InputSource = "upload" | "paste";
type OutputStorageMode = "auto" | "browser" | "disk";
type AirrRow = Record<string, string>;

interface InputData {
  name: string;
  source: string | File;
  count: number | null;
  format: InputFormat;
  formatCode: 1 | 2 | 3;
  size: number;
}

interface DatasetInput extends InputData, DatasetManifestEntry {}

interface ReferenceOverride {
  name: string;
  text: string;
  count: number;
  size: number;
  report: GermlinePreprocessReport;
  sourceKind: "database" | "upload";
  sourceDatabaseId?: string;
}

type ReferenceCellMap = Record<string, ReferenceOverride>;

interface ResultSession {
  id: number;
  store: AirrResultStore;
  total: number;
  seconds: number;
  inputName: string;
  datasets: DatasetManifestEntry[];
  studyDesign: StudyDesign;
  pipeline: PipelinePlan;
  species: string;
  scope: ScopeKey;
  facets: ResultFacets;
  summary: { assigned: number; productive: number; withCdr3: number };
  workers: number;
  outputBytes: number;
  streamedDirectly: boolean;
  inputTotal: number;
  subsampleSize: number | null;
  subsampleSeed: number | null;
  references: CompiledReferences;
  callingProfile: CallingProfile;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  doubleD: DoubleDScreenOptions;
  doubleDCount: number;
  sampleColors: SampleColorMap;
  postAnalysis?: PostAnalysisSessionSnapshot;
  restored?: boolean;
}

const SEGMENTS: SegmentKey[] = ["V", "D", "J", "C"];
const PAGE_SIZE = 50;
const MAX_INLINE_COUNT_BYTES = 2 * 1024 * 1024;
const DIRECT_OUTPUT_COUNT = 100_000;
const DIRECT_OUTPUT_INPUT_BYTES = 32 * 1024 * 1024;

interface SaveFileHandle extends AirrOutputHandle {
  createWritable: () => Promise<AirrOutputWritable>;
}

type SavePicker = (options: {
  suggestedName: string;
  types: Array<{ description: string; accept: Record<string, string[]> }>;
}) => Promise<SaveFileHandle>;

let downloadWorkerRegistration: Promise<ServiceWorkerRegistration> | null = null;

function browserWorkerLimit(): number {
  if (typeof navigator === "undefined") return 4;
  const cores = Math.max(1, navigator.hardwareConcurrency || 4);
  return Math.max(1, Math.min(16, cores - 1 || 1));
}

function recommendedWorkerCount(): number {
  const hardwareLimit = browserWorkerLimit();
  const memory = (navigator as Navigator & { deviceMemory?: number }).deviceMemory;
  const memoryLimit = memory === undefined ? 8 : memory <= 2 ? 2 : memory <= 4 ? 4 : memory <= 8 ? 8 : 12;
  return Math.max(1, Math.min(hardwareLimit, memoryLimit));
}

function outputName(inputName: string): string {
  return `${inputName.replace(/(\.gz)?\.[^.]+$/, "") || "swig"}.airr.tsv`;
}

function formattedOutputName(inputName:string,format:TableExportFormat):string{return `${inputName.replace(/(\.gz)?\.[^.]+$/, "")||"swig"}.airr${tableExtension(format)}`;}

function doubleDOutputName(inputName: string): string {
  return `${inputName.replace(/(\.gz)?\.[^.]+$/, "") || "swig"}.double-d-evidence.tsv`;
}

function formattedDoubleDOutputName(inputName:string,format:TableExportFormat):string{return `${inputName.replace(/(\.gz)?\.[^.]+$/, "")||"swig"}.double-d-evidence${tableExtension(format)}`;}

function likelyLargeInput(input: InputData, selectedCount?: number | null): boolean {
  if (selectedCount) return selectedCount >= DIRECT_OUTPUT_COUNT;
  if (input.count !== null) return input.count >= DIRECT_OUTPUT_COUNT;
  const gzipThreshold = 8 * 1024 * 1024;
  return input.size >= (input.name.toLowerCase().endsWith(".gz") ? gzipThreshold : DIRECT_OUTPUT_INPUT_BYTES);
}

function savePicker(): SavePicker | undefined {
  return (window as Window & { showSaveFilePicker?: SavePicker }).showSaveFilePicker;
}

function registerDownloadWorker(): Promise<ServiceWorkerRegistration> {
  if (!("serviceWorker" in navigator) || !window.isSecureContext) {
    return Promise.reject(new Error("Streaming downloads require a secure browser context."));
  }
  downloadWorkerRegistration ??= navigator.serviceWorker.register(
    `${import.meta.env.BASE_URL}download-worker.js`,
    { scope: import.meta.env.BASE_URL },
  );
  return downloadWorkerRegistration;
}

const LOCUS_LABELS: Record<string, string> = {
  BCR: "All BCR chains",
  TCR: "All TCR chains",
  IGH: "Heavy · IGH",
  IGK: "Kappa light · IGK",
  IGL: "Lambda light · IGL",
  TRA: "Alpha · TRA",
  TRB: "Beta · TRB",
  TRD: "Delta · TRD",
  TRG: "Gamma · TRG",
};

function bytes(value: number): string {
  if (value < 1024) return `${value} B`;
  if (value < 1024 * 1024) return `${(value / 1024).toFixed(1)} KB`;
  if (value < 1024 * 1024 * 1024) return `${(value / (1024 * 1024)).toFixed(1)} MB`;
  return `${(value / (1024 * 1024 * 1024)).toFixed(2)} GB`;
}

function percentage(value: number, total: number): string {
  return total ? `${Math.round(value / total * 100)}%` : "0%";
}

function friendlySpecies(value: string): string {
  return value.replaceAll("_", " · ");
}

function callingProfileLabel(value: CallingProfile): string {
  if (value === "igblast_compatible") return "IgBLAST-agreement";
  if (value === "igblast_balanced") return "IgBLAST-balanced";
  return "Truth-optimized";
}

function favoriteSpecies(species: ReferenceSpecies[]): ReferenceSpecies[] {
  const preferred = [
    "Homo sapiens",
    "Mus musculus_C57BL/6",
    "Macaca mulatta_AG07107",
    "Canis lupus familiaris_boxer",
    "Felis catus_Abyssinian",
    "Bos taurus_Hereford",
    "Danio rerio_Tuebingen",
  ];
  return [...species].sort((a, b) => {
    const aIndex = preferred.indexOf(a.name);
    const bIndex = preferred.indexOf(b.name);
    if (aIndex >= 0 || bIndex >= 0) {
      if (aIndex < 0) return 1;
      if (bIndex < 0) return -1;
      return aIndex - bIndex;
    }
    return a.name.localeCompare(b.name);
  });
}

function templateTiers(
  pack: ReferencePack,
  selected: ReferenceSpecies,
  scope: ScopeKey,
  segment: SegmentKey,
): MetadataAllele[][] {
  const baseTaxon = selected.name.split("_", 1)[0];
  const genus = baseTaxon.split(" ", 1)[0];
  const groups = [
    [selected],
    pack.species.filter((entry) => entry.name !== selected.name && entry.name.split("_", 1)[0] === baseTaxon),
    pack.species.filter((entry) => entry.name.split("_", 1)[0] !== baseTaxon && entry.name.startsWith(`${genus} `)),
    pack.species.filter((entry) => !entry.name.startsWith(`${genus} `)),
  ];
  const seen = new Set<string>();
  return groups.map((group) => group.flatMap((entry) => allelesForScope(entry, scope, segment)).filter((allele) => {
    const key = `${allele[0]}\u0000${allele[1]}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  })).filter((group) => group.length);
}

function formatFromName(name: string): { format: InputFormat; formatCode: 1 | 2 | 3 } | null {
  const plain = name.toLowerCase().replace(/\.gz$/, "");
  if (/\.(fa|fasta|fna|fas)$/.test(plain)) return { format: "FASTA", formatCode: 1 };
  if (/\.(fq|fastq)$/.test(plain)) return { format: "FASTQ", formatCode: 2 };
  if (/\.(tsv|csv)$/.test(plain)) return { format: "AIRR TSV", formatCode: 3 };
  return null;
}

function detectFormat(text: string): { format: InputFormat; formatCode: 1 | 2 | 3 } {
  const first = text.trimStart()[0];
  if (first === ">") return { format: "FASTA", formatCode: 1 };
  if (first === "@") return { format: "FASTQ", formatCode: 2 };
  const header = text.split(/\r?\n/, 1)[0]?.split(/\t|,/g) ?? [];
  if (header.includes("sequence")) return { format: "AIRR TSV", formatCode: 3 };
  throw new Error("Expected FASTA, FASTQ, or an AIRR table with a sequence column.");
}

function countTextRecords(text: string, formatCode: number): number {
  if (formatCode === 1) return (text.match(/^>/gm) ?? []).length;
  if (formatCode === 2) {
    const lines = text.split(/\r?\n/);
    let count = 0;
    let position = 0;
    while (position < lines.length) {
      while (position < lines.length && !lines[position]) position += 1;
      if (position >= lines.length) break;
      if (!lines[position].startsWith("@")) break;
      count += 1;
      position += 1;
      let sequenceLength = 0;
      while (position < lines.length && !lines[position].startsWith("+")) {
        sequenceLength += lines[position].replace(/\s/g, "").length;
        position += 1;
      }
      position += 1;
      let qualityLength = 0;
      while (position < lines.length && qualityLength < sequenceLength) {
        qualityLength += lines[position].length;
        position += 1;
      }
    }
    return count;
  }
  return Math.max(0, text.trimEnd().split(/\r?\n/).filter(Boolean).length - 1);
}

function inspectText(name: string, text: string, size = text.length): InputData {
  const detected = detectFormat(text);
  const count = countTextRecords(text, detected.formatCode);
  if (!count) throw new Error(`No ${detected.format} sequence records were found.`);
  return { name, source: text, count, size, ...detected };
}

async function inspectFile(file: File): Promise<InputData> {
  const namedFormat = formatFromName(file.name);
  if (file.name.toLowerCase().endsWith(".gz")) {
    if (!namedFormat) throw new Error("Name gzip inputs with .fasta.gz, .fastq.gz, or .tsv.gz.");
    return { name: file.name, source: file, count: null, size: file.size, ...namedFormat };
  }
  const sample = await file.slice(0, Math.min(file.size, MAX_INLINE_COUNT_BYTES)).text();
  const detected = namedFormat ?? detectFormat(sample);
  const count = file.size <= MAX_INLINE_COUNT_BYTES ? countTextRecords(sample, detected.formatCode) : null;
  if (count === 0) throw new Error(`No ${detected.format} sequence records were found.`);
  return { name: file.name, source: file, count, size: file.size, ...detected };
}

function inputStem(name: string): string {
  return name.replace(/\.gz$/i, "").replace(/\.(fa|fasta|fna|fas|fq|fastq|tsv|csv|txt)$/i, "").replace(/[^A-Za-z0-9_.-]+/g, "_") || "sample";
}

function datasetInput(input: InputData, ordinal: number): DatasetInput {
  const stem = inputStem(input.name);
  return {
    ...input,
    datasetId: `dataset_${ordinal}`,
    inputName: input.name,
    // Keep independent libraries separated even when files from different
    // directories share the same basename. Technical replicates are merged
    // only by an explicit shared sample ID or the technical-replicate preset.
    sampleId: `${stem}_${ordinal}`,
    subjectId: `subject_${ordinal}`,
    cohort: "cohort_1",
    timepoint: "",
    records: input.count,
  };
}

function copyPipeline(value?: PipelinePlan): PipelinePlan {
  const source = value ?? DEFAULT_PIPELINE_PLAN;
  return {
    ...DEFAULT_PIPELINE_PLAN,
    ...source,
    collapse: { ...DEFAULT_PIPELINE_PLAN.collapse, ...source.collapse },
    chimera: { ...DEFAULT_PIPELINE_PLAN.chimera, ...source.chimera },
    selection: { ...DEFAULT_PIPELINE_PLAN.selection, ...source.selection },
    lineage: { ...DEFAULT_PIPELINE_PLAN.lineage, ...source.lineage },
    shm: { ...DEFAULT_PIPELINE_PLAN.shm, ...source.shm },
    missingAlleles: { ...DEFAULT_PIPELINE_PLAN.missingAlleles, ...source.missingAlleles },
  };
}

async function readUploadedText(file: File): Promise<string> {
  if (!file.name.toLowerCase().endsWith(".gz")) return file.text();
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot decompress gzip files. Decompress the file first.");
  }
  return new Response(file.stream().pipeThrough(new DecompressionStream("gzip"))).text();
}

function downloadBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function Brand() {
  return (
    <span className="brand">
      <span className="brand-glyph" aria-hidden="true"><i /><i /><i /></span>
      <span><b>SWIG</b><small>SwiftIG · WebAssembly</small></span>
    </span>
  );
}

function AppHeader({ page, hasResults, onNavigate, onLoadSession }: {
  page: AppPage;
  hasResults: boolean;
  onNavigate: (page: AppPage) => void;
  onLoadSession: () => void;
}) {
  return (
    <header className="app-header">
      <button className="brand-button" type="button" onClick={() => onNavigate("home")}><Brand /></button>
      <nav aria-label="Primary navigation">
        <button className={page === "home" ? "active" : ""} type="button" onClick={() => onNavigate("home")}>Overview</button>
        <button className={page === "analyze" ? "active" : ""} type="button" onClick={() => onNavigate("analyze")}>Analyze</button>
        {hasResults && <button className={page === "results" ? "active" : ""} type="button" onClick={() => onNavigate("results")}>Results</button>}
        <button type="button" onClick={onLoadSession}>Load session</button>
      </nav>
      <span className="local-badge"><i /> Query data · browser only</span>
    </header>
  );
}

function WorkflowStepper({ active }: { active: 1 | 2 | 3 }) {
  return (
    <ol className="workflow-stepper" aria-label="Analysis workflow">
      {["Choose data", "Run annotation", "Explore results"].map((label, index) => {
        const value = index + 1;
        return <li key={label} className={value === active ? "active" : value < active ? "complete" : ""}><span>{value < active ? "✓" : value}</span><b>{label}</b></li>;
      })}
    </ol>
  );
}

function LandingPage({ references, onStart, onDemo }: {
  references: number | null;
  onStart: () => void;
  onDemo: () => void;
}) {
  return (
    <main className="landing-page">
      <section className="landing-hero">
        <div className="hero-copy">
          <p className="eyebrow"><span>V(D)J sequence annotation</span></p>
          <h1>Local annotation of<br /><em>BCR and TCR sequences.</em></h1>
          <p className="hero-lede">Swig runs the SwiftIG annotation core as WebAssembly for IG and TR loci. Query records, uploaded germlines, and AIRR results are processed in the browser and are not transmitted by Swig.</p>
          <div className="hero-actions">
            <button className="primary-cta" type="button" onClick={onStart}>Start an analysis <span>→</span></button>
            <button className="secondary-cta" type="button" onClick={onDemo}>Explore with demo data</button>
          </div>
          <div className="hero-proof">
            <span><b>FASTA</b> / FASTQ / AIRR · gzip</span>
            <span><b>7</b> IG + TR loci</span>
            <span><b>Query</b> records not transmitted</span>
          </div>
        </div>
        <div className="hero-instrument" aria-label="Example VDJ assignment">
          <div className="instrument-top"><span>example_read_01</span><b>productive</b></div>
          <div className="sequence-ruler"><i /><i /><i /><i /><i /><i /><i /></div>
          <div className="vdj-track"><span className="v">IGHV3-23</span><span className="n">N</span><span className="d">IGHD3-10</span><span className="n">N</span><span className="j">IGHJ4</span></div>
          <div className="instrument-alignment"><code>CAGGTGCAGCTGGTGGAGTCTGGGGGA...</code><code>||||||||||||||||||||||||| ·...</code><code>CAGGTGCAGCTGGTGGAGTCTGGGGGA...</code></div>
          <div className="instrument-callouts"><span><i className="v" />V 98.7%</span><span><i className="d" />D 100%</span><span><i className="j" />J 96.2%</span></div>
        </div>
      </section>

      <section className="workflow-story">
        <div className="section-heading"><p className="eyebrow"><span>Analysis workflow</span></p><h2>Input, annotation, and results.</h2><p>Each run has an explicit configuration, measured progress, and a browser-local result index.</p></div>
        <div className="workflow-cards">
          <article><span>01</span><div className="workflow-icon upload-icon" /><h3>Provide sequences</h3><p>Upload or paste FASTA, FASTQ, or an AIRR table. Gzip inputs are decompressed incrementally.</p></article>
          <article><span>02</span><div className="workflow-icon engine-icon" /><h3>Set references</h3><p>Choose a species and IG/TR search space, then compose IMGT, published, or uploaded V/D/J/C sets by locus.</p></article>
          <article><span>03</span><div className="workflow-icon result-icon" /><h3>Review and post-analyze</h3><p>Download AIRR TSV, inspect calls, deduplicate with abundance, assign or query lineages, and align selected groups on demand.</p></article>
        </div>
      </section>

      <section className="scale-section">
        <div className="scale-copy"><p className="eyebrow"><span>Dataset-size behavior</span></p><h2>Detailed review for small runs;<br />bounded streaming for large runs.</h2><p>Runs of up to three records open at the record view. Larger inputs are decompressed, parsed, annotated, indexed, and optionally written to disk through a bounded queue.</p></div>
        <div className="scale-grid">
          <article><strong>1–3</strong><span>Record detail opens automatically</span><small>Region and alignment layers</small></article>
          <article><strong>10³</strong><span>Facets and paged records</span><small>Indexed categorical filters</small></article>
          <article><strong>10⁶</strong><span>Bounded local batches</span><small>No full-table DOM or full-input buffer</small></article>
          <article><strong>{references ?? "—"}</strong><span>IMGT species/strain sets</span><small>Compatible alternative databases by locus</small></article>
        </div>
      </section>

      <section className="landing-final"><div><span>Annotation setup</span><h2>Configure input, locus, and germline references.</h2></div><button className="primary-cta light" type="button" onClick={onStart}>Open analysis workspace <span>→</span></button></section>
    </main>
  );
}

function ReferenceCellControl({ speciesName, locus, segment, builtInCount, builtInCoverage, reference, busy, onSelect, onFile }: {
  speciesName: string;
  locus: LocusKey;
  segment: SegmentKey;
  builtInCount: number;
  builtInCoverage?: { annotated: number; total: number };
  reference?: ReferenceOverride;
  busy: boolean;
  onSelect: (locus: LocusKey, segment: SegmentKey, sourceId: string) => void;
  onFile: (locus: LocusKey, segment: SegmentKey, file: File) => void;
}) {
  const input = useRef<HTMLInputElement>(null);
  const databases = databasesForCell(speciesName, locus, segment);
  const value = reference?.sourceKind === "upload" ? "upload" : reference?.sourceDatabaseId ?? DEFAULT_DATABASE_ID;
  const count = reference?.count ?? builtInCount;
  const coverage = reference && (segment === "V" || segment === "J")
    ? { annotated: reference.report.annotated, total: reference.report.count }
    : builtInCoverage;
  return (
    <div className={`composition-cell ${reference?.sourceKind ?? "imgt"} ${busy ? "busy" : ""}`}>
      <select aria-label={`${locus} ${segment} reference source`} value={value} disabled={busy} onChange={(event) => onSelect(locus, segment, event.target.value)}>
        <option value={DEFAULT_DATABASE_ID}>IMGT/GENE-DB{builtInCount ? "" : " · no records"}</option>
        {databases.map((database) => <option value={database.id} key={database.id}>{database.name}</option>)}
        {reference?.sourceKind === "upload" && <option value="upload">Uploaded · {reference.name}</option>}
      </select>
      <div className="composition-cell-meta">
        <b>{busy ? "Validating…" : `${count.toLocaleString()} allele${count === 1 ? "" : "s"}`}</b>
        {(segment === "V" || segment === "J") && coverage && <em className={coverage.annotated === coverage.total ? "complete" : coverage.annotated ? "partial" : "missing"} title={reference?.report.warnings.join("\n")}>{coverage.annotated.toLocaleString()}/{coverage.total.toLocaleString()} {segment === "V" ? "regions" : "anchors"}</em>}
      </div>
      <input ref={input} className="visually-hidden" type="file" accept=".fa,.fasta,.fna,.fas,.txt,.gz" onChange={(event) => {
        const file = event.target.files?.[0];
        if (file) onFile(locus, segment, file);
        event.target.value = "";
      }} />
      <button className="cell-upload" type="button" disabled={busy} onClick={() => input.current?.click()}>{reference?.sourceKind === "upload" ? "Replace FASTA" : "Upload FASTA"}</button>
    </div>
  );
}

function ReferenceCompositionMatrix({
  species,
  scope,
  references,
  busyCells,
  onSelect,
  onFile,
}: {
  species: ReferenceSpecies;
  scope: ScopeKey;
  references: ReferenceCellMap;
  busyCells: Set<string>;
  onSelect: (locus: LocusKey, segment: SegmentKey, sourceId: string) => void;
  onFile: (locus: LocusKey, segment: SegmentKey, file: File) => void;
}) {
  const loci = lociForScope(species, scope);
  return (
    <div className="composition-matrix" role="table" aria-label="Germline reference composition">
      <div className="composition-row composition-head" role="row"><span role="columnheader">Locus</span>{SEGMENTS.map((segment) => <span className={`segment-${segment.toLowerCase()}`} role="columnheader" key={segment}>{segment}<small>{segment === "V" ? "variable" : segment === "D" ? "diversity" : segment === "J" ? "joining" : "constant"}</small></span>)}</div>
      {loci.map((locus) => <div className="composition-row" role="row" key={locus}>
        <div className="composition-locus" role="rowheader"><strong>{locus}</strong><small>{LOCUS_LABELS[locus]}</small></div>
        {SEGMENTS.map((segment) => {
          if (!segmentAppliesToLocus(locus, segment)) return <div className="composition-na" role="cell" key={segment}><span>Not used</span></div>;
          const alleles = species.loci[locus]?.[segment] ?? [];
          const builtInCoverage = segment === "V"
            ? { annotated: alleles.filter((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0)).length, total: alleles.length }
            : segment === "J"
              ? { annotated: alleles.filter((allele) => Boolean(allele[2] && allele[2]![0] >= 0 && allele[2]![1] >= 0)).length, total: alleles.length }
              : undefined;
          const key = referenceCellKey(locus, segment);
          return <div role="cell" key={segment}><ReferenceCellControl speciesName={species.name} locus={locus} segment={segment} builtInCount={alleles.length} builtInCoverage={builtInCoverage} reference={references[key]} busy={busyCells.has(key)} onSelect={onSelect} onFile={onFile} /></div>;
        })}
      </div>)}
    </div>
  );
}

function CompositionSummary({
  databases,
  hasUploads,
  busy,
  release,
}: {
  databases: ReferenceDatabase[];
  hasUploads: boolean;
  busy: boolean;
  release: string;
}) {
  const sources = [`IMGT/GENE-DB ${release || "reference pack"}`, ...databases.map((database) => database.name), ...(hasUploads ? ["uploaded FASTA"] : [])];
  return (
    <section className={`database-summary ${databases.length || hasUploads ? "alternative" : ""}`} aria-label="Reference composition summary">
      <div>
        <span>Reference composition</span>
        <strong>{sources.join(" + ")}</strong>
        <p>{busy ? "Downloading and validating published germline FASTA in this browser…" : "Apply a database above, then refine any locus/segment independently in the matrix. Unchanged cells retain IMGT."}</p>
      </div>
      {databases.length > 0 && <div className="database-links">{databases.flatMap((database) => [<a href={database.sourceUrl} target="_blank" rel="noreferrer" key={`${database.id}:source`}>{database.name} source ↗</a>, <a href={database.citationUrl} target="_blank" rel="noreferrer" key={`${database.id}:citation`}>Citation ↗</a>, ...(database.terms ? [<a href={database.terms.url} target="_blank" rel="noreferrer" key={`${database.id}:terms`}>{database.terms.label} ↗</a>] : [])])}</div>}
    </section>
  );
}

function AnalysisProgress({ stage, value, onCancel, lockState, hidden }: {
  stage: string;
  value: number;
  onCancel: () => void;
  lockState: "unsupported" | "waiting" | "held";
  hidden: boolean;
}) {
  const phases = [
    { label: "Load WASM", threshold: 0.08 },
    { label: "Replicate germlines", threshold: 0.14 },
    { label: "Stream + annotate", threshold: 0.96 },
    { label: "Finalize index", threshold: 1 },
  ];
  const percent = Math.min(100, Math.max(0, Math.round(value * 100)));
  return (
    <section className="progress-stage" aria-live="polite">
      <div className="progress-orbit"><span>{percent}<small>%</small></span><i style={{ "--progress": `${percent * 3.6}deg` } as React.CSSProperties} /></div>
      <div className="progress-copy"><p className="eyebrow"><span>SwiftIG is running locally</span></p><h2>{stage}</h2><p>The page will move to Results automatically when the local AIRR index is ready.</p><p className={`background-run-state ${lockState}`}>{lockState === "held" ? `${hidden ? "Background tab · " : ""}active-run Web Lock held; Chrome Energy Saver freeze is inhibited while this analysis is running.` : lockState === "waiting" ? "Waiting for another Swig analysis tab to release the active-run lock." : "This browser does not expose Web Locks; keep the tab active for long runs."}</p><div className="main-progress"><i style={{ width: `${percent}%` }} /></div><ol>{phases.map((phase, index) => {
        const previous = index ? phases[index - 1].threshold : 0;
        const state = value >= phase.threshold ? "complete" : value >= previous ? "active" : "";
        return <li className={state} key={phase.label}><span>{state === "complete" ? "✓" : index + 1}</span>{phase.label}</li>;
      })}</ol><button className="cancel-button" type="button" onClick={onCancel}>Cancel analysis</button></div>
    </section>
  );
}

interface AlternativeHit {
  call: string;
  score: number;
  identity: number;
  queryStart: number;
  queryEnd: number;
  germlineStart: number;
  germlineEnd: number;
}

function splitCalls(value: string): string[] {
  return value.split(",").map((call) => call.trim()).filter(Boolean);
}

function parseAlternatives(value: string): AlternativeHit[] {
  if (!value) return [];
  return value.split(";").flatMap((entry) => {
    const [call, score, identity, queryStart, queryEnd, germlineStart, germlineEnd] = entry.split("|");
    const parsed = { call, score: Number(score), identity: Number(identity), queryStart: Number(queryStart), queryEnd: Number(queryEnd), germlineStart: Number(germlineStart), germlineEnd: Number(germlineEnd) };
    return call && Number.isFinite(parsed.score) ? [parsed] : [];
  });
}

function ResultDetail({ row, onClose }: { row: AirrRow; onClose: () => void }) {
  const [mode, setMode] = useState<"nt" | "aa">("nt");
  const sequenceLength = row.sequence?.length ?? 0;
  const displayedSegments = row.d2_call ? ["v", "d", "d2", "j", "c"] : ["v", "d", "j", "c"];
  const mapSegments = displayedSegments.flatMap((segment) => {
    const start = Number(row[`${segment}_sequence_start`]);
    const end = Number(row[`${segment}_sequence_end`]);
    if (!row[`${segment}_call`] || !start || !end || !sequenceLength) return [];
    return [{ segment, call: row[`${segment}_call`], left: (start - 1) / sequenceLength * 100, width: Math.max(1.5, (end - start + 1) / sequenceLength * 100) }];
  });
  const isotype = row.isotype || inferIsotype(row.c_call, row.c_sequence_alignment, row.c_identity ? Number(row.c_identity) : null);
  const uncertainSegments = ["v", "d", "j", "c"].flatMap((segment) => {
    const selected = splitCalls(row[`${segment}_call`] || "");
    const alternatives = parseAlternatives(row[`${segment}_alternatives`] || "");
    return selected.length > 1 || alternatives.length ? [{ segment, selected, alternatives }] : [];
  });

  return (
    <article className="detail-panel">
      <header className="detail-header">
        <div><span className="section-kicker">Selected rearrangement</span><h2>{row.sequence_id}</h2><div className="detail-tags"><span>{row.locus || "unassigned"}</span>{row.sample_id&&<span>sample · {row.sample_id}</span>}{row.subject_id&&<span>donor · {row.subject_id}</span>}{row.swig_timepoint&&<span>{row.swig_timepoint}</span>}<span className={row.productive === "T" ? "good" : "warn"}>{row.productive === "T" ? "Productive" : "Non-productive"}</span>{row.d2_call && <span>VDDJ screen-supported</span>}{isotype && <span className="isotype-tag">{isotype}</span>}{row.rev_comp === "T" && <span>Reverse complement</span>}</div></div>
        <button className="close-detail" type="button" onClick={onClose}>Close <span>×</span></button>
      </header>

      <section className="rearrangement-overview">
        <div className="overview-heading"><h3>Rearrangement map</h3><span>{sequenceLength.toLocaleString()} nt</span></div>
        <div className="rearrangement-track">
          <i className="baseline" />
          {mapSegments.map((segment) => <span key={segment.segment} className={`map-segment ${segment.segment}`} style={{ left: `${segment.left}%`, width: `${segment.width}%` }} title={segment.call}>{segment.segment.toUpperCase()}</span>)}
        </div>
        <div className="call-grid">
          {displayedSegments.map((segment) => {
            const selected = splitCalls(row[`${segment}_call`] || "");
            const alternatives = parseAlternatives(row[`${segment}_alternatives`] || "");
            return <div key={segment}><span>{segment.toUpperCase()} call</span><strong>{selected[0] || "—"}</strong><small>{row[`${segment}_identity`] ? `${(Number(row[`${segment}_identity`]) * 100).toFixed(1)}% identity` : "not assigned"}{selected.length > 1 ? ` · ${selected.length} co-optimal` : alternatives.length ? ` · ${alternatives.length} near-tied` : ""}</small></div>;
          })}
        </div>
      </section>

      {row.d2_call && <section className="uncertainty-panel double-d-evidence">
        <div className="uncertainty-heading"><div><span className="section-kicker">Opt-in double-D evidence</span><h3>{row.d_call} → {row.d2_call}</h3></div><p>Two ordered, non-overlapping D seeds passed the configured evidence rule and could not be explained by the allowed single-D pseudo-tandem distance. The standard SwiftIG call was {row.standard_d_call || "unassigned"}.</p></div>
        <div className="post-stat-grid compact"><article><span>V–J search span</span><strong>{row.swig_double_d_vj_span} nt</strong></article><article><span>Exact seed</span><strong>{row.swig_double_d_seed_length} nt / D</strong></article><article><span>Two-D score gain</span><strong>{row.swig_double_d_score_gain}</strong></article><article><span>Single-D Δ-distance</span><strong>{row.swig_double_d_pseudo_distance || "not explainable"}</strong></article></div>
        <div className="junction-panel"><div><span className="section-kicker">D1 → D2 insertion</span><code>{row.np2 || "zero length"}</code></div><dl><div><dt>NP2</dt><dd>{row.np2_length || "0"} nt</dd></div><div><dt>D2 → J / NP3</dt><dd>{row.np3_length || "0"} nt</dd></div><div><dt>Screen mode</dt><dd>{row.swig_double_d_mode?.replace("_", " ")}</dd></div></dl></div>
        {row.swig_double_d_alternatives && <p className="scientific-note"><span>i</span>Near-tied ordered pairs: <code>{row.swig_double_d_alternatives}</code></p>}
      </section>}

      {uncertainSegments.length > 0 && <section className="uncertainty-panel">
        <div className="uncertainty-heading"><div><span className="section-kicker">Call uncertainty</span><h3>Co-optimal and near-tied hits</h3></div><p>Comma-separated calls have the same optimal score. Additional candidates are retained within SwiftIG’s configured uncertainty window. Full alignment strings are stored only for the selected call to bound per-record output size.</p></div>
        <div className="uncertainty-stack">{uncertainSegments.map(({ segment, selected, alternatives }) => <article key={segment}>
          <header><span className={`segment-symbol segment-${segment}`}>{segment.toUpperCase()}</span><div><strong>{selected.length > 1 ? `${selected.length} co-optimal calls` : "Selected + alternate evidence"}</strong><small>Ranked alignment evidence for this segment</small></div></header>
          <div className="cooptimal-calls">{selected.map((call, index) => <span key={call} className={index === 0 ? "primary" : ""}>{call}{index === 0 ? " · reported first" : " · tied"}</span>)}</div>
          {alternatives.length > 0 && <div className="alternate-table"><div className="alternate-head"><span>Alternative</span><span>Score</span><span>Identity</span><span>Query</span><span>Germline</span></div>{alternatives.map((hit) => <div key={`${hit.call}-${hit.queryStart}`}><strong>{hit.call}</strong><span>{hit.score}</span><span>{(hit.identity * 100).toFixed(1)}%</span><span>{hit.queryStart}–{hit.queryEnd}</span><span>{hit.germlineStart}–{hit.germlineEnd}</span></div>)}</div>}
        </article>)}</div>
      </section>}

      <section className="junction-panel">
        <div><span className="section-kicker">Junction evidence</span><h3>{row.junction_aa || "No translated junction"}</h3><code>{row.junction || "—"}</code></div>
        <dl><div><dt>CDR3 AA</dt><dd>{row.cdr3_aa || "—"}</dd></div><div><dt>Length</dt><dd>{row.junction_aa_length ? `${row.junction_aa_length} aa` : "—"}</dd></div><div><dt>In frame</dt><dd>{row.vj_in_frame || "—"}</dd></div><div><dt>Stop codon</dt><dd>{row.stop_codon || "—"}</dd></div></dl>
      </section>

      <AlignmentViewer row={row} mode={mode} onMode={setMode} />

      <details className="raw-fields">
        <summary>All AIRR fields <span>{Object.values(row).filter(Boolean).length} populated</span></summary>
        <div>{Object.entries(row).filter(([, value]) => value !== "").map(([key, value]) => <dl key={key}><dt>{key}</dt><dd>{value}</dd></dl>)}</div>
      </details>
    </article>
  );
}

function DoubleDExplorer({session,onInspect}:{session:ResultSession;onInspect:(ordinal:number)=>void}){
  const [records,setRecords]=useState<DoubleDEvidenceRecord[]>([]);const [loading,setLoading]=useState(true);const [d1,setD1]=useState("");const [d2,setD2]=useState("");const [sequenceId,setSequenceId]=useState("");const [minimumGain,setMinimumGain]=useState(0);const [minimumSpan,setMinimumSpan]=useState(0);const [selected,setSelected]=useState<AirrRow|null>(null);const [selectedOrdinal,setSelectedOrdinal]=useState<number|null>(null);const [mode,setMode]=useState<"nt"|"aa">("nt");const detailRef=useRef<HTMLElement>(null);
  useEffect(()=>{let cancelled=false;setLoading(true);session.store.doubleDRecords().then(value=>{if(!cancelled){setRecords(value);setLoading(false);}});return()=>{cancelled=true;};},[session]);
  const filtered=useMemo(()=>records.filter(record=>{const values=record.values;return(!d1||String(values.d_call||"").toUpperCase().includes(d1.toUpperCase()))&&(!d2||String(values.d2_call||"").toUpperCase().includes(d2.toUpperCase()))&&(!sequenceId||String(values.sequence_id||"").toLowerCase().includes(sequenceId.toLowerCase()))&&Number(values.swig_double_d_score_gain||0)>=minimumGain&&Number(values.swig_double_d_vj_span||0)>=minimumSpan;}),[records,d1,d2,sequenceId,minimumGain,minimumSpan]);
  async function open(record:DoubleDEvidenceRecord){setSelectedOrdinal(record.ordinal);const [detail]=await session.store.detailMany([record.ordinal]);if(!detail)return;setSelected(detail.values);window.requestAnimationFrame(()=>detailRef.current?.scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth",block:"start"}));}
  const pairCounts=useMemo(()=>{const map=new Map<string,number>();for(const record of filtered){const key=`${record.values.d_call||"D1 —"} → ${record.values.d2_call||"D2 —"}`;map.set(key,(map.get(key)??0)+1);}return [...map.entries()].sort((a,b)=>b[1]-a[1]).slice(0,12);},[filtered]);
  return <section className="double-d-explorer"><header className="repertoire-heading"><div><span className="section-kicker">Opt-in VDDJ evidence</span><h2>Double-D call explorer</h2><p>Inspect ordered D1/D2 evidence and the full V–D1–D2–J alignment. These calls are separate from the standard AIRR annotation path.</p></div><span className="aggregate-badge">{filtered.length.toLocaleString()} / {records.length.toLocaleString()} supported calls</span></header>
    <div className="double-d-filter-grid"><label><span>Sequence ID</span><CommitTextInput value={sequenceId} onCommit={setSequenceId} placeholder="contains…"/></label><label><span>D1 call</span><CommitTextInput value={d1} onCommit={setD1} placeholder="IGHD…"/></label><label><span>D2 call</span><CommitTextInput value={d2} onCommit={setD2} placeholder="IGHD…"/></label><label><span>Minimum score gain</span><CommitNumberInput min="0" value={minimumGain} onCommit={setMinimumGain}/></label><label><span>Minimum V–J span</span><CommitNumberInput min="0" value={minimumSpan} onCommit={setMinimumSpan}/></label></div>
    {pairCounts.length?<div className="double-d-pair-summary">{pairCounts.map(([pair,count])=><article key={pair}><code>{pair}</code><strong>{count.toLocaleString()}</strong></article>)}</div>:null}
    <div className="lineage-table-wrap"><table><thead><tr><th>Sequence</th><th>Locus</th><th>D1</th><th>D2</th><th>Score gain</th><th>V–J span</th><th>NP2</th><th/></tr></thead><tbody>{filtered.slice(0,1000).map(record=><tr key={record.ordinal} className={selectedOrdinal===record.ordinal?"selected":""} onClick={()=>void open(record)}><td><strong>{record.values.sequence_id||`#${record.ordinal+1}`}</strong></td><td>{record.values.locus||"—"}</td><td><code>{record.values.d_call||"—"}</code></td><td><code>{record.values.d2_call||"—"}</code></td><td>{record.values.swig_double_d_score_gain||"—"}</td><td>{record.values.swig_double_d_vj_span||"—"}</td><td><code>{record.values.np2||"—"}</code></td><td><button type="button">Open alignment →</button></td></tr>)}</tbody></table>{loading?<div className="detail-loading">Loading sparse double-D evidence…</div>:!filtered.length?<div className="empty-results"><span>∅</span><h3>No supported call matches these filters.</h3></div>:null}</div>
    {selected?<section ref={detailRef} className="double-d-alignment-detail" tabIndex={-1}><header><div><span className="section-kicker">Selected VDDJ record</span><h3>{selected.sequence_id}</h3><p>{selected.v_call} → {selected.d_call} → {selected.d2_call} → {selected.j_call}</p></div><div className="result-actions"><button type="button" onClick={()=>onInspect(selectedOrdinal??0)}>Open complete AIRR record</button><button type="button" onClick={()=>setSelected(null)}>Close</button></div></header><div className="post-stat-grid compact"><article><span>Two-D score gain</span><strong>{selected.swig_double_d_score_gain||"—"}</strong></article><article><span>V–J search span</span><strong>{selected.swig_double_d_vj_span||"—"} nt</strong></article><article><span>D1→D2 insertion</span><strong>{selected.np2_length||0} nt</strong></article><article><span>D2→J insertion</span><strong>{selected.np3_length||0} nt</strong></article></div><AlignmentViewer row={selected} mode={mode} onMode={setMode}/></section>:null}
  </section>;
}

function ResultsPage({ session, onNewAnalysis }: { session: ResultSession; onNewAnalysis: () => void }) {
  const [view, setView] = useState<"repertoire" | "sequences" | "double-d" | "post">(session.pipeline.enabled ? "post" : session.total <= 3 ? "sequences" : "repertoire");
  const [postOpened, setPostOpened] = useState(Boolean(session.postAnalysis) || session.pipeline.enabled);
  const [filters, setFilters] = useState<ResultFilters>({ ...EMPTY_FILTERS });
  const [page, setPage] = useState(0);
  const [results, setResults] = useState<ResultPage>({ rows: [], hasMore: false, totalMatches: session.total, scanned: 0 });
  const [searching, setSearching] = useState(true);
  const [scanCount, setScanCount] = useState(0);
  const [selected, setSelected] = useState<AirrIndexRecord | null>(null);
  const [detail, setDetail] = useState<AirrRow | null>(null);
  const [downloading, setDownloading] = useState(false);
  const [doubleDDownloading, setDoubleDDownloading] = useState(false);
  const [downloadError, setDownloadError] = useState("");
  const [downloadFormat,setDownloadFormat]=useState<TableExportFormat>("tsv");
  const [savingSession,setSavingSession]=useState(false);
  const [datasets,setDatasets]=useState<DatasetManifestEntry[]>(()=>session.datasets.map((dataset)=>({...dataset})));
  const [metadataDraft,setMetadataDraft]=useState<DatasetManifestEntry[]>(()=>session.datasets.map((dataset)=>({...dataset})));
  const [facets,setFacets]=useState<ResultFacets>(session.facets);
  const [metadataRevision,setMetadataRevision]=useState(0);
  const [metadataBusy,setMetadataBusy]=useState(false);
  const [metadataProgress,setMetadataProgress]=useState({processed:0,total:session.total});
  const [metadataStatus,setMetadataStatus]=useState("");
  const [sampleColors,setSampleColors]=useState<SampleColorMap>(()=>createSampleColorMap(session.datasets,session.sampleColors));
  const postSessionRef=useRef<PostAnalysisSessionHandle|null>(null);
  const autoOpened = useRef(false);
  const detailRef = useRef<HTMLElement>(null);
  const scrollToDetail = useRef(false);

  useEffect(() => {
    let cancelled = false;
    const controller = new AbortController();
    setSearching(true);
    setScanCount(0);
    const delay = filters.sequenceId.trim() || filters.cdr3.trim() ? 180 : 0;
    const timer = window.setTimeout(() => {
      session.store.page(filters, page * PAGE_SIZE, PAGE_SIZE, (scanned) => {
        if (!cancelled) setScanCount(scanned);
      }, controller.signal).then((next) => {
        if (cancelled) return;
        setResults(next);
        setSearching(false);
      }).catch(() => { if (!cancelled) setSearching(false); });
    }, delay);
    return () => {
      cancelled = true;
      controller.abort();
      window.clearTimeout(timer);
    };
  }, [filters, page, session]);

  useEffect(() => {
    if (session.total <= 3 && !autoOpened.current && !selected && results.rows[0]) {
      autoOpened.current = true;
      setSelected(results.rows[0]);
    }
  }, [results.rows, selected, session.total]);

  useEffect(() => {
    let cancelled = false;
    if (!selected) {
      setDetail(null);
      return;
    }
    session.store.detail(selected).then((row) => { if (!cancelled) setDetail(row); });
    return () => { cancelled = true; };
  }, [selected, session]);

  useEffect(() => {
    if (!selected || !scrollToDetail.current) return;
    scrollToDetail.current = false;
    const frame = window.requestAnimationFrame(() => {
      detailRef.current?.scrollIntoView({
        behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth",
        block: "start",
      });
      detailRef.current?.focus({ preventScroll: true });
    });
    return () => window.cancelAnimationFrame(frame);
  }, [selected]);

  function openRecord(row: AirrIndexRecord) {
    if (selected?.ordinal === row.ordinal) {
      window.requestAnimationFrame(() => detailRef.current?.scrollIntoView({
        behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth",
        block: "start",
      }));
      return;
    }
    setDetail(null);
    scrollToDetail.current = true;
    setSelected(row);
  }

  async function inspectOrdinal(ordinal: number) {
    const [record] = await session.store.indexRecords([ordinal]);
    if (!record) return;
    setView("sequences");
    setDetail(null);
    scrollToDetail.current = true;
    setSelected(record);
  }

  function updateFilter<K extends keyof ResultFilters>(key: K, value: ResultFilters[K]) {
    setFilters((current) => ({ ...current, [key]: value }));
    setPage(0);
    setSelected(null);
  }

  function updateMetadataDraft(datasetId:string,field:"sampleId"|"subjectId"|"cohort"|"timepoint",value:string){
    setMetadataDraft((current)=>current.map((dataset)=>dataset.datasetId===datasetId?{...dataset,[field]:value}:dataset));
    setMetadataStatus("");
  }

  async function applyStudyMetadata(){
    setMetadataBusy(true);setMetadataStatus("");setDownloadError("");setMetadataProgress({processed:0,total:session.total});
    try{
      const updated=metadataDraft.map((dataset)=>({...dataset,sampleId:dataset.sampleId.trim(),subjectId:dataset.subjectId.trim(),cohort:dataset.cohort.trim(),timepoint:dataset.timepoint.trim()}));
      const nextFacets=await session.store.updateStudyMetadata(updated,(processed,total)=>setMetadataProgress({processed,total}));
      setDatasets(updated);setMetadataDraft(updated.map((dataset)=>({...dataset})));setFacets(nextFacets);
      setSampleColors((current)=>createSampleColorMap(updated,current));
      setFilters({...EMPTY_FILTERS});setPage(0);setSelected(null);setDetail(null);
      session.datasets=updated;session.facets=nextFacets;session.postAnalysis=undefined;
      postSessionRef.current=null;setMetadataRevision((value)=>value+1);
      setMetadataStatus("Study metadata applied. V(D)J calls are unchanged; downstream filters, collapse, chimera, lineage, SHM, queries, alignments, and trees were cleared and must be rerun.");
    }catch(error){setDownloadError(error instanceof Error?error.message:String(error));}
    finally{setMetadataBusy(false);}
  }

  async function downloadAll() {
    setDownloading(true);
    setDownloadError("");
    try {
      const name = formattedOutputName(session.inputName,downloadFormat);
      const picker = savePicker();
      if (picker) {
        const handle = await picker.call(window, {
          suggestedName: name,
          types: [{ description: "AIRR rearrangement table", accept: { "text/plain": [tableExtension(downloadFormat)] } }],
        });
        const writable = await handle.createWritable();
        try {
          await session.store.writeAirrFormat(downloadFormat,(part) => writable.write(part));
          await writable.close();
        } catch (error) {
          await writable.abort?.();
          throw error;
        }
      } else if(downloadFormat==="tsv") {
        if (session.streamedDirectly || session.store.hasStudyMetadataOverrides) {
          downloadBlob(await session.store.airrBlob(), name);
          return;
        }
        try {
          await registerDownloadWorker();
          const anchor = document.createElement("a");
          anchor.href = session.store.streamingDownloadUrl(import.meta.env.BASE_URL, name);
          anchor.download = name;
          anchor.click();
        } catch {
          downloadBlob(await session.store.airrBlob(), name);
        }
      } else {
        if(session.outputBytes>256*1024*1024)throw new Error("This converted export is too large for a memory-backed browser download. Use Chrome/Edge on HTTPS so Swig can stream it through Save As, or use AIRR TSV.");
        const parts:BlobPart[]=[];await session.store.writeAirrFormat(downloadFormat,async(part)=>{parts.push(part instanceof Uint8Array?part.slice().buffer:part);});downloadBlob(new Blob(parts,{type:"text/plain;charset=utf-8"}),name);
      }
    } catch (error) {
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setDownloadError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      setDownloading(false);
    }
  }

  async function downloadDoubleD() {
    setDoubleDDownloading(true);
    setDownloadError("");
    try {
      const name = formattedDoubleDOutputName(session.inputName,downloadFormat);
      const picker = savePicker();
      if (picker) {
        const handle = await picker.call(window, {
          suggestedName: name,
          types: [{ description: "Swig double-D evidence table", accept: { "text/plain": [tableExtension(downloadFormat)] } }],
        });
        const writable = await handle.createWritable();
        try {
          await session.store.writeDoubleDFormat(downloadFormat,(part) => writable.write(part));
          await writable.close();
        } catch (error) {
          await writable.abort?.();
          throw error;
        }
      } else {
        const parts: BlobPart[] = [];
        await session.store.writeDoubleDFormat(downloadFormat,async (part) => { parts.push(part as BlobPart); });
        downloadBlob(new Blob(parts, { type: "text/plain;charset=utf-8" }), name);
      }
    } catch (error) {
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setDownloadError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      setDoubleDDownloading(false);
    }
  }

  async function saveAnalysisSession(){
    setSavingSession(true);setDownloadError("");
    try{
      const postAnalysis=postSessionRef.current?await postSessionRef.current.snapshot():session.postAnalysis??{workingStages:[]};
      const artifact:SwigSession={schema:SWIG_SESSION_SCHEMA,application:"Swig",applicationVersion:"0.16.0",savedAt:new Date().toISOString(),linkedAirr:{name:outputName(session.inputName),size:session.outputBytes,lastModified:0,records:session.store.count,headers:[...session.store.airrHeaders],fingerprint:session.store.fingerprint},analysis:{inputName:session.inputName,species:session.species,scope:session.scope,workers:session.workers,callingProfile:session.callingProfile,minimumIdentity:session.minimumIdentity,strand:session.strand,references:session.references,doubleD:{...session.doubleD},datasets,studyDesign:session.studyDesign,pipeline:session.pipeline,sampleColors},doubleD:await session.store.doubleDRecords(),postAnalysis};
      downloadBlob(await encodeSession(artifact),sessionBaseName(session.inputName));
    }catch(error){setDownloadError(error instanceof Error?error.message:String(error));}finally{setSavingSession(false);}
  }

  const filtered = Object.entries(filters).some(([key, value]) => key.startsWith("min") ? Number(value) > 0 : Boolean(value));
  const pageStart = page * PAGE_SIZE + (results.rows.length ? 1 : 0);
  const pageEnd = page * PAGE_SIZE + results.rows.length;

  return (
    <main className="results-page">
      <section className="results-hero">
        <WorkflowStepper active={3} />
        <div className="results-title"><div><p className="eyebrow"><span>{session.restored?"Saved session restored":`Analysis complete · ${session.seconds.toFixed(2)} s · ${Math.round(session.total / Math.max(session.seconds, 0.001)).toLocaleString()} reads/s`}</span></p><h1>{session.total.toLocaleString()} analyzed<br /><em>rearrangements.</em></h1><p>{datasets.length.toLocaleString()} dataset{datasets.length===1?"":"s"} · {friendlySpecies(session.species)} · {LOCUS_LABELS[session.scope]} · {callingProfileLabel(session.callingProfile)} calling · {session.workers} WASM worker{session.workers === 1 ? "" : "s"} · {bytes(session.outputBytes)} AIRR{session.subsampleSize ? ` · exact random sample per dataset from ${session.inputTotal.toLocaleString()} input records (base seed ${session.subsampleSeed})` : ""}{session.doubleD.mode !== "off" ? ` · double-D screen ${session.doubleD.mode === "all" ? "all eligible junctions" : `V–J spans ≥ ${session.doubleD.minimumVjSpan} nt`}` : ""}</p></div><div className="results-actions"><label className="compact-export-format"><span>Table format</span><select value={downloadFormat} onChange={(event)=>setDownloadFormat(event.target.value as TableExportFormat)}><option value="tsv">AIRR TSV</option><option value="csv">CSV</option><option value="jsonl">JSON Lines</option></select></label><button className="download-primary" type="button" onClick={() => void downloadAll()} disabled={downloading}>{downloading ? "Writing results…" : `Download ${downloadFormat.toUpperCase()}`}<span>↓</span></button>{session.doubleDCount > 0 && <button type="button" onClick={() => void downloadDoubleD()} disabled={doubleDDownloading}>{doubleDDownloading ? "Writing double-D evidence…" : `Double-D evidence (${session.doubleDCount.toLocaleString()})`}</button>}<button type="button" onClick={()=>void saveAnalysisSession()} disabled={savingSession}>{savingSession?"Saving session…":"Save session"}</button><button type="button" onClick={onNewAnalysis}>New analysis</button>{downloadError && <small role="alert">{downloadError}</small>}</div></div>
        <div className="result-summary">
          <article><span>V + J assigned</span><strong>{session.summary.assigned.toLocaleString()}</strong><small>{percentage(session.summary.assigned, session.total)} of input</small></article>
          <article><span>Productive</span><strong>{session.summary.productive.toLocaleString()}</strong><small>{percentage(session.summary.productive, session.total)} of input</small></article>
          <article><span>CDR3 called</span><strong>{session.summary.withCdr3.toLocaleString()}</strong><small>{percentage(session.summary.withCdr3, session.total)} of input</small></article>
          <article><span>Loci observed</span><strong>{facets.loci.length}</strong><small>{facets.loci.map((item) => item.value).join(" · ") || "none"}</small></article>
          {session.doubleD.mode !== "off" && <article><span>Supported double-D</span><strong>{session.doubleDCount.toLocaleString()}</strong><small>opt-in screen · separate evidence table</small></article>}
        </div>
        <details className="study-metadata-editor"><summary><span>Study design</span><strong>Edit sample and donor metadata</strong><small>Corrections re-index downstream grouping without rerunning V(D)J assignment.</small></summary><div className="study-metadata-table"><div className="study-metadata-head"><span>Dataset</span><span>Sample</span><span>Donor / subject</span><span>Cohort</span><span>Timepoint</span></div>{metadataDraft.map((dataset)=><div className="study-metadata-row" key={dataset.datasetId}><span><strong>{dataset.inputName}</strong><small>{dataset.datasetId}</small></span><input aria-label={`${dataset.inputName} sample ID`} value={dataset.sampleId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"sampleId",event.target.value)}/><input aria-label={`${dataset.inputName} subject ID`} value={dataset.subjectId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"subjectId",event.target.value)}/><input aria-label={`${dataset.inputName} cohort`} value={dataset.cohort} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"cohort",event.target.value)}/><input aria-label={`${dataset.inputName} timepoint`} value={dataset.timepoint} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"timepoint",event.target.value)}/></div>)}</div><div className="study-metadata-actions"><button className="post-primary" type="button" disabled={metadataBusy} onClick={()=>void applyStudyMetadata()}>{metadataBusy?"Re-indexing metadata…":"Apply metadata + clear downstream state"}</button><button type="button" disabled={metadataBusy} onClick={()=>{setMetadataDraft(datasets.map((dataset)=>({...dataset})));setMetadataStatus("");}}>Discard edits</button>{metadataBusy&&<span>{metadataProgress.processed.toLocaleString()} / {metadataProgress.total.toLocaleString()}</span>}</div>{metadataStatus&&<p className="scientific-note"><span>i</span>{metadataStatus}</p>}</details>
        <details className="sample-palette-editor"><summary><span>Study palette</span><strong>Sample colors</strong><small>One association is reused in sample-level figures and phylogenies.</small></summary><div>{sampleIds(datasets).map((sample)=><label key={sample}><input type="color" value={sampleColor(sample,sampleColors)} onChange={(event)=>setSampleColors((current)=>({...current,[sample]:event.target.value}))}/><span><i style={{background:sampleColor(sample,sampleColors)}}/>{sample}</span></label>)}<button type="button" onClick={()=>setSampleColors(createSampleColorMap(datasets))}>Reset palette</button></div></details>
      </section>

      <nav className="results-view-tabs" aria-label="Results view"><button className={view === "repertoire" ? "active" : ""} type="button" onClick={() => setView("repertoire")}><span>Repertoire</span><small>Figures + composition</small></button><button className={view === "sequences" ? "active" : ""} type="button" onClick={() => setView("sequences")}><span>Sequences</span><small>Filter + inspect calls</small></button>{session.doubleDCount>0?<button className={view === "double-d" ? "active" : ""} type="button" onClick={() => setView("double-d")}><span>Double-D</span><small>{session.doubleDCount.toLocaleString()} VDDJ calls + alignments</small></button>:null}<button className={view === "post" ? "active" : ""} type="button" onClick={() => { setPostOpened(true); setView("post"); }}><span>Post-analysis</span><small>Filter + lineages + SHM + trees</small></button></nav>

      {view === "repertoire" ? <RepertoireDashboard key={metadataRevision} store={session.store} loci={facets.loci} inputName={session.inputName} samples={facets.samples} sampleColors={sampleColors} /> : view === "sequences" ? <>

      <section className="explorer-shell">
        <aside className="filter-panel">
          <div className="filter-heading"><div><span className="section-kicker">Local query</span><h2>Filter results</h2></div>{filtered && <button type="button" onClick={() => { setFilters({ ...EMPTY_FILTERS }); setPage(0); }}>Clear</button>}</div>
          {session.doubleDCount>0?<button type="button" className={`double-d-quick-filter ${filters.hasDoubleD?"active":""}`} onClick={()=>updateFilter("hasDoubleD",!filters.hasDoubleD)}><span>{filters.hasDoubleD?"✓":"DD"}</span><strong>{filters.hasDoubleD?"Showing only double-D positive":"Only double-D positive"}</strong><small>{session.doubleDCount.toLocaleString()} sparse indexed calls</small></button>:null}
          {datasets.length>1&&<div className="study-filter-grid"><label className="filter-field"><span>Dataset</span><select value={filters.datasetId} onChange={(event)=>updateFilter("datasetId",event.target.value)}><option value="">Any dataset</option>{datasets.map((item)=><option value={item.datasetId} key={item.datasetId}>{item.inputName}</option>)}</select></label><label className="filter-field"><span>Sample</span><select value={filters.sampleId} onChange={(event)=>updateFilter("sampleId",event.target.value)}><option value="">Any sample</option>{facets.samples.map((item)=><option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label><label className="filter-field"><span>Donor / subject</span><select value={filters.subjectId} onChange={(event)=>updateFilter("subjectId",event.target.value)}><option value="">Any donor</option>{facets.subjects.map((item)=><option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label><label className="filter-field"><span>Cohort</span><select value={filters.cohort} onChange={(event)=>updateFilter("cohort",event.target.value)}><option value="">Any cohort</option>{facets.cohorts.map((item)=><option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label><label className="filter-field"><span>Timepoint</span><select value={filters.timepoint} onChange={(event)=>updateFilter("timepoint",event.target.value)}><option value="">Any timepoint</option>{facets.timepoints.map((item)=><option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label></div>}
          <label className="filter-field"><span>Sequence ID contains</span><CommitTextInput type="search" value={filters.sequenceId} placeholder="e.g. clonotype_104" onCommit={(value) => updateFilter("sequenceId", value)} /></label>
          <label className="filter-field"><span>CDR3 substring <small>nt or AA</small></span><CommitTextInput className="monospace" type="search" value={filters.cdr3} placeholder="CARDR / TGTGCC…" onCommit={(value) => updateFilter("cdr3", value)} /></label>
          <div className="filter-row">
            <label className="filter-field"><span>Locus</span><select value={filters.locus} onChange={(event) => updateFilter("locus", event.target.value)}><option value="">Any locus</option>{facets.loci.map((item) => <option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>
            <label className="filter-field"><span>Productivity</span><select value={filters.productive} onChange={(event) => updateFilter("productive", event.target.value)}><option value="">Either</option><option value="T">Productive</option><option value="F">Non-productive</option></select></label>
          </div>
          {[{ key: "vCall", label: "V allele", values: facets.vCalls }, { key: "dCall", label: "D allele", values: facets.dCalls }, { key: "jCall", label: "J allele", values: facets.jCalls }, { key: "cCall", label: "C allele", values: facets.cCalls }, { key: "isotype", label: "Isotype / constant class", values: facets.isotypes }].map((field) => <label className="filter-field" key={field.key}><span>{field.label}</span><select value={filters[field.key as keyof ResultFilters] as string} onChange={(event) => updateFilter(field.key as "vCall" | "dCall" | "jCall" | "cCall" | "isotype", event.target.value)}><option value="">Any call</option>{field.values.map((item) => <option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>)}
          <details className="filter-advanced"><summary>Identity, junction + QC</summary><div>
            <label className="identity-filter"><span>Minimum V identity <b>{Math.round(filters.minVIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minVIdentity} onChange={(event) => updateFilter("minVIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum D identity <b>{Math.round(filters.minDIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minDIdentity} onChange={(event) => updateFilter("minDIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum J identity <b>{Math.round(filters.minJIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minJIdentity} onChange={(event) => updateFilter("minJIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum C identity <b>{Math.round(filters.minCIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minCIdentity} onChange={(event) => updateFilter("minCIdentity", Number(event.target.value))} /></label>
            <div className="length-filters">
              <label><span>CDR3 AA min</span><CommitNumberInput min="0" max="250" value={filters.minCdr3AaLength} blankWhenZero placeholder="Any" onCommit={(value) => updateFilter("minCdr3AaLength", value)} /></label>
              <label><span>CDR3 AA max</span><CommitNumberInput min="0" max="250" value={filters.maxCdr3AaLength} blankWhenZero placeholder="Any" onCommit={(value) => updateFilter("maxCdr3AaLength", value)} /></label>
            </div>
            <div className="qc-filter-grid">
              <label><span>VJ frame</span><select value={filters.vjInFrame} onChange={(event) => updateFilter("vjInFrame", event.target.value)}><option value="">Either</option><option value="T">In frame</option><option value="F">Out of frame</option></select></label>
              <label><span>Stop codon</span><select value={filters.stopCodon} onChange={(event) => updateFilter("stopCodon", event.target.value)}><option value="">Either</option><option value="F">Absent</option><option value="T">Present</option></select></label>
              <label><span>Completeness</span><select value={filters.completeVdj} onChange={(event) => updateFilter("completeVdj", event.target.value)}><option value="">Either</option><option value="T">Complete V(D)J</option><option value="F">Partial</option></select></label>
              <label><span>Orientation</span><select value={filters.revComp} onChange={(event) => updateFilter("revComp", event.target.value)}><option value="">Either</option><option value="F">Forward</option><option value="T">Reverse-comp.</option></select></label>
            </div>
            <label className="check-filter"><input type="checkbox" checked={filters.hasD} onChange={(event) => updateFilter("hasD", event.target.checked)} /><span>Require a D assignment</span></label>
            {session.doubleD.mode !== "off" && <label className="check-filter"><input type="checkbox" checked={filters.hasDoubleD} onChange={(event) => updateFilter("hasDoubleD", event.target.checked)} /><span>Require supported double-D evidence</span></label>}
            <label className="check-filter"><input type="checkbox" checked={filters.hasCdr3} onChange={(event) => updateFilter("hasCdr3", event.target.checked)} /><span>Require a CDR3 call</span></label>
          </div></details>
          <p className="index-note"><span>i</span> Exact gene and locus filters use browser-local indexes. Substring filters scan candidate records on demand within the browser.</p>
        </aside>

        <div className="result-browser">
          <header className="browser-heading"><div><span className="section-kicker">AIRR records</span><h2>{searching ? "Searching local index…" : results.totalMatches !== null ? `${results.totalMatches.toLocaleString()} matching records` : `${(pageEnd + (results.hasMore ? 1 : 0)).toLocaleString()}+ matching records`}</h2><p>{searching && scanCount ? `${scanCount.toLocaleString()} candidates scanned` : results.rows.length ? `Showing ${pageStart.toLocaleString()}–${pageEnd.toLocaleString()}` : "Adjust filters to broaden the query."}</p></div><span className="scale-mode">{session.total <= 3 ? "detail mode" : session.total >= 100000 ? "large-run index" : "paged index"}</span></header>
          <div className={`results-table-wrap ${searching ? "loading" : ""}`}>
            <table className="results-table">
              <thead><tr><th>Sequence</th>{datasets.length>1&&<th>Sample</th>}<th>Locus</th><th>V call</th><th>D call</th><th>J call</th><th>Isotype</th><th>CDR3 AA</th><th>Productive</th><th /></tr></thead>
              <tbody>{results.rows.map((row) => <tr className={selected?.ordinal === row.ordinal ? "selected" : ""} key={row.ordinal} tabIndex={0} onClick={() => openRecord(row)} onKeyDown={(event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); openRecord(row); } }}>
                <td><strong title={row.sequenceId}>{row.sequenceId}</strong><small>#{(row.ordinal + 1).toLocaleString()}</small></td>
                {datasets.length>1&&<td><strong className="sample-colored-label"><i style={{background:sampleColor(row.sampleId,sampleColors)}}/>{row.sampleId||"—"}</strong><small>{row.timepoint||row.subjectId||""}</small></td>}
                <td><span className="locus-pill">{row.locus || "—"}</span></td>
                <td title={row.vCall}>{row.vCall || <i>—</i>}</td><td title={[row.dCall, row.d2Call].filter(Boolean).join(" → ")}>{row.dCall || <i>—</i>}{row.d2Call && <small className="d2-table-call">→ {row.d2Call}</small>}</td><td title={row.jCall}>{row.jCall || <i>—</i>}</td><td>{row.isotype || <i>—</i>}</td>
                <td><code title={row.cdr3Aa}>{row.cdr3Aa || "—"}</code></td>
                <td><span className={`productive-dot ${row.productive === "T" ? "yes" : "no"}`} />{row.productive === "T" ? "Yes" : row.productive === "F" ? "No" : "—"}</td>
                <td><button type="button" aria-label={`Open ${row.sequenceId}`}>→</button></td>
              </tr>)}</tbody>
            </table>
            {!searching && !results.rows.length && <div className="empty-results"><span>∅</span><h3>No records match these filters.</h3><p>Try clearing an allele, identity threshold, or substring.</p></div>}
          </div>
          <footer className="table-pagination"><button type="button" disabled={!page || searching} onClick={() => setPage((value) => Math.max(0, value - 1))}>← Previous</button><span>Page {(page + 1).toLocaleString()}</span><button type="button" disabled={!results.hasMore || searching} onClick={() => setPage((value) => value + 1)}>Next →</button></footer>
        </div>
      </section>

      {selected && <section ref={detailRef} className="detail-shell" tabIndex={-1} aria-label={`Details for ${selected.sequenceId}`}>{detail ? <ResultDetail row={detail} onClose={() => setSelected(null)} /> : <div className="detail-loading">Loading selected AIRR record…</div>}</section>}
      </> : view === "double-d" ? <DoubleDExplorer key={metadataRevision} session={session} onInspect={(ordinal)=>void inspectOrdinal(ordinal)}/> : null}
      {postOpened && <div hidden={view !== "post"}><PostAnalysisWorkbench key={metadataRevision} store={session.store} references={session.references} scope={session.scope} loci={facets.loci} inputName={session.inputName} workers={session.workers} callingProfile={session.callingProfile} minimumIdentity={session.minimumIdentity} strand={session.strand} datasets={datasets} sampleColors={sampleColors} defaultCollapseScope={session.pipeline.collapse.scope} defaultLineageScope={session.pipeline.lineage.scope} autoPipeline={metadataRevision===0&&session.pipeline.enabled&&!session.postAnalysis?session.pipeline:null} onInspect={(ordinal) => void inspectOrdinal(ordinal)} sessionHandleRef={postSessionRef} initialSession={session.postAnalysis??null} /></div>}
    </main>
  );
}

export default function SwigApp() {
  const [page, setPage] = useState<AppPage>("home");
  const [pack, setPack] = useState<ReferencePack | null>(null);
  const [packError, setPackError] = useState("");
  const [speciesName, setSpeciesName] = useState("Homo sapiens");
  const [scope, setScope] = useState<ScopeKey>("BCR");
  const [inputSource, setInputSource] = useState<InputSource>("upload");
  const [fileInputs, setFileInputs] = useState<DatasetInput[]>([]);
  const [pasteText, setPasteText] = useState("");
  const [pasteMetadata, setPasteMetadata] = useState<DatasetManifestEntry>({ datasetId: "dataset_paste", inputName: "pasted-sequences.txt", sampleId: "sample_1", subjectId: "subject_1", cohort: "cohort_1", timepoint: "", records: null });
  const [studyName, setStudyName] = useState("swig-study");
  const [studyDesign, setStudyDesign] = useState<StudyDesign>("independent");
  const [pipeline, setPipeline] = useState<PipelinePlan>(() => copyPipeline());
  const [inputError, setInputError] = useState("");
  const [cellReferences, setCellReferences] = useState<ReferenceCellMap>({});
  const [busyCells, setBusyCells] = useState<Set<string>>(new Set());
  const [databaseBusy, setDatabaseBusy] = useState(false);
  const [minimumIdentity, setMinimumIdentity] = useState(0.6);
  const [callingProfile, setCallingProfile] = useState<CallingProfile>("truth_optimized");
  const [strand, setStrand] = useState<0 | 1 | 2>(0);
  const [workerCount, setWorkerCount] = useState(recommendedWorkerCount);
  const [outputStorage, setOutputStorage] = useState<OutputStorageMode>("auto");
  const [subsampleEnabled, setSubsampleEnabled] = useState(false);
  const [subsampleSize, setSubsampleSize] = useState(10_000);
  const [subsampleSeed, setSubsampleSeed] = useState(1);
  const [doubleDMode, setDoubleDMode] = useState<DoubleDScreenMode>("off");
  const [doubleDMinimumSpan, setDoubleDMinimumSpan] = useState(40);
  const [doubleDSeedLength, setDoubleDSeedLength] = useState(11);
  const [doubleDPseudoTrim, setDoubleDPseudoTrim] = useState(5);
  const [doubleDMaximumPseudoMismatches, setDoubleDMaximumPseudoMismatches] = useState(3);
  const [doubleDMinimumScoreGain, setDoubleDMinimumScoreGain] = useState(8);
  const [outputPrompt, setOutputPrompt] = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState({ stage: "Preparing analysis", value: 0 });
  const [analysisLockState, setAnalysisLockState] = useState<"unsupported" | "waiting" | "held">("unsupported");
  const [pageHidden, setPageHidden] = useState(typeof document !== "undefined" && document.hidden);
  const [runError, setRunError] = useState("");
  const [session, setSession] = useState<ResultSession | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  const sessionInputRef=useRef<HTMLInputElement>(null);
  const linkedAirrInputRef=useRef<HTMLInputElement>(null);
  const [pendingLoadedSession,setPendingLoadedSession]=useState<SwigSession|null>(null);
  const [sessionLoadError,setSessionLoadError]=useState("");
  const [sessionLoadProgress,setSessionLoadProgress]=useState({records:0,total:0,stage:""});
  const [loadingSession,setLoadingSession]=useState(false);
  const abortRef = useRef<AbortController | null>(null);
  const databaseRequestRef = useRef(0);
  const cellRequestRef = useRef<Record<string, number>>({});
  const referenceCacheRef = useRef<Map<string, ReferenceOverride>>(new Map());
  const nextDatasetIdRef = useRef(1);

  useEffect(() => {
    loadReferencePack().then(setPack).catch((error) => setPackError(error instanceof Error ? error.message : String(error)));
    if ("serviceWorker" in navigator && window.isSecureContext) void registerDownloadWorker().catch(() => undefined);
  }, []);

  useEffect(() => {
    const updateVisibility = () => setPageHidden(document.hidden);
    document.addEventListener("visibilitychange", updateVisibility);
    return () => document.removeEventListener("visibilitychange", updateVisibility);
  }, []);

  useEffect(() => {
    if (!running) return;
    const warnBeforeLeaving = (event: BeforeUnloadEvent) => {
      event.preventDefault();
      event.returnValue = "";
    };
    window.addEventListener("beforeunload", warnBeforeLeaving);
    return () => window.removeEventListener("beforeunload", warnBeforeLeaving);
  }, [running]);

  useEffect(() => {
    return () => {
      abortRef.current?.abort();
      if (session) void session.store.clear();
    };
  }, [session]);

  const speciesList = useMemo(() => favoriteSpecies(pack?.species ?? []), [pack]);
  const species = useMemo(() => speciesList.find((candidate) => candidate.name === speciesName) ?? speciesList[0], [speciesList, speciesName]);
  const scopes = useMemo(() => species ? availableScopes(species) : [], [species]);
  const activeScope = scopes.includes(scope) ? scope : scopes[0] ?? "BCR";
  const referenceOverrides = useMemo(() => species ? composeReferenceOverrides(
    lociForScope(species, activeScope),
    cellReferences,
    (locus, segment) => allelesToFasta(species.loci[locus]?.[segment] ?? []),
  ) : {}, [activeScope, cellReferences, species]);
  const compiled = useMemo(() => species ? compileReferences(species, activeScope, referenceOverrides) : null, [activeScope, referenceOverrides, species]);
  const databaseOptions = useMemo(() => databaseOptionsFor(species?.name ?? "", pack?.release ?? ""), [pack?.release, species?.name]);
  const activeReferenceEntries = useMemo(() => {
    if (!species) return [] as Array<[string, ReferenceOverride]>;
    const activeLoci = new Set(lociForScope(species, activeScope));
    return Object.entries(cellReferences).filter(([key]) => activeLoci.has(key.split(":", 1)[0] as LocusKey));
  }, [activeScope, cellReferences, species]);
  const usedDatabaseIds = useMemo(() => new Set(activeReferenceEntries.flatMap(([, reference]) => reference.sourceDatabaseId ? [reference.sourceDatabaseId] : [])), [activeReferenceEntries]);
  const usedDatabases = useMemo(() => databaseOptions.flatMap((option) => option.database && usedDatabaseIds.has(option.database.id) ? [option.database] : []), [databaseOptions, usedDatabaseIds]);
  const hasUploadedReferences = activeReferenceEntries.some(([, reference]) => reference.sourceKind === "upload");
  const compositionMode = activeReferenceEntries.length ? "mixed" : DEFAULT_DATABASE_ID;
  const databaseLabel = activeReferenceEntries.length ? ["IMGT", ...usedDatabases.map((database) => database.name), ...(hasUploadedReferences ? ["uploaded FASTA"] : [])].join(" + ") : `IMGT/GENE-DB ${pack?.release ?? "reference pack"}`;
  const receptor = activeScope.startsWith("IG") || activeScope === "BCR" ? "BCR" : "TCR";
  const receptorScopes = scopes.filter((value) => receptor === "BCR" ? value === "BCR" || value.startsWith("IG") : value === "TCR" || value.startsWith("TR"));
  const pasteInput = useMemo(() => {
    if (!pasteText.trim()) return null;
    try { return inspectText("pasted-sequences.txt", pasteText); } catch { return null; }
  }, [pasteText]);
  const activeDatasets = useMemo<DatasetInput[]>(() => {
    if (inputSource === "upload") return fileInputs;
    return pasteInput ? [{ ...pasteInput, ...pasteMetadata, inputName: pasteInput.name, records: pasteInput.count }] : [];
  }, [fileInputs, inputSource, pasteInput, pasteMetadata]);
  const activeInput = activeDatasets[0] ?? null;
  const activeInputName = activeDatasets.length > 1 ? (studyName.trim() || "swig-study") : activeInput?.name ?? "swig";
  const knownInputCount = activeDatasets.every((input) => input.count !== null)
    ? activeDatasets.reduce((sum, input) => sum + (input.count ?? 0), 0)
    : null;

  function navigate(next: AppPage) {
    if (next === "results" && !session) next = "analyze";
    setPage(next);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  async function acceptSessionFile(file:File){
    setSessionLoadError("");
    try{const loaded=await decodeSession(file);setPendingLoadedSession(loaded);setSessionLoadProgress({records:0,total:loaded.linkedAirr.records,stage:"Session manifest read; select the linked AIRR table."});}
    catch(error){setPendingLoadedSession(null);setSessionLoadError(error instanceof Error?error.message:String(error));}
  }

  async function restoreLinkedAirr(file:File){
    const saved=pendingLoadedSession;if(!saved)return;let store:AirrResultStore|undefined;
    setLoadingSession(true);setSessionLoadError("");setSessionLoadProgress({records:0,total:saved.linkedAirr.records,stage:"Streaming linked AIRR table into the local index"});
    try{
      store=new AirrResultStore();let records=0;
      for await(const batch of streamSequenceBatches({source:file,format:3,batchSize:2000,onProgress:(value)=>setSessionLoadProgress({records:value.recordsRead,total:saved.linkedAirr.records,stage:"Reading and indexing linked AIRR records"})})){
        const newline=batch.text.indexOf("\n");const header=batch.text.slice(0,newline).replace(/\r$/,"");const body=batch.text.slice(newline+1);
        if(!header.includes("\t"))throw new Error("Session loading requires the original AIRR TSV or TSV.gz file, not a comma-separated conversion.");
        await store.appendBatch(header,body);records+=batch.count;setSessionLoadProgress({records,total:saved.linkedAirr.records,stage:"Building browser-local AIRR indexes"});
      }
      await store.finalize();
      const mismatches=linkedAirrMatches(saved,{name:file.name,size:file.size,lastModified:file.lastModified,records:store.count,headers:[...store.airrHeaders],fingerprint:store.fingerprint});
      if(mismatches.length)throw new Error(`This is not the AIRR table linked by the session: ${mismatches.join("; ")}.`);
      if(saved.doubleD.length)await store.importDoubleDRecords(saved.doubleD);
      const dd={mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8,...saved.analysis.doubleD} as DoubleDScreenOptions;
      const restoredDatasets=saved.analysis.datasets?.length?saved.analysis.datasets:[{datasetId:"legacy",inputName:saved.analysis.inputName,sampleId:"sample_1",subjectId:"subject_1",cohort:"",timepoint:"",records:store.count}];
      setSessionLoadProgress({records:0,total:store.count,stage:"Applying saved study metadata to local indexes"});
      await store.updateStudyMetadata(restoredDatasets,(processed,total)=>setSessionLoadProgress({records:processed,total,stage:"Applying saved study metadata to local indexes"}));
      setSession({id:Date.now(),store,total:store.count,seconds:0,inputName:saved.analysis.inputName,datasets:restoredDatasets,studyDesign:saved.analysis.studyDesign??"independent",pipeline:copyPipeline(saved.analysis.pipeline),species:saved.analysis.species,scope:saved.analysis.scope,facets:store.facets(),summary:store.summary,workers:saved.analysis.workers,outputBytes:store.outputBytes,streamedDirectly:false,inputTotal:store.count,subsampleSize:null,subsampleSeed:null,references:saved.analysis.references,callingProfile:saved.analysis.callingProfile??"truth_optimized",minimumIdentity:saved.analysis.minimumIdentity,strand:saved.analysis.strand,doubleD:dd,doubleDCount:store.doubleDCount,sampleColors:createSampleColorMap(restoredDatasets,saved.analysis.sampleColors),postAnalysis:saved.postAnalysis,restored:true});
      setPendingLoadedSession(null);setSessionLoadProgress({records:store.count,total:store.count,stage:"Session restored"});setPage("results");window.scrollTo({top:0});
    }catch(error){if(store)await store.clear();setSessionLoadError(error instanceof Error?error.message:String(error));}finally{setLoadingSession(false);}
  }

  async function acceptInputFiles(files: File[]) {
    if (!files.length) return;
    setInputError("");
    setRunError("");
    const accepted: DatasetInput[] = [];
    const failures: string[] = [];
    for (const file of files) {
      try {
        const ordinal = nextDatasetIdRef.current++;
        accepted.push(datasetInput(await inspectFile(file), ordinal));
      } catch (error) {
        failures.push(`${file.name}: ${error instanceof Error ? error.message : String(error)}`);
      }
    }
    if (accepted.length) setFileInputs((current) => [...current, ...accepted]);
    if (failures.length) setInputError(failures.join(" "));
  }

  function updateDataset(datasetId: string, patch: Partial<DatasetManifestEntry>) {
    setFileInputs((current) => current.map((input) => input.datasetId === datasetId ? { ...input, ...patch } : input));
  }

  function applyStudyDesign(next: StudyDesign) {
    setStudyDesign(next);
    const defaults = studyScopeDefaults(next);
    setPipeline((current) => ({
      ...current,
      collapse: { ...current.collapse, scope: defaults.collapse },
      lineage: { ...current.lineage, scope: defaults.lineage },
    }));
    if (next === "technical") {
      setFileInputs((current) => current.map((input) => ({ ...input, sampleId: "sample_1", subjectId: "subject_1" })));
    } else if (next === "independent" || next === "cohort") {
      setFileInputs((current) => current.map((input, index) => ({ ...input, sampleId: `${inputStem(input.inputName)}_${index + 1}`, subjectId: `subject_${index + 1}` })));
    } else if (next === "longitudinal") {
      // A previous technical-replicate preset may have merged sample IDs. Restore
      // file-level specimens while preserving explicitly grouped donor IDs.
      setFileInputs((current) => current.map((input, index) => ({ ...input, sampleId: `${inputStem(input.inputName)}_${index + 1}` })));
    }
  }

  async function prepareReferenceOverride(
    segment: SegmentKey,
    text: string,
    name: string,
    sourceKind: "database" | "upload",
    referenceScope: ScopeKey,
    sourceDatabaseId?: string,
  ): Promise<ReferenceOverride> {
    if (!pack || !species) throw new Error("The baseline reference pack is not ready.");
    const allowedLoci = lociForScope(species, referenceScope);
    if (!allowedLoci.length) throw new Error("Choose a supported locus before adding germlines.");
    const report = await preprocessGermlinesInWorker(
      text,
      segment,
      templateTiers(pack, species, referenceScope, segment),
      allowedLoci,
    );
    return {
      name,
      text: report.fasta,
      count: report.count,
      size: new Blob([report.fasta]).size,
      report,
      sourceKind,
      sourceDatabaseId,
    };
  }

  function markCellsBusy(keys: string[], busy: boolean) {
    setBusyCells((current) => {
      const next = new Set(current);
      for (const key of keys) busy ? next.add(key) : next.delete(key);
      return next;
    });
  }

  async function databaseCellReference(
    database: ReferenceDatabase,
    locus: LocusKey,
    segment: SegmentKey,
  ): Promise<ReferenceOverride> {
    if (!species) throw new Error("Choose a species before loading a database.");
    const cacheKey = `${species.name}:${database.id}:${locus}:${segment}`;
    const cached = referenceCacheRef.current.get(cacheKey);
    if (cached) return cached;
    const collection = collectionsForDatabase(database, locus).find((candidate) => candidate.segments[segment]);
    if (!collection) throw new Error(`${database.name} does not provide ${locus} ${segment} records.`);
    const text = await loadCollectionSegment(collection, segment);
    const prepared = await prepareReferenceOverride(
      segment,
      text,
      `${database.name} · ${locus} ${segment}`,
      "database",
      locus,
      database.id,
    );
    referenceCacheRef.current.set(cacheKey, prepared);
    return prepared;
  }

  async function acceptReferenceFile(locus: LocusKey, segment: SegmentKey, file: File) {
    const contextRequest = databaseRequestRef.current;
    const key = referenceCellKey(locus, segment);
    const cellRequest = (cellRequestRef.current[key] ?? 0) + 1;
    cellRequestRef.current[key] = cellRequest;
    markCellsBusy([key], true);
    try {
      setRunError("");
      const text = await readUploadedText(file);
      if (text.trimStart()[0] !== ">") throw new Error(`${locus} ${segment} references must be FASTA.`);
      const override = await prepareReferenceOverride(segment, text, file.name, "upload", locus);
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setCellReferences((current) => ({ ...current, [key]: override }));
    } catch (error) {
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === contextRequest && cellRequestRef.current[key] === cellRequest) markCellsBusy([key], false);
    }
  }

  async function selectCellSource(locus: LocusKey, segment: SegmentKey, sourceId: string) {
    const key = referenceCellKey(locus, segment);
    const cellRequest = (cellRequestRef.current[key] ?? 0) + 1;
    cellRequestRef.current[key] = cellRequest;
    if (sourceId === "upload") return;
    setRunError("");
    if (sourceId === DEFAULT_DATABASE_ID) {
      setCellReferences((current) => {
        const next = { ...current };
        delete next[key];
        return next;
      });
      markCellsBusy([key], false);
      return;
    }
    const database = databaseOptions.find((option) => option.id === sourceId)?.database;
    if (!database) {
      setRunError("That database is not available for the selected species.");
      return;
    }
    const contextRequest = databaseRequestRef.current;
    markCellsBusy([key], true);
    try {
      const reference = await databaseCellReference(database, locus, segment);
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setCellReferences((current) => ({ ...current, [key]: reference }));
    } catch (error) {
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === contextRequest && cellRequestRef.current[key] === cellRequest) markCellsBusy([key], false);
    }
  }

  async function applyReferenceDatabase(nextId: string) {
    if (nextId === "mixed") return;
    const requestId = ++databaseRequestRef.current;
    cellRequestRef.current = {};
    setRunError("");
    if (nextId === DEFAULT_DATABASE_ID) {
      setCellReferences({});
      setBusyCells(new Set());
      setDatabaseBusy(false);
      return;
    }
    const database = databaseOptions.find((option) => option.id === nextId)?.database;
    if (!database) {
      setRunError("That database is not available for the selected species.");
      return;
    }
    let targetScope = activeScope;
    let collections = collectionsForDatabase(database, targetScope);
    if (!collections.length) {
      targetScope = database.defaultScope;
      collections = collectionsForDatabase(database, targetScope);
    }
    if (!collections.length || !scopes.includes(targetScope)) {
      setRunError(`${database.name} has no records compatible with this species and search space.`);
      return;
    }
    const currentReceptor = activeScope === "BCR" || activeScope.startsWith("IG") ? "BCR" : "TCR";
    const targetReceptor = targetScope === "BCR" || targetScope.startsWith("IG") ? "BCR" : "TCR";
    const replaceComposition = currentReceptor !== targetReceptor;
    if (targetScope !== activeScope) setScope(targetScope);
    const targets = collections.flatMap((collection) => SEGMENTS.flatMap((segment) => collection.segments[segment] ? [{ locus: collection.locus, segment }] : []));
    const keys = targets.map(({ locus, segment }) => referenceCellKey(locus, segment));
    setDatabaseBusy(true);
    markCellsBusy(keys, true);
    try {
      const preparedEntries = await Promise.all(targets.map(async ({ locus, segment }) => ({
        key: referenceCellKey(locus, segment),
        reference: await databaseCellReference(database, locus, segment),
      })));
      if (databaseRequestRef.current !== requestId) return;
      setCellReferences((current) => ({ ...(replaceComposition ? {} : current), ...Object.fromEntries(preparedEntries.map((entry) => [entry.key, entry.reference])) }));
    } catch (error) {
      if (databaseRequestRef.current !== requestId) return;
      setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === requestId) {
        setDatabaseBusy(false);
        markCellsBusy(keys, false);
      }
    }
  }

  function resetReferenceContext() {
    databaseRequestRef.current += 1;
    cellRequestRef.current = {};
    setCellReferences({});
    setBusyCells(new Set());
    setDatabaseBusy(false);
    referenceCacheRef.current.clear();
  }

  function chooseDemo() {
    if (!species) {
      navigate("analyze");
      return;
    }
    const text = makeDemoFasta(species, activeScope);
    if (!text.trim()) {
      setInputError("No complete demo locus is available for this reference selection.");
      navigate("analyze");
      return;
    }
    setPasteText(text);
    setInputSource("paste");
    setInputError("");
    navigate("analyze");
  }

  function requestRun() {
    if (!activeInput || !activeDatasets.length || !compiled || !species) return;
    const invalidDataset = activeDatasets.find((input) => !input.sampleId.trim() || !input.subjectId.trim());
    if (invalidDataset) {
      setRunError(`Dataset ${invalidDataset.inputName} needs both a sample ID and donor/subject ID.`);
      return;
    }
    if (pipeline.enabled && pipeline.lineage.enabled && (pipeline.lineage.identity < 0 || pipeline.lineage.identity > 1)) {
      setRunError("Pipeline lineage identity must be between 0 and 1.");
      return;
    }
    if (pipeline.enabled && pipeline.chimera.enabled && (pipeline.chimera.posteriorThreshold < 0 || pipeline.chimera.posteriorThreshold > 1)) {
      setRunError("Pipeline chimera posterior threshold must be between 0 and 1.");
      return;
    }
    if (pipeline.enabled && pipeline.chimera.enabled && pipeline.chimera.msaSource === "upload" && !pipeline.chimera.uploadedMsa) {
      setRunError("Choose an aligned FASTA reference MSA for the pipeline CHMMAIRRa stage, or build it from the selected references.");
      return;
    }
    if (databaseBusy || busyCells.size) {
      setRunError("Wait for reference validation to finish.");
      return;
    }
    if (!compiled.counts.V || !compiled.counts.J) {
      setRunError("This reference selection requires both V and J germline records.");
      return;
    }
    if (!compiled.annotation.V.annotated) {
      setRunError("None of the selected V records has a validated FWR/CDR delineation. Provide IMGT-gapped records, SWIGMETA/AIRR-C metadata, or a sufficiently close annotated V reference.");
      return;
    }
    if (!compiled.annotation.J.annotated) {
      setRunError("None of the selected J records has a validated coding frame and conserved F/W–G anchor. Junction and CDR3 calls cannot be made from this set.");
      return;
    }
    if (subsampleEnabled && (!Number.isFinite(subsampleSize) || subsampleSize < 1)) {
      setRunError("Random subsample size must be at least one record.");
      return;
    }
    if (doubleDMode !== "off") {
      if (!compiled.counts.D) {
        setRunError("Double-D screening requires at least one D germline in the composed reference set.");
        return;
      }
      const integerParameters = [doubleDMinimumSpan, doubleDSeedLength, doubleDPseudoTrim, doubleDMaximumPseudoMismatches, doubleDMinimumScoreGain];
      if (integerParameters.some((value) => !Number.isFinite(value) || value < 0) || doubleDSeedLength < 4) {
        setRunError("Double-D screen parameters must be non-negative integers; the exact seed must be at least 4 nt.");
        return;
      }
    }
    const selectedCount = subsampleEnabled ? Math.floor(subsampleSize) : null;
    const wantsDisk = outputStorage === "disk" || (outputStorage === "auto" && activeDatasets.some((input) => likelyLargeInput(input, selectedCount)));
    if (wantsDisk && !savePicker() && outputStorage === "disk") {
      setRunError("Direct-to-disk streaming is unavailable in this browser. Use Auto/browser storage or a Chromium-based browser.");
      return;
    }
    if (wantsDisk && savePicker()) {
      setOutputPrompt(true);
      return;
    }
    void run("browser");
  }

  async function run(outputDestination: "browser" | "disk") {
    if (!activeInput || !activeDatasets.length || !compiled || !species) return;
    setOutputPrompt(false);
    const datasetSnapshot = activeDatasets.map((input) => ({ ...input }));
    const runName = activeInputName;

    let directOutput: DirectAirrOutput | undefined;
    if (outputDestination === "disk") {
      const picker = savePicker();
      if (!picker) {
        setRunError("Direct-to-disk streaming is unavailable in this browser. Use Auto/browser storage or a Chromium-based browser.");
        return;
      }
      try {
        const handle = await picker.call(window, {
          suggestedName: outputName(runName),
          types: [{ description: "Swig AIRR analysis output", accept: { "text/tab-separated-values": [".tsv"] } }],
        });
        directOutput = { handle, writable: await handle.createWritable() };
      } catch (error) {
        if (!(error instanceof DOMException && error.name === "AbortError")) {
          setRunError(error instanceof Error ? error.message : String(error));
        }
        return;
      }
    }

    if (session) await session.store.clear();
    setSession(null);
    const store = new AirrResultStore(directOutput);
    const controller = new AbortController();
    abortRef.current = controller;
    setRunning(true);
    setRunError("");
    setProgress({ stage: "Preparing analysis", value: 0.01 });
    const started = performance.now();
    void navigator.storage?.persist?.().catch(() => false);
    const doubleD: DoubleDScreenOptions = {
      mode: doubleDMode,
      minimumVjSpan: Math.round(doubleDMinimumSpan),
      seedLength: Math.round(doubleDSeedLength),
      pseudoTrim: Math.round(doubleDPseudoTrim),
      maximumPseudoMismatches: Math.round(doubleDMaximumPseudoMismatches),
      minimumScoreGain: Math.round(doubleDMinimumScoreGain),
    };
    try {
      const result = await withAnalysisWebLock(controller.signal, setAnalysisLockState, async () => {
        let count = 0;
        let inputRecords = 0;
        let workers = 1;
        const weights = datasetSnapshot.map((input) => input.count ?? 1);
        const totalWeight = weights.reduce((sum, value) => sum + value, 0) || datasetSnapshot.length;
        let completedWeight = 0;
        for (let datasetIndex = 0; datasetIndex < datasetSnapshot.length; datasetIndex += 1) {
          const input = datasetSnapshot[datasetIndex];
          const manifest: DatasetManifestEntry = {
            datasetId: input.datasetId,
            inputName: input.inputName,
            sampleId: input.sampleId,
            subjectId: input.subjectId,
            cohort: input.cohort,
            timepoint: input.timepoint,
            records: input.count,
          };
          const weight = weights[datasetIndex];
          const completed = await runSwiftIg({
            query: input.source,
            format: input.formatCode,
            references: compiled,
            callingProfile,
            minimumIdentity,
            strand,
            workers: workerCount,
            countHint: input.count,
            subsample: subsampleEnabled ? { size: Math.floor(subsampleSize), seed: stableDatasetSeed(subsampleSeed, datasetIndex) } : undefined,
            doubleD,
            signal: controller.signal,
            onProgress: (stage, value) => setProgress({
              stage: datasetSnapshot.length > 1 ? `${input.sampleId} · ${stage}` : stage,
              value: Math.min(0.995, (completedWeight + Math.max(0, Math.min(1, value)) * weight) / totalWeight),
            }),
            onBatch: async (batch) => {
              const annotated = annotateAirrBatch(batch.header, batch.body, manifest);
              const doubleDBatch = batch.doubleDHeader !== undefined && batch.doubleDBody !== undefined
                ? annotateDoubleDBatch(batch.doubleDHeader, batch.doubleDBody, manifest)
                : undefined;
              await store.appendBatch(annotated.header, annotated.body, doubleDBatch);
              setProgress((current) => ({
                ...current,
                stage: `${input.sampleId} · committed ${batch.processed.toLocaleString()}${batch.total === null ? "" : ` / ${batch.total.toLocaleString()}`} AIRR records${store.doubleDCount ? ` · ${store.doubleDCount.toLocaleString()} double-D supported` : ""}`,
              }));
            },
          });
          count += completed.count;
          inputRecords += completed.inputRecords;
          workers = Math.max(workers, completed.workers);
          completedWeight += weight;
        }
        await store.finalize();
        return { count, inputRecords, workers };
      });
      setProgress({ stage: "Results ready", value: 1 });
      setSession({
        id: Date.now(),
        store,
        total: result.count,
        seconds: (performance.now() - started) / 1000,
        inputName: runName,
        datasets: datasetSnapshot.map((input) => ({ datasetId: input.datasetId, inputName: input.inputName, sampleId: input.sampleId, subjectId: input.subjectId, cohort: input.cohort, timepoint: input.timepoint, records: input.count })),
        studyDesign,
        pipeline: copyPipeline(pipeline),
        species: species.name,
        scope: activeScope,
        facets: store.facets(),
        summary: store.summary,
        workers: result.workers,
        outputBytes: store.outputBytes,
        streamedDirectly: store.streamedDirectly,
        inputTotal: result.inputRecords,
        subsampleSize: subsampleEnabled ? result.count : null,
        subsampleSeed: subsampleEnabled ? Math.trunc(subsampleSeed) : null,
        references: compiled,
        callingProfile,
        minimumIdentity,
        strand,
        doubleD,
        doubleDCount: store.doubleDCount,
        sampleColors: createSampleColorMap(datasetSnapshot),
      });
      setPage("results");
      window.scrollTo({ top: 0 });
    } catch (error) {
      await store.abort();
      await store.clear();
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setRunError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      abortRef.current = null;
      setRunning(false);
      setAnalysisLockState("unsupported");
    }
  }

  const hasBcr = scopes.includes("BCR");
  const hasTcr = scopes.includes("TCR");

  return (
    <div className="site-shell">
      <AppHeader page={page} hasResults={Boolean(session)} onNavigate={navigate} onLoadSession={()=>sessionInputRef.current?.click()} />
      <input ref={sessionInputRef} className="visually-hidden" type="file" accept=".swig-session,.json,.gz" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)void acceptSessionFile(file);}}/>
      <input ref={linkedAirrInputRef} className="visually-hidden" type="file" accept=".tsv,.tsv.gz,.gz" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)void restoreLinkedAirr(file);}}/>
      {page === "home" && <LandingPage references={pack?.species.length ?? null} onStart={() => navigate("analyze")} onDemo={chooseDemo} />}

      {page === "analyze" && (
        <main className="analysis-page">
          <section className="analysis-intro"><WorkflowStepper active={running ? 2 : 1} /><div><p className="eyebrow"><span>Analysis workspace</span></p><h1>{running ? "Calling rearrangements…" : "Configure an annotation run."}</h1><p>{running ? "SwiftIG is processing bounded batches and writing AIRR records into a browser-local index." : "Provide sequences, select the biological search space, and specify any germline replacements."}</p></div></section>

          {running ? <AnalysisProgress stage={progress.stage} value={progress.value} onCancel={() => abortRef.current?.abort()} lockState={analysisLockState} hidden={pageHidden} /> : (
            <div className="analysis-layout">
              <div className="analysis-forms">
                <section className="analysis-card input-card">
                  <header><span className="card-number">01</span><div><h2>Datasets and study structure</h2><p>Upload one or more datasets, or paste one dataset directly.</p></div>{activeInput && <span className="ready-tag">{activeDatasets.length} ready</span>}</header>
                  <div className="source-tabs" role="tablist"><button className={inputSource === "upload" ? "active" : ""} type="button" onClick={() => setInputSource("upload")}>Upload dataset(s)</button><button className={inputSource === "paste" ? "active" : ""} type="button" onClick={() => setInputSource("paste")}>Paste one dataset</button></div>
                  <input ref={inputRef} className="visually-hidden" type="file" multiple accept=".fa,.fasta,.fna,.fas,.fq,.fastq,.tsv,.csv,.txt,.gz" onChange={(event: ChangeEvent<HTMLInputElement>) => {
                    const files = [...(event.target.files ?? [])];
                    if (files.length) void acceptInputFiles(files);
                    event.target.value = "";
                  }} />
                  {inputSource === "upload" ? fileInputs.length ? (
                    <div className="dataset-import-stack">
                      <div className="dataset-import-heading"><div><strong>{fileInputs.length.toLocaleString()} dataset{fileInputs.length === 1 ? "" : "s"}</strong><span>{knownInputCount === null ? "Record totals will be counted during analysis" : `${knownInputCount.toLocaleString()} total records`}</span></div><button type="button" onClick={() => inputRef.current?.click()}>＋ Add datasets</button></div>
                      <div className="dataset-manifest-table" role="table" aria-label="Dataset metadata">
                        <div className="dataset-manifest-head" role="row"><span>Input</span><span>Sample ID</span><span>Donor / subject</span><span>Cohort</span><span>Timepoint</span><span /></div>
                        {fileInputs.map((input) => <div className="dataset-manifest-row" role="row" key={input.datasetId}>
                          <div><strong title={input.name}>{input.name}</strong><small>{input.datasetId} · {input.count === null ? "stream counted" : `${input.count.toLocaleString()} records`} · {input.format} · {bytes(input.size)}</small></div>
                          <label><span className="visually-hidden">Sample ID for {input.name}</span><input value={input.sampleId} onChange={(event) => updateDataset(input.datasetId, { sampleId: event.target.value })} /></label>
                          <label><span className="visually-hidden">Donor or subject for {input.name}</span><input value={input.subjectId} onChange={(event) => updateDataset(input.datasetId, { subjectId: event.target.value })} /></label>
                          <label><span className="visually-hidden">Cohort for {input.name}</span><input value={input.cohort} onChange={(event) => updateDataset(input.datasetId, { cohort: event.target.value })} /></label>
                          <label><span className="visually-hidden">Timepoint for {input.name}</span><input value={input.timepoint} placeholder="e.g. day_0" onChange={(event) => updateDataset(input.datasetId, { timepoint: event.target.value })} /></label>
                          <button type="button" aria-label={`Remove ${input.name}`} onClick={() => setFileInputs((current) => current.filter((candidate) => candidate.datasetId !== input.datasetId))}>×</button>
                        </div>)}
                      </div>
                    </div>
                  ) : (
                    <button className="input-dropzone" type="button" onClick={() => inputRef.current?.click()} onDragOver={(event: DragEvent) => event.preventDefault()} onDrop={(event: DragEvent) => {
                      event.preventDefault();
                      const files = [...event.dataTransfer.files];
                      if (files.length) void acceptInputFiles(files);
                    }}><span>＋</span><strong>Drop one or more datasets here</strong><small>.fasta(.gz) · .fastq(.gz) · AIRR .tsv(.gz)</small><i>Choose files</i></button>
                  ) : (
                    <div className="paste-input"><textarea spellCheck={false} value={pasteText} onChange={(event) => setPasteText(event.target.value)} placeholder={">sequence_1\nCAGGTGCAGCTGGTG...\n\n—or—\n\nsequence_id\tsequence\nread_1\tCAGGTGCAGCTGGTG..."} /><footer><span>{pasteInput ? `${pasteInput.count?.toLocaleString()} ${pasteInput.format} records detected` : pasteText.trim() ? "Waiting for valid FASTA, FASTQ, or AIRR…" : "Nothing pasted yet"}</span><button type="button" onClick={chooseDemo}>Insert demo</button></footer>{pasteInput && <div className="paste-metadata"><label><span>Sample ID</span><input value={pasteMetadata.sampleId} onChange={(event)=>setPasteMetadata((current)=>({...current,sampleId:event.target.value}))}/></label><label><span>Donor / subject</span><input value={pasteMetadata.subjectId} onChange={(event)=>setPasteMetadata((current)=>({...current,subjectId:event.target.value}))}/></label><label><span>Cohort</span><input value={pasteMetadata.cohort} onChange={(event)=>setPasteMetadata((current)=>({...current,cohort:event.target.value}))}/></label><label><span>Timepoint</span><input value={pasteMetadata.timepoint} onChange={(event)=>setPasteMetadata((current)=>({...current,timepoint:event.target.value}))}/></label></div>}</div>
                  )}
                  {activeDatasets.length > 0 && <div className="study-design-panel">
                    <header><div><span>Study behavior</span><strong>Define which biological boundaries downstream methods may cross.</strong></div>{activeDatasets.length > 1 && <label><span>Study / output name</span><input value={studyName} onChange={(event)=>setStudyName(event.target.value)} /></label>}</header>
                    <div className="study-design-options">
                      {([
                        ["independent", "Independent samples", "Collapse and lineages stay within each sample."],
                        ["cohort", "Cohort study", "Samples stay independent; cohort metadata supports comparison and plotting."],
                        ["longitudinal", "Longitudinal donors", "Collapse stays within sample; lineages may span timepoints for the same donor."],
                        ["technical", "Technical replicates", "Files sharing a sample ID may collapse together; lineages stay within that sample."],
                        ["custom", "Custom boundaries", "Choose collapse and lineage scopes independently in pipeline or post-analysis."],
                      ] as Array<[StudyDesign,string,string]>).map(([value,label,description])=><label className={studyDesign===value?"selected":""} key={value}><input type="radio" checked={studyDesign===value} onChange={()=>applyStudyDesign(value)}/><span><strong>{label}</strong><small>{description}</small></span></label>)}
                    </div>
                    {studyDesign === "longitudinal" && <p className="scientific-note"><span>i</span>Give all timepoints from one donor the same donor/subject ID. Different donors remain hard-separated during lineage assignment.</p>}
                  </div>}
                  {inputError && <p className="inline-error" role="alert">{inputError}</p>}
                </section>

                <section className="analysis-card reference-card">
                  <header><span className="card-number">02</span><div><h2>Biological search space</h2><p>Compose germline sources by locus and segment.</p></div>{Object.keys(cellReferences).length > 0 && <button className="reset-button" type="button" onClick={resetReferenceContext}>Reset all to IMGT</button>}</header>
                  <div className="reference-selectors">
                    <label><span>Species / strain</span><select value={species?.name ?? ""} disabled={!pack} onChange={(event) => {
                      setSpeciesName(event.target.value);
                      const next = speciesList.find((item) => item.name === event.target.value);
                      const nextScopes = next ? availableScopes(next) : [];
                      if (!nextScopes.includes(activeScope)) setScope(nextScopes[0] ?? "BCR");
                      resetReferenceContext();
                    }}>{!pack && <option>Loading IMGT references…</option>}{speciesList.map((item) => <option value={item.name} key={item.name}>{friendlySpecies(item.name)}</option>)}</select></label>
                    <label><span>Database</span><select value={compositionMode} disabled={!pack || databaseBusy || Boolean(busyCells.size)} onChange={(event) => void applyReferenceDatabase(event.target.value)}>{databaseOptions.map((option) => <option value={option.id} key={option.id}>{option.label}</option>)}{compositionMode === "mixed" && <option value="mixed" disabled>Mixed sources · configured below</option>}</select></label>
                    <div className="receptor-selector"><span>Receptor</span><div><button className={receptor === "BCR" ? "active" : ""} type="button" disabled={!hasBcr} onClick={() => { if (receptor !== "BCR") resetReferenceContext(); setScope("BCR"); }}>BCR <small>IG</small></button><button className={receptor === "TCR" ? "active" : ""} type="button" disabled={!hasTcr} onClick={() => { if (receptor !== "TCR") resetReferenceContext(); setScope("TCR"); }}>TCR <small>TR</small></button></div></div>
                    <label><span>Chain / locus</span><select value={activeScope} disabled={databaseBusy || Boolean(busyCells.size)} onChange={(event) => setScope(event.target.value as ScopeKey)}>{receptorScopes.map((value) => <option value={value} key={value}>{LOCUS_LABELS[value]}</option>)}</select></label>
                  </div>
                  {packError && <p className="inline-error" role="alert">{packError}</p>}
                  <CompositionSummary databases={usedDatabases} hasUploads={hasUploadedReferences} busy={databaseBusy || Boolean(busyCells.size)} release={pack?.release ?? ""} />
                  {species && <ReferenceCompositionMatrix species={species} scope={activeScope} references={cellReferences} busyCells={busyCells} onSelect={(locus, segment, sourceId) => void selectCellSource(locus, segment, sourceId)} onFile={(locus, segment, file) => void acceptReferenceFile(locus, segment, file)} />}
                  <p className="reference-footnote"><span>i</span><b>{compiled?.loci.join(" + ") || "No locus"}</b> · {databaseLabel}. The database selector applies a preset only where compatible; every matrix cell remains independently selectable or uploadable. In-browser preprocessing assigns V FWR/CDR boundaries and J frame/junction-anchor metadata by transfer from validated, locus-matched IMGT relatives; metadata are retained only when mapped intervals and conserved anchors validate. D and C records are validated and indexed but do not have FWR/CDR boundaries.</p>
                </section>

                <section className="analysis-card settings-card">
                  <header><span className="card-number">03</span><div><h2>Analysis parameters</h2><p>Review sampling, strand, identity, worker, and output settings.</p></div><button className="reset-button" type="button" onClick={() => setShowAdvanced((value) => !value)}>{showAdvanced ? "Hide" : "Show"} controls</button></header>
                  <div className="settings-strip"><span><b>Profile</b> {callingProfileLabel(callingProfile)}</span><span><b>Strand</b> {strand === 0 ? "both" : strand === 1 ? "plus" : "minus"}</span><span><b>Identity floor</b> {Math.round(minimumIdentity * 100)}%</span><span><b>Output</b> {outputStorage === "disk" ? "stream to disk" : outputStorage === "browser" ? "compressed browser index" : "adaptive"}</span><span><b>Input</b> {subsampleEnabled ? `random ${Math.floor(subsampleSize || 0).toLocaleString()}` : "all records"}</span><span><b>Double-D</b> {doubleDMode === "off" ? "off" : doubleDMode === "all" ? "screen all" : `span ≥ ${Math.round(doubleDMinimumSpan)} nt`}</span></div>
                  <div className={`subsample-control ${subsampleEnabled ? "active" : ""}`}><label className="subsample-switch"><input type="checkbox" checked={subsampleEnabled} onChange={(event) => setSubsampleEnabled(event.target.checked)} /><span><b>Analyze a random subsample</b><small>Exact reservoir sampling scans the full stream but retains only the requested records {activeDatasets.length > 1 ? "from each dataset" : ""} in memory and output.</small></span></label>{subsampleEnabled && <div className="subsample-fields"><label><span>{activeDatasets.length > 1 ? "Records per dataset" : "Records to analyze"}</span><CommitNumberInput min="1" step="1000" value={subsampleSize} onCommit={setSubsampleSize} /></label><label><span>Base random seed</span><CommitNumberInput step="1" value={subsampleSeed} onCommit={setSubsampleSeed} /></label></div>}</div>
                  {showAdvanced && <div className="advanced-settings">
                    <label><span>Calling profile</span><select value={callingProfile} onChange={(event) => setCallingProfile(event.target.value as CallingProfile)}><option value="truth_optimized">Truth-optimized · default</option><option value="igblast_balanced">IgBLAST-balanced · agreement + truth constraint</option><option value="igblast_compatible">IgBLAST-agreement · agreement only</option></select></label>
                    <label><span>Search strand</span><select value={strand} onChange={(event) => setStrand(Number(event.target.value) as 0 | 1 | 2)}><option value={0}>Both orientations</option><option value={1}>Plus only</option><option value={2}>Minus only</option></select></label>
                    <label><span>Parallel WASM workers</span><select value={workerCount} onChange={(event) => setWorkerCount(Number(event.target.value))}>{Array.from({ length: browserWorkerLimit() }, (_, index) => index + 1).map((count) => <option value={count} key={count}>{count}{count === recommendedWorkerCount() ? " · recommended" : ""}</option>)}</select></label>
                    <label><span>AIRR results destination</span><select value={outputStorage} onChange={(event) => setOutputStorage(event.target.value as OutputStorageMode)}><option value="auto">Auto · ask to save large results</option><option value="browser">Browser · compressed local index</option><option value="disk">File · save while analyzing</option></select></label>
                    <label className="minimum-slider"><span>Minimum alignment identity <b>{Math.round(minimumIdentity * 100)}%</b></span><input type="range" min="0.45" max="0.9" step="0.01" value={minimumIdentity} onChange={(event) => setMinimumIdentity(Number(event.target.value))} /></label>
                    <p className="scientific-note calling-profile-note"><span>i</span>{callingProfile === "igblast_compatible" ? "Selected solely for agreement with the supplied IgBLAST calls: D +2/−4/−11/−1, minimum 5-nt exact run, 3 candidates; J +2/−4/−13/−1, 2 candidates. Calibrated on simulated human IGH; it is not an IgBLAST implementation." : callingProfile === "igblast_balanced" ? "Maximizes IgBLAST agreement subject to combined V/D/J first-call and fair-scored truth accuracy exceeding IgBLAST on the supplied simulation. It uses the IgBLAST-agreement settings, then removes a D call only when its strongest support is exactly five consecutive matches and j_sequence_start − v_sequence_end ≤ 11 nt. Calibrated on simulated human IGH." : "Default profile selected for simulated ground-truth accuracy. Optional agreement-oriented profiles are never selected automatically."}</p>
                    <label><span>Double-D / VDDJ screening</span><select value={doubleDMode} onChange={(event) => setDoubleDMode(event.target.value as DoubleDScreenMode)}><option value="off">Off · standard V(D)J only</option><option value="long_span">Long inter-V/J spans only</option><option value="all">All eligible D-bearing junctions</option></select></label>
                    <section className={`double-d-parameters ${doubleDMode === "off" ? "disabled" : ""}`}>
                      <div><span className="section-kicker">Rare rearrangement screen</span><h3>{doubleDMode === "off" ? "No double-D screening" : "Two ordered D segments"}</h3><p>{doubleDMode === "off" ? "The original SwiftIG annotation entry point is used and ordinary behavior is unchanged." : "This evidence screen runs after the ordinary V(D)J call. It does not rewrite the main AIRR TSV; supported D1/D2 calls are retained in a separate evidence table and merged into the interactive detail viewer."}</p></div>
                      {doubleDMode !== "off" && <div className="double-d-fields">
                        {doubleDMode === "long_span" && <label><span>Minimum inter-V/J span</span><CommitNumberInput min="0" max="10000" step="1" value={doubleDMinimumSpan} onCommit={setDoubleDMinimumSpan} /><small>nt between the baseline V end and J start</small></label>}
                        <label><span>Exact seed per D</span><CommitNumberInput min="6" max="24" step="1" value={doubleDSeedLength} onCommit={setDoubleDSeedLength} /><small>nt; 11 follows the IgScout tandem-D screen</small></label>
                        <label><span>Single-D Δ trim</span><CommitNumberInput min="0" max="24" step="1" value={doubleDPseudoTrim} onCommit={setDoubleDPseudoTrim} /><small>nt removed for the pseudo-tandem test</small></label>
                        <label><span>Maximum pseudo mismatches</span><CommitNumberInput min="0" max="24" step="1" value={doubleDMaximumPseudoMismatches} onCommit={setDoubleDMaximumPseudoMismatches} /><small>reject a pair at or below this distance</small></label>
                        <label><span>Minimum two-D score gain</span><CommitNumberInput min="0" max="1000" step="1" value={doubleDMinimumScoreGain} onCommit={setDoubleDMinimumScoreGain} /><small>over the best supported single D</small></label>
                      </div>}
                    </section>
                  </div>}
                </section>

                <section className={`analysis-card pipeline-card ${pipeline.enabled ? "active" : ""}`}>
                  <header><span className="card-number">04</span><div><h2>Execution mode</h2><p>Stop after annotation, or run the selected repertoire-scale stages unattended.</p></div><div className="mode-toggle"><button className={!pipeline.enabled ? "active" : ""} type="button" onClick={()=>setPipeline((current)=>({...current,enabled:false}))}>Interactive</button><button className={pipeline.enabled ? "active" : ""} type="button" onClick={()=>setPipeline((current)=>({...current,enabled:true}))}>Pipeline</button></div></header>
                  {!pipeline.enabled ? <div className="pipeline-interactive-note"><span>Manual post-analysis</span><strong>Annotation finishes at Results.</strong><p>Collapse, chimera filtering, selection, lineage assignment, SHM, missing-allele diagnostics, queries, alignments, and trees remain available as explicit steps.</p></div> : <>
                    <div className="pipeline-stage-grid">
                      <article className={pipeline.collapse.enabled?"enabled":""}><label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.collapse.enabled} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,enabled:event.target.checked}}))}/><span><b>1 · Collapse / denoise</b><small>Never crosses the selected boundary.</small></span></label><div className="pipeline-fields"><label><span>Method</span><select value={pipeline.collapse.mode} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,mode:event.target.value as PipelinePlan["collapse"]["mode"]}}))}><option value="exact">Exact deduplication</option><option value="fad">FAD-compatible denoising</option><option value="conservative">Conservative indexed model</option><option value="indel">Indel-aware method D</option></select></label>{pipeline.collapse.mode==="exact"&&<label><span>Exact key</span><select value={pipeline.collapse.key} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,key:event.target.value as PipelinePlan["collapse"]["key"]}}))}><option value="sequence">Full input sequence</option><option value="trimmed">V–J-trimmed sequence</option><option value="cdr3">CDR3 nucleotide</option><option value="rearrangement">V/J calls + CDR3</option></select></label>}<label><span>Boundary</span><select value={pipeline.collapse.scope} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,scope:event.target.value as DatasetScope}}))}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>Unusable V/J or trim</span><select value={pipeline.collapse.unresolvedPolicy} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,unresolvedPolicy:event.target.value as "discard"|"retain"}}))}><option value="discard">Discard · default</option><option value="retain">Retain unchanged</option></select></label></div></article>
                      <article className={pipeline.chimera.enabled?"enabled":""}>
                        <label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.chimera.enabled} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,enabled:event.target.checked}}))}/><span><b>2 · CHMMAIRRa filter</b><small>Consumes the collapsed working set.</small></span></label>
                        <div className="pipeline-fields"><label><span>Segment</span><select value={pipeline.chimera.segment} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,segment:event.target.value as "V"|"J"}}))}><option value="V">V</option><option value="J">J</option></select></label><label><span>Model</span><select value={pipeline.chimera.model} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,model:event.target.value as "auto"|"BW"|"DB"}}))}><option value="auto">Auto · IG BW / TR DB</option><option value="BW">Baum–Welch</option><option value="DB">Discretized Bayesian</option></select></label><label><span>Exclude posterior ≥</span><CommitNumberInput min="0" max="1" step="0.01" value={pipeline.chimera.posteriorThreshold} onCommit={(posteriorThreshold)=>setPipeline((current)=>({...current,chimera:{...current.chimera,posteriorThreshold}}))}/></label><label className="check-line"><input type="checkbox" checked={pipeline.chimera.retainUnevaluated} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,retainUnevaluated:event.target.checked}}))}/><span>Retain unevaluated</span></label><label><span>Reference MSA</span><select value={pipeline.chimera.msaSource} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,msaSource:event.target.value as "selected"|"upload"}}))}><option value="selected">Build from selected references</option><option value="upload">Use uploaded aligned FASTA</option></select></label>{pipeline.chimera.msaSource==="upload"&&<label className="pipeline-file-field"><span>Aligned FASTA</span><input type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(!file)return;void file.text().then((text)=>{prepareReferenceMsa(text);setPipeline((current)=>({...current,chimera:{...current.chimera,uploadedMsa:text,uploadedMsaName:file.name}}));setRunError("");}).catch((error)=>setRunError(error instanceof Error?error.message:String(error)));}}/><small>{pipeline.chimera.uploadedMsaName||"No MSA selected"}</small></label>}</div>
                      </article>
                      <article className={pipeline.selection.enabled?"enabled":""}>
                        <label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.selection.enabled} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,enabled:event.target.checked}}))}/><span><b>3 · Repertoire selection</b><small>Commit study, call, CDR3, and QC filters before lineage assignment.</small></span></label>
                        <div className="pipeline-fields pipeline-selection-fields"><label><span>Dataset / library</span><select value={pipeline.selection.datasetId} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,datasetId:event.target.value}}))}><option value="">Any dataset</option>{activeDatasets.map((dataset)=><option value={dataset.datasetId} key={dataset.datasetId}>{dataset.inputName} · {dataset.datasetId}</option>)}</select></label><label><span>Sample ID</span><input value={pipeline.selection.sampleId} placeholder="any" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,sampleId:event.target.value}}))}/></label><label><span>Donor / subject</span><input value={pipeline.selection.subjectId} placeholder="any" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,subjectId:event.target.value}}))}/></label><label><span>Cohort</span><input value={pipeline.selection.cohort} placeholder="any" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,cohort:event.target.value}}))}/></label><label><span>Timepoint</span><input value={pipeline.selection.timepoint} placeholder="any" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,timepoint:event.target.value}}))}/></label><label><span>Locus</span><input value={pipeline.selection.locus} placeholder="IGH, TRB…" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,locus:event.target.value}}))}/></label><label><span>V call</span><input value={pipeline.selection.vCall} placeholder="contains…" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,vCall:event.target.value}}))}/></label><label><span>J call</span><input value={pipeline.selection.jCall} placeholder="contains…" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,jCall:event.target.value}}))}/></label><label><span>CDR3 nt</span><input value={pipeline.selection.cdr3Nt} placeholder="substring" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,cdr3Nt:event.target.value}}))}/></label><label><span>CDR3 aa</span><input value={pipeline.selection.cdr3Aa} placeholder="substring" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,cdr3Aa:event.target.value}}))}/></label><label><span>Productive</span><select value={pipeline.selection.productive} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,productive:event.target.value as "any"|"yes"|"no"}}))}><option value="any">Either</option><option value="yes">Productive</option><option value="no">Non-productive</option></select></label><label><span>CDR3 assigned</span><select value={pipeline.selection.hasCdr3} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,hasCdr3:event.target.value as "any"|"yes"|"no"}}))}><option value="any">Either</option><option value="yes">Required</option><option value="no">Absent</option></select></label><label><span>Double-D</span><select value={pipeline.selection.doubleD} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,doubleD:event.target.value as "any"|"positive"|"negative"}}))}><option value="any">Either</option><option value="positive">Positive only</option><option value="negative">Exclude positive</option></select></label></div>
                      </article>
                      <article className={pipeline.lineage.enabled?"enabled":""}><label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.lineage.enabled} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,enabled:event.target.checked}}))}/><span><b>4 · Lineage assignment</b><small>May span samples only within the selected boundary.</small></span></label><div className="pipeline-fields"><label><span>Boundary</span><select value={pipeline.lineage.scope} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,scope:event.target.value as DatasetScope}}))}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>CDR3 nt identity</span><CommitNumberInput min="0.7" max="1" step="0.01" value={pipeline.lineage.identity} onCommit={(identity)=>setPipeline((current)=>({...current,lineage:{...current.lineage,identity}}))}/></label><label><span>V/J level</span><select value={pipeline.lineage.resolution} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,resolution:event.target.value as "gene"|"allele"}}))}><option value="gene">Gene</option><option value="allele">Allele</option></select></label><label><span>Multiple-call policy</span><select value={pipeline.lineage.ambiguity} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,ambiguity:event.target.value as PipelinePlan["lineage"]["ambiguity"]}}))}><option value="overlap">Any V/J overlap</option><option value="top">Top call only</option><option value="strict">Identical call sets</option></select></label><label className="check-line"><input type="checkbox" checked={pipeline.lineage.productiveOnly} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,productiveOnly:event.target.checked}}))}/><span>Productive only</span></label></div></article>
                      <article className={pipeline.shm.enabled||pipeline.missingAlleles.enabled?"enabled":""}><div className="pipeline-diagnostic-switches"><label><input type="checkbox" checked={pipeline.shm.enabled} onChange={(event)=>setPipeline((current)=>({...current,shm:{...current.shm,enabled:event.target.checked}}))}/><span><b>5 · SHM summary</b><small>On the final working set.</small></span></label><label><input type="checkbox" checked={pipeline.missingAlleles.enabled} onChange={(event)=>setPipeline((current)=>({...current,missingAlleles:{enabled:event.target.checked},lineage:event.target.checked?{...current.lineage,enabled:true}:current.lineage}))}/><span><b>6 · Missing-allele hints</b><small>Enabling this also enables lineage assignment.</small></span></label></div>{pipeline.shm.enabled&&<div className="pipeline-fields"><label><span>SHM metric</span><select value={pipeline.shm.metric} onChange={(event)=>setPipeline((current)=>({...current,shm:{...current.shm,metric:event.target.value as PipelinePlan["shm"]["metric"]}}))}><option value="vNtRate">V nucleotide mutation rate</option><option value="vNtMutations">V nucleotide mutation count</option><option value="vAaRate">V amino-acid replacement rate</option><option value="vAaReplacements">V amino-acid replacement count</option><option value="synonymous">V synonymous mutation count</option><option value="cdrNtRate">CDR nucleotide mutation rate</option><option value="frameworkNtRate">Framework nucleotide mutation rate</option></select></label></div>}</article>
                    </div>
                    <p className="scientific-note"><span>i</span>The pipeline runs annotation → collapse → chimera exclusion → repertoire selection → lineage assignment → diagnostics in that order. Every stage receives the retained set from the previous stage. Targeted queries, lineage alignments, and phylogenies remain on-demand because they require a selected target.</p>
                  </>}
                </section>
              </div>

              <aside className="run-summary">
                <span className="section-kicker">Run manifest</span><h2>{activeInput ? `${activeDatasets.length.toLocaleString()} dataset${activeDatasets.length === 1 ? "" : "s"}` : "Awaiting data"}</h2><p>{activeInput ? activeInputName : "Add an input to configure the run."}</p>
                <dl><div><dt>Records</dt><dd>{subsampleEnabled ? `${Math.floor(subsampleSize).toLocaleString()} per dataset` : knownInputCount === null ? "Count during streaming" : knownInputCount.toLocaleString()}</dd></div><div><dt>Study</dt><dd>{studyDesign.replace(/^./,(value)=>value.toUpperCase())}</dd></div><div><dt>Species</dt><dd>{species ? friendlySpecies(species.name) : "Loading…"}</dd></div><div><dt>Sources</dt><dd>{databaseLabel}</dd></div><div><dt>Search</dt><dd>{LOCUS_LABELS[activeScope]}</dd></div><div><dt>Calling profile</dt><dd>{callingProfileLabel(callingProfile)}</dd></div><div><dt>V / D / J / C</dt><dd>{compiled ? `${compiled.counts.V} / ${compiled.counts.D} / ${compiled.counts.J} / ${compiled.counts.C}` : "—"}</dd></div><div><dt>Non-IMGT cells</dt><dd>{activeReferenceEntries.length ? activeReferenceEntries.map(([key]) => key.replace(":", " ")).join(" · ") : "None"}</dd></div><div><dt>Compute</dt><dd>{workerCount} workers · datasets sequential</dd></div><div><dt>Post-analysis</dt><dd>{pipeline.enabled ? [pipeline.collapse.enabled&&"collapse",pipeline.chimera.enabled&&"chimera",pipeline.selection.enabled&&"selection",pipeline.lineage.enabled&&"lineages",pipeline.shm.enabled&&"SHM",pipeline.missingAlleles.enabled&&"allele hints"].filter(Boolean).join(" → ")||"No stages selected" : "Interactive"}</dd></div><div><dt>Double-D</dt><dd>{doubleDMode === "off" ? "Off · standard path" : doubleDMode === "all" ? "Screen all eligible junctions" : `Screen inter-V/J spans ≥ ${Math.round(doubleDMinimumSpan)} nt`}</dd></div><div><dt>Storage</dt><dd>{outputStorage === "auto" ? "Adaptive streaming" : outputStorage === "disk" ? "Direct AIRR file" : "Compressed local index"}</dd></div></dl>
                {runError && <p className="run-error" role="alert">{runError}</p>}
                <button className="analyze-button" type="button" disabled={!activeInput || !compiled || Boolean(packError) || databaseBusy || Boolean(busyCells.size)} onClick={requestRun}><span>{databaseBusy || busyCells.size ? "Validating references…" : pipeline.enabled ? "Run annotation + pipeline" : "Analyze with SwiftIG"}</span><b>→</b></button>
                <p className="privacy-copy"><span>i</span> Query sequences, uploaded germlines, and AIRR results are processed in this browser; Swig does not transmit them. A remotely hosted alternative database is requested from its named provider only when selected.</p>
              </aside>
            </div>
          )}
        </main>
      )}

      {page === "results" && session && <ResultsPage session={session} onNewAnalysis={() => navigate("analyze")} />}
      {(pendingLoadedSession||sessionLoadError)&&<div className="output-modal-backdrop" role="presentation"><section className="output-modal session-load-modal" role="dialog" aria-modal="true" aria-labelledby="session-load-title"><button className="output-modal-close" type="button" disabled={loadingSession} onClick={()=>{setPendingLoadedSession(null);setSessionLoadError("");}}>×</button><span className="output-direction">SESSION · LINKED AIRR DATA</span><h2 id="session-load-title">{pendingLoadedSession?"Select the AIRR table linked to this session.":"Session could not be loaded."}</h2>{pendingLoadedSession?<><p>The session stores references, options, masks, counts, lineage assignments, plots, and sparse double-D evidence. It deliberately does not duplicate the main AIRR table.</p><div className="output-flow"><div><span>Saved analysis</span><strong>{pendingLoadedSession.analysis.inputName}</strong><small>{pendingLoadedSession.linkedAirr.records.toLocaleString()} records · fingerprint {pendingLoadedSession.linkedAirr.fingerprint.slice(0,12)}…</small></div><b>+</b><div className="destination"><span>Required linked file</span><strong>{pendingLoadedSession.linkedAirr.name}</strong><small>AIRR TSV or TSV.gz; content is verified before restoration</small></div></div>{loadingSession?<div className="post-progress"><div><span>{sessionLoadProgress.stage}</span><strong>{sessionLoadProgress.total?`${Math.min(100,sessionLoadProgress.records/sessionLoadProgress.total*100).toFixed(1)}%`:"working"}</strong></div><progress max={Math.max(1,sessionLoadProgress.total)} value={sessionLoadProgress.records}/><small>{sessionLoadProgress.records.toLocaleString()} / {sessionLoadProgress.total.toLocaleString()} records</small></div>:<button className="output-save-primary" type="button" onClick={()=>linkedAirrInputRef.current?.click()}><span>Choose linked AIRR TSV</span><b>Open →</b></button>}</>:null}{sessionLoadError?<p className="run-error" role="alert">{sessionLoadError}</p>:null}</section></div>}
      {outputPrompt && activeInput && <div className="output-modal-backdrop" role="presentation"><section className="output-modal" role="dialog" aria-modal="true" aria-labelledby="output-dialog-title">
        <button className="output-modal-close" type="button" onClick={() => setOutputPrompt(false)} aria-label="Cancel output selection">×</button>
        <span className="output-direction">OUTPUT · SAVE RESULTS</span>
        <h2 id="output-dialog-title">Choose where Swig will <em>write the AIRR output.</em></h2>
        <p>The next system window is a <b>Save As</b> dialog. You are naming a new results file—this is not another sequence import.</p>
        <div className="output-flow"><div><span>Input already selected</span><strong>{activeDatasets.length === 1 ? activeInput.name : `${activeDatasets.length} datasets · ${activeInputName}`}</strong></div><b>→</b><div className="destination"><span>New output file</span><strong>{outputName(activeInputName)}</strong><small>Written incrementally during analysis</small></div></div>
        <div className="output-modal-actions"><button className="output-save-primary" type="button" onClick={() => void run("disk")}><span>Choose output file &amp; start</span><b>Save AIRR →</b></button><button type="button" onClick={() => void run("browser")}><span>Keep output in browser instead</span><small>Compressed local index; download after the run</small></button></div>
        <p className="output-safety"><span>i</span> Query sequences remain in this browser and are not transmitted by Swig.</p>
      </section></div>}
      <footer className="site-footer"><Brand /><p>Swig 0.16.0 · SwiftIG WebAssembly interface · research software · validate study-critical calls independently.</p><div><a href="https://github.com/MurrellGroup/swiftig" target="_blank" rel="noreferrer">Source ↗</a><a href="https://www.imgt.org/" target="_blank" rel="noreferrer">IMGT ↗</a><a href="https://docs.airr-community.org/" target="_blank" rel="noreferrer">AIRR ↗</a></div></footer>
    </div>
  );
}
