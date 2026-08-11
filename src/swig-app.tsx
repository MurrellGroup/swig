import {
  ChangeEvent,
  DragEvent,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import { AlignmentViewer } from "./alignment-view";
import type { GermlinePreprocessReport, MetadataAllele } from "./germline-preprocess";
import { preprocessGermlinesInWorker } from "./germline-preprocess-client";
import { RepertoireDashboard } from "./repertoire-charts";
import { PostAnalysisWorkbench } from "./post-analysis";
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
  type DirectAirrOutput,
  type ResultFacets,
  type ResultFilters,
  type ResultPage,
} from "./result-store";
import { runSwiftIg } from "./swiftig-runtime";

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
  minimumIdentity: number;
  strand: 0 | 1 | 2;
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

function AppHeader({ page, hasResults, onNavigate }: {
  page: AppPage;
  hasResults: boolean;
  onNavigate: (page: AppPage) => void;
}) {
  return (
    <header className="app-header">
      <button className="brand-button" type="button" onClick={() => onNavigate("home")}><Brand /></button>
      <nav aria-label="Primary navigation">
        <button className={page === "home" ? "active" : ""} type="button" onClick={() => onNavigate("home")}>Overview</button>
        <button className={page === "analyze" ? "active" : ""} type="button" onClick={() => onNavigate("analyze")}>Analyze</button>
        {hasResults && <button className={page === "results" ? "active" : ""} type="button" onClick={() => onNavigate("results")}>Results</button>}
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

function AnalysisProgress({ stage, value, onCancel }: { stage: string; value: number; onCancel: () => void }) {
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
      <div className="progress-copy"><p className="eyebrow"><span>SwiftIG is running locally</span></p><h2>{stage}</h2><p>The page will move to Results automatically when the local AIRR index is ready.</p><div className="main-progress"><i style={{ width: `${percent}%` }} /></div><ol>{phases.map((phase, index) => {
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
  const mapSegments = ["v", "d", "j", "c"].flatMap((segment) => {
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
        <div><span className="section-kicker">Selected rearrangement</span><h2>{row.sequence_id}</h2><div className="detail-tags"><span>{row.locus || "unassigned"}</span><span className={row.productive === "T" ? "good" : "warn"}>{row.productive === "T" ? "Productive" : "Non-productive"}</span>{isotype && <span className="isotype-tag">{isotype}</span>}{row.rev_comp === "T" && <span>Reverse complement</span>}</div></div>
        <button className="close-detail" type="button" onClick={onClose}>Close <span>×</span></button>
      </header>

      <section className="rearrangement-overview">
        <div className="overview-heading"><h3>Rearrangement map</h3><span>{sequenceLength.toLocaleString()} nt</span></div>
        <div className="rearrangement-track">
          <i className="baseline" />
          {mapSegments.map((segment) => <span key={segment.segment} className={`map-segment ${segment.segment}`} style={{ left: `${segment.left}%`, width: `${segment.width}%` }} title={segment.call}>{segment.segment.toUpperCase()}</span>)}
        </div>
        <div className="call-grid">
          {["v", "d", "j", "c"].map((segment) => {
            const selected = splitCalls(row[`${segment}_call`] || "");
            const alternatives = parseAlternatives(row[`${segment}_alternatives`] || "");
            return <div key={segment}><span>{segment.toUpperCase()} call</span><strong>{selected[0] || "—"}</strong><small>{row[`${segment}_identity`] ? `${(Number(row[`${segment}_identity`]) * 100).toFixed(1)}% identity` : "not assigned"}{selected.length > 1 ? ` · ${selected.length} co-optimal` : alternatives.length ? ` · ${alternatives.length} near-tied` : ""}</small></div>;
          })}
        </div>
      </section>

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

function ResultsPage({ session, onNewAnalysis }: { session: ResultSession; onNewAnalysis: () => void }) {
  const [view, setView] = useState<"repertoire" | "sequences" | "post">(session.total <= 3 ? "sequences" : "repertoire");
  const [postOpened, setPostOpened] = useState(false);
  const [filters, setFilters] = useState<ResultFilters>({ ...EMPTY_FILTERS });
  const [page, setPage] = useState(0);
  const [results, setResults] = useState<ResultPage>({ rows: [], hasMore: false, totalMatches: session.total, scanned: 0 });
  const [searching, setSearching] = useState(true);
  const [scanCount, setScanCount] = useState(0);
  const [selected, setSelected] = useState<AirrIndexRecord | null>(null);
  const [detail, setDetail] = useState<AirrRow | null>(null);
  const [downloading, setDownloading] = useState(false);
  const [downloadError, setDownloadError] = useState("");
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

  async function downloadAll() {
    setDownloading(true);
    setDownloadError("");
    try {
      const name = outputName(session.inputName);
      const picker = savePicker();
      if (picker) {
        const handle = await picker.call(window, {
          suggestedName: name,
          types: [{ description: "AIRR rearrangement table", accept: { "text/tab-separated-values": [".tsv"] } }],
        });
        const writable = await handle.createWritable();
        try {
          await session.store.writeAirr((part) => writable.write(part));
          await writable.close();
        } catch (error) {
          await writable.abort?.();
          throw error;
        }
      } else {
        if (session.streamedDirectly) {
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
      }
    } catch (error) {
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setDownloadError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      setDownloading(false);
    }
  }

  const filtered = Object.entries(filters).some(([key, value]) => key.startsWith("min") ? Number(value) > 0 : Boolean(value));
  const pageStart = page * PAGE_SIZE + (results.rows.length ? 1 : 0);
  const pageEnd = page * PAGE_SIZE + results.rows.length;

  return (
    <main className="results-page">
      <section className="results-hero">
        <WorkflowStepper active={3} />
        <div className="results-title"><div><p className="eyebrow"><span>Analysis complete · {session.seconds.toFixed(2)} s · {Math.round(session.total / Math.max(session.seconds, 0.001)).toLocaleString()} reads/s</span></p><h1>{session.total.toLocaleString()} analyzed<br /><em>rearrangements.</em></h1><p>{friendlySpecies(session.species)} · {LOCUS_LABELS[session.scope]} · {session.workers} WASM worker{session.workers === 1 ? "" : "s"} · {bytes(session.outputBytes)} AIRR{session.subsampleSize ? ` · exact random sample from ${session.inputTotal.toLocaleString()} input records (seed ${session.subsampleSeed})` : ""}</p></div><div className="results-actions"><button className="download-primary" type="button" onClick={() => void downloadAll()} disabled={downloading}>{downloading ? "Writing AIRR…" : session.streamedDirectly ? "Save another AIRR copy" : "Download AIRR TSV"}<span>↓</span></button><button type="button" onClick={onNewAnalysis}>New analysis</button>{downloadError && <small role="alert">{downloadError}</small>}</div></div>
        <div className="result-summary">
          <article><span>V + J assigned</span><strong>{session.summary.assigned.toLocaleString()}</strong><small>{percentage(session.summary.assigned, session.total)} of input</small></article>
          <article><span>Productive</span><strong>{session.summary.productive.toLocaleString()}</strong><small>{percentage(session.summary.productive, session.total)} of input</small></article>
          <article><span>CDR3 called</span><strong>{session.summary.withCdr3.toLocaleString()}</strong><small>{percentage(session.summary.withCdr3, session.total)} of input</small></article>
          <article><span>Loci observed</span><strong>{session.facets.loci.length}</strong><small>{session.facets.loci.map((item) => item.value).join(" · ") || "none"}</small></article>
        </div>
      </section>

      <nav className="results-view-tabs" aria-label="Results view"><button className={view === "repertoire" ? "active" : ""} type="button" onClick={() => setView("repertoire")}><span>Repertoire</span><small>Figures + composition</small></button><button className={view === "sequences" ? "active" : ""} type="button" onClick={() => setView("sequences")}><span>Sequences</span><small>Filter + inspect calls</small></button><button className={view === "post" ? "active" : ""} type="button" onClick={() => { setPostOpened(true); setView("post"); }}><span>Post-analysis</span><small>Deduplicate + lineages + trees</small></button></nav>

      {view === "repertoire" ? <RepertoireDashboard store={session.store} loci={session.facets.loci} inputName={session.inputName} /> : view === "sequences" ? <>

      <section className="explorer-shell">
        <aside className="filter-panel">
          <div className="filter-heading"><div><span className="section-kicker">Local query</span><h2>Filter results</h2></div>{filtered && <button type="button" onClick={() => { setFilters({ ...EMPTY_FILTERS }); setPage(0); }}>Clear</button>}</div>
          <label className="filter-field"><span>Sequence ID contains</span><input type="search" value={filters.sequenceId} placeholder="e.g. clonotype_104" onChange={(event) => updateFilter("sequenceId", event.target.value)} /></label>
          <label className="filter-field"><span>CDR3 substring <small>nt or AA</small></span><input className="monospace" type="search" value={filters.cdr3} placeholder="CARDR / TGTGCC…" onChange={(event) => updateFilter("cdr3", event.target.value)} /></label>
          <div className="filter-row">
            <label className="filter-field"><span>Locus</span><select value={filters.locus} onChange={(event) => updateFilter("locus", event.target.value)}><option value="">Any locus</option>{session.facets.loci.map((item) => <option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>
            <label className="filter-field"><span>Productivity</span><select value={filters.productive} onChange={(event) => updateFilter("productive", event.target.value)}><option value="">Either</option><option value="T">Productive</option><option value="F">Non-productive</option></select></label>
          </div>
          {[{ key: "vCall", label: "V allele", values: session.facets.vCalls }, { key: "dCall", label: "D allele", values: session.facets.dCalls }, { key: "jCall", label: "J allele", values: session.facets.jCalls }, { key: "cCall", label: "C allele", values: session.facets.cCalls }, { key: "isotype", label: "Isotype / constant class", values: session.facets.isotypes }].map((field) => <label className="filter-field" key={field.key}><span>{field.label}</span><select value={filters[field.key as keyof ResultFilters] as string} onChange={(event) => updateFilter(field.key as "vCall" | "dCall" | "jCall" | "cCall" | "isotype", event.target.value)}><option value="">Any call</option>{field.values.map((item) => <option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>)}
          <details className="filter-advanced"><summary>Identity, junction + QC</summary><div>
            <label className="identity-filter"><span>Minimum V identity <b>{Math.round(filters.minVIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minVIdentity} onChange={(event) => updateFilter("minVIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum D identity <b>{Math.round(filters.minDIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minDIdentity} onChange={(event) => updateFilter("minDIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum J identity <b>{Math.round(filters.minJIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minJIdentity} onChange={(event) => updateFilter("minJIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum C identity <b>{Math.round(filters.minCIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minCIdentity} onChange={(event) => updateFilter("minCIdentity", Number(event.target.value))} /></label>
            <div className="length-filters">
              <label><span>CDR3 AA min</span><input type="number" min="0" max="250" value={filters.minCdr3AaLength || ""} placeholder="Any" onChange={(event) => updateFilter("minCdr3AaLength", Math.max(0, Number(event.target.value) || 0))} /></label>
              <label><span>CDR3 AA max</span><input type="number" min="0" max="250" value={filters.maxCdr3AaLength || ""} placeholder="Any" onChange={(event) => updateFilter("maxCdr3AaLength", Math.max(0, Number(event.target.value) || 0))} /></label>
            </div>
            <div className="qc-filter-grid">
              <label><span>VJ frame</span><select value={filters.vjInFrame} onChange={(event) => updateFilter("vjInFrame", event.target.value)}><option value="">Either</option><option value="T">In frame</option><option value="F">Out of frame</option></select></label>
              <label><span>Stop codon</span><select value={filters.stopCodon} onChange={(event) => updateFilter("stopCodon", event.target.value)}><option value="">Either</option><option value="F">Absent</option><option value="T">Present</option></select></label>
              <label><span>Completeness</span><select value={filters.completeVdj} onChange={(event) => updateFilter("completeVdj", event.target.value)}><option value="">Either</option><option value="T">Complete V(D)J</option><option value="F">Partial</option></select></label>
              <label><span>Orientation</span><select value={filters.revComp} onChange={(event) => updateFilter("revComp", event.target.value)}><option value="">Either</option><option value="F">Forward</option><option value="T">Reverse-comp.</option></select></label>
            </div>
            <label className="check-filter"><input type="checkbox" checked={filters.hasD} onChange={(event) => updateFilter("hasD", event.target.checked)} /><span>Require a D assignment</span></label>
            <label className="check-filter"><input type="checkbox" checked={filters.hasCdr3} onChange={(event) => updateFilter("hasCdr3", event.target.checked)} /><span>Require a CDR3 call</span></label>
          </div></details>
          <p className="index-note"><span>i</span> Exact gene and locus filters use browser-local indexes. Substring filters scan candidate records on demand within the browser.</p>
        </aside>

        <div className="result-browser">
          <header className="browser-heading"><div><span className="section-kicker">AIRR records</span><h2>{searching ? "Searching local index…" : results.totalMatches !== null ? `${results.totalMatches.toLocaleString()} matching records` : `${(pageEnd + (results.hasMore ? 1 : 0)).toLocaleString()}+ matching records`}</h2><p>{searching && scanCount ? `${scanCount.toLocaleString()} candidates scanned` : results.rows.length ? `Showing ${pageStart.toLocaleString()}–${pageEnd.toLocaleString()}` : "Adjust filters to broaden the query."}</p></div><span className="scale-mode">{session.total <= 3 ? "detail mode" : session.total >= 100000 ? "large-run index" : "paged index"}</span></header>
          <div className={`results-table-wrap ${searching ? "loading" : ""}`}>
            <table className="results-table">
              <thead><tr><th>Sequence</th><th>Locus</th><th>V call</th><th>D call</th><th>J call</th><th>Isotype</th><th>CDR3 AA</th><th>Productive</th><th /></tr></thead>
              <tbody>{results.rows.map((row) => <tr className={selected?.ordinal === row.ordinal ? "selected" : ""} key={row.ordinal} tabIndex={0} onClick={() => openRecord(row)} onKeyDown={(event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); openRecord(row); } }}>
                <td><strong title={row.sequenceId}>{row.sequenceId}</strong><small>#{(row.ordinal + 1).toLocaleString()}</small></td>
                <td><span className="locus-pill">{row.locus || "—"}</span></td>
                <td title={row.vCall}>{row.vCall || <i>—</i>}</td><td title={row.dCall}>{row.dCall || <i>—</i>}</td><td title={row.jCall}>{row.jCall || <i>—</i>}</td><td>{row.isotype || <i>—</i>}</td>
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
      </> : null}
      {postOpened && <div hidden={view !== "post"}><PostAnalysisWorkbench store={session.store} references={session.references} scope={session.scope} loci={session.facets.loci} inputName={session.inputName} workers={session.workers} minimumIdentity={session.minimumIdentity} strand={session.strand} onInspect={(ordinal) => void inspectOrdinal(ordinal)} /></div>}
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
  const [fileInput, setFileInput] = useState<InputData | null>(null);
  const [pasteText, setPasteText] = useState("");
  const [inputError, setInputError] = useState("");
  const [cellReferences, setCellReferences] = useState<ReferenceCellMap>({});
  const [busyCells, setBusyCells] = useState<Set<string>>(new Set());
  const [databaseBusy, setDatabaseBusy] = useState(false);
  const [minimumIdentity, setMinimumIdentity] = useState(0.6);
  const [strand, setStrand] = useState<0 | 1 | 2>(0);
  const [workerCount, setWorkerCount] = useState(recommendedWorkerCount);
  const [outputStorage, setOutputStorage] = useState<OutputStorageMode>("auto");
  const [subsampleEnabled, setSubsampleEnabled] = useState(false);
  const [subsampleSize, setSubsampleSize] = useState(10_000);
  const [subsampleSeed, setSubsampleSeed] = useState(1);
  const [outputPrompt, setOutputPrompt] = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState({ stage: "Preparing analysis", value: 0 });
  const [runError, setRunError] = useState("");
  const [session, setSession] = useState<ResultSession | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  const abortRef = useRef<AbortController | null>(null);
  const databaseRequestRef = useRef(0);
  const cellRequestRef = useRef<Record<string, number>>({});
  const referenceCacheRef = useRef<Map<string, ReferenceOverride>>(new Map());

  useEffect(() => {
    loadReferencePack().then(setPack).catch((error) => setPackError(error instanceof Error ? error.message : String(error)));
    if ("serviceWorker" in navigator && window.isSecureContext) void registerDownloadWorker().catch(() => undefined);
  }, []);

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
  const activeInput = inputSource === "upload" ? fileInput : pasteInput;

  function navigate(next: AppPage) {
    if (next === "results" && !session) next = "analyze";
    setPage(next);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  async function acceptInputFile(file: File) {
    try {
      setInputError("");
      setRunError("");
      setFileInput(await inspectFile(file));
    } catch (error) {
      setFileInput(null);
      setInputError(error instanceof Error ? error.message : String(error));
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
    if (!activeInput || !compiled || !species) return;
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
    const selectedCount = subsampleEnabled ? Math.floor(subsampleSize) : null;
    const wantsDisk = outputStorage === "disk" || (outputStorage === "auto" && likelyLargeInput(activeInput, selectedCount));
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
    if (!activeInput || !compiled || !species) return;
    setOutputPrompt(false);

    let directOutput: DirectAirrOutput | undefined;
    if (outputDestination === "disk") {
      const picker = savePicker();
      if (!picker) {
        setRunError("Direct-to-disk streaming is unavailable in this browser. Use Auto/browser storage or a Chromium-based browser.");
        return;
      }
      try {
        const handle = await picker.call(window, {
          suggestedName: outputName(activeInput.name),
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
    try {
      const result = await runSwiftIg({
        query: activeInput.source,
        format: activeInput.formatCode,
        references: compiled,
        minimumIdentity,
        strand,
        workers: workerCount,
        countHint: activeInput.count,
        subsample: subsampleEnabled ? { size: Math.floor(subsampleSize), seed: Math.trunc(subsampleSeed) } : undefined,
        signal: controller.signal,
        onProgress: (stage, value) => setProgress({ stage, value }),
        onBatch: async (batch) => {
          await store.appendBatch(batch.header, batch.body);
          setProgress((current) => ({
            ...current,
            stage: batch.total === null
              ? `Committed ${batch.processed.toLocaleString()} AIRR records`
              : `Committed ${batch.processed.toLocaleString()} of ${batch.total.toLocaleString()} AIRR records`,
          }));
        },
      });
      await store.finalize();
      setProgress({ stage: "Results ready", value: 1 });
      setSession({
        id: Date.now(),
        store,
        total: result.count,
        seconds: (performance.now() - started) / 1000,
        inputName: activeInput.name,
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
        minimumIdentity,
        strand,
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
    }
  }

  const hasBcr = scopes.includes("BCR");
  const hasTcr = scopes.includes("TCR");

  return (
    <div className="site-shell">
      <AppHeader page={page} hasResults={Boolean(session)} onNavigate={navigate} />
      {page === "home" && <LandingPage references={pack?.species.length ?? null} onStart={() => navigate("analyze")} onDemo={chooseDemo} />}

      {page === "analyze" && (
        <main className="analysis-page">
          <section className="analysis-intro"><WorkflowStepper active={running ? 2 : 1} /><div><p className="eyebrow"><span>Analysis workspace</span></p><h1>{running ? "Calling rearrangements…" : "Configure an annotation run."}</h1><p>{running ? "SwiftIG is processing bounded batches and writing AIRR records into a browser-local index." : "Provide sequences, select the biological search space, and specify any germline replacements."}</p></div></section>

          {running ? <AnalysisProgress stage={progress.stage} value={progress.value} onCancel={() => abortRef.current?.abort()} /> : (
            <div className="analysis-layout">
              <div className="analysis-forms">
                <section className="analysis-card input-card">
                  <header><span className="card-number">01</span><div><h2>Sequence input</h2><p>Upload a file or paste records directly.</p></div>{activeInput && <span className="ready-tag">Ready</span>}</header>
                  <div className="source-tabs" role="tablist"><button className={inputSource === "upload" ? "active" : ""} type="button" onClick={() => setInputSource("upload")}>Upload file</button><button className={inputSource === "paste" ? "active" : ""} type="button" onClick={() => setInputSource("paste")}>Paste sequences</button></div>
                  <input ref={inputRef} className="visually-hidden" type="file" accept=".fa,.fasta,.fna,.fas,.fq,.fastq,.tsv,.csv,.txt,.gz" onChange={(event: ChangeEvent<HTMLInputElement>) => {
                    const file = event.target.files?.[0];
                    if (file) void acceptInputFile(file);
                    event.target.value = "";
                  }} />
                  {inputSource === "upload" ? fileInput ? (
                    <div className="loaded-input"><span className="file-glyph">↳</span><div><strong>{fileInput.name}</strong><span>{fileInput.count === null ? "Counted during analysis" : `${fileInput.count.toLocaleString()} sequences`} · {fileInput.format} · {bytes(fileInput.size)}</span></div><button type="button" onClick={() => setFileInput(null)}>Remove</button></div>
                  ) : (
                    <button className="input-dropzone" type="button" onClick={() => inputRef.current?.click()} onDragOver={(event: DragEvent) => event.preventDefault()} onDrop={(event: DragEvent) => {
                      event.preventDefault();
                      const file = event.dataTransfer.files?.[0];
                      if (file) void acceptInputFile(file);
                    }}><span>＋</span><strong>Drop sequence data here</strong><small>.fasta(.gz) · .fastq(.gz) · AIRR .tsv(.gz)</small><i>Choose file</i></button>
                  ) : (
                    <div className="paste-input"><textarea spellCheck={false} value={pasteText} onChange={(event) => setPasteText(event.target.value)} placeholder={">sequence_1\nCAGGTGCAGCTGGTG...\n\n—or—\n\nsequence_id\tsequence\nread_1\tCAGGTGCAGCTGGTG..."} /><footer><span>{pasteInput ? `${pasteInput.count?.toLocaleString()} ${pasteInput.format} records detected` : pasteText.trim() ? "Waiting for valid FASTA, FASTQ, or AIRR…" : "Nothing pasted yet"}</span><button type="button" onClick={chooseDemo}>Insert demo</button></footer></div>
                  )}
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
                  <div className="settings-strip"><span><b>Strand</b> {strand === 0 ? "both" : strand === 1 ? "plus" : "minus"}</span><span><b>Identity floor</b> {Math.round(minimumIdentity * 100)}%</span><span><b>Output</b> {outputStorage === "disk" ? "stream to disk" : outputStorage === "browser" ? "compressed browser index" : "adaptive"}</span><span><b>Input</b> {subsampleEnabled ? `random ${Math.floor(subsampleSize || 0).toLocaleString()}` : "all records"}</span></div>
                  <div className={`subsample-control ${subsampleEnabled ? "active" : ""}`}><label className="subsample-switch"><input type="checkbox" checked={subsampleEnabled} onChange={(event) => setSubsampleEnabled(event.target.checked)} /><span><b>Analyze a random subsample</b><small>Exact reservoir sampling scans the full stream but retains only the requested records in memory and output.</small></span></label>{subsampleEnabled && <div className="subsample-fields"><label><span>Records to analyze</span><input type="number" min="1" step="1000" value={subsampleSize} onChange={(event) => setSubsampleSize(Math.max(1, Number(event.target.value) || 1))} /></label><label><span>Random seed</span><input type="number" step="1" value={subsampleSeed} onChange={(event) => setSubsampleSeed(Number(event.target.value) || 0)} /></label></div>}</div>
                  {showAdvanced && <div className="advanced-settings">
                    <label><span>Search strand</span><select value={strand} onChange={(event) => setStrand(Number(event.target.value) as 0 | 1 | 2)}><option value={0}>Both orientations</option><option value={1}>Plus only</option><option value={2}>Minus only</option></select></label>
                    <label><span>Parallel WASM workers</span><select value={workerCount} onChange={(event) => setWorkerCount(Number(event.target.value))}>{Array.from({ length: browserWorkerLimit() }, (_, index) => index + 1).map((count) => <option value={count} key={count}>{count}{count === recommendedWorkerCount() ? " · recommended" : ""}</option>)}</select></label>
                    <label><span>AIRR results destination</span><select value={outputStorage} onChange={(event) => setOutputStorage(event.target.value as OutputStorageMode)}><option value="auto">Auto · ask to save large results</option><option value="browser">Browser · compressed local index</option><option value="disk">File · save while analyzing</option></select></label>
                    <label className="minimum-slider"><span>Minimum alignment identity <b>{Math.round(minimumIdentity * 100)}%</b></span><input type="range" min="0.45" max="0.9" step="0.01" value={minimumIdentity} onChange={(event) => setMinimumIdentity(Number(event.target.value))} /></label>
                  </div>}
                </section>
              </div>

              <aside className="run-summary">
                <span className="section-kicker">Run manifest</span><h2>{activeInput ? subsampleEnabled ? `${Math.floor(subsampleSize).toLocaleString()}-read sample` : activeInput.count === null ? "Large file" : `${activeInput.count?.toLocaleString()} sequences` : "Awaiting data"}</h2><p>{activeInput?.name ?? "Add an input to configure the run."}</p>
                <dl><div><dt>Species</dt><dd>{species ? friendlySpecies(species.name) : "Loading…"}</dd></div><div><dt>Sources</dt><dd>{databaseLabel}</dd></div><div><dt>Search</dt><dd>{LOCUS_LABELS[activeScope]}</dd></div><div><dt>V / D / J / C</dt><dd>{compiled ? `${compiled.counts.V} / ${compiled.counts.D} / ${compiled.counts.J} / ${compiled.counts.C}` : "—"}</dd></div><div><dt>Non-IMGT cells</dt><dd>{activeReferenceEntries.length ? activeReferenceEntries.map(([key]) => key.replace(":", " ")).join(" · ") : "None"}</dd></div><div><dt>Compute</dt><dd>{workerCount} workers · bounded queue</dd></div><div><dt>Storage</dt><dd>{outputStorage === "auto" ? "Adaptive streaming" : outputStorage === "disk" ? "Direct AIRR file" : "Compressed local index"}</dd></div></dl>
                {runError && <p className="run-error" role="alert">{runError}</p>}
                <button className="analyze-button" type="button" disabled={!activeInput || !compiled || Boolean(packError) || databaseBusy || Boolean(busyCells.size)} onClick={requestRun}><span>{databaseBusy || busyCells.size ? "Validating references…" : "Analyze with SwiftIG"}</span><b>→</b></button>
                <p className="privacy-copy"><span>i</span> Query sequences, uploaded germlines, and AIRR results are processed in this browser; Swig does not transmit them. A remotely hosted alternative database is requested from its named provider only when selected.</p>
              </aside>
            </div>
          )}
        </main>
      )}

      {page === "results" && session && <ResultsPage session={session} onNewAnalysis={() => navigate("analyze")} />}
      {outputPrompt && activeInput && <div className="output-modal-backdrop" role="presentation"><section className="output-modal" role="dialog" aria-modal="true" aria-labelledby="output-dialog-title">
        <button className="output-modal-close" type="button" onClick={() => setOutputPrompt(false)} aria-label="Cancel output selection">×</button>
        <span className="output-direction">OUTPUT · SAVE RESULTS</span>
        <h2 id="output-dialog-title">Choose where Swig will <em>write the AIRR output.</em></h2>
        <p>The next system window is a <b>Save As</b> dialog. You are naming a new results file—this is not another sequence import.</p>
        <div className="output-flow"><div><span>Input already selected</span><strong>{activeInput.name}</strong></div><b>→</b><div className="destination"><span>New output file</span><strong>{outputName(activeInput.name)}</strong><small>Written incrementally during analysis</small></div></div>
        <div className="output-modal-actions"><button className="output-save-primary" type="button" onClick={() => void run("disk")}><span>Choose output file &amp; start</span><b>Save AIRR →</b></button><button type="button" onClick={() => void run("browser")}><span>Keep output in browser instead</span><small>Compressed local index; download after the run</small></button></div>
        <p className="output-safety"><span>i</span> Query sequences remain in this browser and are not transmitted by Swig.</p>
      </section></div>}
      <footer className="site-footer"><Brand /><p>Swig 0.8.2 · SwiftIG WebAssembly interface · research software · validate study-critical calls independently.</p><div><a href="https://github.com/MurrellGroup/swiftig" target="_blank" rel="noreferrer">Source ↗</a><a href="https://www.imgt.org/" target="_blank" rel="noreferrer">IMGT ↗</a><a href="https://docs.airr-community.org/" target="_blank" rel="noreferrer">AIRR ↗</a></div></footer>
    </div>
  );
}
