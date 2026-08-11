import {
  ChangeEvent,
  DragEvent,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import { AlignmentViewer } from "./alignment-view";
import { RepertoireDashboard } from "./repertoire-charts";
import {
  availableScopes,
  compileReferences,
  countFastaRecords,
  loadReferencePack,
  makeDemoFasta,
  type ReferencePack,
  type ReferenceSpecies,
  type ScopeKey,
  type SegmentKey,
} from "./reference-pack";
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
}

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
      <span className="local-badge"><i /> Local-only</span>
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
          <p className="eyebrow"><span>Browser-native immunogenetics</span></p>
          <h1>From rearranged reads<br />to <em>inspectable evidence.</em></h1>
          <p className="hero-lede">Call V(D)J genes across BCR and TCR loci, then interrogate every assignment down to the nucleotide or translated alignment—without uploading a base.</p>
          <div className="hero-actions">
            <button className="primary-cta" type="button" onClick={onStart}>Start an analysis <span>→</span></button>
            <button className="secondary-cta" type="button" onClick={onDemo}>Explore with demo data</button>
          </div>
          <div className="hero-proof">
            <span><b>FASTA</b> / FASTQ / AIRR · gzip</span>
            <span><b>7</b> IG + TR loci</span>
            <span><b>0</b> bases uploaded</span>
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
        <div className="section-heading"><p className="eyebrow"><span>A clear path through the data</span></p><h2>Three stages, no mystery state.</h2><p>Each run moves from configuration to measured progress to a durable local result index.</p></div>
        <div className="workflow-cards">
          <article><span>01</span><div className="workflow-icon upload-icon" /><h3>Bring sequences</h3><p>Upload or paste FASTA, FASTQ, or AIRR. Gzip is accepted, and large files remain out of the main UI thread.</p></article>
          <article><span>02</span><div className="workflow-icon engine-icon" /><h3>Annotate locally</h3><p>Choose species, BCR/TCR locus, and germlines. Swap V, D, J, or C independently. Use multiple isolated WASM workers across CPU cores.</p></article>
          <article><span>03</span><div className="workflow-icon result-icon" /><h3>Interrogate calls</h3><p>Filter a million records without rendering them all, then open any read for nucleotide or on-demand protein alignments.</p></article>
        </div>
      </section>

      <section className="scale-section">
        <div className="scale-copy"><p className="eyebrow"><span>Designed across scales</span></p><h2>One read gets a microscope.<br />A million reads get a pipeline.</h2><p>Small runs open directly into detailed evidence. Large runs are decompressed, parsed, annotated, indexed, and optionally written to disk under one bounded backpressure queue.</p></div>
        <div className="scale-grid">
          <article><strong>1–3</strong><span>Auto-open individual calls</span><small>Alignment-first review</small></article>
          <article><strong>10³</strong><span>Facets + paged table</span><small>Instant categorical filters</small></article>
          <article><strong>10⁶</strong><span>Streaming local batches</span><small>No million-row DOM or full-file buffer</small></article>
          <article><strong>{references ?? "—"}</strong><span>Reference sets</span><small>Plus arbitrary custom FASTA</small></article>
        </div>
      </section>

      <section className="landing-final"><div><span>Ready when your sequences are.</span><h2>Keep the data. See the evidence.</h2></div><button className="primary-cta light" type="button" onClick={onStart}>Open the analysis workspace <span>→</span></button></section>
    </main>
  );
}

function ReferenceSegment({ segment, count, override, optional, onFile, onClear }: {
  segment: SegmentKey;
  count: number;
  override?: ReferenceOverride;
  optional?: boolean;
  onFile: (segment: SegmentKey, file: File) => void;
  onClear: (segment: SegmentKey) => void;
}) {
  const input = useRef<HTMLInputElement>(null);
  return (
    <article className={`reference-segment ${override ? "custom" : ""}`}>
      <span className={`segment-symbol segment-${segment.toLowerCase()}`}>{segment}</span>
      <div><strong>{segment} germlines {optional && <small>optional</small>}</strong><span title={override?.name}>{override?.name ?? (count ? "Built-in IMGT set" : segment === "C" ? "Upload to enable constant calls" : "Not used for this locus")}</span><b>{count.toLocaleString()} alleles</b></div>
      <input ref={input} className="visually-hidden" type="file" accept=".fa,.fasta,.fna,.fas,.txt,.gz" onChange={(event) => {
        const file = event.target.files?.[0];
        if (file) onFile(segment, file);
        event.target.value = "";
      }} />
      {override ? <button className="segment-action remove" type="button" onClick={() => onClear(segment)} aria-label={`Remove custom ${segment} set`}>×</button> : <button className="segment-action" type="button" onClick={() => input.current?.click()}>{count ? "Replace" : "Upload"}</button>}
    </article>
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
  const regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"];
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
        <div className="uncertainty-heading"><div><span className="section-kicker">Call uncertainty</span><h3>Co-optimal and near-tied hits</h3></div><p>Comma-separated calls are exact co-optima. Alternate rows are retained only within SwiftIG’s uncertainty window; full alignment strings stay selected-call-only so million-read AIRR output remains bounded.</p></div>
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

      <section className="region-section">
        <div className="region-heading"><div><span className="section-kicker">Translated regions</span><h3>Framework + CDR sequence</h3></div><button type="button" onClick={() => void navigator.clipboard.writeText(`>${row.sequence_id}\n${row.sequence}\n`)}>Copy FASTA</button></div>
        <div className="region-strip">{regions.map((region) => <article className={region.startsWith("cdr") ? "cdr" : "fwr"} key={region}><span>{region.toUpperCase()}</span><code>{row[`${region}_aa`] || "—"}</code><small>{row[region]?.length ? `${row[region].length} nt` : "not called"}</small></article>)}</div>
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
  const [view, setView] = useState<"repertoire" | "sequences">(session.total <= 3 ? "sequences" : "repertoire");
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
        <div className="results-title"><div><p className="eyebrow"><span>Analysis complete · {session.seconds.toFixed(2)} s · {Math.round(session.total / Math.max(session.seconds, 0.001)).toLocaleString()} reads/s</span></p><h1>{session.total.toLocaleString()} analyzed rearrangements,<br /><em>ready to interrogate.</em></h1><p>{friendlySpecies(session.species)} · {LOCUS_LABELS[session.scope]} · {session.workers} WASM worker{session.workers === 1 ? "" : "s"} · {bytes(session.outputBytes)} AIRR{session.subsampleSize ? ` · exact random sample from ${session.inputTotal.toLocaleString()} input records (seed ${session.subsampleSeed})` : ""}</p></div><div className="results-actions"><button className="download-primary" type="button" onClick={() => void downloadAll()} disabled={downloading}>{downloading ? "Writing AIRR…" : session.streamedDirectly ? "Save another AIRR copy" : "Download AIRR TSV"}<span>↓</span></button><button type="button" onClick={onNewAnalysis}>New analysis</button>{downloadError && <small role="alert">{downloadError}</small>}</div></div>
        <div className="result-summary">
          <article><span>V + J assigned</span><strong>{session.summary.assigned.toLocaleString()}</strong><small>{percentage(session.summary.assigned, session.total)} of input</small></article>
          <article><span>Productive</span><strong>{session.summary.productive.toLocaleString()}</strong><small>{percentage(session.summary.productive, session.total)} of input</small></article>
          <article><span>CDR3 called</span><strong>{session.summary.withCdr3.toLocaleString()}</strong><small>{percentage(session.summary.withCdr3, session.total)} of input</small></article>
          <article><span>Loci observed</span><strong>{session.facets.loci.length}</strong><small>{session.facets.loci.map((item) => item.value).join(" · ") || "none"}</small></article>
        </div>
      </section>

      <nav className="results-view-tabs" aria-label="Results view"><button className={view === "repertoire" ? "active" : ""} type="button" onClick={() => setView("repertoire")}><span>Repertoire</span><small>Figures + composition</small></button><button className={view === "sequences" ? "active" : ""} type="button" onClick={() => setView("sequences")}><span>Sequences</span><small>Filter + inspect calls</small></button></nav>

      {view === "repertoire" ? <RepertoireDashboard store={session.store} loci={session.facets.loci} inputName={session.inputName} /> : <>

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
          <p className="index-note"><span>⌁</span> Exact gene and locus filters use local indexes. Substring filters scan candidates on demand and never leave this device.</p>
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
      </>}
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
  const [overrides, setOverrides] = useState<Partial<Record<SegmentKey, ReferenceOverride>>>({});
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
  const compiled = useMemo(() => species ? compileReferences(species, activeScope, Object.fromEntries(Object.entries(overrides).map(([key, value]) => [key, value?.text]))) : null, [activeScope, overrides, species]);
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

  async function acceptReferenceFile(segment: SegmentKey, file: File) {
    try {
      setRunError("");
      const text = await readUploadedText(file);
      if (text.trimStart()[0] !== ">") throw new Error(`${segment} references must be FASTA.`);
      const count = countFastaRecords(text);
      if (!count) throw new Error(`No ${segment} germline records were found.`);
      setOverrides((current) => ({ ...current, [segment]: { name: file.name, text, count, size: file.size } }));
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    }
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
    if (!compiled.counts.V || !compiled.counts.J) {
      setRunError("This reference selection requires both V and J germline records.");
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
          <section className="analysis-intro"><WorkflowStepper active={running ? 2 : 1} /><div><p className="eyebrow"><span>Analysis workspace</span></p><h1>{running ? "Calling rearrangements…" : "Configure one local run."}</h1><p>{running ? "SwiftIG is processing bounded batches and writing AIRR records into a local index." : "Add sequences, choose the biological search space, and make any germline substitutions explicit."}</p></div></section>

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
                  <header><span className="card-number">02</span><div><h2>Biological search space</h2><p>Choose a baseline, then replace any segment independently.</p></div>{Object.keys(overrides).length > 0 && <button className="reset-button" type="button" onClick={() => setOverrides({})}>Reset custom sets</button>}</header>
                  <div className="reference-selectors">
                    <label><span>Species / strain</span><select value={species?.name ?? ""} disabled={!pack} onChange={(event) => {
                      setSpeciesName(event.target.value);
                      const next = speciesList.find((item) => item.name === event.target.value);
                      const nextScopes = next ? availableScopes(next) : [];
                      if (!nextScopes.includes(activeScope)) setScope(nextScopes[0] ?? "BCR");
                      setOverrides({});
                    }}>{!pack && <option>Loading IMGT references…</option>}{speciesList.map((item) => <option value={item.name} key={item.name}>{friendlySpecies(item.name)}</option>)}</select></label>
                    <div className="receptor-selector"><span>Receptor</span><div><button className={receptor === "BCR" ? "active" : ""} type="button" disabled={!hasBcr} onClick={() => setScope("BCR")}>BCR <small>IG</small></button><button className={receptor === "TCR" ? "active" : ""} type="button" disabled={!hasTcr} onClick={() => setScope("TCR")}>TCR <small>TR</small></button></div></div>
                    <label><span>Chain / locus</span><select value={activeScope} onChange={(event) => setScope(event.target.value as ScopeKey)}>{receptorScopes.map((value) => <option value={value} key={value}>{LOCUS_LABELS[value]}</option>)}</select></label>
                  </div>
                  {packError && <p className="inline-error" role="alert">{packError}</p>}
                  <div className="reference-segments">{SEGMENTS.map((segment) => <ReferenceSegment key={segment} segment={segment} count={compiled?.counts[segment] ?? 0} override={overrides[segment]} optional={segment === "D" || segment === "C"} onFile={(key, file) => void acceptReferenceFile(key, file)} onClear={(key) => setOverrides((current) => { const next = { ...current }; delete next[key]; return next; })} />)}</div>
                  <p className="reference-footnote"><span>i</span><b>{compiled?.loci.join(" + ") || "No locus"}</b> · IMGT/GENE-DB {pack?.release ?? "…"}. Each upload replaces only that segment. C references enable constant-region and isotype calls when ≥30 nt align.</p>
                </section>

                <section className="analysis-card settings-card">
                  <header><span className="card-number">03</span><div><h2>Analysis behavior</h2><p>Defaults are tuned for rearranged repertoire sequences.</p></div><button className="reset-button" type="button" onClick={() => setShowAdvanced((value) => !value)}>{showAdvanced ? "Hide" : "Show"} controls</button></header>
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
                <span className="section-kicker">Run manifest</span><h2>{activeInput ? subsampleEnabled ? `${Math.floor(subsampleSize).toLocaleString()}-read sample` : activeInput.count === null ? "Large file" : `${activeInput.count?.toLocaleString()} sequences` : "Awaiting data"}</h2><p>{activeInput?.name ?? "Add an input to unlock analysis."}</p>
                <dl><div><dt>Reference</dt><dd>{species ? friendlySpecies(species.name) : "Loading…"}</dd></div><div><dt>Search</dt><dd>{LOCUS_LABELS[activeScope]}</dd></div><div><dt>V / D / J / C</dt><dd>{compiled ? `${compiled.counts.V} / ${compiled.counts.D} / ${compiled.counts.J} / ${compiled.counts.C}` : "—"}</dd></div><div><dt>Custom</dt><dd>{Object.keys(overrides).length ? Object.keys(overrides).join(" + ") : "None"}</dd></div><div><dt>Compute</dt><dd>{workerCount} workers · bounded queue</dd></div><div><dt>Storage</dt><dd>{outputStorage === "auto" ? "Adaptive streaming" : outputStorage === "disk" ? "Direct AIRR file" : "Compressed local index"}</dd></div></dl>
                {runError && <p className="run-error" role="alert">{runError}</p>}
                <button className="analyze-button" type="button" disabled={!activeInput || !compiled || Boolean(packError)} onClick={requestRun}><span>Analyze with SwiftIG</span><b>→</b></button>
                <p className="privacy-copy"><span>⌁</span> Query data, custom germlines, and AIRR output remain on this device. For large runs, Swig clearly asks where to save a new output file before analysis starts.</p>
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
        <p className="output-safety"><span>i</span> Your sequence data still never leave this device.</p>
      </section></div>}
      <footer className="site-footer"><Brand /><p>Swig 0.4 · parallel streaming WebAssembly · research software · validate study-critical calls independently.</p><div><a href="https://github.com/MurrellGroup/swiftig" target="_blank" rel="noreferrer">Source ↗</a><a href="https://www.imgt.org/" target="_blank" rel="noreferrer">IMGT ↗</a><a href="https://docs.airr-community.org/" target="_blank" rel="noreferrer">AIRR ↗</a></div></footer>
    </div>
  );
}
