import { ChangeEvent, DragEvent, useEffect, useMemo, useRef, useState } from "react";

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
import { runSwiftIg } from "./swiftig-runtime";

type InputFormat = "FASTA" | "FASTQ" | "AIRR TSV";

interface InputData {
  name: string;
  text: string;
  count: number;
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

interface ResultRow {
  [key: string]: string;
}

const SEGMENTS: SegmentKey[] = ["V", "D", "J"];
const RESULT_COLUMNS = [
  "sequence_id",
  "locus",
  "v_call",
  "d_call",
  "j_call",
  "productive",
  "junction_aa",
  "rev_comp",
];

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
  return `${(value / (1024 * 1024)).toFixed(1)} MB`;
}

async function readUploadedText(file: File): Promise<string> {
  if (!file.name.toLowerCase().endsWith(".gz")) return file.text();
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot decompress gzip files. Decompress the file first.");
  }
  const decompressed = file.stream().pipeThrough(new DecompressionStream("gzip"));
  return new Response(decompressed).text();
}

function inspectInput(name: string, text: string, size = text.length): InputData {
  const first = text.trimStart()[0];
  if (first === ">") {
    const count = text.split(/\r?\n/).filter((line) => line.startsWith(">")).length;
    if (!count) throw new Error("No FASTA records were found.");
    return { name, text, count, format: "FASTA", formatCode: 1, size };
  }
  if (first === "@") {
    const lines = text.split(/\r?\n/);
    let count = 0;
    for (let index = 0; index < lines.length; index += 4) {
      if (lines[index]?.startsWith("@")) count += 1;
    }
    if (!count) throw new Error("No FASTQ records were found.");
    return { name, text, count, format: "FASTQ", formatCode: 2, size };
  }
  const header = text.split(/\r?\n/, 1)[0]?.split("\t") ?? [];
  if (header.includes("sequence")) {
    const count = Math.max(0, text.trimEnd().split(/\r?\n/).length - 1);
    if (!count) throw new Error("The AIRR table has no data rows.");
    return { name, text, count, format: "AIRR TSV", formatCode: 3, size };
  }
  throw new Error("Expected FASTA, FASTQ, or an AIRR TSV with a sequence column.");
}

function parseResults(tsv: string): { headers: string[]; rows: ResultRow[] } {
  const lines = tsv.trimEnd().split("\n");
  const headers = lines.shift()?.split("\t") ?? [];
  const rows = lines.map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
  });
  return { headers, rows };
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

function friendlySpecies(value: string): string {
  return value.replaceAll("_", " · ");
}

function SegmentCard({
  segment,
  count,
  override,
  optional,
  onFile,
  onClear,
}: {
  segment: SegmentKey;
  count: number;
  override?: ReferenceOverride;
  optional?: boolean;
  onFile: (segment: SegmentKey, file: File) => void;
  onClear: (segment: SegmentKey) => void;
}) {
  const inputRef = useRef<HTMLInputElement>(null);
  return (
    <article className={`segment-card ${override ? "segment-card--custom" : ""}`}>
      <div className="segment-letter" aria-hidden="true">{segment}</div>
      <div className="segment-copy">
        <div className="segment-heading">
          <strong>{segment} genes</strong>
          {optional && <span className="optional-tag">optional</span>}
        </div>
        <span className="segment-source" title={override?.name}>
          {override ? override.name : count ? "Built-in reference" : "Not required for this chain"}
        </span>
        <span className="segment-count">
          {count.toLocaleString()} {count === 1 ? "allele" : "alleles"}
        </span>
      </div>
      <input
        ref={inputRef}
        className="visually-hidden"
        type="file"
        accept=".fa,.fasta,.fna,.fas,.txt,.gz"
        aria-label={`Upload custom ${segment} germline FASTA`}
        onChange={(event) => {
          const file = event.target.files?.[0];
          if (file) onFile(segment, file);
          event.target.value = "";
        }}
      />
      {override ? (
        <button className="icon-button" type="button" onClick={() => onClear(segment)} aria-label={`Remove custom ${segment} database`}>
          ×
        </button>
      ) : (
        <button className="tiny-button" type="button" onClick={() => inputRef.current?.click()}>
          Swap
        </button>
      )}
    </article>
  );
}

export default function SwigApp() {
  const [pack, setPack] = useState<ReferencePack | null>(null);
  const [packError, setPackError] = useState("");
  const [speciesName, setSpeciesName] = useState("Homo sapiens");
  const [scope, setScope] = useState<ScopeKey>("BCR");
  const [input, setInput] = useState<InputData | null>(null);
  const [inputError, setInputError] = useState("");
  const [overrides, setOverrides] = useState<Partial<Record<SegmentKey, ReferenceOverride>>>({});
  const [minimumIdentity, setMinimumIdentity] = useState(0.6);
  const [strand, setStrand] = useState<0 | 1 | 2>(0);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState({ stage: "Ready", value: 0 });
  const [resultTsv, setResultTsv] = useState("");
  const [resultSeconds, setResultSeconds] = useState(0);
  const [runError, setRunError] = useState("");
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    loadReferencePack().then(setPack).catch((error) => {
      setPackError(error instanceof Error ? error.message : String(error));
    });
  }, []);

  const speciesList = useMemo(() => favoriteSpecies(pack?.species ?? []), [pack]);
  const species = useMemo(
    () => speciesList.find((candidate) => candidate.name === speciesName) ?? speciesList[0],
    [speciesList, speciesName],
  );
  const scopes = useMemo(() => species ? availableScopes(species) : [], [species]);
  const activeScope = scopes.includes(scope) ? scope : scopes[0] ?? "BCR";

  const compiled = useMemo(() => {
    if (!species) return null;
    return compileReferences(
      species,
      activeScope,
      Object.fromEntries(Object.entries(overrides).map(([key, value]) => [key, value?.text])),
    );
  }, [activeScope, overrides, species]);

  const parsedResults = useMemo(() => resultTsv ? parseResults(resultTsv) : null, [resultTsv]);
  const previewRows = parsedResults?.rows.slice(0, 100) ?? [];
  const resultSummary = useMemo(() => {
    const rows = parsedResults?.rows ?? [];
    const called = rows.filter((row) => row.v_call && row.j_call).length;
    const productive = rows.filter((row) => row.productive === "T").length;
    const loci = new Set(rows.map((row) => row.locus).filter(Boolean)).size;
    return { total: rows.length, called, productive, loci };
  }, [parsedResults]);

  const receptor = activeScope.startsWith("IG") || activeScope === "BCR" ? "BCR" : "TCR";
  const receptorScopes = scopes.filter((value) => receptor === "BCR"
    ? value === "BCR" || value.startsWith("IG")
    : value === "TCR" || value.startsWith("TR"));

  async function acceptInputFile(file: File) {
    try {
      setInputError("");
      setRunError("");
      setResultTsv("");
      const text = await readUploadedText(file);
      setInput(inspectInput(file.name, text, file.size));
    } catch (error) {
      setInput(null);
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
      setOverrides((current) => ({
        ...current,
        [segment]: { name: file.name, text, count, size: file.size },
      }));
      setResultTsv("");
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    }
  }

  function chooseDemo() {
    if (!species) return;
    const text = makeDemoFasta(species, activeScope);
    if (!text.trim()) {
      setInputError("No complete demo locus is available for this reference selection.");
      return;
    }
    setInputError("");
    setResultTsv("");
    setInput(inspectInput("swig-demo.fasta", text));
  }

  async function run() {
    if (!input || !compiled) return;
    if (!compiled.counts.V || !compiled.counts.J) {
      setRunError("This reference selection requires both V and J germline records.");
      return;
    }
    setRunning(true);
    setRunError("");
    setResultTsv("");
    setProgress({ stage: "Preparing annotation", value: 0.02 });
    const started = performance.now();
    try {
      const result = await runSwiftIg({
        query: input.text,
        format: input.formatCode,
        references: compiled,
        minimumIdentity,
        strand,
        onProgress: (stage, value) => setProgress({ stage, value }),
      });
      setResultSeconds((performance.now() - started) / 1000);
      setResultTsv(result.tsv);
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      setRunning(false);
    }
  }

  function downloadResults() {
    const blob = new Blob([resultTsv], { type: "text/tab-separated-values;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = `${input?.name.replace(/\.[^.]+$/, "") || "swig"}.airr.tsv`;
    anchor.click();
    URL.revokeObjectURL(url);
  }

  const hasBcr = scopes.includes("BCR");
  const hasTcr = scopes.includes("TCR");

  return (
    <div className="site-shell">
      <header className="topbar">
        <a className="brand" href="#top" aria-label="Swig home">
          <span className="brand-mark" aria-hidden="true"><i /><i /><i /></span>
          <span className="brand-name">SWIG</span>
          <span className="brand-subtitle">SwiftIG · browser edition</span>
        </a>
        <nav className="nav-links" aria-label="Primary navigation">
          <a href="#annotate">Annotate</a>
          <a href="#references">References</a>
          <a href="#about">About</a>
        </nav>
        <div className="privacy-pill"><span /> Local · WASM</div>
      </header>

      <main id="top">
        <section className="hero">
          <div className="hero-copy">
            <p className="eyebrow"><span>01</span> Private immunogenetics</p>
            <h1>V(D)J calls at<br /><em>browser speed.</em></h1>
            <p className="hero-lede">
              Annotate BCR and TCR reads with SwiftIG—without sending a single base off your machine.
              FASTA, FASTQ, and AIRR in; AIRR out.
            </p>
            <div className="hero-actions">
              <a className="primary-link" href="#annotate">Start annotating <span>↓</span></a>
              <button className="text-button" type="button" onClick={chooseDemo} disabled={!species}>Load a demo</button>
            </div>
          </div>
          <div className="hero-visual" aria-label="Swig capabilities">
            <div className="orbit orbit-one" />
            <div className="orbit orbit-two" />
            <div className="helix" aria-hidden="true">
              {Array.from({ length: 7 }, (_, index) => <i key={index} />)}
            </div>
            <div className="metric metric-top"><strong>{pack?.species.length ?? "—"}</strong><span>animal reference sets</span></div>
            <div className="metric metric-bottom"><strong>7</strong><span>IG + TR loci</span></div>
            <div className="metric metric-side"><strong>0</strong><span>bases uploaded</span></div>
          </div>
        </section>

        <section className="trust-strip" aria-label="Supported capabilities">
          <span><b>Formats</b> FASTA · FASTQ · AIRR TSV</span>
          <span><b>Receptors</b> IGH · IGK · IGL · TRA · TRB · TRD · TRG</span>
          <span><b>Output</b> AIRR Rearrangement</span>
        </section>

        <section className="workspace" id="annotate">
          <div className="workspace-heading">
            <div>
              <p className="eyebrow"><span>02</span> Annotation workspace</p>
              <h2>Three choices. One local run.</h2>
            </div>
            <p>Your data and custom references stay inside this browser tab.</p>
          </div>

          <div className="workspace-grid">
            <div className="workflow-column">
              <article className="workflow-card data-card">
                <div className="card-index">1</div>
                <div className="card-title-row">
                  <div>
                    <h3>Add sequences</h3>
                    <p>Assemblies, reads, or an AIRR export with a <code>sequence</code> column.</p>
                  </div>
                  {input && <button className="text-button subtle" type="button" onClick={() => setInput(null)}>Clear</button>}
                </div>
                <input
                  ref={inputRef}
                  className="visually-hidden"
                  type="file"
                  accept=".fa,.fasta,.fna,.fas,.fq,.fastq,.tsv,.txt,.gz"
                  onChange={(event: ChangeEvent<HTMLInputElement>) => {
                    const file = event.target.files?.[0];
                    if (file) void acceptInputFile(file);
                    event.target.value = "";
                  }}
                />
                {input ? (
                  <div className="file-loaded">
                    <div className="file-icon" aria-hidden="true">↳</div>
                    <div className="file-details">
                      <strong title={input.name}>{input.name}</strong>
                      <span>{input.count.toLocaleString()} sequences · {input.format} · {bytes(input.size)}</span>
                    </div>
                    <span className="checkmark">Ready</span>
                  </div>
                ) : (
                  <button
                    className="dropzone"
                    type="button"
                    onClick={() => inputRef.current?.click()}
                    onDragOver={(event: DragEvent) => event.preventDefault()}
                    onDrop={(event: DragEvent) => {
                      event.preventDefault();
                      const file = event.dataTransfer.files?.[0];
                      if (file) void acceptInputFile(file);
                    }}
                  >
                    <span className="drop-icon" aria-hidden="true">＋</span>
                    <strong>Drop sequence data here</strong>
                    <span>or choose a file · FASTA, FASTQ, AIRR TSV · gzip is fine</span>
                  </button>
                )}
                {inputError && <p className="inline-error" role="alert">{inputError}</p>}
                {!input && <button className="demo-link" type="button" onClick={chooseDemo} disabled={!species}>No file handy? Use a generated {activeScope} demo →</button>}
              </article>

              <article className="workflow-card" id="references">
                <div className="card-index">2</div>
                <div className="card-title-row">
                  <div>
                    <h3>Choose references</h3>
                    <p>Start broad, then replace any one segment independently.</p>
                  </div>
                  {Object.keys(overrides).length > 0 && (
                    <button className="text-button subtle" type="button" onClick={() => setOverrides({})}>Reset swaps</button>
                  )}
                </div>

                <div className="reference-controls">
                  <label className="field-label">
                    <span>Species / strain</span>
                    <select
                      value={species?.name ?? ""}
                      disabled={!pack}
                      onChange={(event) => {
                        setSpeciesName(event.target.value);
                        const nextSpecies = speciesList.find((item) => item.name === event.target.value);
                        const nextScopes = nextSpecies ? availableScopes(nextSpecies) : [];
                        if (!nextScopes.includes(activeScope)) setScope(nextScopes[0] ?? "BCR");
                        setOverrides({});
                        setResultTsv("");
                      }}
                    >
                      {!pack && <option>Loading IMGT reference library…</option>}
                      {speciesList.map((item) => <option key={item.name} value={item.name}>{friendlySpecies(item.name)}</option>)}
                    </select>
                  </label>

                  <div className="field-label">
                    <span>Receptor system</span>
                    <div className="segmented-control" role="group" aria-label="Receptor system">
                      <button
                        type="button"
                        className={receptor === "BCR" ? "active" : ""}
                        disabled={!hasBcr}
                        onClick={() => { setScope("BCR"); setResultTsv(""); }}
                      >BCR <small>IG</small></button>
                      <button
                        type="button"
                        className={receptor === "TCR" ? "active" : ""}
                        disabled={!hasTcr}
                        onClick={() => { setScope("TCR"); setResultTsv(""); }}
                      >TCR <small>TR</small></button>
                    </div>
                  </div>

                  <label className="field-label">
                    <span>Chain / locus</span>
                    <select value={activeScope} onChange={(event) => { setScope(event.target.value as ScopeKey); setResultTsv(""); }}>
                      {receptorScopes.map((value) => <option key={value} value={value}>{LOCUS_LABELS[value]}</option>)}
                    </select>
                  </label>
                </div>

                {packError && <p className="inline-error" role="alert">{packError}</p>}
                <div className="segment-grid">
                  {SEGMENTS.map((segment) => (
                    <SegmentCard
                      key={segment}
                      segment={segment}
                      count={compiled?.counts[segment] ?? 0}
                      override={overrides[segment]}
                      optional={segment === "D"}
                      onFile={(key, file) => void acceptReferenceFile(key, file)}
                      onClear={(key) => {
                        setOverrides((current) => {
                          const next = { ...current };
                          delete next[key];
                          return next;
                        });
                        setResultTsv("");
                      }}
                    />
                  ))}
                </div>
                <p className="reference-note">
                  <span>i</span> {compiled?.loci.join(" + ") || "No locus selected"} · IMGT/GENE-DB release {pack?.release ?? "…"}.
                  Uploaded FASTA replaces only its matching segment.
                </p>
              </article>

              <article className="workflow-card settings-card">
                <div className="card-index">3</div>
                <div className="card-title-row">
                  <div>
                    <h3>Run settings</h3>
                    <p>Sensible repertoire defaults, with explicit strand and identity control.</p>
                  </div>
                  <button className="text-button subtle" type="button" onClick={() => setShowAdvanced((value) => !value)}>
                    {showAdvanced ? "Hide" : "Advanced"} {showAdvanced ? "↑" : "↓"}
                  </button>
                </div>
                <div className="settings-summary">
                  <span><b>Mode</b> nucleotide V(D)J</span>
                  <span><b>Domain</b> IMGT</span>
                  <span><b>Output</b> AIRR TSV</span>
                </div>
                {showAdvanced && (
                  <div className="advanced-panel">
                    <label className="field-label">
                      <span>Search strand</span>
                      <select value={strand} onChange={(event) => setStrand(Number(event.target.value) as 0 | 1 | 2)}>
                        <option value={0}>Both</option>
                        <option value={1}>Plus only</option>
                        <option value={2}>Minus only</option>
                      </select>
                    </label>
                    <label className="range-field">
                      <span><b>Minimum identity</b><output>{Math.round(minimumIdentity * 100)}%</output></span>
                      <input
                        type="range"
                        min="0.45"
                        max="0.9"
                        step="0.01"
                        value={minimumIdentity}
                        onChange={(event) => setMinimumIdentity(Number(event.target.value))}
                      />
                    </label>
                    <div className="constant-upload">
                      <div><b>Constant-region database</b><span>Optional C-gene calls</span></div>
                      <label className="tiny-button">
                        {overrides.C ? overrides.C.name : "Upload C FASTA"}
                        <input
                          className="visually-hidden"
                          type="file"
                          accept=".fa,.fasta,.fna,.txt,.gz"
                          onChange={(event) => {
                            const file = event.target.files?.[0];
                            if (file) void acceptReferenceFile("C", file);
                            event.target.value = "";
                          }}
                        />
                      </label>
                    </div>
                  </div>
                )}
              </article>
            </div>

            <aside className="run-panel">
              <div className="run-panel-top">
                <span className="run-kicker">Run manifest</span>
                <span className={`engine-state ${running ? "engine-state--running" : ""}`}><i /> {running ? "Working" : "WASM ready"}</span>
              </div>
              <h3>{input ? `${input.count.toLocaleString()} sequences` : "Awaiting data"}</h3>
              <dl className="manifest-list">
                <div><dt>Reference</dt><dd>{species ? friendlySpecies(species.name) : "Loading…"}</dd></div>
                <div><dt>Chains</dt><dd>{LOCUS_LABELS[activeScope]}</dd></div>
                <div><dt>Segments</dt><dd>{compiled ? `${compiled.counts.V} V · ${compiled.counts.D} D · ${compiled.counts.J} J` : "—"}</dd></div>
                <div><dt>Custom</dt><dd>{Object.keys(overrides).length ? Object.keys(overrides).join(" + ") : "None"}</dd></div>
              </dl>
              {running && (
                <div className="progress-wrap" aria-live="polite">
                  <div className="progress-label"><span>{progress.stage}</span><b>{Math.round(progress.value * 100)}%</b></div>
                  <div className="progress-track"><i style={{ width: `${progress.value * 100}%` }} /></div>
                </div>
              )}
              {runError && <p className="run-error" role="alert">{runError}</p>}
              <button className="run-button" type="button" disabled={!input || !compiled || running || Boolean(packError)} onClick={() => void run()}>
                <span>{running ? "Annotating…" : "Run SwiftIG"}</span>
                <b aria-hidden="true">→</b>
              </button>
              <p className="local-note"><span aria-hidden="true">⌁</span> Runs locally in an isolated worker. Closing the tab clears sequence data.</p>
            </aside>
          </div>
        </section>

        {parsedResults && (
          <section className="results-section" aria-live="polite">
            <div className="results-heading">
              <div>
                <p className="eyebrow"><span>03</span> Results</p>
                <h2>{resultSummary.total.toLocaleString()} AIRR records, ready.</h2>
                <p>Completed locally in {resultSeconds.toFixed(2)} seconds.</p>
              </div>
              <button className="download-button" type="button" onClick={downloadResults}>Download AIRR TSV <span>↓</span></button>
            </div>
            <div className="result-metrics">
              <article><span>V + J assigned</span><strong>{resultSummary.called.toLocaleString()}</strong><small>{resultSummary.total ? Math.round(resultSummary.called / resultSummary.total * 100) : 0}% of input</small></article>
              <article><span>Productive</span><strong>{resultSummary.productive.toLocaleString()}</strong><small>{resultSummary.total ? Math.round(resultSummary.productive / resultSummary.total * 100) : 0}% of input</small></article>
              <article><span>Loci observed</span><strong>{resultSummary.loci}</strong><small>from selected references</small></article>
              <article><span>Throughput</span><strong>{Math.round(resultSummary.total / Math.max(resultSeconds, 0.001)).toLocaleString()}</strong><small>sequences / second</small></article>
            </div>
            <div className="table-shell">
              <div className="table-caption"><b>AIRR preview</b><span>Showing {previewRows.length.toLocaleString()} of {resultSummary.total.toLocaleString()} rows · full schema in download</span></div>
              <div className="table-scroll">
                <table>
                  <thead><tr>{RESULT_COLUMNS.map((column) => <th key={column}>{column}</th>)}</tr></thead>
                  <tbody>
                    {previewRows.map((row, rowIndex) => (
                      <tr key={`${row.sequence_id}-${rowIndex}`}>
                        {RESULT_COLUMNS.map((column) => (
                          <td key={column} className={column.endsWith("_call") ? "gene-cell" : ""} title={row[column]}>
                            {row[column] || <span className="null-value">—</span>}
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </section>
        )}

        <section className="about-section" id="about">
          <div>
            <p className="eyebrow"><span>04</span> What is running?</p>
            <h2>SwiftIG’s C++20 annotation core, compiled for the web.</h2>
          </div>
          <div className="about-grid">
            <article><span>01</span><h3>Private by construction</h3><p>Parsing, reference indexing, alignment, traceback, and AIRR serialization all execute in-browser.</p></article>
            <article><span>02</span><h3>BCR and TCR aware</h3><p>Reference identifiers determine IGH, IGK, IGL, TRA, TRB, TRD, or TRG. D search is automatically limited to heavy-like loci.</p></article>
            <article><span>03</span><h3>Easy to perturb</h3><p>Keep the standard V and D sets, swap one J FASTA, and rerun—no database build step and no BLAST indexes.</p></article>
          </div>
          <p className="science-note">
            SwiftIG 0.2 is research-grade. Validate calls against IgBLAST for the organism, reference set,
            protocol, and mutation regime that matter to your study. Built-in references are a compact,
            segment-organized derivative of IMGT/GENE-DB release {pack?.release ?? "202632-7"}; IMGT terms and attribution apply.
          </p>
        </section>
      </main>

      <footer>
        <a className="brand footer-brand" href="#top"><span className="brand-name">SWIG</span></a>
        <p>SwiftIG in WebAssembly · local-first V(D)J annotation</p>
        <div><a href="https://github.com/MurrellGroup/swiftig" target="_blank" rel="noreferrer">Source ↗</a><a href="https://www.imgt.org/" target="_blank" rel="noreferrer">IMGT ↗</a><a href="https://docs.airr-community.org/" target="_blank" rel="noreferrer">AIRR ↗</a></div>
      </footer>
    </div>
  );
}
