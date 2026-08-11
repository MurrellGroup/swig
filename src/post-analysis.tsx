import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type ChangeEvent,
  type ReactNode,
} from "react";

import { runCodonAwareKalign, runFastTree, runKalign } from "./biowasm-runtime";
import {
  runChmmairra,
  writeChmmairraTsv,
  type ChmmDashboard,
  type ChmmSegment,
} from "./chmmairra-runtime";
import { GERMLINE_OUTGROUP, lineageInputFasta, quickAirrAlignment } from "./lineage-alignment";
import { layoutTree, parseNewick, rootOnOutgroup, serializeNewick } from "./phylogeny";
import {
  parseFasta,
  prepareReferenceMsa,
  type AmbiguityPolicy,
  type CallResolution,
  type DedupKey,
  type LineageSummary,
  type QueryHit,
  type QueryMetric,
  type QueryTarget,
} from "./post-analysis-core";
import {
  PostAnalysisRuntime,
  type DedupDashboard,
  type LineageDashboard,
} from "./post-analysis-runtime";
import type { CompiledReferences, ScopeKey } from "./reference-pack";
import type { AirrDetailRow, AirrIndexRecord, AirrResultStore, FacetValue } from "./result-store";

interface Props {
  store: AirrResultStore;
  references: CompiledReferences;
  scope: ScopeKey;
  loci: FacetValue[];
  inputName: string;
  workers: number;
  onInspect: (ordinal: number) => void;
}

interface SaveFileHandle {
  createWritable: () => Promise<{ write: (value: string | Blob | Uint8Array) => Promise<void>; close: () => Promise<void>; abort?: () => Promise<void> }>;
}

interface ChartDatum {
  label: string;
  value: number;
}

function baseName(name: string): string {
  return name.replace(/(\.airr)?\.(tsv|csv|txt|fa|fasta|fastq)(\.gz)?$/i, "") || "swig";
}

function downloadBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1_000);
}

async function saveStream(
  name: string,
  description: string,
  extension: string,
  produce: (writer: { write: (value: string | Blob | Uint8Array) => Promise<void> }) => Promise<void>,
) {
  const picker = (window as Window & { showSaveFilePicker?: (options: unknown) => Promise<SaveFileHandle> }).showSaveFilePicker;
  if (picker) {
    const handle = await picker.call(window, {
      suggestedName: name,
      types: [{ description, accept: { "text/plain": [extension] } }],
    });
    const writable = await handle.createWritable();
    try {
      await produce(writable);
      await writable.close();
    } catch (error) {
      await writable.abort?.();
      throw error;
    }
    return;
  }
  const parts: BlobPart[] = [];
  await produce({ write: async (value) => { parts.push(value instanceof Uint8Array ? value.slice().buffer : value); } });
  downloadBlob(new Blob(parts, { type: "text/plain;charset=utf-8" }), name);
}

function downloadText(value: string, name: string, type = "text/plain;charset=utf-8") {
  downloadBlob(new Blob([value], { type }), name);
}

function saveSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  downloadText(new XMLSerializer().serializeToString(clone), name, "image/svg+xml;charset=utf-8");
}

function formatNumber(value: number): string {
  return value >= 1_000_000 ? `${(value / 1_000_000).toFixed(1)}m` : value >= 1_000 ? `${(value / 1_000).toFixed(1)}k` : value.toLocaleString();
}

function BarChart({ title, subtitle, data, color, name, controls }: {
  title: string;
  subtitle: string;
  data: ChartDatum[];
  color: string;
  name: string;
  controls?: ReactNode;
}) {
  const svgRef = useRef<SVGSVGElement>(null);
  const visible = data.slice(0, 24);
  const maximum = Math.max(1, ...visible.map((item) => item.value));
  const height = Math.max(210, visible.length * 25 + 52);
  return <article className="post-chart-card">
    <header><div><span className="section-kicker">Repertoire summary</span><h3>{title}</h3><p>{subtitle}</p></div><div className="post-chart-actions">{controls}<button type="button" onClick={() => saveSvg(svgRef.current, name)}>SVG ↓</button></div></header>
    <div className="post-chart-scroll"><svg ref={svgRef} role="img" aria-label={title} viewBox={`0 0 760 ${height}`}>
      <rect width="760" height={height} fill="#fffdf8" />
      {visible.map((item, index) => {
        const y = 34 + index * 25;
        const width = item.value / maximum * 450;
        return <g key={`${item.label}-${index}`}>
          <text x="8" y={y + 14} fontFamily="Inter,Arial,sans-serif" fontSize="11" fill="#31413c">{item.label.length > 25 ? `${item.label.slice(0, 23)}…` : item.label}</text>
          <rect x="190" y={y} width="450" height="16" rx="3" fill="#e7ebe5" />
          <rect x="190" y={y} width={Math.max(1, width)} height="16" rx="3" fill={color} />
          <text x={Math.min(704, 200 + width)} y={y + 13} fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700" fill="#16231f">{formatNumber(item.value)}</text>
        </g>;
      })}
      {!visible.length && <text x="380" y="110" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="13" fill="#6d7975">No values</text>}
    </svg></div>
  </article>;
}

function TreeView({ newick, name }: { newick: string; name: string }) {
  const svgRef = useRef<SVGSVGElement>(null);
  const tree = useMemo(() => layoutTree(parseNewick(newick), 980, 24), [newick]);
  return <section className="tree-viewer">
    <header><div><span className="section-kicker">Rooted phylogeny</span><h4>{tree.leaves.toLocaleString()} tips</h4></div><div><button type="button" onClick={() => saveSvg(svgRef.current, `${name}.svg`)}>Download SVG</button><button type="button" onClick={() => downloadText(`${newick.replace(/;?$/, ";")}\n`, `${name}.nwk`)}>Newick</button></div></header>
    <div className="tree-scroll"><svg ref={svgRef} viewBox={`0 0 ${tree.width} ${tree.height}`} role="img" aria-label="N-masked germline rooted lineage tree">
      <rect width={tree.width} height={tree.height} fill="#fbfaf5" />
      {tree.edges.map((edge, index) => <g key={index} stroke="#49605a" strokeWidth="1.2" fill="none"><path d={`M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}`} /></g>)}
      {tree.nodes.map((node, index) => <g key={index}><circle cx={node.x} cy={node.y} r={node.children.length ? 2.2 : 2.8} fill={node.name === GERMLINE_OUTGROUP ? "#d49a19" : "#08796f"} />{!node.children.length && <text x={node.x + 7} y={node.y + 4} fontFamily="ui-monospace,monospace" fontSize="10" fontWeight={node.name === GERMLINE_OUTGROUP ? "700" : "500"} fill="#1a2925">{node.name}</text>}</g>)}
    </svg></div>
    <p className="scientific-note"><span>i</span> FastTree is inferred from nucleotide alignment. The reconstructed germline has N at non-templated junction positions and is placed as the outgroup; its incident edge defines the displayed root.</p>
  </section>;
}

function translateAligned(sequence: string): string {
  const codons: Record<string, string> = {
    TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S", TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
    CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P", CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
    ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T", AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
    GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A", GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
  };
  let result = "";
  for (let index = 0; index < sequence.length; index += 3) {
    const codon = sequence.slice(index, index + 3);
    result += codon === "---" ? "-" : codon.includes("-") || codon.includes("N") || codon.length < 3 ? "X" : codons[codon] ?? "X";
  }
  return result;
}

function AlignmentPreview({ fasta, mode }: { fasta: string; mode: "nt" | "aa" }) {
  const records = parseFasta(fasta, true);
  return <div className="lineage-alignment-preview">
    <div className="alignment-ruler"><span>Name</span><span>{mode === "nt" ? "Nucleotide alignment" : "Codon translation"} · showing {Math.min(80, records.length)} of {records.length}</span></div>
    {records.slice(0, 80).map((record) => <div className={record.name === GERMLINE_OUTGROUP ? "germline-row" : ""} key={record.name}><strong title={record.name}>{record.name}</strong><code>{mode === "nt" ? record.sequence : translateAligned(record.sequence)}</code></div>)}
  </div>;
}

function parseQueries(text: string): string[] {
  if (text.trimStart().startsWith(">")) return parseFasta(text).map((record) => record.sequence).filter(Boolean);
  return text.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith("#")).map((line) => line.replace(/\s/g, ""));
}

export function PostAnalysisWorkbench({ store, references, scope, loci, inputName, workers, onInspect }: Props) {
  const runtime = useMemo(() => new PostAnalysisRuntime(store), [store]);
  useEffect(() => () => runtime.terminate(), [runtime]);
  const [busy, setBusy] = useState("");
  const [progress, setProgress] = useState({ processed: 0, total: store.count });
  const [error, setError] = useState("");

  const [dedupKey, setDedupKey] = useState<DedupKey>("sequence");
  const [dedup, setDedup] = useState<DedupDashboard | null>(null);

  const [identity, setIdentity] = useState(0.85);
  const [resolution, setResolution] = useState<CallResolution>("gene");
  const [ambiguity, setAmbiguity] = useState<AmbiguityPolicy>("overlap");
  const [productiveOnly, setProductiveOnly] = useState(true);
  const [candidateCap, setCandidateCap] = useState(50_000);
  const [useDedup, setUseDedup] = useState(true);
  const [lineages, setLineages] = useState<LineageDashboard | null>(null);
  const [geneMetric, setGeneMetric] = useState<"abundance" | "lineages">("abundance");
  const [chartColor, setChartColor] = useState("#08796f");
  const [topGenes, setTopGenes] = useState(15);

  const [selectedLineage, setSelectedLineage] = useState<LineageSummary | null>(null);
  const [lineageRows, setLineageRows] = useState<AirrDetailRow[]>([]);
  const [lineageTotal, setLineageTotal] = useState(0);
  const workbenchRef = useRef<HTMLElement>(null);
  const [alignment, setAlignment] = useState("");
  const [alignmentMode, setAlignmentMode] = useState<"nt" | "aa">("nt");
  const [alignmentMethod, setAlignmentMethod] = useState<"quick" | "kalign" | "codon">("codon");
  const [alignmentLimit, setAlignmentLimit] = useState(200);
  const [treeNewick, setTreeNewick] = useState("");
  const [treeModel, setTreeModel] = useState<"gtr" | "jc">("gtr");
  const [treeFast, setTreeFast] = useState(false);

  const isTcr = scope === "TCR" || String(scope).startsWith("TR");
  const [chmmSegment, setChmmSegment] = useState<ChmmSegment>("V");
  const [chmmMethod, setChmmMethod] = useState<"BW" | "DB">(isTcr ? "DB" : "BW");
  const [chmmSource, setChmmSource] = useState<"selected" | "upload">("selected");
  const [uploadedMsa, setUploadedMsa] = useState("");
  const [uploadedMsaName, setUploadedMsaName] = useState("");
  const [preparedMsa, setPreparedMsa] = useState("");
  const [chmmThreshold, setChmmThreshold] = useState(0.95);
  const [chmmPrior, setChmmPrior] = useState(0.05);
  const [chmmMinDfr, setChmmMinDfr] = useState(1);
  const [chmmDetailed, setChmmDetailed] = useState(false);
  const [mutationRates, setMutationRates] = useState(isTcr ? "0.005" : "0,0.0179,0.0357,0.0536,0.0714,0.0893,0.1071,0.125,0.1429,0.1607,0.1786,0.1964,0.2143,0.2321,0.25");
  const [chmm, setChmm] = useState<ChmmDashboard | null>(null);

  const [queryText, setQueryText] = useState("");
  const [queryTarget, setQueryTarget] = useState<QueryTarget>("cdr3_nt");
  const [queryMetric, setQueryMetric] = useState<QueryMetric>("hamming");
  const [queryIdentity, setQueryIdentity] = useState(0.85);
  const [queryLimit, setQueryLimit] = useState(250);
  const [queryLocus, setQueryLocus] = useState("");
  const [queryV, setQueryV] = useState("");
  const [queryJ, setQueryJ] = useState("");
  const [queryHits, setQueryHits] = useState<QueryHit[]>([]);
  const [queryRecords, setQueryRecords] = useState<Map<number, AirrIndexRecord>>(new Map());
  const [expanded, setExpanded] = useState<{ ordinals: number[]; comparisons: number; capped: boolean } | null>(null);

  async function operation<T>(label: string, action: () => Promise<T>): Promise<T | undefined> {
    setBusy(label);
    setError("");
    try {
      await runtime.ensureIndexed((processed, total) => setProgress({ processed, total }));
      return await action();
    } catch (operationError) {
      if (!(operationError instanceof DOMException && operationError.name === "AbortError")) setError(operationError instanceof Error ? operationError.message : String(operationError));
      return undefined;
    } finally {
      setBusy("");
    }
  }

  async function runDedup() {
    const result = await operation("Deduplicating AIRR records", () => runtime.deduplicate(dedupKey));
    if (result) {
      setDedup(result);
      setUseDedup(true);
      setLineages(null);
    }
  }

  async function downloadDeduplicated() {
    setBusy("Writing deduplicated AIRR table");
    setError("");
    try {
      const counts = await runtime.dedupCounts();
      await saveStream(`${baseName(inputName)}.deduplicated.airr.tsv`, "Deduplicated AIRR rearrangement table", ".tsv", async (writer) => store.writeDeduplicatedAirr(counts, writer.write));
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function runLineages() {
    const result = await operation("Assigning lineages", () => runtime.assignLineages({
      identity,
      callResolution: resolution,
      ambiguity,
      productiveOnly,
      requireSameLocus: true,
      maxCandidateComparisons: candidateCap,
    }, Boolean(dedup && useDedup)));
    if (result) {
      setLineages(result);
      setSelectedLineage(null);
      setLineageRows([]);
      setAlignment("");
      setTreeNewick("");
    }
  }

  async function downloadLineages() {
    setBusy("Writing AIRR table with clone identifiers");
    setError("");
    try {
      const assignments = await runtime.lineageAssignments();
      await saveStream(`${baseName(inputName)}.lineages.airr.tsv`, "AIRR rearrangement table with clone identifiers", ".tsv", async (writer) => store.writeLineageAirr(assignments, writer.write));
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function openLineage(summary: LineageSummary) {
    setBusy(`Loading lineage ${summary.id}`);
    setError("");
    try {
      const members = await runtime.lineageMembers(summary.id, 0, 500);
      const rows = await store.detailMany(members.ordinals);
      setSelectedLineage(summary);
      setLineageRows(rows);
      setLineageTotal(members.total);
      setAlignment("");
      setTreeNewick("");
      window.requestAnimationFrame(() => workbenchRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function runAlignment() {
    if (!lineageRows.length) return;
    setBusy(alignmentMethod === "quick" ? "Preparing AIRR-anchored alignment" : alignmentMethod === "codon" ? "Running codon-aware Kalign WASM" : "Running Kalign WASM");
    setError("");
    try {
      const rows = lineageRows.slice(0, Math.max(2, alignmentLimit));
      let next = "";
      if (alignmentMethod === "quick") next = quickAirrAlignment(rows);
      else {
        const input = lineageInputFasta(rows);
        next = alignmentMethod === "codon" ? await runCodonAwareKalign(input.fasta, input.frames) : await runKalign(input.fasta);
      }
      setAlignment(next);
      setTreeNewick("");
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function inferTree() {
    if (!alignment) return;
    setBusy("Running FastTree WASM and rooting on the N-masked germline");
    setError("");
    try {
      const unrooted = await runFastTree(alignment, treeModel, treeFast);
      const rooted = rootOnOutgroup(parseNewick(unrooted), GERMLINE_OUTGROUP);
      setTreeNewick(`${serializeNewick(rooted)};`);
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function acceptMsa(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    try {
      const text = await file.text();
      prepareReferenceMsa(text);
      setUploadedMsa(text);
      setUploadedMsaName(file.name);
      setChmmSource("upload");
      setPreparedMsa("");
      setError("");
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    }
  }

  async function runChmmAnalysis() {
    setBusy(chmmSource === "selected" ? `Building ${chmmSegment} reference MSA with Kalign WASM` : "Validating uploaded reference MSA");
    setError("");
    try {
      let msa = chmmSource === "upload" ? uploadedMsa : preparedMsa;
      if (!msa) msa = await runKalign(references[chmmSegment]);
      prepareReferenceMsa(msa);
      setPreparedMsa(msa);
      setBusy(`Running CHMMAIRRa ${chmmSegment} model with ${Math.min(workers, 16)} workers`);
      const rates = mutationRates.split(/[\s,;]+/).map(Number).filter((value) => Number.isFinite(value) && value >= 0 && value < 1);
      if (chmmMethod === "DB" && !rates.length) throw new Error("Provide at least one mutation-rate state for the discretized Bayesian model.");
      const result = await runChmmairra(store, msa, {
        segment: chmmSegment,
        method: chmmMethod,
        priorProbability: chmmPrior,
        baseMutationProbability: 0.05,
        mutationRates: rates,
        mutationSwitchProbability: 0,
        detailed: chmmDetailed,
        minDfr: chmmMinDfr,
        threshold: chmmThreshold,
        workers,
      }, (processed, total) => setProgress({ processed, total }));
      setChmm(result);
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function downloadChmm() {
    if (!chmm) return;
    setBusy("Writing CHMMAIRRa result table");
    try {
      await saveStream(`${baseName(inputName)}.chmmairra-${chmm.segment.toLowerCase()}.tsv`, "CHMMAIRRa result table", ".tsv", (writer) => writeChmmairraTsv(store, chmm, writer));
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function runQuery() {
    const queries = parseQueries(queryText);
    if (!queries.length) {
      setError("Paste at least one nucleotide/amino-acid sequence or FASTA record.");
      return;
    }
    const metric = queryTarget === "trimmed" ? "sketch" : queryMetric;
    const hits = await operation("Searching the assigned repertoire", () => runtime.query(queries, {
      target: queryTarget,
      metric,
      identity: queryIdentity,
      maxResults: queryLimit,
      locus: queryLocus || undefined,
      vCall: queryV || undefined,
      jCall: queryJ || undefined,
      callResolution: resolution,
      ambiguity,
      productiveOnly,
    }));
    if (!hits) return;
    setQueryHits(hits);
    setExpanded(null);
    const indexRows = await store.indexRecords([...new Set(hits.map((hit) => hit.ordinal))]);
    setQueryRecords(new Map(indexRows.map((record) => [record.ordinal, record])));
  }

  async function expandMatches() {
    const seeds = expanded?.ordinals.length ? expanded.ordinals : [...new Set(queryHits.map((hit) => hit.ordinal))];
    if (!seeds.length) {
      setError("Run a sequence query before expanding its matched set.");
      return;
    }
    const result = await operation("Expanding a single-linkage CDR3 neighborhood", () => runtime.expand(seeds, {
      identity: queryIdentity,
      callResolution: resolution,
      ambiguity,
      productiveOnly,
      requireSameLocus: true,
      maxCandidateComparisons: candidateCap,
      maxResults: Math.max(queryLimit, 10_000),
    }));
    if (!result) return;
    setExpanded(result);
    const indexRows = await store.indexRecords(result.ordinals.slice(0, 500));
    setQueryRecords(new Map(indexRows.map((record) => [record.ordinal, record])));
  }

  const displayedQueryOrdinals = expanded?.ordinals ?? queryHits.map((hit) => hit.ordinal);
  const displayedQueryRows = [...new Set(displayedQueryOrdinals)].slice(0, 500).flatMap((ordinal) => queryRecords.get(ordinal) ?? []);
  const vChart = lineages?.vUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];
  const jChart = lineages?.jUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];

  return <section className="post-analysis-shell">
    <header className="post-analysis-heading"><div><span className="section-kicker">Post-assignment methods</span><h2>Repertoire structure and targeted phylogenetics</h2><p>Each method is opt-in. AIRR chunks are scanned in bounded batches; lineage candidate generation uses V/J and CDR3-length partitions before exact nucleotide-distance checks.</p></div><div className="local-method-note"><span>Execution</span><strong>Browser-local</strong><small>Input, germlines, and results are not submitted to an analysis server.</small></div></header>

    <div className="post-method-map"><article><b>01</b><span>Collapse</span><strong>Deduplicate + retain count</strong></article><article><b>02</b><span>QC</span><strong>Optional CHMMAIRRa</strong></article><article><b>03</b><span>Repertoire</span><strong>Assign lineages</strong></article><article><b>04</b><span>Targeted</span><strong>Query + expand</strong></article><article><b>05</b><span>On demand</span><strong>Align + infer tree</strong></article></div>

    {busy && <div className="post-progress" role="status"><div><span>{busy}</span><strong>{progress.total ? `${Math.min(100, progress.processed / progress.total * 100).toFixed(1)}%` : "working"}</strong></div><progress max={Math.max(1, progress.total)} value={progress.processed} /><small>{progress.processed.toLocaleString()} / {progress.total.toLocaleString()} AIRR records indexed or scanned</small></div>}
    {error && <div className="post-error" role="alert"><strong>Post-analysis stopped</strong><span>{error}</span><button type="button" onClick={() => setError("")}>Dismiss</button></div>}

    <section className="post-module dedup-module">
      <header><div className="module-number">01</div><div><span className="section-kicker">Abundance preservation</span><h3>Deduplicate records</h3><p>One representative is retained for each key; <code>duplicate_count</code> carries the collapsed abundance into lineage sizes and export.</p></div></header>
      <div className="module-controls">
        <label><span>Identity key</span><select value={dedupKey} onChange={(event) => setDedupKey(event.target.value as DedupKey)}><option value="sequence">Full input sequence</option><option value="trimmed">VDJ-aligned sequence</option><option value="cdr3">Locus + CDR3 nucleotide</option><option value="rearrangement">Locus + V/J calls + CDR3</option></select></label>
        <button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runDedup()}>Run deduplication</button>
      </div>
      {(dedupKey === "sequence" || dedupKey === "trimmed") && <p className="scientific-note"><span>i</span>Sequence-key modes compare normalized length plus a 128-bit fingerprint so complete sequence payloads do not remain in memory. CDR3 and rearrangement modes retain and compare their exact key strings.</p>}
      {dedup && <div className="module-result"><div className="post-stat-grid"><article><span>Input records</span><strong>{dedup.inputRecords.toLocaleString()}</strong></article><article><span>Unique representatives</span><strong>{dedup.uniqueRecords.toLocaleString()}</strong></article><article><span>Collapsed duplicates</span><strong>{dedup.collapsedRecords.toLocaleString()}</strong></article><article><span>Largest abundance</span><strong>{(dedup.largestGroups[0]?.count ?? 1).toLocaleString()}</strong></article></div><div className="result-actions"><button type="button" onClick={() => void downloadDeduplicated()}>Download deduplicated AIRR + counts</button><label><input type="checkbox" checked={useDedup} onChange={(event) => setUseDedup(event.target.checked)} /> Use representatives and abundance for lineage assignment</label></div></div>}
    </section>

    <section className="post-module chmm-module">
      <header><div className="module-number amber">02</div><div><span className="section-kicker">Optional PCR-chimera model</span><h3>CHMMAIRRa after V(D)J assignment</h3><p>The browser port threads each AIRR local V or J alignment onto a reference MSA, then evaluates the CHMMera posterior. V is the manuscript default; D is not modeled.</p></div><a href="https://github.com/MurrellGroup/CHMMAIRRa.jl" target="_blank" rel="noreferrer">Method source ↗</a></header>
      <div className="chmm-grid">
        <div className="chmm-config">
          <div className="control-grid three"><label><span>Segment</span><select value={chmmSegment} onChange={(event) => { setChmmSegment(event.target.value as ChmmSegment); setPreparedMsa(""); setChmm(null); }}><option value="V">V (recommended)</option><option value="J">J (optional)</option></select></label><label><span>Model</span><select value={chmmMethod} onChange={(event) => setChmmMethod(event.target.value as "BW" | "DB")}><option value="BW">Baum–Welch · IG default</option><option value="DB">Discretized Bayesian · TCR default</option></select></label><label><span>Posterior threshold</span><input type="number" min="0" max="1" step="0.01" value={chmmThreshold} onChange={(event) => setChmmThreshold(Number(event.target.value))} /></label></div>
          <fieldset className="msa-source"><legend>Reference multiple-sequence alignment</legend><label className={chmmSource === "selected" ? "selected" : ""}><input type="radio" checked={chmmSource === "selected"} onChange={() => setChmmSource("selected")} /><span><strong>Build from this run’s {chmmSegment} references</strong><small>Kalign 3.3.1 WASM; preserves selected IMGT/KI/uploaded composition.</small></span></label><label className={chmmSource === "upload" ? "selected" : ""}><input type="radio" checked={chmmSource === "upload"} onChange={() => setChmmSource("upload")} /><span><strong>Use an aligned FASTA MSA</strong><small>{uploadedMsaName || "Every record must have equal aligned length and names matching AIRR calls."}</small></span><input className="file-inline" type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event) => void acceptMsa(event)} /></label></fieldset>
          <details className="post-advanced"><summary>Model parameters</summary><div className="control-grid three"><label><span>Chimera prior</span><input type="number" min="0.00001" max="0.5" step="0.01" value={chmmPrior} onChange={(event) => setChmmPrior(Number(event.target.value))} /></label><label><span>Minimum DFR</span><input type="number" min="0" max="100" step="1" value={chmmMinDfr} onChange={(event) => setChmmMinDfr(Number(event.target.value))} /></label><label><span>DB mutation rates</span><input type="text" value={mutationRates} onChange={(event) => setMutationRates(event.target.value)} /></label><label className="check-line"><input type="checkbox" checked={chmmDetailed} onChange={(event) => setChmmDetailed(event.target.checked)} /><span>Infer Viterbi parents and breakpoints for evaluated records</span></label></div></details>
          <div className="scientific-note warning"><span>!</span><p>Reference completeness matters: an absent true V/J allele can produce a false switch signal. Uploaded MSAs are never silently supplemented. Low-DFR records are reported as unevaluated rather than forced through the model.</p></div>
          <div className="result-actions"><button className="post-primary amber" type="button" disabled={Boolean(busy) || (chmmSource === "upload" && !uploadedMsa)} onClick={() => void runChmmAnalysis()}>Run CHMMAIRRa</button>{preparedMsa && <button type="button" onClick={() => downloadText(preparedMsa, `${baseName(inputName)}.${chmmSegment.toLowerCase()}-reference-msa.fasta`)}>Download reference MSA</button>}</div>
        </div>
        <div className="chmm-result-panel">
          {chmm ? <><div className="post-stat-grid compact"><article><span>Evaluated</span><strong>{chmm.evaluated.toLocaleString()}</strong></article><article><span>Posterior ≥ {chmm.threshold}</span><strong>{chmm.flagged.toLocaleString()}</strong></article><article><span>Below DFR</span><strong>{chmm.lowDfr.toLocaleString()}</strong></article><article><span>Missing reference</span><strong>{chmm.missingReference.toLocaleString()}</strong></article></div><BarChart title={`${chmm.segment} chimera posterior`} subtitle="Evaluated rearrangements by posterior interval" data={chmm.histogram.map((item) => ({ label: item.label, value: item.count }))} color="#d49a19" name={`${baseName(inputName)}.chmmairra-posterior.svg`} />{chmm.top.length > 0 && <div className="chmm-top-list"><header><strong>Highest posteriors</strong><span>Click a record for AIRR evidence</span></header>{chmm.top.slice(0, 12).map((record) => <button type="button" key={record.ordinal} onClick={() => onInspect(record.ordinal)}><b>#{(record.ordinal + 1).toLocaleString()}</b><span>{(record.probability * 100).toFixed(2)}%</span><small>DFR {record.dfr}{record.recombinations.length ? ` · ${record.recombinations.map((event) => `${event.left}→${event.right}@${event.position}`).join("; ")}` : ""}</small></button>)}</div>}<button type="button" onClick={() => void downloadChmm()}>Download CHMMAIRRa TSV</button></> : <div className="method-placeholder"><span>HMM</span><h4>No CHMMAIRRa run</h4><p>{isTcr ? "TCR mode defaults to a fixed 0.005 mutation-rate state (DB)." : "IG mode defaults to per-reference Baum–Welch mutation estimates."}</p></div>}
        </div>
      </div>
    </section>

    <section className="post-module lineage-module">
      <header><div className="module-number">03</div><div><span className="section-kicker">Repertoire-scale clonal grouping</span><h3>Assign lineages from CDR3 nucleotide distance</h3><p>Default: same locus, overlapping V/J gene assignments, exact CDR3 nucleotide length, and single-linkage at ≥85% identity. The threshold is a starting point and remains dataset-adjustable.</p></div><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5340603/" target="_blank" rel="noreferrer">Clonal threshold literature ↗</a></header>
      <div className="lineage-config"><div className="control-grid five"><label><span>CDR3 identity</span><div className="range-number"><input type="range" min="0.7" max="1" step="0.01" value={identity} onChange={(event) => setIdentity(Number(event.target.value))} /><b>{Math.round(identity * 100)}%</b></div></label><label><span>Call level</span><select value={resolution} onChange={(event) => setResolution(event.target.value as CallResolution)}><option value="gene">Gene</option><option value="allele">Allele</option></select></label><label><span>Ambiguous calls</span><select value={ambiguity} onChange={(event) => setAmbiguity(event.target.value as AmbiguityPolicy)}><option value="overlap">Any assignment overlaps</option><option value="top">Top call only</option><option value="strict">Exact call sets</option></select></label><label><span>Candidate cap / record</span><input type="number" min="100" max="1000000" step="1000" value={candidateCap} onChange={(event) => setCandidateCap(Number(event.target.value))} /></label><label className="check-line"><input type="checkbox" checked={productiveOnly} onChange={(event) => setProductiveOnly(event.target.checked)} /><span>Productive only</span></label></div><div className="algorithm-note"><strong>Exact accelerated single-linkage</strong><span>Partition by locus → V/J calls → CDR3 length → d+1 exact blocks → verify normalized Hamming distance → union-find components.</span></div><button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runLineages()}>Assign lineages</button></div>
      {lineages && <div className="lineage-results"><div className="post-stat-grid"><article><span>Lineages</span><strong>{lineages.lineageCount.toLocaleString()}</strong></article><article><span>Assigned records</span><strong>{lineages.assignedRecords.toLocaleString()}</strong></article><article><span>Largest lineage</span><strong>{(lineages.summaries[0]?.abundance ?? 0).toLocaleString()}</strong></article><article><span>Exact comparisons</span><strong>{lineages.candidateComparisons.toLocaleString()}</strong></article></div>{lineages.truncatedCandidates > 0 && <div className="scientific-note warning"><span>!</span><p>{lineages.truncatedCandidates.toLocaleString()} records reached the candidate cap. Increase it and rerun before treating components as complete.</p></div>}<div className="result-actions"><button type="button" onClick={() => void downloadLineages()}>Download AIRR + clone_id</button></div><div className="chart-customizer"><label><span>Gene chart metric</span><select value={geneMetric} onChange={(event) => setGeneMetric(event.target.value as "abundance" | "lineages")}><option value="abundance">Sequence abundance</option><option value="lineages">Lineage count</option></select></label><label><span>Top genes</span><input type="number" min="5" max="24" value={topGenes} onChange={(event) => setTopGenes(Number(event.target.value))} /></label><label><span>Figure color</span><input type="color" value={chartColor} onChange={(event) => setChartColor(event.target.value)} /></label></div><div className="post-chart-grid"><BarChart title="Lineage abundance distribution" subtitle="Lineage count in each abundance interval" data={lineages.sizeHistogram.map((item) => ({ label: item.label, value: item.count }))} color={chartColor} name={`${baseName(inputName)}.lineage-size-distribution.svg`} /><BarChart title="Largest lineages" subtitle="Abundance retained after deduplication" data={lineages.summaries.slice(0, 20).map((item) => ({ label: `Lineage ${item.id}`, value: item.abundance }))} color={chartColor} name={`${baseName(inputName)}.largest-lineages.svg`} /><BarChart title="V germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={vChart} color={chartColor} name={`${baseName(inputName)}.lineage-v-use.svg`} /><BarChart title="J germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={jChart} color={chartColor} name={`${baseName(inputName)}.lineage-j-use.svg`} /></div><div className="lineage-table-wrap"><table><thead><tr><th>Lineage</th><th>Abundance</th><th>Unique</th><th>Locus</th><th>V calls</th><th>J calls</th><th>CDR3 nt</th><th /></tr></thead><tbody>{lineages.summaries.slice(0, 250).map((summary) => <tr key={summary.id} className={selectedLineage?.id === summary.id ? "selected" : ""} onClick={() => void openLineage(summary)}><td><strong>{summary.id}</strong></td><td>{summary.abundance.toLocaleString()}</td><td>{summary.uniqueMembers.toLocaleString()}</td><td>{summary.locus}</td><td>{summary.vCalls.join(", ")}</td><td>{summary.jCalls.join(", ")}</td><td>{summary.cdr3Length} nt</td><td><button type="button">Open →</button></td></tr>)}</tbody></table></div></div>}
    </section>

    <section className="post-module query-module">
      <header><div className="module-number dark">04</div><div><span className="section-kicker">Targeted repertoire search</span><h3>Query sequences, then expand the matched set</h3><p>Paste one or more sequences or FASTA records. Search CDR3 directly or use a compact k-mer sketch of the VDJ-aligned sequence; expansion follows exact CDR3 single-linkage edges from the current seed set.</p></div></header>
      <div className="query-layout"><div className="query-input"><label><span>Query sequence(s) or FASTA</span><textarea value={queryText} onChange={(event) => setQueryText(event.target.value)} placeholder=">seed_1\nTGTGCGAGAGAT…\n>seed_2\nTGTGCGAGGGAT…" /></label><div className="control-grid three"><label><span>Target</span><select value={queryTarget} onChange={(event) => setQueryTarget(event.target.value as QueryTarget)}><option value="cdr3_nt">CDR3 nucleotide</option><option value="cdr3_aa">CDR3 amino acid</option><option value="trimmed">VDJ-aligned k-mer sketch</option></select></label><label><span>Metric</span><select disabled={queryTarget === "trimmed"} value={queryTarget === "trimmed" ? "sketch" : queryMetric} onChange={(event) => setQueryMetric(event.target.value as QueryMetric)}><option value="exact">Exact</option><option value="substring">Substring</option><option value="hamming">Hamming</option><option value="edit">Banded edit distance</option><option value="sketch">MinHash k-mer estimate</option></select></label><label><span>Identity / similarity</span><input type="number" min="0" max="1" step="0.01" value={queryIdentity} onChange={(event) => setQueryIdentity(Number(event.target.value))} /></label><label><span>Locus constraint</span><select value={queryLocus} onChange={(event) => setQueryLocus(event.target.value)}><option value="">Any locus</option>{loci.map((item) => <option value={item.value} key={item.value}>{item.value}</option>)}</select></label><label><span>V constraint</span><input value={queryV} onChange={(event) => setQueryV(event.target.value)} placeholder="optional IGHV…" /></label><label><span>J constraint</span><input value={queryJ} onChange={(event) => setQueryJ(event.target.value)} placeholder="optional IGHJ…" /></label><label><span>Maximum initial hits</span><input type="number" min="1" max="10000" value={queryLimit} onChange={(event) => setQueryLimit(Number(event.target.value))} /></label></div><div className="result-actions"><button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void runQuery()}>Search repertoire</button><button type="button" disabled={Boolean(busy) || !queryHits.length} onClick={() => void expandMatches()}>Single-linkage expand set</button></div><p className="scientific-note"><span>i</span>{queryTarget === "trimmed" ? "VDJ searches lazily build a packed index of eight independent 7-mer MinHash values per record. Scores are approximate and intended for candidate retrieval." : "Hamming search requires equal length. Edit distance is banded by the selected identity. Expansion always uses equal-length CDR3 nucleotide Hamming edges."}</p></div><div className="query-results"><header><div><span className="section-kicker">Matched set</span><h4>{expanded ? `${expanded.ordinals.length.toLocaleString()} expanded records` : `${queryHits.length.toLocaleString()} initial hits`}</h4><p>{expanded ? `${expanded.comparisons.toLocaleString()} exact edge checks${expanded.capped ? " · result cap reached" : " · fixed point reached"}` : "Ranked by similarity, then distance."}</p></div></header><div className="query-result-list">{displayedQueryRows.slice(0, 100).map((record) => { const hit = queryHits.find((value) => value.ordinal === record.ordinal); return <button type="button" key={record.ordinal} onClick={() => onInspect(record.ordinal)}><span><strong>{record.sequenceId}</strong><small>{record.locus} · {record.vCall || "V—"} · {record.jCall || "J—"}</small></span><code>{record.cdr3 || record.cdr3Aa || "—"}</code><b>{hit ? `${(hit.score * 100).toFixed(1)}%` : "expanded"}</b></button>; })}{!displayedQueryRows.length && <div className="method-placeholder small"><span>⌕</span><h4>No query results</h4><p>Provide a sequence and search the assigned repertoire.</p></div>}</div></div></div>
    </section>

    {selectedLineage && <section ref={workbenchRef} className="post-module lineage-workbench" tabIndex={-1}>
      <header><div className="module-number dark">05</div><div><span className="section-kicker">Selected lineage {selectedLineage.id}</span><h3>Alignment and rooted phylogeny</h3><p>{lineageTotal.toLocaleString()} records · {selectedLineage.uniqueMembers.toLocaleString()} unique representatives · {selectedLineage.locus} · {selectedLineage.vCalls.join(", ")} / {selectedLineage.jCalls.join(", ")}</p></div><button type="button" onClick={() => { setSelectedLineage(null); setLineageRows([]); }}>Close lineage</button></header>
      <div className="member-strip"><div><strong>{lineageRows.length.toLocaleString()}</strong><span>records loaded</span><small>{lineageTotal > lineageRows.length ? `first ${lineageRows.length}; alignment remains on-demand` : "complete selected lineage"}</small></div><div className="member-pills">{lineageRows.slice(0, 18).map((row) => <button type="button" key={row.record.ordinal} onClick={() => onInspect(row.record.ordinal)}>{row.record.sequenceId}</button>)}</div></div>
      <div className="alignment-controls"><label><span>Alignment method</span><select value={alignmentMethod} onChange={(event) => setAlignmentMethod(event.target.value as "quick" | "kalign" | "codon")}><option value="quick">AIRR-anchored quick view</option><option value="kalign">Kalign 3.3.1 WASM · nucleotide</option><option value="codon">Kalign 3.3.1 WASM · codon-aware</option></select></label><label><span>Maximum sequences</span><input type="number" min="2" max="500" value={alignmentLimit} onChange={(event) => setAlignmentLimit(Number(event.target.value))} /></label><button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runAlignment()}>Align selected lineage</button>{alignment && <button type="button" onClick={() => downloadText(alignment, `${baseName(inputName)}.lineage-${selectedLineage.id}.alignment.fasta`)}>Alignment FASTA ↓</button>}</div>
      {alignment && <><div className="alignment-view-controls"><div className="mode-toggle"><button className={alignmentMode === "nt" ? "active" : ""} type="button" onClick={() => setAlignmentMode("nt")}>Nucleotide</button><button className={alignmentMode === "aa" ? "active" : ""} type="button" onClick={() => setAlignmentMode("aa")}>Amino acid</button></div><span>{parseFasta(alignment, true).length.toLocaleString()} aligned rows · N-masked germline included</span></div><AlignmentPreview fasta={alignment} mode={alignmentMode} /><div className="tree-controls"><div><span className="section-kicker">On-demand tree</span><h4>FastTree 2.1.11 WASM</h4><p>Nucleotide inference only; the N-masked germline is used as the outgroup after inference.</p></div><label><span>Model</span><select value={treeModel} onChange={(event) => setTreeModel(event.target.value as "gtr" | "jc")}><option value="gtr">GTR</option><option value="jc">Jukes–Cantor</option></select></label><label className="check-line"><input type="checkbox" checked={treeFast} onChange={(event) => setTreeFast(event.target.checked)} /><span>Fastest heuristic</span></label><button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void inferTree()}>Infer + root tree</button></div>{treeNewick && <TreeView newick={treeNewick} name={`${baseName(inputName)}.lineage-${selectedLineage.id}.rooted-tree`} />}</>}
    </section>}
  </section>;
}
