import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type ChangeEvent,
  type ReactNode,
} from "react";

import { runCodonAwareKalign, runFastTree, runKalign, type FastTreeRun } from "./biowasm-runtime";
import { inspectAlignment, validateCorrectedAlignment } from "./alignment-provenance";
import { chimeraVisiblePositions, classifyChimeraQuerySite } from "./chimera-view-model";
import { CommitNumberInput } from "./commit-number-input";
import { CommitTextInput, CommitTextarea } from "./commit-text-input";
import {
  runChmmairra,
  runChmmairraDetail,
  writeChmmairraTsv,
  type ChmmDashboard,
  type ChmmDetail,
  type ChmmRunOptions,
  type ChmmSegment,
} from "./chmmairra-runtime";
import { GERMLINE_OUTGROUP, lineageInputFasta, quickAirrAlignment } from "./lineage-alignment";
import { LineageTreeViewer } from "./lineage-tree-viewer";
import {
  canonicalizeTree,
  collapseShortInternalBranches,
  parseNewick,
  rootOnOutgroup,
  serializeNewick,
} from "./phylogeny";
import {
  parseFasta,
  prepareReferenceMsa,
  type AmbiguityPolicy,
  type CallResolution,
  type CollapseMode,
  type DedupKey,
  type LineageSummary,
  type QueryHit,
  type QueryConstraint,
  type QueryMetric,
  type QueryTarget,
} from "./post-analysis-core";
import {
  PostAnalysisRuntime,
  type DedupDashboard,
  type LineageDashboard,
} from "./post-analysis-runtime";
import { inferQueryAssignments, type InferredQueryAssignment } from "./query-inference-runtime";
import type { CompiledReferences, ScopeKey } from "./reference-pack";
import type { AirrDetailRow, AirrIndexRecord, AirrResultStore, FacetValue } from "./result-store";
import { ColoredSequence, sequenceColor } from "./sequence-colors";

interface Props {
  store: AirrResultStore;
  references: CompiledReferences;
  scope: ScopeKey;
  loci: FacetValue[];
  inputName: string;
  workers: number;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  onInspect: (ordinal: number) => void;
}

interface SaveFileHandle {
  createWritable: () => Promise<{ write: (value: string | Blob | Uint8Array) => Promise<void>; close: () => Promise<void>; abort?: () => Promise<void> }>;
}

interface ChartDatum {
  label: string;
  value: number;
}

interface WorkingSetStage {
  id: "dedup" | "chimera";
  label: string;
  input: number;
  retained: number;
  discarded: number;
  detail: string;
}

interface TreeSnapshot extends FastTreeRun {
  rootedNewick: string;
  stableNewick: string;
  collapsedEdges: number;
  collapseThreshold: number;
  source: string;
}

type TreeViewMode = "stable" | "rooted" | "raw";

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
    {records.slice(0, 80).map((record) => <div className={record.name === GERMLINE_OUTGROUP ? "germline-row" : ""} key={record.name}><strong title={record.name}>{record.name}</strong><ColoredSequence sequence={mode === "nt" ? record.sequence : translateAligned(record.sequence)} alphabet={mode} /></div>)}
  </div>;
}

function ChimeraHighlighter({ detail, name, onInspect }: { detail: ChmmDetail; name: string; onInspect: (ordinal: number) => void }) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [colorMode, setColorMode] = useState<"nucleotide" | "highlighter">("highlighter");
  const [eventIndex, setEventIndex] = useState(0);
  useEffect(() => setEventIndex(0), [detail.ordinal]);
  const cellWidth = 12;
  const labelWidth = 210;
  const rowHeight = 27;
  const top = 45;
  const selectedEvent = detail.recombinations[Math.min(eventIndex, Math.max(0, detail.recombinations.length - 1))];
  const parentIndex = (parentName: string | undefined, fallback: number) => {
    const found = detail.parents.findIndex((parent) => parent.name === parentName);
    return found >= 0 ? found : Math.min(fallback, Math.max(0, detail.parents.length - 1));
  };
  const parentAIndex = parentIndex(selectedEvent?.left, 0);
  let parentBIndex = parentIndex(selectedEvent?.right, parentAIndex === 0 ? 1 : 0);
  if (parentBIndex === parentAIndex && detail.parents.length > 1) parentBIndex = parentAIndex === 0 ? 1 : 0;
  const parentA = detail.parents[parentAIndex] ?? { name: "Parent A unavailable", sequence: "" };
  const parentB = detail.parents[parentBIndex] ?? { name: "Parent B unavailable", sequence: "" };
  const query = detail.threadedObservation;
  const visiblePositions = chimeraVisiblePositions(query, parentA.sequence, parentB.sequence);
  const hiddenTripleGaps = query.length - visiblePositions.length;
  const rows = [
    { kind: "parent" as const, label: `PARENT A · ${parentA.name}`, sequence: parentA.sequence },
    { kind: "query" as const, label: `QUERY · ${detail.sequenceId}`, sequence: query },
    { kind: "parent" as const, label: `PARENT B · ${parentB.name}`, sequence: parentB.sequence },
  ];
  const alignmentBottom = top + (rows.length - 1) * rowHeight + rowHeight - 16;
  const breakpointLabelStart = alignmentBottom + 22;
  const displayBoundary = (position: number) => {
    const original = Math.max(0, position - 1);
    const index = visiblePositions.findIndex((value) => value >= original);
    return index < 0 ? visiblePositions.length : index;
  };
  const highlighterColor = (position: number) => {
    const category = classifyChimeraQuerySite(query[position] ?? "-", parentA.sequence[position] ?? "-", parentB.sequence[position] ?? "-");
    return category === "parent_a" ? "#d84a4a" : category === "parent_b" ? "#4277c7" : "#cfd3d1";
  };
  const width = labelWidth + visiblePositions.length * cellWidth + 18;
  const height = breakpointLabelStart + Math.max(1, detail.recombinations.length) * 10 + 24;
  const informative = [...detail.threadedObservation].filter((base) => base !== "-" && base !== "N").length;
  return <section className="chimera-detail" tabIndex={-1}>
    <header><div><span className="section-kicker">On-demand Viterbi reconstruction · record #{(detail.ordinal + 1).toLocaleString()}</span><h4>{detail.sequenceId}</h4><p>{detail.segment} call {detail.call || "unassigned"} · posterior {(detail.probability * 100).toFixed(3)}% · DFR {detail.dfr} · {informative} informative MSA columns · {hiddenTripleGaps} triple-gap sites hidden</p></div><div className="post-chart-actions"><button type="button" onClick={() => onInspect(detail.ordinal)}>Open AIRR record</button><button type="button" onClick={() => saveSvg(svgRef.current, name)}>SVG ↓</button></div></header>
    <div className="chimera-viz-controls"><div className="mode-toggle"><button className={colorMode === "nucleotide" ? "active" : ""} type="button" onClick={() => setColorMode("nucleotide")}>Nucleotide colors</button><button className={colorMode === "highlighter" ? "active" : ""} type="button" onClick={() => setColorMode("highlighter")}>Parent-match highlighter</button></div>{detail.recombinations.length > 1 && <label><span>Parent pair / breakpoint</span><select value={eventIndex} onChange={(event) => setEventIndex(Number(event.target.value))}>{detail.recombinations.map((event, index) => <option value={index} key={`${event.position}-${index}`}>{index + 1}. {event.left} → {event.right} at MSA {event.position}</option>)}</select></label>}<p>Tile colors are computed from nucleotide identity only; the Viterbi path is never used for coloring.</p></div>
    {colorMode === "highlighter" && <div className="chimera-breakpoint-strip"><span className="parent-a-key"><i />red · query matches parent A only</span><span className="parent-b-key"><i />blue · query matches parent B only</span><span className="neutral-key"><i />gray · matches both or neither</span></div>}
    <div className="chimera-alignment-scroll"><svg ref={svgRef} width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label="Viterbi parent highlighter nucleotide alignment">
      <rect width={width} height={height} fill="#fbfaf5" />
      <text x="12" y="22" fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700" fill="#34433e">CHMMAIRRa parent-match alignment · {colorMode === "nucleotide" ? "standard nucleotide palette" : "literal query/parent matches"}</text>
      {visiblePositions.filter((_, index) => index % 10 === 0).map((position, tick) => <g key={position}><line x1={labelWidth + tick * 10 * cellWidth} x2={labelWidth + tick * 10 * cellWidth} y1={top - 11} y2={alignmentBottom} stroke={tick % 5 === 0 ? "#9aa7a1" : "#d9ddd8"} strokeWidth={tick % 5 === 0 ? 1 : 0.5} /><text x={labelWidth + tick * 10 * cellWidth + 2} y={top - 16} fontFamily="ui-monospace,monospace" fontSize="8" fill="#6b7974">{position + 1}</text></g>)}
      {rows.map((row, rowIndex) => {
        const y = top + rowIndex * rowHeight;
        return <g key={`${row.kind}-${row.label}`}>
          <rect x="0" y={y - 14} width={labelWidth - 5} height={rowHeight - 2} fill={row.kind === "query" ? "#172622" : "#e9ece7"} />
          <text x="10" y={y + 3} fontFamily="ui-monospace,monospace" fontSize="8" fontWeight="700" fill={row.kind === "query" ? "#f4f6f3" : "#384842"}>{row.label.length > 31 ? `${row.label.slice(0, 29)}…` : row.label}<title>{row.label}</title></text>
          {visiblePositions.map((position, displayIndex) => {
            const base = row.sequence[position] ?? "-";
            const fill = colorMode === "nucleotide" ? sequenceColor(base, "nt") : row.kind === "query" ? highlighterColor(position) : "#f7f6f1";
            const dark = fill === "#d84a4a" || fill === "#4277c7";
            return <g key={position}>
              <rect x={labelWidth + displayIndex * cellWidth} y={y - 14} width={cellWidth} height={rowHeight - 2} fill={fill} />
              <text x={labelWidth + displayIndex * cellWidth + cellWidth / 2} y={y + 4} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize="10" fontWeight="700" fill={dark ? "#ffffff" : "#17231f"}>{base}</text>
            </g>;
          })}
        </g>;
      })}
      {detail.recombinations.map((event, index) => {
        const x = labelWidth + displayBoundary(event.position) * cellWidth;
        const selected = index === eventIndex || detail.recombinations.length === 1;
        return <g key={`${event.position}-${index}`}><line x1={x} x2={x} y1={top - 12} y2={alignmentBottom + 7} stroke={selected ? "#7c2435" : "#8c9692"} strokeWidth={selected ? 2.4 : 1.2} strokeDasharray="4 3" /><text x={x + 4} y={breakpointLabelStart + index * 10} fontFamily="ui-monospace,monospace" fontSize="8" fontWeight={selected ? "700" : "500"} fill={selected ? "#84253a" : "#66736e"}>Viterbi break {event.position}: {event.left.replace(/\s.*/, "")} → {event.right.replace(/\s.*/, "")}</text></g>;
      })}
      <text x="12" y={height - 11} fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#65736e">Pure A/query/B gap columns are hidden. Dashed lines are Viterbi breakpoints; tile colors are independent of the inferred path.</text>
    </svg></div>
  </section>;
}

function parseQueries(text: string): string[] {
  if (text.trimStart().startsWith(">")) return parseFasta(text).map((record) => record.sequence).filter(Boolean);
  return text.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith("#")).map((line) => line.replace(/\s/g, ""));
}

export function PostAnalysisWorkbench({ store, references, scope, loci, inputName, workers, minimumIdentity, strand, onInspect }: Props) {
  const runtime = useMemo(() => new PostAnalysisRuntime(store), [store]);
  useEffect(() => () => runtime.terminate(), [runtime]);
  const [busy, setBusy] = useState("");
  const [progress, setProgress] = useState({ processed: 0, total: store.count });
  const [error, setError] = useState("");

  const [dedupKey, setDedupKey] = useState<DedupKey>("sequence");
  const [collapseMode, setCollapseMode] = useState<CollapseMode>("exact");
  const [dedup, setDedup] = useState<DedupDashboard | null>(null);
  const [denoiseErrorRate, setDenoiseErrorRate] = useState(0.00473);
  const [denoiseAlpha, setDenoiseAlpha] = useState(0.01);
  const [denoiseResolution, setDenoiseResolution] = useState<CallResolution>("allele");
  const [denoiseAmbiguity, setDenoiseAmbiguity] = useState<"top" | "strict">("strict");
  const [minimumParentCount, setMinimumParentCount] = useState(2);
  const [denoiseAmbiguousPolicy, setDenoiseAmbiguousPolicy] = useState<"exclude" | "retain">("exclude");
  const [fadNeighborThreshold, setFadNeighborThreshold] = useState(1);
  const [fadMethod, setFadMethod] = useState<1 | 2>(2);
  const [expectedZeroErrorFraction, setExpectedZeroErrorFraction] = useState(1);
  const [maximumDenoiseDistance, setMaximumDenoiseDistance] = useState(1);
  const [maximumEditDistance, setMaximumEditDistance] = useState(2);
  const [minimumIndelParentRatio, setMinimumIndelParentRatio] = useState(2);
  const [denoiseCandidateCap, setDenoiseCandidateCap] = useState(50_000);
  const [workingMask, setWorkingMask] = useState<Uint8Array | null>(null);
  const [workingStages, setWorkingStages] = useState<WorkingSetStage[]>([]);

  const [identity, setIdentity] = useState(0.85);
  const [resolution, setResolution] = useState<CallResolution>("gene");
  const [ambiguity, setAmbiguity] = useState<AmbiguityPolicy>("overlap");
  const [productiveOnly, setProductiveOnly] = useState(true);
  const [candidateCap, setCandidateCap] = useState(50_000);
  const [lineages, setLineages] = useState<LineageDashboard | null>(null);
  const [geneMetric, setGeneMetric] = useState<"abundance" | "lineages">("abundance");
  const [chartColor, setChartColor] = useState("#08796f");
  const [topGenes, setTopGenes] = useState(15);

  const [selectedLineage, setSelectedLineage] = useState<LineageSummary | null>(null);
  const [lineageRows, setLineageRows] = useState<AirrDetailRow[]>([]);
  const [lineageMultiplicity, setLineageMultiplicity] = useState<Map<number, number>>(new Map());
  const [lineageTotal, setLineageTotal] = useState(0);
  const workbenchRef = useRef<HTMLElement>(null);
  const [alignment, setAlignment] = useState("");
  const [alignmentMode, setAlignmentMode] = useState<"nt" | "aa">("nt");
  const [alignmentMethod, setAlignmentMethod] = useState<"quick" | "kalign" | "codon">("codon");
  const [alignmentLimit, setAlignmentLimit] = useState(200);
  const [alignmentSource, setAlignmentSource] = useState("");
  const alignmentRevisionRef = useRef(0);
  const [treeRun, setTreeRun] = useState<TreeSnapshot | null>(null);
  const [treeViewMode, setTreeViewMode] = useState<TreeViewMode>("rooted");
  const [treeError, setTreeError] = useState("");
  const [treeModel, setTreeModel] = useState<"gtr" | "jc">("gtr");
  const [treeFast, setTreeFast] = useState(false);
  const treeResultRef = useRef<HTMLDivElement>(null);
  const alivibeWindowRef = useRef<Window | null>(null);
  const [alignmentEditorStatus, setAlignmentEditorStatus] = useState("");
  const [alignmentEditorError, setAlignmentEditorError] = useState("");
  const alignmentInfo = useMemo(() => {
    if (!alignment) return null;
    try { return inspectAlignment(alignment); } catch { return null; }
  }, [alignment]);

  function clearAlignmentArtifacts() {
    alignmentRevisionRef.current += 1;
    setAlignment("");
    setAlignmentSource("");
    setTreeRun(null);
    setTreeError("");
  }

  function installAlignment(next: string, source: string) {
    const inspected = inspectAlignment(next);
    alignmentRevisionRef.current += 1;
    setAlignment(inspected.fasta);
    setAlignmentSource(source);
    setTreeRun(null);
    setTreeError("");
    return inspected;
  }

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
  const [chmmRun, setChmmRun] = useState<{ msa: string; options: ChmmRunOptions; inputMask: Uint8Array | null } | null>(null);
  const [chmmFilterThreshold, setChmmFilterThreshold] = useState(0.95);
  const [retainUnevaluated, setRetainUnevaluated] = useState(true);
  const [chimeraDetail, setChimeraDetail] = useState<ChmmDetail | null>(null);
  const chimeraDetailRef = useRef<HTMLDivElement>(null);

  const [queryText, setQueryText] = useState("");
  const [queryTarget, setQueryTarget] = useState<QueryTarget>("cdr3_nt");
  const [queryMetric, setQueryMetric] = useState<QueryMetric>("hamming");
  const [queryIdentity, setQueryIdentity] = useState(0.85);
  const [queryLimit, setQueryLimit] = useState(250);
  const [queryLocus, setQueryLocus] = useState("");
  const [queryV, setQueryV] = useState("");
  const [queryJ, setQueryJ] = useState("");
  const [queryConstraintMode, setQueryConstraintMode] = useState<"manual" | "infer">("manual");
  const [queryInference, setQueryInference] = useState<InferredQueryAssignment[]>([]);
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
    setProgress({ processed: 0, total: store.count });
    const label = collapseMode === "exact" ? "Deduplicating AIRR records" : collapseMode === "fad" ? "Running FAD-compatible denoising" : collapseMode === "indel" ? "Running indel-aware denoising" : "Running conservative error-model denoising";
    const result = await operation(label, () => collapseMode === "exact" ? runtime.deduplicate(dedupKey) : runtime.denoise({
      mode: collapseMode,
      errorRate: denoiseErrorRate,
      alpha: denoiseAlpha,
      callResolution: denoiseResolution,
      ambiguity: denoiseAmbiguity,
      minimumParentCount,
      ambiguousPolicy: denoiseAmbiguousPolicy,
      fadNeighborThreshold,
      fadMethod,
      expectedZeroErrorFraction,
      maximumHammingDistance: maximumDenoiseDistance,
      maximumEditDistance,
      minimumIndelParentRatio,
      maxCandidatesPerVariant: denoiseCandidateCap,
    }, (processed, total) => setProgress({ processed, total })));
    if (result) {
      setDedup(result);
      setWorkingMask(null);
      setWorkingStages([]);
      setChmm(null);
      setChmmRun(null);
      setChimeraDetail(null);
      setLineages(null);
      setSelectedLineage(null);
      setLineageMultiplicity(new Map());
      setQueryHits([]);
      setExpanded(null);
    }
  }

  async function applyDedupFilter() {
    if (!dedup) return;
    const result = await operation("Applying deduplicated representatives to the downstream working set", () => runtime.applyDedupFilter());
    if (!result) return;
    setWorkingMask(result.mask);
    setWorkingStages([{
      id: "dedup",
      label: dedup.mode === "exact" ? `Exact collapse · ${dedup.key}` : dedup.mode === "fad" ? "FAD denoising" : dedup.mode === "indel" ? "Indel-aware denoising" : "Conservative denoising",
      input: dedup.inputRecords,
      retained: result.retained,
      discarded: dedup.collapsedRecords,
      detail: `${dedup.algorithm}. Collapsed abundance remains in duplicate_count and lineage weights.`,
    }]);
    setChmm(null);
    setChmmRun(null);
    setChimeraDetail(null);
    setLineages(null);
    setSelectedLineage(null);
    setLineageMultiplicity(new Map());
    clearAlignmentArtifacts();
    setQueryHits([]);
    setExpanded(null);
  }

  async function resetWorkingSet() {
    const result = await operation("Restoring the complete downstream working set", () => runtime.setActiveMask(null));
    if (!result) return;
    setWorkingMask(null);
    setWorkingStages([]);
    setChmm(null);
    setChmmRun(null);
    setChimeraDetail(null);
    setLineages(null);
    setSelectedLineage(null);
    setLineageMultiplicity(new Map());
    clearAlignmentArtifacts();
    setQueryHits([]);
    setExpanded(null);
  }

  async function downloadDeduplicated() {
    setBusy("Writing deduplicated AIRR table");
    setError("");
    try {
      const counts = await runtime.dedupCounts();
      const suffix = dedup?.mode === "exact" ? "deduplicated" : "denoised";
      await saveStream(`${baseName(inputName)}.${suffix}.airr.tsv`, "Collapsed AIRR rearrangement table with multiplicity", ".tsv", async (writer) => store.writeDeduplicatedAirr(counts, writer.write));
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
    }, workingStages.some((stage) => stage.id === "dedup")));
    if (result) {
      setLineages(result);
      setSelectedLineage(null);
      setLineageRows([]);
      setLineageMultiplicity(new Map());
      clearAlignmentArtifacts();
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
      const deduplicationApplied = workingStages.some((stage) => stage.id === "dedup");
      const counts = deduplicationApplied ? await runtime.dedupCounts() : null;
      const multiplicity = new Map(rows.map((row) => {
        const imported = Number(row.values.duplicate_count);
        const value = counts?.[row.record.ordinal] || (Number.isFinite(imported) && imported > 0 ? imported : 1);
        return [row.record.ordinal, Math.max(1, Math.floor(value))] as const;
      }));
      setSelectedLineage(summary);
      setLineageRows(rows);
      setLineageMultiplicity(multiplicity);
      setLineageTotal(members.total);
      clearAlignmentArtifacts();
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
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
      installAlignment(next, alignmentMethod === "quick" ? "AIRR-anchored alignment" : alignmentMethod === "codon" ? "Codon-aware Kalign 3.3.1" : "Nucleotide Kalign 3.3.1");
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function inferTree() {
    if (!alignment) return;
    setBusy("Running FastTree WASM and rooting on the N-masked germline");
    setTreeError("");
    const revision = alignmentRevisionRef.current;
    const alignmentSnapshot = alignment;
    const sourceSnapshot = alignmentSource;
    try {
      const execution = await runFastTree(alignmentSnapshot, treeModel, treeFast);
      if (revision !== alignmentRevisionRef.current) throw new Error("The alignment changed while FastTree was running. Run the tree again on the current alignment.");
      const rooted = rootOnOutgroup(parseNewick(execution.newick), GERMLINE_OUTGROUP);
      const canonicalRooted = canonicalizeTree(rooted);
      const stable = collapseShortInternalBranches(canonicalRooted);
      setTreeRun({
        ...execution,
        rootedNewick: `${serializeNewick(canonicalRooted)};`,
        stableNewick: `${serializeNewick(stable.root)};`,
        collapsedEdges: stable.collapsedEdges,
        collapseThreshold: stable.threshold,
        source: sourceSnapshot,
      });
      setTreeViewMode("rooted");
      window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (operationError) {
      setTreeError(operationError instanceof Error ? operationError.message : String(operationError));
      window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "center" }));
    } finally {
      setBusy("");
    }
  }

  function importEditedAlignment(text: string, source = "Alivibe-corrected alignment") {
    if (!alignment) throw new Error("Create a lineage alignment before importing a correction.");
    const inspected = validateCorrectedAlignment(alignment, text);
    installAlignment(inspected.fasta, source);
    setAlignmentEditorError("");
    setAlignmentEditorStatus(`Accepted the complete corrected alignment: ${inspected.rows.toLocaleString()} rows × ${inspected.columns.toLocaleString()} columns · fingerprint ${inspected.fingerprint}. FastTree will use these exact aligned rows.`);
  }

  function openAlivibeEditor() {
    if (!alignment) return;
    setAlignmentEditorError("");
    setAlignmentEditorStatus("Opening Alivibe…");
    const popup = window.open("https://murrellgroup.github.io/WebWidgets/alivibe.html", "swig-alivibe", "popup,width=1500,height=920");
    if (!popup) {
      setAlignmentEditorError("The browser blocked the Alivibe window. Allow pop-ups for this page and try again.");
      return;
    }
    alivibeWindowRef.current = popup;
    if (navigator.clipboard?.writeText) {
      void navigator.clipboard.writeText(alignment).catch(() => undefined);
    }
    let attempts = 0;
    const timer = window.setInterval(() => {
      attempts += 1;
      if (popup.closed) {
        window.clearInterval(timer);
        setAlignmentEditorStatus("Alivibe closed. Import its downloaded full-alignment FASTA; clipboard import is also accepted when it contains every complete row.");
        return;
      }
      try {
        const editor = popup as Window & {
          parseFasta?: (text: string) => void;
          getClipboardContent?: (preferSelection?: boolean) => string;
        };
        if (typeof editor.parseFasta !== "function") return;
        const controls = editor.document.getElementById("controls");
        if (!controls) return;
        editor.parseFasta(alignment);
        if (controls && !editor.document.getElementById("swig-return-control")) {
          const group = editor.document.createElement("div");
          group.id = "swig-return-control";
          group.className = "control-group";
          const label = editor.document.createElement("label");
          label.textContent = "Swig round trip";
          const button = editor.document.createElement("button");
          button.type = "button";
          button.textContent = "Return alignment to Swig";
          button.className = "active";
          button.onclick = () => {
            try {
              const corrected = editor.getClipboardContent?.(false) ?? "";
              importEditedAlignment(corrected);
              window.focus();
              popup.close();
            } catch (importError) {
              setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
              window.focus();
            }
          };
          group.append(label, button);
          controls.prepend(group);
        }
        setAlignmentEditorStatus("Alignment loaded locally in Alivibe. Edit it, then press Return alignment to Swig in Alivibe’s toolbar.");
        window.clearInterval(timer);
      } catch {
        if (attempts < 20) return;
        setAlignmentEditorStatus("Alivibe is on a different origin. The FASTA was copied locally: paste and edit it there, download the complete alignment, then use Import corrected FASTA. Alivibe’s Copy action may contain only a selection.");
        window.clearInterval(timer);
      }
    }, 250);
  }

  async function importFromAlivibe() {
    setAlignmentEditorError("");
    try {
      const editor = alivibeWindowRef.current as (Window & { getClipboardContent?: (preferSelection?: boolean) => string }) | null;
      let corrected = "";
      if (editor && !editor.closed) {
        try { corrected = editor.getClipboardContent?.(false) ?? ""; } catch { /* cross-origin fallback below */ }
      }
      if (!corrected && navigator.clipboard?.readText) corrected = await navigator.clipboard.readText();
      if (!corrected) throw new Error("No complete corrected FASTA was available. Download the full alignment from Alivibe, then use Import corrected FASTA.");
      importEditedAlignment(corrected);
    } catch (importError) {
      setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
    }
  }

  async function acceptEditedAlignment(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    try {
      importEditedAlignment(await file.text(), `Corrected alignment imported from ${file.name}`);
    } catch (importError) {
      setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
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
      const options: ChmmRunOptions = {
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
      };
      const inputMask = workingMask?.slice() ?? null;
      const result = await runChmmairra(store, msa, options, inputMask ?? undefined, (processed, total) => setProgress({ processed, total }));
      setChmm(result);
      setChmmRun({ msa, options, inputMask });
      setChmmFilterThreshold(chmmThreshold);
      setChimeraDetail(null);
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
  }

  async function applyChimeraFilter() {
    if (!chmm || !chmmRun) return;
    const next = new Uint8Array(store.count);
    let retained = 0;
    let unevaluatedRetained = 0;
    for (let ordinal = 0; ordinal < store.count; ordinal += 1) {
      if (chmmRun.inputMask && !chmmRun.inputMask[ordinal]) continue;
      const probability = chmm.probabilities[ordinal];
      const keep = Number.isFinite(probability) ? probability < chmmFilterThreshold : retainUnevaluated;
      if (!keep) continue;
      next[ordinal] = 1;
      retained += 1;
      if (!Number.isFinite(probability)) unevaluatedRetained += 1;
    }
    const result = await operation("Applying CHMMAIRRa exclusion to the downstream working set", () => runtime.setActiveMask(next));
    if (!result) return;
    setWorkingMask(next);
    setWorkingStages((stages) => [...stages.filter((stage) => stage.id !== "chimera"), {
      id: "chimera",
      label: `${chmm.segment} chimera posterior < ${chmmFilterThreshold}`,
      input: chmm.inputRecords,
      retained,
      discarded: chmm.inputRecords - retained,
      detail: retainUnevaluated ? `${unevaluatedRetained.toLocaleString()} unevaluated records retained.` : "Unevaluated records excluded.",
    }]);
    setLineages(null);
    setSelectedLineage(null);
    setLineageMultiplicity(new Map());
    clearAlignmentArtifacts();
    setTreeError("");
    setQueryHits([]);
    setExpanded(null);
  }

  async function openChimera(ordinal: number) {
    if (!chmmRun) {
      setError("Run CHMMAIRRa before opening a Viterbi parent reconstruction.");
      return;
    }
    setBusy("Running the selected record's detailed CHMMAIRRa Viterbi path");
    setError("");
    try {
      const detail = await runChmmairraDetail(store, chmmRun.msa, chmmRun.options, ordinal);
      setChimeraDetail(detail);
      window.requestAnimationFrame(() => chimeraDetailRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
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
    const hits = await operation(queryConstraintMode === "infer" ? "Assigning query V/J calls with SwiftIG, then searching" : "Searching the assigned repertoire", async () => {
      let searchQueries = queries;
      let queryConstraints: QueryConstraint[] | undefined;
      if (queryConstraintMode === "infer") {
        const inferred = await inferQueryAssignments(queries, queryTarget, references, minimumIdentity, strand, workers);
        setQueryInference(inferred);
        const usable = inferred.filter((assignment) => assignment.searchSequence && (queryV || assignment.vCall) && (queryJ || assignment.jCall));
        if (!usable.length) {
          throw new Error("SwiftIG could not infer both V and J for any query. Use complete rearranged nucleotide sequences, lower the main identity floor in a new run, or provide manual V/J overrides.");
        }
        searchQueries = usable.map((assignment) => assignment.searchSequence);
        queryConstraints = usable.map((assignment) => ({
          locus: queryLocus || assignment.locus || undefined,
          vCall: queryV || assignment.vCall || undefined,
          jCall: queryJ || assignment.jCall || undefined,
        }));
      } else {
        setQueryInference([]);
      }
      return runtime.query(searchQueries, {
        target: queryTarget,
        metric,
        identity: queryIdentity,
        maxResults: queryLimit,
        locus: queryConstraintMode === "manual" ? queryLocus || undefined : undefined,
        vCall: queryConstraintMode === "manual" ? queryV || undefined : undefined,
        jCall: queryConstraintMode === "manual" ? queryJ || undefined : undefined,
        queryConstraints,
        callResolution: resolution,
        ambiguity,
        productiveOnly,
      });
    });
    if (!hits) return;
    setQueryHits(hits);
    setExpanded(null);
    const indexRows = await store.indexRecords([...new Set(hits.map((hit) => hit.ordinal))]);
    setQueryRecords(new Map(indexRows.map((record) => [record.ordinal, record])));
  }

  async function previewQueryInference() {
    const queries = parseQueries(queryText);
    if (!queries.length) {
      setError("Paste at least one complete rearranged nucleotide sequence or FASTA record.");
      return;
    }
    setBusy("Assigning seed V/J calls with the main SwiftIG configuration");
    setError("");
    try {
      setQueryInference(await inferQueryAssignments(queries, queryTarget, references, minimumIdentity, strand, workers));
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      setBusy("");
    }
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

  const workingCount = workingStages.length ? workingStages[workingStages.length - 1].retained : store.count;
  const chmmFilterThresholdValid = Number.isFinite(chmmFilterThreshold) && chmmFilterThreshold >= 0 && chmmFilterThreshold <= 1;
  const chmmFilterPreview = useMemo(() => {
    if (!chmm || !chmmRun || !Number.isFinite(chmmFilterThreshold) || chmmFilterThreshold < 0 || chmmFilterThreshold > 1) return null;
    let retained = 0;
    let excluded = 0;
    let unevaluated = 0;
    for (let ordinal = 0; ordinal < store.count; ordinal += 1) {
      if (chmmRun.inputMask && !chmmRun.inputMask[ordinal]) continue;
      const probability = chmm.probabilities[ordinal];
      if (!Number.isFinite(probability)) {
        unevaluated += 1;
        if (retainUnevaluated) retained += 1;
        else excluded += 1;
      } else if (probability >= chmmFilterThreshold) excluded += 1;
      else retained += 1;
    }
    return { retained, excluded, unevaluated };
  }, [chmm, chmmFilterThreshold, chmmRun, retainUnevaluated, store.count]);
  const displayedQueryOrdinals = expanded?.ordinals ?? queryHits.map((hit) => hit.ordinal);
  const displayedQueryRows = [...new Set(displayedQueryOrdinals)].slice(0, 500).flatMap((ordinal) => queryRecords.get(ordinal) ?? []);
  const vChart = lineages?.vUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];
  const jChart = lineages?.jUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];

  return <section className="post-analysis-shell">
    <header className="post-analysis-heading"><div><span className="section-kicker">Post-assignment methods</span><h2>Repertoire structure and targeted phylogenetics</h2><p>Exact collapse or denoising and chimera exclusion modify an explicit cumulative working set. CHMMAIRRa, lineage assignment, repertoire querying, and expansion consume that set; alignment and tree inference consume the selected lineage.</p></div><div className="local-method-note"><span>Execution</span><strong>Browser-local</strong><small>Input, germlines, and results are not submitted to an analysis server.</small></div></header>

    <div className="post-method-map"><article><b>01</b><span>Collapse</span><strong>Exact deduplication or denoising</strong></article><article><b>02</b><span>QC</span><strong>Optional CHMMAIRRa</strong></article><article><b>03</b><span>Repertoire</span><strong>Assign lineages</strong></article><article><b>04</b><span>Targeted</span><strong>Query + expand</strong></article><article><b>05</b><span>On demand</span><strong>Align + infer tree</strong></article></div>

    <section className="working-set-panel" aria-label="Downstream working-set pipeline">
      <header><div><span className="section-kicker">Cumulative downstream filter</span><h3>{workingCount.toLocaleString()} of {store.count.toLocaleString()} records active</h3><p>Nothing is deleted from the AIRR result. Applying a stage excludes rows from later computation; resetting restores the complete input.</p></div><button type="button" disabled={Boolean(busy) || !workingStages.length} onClick={() => void resetWorkingSet()}>Reset to all records</button></header>
      <div className="working-set-flow"><article className="source"><span>Assigned input</span><strong>{store.count.toLocaleString()}</strong><small>records</small></article>{workingStages.map((stage) => <article key={stage.id} className={stage.id}><span>{stage.label}</span><strong>{stage.retained.toLocaleString()}</strong><small>retained · {stage.discarded.toLocaleString()} excluded at this step</small><p>{stage.detail}</p></article>)}{!workingStages.length && <article className="pass-through"><span>No applied filter</span><strong>All records</strong><small>pass downstream</small></article>}<article className="consumers"><span>Current consumers</span><strong>CHMMAIRRa · lineages · query</strong><small>alignment/tree follow the selected lineage</small></article></div>
    </section>

    {busy && <div className="post-progress" role="status"><div><span>{busy}</span><strong>{progress.total ? `${Math.min(100, progress.processed / progress.total * 100).toFixed(1)}%` : "working"}</strong></div><progress max={Math.max(1, progress.total)} value={progress.processed} /><small>{progress.processed.toLocaleString()} / {progress.total.toLocaleString()} AIRR records indexed or scanned</small></div>}
    {error && <div className="post-error" role="alert"><strong>Post-analysis stopped</strong><span>{error}</span><button type="button" onClick={() => setError("")}>Dismiss</button></div>}

    <section className="post-module dedup-module">
      <header><div className="module-number">01</div><div><span className="section-kicker">Abundance preservation</span><h3>Collapse exact duplicates or denoise read errors</h3><p>Select one method explicitly. Every retained representative carries the sum of its source multiplicities in <code>duplicate_count</code>; lineage abundance and phylogeny bubbles use that value.</p></div><a href="https://academic.oup.com/nar/article/47/18/e104/5550323" target="_blank" rel="noreferrer">FAD paper ↗</a></header>
      <div className="collapse-mode-grid" role="radiogroup" aria-label="Collapse method">
        <button type="button" role="radio" aria-checked={collapseMode === "exact"} className={collapseMode === "exact" ? "selected" : ""} onClick={() => { setCollapseMode("exact"); setDedup(null); }}><b>A</b><span><strong>Exact deduplication</strong><small>Collapse identical keys only. No error model.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "fad"} className={collapseMode === "fad" ? "selected" : ""} onClick={() => { setCollapseMode("fad"); setDenoiseAmbiguousPolicy("exclude"); setDedup(null); }}><b>B</b><span><strong>FAD-compatible denoising</strong><small>Published 6-mer distance and abundance/Poisson rule.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "conservative"} className={collapseMode === "conservative" ? "selected" : ""} onClick={() => { setCollapseMode("conservative"); setDenoiseAmbiguousPolicy("retain"); setDedup(null); }}><b>C</b><span><strong>Exact-neighbor error model</strong><small>Experimental; conservative, indexed Hamming candidates.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "indel"} className={collapseMode === "indel" ? "selected" : ""} onClick={() => { setCollapseMode("indel"); setDenoiseAmbiguousPolicy("retain"); setDedup(null); }}><b>D</b><span><strong>Indel-aware error model</strong><small>Complete 1–2 edit index; abundance-directed indel collapse.</small></span></button>
      </div>
      {collapseMode === "exact" ? <>
        <div className="module-controls">
          <label><span>Identity key</span><select value={dedupKey} onChange={(event) => setDedupKey(event.target.value as DedupKey)}><option value="sequence">Full input sequence</option><option value="trimmed">VDJ-aligned sequence</option><option value="cdr3">Locus + CDR3 nucleotide</option><option value="rearrangement">Locus + V/J calls + CDR3</option></select></label>
          <button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runDedup()}>Run exact deduplication</button>
        </div>
        {(dedupKey === "sequence" || dedupKey === "trimmed") && <p className="scientific-note"><span>i</span>Sequence-key modes compare normalized length plus a 128-bit fingerprint so complete sequence payloads do not remain in memory. Existing <code>duplicate_count</code> values are summed rather than reset.</p>}
      </> : <div className="denoise-config">
        <div className="control-grid three">
          <label><span>V/J call level</span><select value={denoiseResolution} onChange={(event) => setDenoiseResolution(event.target.value as CallResolution)}><option value="allele">Allele</option><option value="gene">Gene</option></select></label>
          <label><span>Uncertain assignments</span><select value={denoiseAmbiguity} onChange={(event) => setDenoiseAmbiguity(event.target.value as "top" | "strict")}><option value="strict">Exact sorted call set</option><option value="top">Top V and J calls</option></select></label>
          <label><span>Ambiguous N bases</span><select value={denoiseAmbiguousPolicy} onChange={(event) => setDenoiseAmbiguousPolicy(event.target.value as "exclude" | "retain")}><option value="exclude">Exclude (published FAD behavior)</option><option value="retain">Keep as exact-only representatives</option></select></label>
          <label><span>Per-base error rate</span><CommitNumberInput min="0.000001" max="0.2" step="0.00001" value={denoiseErrorRate} onCommit={setDenoiseErrorRate} /></label>
          <label><span>Significance α</span><CommitNumberInput min="0.000001" max="0.5" step="0.001" value={denoiseAlpha} onCommit={setDenoiseAlpha} /></label>
          <label><span>Minimum parent abundance</span><CommitNumberInput min="1" max="1000000" step="1" value={minimumParentCount} onCommit={setMinimumParentCount} /></label>
        </div>
        {collapseMode === "fad" ? <div className="control-grid three denoise-specific">
          <label><span>FAD decision</span><select value={fadMethod} onChange={(event) => setFadMethod(Number(event.target.value) as 1 | 2)}><option value={2}>Method 2 · Poisson</option><option value={1}>Method 1 · abundance only</option></select></label>
          <label><span>Corrected 6-mer radius</span><CommitNumberInput min="0" max="5" step="0.25" value={fadNeighborThreshold} onCommit={setFadNeighborThreshold} /></label>
          <label><span>Expected zero-error fraction</span><CommitNumberInput min="0.0001" max="1" step="0.01" value={expectedZeroErrorFraction} onCommit={setExpectedZeroErrorFraction} /></label>
        </div> : collapseMode === "conservative" ? <div className="control-grid three denoise-specific"><label><span>Maximum Hamming errors</span><select value={maximumDenoiseDistance} onChange={(event) => setMaximumDenoiseDistance(Number(event.target.value))}><option value={1}>1 substitution</option><option value={2}>2 substitutions</option><option value={3}>3 substitutions</option></select></label></div> : <div className="control-grid three denoise-specific">
          <label><span>Maximum edit distance</span><select value={maximumEditDistance} onChange={(event) => setMaximumEditDistance(Number(event.target.value))}><option value={1}>1 edit</option><option value={2}>2 edits</option></select></label>
          <label><span>Minimum indel parent:child ratio</span><CommitNumberInput min="1.01" max="1000000" step="0.1" value={minimumIndelParentRatio} onCommit={setMinimumIndelParentRatio} /></label>
        </div>}
        <details className="post-advanced"><summary>Performance guard</summary><div className="control-grid three"><label><span>Candidate cap / variant</span><CommitNumberInput min="100" max="1000000" step="1000" value={denoiseCandidateCap} onCommit={setDenoiseCandidateCap} /></label></div></details>
        <div className="algorithm-note"><strong>{collapseMode === "fad" ? "FAD-compatible indexed implementation" : collapseMode === "indel" ? "Complete bounded-edit implementation" : "Conservative exact-neighbor implementation"}</strong><span>{collapseMode === "fad" ? "Partition locus/V/J → exact dereplication → complete corrected-6-mer radius index → published template test → exact nearest accepted centroid within each partition." : collapseMode === "indel" ? "Partition locus/V/J → exact dereplication → complete length-aware d+1 segment join → exact banded edit profiling. Low-abundance indel paths use the parent:child ratio; substitution-only paths retain method C’s Poisson test. No distant-centroid assignment is performed." : "Partition locus/V/J → exact dereplication → d+1 block index → exact Hamming verification → sequence-specific Poisson error test. Isolated singletons are retained; no read is assigned to a distant centroid."}</span></div>
        {collapseMode === "indel" && <p className="scientific-note warning"><span>!</span>Method D deliberately treats a one- or two-edit path containing an indel as likely technical error once the parent:child abundance ratio is met. True biological length variants can therefore be merged; set a stricter ratio or use method C when that assumption is inappropriate.</p>}
        <p className="scientific-note"><span>i</span>Trimmed VDJ sequences are streamed from the AIRR store into a two-bit packed arena. Temporary neighbor profiles are released one V/J partition at a time. The 0.00473 default reproduces the linked workflow’s MiSeq median error-rate setting; it should be changed for a different assay.</p>
        <button className="post-primary denoise-run" type="button" disabled={Boolean(busy)} onClick={() => void runDedup()}>{collapseMode === "fad" ? "Run FAD-compatible denoising" : collapseMode === "indel" ? "Run indel-aware denoising" : "Run exact-neighbor denoising"}</button>
      </div>}
      {dedup && <div className="module-result"><div className="post-stat-grid"><article><span>Input rows / abundance</span><strong>{dedup.inputRecords.toLocaleString()} / {dedup.inputAbundance.toLocaleString()}</strong></article><article><span>Retained representatives</span><strong>{dedup.uniqueRecords.toLocaleString()}</strong></article><article><span>Collapsed or excluded rows</span><strong>{dedup.collapsedRecords.toLocaleString()}</strong></article><article><span>Largest multiplicity</span><strong>{(dedup.largestGroups[0]?.count ?? 1).toLocaleString()}</strong></article>{dedup.mode !== "exact" && <><article><span>V/J partitions</span><strong>{dedup.partitions.toLocaleString()}</strong></article><article><span>Verified candidates</span><strong>{dedup.candidateComparisons.toLocaleString()}</strong></article></>}{dedup.mode === "indel" && <><article><span>Indel variants merged</span><strong>{dedup.indelMergedVariants.toLocaleString()}</strong></article><article><span>Substitution variants merged</span><strong>{dedup.substitutionMergedVariants.toLocaleString()}</strong></article></>}</div><div className="denoise-provenance"><strong>{dedup.algorithm}</strong><span>{dedup.mode === "exact" ? `Key: ${dedup.key}` : `${denoiseResolution}-level V/J partitioning · ${denoiseAmbiguity} assignment policy`}</span></div>{dedup.warnings.map((warning) => <div key={warning} className="scientific-note warning"><span>!</span><p>{warning}</p></div>)}<div className="filter-commit"><div><span>Downstream action</span><strong>Retain representatives; preserve collapsed abundance as counts</strong><p>This action changes the working set from {dedup.inputRecords.toLocaleString()} to {dedup.uniqueRecords.toLocaleString()} rows and invalidates downstream results.</p></div><button className="post-primary" type="button" disabled={Boolean(busy) || workingStages.some((stage) => stage.id === "dedup")} onClick={() => void applyDedupFilter()}>{workingStages.some((stage) => stage.id === "dedup") ? "Applied downstream" : `Apply ${dedup.uniqueRecords.toLocaleString()} representatives`}</button></div><div className="result-actions"><button type="button" onClick={() => void downloadDeduplicated()}>Download collapsed AIRR + multiplicity</button></div></div>}
    </section>

    <section className="post-module chmm-module">
      <header><div className="module-number amber">02</div><div><span className="section-kicker">Optional PCR-chimera model</span><h3>CHMMAIRRa after V(D)J assignment</h3><p>The browser port threads each AIRR local V or J alignment onto a reference MSA, then evaluates the CHMMera posterior. V is the manuscript default; D is not modeled.</p></div><a href="https://github.com/MurrellGroup/CHMMAIRRa.jl" target="_blank" rel="noreferrer">Method source ↗</a></header>
      <div className="chmm-grid">
        <div className="chmm-config">
          <div className="control-grid three"><label><span>Segment</span><select value={chmmSegment} onChange={(event) => { setChmmSegment(event.target.value as ChmmSegment); setPreparedMsa(""); setChmm(null); setChmmRun(null); setChimeraDetail(null); }}><option value="V">V (recommended)</option><option value="J">J (optional)</option></select></label><label><span>Model</span><select value={chmmMethod} onChange={(event) => setChmmMethod(event.target.value as "BW" | "DB")}><option value="BW">Baum–Welch · IG default</option><option value="DB">Discretized Bayesian · TCR default</option></select></label><label><span>Posterior threshold</span><CommitNumberInput min="0" max="1" step="0.01" value={chmmThreshold} onCommit={setChmmThreshold} /></label></div>
          <fieldset className="msa-source"><legend>Reference multiple-sequence alignment</legend><label className={chmmSource === "selected" ? "selected" : ""}><input type="radio" checked={chmmSource === "selected"} onChange={() => setChmmSource("selected")} /><span><strong>Build from this run’s {chmmSegment} references</strong><small>Kalign 3.3.1 WASM; preserves selected IMGT/KI/uploaded composition.</small></span></label><label className={chmmSource === "upload" ? "selected" : ""}><input type="radio" checked={chmmSource === "upload"} onChange={() => setChmmSource("upload")} /><span><strong>Use an aligned FASTA MSA</strong><small>{uploadedMsaName || "Every record must have equal aligned length and names matching AIRR calls."}</small></span><input className="file-inline" type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event) => void acceptMsa(event)} /></label></fieldset>
          <details className="post-advanced"><summary>Model parameters</summary><div className="control-grid three"><label><span>Chimera prior</span><CommitNumberInput min="0.00001" max="0.5" step="0.01" value={chmmPrior} onCommit={setChmmPrior} /></label><label><span>Minimum DFR</span><CommitNumberInput min="0" max="100" step="1" value={chmmMinDfr} onCommit={setChmmMinDfr} /></label><label><span>DB mutation rates</span><CommitTextInput value={mutationRates} onCommit={setMutationRates} /></label><label className="check-line"><input type="checkbox" checked={chmmDetailed} onChange={(event) => setChmmDetailed(event.target.checked)} /><span>Precompute breakpoint labels during the repertoire scan (the full path remains on-demand)</span></label></div></details>
          <div className="scientific-note warning"><span>!</span><p>Reference completeness matters: an absent true V/J allele can produce a false switch signal. Uploaded MSAs are never silently supplemented. Low-DFR records are reported as unevaluated rather than forced through the model.</p></div>
          <div className="result-actions"><button className="post-primary amber" type="button" disabled={Boolean(busy) || (chmmSource === "upload" && !uploadedMsa)} onClick={() => void runChmmAnalysis()}>Run CHMMAIRRa on {workingCount.toLocaleString()}</button>{preparedMsa && <button type="button" onClick={() => downloadText(preparedMsa, `${baseName(inputName)}.${chmmSegment.toLowerCase()}-reference-msa.fasta`)}>Download reference MSA</button>}</div>
        </div>
        <div className="chmm-result-panel">
          {chmm ? <>
            <div className="post-stat-grid compact"><article><span>Working-set input</span><strong>{chmm.inputRecords.toLocaleString()}</strong></article><article><span>Evaluated</span><strong>{chmm.evaluated.toLocaleString()}</strong></article><article><span>Posterior ≥ {chmm.threshold}</span><strong>{chmm.flagged.toLocaleString()}</strong></article><article><span>Below DFR</span><strong>{chmm.lowDfr.toLocaleString()}</strong></article><article><span>Missing reference</span><strong>{chmm.missingReference.toLocaleString()}</strong></article>{chmm.upstreamExcluded > 0 && <article><span>Excluded upstream</span><strong>{chmm.upstreamExcluded.toLocaleString()}</strong></article>}</div>
            <div className="chmm-filter-commit">
              <header><div><span className="section-kicker">Explicit downstream exclusion</span><h4>Filter chimeras from later steps</h4><p>Change the threshold without rerunning the HMM, inspect the resulting counts, then apply the filter.</p></div></header>
              <div className="chmm-filter-controls"><label><span>Exclude posterior ≥</span><CommitNumberInput min="0" max="1" step="0.01" value={chmmFilterThreshold} onCommit={setChmmFilterThreshold} /></label><label className="check-line"><input type="checkbox" checked={retainUnevaluated} onChange={(event) => setRetainUnevaluated(event.target.checked)} /><span>Retain unevaluated records</span></label></div>
              {chmmFilterPreview && <div className="filter-preview"><div><span>Retained downstream</span><strong>{chmmFilterPreview.retained.toLocaleString()}</strong></div><div><span>Excluded as chimera</span><strong>{chmmFilterPreview.excluded.toLocaleString()}</strong></div><div><span>Unevaluated in input</span><strong>{chmmFilterPreview.unevaluated.toLocaleString()}</strong></div></div>}
              {!chmmFilterThresholdValid && <div className="inline-field-error">Enter a posterior threshold from 0 to 1.</div>}
              <button className="post-primary amber" type="button" disabled={Boolean(busy) || !chmmFilterThresholdValid} onClick={() => void applyChimeraFilter()}>Apply chimera filter downstream</button>
            </div>
            <BarChart title={`${chmm.segment} chimera posterior`} subtitle="Evaluated rearrangements by posterior interval" data={chmm.histogram.map((item) => ({ label: item.label, value: item.count }))} color="#d49a19" name={`${baseName(inputName)}.chmmairra-posterior.svg`} />
            {chmm.top.length > 0 && <div className="chmm-top-list"><header><strong>Highest posteriors</strong><span>Click to run the detailed Viterbi parent view</span></header>{chmm.top.slice(0, 12).map((record) => <button type="button" key={record.ordinal} onClick={() => void openChimera(record.ordinal)}><b>#{(record.ordinal + 1).toLocaleString()}</b><span>{(record.probability * 100).toFixed(2)}%</span><small>DFR {record.dfr}{record.recombinations.length ? ` · ${record.recombinations.map((event) => `${event.left}→${event.right}@${event.position}`).join("; ")}` : " · Viterbi path on click"}</small></button>)}</div>}
            <button type="button" onClick={() => void downloadChmm()}>Download CHMMAIRRa TSV</button>
          </> : <div className="method-placeholder"><span>HMM</span><h4>No CHMMAIRRa run</h4><p>{isTcr ? "TCR mode defaults to a fixed 0.005 mutation-rate state (DB)." : "IG mode defaults to per-reference Baum–Welch mutation estimates."}</p></div>}
        </div>
      </div>
      {chimeraDetail && <div ref={chimeraDetailRef}><ChimeraHighlighter detail={chimeraDetail} name={`${baseName(inputName)}.record-${chimeraDetail.ordinal + 1}.chmmairra-viterbi.svg`} onInspect={onInspect} /></div>}
    </section>

    <section className="post-module lineage-module">
      <header><div className="module-number">03</div><div><span className="section-kicker">Repertoire-scale clonal grouping</span><h3>Assign lineages from CDR3 nucleotide distance</h3><p>Default: same locus, overlapping V/J gene assignments, exact CDR3 nucleotide length, and single-linkage at ≥85% identity. The threshold is a starting point and remains dataset-adjustable.</p></div><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5340603/" target="_blank" rel="noreferrer">Clonal threshold literature ↗</a></header>
      <div className="lineage-config"><div className="control-grid five"><label><span>CDR3 identity</span><div className="range-number"><input type="range" min="0.7" max="1" step="0.01" value={identity} onChange={(event) => setIdentity(Number(event.target.value))} /><b>{Math.round(identity * 100)}%</b></div></label><label><span>Call level</span><select value={resolution} onChange={(event) => setResolution(event.target.value as CallResolution)}><option value="gene">Gene</option><option value="allele">Allele</option></select></label><label><span>Ambiguous calls</span><select value={ambiguity} onChange={(event) => setAmbiguity(event.target.value as AmbiguityPolicy)}><option value="overlap">Any assignment overlaps</option><option value="top">Top call only</option><option value="strict">Exact call sets</option></select></label><label><span>Candidate cap / record</span><CommitNumberInput min="100" max="1000000" step="1000" value={candidateCap} onCommit={setCandidateCap} /></label><label className="check-line"><input type="checkbox" checked={productiveOnly} onChange={(event) => setProductiveOnly(event.target.checked)} /><span>Productive only</span></label></div><div className="algorithm-note"><strong>Exact accelerated single-linkage</strong><span>Partition by locus → V/J calls → CDR3 length → d+1 exact blocks → verify normalized Hamming distance → union-find components.</span></div><div className="current-step-input"><span>Input inherited from applied filters</span><strong>{workingCount.toLocaleString()} active records</strong><small>{workingStages.length ? workingStages.map((stage) => stage.label).join(" → ") : "No upstream exclusion applied"}</small></div><button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runLineages()}>Assign lineages on current set</button></div>
      {lineages && <div className="lineage-results"><div className="post-stat-grid"><article><span>Lineages</span><strong>{lineages.lineageCount.toLocaleString()}</strong></article><article><span>Assigned records</span><strong>{lineages.assignedRecords.toLocaleString()}</strong></article><article><span>Largest lineage</span><strong>{(lineages.summaries[0]?.abundance ?? 0).toLocaleString()}</strong></article><article><span>Exact comparisons</span><strong>{lineages.candidateComparisons.toLocaleString()}</strong></article></div>{lineages.truncatedCandidates > 0 && <div className="scientific-note warning"><span>!</span><p>{lineages.truncatedCandidates.toLocaleString()} records reached the candidate cap. Increase it and rerun before treating components as complete.</p></div>}<div className="result-actions"><button type="button" onClick={() => void downloadLineages()}>Download AIRR + clone_id</button></div><div className="chart-customizer"><label><span>Gene chart metric</span><select value={geneMetric} onChange={(event) => setGeneMetric(event.target.value as "abundance" | "lineages")}><option value="abundance">Sequence abundance</option><option value="lineages">Lineage count</option></select></label><label><span>Top genes</span><CommitNumberInput min="5" max="24" value={topGenes} onCommit={setTopGenes} /></label><label><span>Figure color</span><input type="color" value={chartColor} onChange={(event) => setChartColor(event.target.value)} /></label></div><div className="post-chart-grid"><BarChart title="Lineage abundance distribution" subtitle="Lineage count in each abundance interval" data={lineages.sizeHistogram.map((item) => ({ label: item.label, value: item.count }))} color={chartColor} name={`${baseName(inputName)}.lineage-size-distribution.svg`} /><BarChart title="Largest lineages" subtitle="Abundance retained after deduplication" data={lineages.summaries.slice(0, 20).map((item) => ({ label: `Lineage ${item.id}`, value: item.abundance }))} color={chartColor} name={`${baseName(inputName)}.largest-lineages.svg`} /><BarChart title="V germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={vChart} color={chartColor} name={`${baseName(inputName)}.lineage-v-use.svg`} /><BarChart title="J germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={jChart} color={chartColor} name={`${baseName(inputName)}.lineage-j-use.svg`} /></div><div className="lineage-table-wrap"><table><thead><tr><th>Lineage</th><th>Abundance</th><th>Unique</th><th>Locus</th><th>V calls</th><th>J calls</th><th>CDR3 nt</th><th /></tr></thead><tbody>{lineages.summaries.slice(0, 250).map((summary) => <tr key={summary.id} className={selectedLineage?.id === summary.id ? "selected" : ""} onClick={() => void openLineage(summary)}><td><strong>{summary.id}</strong></td><td>{summary.abundance.toLocaleString()}</td><td>{summary.uniqueMembers.toLocaleString()}</td><td>{summary.locus}</td><td>{summary.vCalls.join(", ")}</td><td>{summary.jCalls.join(", ")}</td><td>{summary.cdr3Length} nt</td><td><button type="button">Open →</button></td></tr>)}</tbody></table></div></div>}
    </section>

    <section className="post-module query-module">
      <header><div className="module-number dark">04</div><div><span className="section-kicker">Targeted repertoire search</span><h3>Query sequences, then expand the matched set</h3><p>Paste one or more sequences or FASTA records. V/J constraints can be supplied directly or inferred per seed with the same SwiftIG references and assignment parameters as the main analysis.</p></div></header>
      <div className="query-layout">
        <div className="query-input">
          <label><span>Query sequence(s) or FASTA</span><CommitTextarea value={queryText} onCommit={setQueryText} placeholder=">seed_1\nTGTGCGAGAGAT…\n>seed_2\nTGTGCGAGGGAT…" /></label>
          <div className="control-grid three">
            <label><span>Target</span><select value={queryTarget} onChange={(event) => setQueryTarget(event.target.value as QueryTarget)}><option value="cdr3_nt">CDR3 nucleotide</option><option value="cdr3_aa">CDR3 amino acid</option><option value="trimmed">VDJ-aligned k-mer sketch</option></select></label>
            <label><span>Metric</span><select disabled={queryTarget === "trimmed"} value={queryTarget === "trimmed" ? "sketch" : queryMetric} onChange={(event) => setQueryMetric(event.target.value as QueryMetric)}><option value="exact">Exact</option><option value="substring">Substring</option><option value="hamming">Hamming</option><option value="edit">Banded edit distance</option><option value="sketch">MinHash k-mer estimate</option></select></label>
            <label><span>Identity / similarity</span><CommitNumberInput min="0" max="1" step="0.01" value={queryIdentity} onCommit={setQueryIdentity} /></label>
            <label><span>V/J constraint source</span><select value={queryConstraintMode} onChange={(event) => { setQueryConstraintMode(event.target.value as "manual" | "infer"); setQueryInference([]); }}><option value="manual">Manual / unconstrained</option><option value="infer">Infer per query with SwiftIG</option></select></label>
            <label><span>{queryConstraintMode === "infer" ? "Locus override" : "Locus constraint"}</span><select value={queryLocus} onChange={(event) => setQueryLocus(event.target.value)}><option value="">{queryConstraintMode === "infer" ? "Use each inferred locus" : "Any locus"}</option>{loci.map((item) => <option value={item.value} key={item.value}>{item.value}</option>)}</select></label>
            <label><span>{queryConstraintMode === "infer" ? "V override" : "V constraint"}</span><CommitTextInput value={queryV} onCommit={setQueryV} placeholder={queryConstraintMode === "infer" ? "otherwise inferred per seed" : "optional IGHV…"} /></label>
            <label><span>{queryConstraintMode === "infer" ? "J override" : "J constraint"}</span><CommitTextInput value={queryJ} onCommit={setQueryJ} placeholder={queryConstraintMode === "infer" ? "otherwise inferred per seed" : "optional IGHJ…"} /></label>
            <label><span>Maximum initial hits</span><CommitNumberInput min="1" max="10000" value={queryLimit} onCommit={setQueryLimit} /></label>
          </div>
          {queryInference.length > 0 && <div className="query-inference"><header><strong>SwiftIG seed assignments</strong><span>{queryInference.filter((item) => item.assigned).length} / {queryInference.length} with both V and J</span></header>{queryInference.slice(0, 20).map((assignment) => <div className={assignment.assigned ? "assigned" : "unassigned"} key={assignment.queryIndex}><b>Seed {assignment.queryIndex + 1}</b><span>{assignment.locus || "locus —"}</span><code>{assignment.vCall || "V —"}</code><code>{assignment.jCall || "J —"}</code><small>{assignment.searchSequence ? `${assignment.searchSequence.length} ${queryTarget === "cdr3_aa" ? "aa" : "nt"} search target` : "no target derived"}</small></div>)}</div>}
          <div className="current-step-input"><span>Search universe inherited from applied filters</span><strong>{workingCount.toLocaleString()} active records</strong><small>Initial search and every single-linkage expansion both use this same working set.</small></div>
          <div className="result-actions">{queryConstraintMode === "infer" && <button type="button" disabled={Boolean(busy)} onClick={() => void previewQueryInference()}>Preview inferred V/J</button>}<button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void runQuery()}>{queryConstraintMode === "infer" ? "Assign seeds + search" : "Search repertoire"}</button><button type="button" disabled={Boolean(busy) || !queryHits.length} onClick={() => void expandMatches()}>Single-linkage expand set</button></div>
          <p className="scientific-note"><span>i</span>{queryConstraintMode === "infer" ? `Input is interpreted as rearranged nucleotide sequence. SwiftIG uses this run’s composed references, ${Math.round(minimumIdentity * 100)}% identity floor, and ${strand === 0 ? "both strands" : strand === 1 ? "plus strand" : "minus strand"}; ambiguous V/J calls remain ambiguity-aware and are applied per seed. Non-empty override fields replace the inferred value.` : queryTarget === "trimmed" ? "VDJ searches lazily build a packed index of eight independent 7-mer MinHash values per record. Scores are approximate and intended for candidate retrieval." : "Hamming search requires equal length. Edit distance is banded by the selected identity. Expansion always uses equal-length CDR3 nucleotide Hamming edges."}</p>
        </div>
        <div className="query-results"><header><div><span className="section-kicker">Matched set</span><h4>{expanded ? `${expanded.ordinals.length.toLocaleString()} expanded records` : `${queryHits.length.toLocaleString()} initial hits`}</h4><p>{expanded ? `${expanded.comparisons.toLocaleString()} exact edge checks${expanded.capped ? " · result cap reached" : " · fixed point reached"}` : "Ranked by similarity, then distance."}</p></div></header><div className="query-result-list">{displayedQueryRows.slice(0, 100).map((record) => { const hit = queryHits.find((value) => value.ordinal === record.ordinal); return <button type="button" key={record.ordinal} onClick={() => onInspect(record.ordinal)}><span><strong>{record.sequenceId}</strong><small>{record.locus} · {record.vCall || "V—"} · {record.jCall || "J—"}</small></span><code>{record.cdr3 || record.cdr3Aa || "—"}</code><b>{hit ? `${(hit.score * 100).toFixed(1)}%` : "expanded"}</b></button>; })}{!displayedQueryRows.length && <div className="method-placeholder small"><span>⌕</span><h4>No query results</h4><p>Provide a sequence and search the assigned repertoire.</p></div>}</div></div>
      </div>
    </section>

    {selectedLineage && <section ref={workbenchRef} className="post-module lineage-workbench" tabIndex={-1}>
      <header><div className="module-number dark">05</div><div><span className="section-kicker">Selected lineage {selectedLineage.id}</span><h3>Alignment and rooted phylogeny</h3><p>{lineageTotal.toLocaleString()} active rows · abundance {selectedLineage.abundance.toLocaleString()} · {selectedLineage.locus} · {selectedLineage.vCalls.join(", ")} / {selectedLineage.jCalls.join(", ")}</p></div><button type="button" onClick={() => { setSelectedLineage(null); setLineageRows([]); setLineageMultiplicity(new Map()); }}>Close lineage</button></header>
      <div className="member-strip"><div><strong>{lineageRows.length.toLocaleString()}</strong><span>active rows loaded</span><small>{lineageTotal > lineageRows.length ? `first ${lineageRows.length}; alignment remains on-demand` : "complete selected lineage working set"}</small></div><div className="member-pills">{lineageRows.slice(0, 18).map((row) => <button type="button" key={row.record.ordinal} onClick={() => onInspect(row.record.ordinal)}>{row.record.sequenceId}</button>)}</div></div>
      <div className="alignment-controls"><label><span>Alignment method</span><select value={alignmentMethod} onChange={(event) => setAlignmentMethod(event.target.value as "quick" | "kalign" | "codon")}><option value="quick">AIRR-anchored quick view</option><option value="kalign">Kalign 3.3.1 WASM · nucleotide</option><option value="codon">Kalign 3.3.1 WASM · codon-aware</option></select></label><label><span>Maximum sequences</span><CommitNumberInput min="2" max="500" value={alignmentLimit} onCommit={setAlignmentLimit} /></label><button className="post-primary" type="button" disabled={Boolean(busy)} onClick={() => void runAlignment()}>Align selected lineage</button>{alignment && <button type="button" onClick={() => downloadText(alignment, `${baseName(inputName)}.lineage-${selectedLineage.id}.alignment.fasta`)}>Alignment FASTA ↓</button>}</div>
      {alignment && <>
        <div className="alignment-editor-transfer"><div><span className="section-kicker">Manual correction</span><h4>Round trip through Alivibe</h4><p>Swig keeps its internal aligners. A returned alignment is accepted only if it contains every original row and every ungapped nucleotide sequence is unchanged; selected or truncated exports are rejected.</p></div><div className="result-actions"><button type="button" onClick={openAlivibeEditor}>Open + load in Alivibe ↗</button><button type="button" onClick={() => void importFromAlivibe()}>Import full alignment</button><label className="file-button">Import corrected FASTA<input type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event) => void acceptEditedAlignment(event)} /></label></div>{alignmentEditorStatus && <p className="editor-status">{alignmentEditorStatus}</p>}{alignmentEditorError && <div className="inline-method-error" role="alert">{alignmentEditorError}</div>}<small>For the cross-origin fallback, use Alivibe’s full-alignment download. Its Copy action can be selection-sensitive. Sequence data remain in the browser.</small></div>
        <div className="alignment-view-controls"><div className="mode-toggle"><button className={alignmentMode === "nt" ? "active" : ""} type="button" onClick={() => setAlignmentMode("nt")}>Nucleotide</button><button className={alignmentMode === "aa" ? "active" : ""} type="button" onClick={() => setAlignmentMode("aa")}>Amino acid</button></div><span>{alignmentInfo?.rows.toLocaleString() ?? parseFasta(alignment, true).length.toLocaleString()} aligned rows · {alignmentInfo?.columns.toLocaleString() ?? "—"} columns · {alignmentSource || "alignment"} · fingerprint {alignmentInfo?.fingerprint ?? "—"}</span></div>
        <AlignmentPreview fasta={alignment} mode={alignmentMode} />
        <div ref={treeResultRef} className="tree-operation-region">
          <div className="tree-controls"><div><span className="section-kicker">On-demand tree</span><h4>FastTree 2.1.11 double-precision WASM</h4><p>The exact current nucleotide alignment is rewritten into the WASM filesystem before every run. Rooting is a separate post-inference operation.</p></div><label><span>Model</span><select value={treeModel} onChange={(event) => setTreeModel(event.target.value as "gtr" | "jc")}><option value="gtr">GTR</option><option value="jc">Jukes–Cantor</option></select></label><label className="check-line"><input type="checkbox" checked={treeFast} onChange={(event) => setTreeFast(event.target.checked)} /><span>Fastest heuristic</span></label><button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void inferTree()}>{busy.startsWith("Running FastTree") ? "Inferring tree…" : "Infer tree"}</button></div>
          {treeError && <div className="inline-method-error tree-error" role="alert"><strong>Tree inference stopped</strong><span>{treeError}</span><button type="button" onClick={() => setTreeError("")}>Dismiss</button></div>}
          {treeRun && <>
            <div className="tree-provenance"><div><span className="section-kicker">Executed input</span><strong>{treeRun.rows.toLocaleString()} rows × {treeRun.columns.toLocaleString()} columns</strong><small>{treeRun.source} · alignment fingerprint {treeRun.fingerprint}</small><code>{treeRun.command}</code></div><div className="result-actions"><button type="button" onClick={() => downloadText(treeRun.alignmentFasta, `${baseName(inputName)}.lineage-${selectedLineage.id}.fasttree-alignment.fasta`)}>Named input FASTA ↓</button><button type="button" onClick={() => downloadText(treeRun.inputFasta, `${baseName(inputName)}.lineage-${selectedLineage.id}.fasttree-exact-input.fasta`)}>Exact numeric input ↓</button></div></div>
            <div className="tree-output-switch"><div className="mode-toggle"><button className={treeViewMode === "rooted" ? "active" : ""} type="button" onClick={() => setTreeViewMode("rooted")}>Rooted · resolved</button><button className={treeViewMode === "stable" ? "active" : ""} type="button" onClick={() => setTreeViewMode("stable")}>Rooted · floor-collapsed</button><button className={treeViewMode === "raw" ? "active" : ""} type="button" onClick={() => setTreeViewMode("raw")}>Raw FastTree</button></div><span>{treeRun.collapsedEdges.toLocaleString()} numerical-floor internal edges can optionally be displayed as polytomies; the complete resolved tree is the default.</span></div>
            <LineageTreeViewer
              newick={treeViewMode === "stable" ? treeRun.stableNewick : treeViewMode === "rooted" ? treeRun.rootedNewick : treeRun.newick}
              alignmentFasta={treeRun.alignmentFasta}
              rows={lineageRows}
              multiplicityByOrdinal={lineageMultiplicity}
              variant={treeViewMode}
              collapsedEdges={treeRun.collapsedEdges}
              collapseThreshold={treeRun.collapseThreshold}
              mode={alignmentMode}
              onModeChange={setAlignmentMode}
              isTcr={isTcr}
              name={`${baseName(inputName)}.lineage-${selectedLineage.id}.${treeViewMode === "stable" ? "rooted-floor-collapsed" : treeViewMode === "rooted" ? "rooted-resolved" : "raw-fasttree"}-tree`}
            />
          </>}
        </div>
      </>}
    </section>}
  </section>;
}
