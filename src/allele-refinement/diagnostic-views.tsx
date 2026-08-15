import { useEffect, useMemo, useRef, useState, type CSSProperties } from "react";

import { CommitNumberInput } from "../commit-number-input.tsx";
import { FacetPicker } from "../facet-picker.tsx";
import type { CompiledReferences } from "../reference-pack.ts";
import { sequenceColor } from "../sequence-colors.tsx";
import {
  alignReferenceKernelInspection,
  assignmentShiftData,
  inspectReferenceEvidenceKernel,
} from "./diagnostics.ts";
import { adaptiveNeighbourOdds } from "./evidence.ts";
import { buildReferenceAlleleGraph } from "./reference-graph.ts";
import type {
  AlleleRefinementOptions,
  RefinementModelSummary,
  RefinementSegment,
} from "./types.ts";

const SEGMENTS: RefinementSegment[] = ["V", "D", "J"];

function formatProbability(value: number): string {
  const percent = value * 100;
  if (percent === 0) return "0%";
  if (percent < 0.0001) return "<0.0001%";
  if (percent < 0.01) return `${percent.toFixed(4)}%`;
  if (percent < 1) return `${percent.toFixed(3)}%`;
  return `${percent.toFixed(2)}%`;
}

function download(value: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([value], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

function safeName(value: string): string {
  return value.replace(/[^a-z0-9_.-]+/gi, "-").replace(/^-+|-+$/g, "").slice(0, 120) || "allele-pool";
}

function csvCell(value: string | number): string {
  const text = String(value);
  return /[\n\r,"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function kernelCellStyle(base: string, primaryBase: string, primary: boolean, mode: "nucleotide" | "highlighter"): CSSProperties {
  if (mode === "nucleotide") return { backgroundColor: sequenceColor(base, "nt"), color: "#10201c" };
  if (primary) return { backgroundColor: base === "-" ? "#ffffff" : "#e8ece9", color: "#26332f" };
  if (base === primaryBase) return { backgroundColor: base === "-" ? "#ffffff" : "#f0f1ee", color: "#8a938f" };
  if (base === "-" || primaryBase === "-") return { backgroundColor: "#efc267", color: "#302814" };
  return { backgroundColor: "#df735f", color: "#ffffff" };
}

export function ReferenceKernelInspector({
  references,
  options,
}: {
  references: CompiledReferences;
  options: AlleleRefinementOptions;
}) {
  const availableSegments = useMemo(
    () => SEGMENTS.filter((segment) => references[segment].trim().length > 0),
    [references],
  );
  const [segment, setSegment] = useState<RefinementSegment>(availableSegments[0] ?? "V");
  const [selectedName, setSelectedName] = useState("");
  const [assumedShm, setAssumedShm] = useState(0.05);
  const [mode, setMode] = useState<"nucleotide" | "highlighter">("highlighter");

  useEffect(() => {
    if (!availableSegments.includes(segment)) setSegment(availableSegments[0] ?? "V");
  }, [availableSegments, segment]);

  const graph = useMemo(
    () => buildReferenceAlleleGraph(references[segment], segment, options.neighbourRadius),
    [options.neighbourRadius, references, segment],
  );
  const alleleNames = useMemo(
    () => [...graph.callToNode.keys()].sort((left, right) => left.localeCompare(right, undefined, { numeric: true })),
    [graph],
  );
  const effectiveName = graph.callToNode.has(selectedName) ? selectedName : alleleNames[0] ?? "";
  const selectedNode = graph.callToNode.get(effectiveName) ?? -1;
  const inspection = useMemo(
    () => inspectReferenceEvidenceKernel(graph, selectedNode, assumedShm, options),
    [assumedShm, graph, options, selectedNode],
  );
  const alignment = useMemo(
    () => inspection ? alignReferenceKernelInspection(inspection) : null,
    [inspection],
  );
  const maximumAlternative = Math.max(0, ...(inspection?.alternatives.map((alternative) => alternative.probability) ?? []));
  const oneSnpOdds = adaptiveNeighbourOdds(assumedShm, options);
  const oneSnpProbability = oneSnpOdds / (1 + oneSnpOdds);

  if (!availableSegments.length) return <div className="scientific-note warning"><span>!</span>No V, D, or J reference sequences are available for this analysis.</div>;

  return <section className="allele-kernel-inspector" aria-label="Reference-neighbour error model inspector">
    <header>
      <div><span className="section-kicker">Local evidence diagnostic</span><h4>Inspect the reference-neighbour error model</h4><p>Select a database sequence and an assumed read SHM level. The table recomputes the same sparse, pre-repertoire evidence kernel used during fitting.</p></div>
      <div className="mode-toggle" aria-label="Reference alignment coloring">
        <button type="button" className={mode === "nucleotide" ? "active" : ""} onClick={() => setMode("nucleotide")}>Nucleotide</button>
        <button type="button" className={mode === "highlighter" ? "active" : ""} onClick={() => setMode("highlighter")}>Differences</button>
      </div>
    </header>
    <div className="allele-kernel-controls">
      <label title="Choose which reference segment database to inspect. This does not enable or disable that segment in the repertoire model."><span>Reference segment</span><select value={segment} onChange={(event) => { setSegment(event.target.value as RefinementSegment); setSelectedName(""); }}>{availableSegments.map((value) => <option key={value} value={value}>{value} database</option>)}</select></label>
      <FacetPicker label="Reference allele" value={effectiveName} items={alleleNames.map((value) => ({ value }))} placeholder="Choose an allele" onChange={(value) => setSelectedName(value || alleleNames[0] || "")} help="Search and select one allele label from the active reference database. Sequence-identical labels resolve to the same reference node." />
      <label title="The assumed best-reference mutation fraction used by the SHM-adaptive neighbour leakage term. The model clamps it at the configured SHM estimate cap."><span>Assumed read SHM</span><div className="range-number"><input type="range" min="0" max="0.4" step="0.005" value={assumedShm} onChange={(event) => setAssumedShm(Number(event.target.value))} /><b>{(assumedShm * 100).toFixed(1)}%</b></div></label>
    </div>
    {inspection && alignment ? <>
      <div className="allele-kernel-summary">
        <div><span>Selected-reference mass</span><strong>{formatProbability(inspection.primaryProbability)}</strong></div>
        <div><span>Alternative mass</span><strong>{formatProbability(inspection.alternativeProbability)}</strong></div>
        <div><span>One-SNP pairwise probability</span><strong>{formatProbability(oneSnpProbability)}</strong></div>
        <div><span>Non-primary candidates</span><strong>{inspection.alternatives.length.toLocaleString()}</strong></div>
      </div>
      <div className="allele-kernel-legend"><span><i className="kernel-same" />same as selected</span><span><i className="kernel-substitution" />substitution</span><span><i className="kernel-indel" />indel column</span><small>Alternative bar lengths are scaled to the largest non-primary probability; printed percentages use the complete normalized row, including the primary.</small></div>
      <div className="allele-kernel-table-wrap">
        <table style={{ minWidth: `${520 + alignment.columns * 11}px` }}>
          <thead><tr><th>Reference</th><th>Probability</th><th>Relative alternative mass</th><th>Aligned nucleotide sequence</th></tr></thead>
          <tbody>{alignment.rows.map((row) => {
            const label = row.names.join(", ");
            return <tr key={row.nodeIndex} className={row.primary ? "primary" : ""}>
              <th scope="row" title={label}><strong>{label}</strong><small>{row.primary ? "selected reference" : `edit distance ${row.distance} · ${row.substitutionOnly ? "substitution-only" : "includes indel"}`}</small></th>
              <td>{formatProbability(row.probability)}</td>
              <td>{row.primary ? <small>bar omitted</small> : <div className="allele-kernel-bar" role="img" aria-label={`${formatProbability(row.probability)} local probability`}><i style={{ width: `${maximumAlternative > 0 ? row.probability / maximumAlternative * 100 : 0}%` }} /></div>}</td>
              <td><code className="allele-kernel-sequence" aria-label={row.sequence}>{[...row.sequence].map((base, column) => <span key={column} title={`Alignment column ${column + 1}: ${base}`} style={kernelCellStyle(base, alignment.rows[0].sequence[column], row.primary, mode)}>{base}</span>)}</code></td>
            </tr>;
          })}</tbody>
        </table>
      </div>
      <p className="scientific-note"><span>i</span>{alignment.strippedAllGapColumns.toLocaleString()} all-gap alignment columns removed. Sequence-identical labels are shown as one unresolved reference node. The selected-reference row remains visible for comparison but deliberately has no bar.</p>
      {inspection.truncated && <p className="scientific-note warning"><span>!</span>The candidate cap truncated this diagnostic row. Increase “Candidate cap / record” to inspect the complete configured neighbourhood.</p>}
    </> : <p className="scientific-note warning"><span>!</span>The selected database sequence could not be resolved in the active reference graph.</p>}
  </section>;
}

function modelLabel(model: RefinementModelSummary): string {
  return `${model.scopeValue || "all data"} · ${model.locus || "unknown locus"} · ${model.segment} · ${model.rows.toLocaleString()} rows`;
}

function niceFrequencyMaximum(value: number): number {
  if (value <= 0.01) return 0.01;
  if (value <= 0.025) return 0.025;
  if (value <= 0.05) return 0.05;
  if (value <= 0.1) return 0.1;
  if (value <= 0.25) return 0.25;
  if (value <= 0.5) return 0.5;
  return 1;
}

export function AlleleAssignmentShiftChart({ models }: { models: RefinementModelSummary[] }) {
  const [modelKey, setModelKey] = useState(models[0]?.key ?? "");
  const [allelesShown, setAllelesShown] = useState(20);
  const svgRef = useRef<SVGSVGElement>(null);
  useEffect(() => {
    if (!models.some((model) => model.key === modelKey)) setModelKey(models[0]?.key ?? "");
  }, [modelKey, models]);
  const model = models.find((candidate) => candidate.key === modelKey) ?? models[0];
  const allRows = useMemo(() => model ? assignmentShiftData(model) : [], [model]);
  const visible = allRows.slice(0, Math.max(1, Math.floor(allelesShown)));
  const maximum = niceFrequencyMaximum(Math.max(0, ...visible.flatMap((row) => [row.before, row.after])));
  const width = 920;
  const labelWidth = 250;
  const plotWidth = 500;
  const valueX = 774;
  const top = 66;
  const rowHeight = 40;
  const height = Math.max(180, top + visible.length * rowHeight + 28);
  const exportStem = safeName(`${model?.scopeValue ?? "pool"}-${model?.locus ?? "locus"}-${model?.segment ?? "segment"}-assignment-shift`);

  if (!model) return null;
  const saveSvg = () => {
    if (!svgRef.current) return;
    const clone = svgRef.current.cloneNode(true) as SVGSVGElement;
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("width", String(width));
    clone.setAttribute("height", String(height));
    download(`<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}`, `${exportStem}.svg`, "image/svg+xml;charset=utf-8");
  };
  const saveCsv = () => {
    const header = ["pool", "locus", "segment", "allele", "before_frequency", "after_frequency", "delta", "before_expected_assignments", "after_expected_assignments"];
    const rows = allRows.map((row) => [model.scopeValue, model.locus, model.segment, row.label, row.before, row.after, row.delta, row.beforeAssignments, row.afterAssignments]);
    download([header, ...rows].map((row) => row.map(csvCell).join(",")).join("\n") + "\n", `${exportStem}.csv`, "text/csv;charset=utf-8");
  };

  return <article className="allele-shift-chart">
    <header>
      <div><span className="section-kicker">Assignment redistribution</span><h4>Allele frequencies before and after repertoire pooling</h4><p>Before uses normalized local evidence. After uses normalized fitted assignment responsibilities. Dirichlet prior-only mass is excluded from both series.</p></div>
      <div className="result-actions"><button type="button" onClick={saveCsv}>Data CSV ↓</button><button type="button" onClick={saveSvg}>SVG ↓</button></div>
    </header>
    <div className="allele-shift-controls">
      <label title="Each donor/study-boundary, locus, and reference segment is fitted independently. Choose the fitted pool to display."><span>Fitted pool</span><select value={model.key} onChange={(event) => setModelKey(event.target.value)}>{models.map((candidate) => <option key={candidate.key} value={candidate.key}>{modelLabel(candidate)}</option>)}</select></label>
      <label title="The chart is sorted by post-pooling frequency. Increase this display limit to include more of the fitted allele tail; CSV export always contains every modeled allele."><span>Alleles shown</span><CommitNumberInput min="1" max="500" step="1" value={allelesShown} onCommit={(value) => setAllelesShown(Math.max(1, Math.floor(value)))} /></label>
      <div className="allele-shift-legend" aria-label="Chart series"><span><i className="before" />Before</span><span><i className="after" />After</span></div>
    </div>
    <div className="allele-shift-scroll"><svg ref={svgRef} role="img" aria-label={`Allele assignment frequency before and after repertoire pooling for ${modelLabel(model)}`} width={width} height={height} viewBox={`0 0 ${width} ${height}`}>
      <rect width={width} height={height} fill="#fffdf8" />
      {[0, 0.5, 1].map((fraction) => {
        const x = labelWidth + fraction * plotWidth;
        return <g key={fraction}><line x1={x} x2={x} y1="42" y2={height - 12} stroke={fraction === 0 ? "#84918c" : "#dfe4e0"} strokeWidth={fraction === 0 ? 1.2 : 1} /><text x={x} y="32" textAnchor={fraction === 0 ? "start" : fraction === 1 ? "end" : "middle"} fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#68746f">{formatProbability(maximum * fraction)}</text></g>;
      })}
      <g transform="translate(620 15)"><rect width="12" height="8" rx="2" fill="#aeb8b2" /><text x="18" y="8" fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#46534e">Before</text><rect x="76" width="12" height="8" rx="2" fill="#2f7767" /><text x="94" y="8" fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#46534e">After</text></g>
      {visible.map((row, index) => {
        const y = top + index * rowHeight;
        const beforeWidth = row.before / maximum * plotWidth;
        const afterWidth = row.after / maximum * plotWidth;
        const shownLabel = row.label.length > 34 ? `${row.label.slice(0, 32)}…` : row.label;
        return <g key={row.nodeIndex}>
          <title>{`${row.label}\nBefore ${formatProbability(row.before)} (${row.beforeAssignments.toFixed(3)} expected assignments)\nAfter ${formatProbability(row.after)} (${row.afterAssignments.toFixed(3)} expected assignments)\nChange ${row.delta >= 0 ? "+" : ""}${formatProbability(row.delta)}`}</title>
          <text x="10" y={y + 14} fontFamily="Inter, Arial, sans-serif" fontSize="11" fontWeight="600" fill="#26332f">{shownLabel}</text>
          <rect x={labelWidth} y={y + 2} width={Math.max(0, beforeWidth)} height="9" rx="2" fill="#aeb8b2" />
          <rect x={labelWidth} y={y + 16} width={Math.max(0, afterWidth)} height="9" rx="2" fill="#2f7767" />
          <text x={valueX} y={y + 10} fontFamily="Inter, Arial, sans-serif" fontSize="10" fill="#596661">{formatProbability(row.before)}</text>
          <text x={valueX} y={y + 24} fontFamily="Inter, Arial, sans-serif" fontSize="10" fontWeight="700" fill="#245f53">{formatProbability(row.after)}</text>
          <line x1="8" x2="864" y1={y + 32} y2={y + 32} stroke="#edf0ed" />
        </g>;
      })}
      {!visible.length && <text x={width / 2} y="110" textAnchor="middle" fontFamily="Inter, Arial, sans-serif" fontSize="13" fill="#68746f">No modeled alleles in this pool</text>}
    </svg></div>
    <p className="scientific-note"><span>i</span>Rows are sorted by post-pooling assignment frequency. {allRows.length > visible.length ? `${(allRows.length - visible.length).toLocaleString()} lower-frequency modeled alleles are omitted from the figure but retained in the CSV.` : "Every modeled allele in this pool is shown."}</p>
  </article>;
}
