import { useEffect, useMemo, useRef, useState, type CSSProperties } from "react";

import { CommitNumberInput } from "../commit-number-input.tsx";
import { FacetPicker } from "../facet-picker.tsx";
import type { CompiledReferences } from "../reference-pack.ts";
import { sequenceColor } from "../sequence-colors.tsx";
import {
  alignReferenceKernelInspection,
  hardAssignmentShiftData,
  inspectReferenceEvidenceKernel,
  survivingAlleleReference,
} from "./diagnostics.ts";
import { adaptiveNeighbourOdds } from "./evidence.ts";
import { buildReferenceAlleleGraph } from "./reference-graph.ts";
import type {
  AlleleReassignmentPolicy,
  AlleleRefinementOptions,
  RefinementModelSummary,
  RefinementSegment,
  SegmentRefinementResult,
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

function csvCell(value: string | number | boolean): string {
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

function niceCountMaximum(value: number): number {
  if (!(value > 0)) return 1;
  const magnitude = 10 ** Math.floor(Math.log10(value));
  const scaled = value / magnitude;
  return (scaled <= 1 ? 1 : scaled <= 2 ? 2 : scaled <= 5 ? 5 : 10) * magnitude;
}

function formatAssignmentCount(value: number): string {
  if (Math.abs(value - Math.round(value)) < 1e-6) return Math.round(value).toLocaleString();
  return value.toLocaleString(undefined, { maximumFractionDigits: 2 });
}

type HardAssignmentFilter = "all" | "changed" | "vanished";

export function AlleleAssignmentShiftChart({
  results,
  reassignmentPolicy,
  minimumPosterior,
  weighting,
}: {
  results: SegmentRefinementResult[];
  reassignmentPolicy: AlleleReassignmentPolicy;
  minimumPosterior: number;
  weighting: AlleleRefinementOptions["weighting"];
}) {
  const models = useMemo(() => results.flatMap((result) => result.models.map((model, modelIndex) => ({ result, model, modelIndex }))), [results]);
  const [modelKey, setModelKey] = useState(models[0]?.model.key ?? "");
  const [allelesShown, setAllelesShown] = useState(20);
  const [filter, setFilter] = useState<HardAssignmentFilter>("all");
  const [minimumReferenceReads, setMinimumReferenceReads] = useState(0);
  const svgRef = useRef<SVGSVGElement>(null);
  useEffect(() => {
    if (!models.some((entry) => entry.model.key === modelKey)) setModelKey(models[0]?.model.key ?? "");
  }, [modelKey, models]);
  const selected = models.find((entry) => entry.model.key === modelKey) ?? models[0];
  const shift = useMemo(() => selected ? hardAssignmentShiftData(selected.result, selected.modelIndex, reassignmentPolicy, minimumPosterior) : null, [minimumPosterior, reassignmentPolicy, selected]);
  const allRows = shift?.rows ?? [];
  const filtered = allRows.filter((row) => filter === "all" || (filter === "changed" ? row.before !== row.after : row.vanishes));
  const visible = filtered.slice(0, Math.max(1, Math.floor(allelesShown)));
  const maximum = niceCountMaximum(Math.max(0, ...visible.flatMap((row) => [row.before, row.after])));
  const width = 920;
  const labelWidth = 250;
  const plotWidth = 500;
  const valueX = 774;
  const top = 66;
  const rowHeight = 40;
  const height = Math.max(180, top + visible.length * rowHeight + 28);
  const model = selected?.model;
  const effectiveMinimumPosterior = Math.max(0, Math.min(1, minimumPosterior));
  const exportStem = safeName(`${model?.scopeValue ?? "pool"}-${model?.locus ?? "locus"}-${model?.segment ?? "segment"}-hard-assignment-shift`);

  if (!model) return null;
  const modelAlleleByNode = new Map(model.alleles.map((allele) => [allele.nodeIndex, allele] as const));
  const saveSvg = () => {
    if (!svgRef.current) return;
    const clone = svgRef.current.cloneNode(true) as SVGSVGElement;
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("width", String(width));
    clone.setAttribute("height", String(height));
    download(`<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}`, `${exportStem}.svg`, "image/svg+xml;charset=utf-8");
  };
  const saveCsv = () => {
    const header = ["pool", "locus", "segment", "inference_model", "allele", "model_active", "inclusion_probability", "local_best_count", "after_policy_count", "delta", "vanishes", "appears", "record_weighting", "reassignment_policy", "minimum_posterior"];
    const summaryByNode = new Map(model.alleles.map((allele) => [allele.nodeIndex, allele] as const));
    const rows = allRows.map((row) => {
      const allele = summaryByNode.get(row.nodeIndex);
      return [model.scopeValue, model.locus, model.segment, model.inferenceModel ?? "dirichlet", row.label, allele?.active ?? "", allele?.inclusionProbability ?? "", row.before, row.after, row.delta, row.vanishes, row.appears, weighting, reassignmentPolicy, reassignmentPolicy === "confidence" ? effectiveMinimumPosterior : ""];
    });
    download([header, ...rows].map((row) => row.map(csvCell).join(",")).join("\n") + "\n", `${exportStem}.csv`, "text/csv;charset=utf-8");
  };
  const saveSurvivingReference = () => {
    if (!selected) return;
    const reference = survivingAlleleReference(selected.result, selected.modelIndex, reassignmentPolicy, effectiveMinimumPosterior, minimumReferenceReads);
    if (!reference) return;
    download(reference.fasta, `${safeName(`${model.scopeValue}-${model.locus}-${model.segment}`)}-surviving-alleles.fasta`, "text/plain;charset=utf-8");
  };
  const policyLabel = reassignmentPolicy === "best" ? "posterior MAP for every modeled record" : `posterior MAP at ≥ ${(effectiveMinimumPosterior * 100).toFixed(0)}% confidence, otherwise local best`;
  const fittedModelLabel = (model.inferenceModel ?? "dirichlet") === "active-set" ? "fast hurdle active-set model" : "continuous Dirichlet model";

  return <article className="allele-shift-chart">
    <header>
      <div><span className="section-kicker">Hard assignment projection</span><h4>Best-match allele counts before and after reassignment</h4><p>Each modeled record contributes its complete configured weight to one allele. Local best is the pre-pooling evidence argmax; after policy applies {policyLabel}. Counts are recomputed from the fitted {fittedModelLabel} and the currently selected decision policy.</p></div>
      <div className="result-actions"><button type="button" disabled={!shift} onClick={saveCsv}>Data CSV ↓</button><button type="button" disabled={!shift} onClick={saveSvg}>SVG ↓</button></div>
    </header>
    <div className="allele-shift-controls">
      <label title="Each donor/study-boundary, locus, and reference segment is fitted independently. Choose the fitted pool to display."><span>Fitted pool</span><select value={model.key} onChange={(event) => setModelKey(event.target.value)}>{models.map((entry) => <option key={entry.model.key} value={entry.model.key}>{modelLabel(entry.model)}</option>)}</select></label>
      <label title="Restrict the display without changing the model or CSV export."><span>Show</span><select value={filter} onChange={(event) => setFilter(event.target.value as HardAssignmentFilter)}><option value="all">All hard-assigned alleles</option><option value="changed">Changed counts only</option><option value="vanished">Vanished alleles only</option></select></label>
      <label title="Increase this display limit to include more rows; CSV export always contains every hard-assigned allele in the selected pool."><span>Alleles shown</span><CommitNumberInput min="1" max="500" step="1" value={allelesShown} onCommit={(value) => setAllelesShown(Math.max(1, Math.floor(value)))} /></label>
      <div className="allele-shift-legend" aria-label="Chart series"><span><i className="before" />Local best</span><span><i className="after" />After policy</span></div>
    </div>
    {!shift ? <p className="scientific-note warning"><span>!</span>This saved result predates hard-assignment tracking. Refit the repertoire allele model to generate the count projection.</p> : <div className="allele-shift-scroll"><svg ref={svgRef} role="img" aria-label={`Hard best-match allele counts before and after reassignment for ${modelLabel(model)}`} width={width} height={height} viewBox={`0 0 ${width} ${height}`}>
      <rect width={width} height={height} fill="#fffdf8" />
      {[0, 0.5, 1].map((fraction) => {
        const x = labelWidth + fraction * plotWidth;
        return <g key={fraction}><line x1={x} x2={x} y1="42" y2={height - 12} stroke={fraction === 0 ? "#84918c" : "#dfe4e0"} strokeWidth={fraction === 0 ? 1.2 : 1} /><text x={x} y="32" textAnchor={fraction === 0 ? "start" : fraction === 1 ? "end" : "middle"} fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#68746f">{formatAssignmentCount(maximum * fraction)}</text></g>;
      })}
      <g transform="translate(590 15)"><rect width="12" height="8" rx="2" fill="#aeb8b2" /><text x="18" y="8" fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#46534e">Local best</text><rect x="104" width="12" height="8" rx="2" fill="#2f7767" /><text x="122" y="8" fontFamily="Inter, Arial, sans-serif" fontSize="11" fill="#46534e">After policy</text></g>
      {visible.map((row, index) => {
        const y = top + index * rowHeight;
        const beforeWidth = row.before / maximum * plotWidth;
        const afterWidth = row.after / maximum * plotWidth;
        const shownLabel = row.label.length > 34 ? `${row.label.slice(0, 32)}…` : row.label;
        const allele = modelAlleleByNode.get(row.nodeIndex);
        const inclusionDetail = allele?.inclusionProbability === undefined ? "" : `\nModel active ${allele.active ? "yes" : "no"}\nInclusion probability ${(allele.inclusionProbability * 100).toFixed(2)}%`;
        return <g key={row.nodeIndex}>
          <title>{`${row.label}${inclusionDetail}\nLocal best ${formatAssignmentCount(row.before)}\nAfter policy ${formatAssignmentCount(row.after)}\nChange ${row.delta >= 0 ? "+" : ""}${formatAssignmentCount(row.delta)}${row.vanishes ? "\nVanishes from this hard-assigned pool" : row.appears ? "\nAppears after reassignment" : ""}`}</title>
          <text x="10" y={y + 14} fontFamily="Inter, Arial, sans-serif" fontSize="11" fontWeight="600" fill={row.vanishes ? "#a33b32" : "#26332f"}>{shownLabel}</text>
          <rect x={labelWidth} y={y + 2} width={Math.max(0, beforeWidth)} height="9" rx="2" fill="#aeb8b2" />
          <rect x={labelWidth} y={y + 16} width={Math.max(0, afterWidth)} height="9" rx="2" fill="#2f7767" />
          <text x={valueX} y={y + 10} fontFamily="Inter, Arial, sans-serif" fontSize="10" fill="#596661">{formatAssignmentCount(row.before)}</text>
          <text x={valueX} y={y + 24} fontFamily="Inter, Arial, sans-serif" fontSize="10" fontWeight="700" fill={row.vanishes ? "#a33b32" : "#245f53"}>{row.vanishes ? "0 · vanishes" : formatAssignmentCount(row.after)}</text>
          <line x1="8" x2="864" y1={y + 32} y2={y + 32} stroke="#edf0ed" />
        </g>;
      })}
      {!visible.length && <text x={width / 2} y="110" textAnchor="middle" fontFamily="Inter, Arial, sans-serif" fontSize="13" fill="#68746f">No alleles match this display filter</text>}
    </svg></div>}
    {shift && <div className="allele-surviving-reference"><div><strong>Surviving allele reference</strong><small>Uses this fitted pool, the current reassignment policy, confidence gate, and record weighting. A threshold of 0 retains every candidate reference; use 1 to remove zero-count alleles.</small></div><label title="Exclude reference nodes with fewer hard post-reassignment reads than this value."><span>Minimum post-reassignment reads</span><CommitNumberInput min="0" step="1" value={minimumReferenceReads} onCommit={(value) => setMinimumReferenceReads(Math.max(0, value))} /></label><button type="button" onClick={saveSurvivingReference}>Download surviving allele reference ↓</button></div>}
    {shift && <p className="scientific-note"><span>i</span>{formatAssignmentCount(shift.totalAssignments)} {weighting === "abundance" ? "duplicate-count-weighted" : "unique-record"} hard assignments in this pool; {formatAssignmentCount(shift.changedAssignments)} change under this policy. {reassignmentPolicy === "confidence" ? `${formatAssignmentCount(shift.heldBelowConfidence)} differing posterior MAP assignments are held below confidence and remain at local best in this projection; the applied AIRR overlay preserves their original call strings.` : "Every modeled record is projected to its posterior MAP."} {shift.vanishedAlleles.toLocaleString()} allele class{shift.vanishedAlleles === 1 ? "" : "es"} vanish and {shift.appearedAlleles.toLocaleString()} appear. {filtered.length > visible.length ? `${(filtered.length - visible.length).toLocaleString()} matching rows are omitted from the figure but retained in the CSV.` : "Every matching row is shown."}</p>}
  </article>;
}
