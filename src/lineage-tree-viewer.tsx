import { useEffect, useMemo, useRef, useState } from "react";

import { inferKabatColumns, type KabatColumnMap } from "./kabat-numbering";
import { CommitTextInput } from "./commit-text-input";
import { GERMLINE_OUTGROUP } from "./lineage-alignment";
import {
  aminoAcidBranchMutations,
  alignmentRegionMap,
  aminoAcidRegionMap,
  cladeSignature,
  columnsForRegionPreset,
  mapParsimonyMutations,
  motifCellMap,
  ordinalFromAlignmentName,
  parseColumnSelection,
  spacedColumnOffsets,
  translateAlignedNucleotides,
  type RegionPreset,
  type VariableRegion,
} from "./lineage-phylogeny";
import { ladderizeTree, layoutTree, parseNewick, serializeNewick } from "./phylogeny";
import { parseFasta } from "./post-analysis-core";
import type { AirrDetailRow } from "./result-store";
import { sequenceColor } from "./sequence-colors";
import { categoricalLineageColor, categoricalValueColor, sampleColor, type SampleColorMap } from "./sample-colors";

export type CoordinatedTreeVariant = "stable" | "rooted" | "raw";

interface Props {
  newick: string;
  alignmentFasta: string;
  rows: AirrDetailRow[];
  multiplicityByOrdinal: Map<number, number>;
  name: string;
  variant: CoordinatedTreeVariant;
  collapsedEdges?: number;
  collapseThreshold?: number;
  mode: "nt" | "aa";
  onModeChange: (mode: "nt" | "aa") => void;
  frameOffset: 0 | 1 | 2;
  isTcr: boolean;
  sampleColors: SampleColorMap;
  lineageByOrdinal: Map<number, number>;
}

const REGION_LABELS: Record<VariableRegion, string> = {
  fwr1: "FWR1", cdr1: "CDR1", fwr2: "FWR2", cdr2: "CDR2", fwr3: "FWR3", cdr3: "CDR3", fwr4: "FWR4",
};
const REGION_COLORS: Record<VariableRegion, string> = {
  fwr1: "#dfe9e4", cdr1: "#f5d9a8", fwr2: "#dfe9e4", cdr2: "#edb9ae", fwr3: "#dfe9e4", cdr3: "#d4c4eb", fwr4: "#dfe9e4",
};
const MOTIF_COLORS = ["#d64b64", "#3f78c5", "#da8d1e", "#128a79", "#8b5cc4", "#a1633f", "#d85fa6", "#55713a"];

type TipColorMode = "sample" | "lineage" | "isotype" | "constant" | "subject" | "cohort" | "timepoint" | "compartment" | "v_gene" | "j_gene" | "productive" | "double_d" | "uniform";

const TIP_COLOR_LABELS: Record<TipColorMode, string> = {
  sample: "Sample", lineage: "Original lineage", isotype: "Isotype", constant: "Constant gene",
  subject: "Donor / subject", cohort: "Cohort", timepoint: "Timepoint", compartment: "Compartment",
  v_gene: "V gene", j_gene: "J gene", productive: "Productivity", double_d: "Double-D status", uniform: "Uniform",
};

function topGene(value: string): string {
  return value.split(",", 1)[0]?.trim().replace(/\*.*$/, "") || "Unassigned";
}

function truthLabel(value: string): string {
  const normalized = value.trim().toLowerCase();
  if (normalized === "t" || normalized === "true" || normalized === "yes" || normalized === "1") return "Productive";
  if (normalized === "f" || normalized === "false" || normalized === "no" || normalized === "0") return "Non-productive";
  return "Unassigned";
}

function tipCategory(row: AirrDetailRow | undefined, mode: TipColorMode, originalLineage: number): string {
  if (mode === "uniform") return "Uniform";
  if (mode === "lineage") return originalLineage > 0 ? `Lineage ${originalLineage}` : "Unassigned";
  if (!row) return "Unassigned";
  if (mode === "sample") return row.values.sample_id || row.record.sampleId || "Unassigned";
  if (mode === "isotype") return row.values.isotype || row.record.isotype || "Unassigned";
  if (mode === "constant") return topGene(row.values.c_call || row.record.cCall || "");
  if (mode === "subject") return row.values.subject_id || row.record.subjectId || "Unassigned";
  if (mode === "cohort") return row.values.swig_cohort || row.record.cohort || "Unassigned";
  if (mode === "timepoint") return row.values.swig_timepoint || row.record.timepoint || "Unassigned";
  if (mode === "compartment") return row.values.swig_compartment || row.record.compartment || "Unassigned";
  if (mode === "v_gene") return topGene(row.values.v_call || row.record.vCall || "");
  if (mode === "j_gene") return topGene(row.values.j_call || row.record.jCall || "");
  if (mode === "productive") return truthLabel(row.values.productive || row.record.productive || "");
  return row.values.d2_call || row.record.d2Call ? "Double-D positive" : "No double-D call";
}

function downloadText(value: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([value], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1_000);
}

function saveSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("width", svg.viewBox.baseVal.width.toString());
  clone.setAttribute("height", svg.viewBox.baseVal.height.toString());
  downloadText(new XMLSerializer().serializeToString(clone), name, "image/svg+xml;charset=utf-8");
}

function textColor(background: string): string {
  const red = Number.parseInt(background.slice(1, 3), 16);
  const green = Number.parseInt(background.slice(3, 5), 16);
  const blue = Number.parseInt(background.slice(5, 7), 16);
  return red * 0.299 + green * 0.587 + blue * 0.114 < 135 ? "#ffffff" : "#17231f";
}

function shortName(name: string): string {
  return name.length > 28 ? `${name.slice(0, 26)}…` : name;
}

function niceBranchScale(maximum: number): number {
  if (!(maximum > 0)) return 0;
  const target = maximum / 4;
  const power = 10 ** Math.floor(Math.log10(target));
  const normalized = target / power;
  const unit = normalized >= 5 ? 5 : normalized >= 2 ? 2 : 1;
  return unit * power;
}

function formatBranchScale(value: number): string {
  if (value >= 0.01) return value.toPrecision(2).replace(/\.0$/, "");
  if (value >= 0.0001) return value.toFixed(4).replace(/0+$/, "").replace(/\.$/, "");
  return value.toExponential(1);
}

/** Radius is proportional to sqrt(count), so circle area represents count. */
function abundanceBubbleRadius(count: number): { radius: number; capped: boolean } {
  const uncapped = 2.8 * Math.sqrt(Math.max(1, count));
  return { radius: Math.min(14, uncapped), capped: uncapped > 14 };
}

function regionRuns(columns: number[], regions: Array<VariableRegion | null>) {
  const runs: Array<{ start: number; end: number; region: VariableRegion | null }> = [];
  columns.forEach((column, displayIndex) => {
    const region = regions[column] ?? null;
    const previous = runs.at(-1);
    if (previous && previous.region === region && columns[displayIndex - 1] + 1 === column) previous.end = displayIndex;
    else runs.push({ start: displayIndex, end: displayIndex, region });
  });
  return runs;
}

export function LineageTreeViewer({
  newick,
  alignmentFasta,
  rows,
  multiplicityByOrdinal,
  name,
  variant,
  collapsedEdges = 0,
  collapseThreshold = 1e-8,
  mode,
  onModeChange,
  frameOffset,
  isTcr,
  sampleColors,
  lineageByOrdinal,
}: Props) {
  const svgRef = useRef<SVGSVGElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const [layoutMode, setLayoutMode] = useState<"phylogram" | "cladogram">("phylogram");
  const [ladderization, setLadderization] = useState<"none" | "large-first" | "small-first">("none");
  const [treeWidth, setTreeWidth] = useState(620);
  const [rowHeight, setRowHeight] = useState(26);
  const [cellWidth, setCellWidth] = useState(13);
  const [showMutations, setShowMutations] = useState(false);
  const [mutationLabelLimit, setMutationLabelLimit] = useState(2);
  const [mutationFontSize, setMutationFontSize] = useState(8);
  const [numbering, setNumbering] = useState<"alignment" | "kabat">("alignment");
  const [kabat, setKabat] = useState<KabatColumnMap | null>(null);
  const [kabatStatus, setKabatStatus] = useState("");
  const [preset, setPreset] = useState<RegionPreset | "motifs">("full");
  const [customColumns, setCustomColumns] = useState("");
  const [colorMode, setColorMode] = useState<"residue" | "motif">("residue");
  const [motifText, setMotifText] = useState("");
  const [fullscreen, setFullscreen] = useState(false);
  const [tipColorMode, setTipColorMode] = useState<TipColorMode>(() => new Set(lineageByOrdinal.values()).size > 1 ? "lineage" : Object.keys(sampleColors).length > 1 ? "sample" : "uniform");

  const records = useMemo(() => parseFasta(alignmentFasta, true), [alignmentFasta]);
  const nucleotideByName = useMemo(() => new Map(records.map((record) => [record.name, record.sequence])), [records]);
  const displayedByName = useMemo(() => new Map(records.map((record) => [record.name, mode === "nt" ? record.sequence : translateAlignedNucleotides(record.sequence, frameOffset)])), [records, mode, frameOffset]);
  const tree = useMemo(() => parseNewick(newick), [newick]);
  const displayTree = useMemo(() => ladderization === "none" ? tree : ladderizeTree(tree, ladderization), [tree, ladderization]);
  const layout = useMemo(() => layoutTree(displayTree, treeWidth, rowHeight, layoutMode, 24, 74), [displayTree, treeWidth, rowHeight, layoutMode]);
  const nucleotideRegions = useMemo(() => alignmentRegionMap(records, rows), [records, rows]);
  const rowByOrdinal = useMemo(() => new Map(rows.map((row) => [row.record.ordinal, row])), [rows]);
  const tipLegendEntries = useMemo(() => {
    if (tipColorMode === "uniform") return [];
    const categories = new Map<string,string>();
    for (const row of rows) {
      const originalLineage = lineageByOrdinal.get(row.record.ordinal) ?? 0;
      const category = tipCategory(row, tipColorMode, originalLineage);
      const color = tipColorMode === "sample" ? sampleColor(category, sampleColors) : tipColorMode === "lineage" ? categoricalLineageColor(originalLineage) : categoricalValueColor(category);
      categories.set(category, color);
    }
    return [...categories].sort((left,right)=>left[0].localeCompare(right[0],undefined,{numeric:true})).slice(0,32).map(([label,color])=>({label,color}));
  }, [rows, lineageByOrdinal, tipColorMode, sampleColors]);
  const regions = useMemo(() => mode === "nt" ? nucleotideRegions : aminoAcidRegionMap(nucleotideRegions, frameOffset), [mode, nucleotideRegions, frameOffset]);
  const alignmentColumns = displayedByName.values().next().value?.length ?? 0;
  const alignmentLabels = useMemo(() => Array.from({ length: alignmentColumns }, (_, index) => String(index + 1)), [alignmentColumns]);
  const coordinateLabels = mode === "aa" && numbering === "kabat" && kabat ? kabat.labels : alignmentLabels;
  const motifs = useMemo(() => motifText.split(/[\n,;]+/).map((value) => value.trim()).filter(Boolean), [motifText]);
  const motifByName = useMemo(() => new Map([...displayedByName].map(([recordName, sequence]) => [recordName, motifCellMap(sequence, motifs, mode)])), [displayedByName, motifs, mode]);
  const selectedColumns = useMemo(() => {
    if (preset === "motifs") {
      const matched = new Set<number>();
      motifByName.forEach((cells) => cells.forEach((motif, column) => { if (motif) matched.add(column); }));
      return [...matched].sort((left, right) => left - right);
    }
    if (preset === "custom") return parseColumnSelection(customColumns, coordinateLabels, alignmentColumns);
    const selected = columnsForRegionPreset(regions, preset);
    return selected.length || preset !== "full" ? selected : Array.from({ length: alignmentColumns }, (_, index) => index);
  }, [preset, customColumns, coordinateLabels, alignmentColumns, regions, motifByName]);
  const selectedColumnSet = useMemo(() => new Set(selectedColumns), [selectedColumns]);
  const parsimony = useMemo(() => {
    if (variant === "raw") return null;
    try { return mapParsimonyMutations(displayTree, nucleotideByName, GERMLINE_OUTGROUP); } catch { return null; }
  }, [displayTree, nucleotideByName, variant]);
  const aminoMutationsByClade = useMemo(() => {
    const result = new Map<string, ReturnType<typeof aminoAcidBranchMutations>>();
    if (!parsimony) return result;
    for (const edge of layout.edges) {
      const childSignature = cladeSignature(edge.child);
      const parentSequence = parsimony.sequencesByClade.get(cladeSignature(edge.parent));
      const childSequence = parsimony.sequencesByClade.get(childSignature);
      if (parentSequence && childSequence) result.set(childSignature, aminoAcidBranchMutations(parentSequence, childSequence, childSignature, frameOffset));
    }
    return result;
  }, [parsimony, layout, frameOffset]);

  useEffect(() => {
    setKabat(null);
    setKabatStatus("");
  }, [frameOffset]);

  useEffect(() => {
    if (numbering !== "kabat" || mode !== "aa" || isTcr || kabat) return;
    let cancelled = false;
    setKabatStatus("Numbering aligned IG variable domains in-browser…");
    void inferKabatColumns(records, 24, frameOffset).then((result) => {
      if (cancelled) return;
      setKabat(result);
      setKabatStatus("");
    }).catch((error) => {
      if (cancelled) return;
      setNumbering("alignment");
      setKabatStatus(error instanceof Error ? error.message : String(error));
    });
    return () => { cancelled = true; };
  }, [numbering, mode, isTcr, kabat, records, frameOffset]);

  useEffect(() => {
    const changed = () => setFullscreen(document.fullscreenElement === containerRef.current);
    document.addEventListener("fullscreenchange", changed);
    return () => document.removeEventListener("fullscreenchange", changed);
  }, []);

  const labelWidth = 208;
  const matrixX = treeWidth + labelWidth;
  const displayColumnX = spacedColumnOffsets(selectedColumns, cellWidth);
  const matrixWidth = selectedColumns.length ? displayColumnX.at(-1)! + cellWidth : cellWidth;
  const svgWidth = matrixX + matrixWidth + 34;
  const svgLegendColumns = Math.max(1, Math.min(4, Math.floor(svgWidth / 230)));
  const svgLegendRows = Math.ceil(tipLegendEntries.length / svgLegendColumns);
  const svgLegendHeight = tipLegendEntries.length ? 34 + svgLegendRows * 19 : 0;
  const svgHeight = layout.height + svgLegendHeight;
  const runs = regionRuns(selectedColumns, regions);
  const leaves = layout.nodes.filter((node) => !node.children.length);
  const maximumMultiplicity = Math.max(1, ...leaves.map((leaf) => {
    const ordinal = ordinalFromAlignmentName(leaf.name);
    return ordinal === null ? 1 : multiplicityByOrdinal.get(ordinal) ?? 1;
  }));
  const variantLabel = variant === "stable" ? "Germline-rooted · numerical-floor polytomies" : variant === "rooted" ? "Germline-rooted · complete FastTree resolution" : "Raw unrooted FastTree output";
  const inferredN = records.find((record) => record.name === GERMLINE_OUTGROUP)?.sequence.split("").filter((value) => value === "N").length ?? 0;
  const scaleDistance = layout.mode === "phylogram" ? niceBranchScale(layout.maximumDistance) : 0;
  const horizontalPadding = Math.min(24, Math.max(0, treeWidth * 0.1));
  const scalePixels = scaleDistance && layout.maximumDistance ? scaleDistance / layout.maximumDistance * Math.max(0, treeWidth - 2 * horizontalPadding) : 0;

  const mutationLabel = (column: number, from: string, to: string) => {
    if (mode === "nt") return `${column + 1} ${from}→${to}`;
    const position = numbering === "kabat" && kabat?.labels[column] ? `${kabat.chain}${kabat.labels[column]}` : String(column + 1);
    return `${position} ${from}→${to}`;
  };

  const toggleFullscreen = async () => {
    if (document.fullscreenElement) await document.exitFullscreen();
    else await containerRef.current?.requestFullscreen();
  };

  return <div ref={containerRef} className={`coordinated-tree-viewer${fullscreen ? " is-fullscreen" : ""}`}>
    <header className="coordinated-tree-header">
      <div><span className="section-kicker">{variantLabel}</span><h4>Tree, abundance and aligned leaf sequences</h4><p>{layout.leaves.toLocaleString()} tips · {layout.mode} · {selectedColumns.length.toLocaleString()} of {alignmentColumns.toLocaleString()} {mode === "nt" ? "nucleotide" : `amino-acid · frame starts at nucleotide column ${frameOffset + 1}`} columns shown</p></div>
      <div className="tree-header-actions"><button type="button" onClick={() => void toggleFullscreen()}>{fullscreen ? "Exit full screen" : "Full screen"}</button><button type="button" onClick={() => saveSvg(svgRef.current, `${name}.svg`)}>SVG with alignment ↓</button><button type="button" onClick={() => downloadText(`${serializeNewick(displayTree)};\n`, `${name}.nwk`, "text/plain;charset=utf-8")}>Displayed Newick</button></div>
    </header>
    <div className="coordinated-tree-controls">
      <div className="mode-toggle"><button className={mode === "nt" ? "active" : ""} type="button" onClick={() => onModeChange("nt")}>Nucleotide</button><button className={mode === "aa" ? "active" : ""} type="button" onClick={() => onModeChange("aa")}>Amino acid</button></div>
      <div className="mode-toggle"><button className={layoutMode === "phylogram" ? "active" : ""} type="button" onClick={() => setLayoutMode("phylogram")}>Branch lengths</button><button className={layoutMode === "cladogram" ? "active" : ""} type="button" onClick={() => setLayoutMode("cladogram")}>Topology</button></div>
      <label className="check-line"><input type="checkbox" checked={showMutations} disabled={!parsimony} onChange={(event) => setShowMutations(event.target.checked)} /><span>{mode === "nt" ? "Nucleotide branch mutations" : "AA replacements only"}</span></label>
      <label className="tip-color-control"><span>Tip circles</span><select value={tipColorMode} onChange={(event)=>setTipColorMode(event.target.value as TipColorMode)}><option value="sample">Color by sample</option><option value="lineage">Color by original lineage</option><option value="isotype">Color by isotype</option><option value="constant">Color by constant gene</option><option value="subject">Color by donor / subject</option><option value="cohort">Color by cohort</option><option value="timepoint">Color by timepoint</option><option value="compartment">Color by compartment</option><option value="v_gene">Color by V gene</option><option value="j_gene">Color by J gene</option><option value="productive">Color by productivity</option><option value="double_d">Color by double-D status</option><option value="uniform">Uniform</option></select></label>
      {mode === "aa" && <div className="mode-toggle"><button className={numbering === "alignment" ? "active" : ""} type="button" onClick={() => setNumbering("alignment")}>Alignment positions</button><button className={numbering === "kabat" ? "active" : ""} type="button" disabled={isTcr} title={isTcr ? "Kabat numbering is defined for IGH, IGK and IGL, not TCR chains." : "Number IG variable domains with Kabat positions"} onClick={() => setNumbering("kabat")}>Kabat</button></div>}
    </div>
    <div className="coordinated-tree-options">
      <label><span>Alignment region</span><select value={preset} onChange={(event) => { const value = event.target.value as RegionPreset | "motifs"; setPreset(value); if (value === "motifs") setColorMode("motif"); }}><option value="full">Full alignment</option><option value="variable">Annotated variable domain</option><option value="cdrs">All CDRs</option><option value="fwr1">FWR1</option><option value="cdr1">CDR1</option><option value="fwr2">FWR2</option><option value="cdr2">CDR2</option><option value="fwr3">FWR3</option><option value="cdr3">CDR3</option><option value="fwr4">FWR4</option><option value="motifs">Motif-matched columns</option><option value="custom">Positions / ranges</option></select></label>
      {preset === "custom" && <label className="wide-option"><span>{numbering === "kabat" && mode === "aa" ? "Kabat positions or ranges" : "1-based columns or ranges"}</span><CommitTextInput value={customColumns} onCommit={setCustomColumns} placeholder={numbering === "kabat" && mode === "aa" ? "31-35B, 50, 95-102" : "95-110, 145, 220-235"} /><small>Applied on Enter or when focus leaves the field.</small></label>}
      <label><span>Cell colors</span><select value={colorMode} onChange={(event) => setColorMode(event.target.value as "residue" | "motif")}><option value="residue">Standard residue colors</option><option value="motif">Motif colors</option></select></label>
      <label><span>Tip order</span><select value={ladderization} onChange={(event) => setLadderization(event.target.value as "none" | "large-first" | "small-first")}><option value="none">Newick order</option><option value="large-first">Ladderize · large clades first</option><option value="small-first">Ladderize · small clades first</option></select></label>
      {(colorMode === "motif" || preset === "motifs") && <label className="wide-option"><span>Motifs · comma/newline separated</span><CommitTextInput value={motifText} onCommit={setMotifText} placeholder={mode === "aa" ? "YYC, AR..Y, WGQ" : "TGT, AAR, TGG"} /><small>Applied on Enter or when focus leaves the field.</small></label>}
      <label><span>Tree width · {treeWidth}px</span><input type="range" min="0" max="1400" step="10" value={treeWidth} onChange={(event) => setTreeWidth(Number(event.target.value))} /></label>
      <label><span>Tip spacing · {rowHeight}px</span><input type="range" min="0" max="64" step="1" value={rowHeight} onChange={(event) => setRowHeight(Number(event.target.value))} /></label>
      <label><span>Residue width · {cellWidth}px</span><input type="range" min="9" max="24" step="1" value={cellWidth} onChange={(event) => setCellWidth(Number(event.target.value))} /></label>
      {showMutations&&<><label><span>Mutations before +N · {mutationLabelLimit}</span><input type="range" min="1" max="20" step="1" value={mutationLabelLimit} onChange={(event)=>setMutationLabelLimit(Number(event.target.value))}/></label><label><span>Mutation font · {mutationFontSize}px</span><input type="range" min="5" max="24" step="1" value={mutationFontSize} onChange={(event)=>setMutationFontSize(Number(event.target.value))}/></label></>}
    </div>
    {kabatStatus && <div className={numbering === "kabat" ? "viewer-numbering-status" : "viewer-numbering-status warning"}>{kabatStatus}</div>}
    {mode === "aa" && numbering === "kabat" && kabat && <div className="viewer-numbering-status">Kabat {kabat.chain} consensus · {kabat.numberedColumns.toLocaleString()} columns · {kabat.contributingSequences} members · mean confidence {kabat.confidence.toFixed(2)}</div>}
    {mode === "aa" && numbering === "kabat" && kabat?.warnings.map((warning) => <div key={warning} className="viewer-numbering-status warning">{warning}</div>)}
    <div className="coordinated-tree-legend">
      <span><i className="bubble-key" />tip bubble area is proportional to multiplicity (maximum {maximumMultiplicity.toLocaleString()}; radius capped at 14 px)</span>
      {tipLegendEntries.map((entry)=><span key={`tip-${tipColorMode}-${entry.label}`}><i style={{background:entry.color}}/>{entry.label}</span>)}
      {parsimony && <span><i className="mutation-key" />{mode === "nt" ? "germline-constrained nucleotide mutations" : "nonsynonymous amino-acid replacements"} · {parsimony.score.toLocaleString()} nucleotide parsimony steps · {inferredN.toLocaleString()} germline N sites inferred</span>}
      {variant === "stable" && <span>{collapsedEdges.toLocaleString()} internal edges ≤ {collapseThreshold.toExponential(0)} shown as polytomies</span>}
      {colorMode === "motif" && motifs.map((motif, index) => <span key={`${motif}-${index}`}><i style={{ background: MOTIF_COLORS[index % MOTIF_COLORS.length] }} />{index + 1}. {motif}</span>)}
    </div>
    <div className="coordinated-tree-scroll"><svg ref={svgRef} width={svgWidth} height={svgHeight} viewBox={`0 0 ${svgWidth} ${svgHeight}`} role="img" aria-label="Rooted lineage phylogeny coordinated with aligned leaf sequences">
      <rect width={svgWidth} height={svgHeight} fill="#fbfaf5" />
      <text x="12" y="18" fontFamily="Inter,Arial,sans-serif" fontSize="11" fontWeight="700" fill="#31413c">{variantLabel} · {mode === "nt" ? "nucleotide" : numbering === "kabat" ? `amino acid · Kabat ${kabat?.chain ?? ""}` : "amino acid"}</text>
      {scaleDistance > 0 && scalePixels >= 8 && <g aria-label={`${formatBranchScale(scaleDistance)} substitutions per site scale`}>
        <line x1={horizontalPadding} x2={horizontalPadding + scalePixels} y1="48" y2="48" stroke="#273d37" strokeWidth="1.25" />
        <line x1={horizontalPadding} x2={horizontalPadding} y1="44" y2="52" stroke="#273d37" strokeWidth="1.25" />
        <line x1={horizontalPadding + scalePixels} x2={horizontalPadding + scalePixels} y1="44" y2="52" stroke="#273d37" strokeWidth="1.25" />
        <text x={horizontalPadding + scalePixels / 2} y="41" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#4e5d58">{formatBranchScale(scaleDistance)} substitutions/site</text>
      </g>}
      {runs.map((run, index) => {
        const x = matrixX + displayColumnX[run.start];
        const width = displayColumnX[run.end] - displayColumnX[run.start] + cellWidth;
        return <g key={`${run.start}-${index}`}><rect x={x} y="24" width={width} height="17" fill={run.region ? REGION_COLORS[run.region] : "#eceeea"} /><text x={x + width / 2} y="36" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="7" fontWeight="700" fill="#364640">{run.region && width >= 25 ? REGION_LABELS[run.region] : ""}</text></g>;
      })}
      {selectedColumns.map((column, displayIndex) => {
        const label = coordinateLabels[column] || String(column + 1);
        const show = displayIndex === 0 || displayIndex === selectedColumns.length - 1 || displayIndex % Math.max(1, Math.ceil(45 / cellWidth)) === 0 || selectedColumns[displayIndex - 1] + 1 !== column;
        return show ? <text key={`tick-${column}`} transform={`translate(${matrixX + displayColumnX[displayIndex] + cellWidth / 2},67) rotate(-55)`} textAnchor="end" fontFamily="ui-monospace,monospace" fontSize="7" fill="#65736e">{label}</text> : null;
      })}
      <line x1={treeWidth} x2={treeWidth} y1="23" y2={layout.height - 10} stroke="#c9d0ca" />
      {layout.edges.map((edge, index) => {
        const childSignature = cladeSignature(edge.child);
        const mutations = (mode === "nt" ? parsimony?.mutationsByClade.get(childSignature) : aminoMutationsByClade.get(childSignature))?.filter((mutation) => selectedColumnSet.has(mutation.column)) ?? [];
        const label = mutations.slice(0, mutationLabelLimit).map((mutation) => mutationLabel(mutation.column, mutation.from, mutation.to)).join(" · ") + (mutations.length > mutationLabelLimit ? ` · +${mutations.length - mutationLabelLimit}` : "");
        const labelX = edge.parent.x + (edge.child.x - edge.parent.x) * 0.5;
        return <g key={`edge-${index}`}>
          <path d={`M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}`} stroke="#3f5650" strokeWidth="0.95" strokeLinecap="square" strokeLinejoin="miter" fill="none" />
          {showMutations && label && <g>
            <line x1={labelX} x2={labelX} y1={edge.child.y - 2.5} y2={edge.child.y + 2.5} stroke="#88465b" strokeWidth="0.8" />
            <text x={labelX} y={edge.child.y - Math.max(4.5, mutationFontSize * 0.58)} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize={mutationFontSize} fontWeight="650" fill="#71374a" stroke="#fbfaf5" strokeWidth={Math.max(2.5,mutationFontSize*0.36)} paintOrder="stroke fill">{label}</text>
          </g>}
        </g>;
      })}
      {layout.nodes.map((node, index) => {
        const leaf = !node.children.length;
        const ordinal = ordinalFromAlignmentName(node.name);
        const multiplicity = ordinal === null ? 1 : Math.max(1, multiplicityByOrdinal.get(ordinal) ?? 1);
        const bubble = abundanceBubbleRadius(multiplicity);
        const radius = leaf ? bubble.radius : 1.7;
        const recordSequence = displayedByName.get(node.name) ?? "";
        const motifsForRecord = motifByName.get(node.name);
        const row = ordinal === null ? undefined : rowByOrdinal.get(ordinal);
        const sample = row?.values.sample_id || row?.record.sampleId || "";
        const originalLineage = ordinal === null ? 0 : lineageByOrdinal.get(ordinal) ?? 0;
        const category = tipCategory(row,tipColorMode,originalLineage);
        const tipFill = node.name === GERMLINE_OUTGROUP ? "#d49a19" : tipColorMode === "sample" ? sampleColor(category,sampleColors) : tipColorMode === "lineage" ? categoricalLineageColor(originalLineage) : tipColorMode === "uniform" ? "#08796f" : categoricalValueColor(category);
        return <g key={`node-${index}`}>
          {leaf && <><line x1={node.x + radius + 2} x2={treeWidth + 6} y1={node.y} y2={node.y} stroke="#b3bdb8" strokeWidth="0.65" strokeDasharray="2 3" /><text x={treeWidth + 10} y={node.y + 3.4} fontFamily="ui-monospace,monospace" fontSize="9" fontWeight={node.name === GERMLINE_OUTGROUP ? "700" : "500"} fill="#263630">{shortName(node.name)}<title>{node.name} · multiplicity {multiplicity}</title></text></>}
          <circle cx={node.x} cy={node.y} r={radius} fill={leaf ? tipFill : "#fbfaf5"} fillOpacity={leaf ? 0.9 : 1} stroke={leaf ? "#244d45" : "#3f5650"} strokeWidth={leaf ? 0.8 : 0.75}><title>{leaf ? `${node.name} · multiplicity ${multiplicity} · ${TIP_COLOR_LABELS[tipColorMode]}: ${category}${sample?` · sample ${sample}`:""}${originalLineage?` · original lineage ${originalLineage}`:""}${bubble.capped ? " · display radius capped" : ""}` : `Internal node · ${cladeSignature(node).split("\u0000").length} descendants`}</title></circle>
          {leaf && selectedColumns.map((column, displayIndex) => {
            const value = recordSequence[column] ?? "-";
            const motif = motifsForRecord?.[column] ?? 0;
            const background = colorMode === "motif" ? motif ? MOTIF_COLORS[(motif - 1) % MOTIF_COLORS.length] : value === "-" ? "#ffffff" : "#eceeea" : sequenceColor(value, mode);
            const x = matrixX + displayColumnX[displayIndex];
            return <g key={`${node.name}-${column}`}><rect x={x} y={node.y - rowHeight / 2} width={cellWidth} height={rowHeight} fill={background} /><text x={x + cellWidth / 2} y={node.y + Math.min(4, rowHeight * 0.16)} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize={Math.min(10.5, cellWidth - 2)} fontWeight="700" fill={textColor(background)}>{value}</text><title>{node.name} · {coordinateLabels[column] || column + 1} · {value}{motif ? ` · motif ${motif}` : ""}</title></g>;
          })}
        </g>;
      })}
      {tipLegendEntries.length>0&&<g aria-label={`${TIP_COLOR_LABELS[tipColorMode]} tip-color legend`}>
        <line x1="12" x2={svgWidth-12} y1={layout.height+4} y2={layout.height+4} stroke="#c9d0ca"/>
        <text x="14" y={layout.height+20} fontFamily="Inter,Arial,sans-serif" fontSize="9" fontWeight="700" fill="#31413c">Tip color · {TIP_COLOR_LABELS[tipColorMode]}</text>
        {tipLegendEntries.map((entry,index)=>{const column=index%svgLegendColumns,row=Math.floor(index/svgLegendColumns),cell=(svgWidth-28)/svgLegendColumns,x=18+column*cell,y=layout.height+35+row*19;return <g key={`svg-legend-${entry.label}`}><circle cx={x} cy={y-3} r="4" fill={entry.color} stroke="#244d45" strokeWidth=".6"/><text x={x+9} y={y} fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#34453f">{shortName(entry.label)}<title>{entry.label}</title></text></g>;})}
      </g>}
    </svg></div>
    <p className="scientific-note"><span>i</span>{variant === "raw" ? "Mutation mapping is disabled for the unrooted raw tree. Choose a germline-rooted view to reconstruct the UCA and branch changes." : mode === "nt" ? "Equal-cost nucleotide parsimony is reconstructed on demand. Known germline bases constrain the UCA; N bases are unknown A/C/G/T states inferred from descendants, and gaps are explicit indel states." : "Amino-acid branch labels are translated from the reconstructed parent and child codons. Synonymous nucleotide changes and unresolved X codons are omitted; multiple nucleotide changes in one codon are shown as one replacement."} {isTcr ? "Kabat is not defined for TCR chains, so TCR amino-acid views use alignment positions." : "Kabat numbering is computed in-browser from a multi-member consensus of IG variable-domain assignments."}</p>
  </div>;
}
