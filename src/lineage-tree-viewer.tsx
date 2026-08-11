import { useEffect, useMemo, useRef, useState } from "react";

import { inferKabatColumns, type KabatColumnMap } from "./kabat-numbering";
import { GERMLINE_OUTGROUP } from "./lineage-alignment";
import {
  alignmentRegionMap,
  aminoAcidRegionMap,
  cladeSignature,
  columnsForRegionPreset,
  mapParsimonyMutations,
  motifCellMap,
  ordinalFromAlignmentName,
  parseColumnSelection,
  translateAlignedNucleotides,
  type RegionPreset,
  type VariableRegion,
} from "./lineage-phylogeny";
import { layoutTree, parseNewick } from "./phylogeny";
import { parseFasta } from "./post-analysis-core";
import type { AirrDetailRow } from "./result-store";
import { sequenceColor } from "./sequence-colors";

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
  isTcr: boolean;
}

const REGION_LABELS: Record<VariableRegion, string> = {
  fwr1: "FWR1", cdr1: "CDR1", fwr2: "FWR2", cdr2: "CDR2", fwr3: "FWR3", cdr3: "CDR3", fwr4: "FWR4",
};
const REGION_COLORS: Record<VariableRegion, string> = {
  fwr1: "#dfe9e4", cdr1: "#f5d9a8", fwr2: "#dfe9e4", cdr2: "#edb9ae", fwr3: "#dfe9e4", cdr3: "#d4c4eb", fwr4: "#dfe9e4",
};
const MOTIF_COLORS = ["#d64b64", "#3f78c5", "#da8d1e", "#128a79", "#8b5cc4", "#a1633f", "#d85fa6", "#55713a"];

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
  isTcr,
}: Props) {
  const svgRef = useRef<SVGSVGElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const [layoutMode, setLayoutMode] = useState<"phylogram" | "cladogram">("phylogram");
  const [treeWidth, setTreeWidth] = useState(620);
  const [rowHeight, setRowHeight] = useState(26);
  const [cellWidth, setCellWidth] = useState(13);
  const [showMutations, setShowMutations] = useState(false);
  const [numbering, setNumbering] = useState<"alignment" | "kabat">("alignment");
  const [kabat, setKabat] = useState<KabatColumnMap | null>(null);
  const [kabatStatus, setKabatStatus] = useState("");
  const [preset, setPreset] = useState<RegionPreset | "motifs">("full");
  const [customColumns, setCustomColumns] = useState("");
  const [colorMode, setColorMode] = useState<"residue" | "motif">("residue");
  const [motifText, setMotifText] = useState("");
  const [fullscreen, setFullscreen] = useState(false);

  const records = useMemo(() => parseFasta(alignmentFasta, true), [alignmentFasta]);
  const nucleotideByName = useMemo(() => new Map(records.map((record) => [record.name, record.sequence])), [records]);
  const displayedByName = useMemo(() => new Map(records.map((record) => [record.name, mode === "nt" ? record.sequence : translateAlignedNucleotides(record.sequence)])), [records, mode]);
  const tree = useMemo(() => parseNewick(newick), [newick]);
  const layout = useMemo(() => layoutTree(tree, treeWidth, rowHeight, layoutMode, 24, 74), [tree, treeWidth, rowHeight, layoutMode]);
  const nucleotideRegions = useMemo(() => alignmentRegionMap(records, rows), [records, rows]);
  const regions = useMemo(() => mode === "nt" ? nucleotideRegions : aminoAcidRegionMap(nucleotideRegions), [mode, nucleotideRegions]);
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
    try { return mapParsimonyMutations(tree, nucleotideByName, GERMLINE_OUTGROUP); } catch { return null; }
  }, [tree, nucleotideByName, variant]);

  useEffect(() => {
    if (numbering !== "kabat" || mode !== "aa" || isTcr || kabat) return;
    let cancelled = false;
    setKabatStatus("Numbering aligned IG variable domains in-browser…");
    void inferKabatColumns(records).then((result) => {
      if (cancelled) return;
      setKabat(result);
      setKabatStatus("");
    }).catch((error) => {
      if (cancelled) return;
      setNumbering("alignment");
      setKabatStatus(error instanceof Error ? error.message : String(error));
    });
    return () => { cancelled = true; };
  }, [numbering, mode, isTcr, kabat, records]);

  useEffect(() => {
    const changed = () => setFullscreen(document.fullscreenElement === containerRef.current);
    document.addEventListener("fullscreenchange", changed);
    return () => document.removeEventListener("fullscreenchange", changed);
  }, []);

  const labelWidth = 208;
  const matrixX = treeWidth + labelWidth;
  const svgWidth = matrixX + Math.max(1, selectedColumns.length) * cellWidth + 34;
  const svgHeight = layout.height;
  const runs = regionRuns(selectedColumns, regions);
  const leaves = layout.nodes.filter((node) => !node.children.length);
  const maximumMultiplicity = Math.max(1, ...leaves.map((leaf) => {
    const ordinal = ordinalFromAlignmentName(leaf.name);
    return ordinal === null ? 1 : multiplicityByOrdinal.get(ordinal) ?? 1;
  }));
  const variantLabel = variant === "stable" ? "Germline-rooted · numerical-floor polytomies" : variant === "rooted" ? "Germline-rooted · complete FastTree resolution" : "Raw unrooted FastTree output";
  const inferredN = records.find((record) => record.name === GERMLINE_OUTGROUP)?.sequence.split("").filter((value) => value === "N").length ?? 0;

  const mutationLabel = (column: number, from: string, to: string) => {
    if (mode === "nt") return `${column + 1} ${from}→${to}`;
    const aaColumn = Math.floor(column / 3);
    const codonOffset = column % 3 + 1;
    const position = numbering === "kabat" && kabat?.labels[aaColumn] ? `${kabat.chain}${kabat.labels[aaColumn]}` : String(aaColumn + 1);
    return `${position}.${codonOffset} ${from}→${to}`;
  };

  const toggleFullscreen = async () => {
    if (document.fullscreenElement) await document.exitFullscreen();
    else await containerRef.current?.requestFullscreen();
  };

  return <div ref={containerRef} className={`coordinated-tree-viewer${fullscreen ? " is-fullscreen" : ""}`}>
    <header className="coordinated-tree-header">
      <div><span className="section-kicker">{variantLabel}</span><h4>Tree, abundance and aligned leaf sequences</h4><p>{layout.leaves.toLocaleString()} tips · {layout.mode} · {selectedColumns.length.toLocaleString()} of {alignmentColumns.toLocaleString()} {mode === "nt" ? "nucleotide" : "amino-acid"} columns shown</p></div>
      <div className="tree-header-actions"><button type="button" onClick={() => void toggleFullscreen()}>{fullscreen ? "Exit full screen" : "Full screen"}</button><button type="button" onClick={() => saveSvg(svgRef.current, `${name}.svg`)}>SVG with alignment ↓</button><button type="button" onClick={() => downloadText(`${newick.replace(/;?$/, ";")}\n`, `${name}.nwk`, "text/plain;charset=utf-8")}>Newick</button></div>
    </header>
    <div className="coordinated-tree-controls">
      <div className="mode-toggle"><button className={mode === "nt" ? "active" : ""} type="button" onClick={() => onModeChange("nt")}>Nucleotide</button><button className={mode === "aa" ? "active" : ""} type="button" onClick={() => onModeChange("aa")}>Amino acid</button></div>
      <div className="mode-toggle"><button className={layoutMode === "phylogram" ? "active" : ""} type="button" onClick={() => setLayoutMode("phylogram")}>Branch lengths</button><button className={layoutMode === "cladogram" ? "active" : ""} type="button" onClick={() => setLayoutMode("cladogram")}>Topology</button></div>
      <label className="check-line"><input type="checkbox" checked={showMutations} disabled={!parsimony} onChange={(event) => setShowMutations(event.target.checked)} /><span>Branch mutations</span></label>
      {mode === "aa" && <div className="mode-toggle"><button className={numbering === "alignment" ? "active" : ""} type="button" onClick={() => setNumbering("alignment")}>Alignment positions</button><button className={numbering === "kabat" ? "active" : ""} type="button" disabled={isTcr} title={isTcr ? "Kabat numbering is defined for IGH, IGK and IGL, not TCR chains." : "Number IG variable domains with Kabat positions"} onClick={() => setNumbering("kabat")}>Kabat</button></div>}
    </div>
    <div className="coordinated-tree-options">
      <label><span>Alignment region</span><select value={preset} onChange={(event) => { const value = event.target.value as RegionPreset | "motifs"; setPreset(value); if (value === "motifs") setColorMode("motif"); }}><option value="full">Full alignment</option><option value="variable">Annotated variable domain</option><option value="cdrs">All CDRs</option><option value="fwr1">FWR1</option><option value="cdr1">CDR1</option><option value="fwr2">FWR2</option><option value="cdr2">CDR2</option><option value="fwr3">FWR3</option><option value="cdr3">CDR3</option><option value="fwr4">FWR4</option><option value="motifs">Motif-matched columns</option><option value="custom">Positions / ranges</option></select></label>
      {preset === "custom" && <label className="wide-option"><span>{numbering === "kabat" && mode === "aa" ? "Kabat positions or ranges" : "1-based columns or ranges"}</span><input value={customColumns} onChange={(event) => setCustomColumns(event.target.value)} placeholder={numbering === "kabat" && mode === "aa" ? "31-35B, 50, 95-102" : "95-110, 145, 220-235"} /></label>}
      <label><span>Cell colors</span><select value={colorMode} onChange={(event) => setColorMode(event.target.value as "residue" | "motif")}><option value="residue">Standard residue colors</option><option value="motif">Motif colors</option></select></label>
      {(colorMode === "motif" || preset === "motifs") && <label className="wide-option"><span>Motifs · comma/newline separated</span><input value={motifText} onChange={(event) => setMotifText(event.target.value)} placeholder={mode === "aa" ? "YYC, AR..Y, WGQ" : "TGT, AAR, TGG"} /></label>}
      <label><span>Tree width · {treeWidth}px</span><input type="range" min="320" max="1400" step="20" value={treeWidth} onChange={(event) => setTreeWidth(Number(event.target.value))} /></label>
      <label><span>Vertical scale · {rowHeight}px</span><input type="range" min="16" max="64" step="2" value={rowHeight} onChange={(event) => setRowHeight(Number(event.target.value))} /></label>
      <label><span>Residue width · {cellWidth}px</span><input type="range" min="9" max="24" step="1" value={cellWidth} onChange={(event) => setCellWidth(Number(event.target.value))} /></label>
    </div>
    {kabatStatus && <div className={numbering === "kabat" ? "viewer-numbering-status" : "viewer-numbering-status warning"}>{kabatStatus}</div>}
    {mode === "aa" && numbering === "kabat" && kabat && <div className="viewer-numbering-status">Kabat {kabat.chain} consensus · {kabat.numberedColumns.toLocaleString()} columns · {kabat.contributingSequences} members · mean confidence {kabat.confidence.toFixed(2)}</div>}
    <div className="coordinated-tree-legend">
      <span><i className="bubble-key" />tip bubble area increases with duplicate_count (maximum {maximumMultiplicity.toLocaleString()})</span>
      {parsimony && <span><i className="mutation-key" />germline-constrained nucleotide parsimony · {parsimony.score.toLocaleString()} steps · {inferredN.toLocaleString()} germline N sites inferred</span>}
      {variant === "stable" && <span>{collapsedEdges.toLocaleString()} internal edges ≤ {collapseThreshold.toExponential(0)} shown as polytomies</span>}
      {colorMode === "motif" && motifs.map((motif, index) => <span key={`${motif}-${index}`}><i style={{ background: MOTIF_COLORS[index % MOTIF_COLORS.length] }} />{index + 1}. {motif}</span>)}
    </div>
    <div className="coordinated-tree-scroll"><svg ref={svgRef} width={svgWidth} height={svgHeight} viewBox={`0 0 ${svgWidth} ${svgHeight}`} role="img" aria-label="Rooted lineage phylogeny coordinated with aligned leaf sequences">
      <rect width={svgWidth} height={svgHeight} fill="#fbfaf5" />
      <text x="12" y="18" fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700" fill="#31413c">{variantLabel} · {mode === "nt" ? "nucleotide" : numbering === "kabat" ? `amino acid · Kabat ${kabat?.chain ?? ""}` : "amino acid"}</text>
      {runs.map((run, index) => {
        const x = matrixX + run.start * cellWidth;
        const width = (run.end - run.start + 1) * cellWidth;
        return <g key={`${run.start}-${index}`}><rect x={x} y="24" width={width} height="17" fill={run.region ? REGION_COLORS[run.region] : "#eceeea"} /><text x={x + width / 2} y="36" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="7" fontWeight="700" fill="#364640">{run.region && width >= 25 ? REGION_LABELS[run.region] : ""}</text></g>;
      })}
      {selectedColumns.map((column, displayIndex) => {
        const label = coordinateLabels[column] || String(column + 1);
        const show = displayIndex === 0 || displayIndex === selectedColumns.length - 1 || displayIndex % Math.max(1, Math.ceil(45 / cellWidth)) === 0 || selectedColumns[displayIndex - 1] + 1 !== column;
        return show ? <text key={`tick-${column}`} transform={`translate(${matrixX + displayIndex * cellWidth + cellWidth / 2},67) rotate(-55)`} textAnchor="end" fontFamily="ui-monospace,monospace" fontSize="7" fill="#65736e">{label}</text> : null;
      })}
      <line x1={treeWidth} x2={treeWidth} y1="23" y2={svgHeight - 10} stroke="#c9d0ca" />
      {layout.edges.map((edge, index) => {
        const mutations = parsimony?.mutationsByClade.get(cladeSignature(edge.child))?.filter((mutation) => mode === "nt" ? selectedColumnSet.has(mutation.column) : selectedColumnSet.has(Math.floor(mutation.column / 3))) ?? [];
        const label = mutations.slice(0, 2).map((mutation) => mutationLabel(mutation.column, mutation.from, mutation.to)).join(" · ") + (mutations.length > 2 ? ` · +${mutations.length - 2}` : "");
        const labelX = edge.parent.x + Math.max(2, edge.child.x - edge.parent.x) / 2;
        return <g key={`edge-${index}`}><path d={`M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}`} stroke="#49605a" strokeWidth="1.2" fill="none" />{showMutations && label && <g><rect x={labelX - 3} y={edge.child.y - 14} width={Math.max(28, label.length * 5.2)} height="12" rx="2" fill="#fffaf0" stroke="#d9cba9" /><text x={labelX} y={edge.child.y - 5} fontFamily="ui-monospace,monospace" fontSize="7" fontWeight="700" fill="#71374a">{label}</text></g>}</g>;
      })}
      {layout.nodes.map((node, index) => {
        const leaf = !node.children.length;
        const ordinal = ordinalFromAlignmentName(node.name);
        const multiplicity = ordinal === null ? 1 : Math.max(1, multiplicityByOrdinal.get(ordinal) ?? 1);
        const radius = leaf ? 3 + Math.min(10, Math.log2(multiplicity + 1) * 1.45) : 2.2;
        const recordSequence = displayedByName.get(node.name) ?? "";
        const motifsForRecord = motifByName.get(node.name);
        return <g key={`node-${index}`}>
          {leaf && <><line x1={node.x + radius + 2} x2={treeWidth + 6} y1={node.y} y2={node.y} stroke="#a7b2ad" strokeWidth="0.8" strokeDasharray="2 3" /><text x={treeWidth + 10} y={node.y + 3.2} fontFamily="ui-monospace,monospace" fontSize="8" fontWeight={node.name === GERMLINE_OUTGROUP ? "700" : "500"} fill="#263630">{shortName(node.name)}<title>{node.name} · multiplicity {multiplicity}</title></text></>}
          <circle cx={node.x} cy={node.y} r={radius} fill={node.name === GERMLINE_OUTGROUP ? "#d49a19" : leaf ? "#08796f" : "#49605a"} fillOpacity={leaf ? 0.84 : 1} stroke={leaf && multiplicity > 1 ? "#133f39" : "none"}><title>{leaf ? `${node.name} · multiplicity ${multiplicity}` : `Internal node · ${cladeSignature(node).split("\u0000").length} descendants`}</title></circle>
          {leaf && selectedColumns.map((column, displayIndex) => {
            const value = recordSequence[column] ?? "-";
            const motif = motifsForRecord?.[column] ?? 0;
            const background = colorMode === "motif" ? motif ? MOTIF_COLORS[(motif - 1) % MOTIF_COLORS.length] : value === "-" ? "#ffffff" : "#eceeea" : sequenceColor(value, mode);
            const x = matrixX + displayIndex * cellWidth;
            return <g key={`${node.name}-${column}`}><rect x={x} y={node.y - rowHeight * 0.4} width={cellWidth} height={rowHeight * 0.8} fill={background} stroke="#ffffff" strokeWidth="0.35" /><text x={x + cellWidth / 2} y={node.y + Math.min(4, rowHeight * 0.16)} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize={Math.min(10, cellWidth - 2)} fontWeight="700" fill={textColor(background)}>{value}</text><title>{node.name} · {coordinateLabels[column] || column + 1} · {value}{motif ? ` · motif ${motif}` : ""}</title></g>;
          })}
        </g>;
      })}
    </svg></div>
    <p className="scientific-note"><span>i</span>{variant === "raw" ? "Mutation mapping is disabled for the unrooted raw tree. Choose a germline-rooted view to reconstruct the UCA and branch changes." : "Equal-cost nucleotide parsimony is reconstructed on demand. Known germline bases constrain the UCA; N bases are unknown A/C/G/T states inferred from descendants, and gaps are explicit indel states."} {isTcr ? "Kabat is not defined for TCR chains, so TCR amino-acid views use alignment positions." : "Kabat numbering is computed in-browser from a multi-member consensus of IG variable-domain assignments."}</p>
  </div>;
}
