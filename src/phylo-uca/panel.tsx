import { useEffect, useMemo, useRef, useState } from "react";

import { runFastTree } from "../biowasm-runtime.ts";
import { CommitNumberInput } from "../commit-number-input.tsx";
import { inspectAlignment } from "../alignment-provenance.ts";
import { GERMLINE_OUTGROUP } from "../lineage-alignment.ts";
import { LineageTreeViewer } from "../lineage-tree-viewer.tsx";
import { ProbabilityLogo, serializeProbabilityLogoSvg, type ProbabilityLogoAnnotation } from "../probability-logo.tsx";
import type { CompiledReferences } from "../reference-pack.ts";
import type { AirrDetailRow } from "../result-store.ts";
import { alignmentRegionMap, type AlignmentFrameOffset, type VariableRegion } from "../lineage-phylogeny.ts";
import { parseFasta } from "../post-analysis-core.ts";
import type { SampleColorMap } from "../sample-colors.ts";
import { defaultPhyloUcaOptions } from "./defaults.ts";
import { PhyloUcaHmmAnnotationTracks, serializePhyloUcaTrackLogoSvg, type PhyloUcaAnnotationColumnLayout } from "./hmm-annotation.tsx";
import { collapseAndOrderHmmAnnotationTracks } from "./hmm-annotation-model.ts";
import { aminoAcidUcaLogoColumns, codonUcaLogoColumns, nucleotideUcaLogoColumns } from "./logo.ts";
import { PHYLO_UCA_CODON_SYMBOLS, translatePhyloUcaCodonState } from "./codons.ts";
import { prepareObservedOnlyAlignment } from "./references.ts";
import { runPhyloUca } from "./runtime.ts";
import { PhyloUcaPlacementMap } from "./placement-map.tsx";
import type { PhyloUcaOptions, PhyloUcaProgress, PhyloUcaResult, PhyloUcaSavedState } from "./types.ts";

export type PhyloUcaPanelState = PhyloUcaSavedState;

interface Props {
  alignment: string;
  lineageRows: AirrDetailRow[];
  lineageIds: number[];
  lineageLabel: string;
  locus: string;
  references: CompiledReferences;
  inputName: string;
  frameOffset: AlignmentFrameOffset;
  isTcr: boolean;
  sampleColors: SampleColorMap;
  multiplicityByOrdinal: Map<number, number>;
  lineageByOrdinal: Map<number, number>;
  initialState?: PhyloUcaPanelState | null;
  onStateChange?: (state: PhyloUcaPanelState | null) => void;
}

function download(value: string, name: string, type = "text/plain;charset=utf-8") {
  const url = URL.createObjectURL(new Blob([value], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 500);
}

function stem(name: string): string {
  return name.replace(/(?:\.airr)?\.(?:tsv|csv|txt|fasta|fastq)(?:\.gz)?$/i, "") || "swig";
}

function phaseLabel(progress: PhyloUcaProgress | null): string {
  if (!progress) return "";
  const labels: Record<PhyloUcaProgress["phase"], string> = {
    references: "Reference hypotheses",
    "tree-messages": "Tree messages",
    "edge-screen": "Edge screen",
    "hmm-search": "Placement + HMM search",
    posterior: "UCA posterior",
    finalize: "Finalizing",
  };
  return labels[progress.phase];
}

function posteriorTsv(result: PhyloUcaResult): string {
  const header = ["alignment_column", "map_character", "map_probability", "entropy_bits", "segment", "call", "p_A", "p_C", "p_G", "p_T", "p_gap"];
  const rows = result.posterior.map((site) => [site.alignmentColumn, site.mapCharacter, site.mapProbability, site.entropyBits, site.segment, site.call ?? "", ...site.probabilities]);
  return [header, ...rows].map((row) => row.join("\t")).join("\n") + "\n";
}

function codonPosteriorTsv(result: PhyloUcaResult): string {
  const header = ["codon_index", "alignment_columns", "codon", "amino_acid", "probability", "is_map"];
  const rows = (result.codonPosterior ?? []).flatMap((codon) => codon.probabilities.map((probability, state) => [
    codon.codonIndex,
    codon.alignmentColumns.join(","),
    PHYLO_UCA_CODON_SYMBOLS[state],
    translatePhyloUcaCodonState(state),
    probability,
    Number(PHYLO_UCA_CODON_SYMBOLS[state] === codon.mapCodon),
  ]));
  return [header, ...rows].map((row) => row.join("\t")).join("\n") + "\n";
}

export function PhyloUcaPanel({ alignment, lineageRows, lineageIds, lineageLabel, locus, references, inputName, frameOffset, isTcr, sampleColors, multiplicityByOrdinal, lineageByOrdinal, initialState, onStateChange }: Props) {
  const fingerprint = useMemo(() => inspectAlignment(alignment).fingerprint, [alignment]);
  const initialFrame = initialState?.frameOffset ?? initialState?.result?.frameOffset ?? 0;
  const matchingInitial = initialState?.alignmentFingerprint === fingerprint && initialState.lineageIds.join(",") === lineageIds.join(",") && initialFrame === frameOffset ? initialState : null;
  const [options, setOptions] = useState<PhyloUcaOptions>(() => matchingInitial?.options ?? defaultPhyloUcaOptions());
  const [result, setResult] = useState<PhyloUcaResult | null>(() => matchingInitial?.result ?? null);
  const [progress, setProgress] = useState<PhyloUcaProgress | null>(null);
  const [treeStage, setTreeStage] = useState(false);
  const [error, setError] = useState("");
  const [running, setRunning] = useState(false);
  const [treeMode, setTreeMode] = useState<"nt" | "aa">("nt");
  const [logoMode, setLogoMode] = useState<"nt" | "codon" | "aa">("nt");
  const [annotationMode, setAnnotationMode] = useState<"viterbi" | "marginalized">("viterbi");
  const [annotationScrollLeft, setAnnotationScrollLeft] = useState(0);
  const abortRef = useRef<AbortController | null>(null);
  const resultRef = useRef<HTMLDivElement>(null);
  const logoRef = useRef<SVGSVGElement>(null);
  const annotationSvgRef = useRef<SVGSVGElement>(null);
  const annotationViewportRef = useRef<HTMLDivElement>(null);
  const annotationScrollRef = useRef<HTMLDivElement>(null);
  const adoptedResultSnapshotRef = useRef(matchingInitial?.result?.generatedAt ?? null);
  const identityKey = `${lineageIds.join(",")}|${fingerprint}|${frameOffset}`;
  const identityKeyRef = useRef(identityKey);
  useEffect(() => () => abortRef.current?.abort(), []);
  useEffect(() => {
    if (identityKeyRef.current === identityKey) return;
    identityKeyRef.current = identityKey;
    abortRef.current?.abort();
    const restoredFrame = initialState?.frameOffset ?? initialState?.result?.frameOffset ?? 0;
    const restored = initialState?.alignmentFingerprint === fingerprint && initialState.lineageIds.join(",") === lineageIds.join(",") && restoredFrame === frameOffset ? initialState : null;
    adoptedResultSnapshotRef.current = restored?.result?.generatedAt ?? null;
    setOptions(restored?.options ?? defaultPhyloUcaOptions());
    setResult(restored?.result ?? null);
    setProgress(null);
    setError("");
  }, [fingerprint, identityKey, initialState, lineageIds]);
  useEffect(() => {
    const incoming = matchingInitial;
    if (!incoming?.result) return;
    const incomingGeneratedAt = incoming.result.generatedAt;
    if (adoptedResultSnapshotRef.current === incomingGeneratedAt) return;
    adoptedResultSnapshotRef.current = incomingGeneratedAt;
    setOptions(incoming.options);
    setResult(incoming.result);
  }, [matchingInitial]);
  useEffect(() => {
    const currentResult = result && result.alignmentFingerprint === fingerprint && (result.frameOffset ?? 0) === frameOffset ? result : undefined;
    onStateChange?.({ lineageIds: [...lineageIds], alignmentFingerprint: fingerprint, frameOffset, options, result: currentResult });
  }, [fingerprint, frameOffset, lineageIds, onStateChange, options, result]);
  useEffect(() => {
    if (result && !(result.codonPosterior?.length) && logoMode !== "nt") setLogoMode("nt");
  }, [logoMode, result]);

  const updateHmm = <K extends keyof PhyloUcaOptions["hmm"]>(key: K, value: PhyloUcaOptions["hmm"][K]) => {
    setOptions((current) => ({ ...current, hmm: { ...current.hmm, [key]: value } }));
    setResult(null);
  };
  const updateCandidates = <K extends keyof PhyloUcaOptions["candidates"]>(key: K, value: PhyloUcaOptions["candidates"][K]) => {
    setOptions((current) => ({ ...current, candidates: { ...current.candidates, [key]: value } }));
    setResult(null);
  };
  const updateSearch = <K extends keyof PhyloUcaOptions["search"]>(key: K, value: PhyloUcaOptions["search"][K]) => {
    setOptions((current) => ({ ...current, search: { ...current.search, [key]: value } }));
    setResult(null);
  };
  const resetOptions = () => {
    setOptions(defaultPhyloUcaOptions());
    setResult(null);
    setError("");
  };

  async function run() {
    const controller = new AbortController();
    abortRef.current = controller;
    setRunning(true);
    setTreeStage(true);
    setProgress(null);
    setError("");
    setResult(null);
    try {
      const observed = prepareObservedOnlyAlignment(alignment, GERMLINE_OUTGROUP);
      const observedTree = await runFastTree(observed.fasta, "gtr", false);
      if (controller.signal.aborted) throw new DOMException("Cancelled", "AbortError");
      setTreeStage(false);
      const inference = await runPhyloUca({
        curatedAlignmentFasta: alignment,
        observedTreeNewick: observedTree.newick,
        observedAlignmentFasta: observed.posteriorFasta,
        retainedColumns: observed.posteriorColumns,
        germlineGuideName: GERMLINE_OUTGROUP,
        lineageRows: lineageRows.map((row) => ({ ordinal: row.record.ordinal, sequenceId: row.record.sequenceId, locus: row.values.locus || row.record.locus, values: { ...row.values } })),
        references: { V: references.V, D: references.D, J: references.J },
        locus,
        lineageLabel,
        alignmentFingerprint: fingerprint,
        frameOffset,
        options,
      }, setProgress, controller.signal);
      adoptedResultSnapshotRef.current = inference.generatedAt;
      setResult(inference);
      window.requestAnimationFrame(() => resultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (runError) {
      if ((runError as Error).name !== "AbortError") setError(runError instanceof Error ? runError.message : String(runError));
    } finally {
      if (abortRef.current === controller) abortRef.current = null;
      setRunning(false);
      setTreeStage(false);
    }
  }

  function downloadAnnotationSvg(visibleOnly: boolean) {
    const svg = annotationSvgRef.current;
    const logo = logoRef.current;
    if (!svg || !logo) return;
    let serialized: string;
    if (visibleOnly) {
      const viewport = annotationViewportRef.current;
      const trackScroll = annotationScrollRef.current;
      if (!viewport || !trackScroll) return;
      serialized = serializePhyloUcaTrackLogoSvg(svg, logo, {
        x: viewport.scrollLeft,
        y: trackScroll.scrollTop,
        width: viewport.clientWidth,
        height: trackScroll.clientHeight,
      });
    } else {
      serialized = serializePhyloUcaTrackLogoSvg(svg, logo);
    }
    download(serialized, `${base}.${annotationMode}.${visibleOnly ? "visible" : "full"}-hmm-tracks-and-logo.svg`, "image/svg+xml;charset=utf-8");
  }

  const progressFraction = treeStage ? 0.03 : progress?.total ? Math.max(0.04, Math.min(1, progress.processed / progress.total)) : 0;
  const base = `${stem(inputName)}.${lineageLabel.replace(/[^A-Za-z0-9_.-]+/g, "-")}.phylo-uca`;
  const displayAlignment = result ? `${result.observedAlignmentFasta.trim()}\n>phylo_UCA\n${result.retainedColumns.map((column) => result.mapAlignedSequence[column] ?? "-").join("")}\n` : "";
  const logoColumns = useMemo(() => !result ? [] : logoMode === "nt"
    ? nucleotideUcaLogoColumns(result.posterior)
    : logoMode === "codon"
      ? codonUcaLogoColumns(result.codonPosterior ?? [])
      : aminoAcidUcaLogoColumns(result.codonPosterior ?? []), [logoMode, result]);
  const logoColumnWidth = logoMode === "codon" ? 48 : 22;
  const logoLeftInset = result?.hmmAnnotations ? 220 : 38;
  const annotationColumns = useMemo<PhyloUcaAnnotationColumnLayout[]>(() => {
    if (!result) return [];
    if (logoMode === "nt") return result.posterior.map((site, index) => ({ alignmentColumn: site.alignmentColumn, x: index * logoColumnWidth, width: logoColumnWidth }));
    return (result.codonPosterior ?? []).flatMap((codon, codonIndex) => codon.alignmentColumns.map((alignmentColumn, position) => ({
      alignmentColumn,
      x: codonIndex * logoColumnWidth + position * logoColumnWidth / 3,
      width: logoColumnWidth / 3,
    })));
  }, [logoColumnWidth, logoMode, result]);
  const annotationTracks = useMemo(() => collapseAndOrderHmmAnnotationTracks(result?.hmmAnnotations?.[annotationMode] ?? [], annotationMode), [annotationMode, result]);
  const logoBottomAnnotations = useMemo<ProbabilityLogoAnnotation[]>(() => {
    if (!result) return [];
    const nucleotideRegions = alignmentRegionMap(parseFasta(alignment, true), lineageRows);
    const displayed: Array<VariableRegion | null> = logoMode === "nt"
      ? result.posterior.map((site) => nucleotideRegions[site.alignmentColumn - 1] ?? null)
      : (result.codonPosterior ?? []).map((codon) => {
        const votes = new Map<VariableRegion, number>();
        codon.alignmentColumns.forEach((column) => {
          const region = nucleotideRegions[column - 1];
          if (region) votes.set(region, (votes.get(region) ?? 0) + 1);
        });
        return [...votes.entries()].sort((left, right) => right[1] - left[1])[0]?.[0] ?? null;
      });
    const colors: Record<string, string> = { cdr1: "#d8c3eb", cdr2: "#f0cd78", cdr3: "#ee9a87" };
    const spans: ProbabilityLogoAnnotation[] = [];
    let start = 0;
    while (start < displayed.length) {
      const region = displayed[start];
      let end = start;
      while (end + 1 < displayed.length && displayed[end + 1] === region) end += 1;
      if (region?.startsWith("cdr")) spans.push({ startColumn: start, endColumn: end, label: region.toUpperCase(), color: colors[region], title: `${region.toUpperCase()} AIRR annotation; display columns ${start + 1}–${end + 1}` });
      start = end + 1;
    }
    return spans;
  }, [alignment, lineageRows, logoMode, result]);
  const logoContentWidth = logoColumns.length * logoColumnWidth;
  return <section className="phylo-uca-panel">
    <header>
      <div><span className="section-kicker">Fixed-tree empirical Bayes</span><h4>Phylogenetic UCA inference</h4><p>Infer the UCA attachment, pendant length, recombination path, and nucleotide posterior from the exact curated lineage alignment. The ordinary germline guide is removed before the observed tree is inferred.</p></div>
      <a href="./PHYLO_UCA_INFERENCE.md" target="_blank" rel="noreferrer">Method details ↗</a>
    </header>
    <div className="phylo-uca-run-row">
      <div><strong>{options.characterMode === "auto" ? "Automatic GTR4 / internal-gap GTR5" : options.characterMode === "nucleotide-gtr4" ? "Forced nucleotide GTR4" : "Forced gap-aware GTR5"}</strong><span>Terminal tip gaps are missing coverage · broad V/J hypotheses · all {references.counts.D.toLocaleString()} active D records · up to {options.hmm.maximumDSegments} D segments · local placement marginalization {options.search.marginalizeLocally ? "on" : "off"}</span></div>
      {!running ? <button className="post-primary dark" type="button" onClick={() => void run()}>Infer phylogenetic UCA</button> : <button type="button" onClick={() => abortRef.current?.abort()}>Cancel</button>}
    </div>
    <details className="post-advanced phylo-uca-advanced">
      <summary><span>Advanced UCA model and search settings</span><small>Candidate breadth, recombination priors, repeated D, gap model, and placement integration</small></summary>
      <div className="phylo-uca-advanced-body">
        <fieldset><legend>Character model</legend><label title="Auto uses four states unless a gap occurs between a tip's first and last observed nucleotide. Leading and trailing gap padding is always missing data."><span>Alignment character model</span><select value={options.characterMode} onChange={(event) => { setOptions((current) => ({ ...current, characterMode: event.target.value as PhyloUcaOptions["characterMode"] })); setResult(null); }}><option value="auto">Automatic · GTR4 unless internal gaps occur</option><option value="nucleotide-gtr4">Force nucleotide GTR4 · all gaps missing</option><option value="gap-aware-gtr5">Force GTR5 · internal gaps only</option></select></label><label title="Stationary frequency of the explicit internal alignment-gap character. Leading and trailing tip gaps remain missing. Ignored by GTR4."><span>Gap equilibrium frequency</span><CommitNumberInput min="0.0001" max="0.5" step="0.005" value={options.model.frequencies[4]} onCommit={(value) => { setOptions((current) => ({ ...current, model: { ...current.model, frequencies: [current.model.frequencies[0], current.model.frequencies[1], current.model.frequencies[2], current.model.frequencies[3], value] } })); setResult(null); }} /></label></fieldset>
        <fieldset><legend>Candidate hypotheses</legend><label title="Retain V alleles no more than this many fixed-alignment differences beyond the best candidate, in addition to every observed call hypothesis."><span>V extra-difference window</span><CommitNumberInput min="0" max="30" value={options.candidates.vMaximumExtraDifferences} onCommit={(value) => updateCandidates("vMaximumExtraDifferences", value)} /></label><label title="Equivalent broad-screen window for J candidates."><span>J extra-difference window</span><CommitNumberInput min="0" max="20" value={options.candidates.jMaximumExtraDifferences} onCommit={(value) => updateCandidates("jMaximumExtraDifferences", value)} /></label><label title="Computational guard after observed V hypotheses have been retained."><span>Maximum V candidates</span><CommitNumberInput min="1" max="250" value={options.candidates.maximumVCandidates} onCommit={(value) => updateCandidates("maximumVCandidates", value)} /></label><label title="Computational guard after observed J hypotheses have been retained."><span>Maximum J candidates</span><CommitNumberInput min="1" max="100" value={options.candidates.maximumJCandidates} onCommit={(value) => updateCandidates("maximumJCandidates", value)} /></label></fieldset>
        <fieldset><legend>Recombination HMM</legend><label title="The automaton can use zero through this many D segments. Values above one admit rare VDDJ and higher-order hypotheses."><span>Maximum D segments</span><CommitNumberInput min="0" max="5" value={options.hmm.maximumDSegments} onCommit={(value) => updateHmm("maximumDSegments", value)} /></label><label title="Minimum number of consecutive templated D nucleotides before a D path may exit."><span>Minimum D match</span><CommitNumberInput min="1" max="12" value={options.hmm.minimumDMatch} onCommit={(value) => updateHmm("minimumDMatch", value)} /></label><label title="Prior probability weight for entering another D after one D has ended. Any probability from 0 through 1 is accepted; the default is 0.015."><span>Additional-D probability</span><CommitNumberInput min="0" max="1" step="0.005" value={options.hmm.additionalDProbability} onCommit={(value) => updateHmm("additionalDProbability", value)} /></label><label title="Expected run length of non-templated nucleotides before another transition is attempted."><span>Mean N length</span><CommitNumberInput min="0" max="40" step="0.5" value={options.hmm.meanNLength} onCommit={(value) => updateHmm("meanNLength", value)} /></label><label title="Small leakage probability away from a deterministic reference nucleotide; this is not an SHM context model."><span>Template leakage</span><CommitNumberInput min="0.000001" max="0.25" step="0.001" value={options.hmm.templateMismatchProbability} onCommit={(value) => updateHmm("templateMismatchProbability", value)} /></label><label title="Prior gap probability in non-templated junction states; ignored by GTR4."><span>Junction gap probability</span><CommitNumberInput min="0.000001" max="0.5" step="0.005" value={options.hmm.junctionGapProbability} onCommit={(value) => updateHmm("junctionGapProbability", value)} /></label></fieldset>
        <fieldset><legend>Placement search</legend>
          <label title="The default independently mixes retained V/J allele nucleotides at each fixed-alignment column. This is only a fast starting-position screen; retained points are always recomputed with the complete recombination HMM."><span>Starting-position screen</span><select value={options.search.screenMode ?? "vj-mixture"} onChange={(event) => updateSearch("screenMode", event.target.value as PhyloUcaOptions["search"]["screenMode"])}><option value="vj-mixture">V/J nucleotide mixture · default</option><option value="germline-guide">Single N-masked germline guide</option></select></label>
          <label title="Number of attachment positions screened along every observed-tree edge before full-HMM refinement. Endpoints and at least one branch-interior point are always included."><span>Screen points / edge</span><CommitNumberInput min="3" max="101" value={options.search.screenEdgeGridPoints ?? 5} onCommit={(value) => updateSearch("screenEdgeGridPoints", Math.max(3, Math.floor(value)))} /></label>
          <label title="Top screen-ranked edges receiving the full recombination HMM. Set to zero to search every edge."><span>Full-HMM edges (0 = all)</span><CommitNumberInput min="0" max="10000" value={options.search.fullHmmEdges} onCommit={(value) => updateSearch("fullHmmEdges", Math.max(0, Math.floor(value)))} /></label>
          <label title="Full-HMM attachment grid size on each retained edge. Endpoints and at least one interior branch point are always included."><span>Full-HMM points / edge</span><CommitNumberInput min="3" max="101" value={options.search.edgeGridPoints} onCommit={(value) => updateSearch("edgeGridPoints", Math.max(3, Math.floor(value)))} /></label>
          <label title="Full-HMM grid size for the UCA pendant branch length."><span>UCA branch grid points</span><CommitNumberInput min="2" max="101" value={options.search.branchGridPoints} onCommit={(value) => updateSearch("branchGridPoints", Math.max(2, Math.floor(value)))} /></label>
          <label title="Full-HMM local refinement rounds around the leading placement points."><span>Full-HMM refinement rounds</span><CommitNumberInput min="0" max="10" value={options.search.localRefinementRounds} onCommit={(value) => updateSearch("localRefinementRounds", Math.max(0, Math.floor(value)))} /></label>
          <label title="Maximum UCA-to-observed-tree branch length in expected substitutions per character."><span>Maximum UCA branch</span><CommitNumberInput min="0.001" max="3" step="0.01" value={options.search.maximumUcaBranchLength} onCommit={(value) => updateSearch("maximumUcaBranchLength", value)} /></label>
          <label title="Mean of the fixed exponential prior on UCA pendant length."><span>Branch-prior mean</span><CommitNumberInput min="0.0001" max="2" step="0.01" value={options.search.branchPriorMean} onCommit={(value) => updateSearch("branchPriorMean", value)} /></label>
          <label title="Number of high-scoring full-HMM branch points integrated in the reported marginal nucleotide probabilities."><span>Local posterior points</span><CommitNumberInput min="1" max="500" value={options.search.localPosteriorPoints} onCommit={(value) => updateSearch("localPosteriorPoints", Math.max(1, Math.floor(value)))} /></label>
          <label className="check-line" title="Average nucleotide probabilities over nearby attachment/length hypotheses rather than conditioning only on the optimum."><input type="checkbox" checked={options.search.marginalizeLocally} onChange={(event) => updateSearch("marginalizeLocally", event.target.checked)} /><span>Marginalize nearby full-HMM placements</span></label>
        </fieldset>
      </div>
      <div className="phylo-uca-default-reset"><div><strong>Restore calibrated starting values</strong><span>Resets every UCA character, candidate, HMM, and placement option; it does not alter the curated alignment.</span></div><button type="button" onClick={resetOptions}>Reset all UCA settings to defaults</button></div>
    </details>
    {running && <div className="phylo-uca-progress" role="status"><div><strong>{treeStage ? "Observed-only FastTree" : phaseLabel(progress)}</strong><span>{treeStage ? "Removing the germline guide and inferring the fixed tree" : progress?.detail}</span></div><progress max="1" value={progressFraction} /><b>{treeStage ? "tree" : progress?.total ? `${progress.processed.toLocaleString()} / ${progress.total.toLocaleString()}` : "working"}</b></div>}
    {error && <div className="inline-method-error" role="alert"><strong>Phylogenetic UCA inference stopped</strong><span>{error}</span><button type="button" onClick={() => setError("")}>Dismiss</button></div>}
    {result && <div ref={resultRef} className="phylo-uca-results">
      <div className="phylo-uca-result-head"><div><span className="section-kicker">Maximum joint path + exact marginal characters</span><h4>{result.mapVCall || "V?"} · {result.mapDCalls.length ? result.mapDCalls.join(" → ") : "no D"} · {result.mapJCall || "J?"}</h4><p>{result.characterModel === "nucleotide-gtr4" ? "Four-state nucleotide GTR" : "Five-state A/C/G/T/gap GTR"} · edge {result.bestPlacement.endpointA} ↔ {result.bestPlacement.endpointB} at {(result.bestPlacement.edgeFraction * 100).toFixed(1)}% · UCA branch {result.bestPlacement.ucaBranchLength.toFixed(5)} · effective local placements {result.effectivePlacementCount.toFixed(2)}</p></div><div className="result-actions"><button type="button" onClick={() => download(`>phylo_UCA_joint_MAP_aligned\n${result.mapAlignedSequence}\n>phylo_UCA_marginal_consensus_aligned\n${result.posteriorConsensusAligned}\n`, `${base}.fasta`)}>UCA FASTA ↓</button><button type="button" onClick={() => download(posteriorTsv(result), `${base}.nucleotide-posterior.tsv`)}>Nucleotide TSV ↓</button>{Boolean(result.codonPosterior?.length) && <button type="button" onClick={() => download(codonPosteriorTsv(result), `${base}.codon-posterior.tsv`)}>Codon TSV ↓</button>}<button type="button" onClick={() => download(result.placedTreeNewick, `${base}.placed-tree.nwk`)}>Placed tree ↓</button><button type="button" onClick={() => download(JSON.stringify(result, null, 2), `${base}.json`, "application/json")}>Complete JSON ↓</button></div></div>
      <div className="phylo-uca-stats"><article><span>Observed tree</span><strong>{result.observedAlignmentFasta.split("\n>").length.toLocaleString()} tips</strong></article><article><span>Candidate set</span><strong>{result.candidateReport.v.length} V · {result.candidateReport.d.length} D · {result.candidateReport.j.length} J</strong></article><article><span>Placement log marginal</span><strong>{result.logMarginalLikelihood.toFixed(2)}</strong></article><article><span>Runtime</span><strong>{(result.elapsedMs / 1000).toFixed(2)} s</strong></article></div>
      <PhyloUcaPlacementMap newick={result.observedTreeNewick} placements={result.placements} title={`${lineageLabel} UCA local placement likelihood surface`} />
      <div className="phylo-uca-sequence"><strong>Joint MAP aligned UCA</strong><code>{result.mapAlignedSequence}</code><strong>Marginal consensus</strong><code>{result.posteriorConsensusAligned}</code></div>
      <section className="phylo-uca-logo">
        <header><div><span className="section-kicker">Marginal character probabilities</span><h5>UCA posterior frequency logo</h5><p>Every stack has height 1. Letter height is posterior frequency; entropy does not rescale the column.</p></div><div><div className="mode-toggle"><button className={logoMode === "nt" ? "active" : ""} type="button" onClick={() => setLogoMode("nt")}>Nucleotide</button><button className={logoMode === "codon" ? "active" : ""} type="button" disabled={!result.codonPosterior?.length} onClick={() => setLogoMode("codon")}>Codon</button><button className={logoMode === "aa" ? "active" : ""} type="button" disabled={!result.codonPosterior?.length} onClick={() => setLogoMode("aa")}>Amino acid</button></div><button type="button" onClick={() => logoRef.current && download(serializeProbabilityLogoSvg(logoRef.current), `${base}.${logoMode}-posterior-logo.svg`, "image/svg+xml;charset=utf-8")}>Logo SVG ↓</button></div></header>
        {result.hmmAnnotations && <div className="phylo-uca-annotation-toolbar"><div className="mode-toggle"><button className={annotationMode === "viterbi" ? "active" : ""} type="button" onClick={() => setAnnotationMode("viterbi")}>Best path</button><button className={annotationMode === "marginalized" ? "active" : ""} type="button" onClick={() => setAnnotationMode("marginalized")}>Marginalized</button></div><div className="phylo-uca-track-legend" aria-label="Track background colors"><span className="v">V</span><span className="n">NT</span><span className="d">D</span><span className="j">J</span></div><p>{annotationMode === "viterbi"
          ? `One Viterbi recombination path at the single highest-scoring tree placement. Identical D alleles are combined across alignment registers and D-use ordinals; disagreements become nucleotide mixtures. Rows follow weighted left-to-right source position. ${annotationTracks.length} display tracks.`
          : `Forward–backward source occupancy summed over HMM paths and retained tree placements. Every non-template and unresolved-boundary route is summed into the single NT mixture at the top. V, D, and J allele rows then follow in segment order and are sorted by their posterior-mass center from left to right. ${annotationTracks.length} display tracks; allele groups appear at ≥ ${(result.hmmAnnotations.minimumDisplayedWeight * 100).toFixed(1)}%.${result.hmmAnnotations.omittedMarginalTrackCount ? ` ${result.hmmAnnotations.omittedMarginalTrackCount} subthreshold source track${result.hmmAnnotations.omittedMarginalTrackCount === 1 ? " is" : "s are"} hidden before display aggregation.` : ""}`}</p><div className="phylo-uca-track-downloads"><button type="button" onClick={() => downloadAnnotationSvg(false)}>Full tracks + logo SVG ↓</button><button type="button" onClick={() => downloadAnnotationSvg(true)}>Visible tracks + logo SVG ↓</button></div></div>}
        <div ref={annotationViewportRef} className="phylo-uca-posterior-viewport" onScroll={(event) => setAnnotationScrollLeft(event.currentTarget.scrollLeft)}>
          {result.hmmAnnotations && <div ref={annotationScrollRef} className="phylo-uca-annotation-scroll" style={{ width: logoLeftInset + logoContentWidth }}><PhyloUcaHmmAnnotationTracks ref={annotationSvgRef} tracks={annotationTracks} columns={annotationColumns} leftInset={logoLeftInset} contentWidth={logoContentWidth} labelOffset={annotationScrollLeft} title={`${lineageLabel} ${annotationMode === "viterbi" ? "best-path" : "marginalized"} phylo-HMM source annotation`} /></div>}
          <ProbabilityLogo ref={logoRef} embedded leftInset={logoLeftInset} columns={logoColumns} alphabet={logoMode === "nt" ? "nucleotide" : logoMode === "aa" ? "amino-acid" : "custom"} title={`${lineageLabel} phylogenetic UCA ${logoMode === "nt" ? "nucleotide" : logoMode === "codon" ? "codon" : "amino-acid"} posterior`} labelEvery={logoMode === "nt" ? 5 : 2} columnWidth={logoColumnWidth} bottomAnnotations={logoBottomAnnotations} />
        </div>
        {logoMode !== "nt" && <p className="phylo-uca-logo-note">Exact joint codon probabilities are computed in reading frame {(result.frameOffset ?? frameOffset) + 1} by summing over HMM state paths, V/D/J candidate identities, and retained placement hypotheses. The amino-acid view sums synonymous codon states; it does not multiply nucleotide marginals. A complete gap codon is “-” and a partly gapped codon is “X”.</p>}
        {!result.codonPosterior?.length && <p className="phylo-uca-logo-note">This result predates exact codon posteriors. Re-run phylogenetic UCA inference to enable codon and amino-acid marginals.</p>}
        {!result.hmmAnnotations && <p className="phylo-uca-logo-note">This saved result predates HMM-source annotation tracks. Re-run phylogenetic UCA inference to add best-path and marginalized V/D/J/N tracks.</p>}
      </section>
      {result.warnings.map((warning, index) => <div className="scientific-note" key={index}><span>i</span><p>{warning}</p></div>)}
      <div className="tree-output-switch"><div className="mode-toggle"><button className={treeMode === "nt" ? "active" : ""} type="button" onClick={() => setTreeMode("nt")}>Nucleotide</button><button className={treeMode === "aa" ? "active" : ""} type="button" onClick={() => setTreeMode("aa")}>Amino acid</button></div><span>Tree rooted at the inferred UCA; its named sequence carrier is zero length and the complete inferred branch lies on the observed-tree side.</span></div>
      <LineageTreeViewer newick={result.placedTreeNewick} alignmentFasta={displayAlignment} rows={lineageRows} multiplicityByOrdinal={multiplicityByOrdinal} sampleColors={sampleColors} lineageByOrdinal={lineageByOrdinal} variant="rooted" collapsedEdges={0} collapseThreshold={0} mode={treeMode} onModeChange={setTreeMode} frameOffset={frameOffset} isTcr={isTcr} name={`${base}.tree`} />
    </div>}
  </section>;
}
