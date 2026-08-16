import { useEffect, useMemo, useRef, useState } from "react";

import { runFastTreeTask } from "../biowasm-task-runtime.ts";
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
import { parseReferenceFasta, prepareObservedOnlyAlignment } from "./references.ts";
import { runPhyloUca } from "./runtime.ts";
import { PhyloUcaPlacementMap } from "./placement-map.tsx";
import { phyloUcaBranchLengthGrid } from "./search-grid.ts";
import { phyloUcaPriorPredictiveSummary } from "./prior-predictive.ts";
import type { PhyloUcaMcmcDiagnostics, PhyloUcaOptions, PhyloUcaProgress, PhyloUcaResult, PhyloUcaSavedState } from "./types.ts";

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

function completePhyloUcaOptions(saved?: PhyloUcaOptions): PhyloUcaOptions {
  const defaults = defaultPhyloUcaOptions();
  if (!saved) return defaults;
  const savedHmm = (saved.hmm ?? {}) as Partial<PhyloUcaOptions["hmm"]>;
  const migratedHmm: PhyloUcaOptions["hmm"] = {
    ...defaults.hmm,
    ...savedHmm,
    dThreePrimeTrimContinuation: savedHmm.dThreePrimeTrimContinuation
      ?? (savedHmm.dExitProbability === undefined ? defaults.hmm.dThreePrimeTrimContinuation : 1 - savedHmm.dExitProbability),
    terminalPaddingGapProbability: savedHmm.terminalPaddingGapProbability
      ?? savedHmm.unknownTemplateGapProbability
      ?? defaults.hmm.terminalPaddingGapProbability,
  };
  return {
    ...defaults,
    ...saved,
    model: { ...defaults.model, ...saved.model },
    candidates: { ...defaults.candidates, ...saved.candidates },
    hmm: migratedHmm,
    search: { ...defaults.search, ...saved.search },
  };
}

function phaseLabel(progress: PhyloUcaProgress | null): string {
  if (!progress) return "";
  const labels: Record<PhyloUcaProgress["phase"], string> = {
    references: "Reference hypotheses",
    "tree-messages": "Tree messages",
    "edge-screen": "Edge screen",
    "hmm-search": "Placement + HMM search",
    mcmc: "Gibbs/MH posterior sampling",
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

function acceptance(accepted: number, proposed: number): string {
  return proposed ? `${(100 * accepted / proposed).toFixed(1)}%` : "—";
}

function tracePath(values: readonly number[], left: number, top: number, width: number, height: number): string {
  if (!values.length) return "";
  const minimum = Math.min(...values);
  const maximum = Math.max(...values);
  const span = Math.max(1e-12, maximum - minimum);
  return values.map((value, index) => `${index ? "L" : "M"} ${left + width * index / Math.max(1, values.length - 1)} ${top + height * (1 - (value - minimum) / span)}`).join(" ");
}

function PhyloUcaMcmcMixing({ diagnostics }: { diagnostics: PhyloUcaMcmcDiagnostics }) {
  const trace = diagnostics.trace;
  if (!trace.length) return null;
  const branchValues = trace.map((point) => point.ucaBranchLength);
  const marginalValues = trace.map((point) => point.logMarginalLikelihood);
  const width = 920;
  const left = 46;
  const plotWidth = width - left - 18;
  const burnX = left + plotWidth * Math.max(0, diagnostics.burnIn - 1) / Math.max(1, trace.length - 1);
  const samplingSeconds = Math.max(1e-9, (diagnostics.samplingMilliseconds ?? 0) / 1000);
  const average = (milliseconds: number | undefined, count: number | undefined) => count && milliseconds !== undefined ? `${(milliseconds / count).toFixed(1)} ms` : "—";
  const conditionalProposals = diagnostics.branchProposals + diagnostics.positionProposals + diagnostics.globalProposals;
  return <section className="phylo-uca-mcmc">
    <header><div><span className="section-kicker">Sampler diagnostics</span><h5>Gibbs/MH mixing</h5><p>Every iteration draws one exact joint recombination path and UCA sequence, then updates continuous attachment fraction and continuous pendant length. The vertical line marks burn-in.</p></div><div className="phylo-uca-mcmc-stats"><span>Branch MH <b>{acceptance(diagnostics.branchAccepted, diagnostics.branchProposals)}</b></span><span>Position MH <b>{acceptance(diagnostics.positionAccepted, diagnostics.positionProposals)}</b></span><span>Global jumps <b>{acceptance(diagnostics.globalAccepted, diagnostics.globalProposals)}</b></span><span>Collapsed refresh <b>{acceptance(diagnostics.collapsedAccepted ?? 0, diagnostics.collapsedProposals ?? 0)}</b></span><span>Edge switches <b>{diagnostics.edgeSwitches}</b></span><span>Branch ESS <b>{diagnostics.branchEffectiveSampleSize.toFixed(1)} / {diagnostics.retainedSamples}</b></span><span>Log-target ESS <b>{diagnostics.logTargetEffectiveSampleSize.toFixed(1)} / {diagnostics.retainedSamples}</b></span>{diagnostics.samplingMilliseconds !== undefined && <><span>Branch ESS/s <b>{(diagnostics.branchEffectiveSampleSize / samplingSeconds).toFixed(2)}</b></span><span>Log-target ESS/s <b>{(diagnostics.logTargetEffectiveSampleSize / samplingSeconds).toFixed(2)}</b></span><span>Full-HMM draw <b>{average(diagnostics.gibbsMilliseconds, diagnostics.gibbsDraws)}</b></span><span>Collapsed marginal <b>{average(diagnostics.collapsedMarginalMilliseconds, diagnostics.collapsedProposals)}</b></span><span>Fixed-UCA proposal <b>{average(diagnostics.conditionalMhMilliseconds, conditionalProposals)}</b></span></>}</div></header>
    <div className="phylo-uca-mcmc-chart"><svg viewBox={`0 0 ${width} 214`} role="img" aria-label="MCMC branch-length and full-HMM marginal-likelihood traces">
      <title>Continuous Gibbs/MH trace; branch length above and full-HMM marginal likelihood below</title>
      <rect x={left} y="16" width={plotWidth} height="78" fill="#f5f8f6" />
      <rect x={left} y="116" width={plotWidth} height="78" fill="#f5f8f6" />
      <line x1={burnX} x2={burnX} y1="10" y2="200" stroke="#d36d3f" strokeDasharray="4 3" />
      <path d={tracePath(branchValues, left, 20, plotWidth, 70)} fill="none" stroke="#08796f" strokeWidth="1.4" />
      <path d={tracePath(marginalValues, left, 120, plotWidth, 70)} fill="none" stroke="#845a9e" strokeWidth="1.4" />
      <text x="7" y="55">branch</text><text x="7" y="66">length</text>
      <text x="7" y="153">full-HMM</text><text x="7" y="164">log L</text>
      <text x={left} y="210">1</text><text x={left + plotWidth} y="210" textAnchor="end">{trace.length}</text>
      <text x={Math.min(left + plotWidth - 4, burnX + 4)} y="14" fill="#a14b2d">burn-in</text>
    </svg></div>
  </section>;
}

export function PhyloUcaPanel({ alignment, lineageRows, lineageIds, lineageLabel, locus, references, inputName, frameOffset, isTcr, sampleColors, multiplicityByOrdinal, lineageByOrdinal, initialState, onStateChange }: Props) {
  const fingerprint = useMemo(() => inspectAlignment(alignment).fingerprint, [alignment]);
  const initialFrame = initialState?.frameOffset ?? initialState?.result?.frameOffset ?? 0;
  const matchingInitial = initialState?.alignmentFingerprint === fingerprint && initialState.lineageIds.join(",") === lineageIds.join(",") && initialFrame === frameOffset ? initialState : null;
  const [options, setOptions] = useState<PhyloUcaOptions>(() => completePhyloUcaOptions(matchingInitial?.options));
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
    setOptions(completePhyloUcaOptions(restored?.options));
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
    setOptions(completePhyloUcaOptions(incoming.options));
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
    try {
      const observed = prepareObservedOnlyAlignment(alignment, GERMLINE_OUTGROUP);
      const observedTree = await runFastTreeTask(observed.fasta, "gtr", false, controller.signal);
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
  const displayedBranchGrid = useMemo(() => phyloUcaBranchLengthGrid(options.search), [options.search.branchGridPoints, options.search.maximumUcaBranchLength, options.search.minimumPositiveUcaBranchLength]);
  const priorPredictive = useMemo(() => phyloUcaPriorPredictiveSummary(
    options.hmm,
    parseReferenceFasta(references.D, locus).map((record) => record.sequence.length),
  ), [locus, options.hmm, references.D]);
  const resultInferenceMode = result?.options.search.inferenceMode ?? (result?.options.search.marginalizeLocally ? "grid-marginalization" : "maximum-likelihood");
  return <section className="phylo-uca-panel">
    <header>
      <div><span className="section-kicker">Fixed-tree empirical Bayes</span><h4>Phylogenetic UCA inference</h4><p>Infer the UCA attachment, pendant length, recombination path, and nucleotide posterior from the exact curated lineage alignment. The ordinary germline guide is removed before the observed tree is inferred.</p></div>
      <a href="./PHYLO_UCA_INFERENCE.md" target="_blank" rel="noreferrer">Method details ↗</a>
    </header>
    <div className="phylo-uca-run-row">
      <div><strong>{options.characterMode === "auto" ? "Automatic GTR4 / internal-gap GTR5" : options.characterMode === "nucleotide-gtr4" ? "Forced nucleotide GTR4" : "Forced gap-aware GTR5"}</strong><span>Terminal tip gaps are missing coverage · broad V/J hypotheses · all {references.counts.D.toLocaleString()} active D records · up to {options.hmm.maximumDSegments} D segments · {options.search.inferenceMode === "maximum-likelihood" ? "conditional ML placement" : options.search.inferenceMode === "grid-marginalization" ? "explicit grid marginalization" : "continuous Gibbs/MH"}</span></div>
      {!running ? <button className="post-primary dark" type="button" onClick={() => void run()}>Infer phylogenetic UCA</button> : <button type="button" onClick={() => abortRef.current?.abort()}>Cancel</button>}
    </div>
    <details className="post-advanced phylo-uca-advanced">
      <summary><span>Advanced UCA model and search settings</span><small>Candidate breadth, recombination priors, repeated D, gap model, and placement integration</small></summary>
      <div className="phylo-uca-advanced-body">
        <fieldset><legend>Character model</legend><label title="Auto uses four states unless a gap occurs between a tip's first and last observed nucleotide. Leading and trailing gap padding is always missing data."><span>Alignment character model</span><select value={options.characterMode} onChange={(event) => { setOptions((current) => ({ ...current, characterMode: event.target.value as PhyloUcaOptions["characterMode"] })); setResult(null); }}><option value="auto">Automatic · GTR4 unless internal gaps occur</option><option value="nucleotide-gtr4">Force nucleotide GTR4 · all gaps missing</option><option value="gap-aware-gtr5">Force GTR5 · internal gaps only</option></select></label><label title="Stationary frequency of the explicit internal alignment-gap character. Leading and trailing tip gaps remain missing. Ignored by GTR4."><span>Gap equilibrium frequency</span><CommitNumberInput min="0.0001" max="0.5" step="0.005" value={options.model.frequencies[4]} onCommit={(value) => { setOptions((current) => ({ ...current, model: { ...current.model, frequencies: [current.model.frequencies[0], current.model.frequencies[1], current.model.frequencies[2], current.model.frequencies[3], value] } })); setResult(null); }} /></label></fieldset>
        <fieldset><legend>Candidate hypotheses</legend><label title="Retain V alleles no more than this many fixed-alignment differences beyond the best candidate, in addition to every observed call hypothesis."><span>V extra-difference window</span><CommitNumberInput min="0" max="30" value={options.candidates.vMaximumExtraDifferences} onCommit={(value) => updateCandidates("vMaximumExtraDifferences", value)} /></label><label title="Equivalent broad-screen window for J candidates."><span>J extra-difference window</span><CommitNumberInput min="0" max="20" value={options.candidates.jMaximumExtraDifferences} onCommit={(value) => updateCandidates("jMaximumExtraDifferences", value)} /></label><label title="Computational guard after observed V hypotheses have been retained."><span>Maximum V candidates</span><CommitNumberInput min="1" max="250" value={options.candidates.maximumVCandidates} onCommit={(value) => updateCandidates("maximumVCandidates", value)} /></label><label title="Computational guard after observed J hypotheses have been retained."><span>Maximum J candidates</span><CommitNumberInput min="1" max="100" value={options.candidates.maximumJCandidates} onCommit={(value) => updateCandidates("maximumJCandidates", value)} /></label></fieldset>
        <fieldset className="phylo-uca-hmm-settings"><legend>Recombination HMM · every parameter</legend>
          <label title="The automaton can use zero through this many D segments. Values above one admit rare VDDJ and higher-order hypotheses."><span>Maximum D segments</span><CommitNumberInput min="0" max="5" value={options.hmm.maximumDSegments} onCommit={(value) => updateHmm("maximumDSegments", Math.max(0, Math.floor(value)))} /></label>
          <label title="Identifiability threshold: a D path must emit at least this many consecutive templated nucleotides before it may end."><span>Minimum D match</span><CommitNumberInput min="1" max="12" value={options.hmm.minimumDMatch} onCommit={(value) => updateHmm("minimumDMatch", Math.max(1, Math.floor(value)))} /></label>
          <label title="Prior probability that an identifiable first D segment is used. The default 0.934 is the human-IGH IGoR probability that at least five D nucleotides survive trimming."><span>First-D probability</span><CommitNumberInput min="0" max="1" step="0.001" value={options.hmm.initialDProbability} onCommit={(value) => updateHmm("initialDProbability", value)} /></label>
          <label title="P(add another D | at least one D was used). Endpoints 0 and 1 are accepted and remain stable after rerunning."><span>Additional-D probability</span><CommitNumberInput min="0" max="1" step="0.001" value={options.hmm.additionalDProbability} onCommit={(value) => updateHmm("additionalDProbability", value)} /></label>
          <label title="Probability that each V–D, D–D, D–J, or direct V–J boundary contains at least one non-templated base."><span>Non-empty N probability</span><CommitNumberInput min="0" max="1" step="0.001" value={options.hmm.junctionNProbability} onCommit={(value) => updateHmm("junctionNProbability", value)} /></label>
          <label title="Expected N-run length conditional on the run being non-empty. This is now the actual mean, not the previous off-by-one parameterization."><span>Mean positive N length</span><CommitNumberInput min="1" max="100" step="0.1" value={options.hmm.meanNLength} onCommit={(value) => updateHmm("meanNLength", value)} /></label>
          <label title="Conditional probability that a non-empty N insertion is exactly one nucleotide. Separating this mass prevents a geometric tail from grossly overproducing 1-nt N calls."><span>One-nt N probability</span><CommitNumberInput min="0" max="1" step="0.001" value={options.hmm.singleNProbability} onCommit={(value) => updateHmm("singleNProbability", value)} /></label>
          <label title="Number of geometric phases in the longer N-duration tail. Two gives a peaked, non-geometric distribution with only a few extra HMM states."><span>Long-N tail phases</span><CommitNumberInput min="1" max="4" step="1" value={options.hmm.nLengthPhases} onCommit={(value) => updateHmm("nLengthPhases", Math.max(1, Math.floor(value)))} /></label>
          <label title="Geometric deletion-tail ratio P(k+1)/P(k) at the V 3′ end. Larger values permit more V trimming."><span>V 3′ trim tail ratio</span><CommitNumberInput min="0" max="0.999999" step="0.001" value={options.hmm.vThreePrimeTrimContinuation} onCommit={(value) => updateHmm("vThreePrimeTrimContinuation", value)} /></label>
          <label title="Geometric deletion-tail ratio P(k+1)/P(k) at the D 5′ end. D alleles receive equal prior weight before this within-allele trim prior."><span>D 5′ trim tail ratio</span><CommitNumberInput min="0" max="0.999999" step="0.001" value={options.hmm.dFivePrimeTrimContinuation} onCommit={(value) => updateHmm("dFivePrimeTrimContinuation", value)} /></label>
          <label title="Geometric deletion-tail ratio P(k+1)/P(k) at the D 3′ end. This replaces the old per-retained-base D-exit penalty."><span>D 3′ trim tail ratio</span><CommitNumberInput min="0" max="0.999999" step="0.001" value={options.hmm.dThreePrimeTrimContinuation} onCommit={(value) => updateHmm("dThreePrimeTrimContinuation", value)} /></label>
          <label title="Geometric deletion-tail ratio P(k+1)/P(k) at the J 5′ end. J entry is possible only on a concrete projected J nucleotide."><span>J 5′ trim tail ratio</span><CommitNumberInput min="0" max="0.999999" step="0.001" value={options.hmm.jFivePrimeTrimContinuation} onCommit={(value) => updateHmm("jFivePrimeTrimContinuation", value)} /></label>
          {(["A", "C", "G", "T"] as const).map((base, index) => <label key={base} title={`Unnormalized ${base} weight in non-templated N sequence; all four weights are normalized together.`}><span>N-base {base} weight</span><CommitNumberInput min="0" max="1" step="0.001" value={options.hmm.nBaseFrequencies[index]} onCommit={(value) => updateHmm("nBaseFrequencies", options.hmm.nBaseFrequencies.map((current, entry) => entry === index ? value : current) as [number, number, number, number])} /></label>)}
          <label title="Direct mismatch leakage from a deterministic germline-template nucleotide. Zero is now exact: mutations between the UCA and observed tree are still modeled by the GTR branch."><span>Template leakage</span><CommitNumberInput min="0" max="0.25" step="0.001" value={options.hmm.templateMismatchProbability} onCommit={(value) => updateHmm("templateMismatchProbability", value)} /></label>
          <label title="Prior gap probability in a non-templated junction state; ignored by GTR4."><span>Junction gap probability</span><CommitNumberInput min="0" max="0.5" step="0.001" value={options.hmm.junctionGapProbability} onCommit={(value) => updateHmm("junctionGapProbability", value)} /></label>
          <label title="Gap prior only for leading V or trailing J alignment padding. It is never used inside the V–D–J junction."><span>Terminal-padding gap</span><CommitNumberInput min="0" max="0.5" step="0.001" value={options.hmm.terminalPaddingGapProbability} onCommit={(value) => updateHmm("terminalPaddingGapProbability", value)} /></label>
          <label title="Performance support around the observed V/J anchors in which D states are evaluated. This does not make a J state emit unresolved bases."><span>D-search flank columns</span><CommitNumberInput min="0" max="100" step="1" value={options.hmm.junctionSearchFlankColumns} onCommit={(value) => updateHmm("junctionSearchFlankColumns", Math.max(0, Math.floor(value)))} /></label>
        </fieldset>
        <section className="phylo-uca-prior-audit">
          <header><div><span className="section-kicker">Generative audit</span><h5>Prior-predictive recombination statistics</h5></div><p>{priorPredictive.draws.toLocaleString()} deterministic route/N draws with exact finite D trims against active D lengths; V/J rows show the exposed geometric-tail moments.</p></header>
          <div className="phylo-uca-prior-table"><table><thead><tr><th>Quantity</th><th>Mean</th><th>Median</th><th>90%</th><th>95%</th></tr></thead><tbody>{priorPredictive.metrics.map((metric) => <tr key={metric.id}><th>{metric.label}</th><td>{metric.mean.toFixed(2)}</td><td>{metric.median}</td><td>{metric.p90}</td><td>{metric.p95}</td></tr>)}</tbody></table></div>
          <p>The default shared trim tails moment-match the public human-IGH model learned by IGoR/OLGA: mean deletions V3 3.04, D5 6.01, D3 5.54, J5 6.74 nt; positive N runs average about 8.8 nt. Partis independently shows that real deletion distributions vary by allele and need categorical/tiered estimates when enough repertoire data exist; Swig deliberately does not fit that larger model here. These are human starting values, not universal constants. <a href="https://www.nature.com/articles/s41467-018-02832-w" target="_blank" rel="noreferrer">IGoR study ↗</a> · <a href="https://github.com/statbiophys/OLGA/tree/master/olga/default_models/human_B_heavy" target="_blank" rel="noreferrer">model files ↗</a> · <a href="https://doi.org/10.1371/journal.pcbi.1004409" target="_blank" rel="noreferrer">partis/ham study ↗</a></p>
        </section>
        <fieldset className="phylo-uca-placement-settings"><legend>Placement inference</legend>
          <label title="Gibbs/MH is the default and samples placement uncertainty with continuous—not discretized—branch lengths. Conditional ML is the fastest route. Grid mode explicitly integrates a displayed branch-length grid."><span>Inference route</span><select value={options.search.inferenceMode} onChange={(event) => updateSearch("inferenceMode", event.target.value as PhyloUcaOptions["search"]["inferenceMode"])}><option value="maximum-likelihood">Conditional ML · fastest</option><option value="grid-marginalization">Explicit grid marginalization</option><option value="gibbs-mh">Gibbs/MH · default · continuous placement</option></select></label>
          <label title="The default independently mixes retained V/J allele nucleotides at each fixed-alignment column. This is only a fast starting-position screen; retained points are always recomputed with the complete recombination HMM."><span>Starting-position screen</span><select value={options.search.screenMode ?? "vj-mixture"} onChange={(event) => updateSearch("screenMode", event.target.value as PhyloUcaOptions["search"]["screenMode"])}><option value="vj-mixture">V/J nucleotide mixture · default</option><option value="germline-guide">Single N-masked germline guide</option></select></label>
          <label title="Number of attachment positions screened along every observed-tree edge before full-HMM refinement. Endpoints and at least one branch-interior point are always included."><span>Screen points / edge</span><CommitNumberInput min="3" max="101" value={options.search.screenEdgeGridPoints ?? 5} onCommit={(value) => updateSearch("screenEdgeGridPoints", Math.max(3, Math.floor(value)))} /></label>
          <label title={options.search.inferenceMode === "gibbs-mh" ? "Number of leading V/J-screen candidates refined for initialization; zero refines every edge. Global MH proposals can still visit every edge." : "Top screen-ranked edges receiving the full recombination HMM. Set to zero to search every edge."}><span>{options.search.inferenceMode === "gibbs-mh" ? "Refined initializer edges" : "Full-HMM edges"} (0 = all)</span><CommitNumberInput min="0" max="10000" value={options.search.fullHmmEdges} onCommit={(value) => updateSearch("fullHmmEdges", Math.max(0, Math.floor(value)))} /></label>
          <label title="Maximum UCA-to-observed-tree branch length in expected substitutions per character."><span>Maximum UCA branch</span><CommitNumberInput min="0.001" max="3" step="0.01" value={options.search.maximumUcaBranchLength} onCommit={(value) => updateSearch("maximumUcaBranchLength", value)} /></label>
          {options.search.inferenceMode === "maximum-likelihood" && <>
            <label title="Alternating scalar optimizations of continuous pendant length and continuous within-edge attachment position under the full HMM."><span>ML coordinate rounds</span><CommitNumberInput min="1" max="10" value={options.search.mlOptimizationRounds} onCommit={(value) => updateSearch("mlOptimizationRounds", Math.max(1, Math.floor(value)))} /></label>
            <label title="Stopping tolerance in the transformed unit coordinate. Pendant length itself remains continuous and is concentrated near zero by a fourth-power coordinate transform."><span>ML search tolerance</span><CommitNumberInput min="0.00001" max="0.1" step="0.0005" value={options.search.mlOptimizationTolerance} onCommit={(value) => updateSearch("mlOptimizationTolerance", value)} /></label>
            <p className="phylo-uca-setting-note">The V/J approximation only selects starting edges. The reported point maximizes the complete HMM marginal likelihood; neither attachment nor branch length is marginalized.</p>
          </>}
          {options.search.inferenceMode === "grid-marginalization" && <>
            <label title="Attachment fractions evaluated on each retained edge. Endpoints and interior positions are included."><span>Attachment points / edge</span><CommitNumberInput min="3" max="101" value={options.search.edgeGridPoints} onCommit={(value) => updateSearch("edgeGridPoints", Math.max(3, Math.floor(value)))} /></label>
            <label title="Number of points in the explicit zero-plus-logarithmic pendant-length grid."><span>UCA branch grid points</span><CommitNumberInput min="2" max="101" value={options.search.branchGridPoints} onCommit={(value) => updateSearch("branchGridPoints", Math.max(2, Math.floor(value)))} /></label>
            <label title="Smallest positive pendant length in the explicit logarithmic grid; zero is always included separately."><span>Smallest positive branch</span><CommitNumberInput min="0.000000001" max={String(options.search.maximumUcaBranchLength)} step="0.00001" value={options.search.minimumPositiveUcaBranchLength} onCommit={(value) => updateSearch("minimumPositiveUcaBranchLength", value)} /></label>
            <label title="Number of highest-mass quadrature points receiving exact nucleotide, codon, and HMM-track posterior calculations."><span>Retained posterior points</span><CommitNumberInput min="1" max="10000" value={options.search.localPosteriorPoints} onCommit={(value) => updateSearch("localPosteriorPoints", Math.max(1, Math.floor(value)))} /></label>
            <div className="phylo-uca-grid-list"><strong>Exact pendant-length grid</strong><code>{displayedBranchGrid.map((value) => value.toPrecision(6)).join(", ")}</code></div>
          </>}
          {options.search.inferenceMode === "gibbs-mh" && <>
            <label title="Total exact HMM Gibbs updates, including burn-in."><span>MCMC iterations</span><CommitNumberInput min="2" max="100000" value={options.search.mcmcIterations} onCommit={(value) => updateSearch("mcmcIterations", Math.max(2, Math.floor(value)))} /></label>
            <label title="Initial Gibbs/MH iterations discarded from posterior output."><span>Burn-in iterations</span><CommitNumberInput min="0" max={String(Math.max(0, options.search.mcmcIterations - 1))} value={options.search.mcmcBurnIn} onCommit={(value) => updateSearch("mcmcBurnIn", Math.max(0, Math.floor(value)))} /></label>
            <label title="Retain one joint HMM/UCA draw after this many post-burn-in iterations."><span>Thinning interval</span><CommitNumberInput min="1" max="1000" value={options.search.mcmcThin} onCommit={(value) => updateSearch("mcmcThin", Math.max(1, Math.floor(value)))} /></label>
            <label title="Cheap continuous placement/length MH updates after each exact HMM Gibbs draw."><span>MH steps / Gibbs draw</span><CommitNumberInput min="1" max="100" value={options.search.mcmcMhStepsPerIteration} onCommit={(value) => updateSearch("mcmcMhStepsPerIteration", Math.max(1, Math.floor(value)))} /></label>
            <label title="Half-width of the reflected continuous pendant-length random-walk proposal, in substitutions/site."><span>Branch proposal scale</span><CommitNumberInput min="0.0000001" max="1" step="0.001" value={options.search.mcmcBranchProposalScale} onCommit={(value) => updateSearch("mcmcBranchProposalScale", value)} /></label>
            <label title="Half-width of the reflected continuous within-edge attachment-fraction proposal."><span>Position proposal scale</span><CommitNumberInput min="0.000001" max="1" step="0.01" value={options.search.mcmcPositionProposalScale} onCommit={(value) => updateSearch("mcmcPositionProposalScale", value)} /></label>
            <label title="Probability that an MH step proposes an attachment on another edge. Candidate probabilities use the V/J-only screen, with a nonzero floor for every edge."><span>Global-jump probability</span><CommitNumberInput min="0" max="1" step="0.01" value={options.search.mcmcGlobalJumpProbability} onCommit={(value) => updateSearch("mcmcGlobalJumpProbability", value)} /></label>
            <label title="Fraction of global edge jumps whose within-edge position is proposed near that edge's independently optimized V/J-screen position. The remaining fraction is uniform, preserving support everywhere."><span>Focused global-position mix</span><CommitNumberInput min="0" max="1" step="0.01" value={options.search.mcmcGlobalPositionMixture} onCommit={(value) => updateSearch("mcmcGlobalPositionMixture", value)} /></label>
            <label title="Circular half-width around the V/J-screen optimum used by the focused component of a global within-edge proposal. The exact proposal density is included in the Hastings ratio."><span>Focused global-position scale</span><CommitNumberInput min="0.001" max="0.5" step="0.01" value={options.search.mcmcGlobalPositionScale} onCommit={(value) => updateSearch("mcmcGlobalPositionScale", value)} /></label>
            <label title="Fraction of exact collapsed-refresh branch proposals drawn from the short-branch interval below. The remainder is uniform over the complete allowed branch range, so no continuous branch length is excluded."><span>Focused branch-proposal mix</span><CommitNumberInput min="0" max="1" step="0.01" value={options.search.mcmcGlobalBranchMixture} onCommit={(value) => updateSearch("mcmcGlobalBranchMixture", value)} /></label>
            <label title="Upper bound, in substitutions/site, of the short-branch component used only to propose collapsed refreshes. The exact mixture density is included in the Hastings ratio."><span>Focused branch-proposal max</span><CommitNumberInput min="0.000001" max={String(options.search.maximumUcaBranchLength)} step="0.005" value={options.search.mcmcGlobalBranchMaximum} onCommit={(value) => updateSearch("mcmcGlobalBranchMaximum", value)} /></label>
            <label title="Fraction of exact collapsed edge proposals weighted by the full-HMM initializer evaluations that were already computed. The remainder uses the broad V/J screen, leaving every screened edge reachable; the exact mixture density is included in the Hastings ratio."><span>Full-HMM edge-proposal mix</span><CommitNumberInput min="0" max="1" step="0.01" value={options.search.mcmcCollapsedInitializerMixture} onCommit={(value) => updateSearch("mcmcCollapsedInitializerMixture", value)} /></label>
            <label title="Run an exact collapsed independence update for edge, continuous within-edge position, and continuous UCA branch length after this many iterations. This breaks sticky UCA/placement coupling. Zero disables the refresh."><span>Collapsed refresh interval</span><CommitNumberInput min="0" max="10000" step="1" value={options.search.mcmcCollapsedRefreshInterval} onCommit={(value) => updateSearch("mcmcCollapsedRefreshInterval", Math.max(0, Math.floor(value)))} /></label>
            <label title="Fixed 32-bit seed for a reproducible chain."><span>MCMC seed</span><CommitNumberInput min="0" max="4294967295" value={options.search.mcmcSeed} onCommit={(value) => updateSearch("mcmcSeed", Math.max(0, Math.floor(value)))} /></label>
            <p className="phylo-uca-setting-note">Pendant length and within-edge position are continuous floating-point states. The sampler does not use the grid settings above; each proposal evaluates the exact GTR transition against cached directed tree messages.</p>
          </>}
          {options.search.inferenceMode !== "maximum-likelihood" && <>
            <label title="Mean of the exponential prior on UCA pendant length."><span>Branch-prior mean</span><CommitNumberInput min="0.0001" max="2" step="0.01" value={options.search.branchPriorMean} onCommit={(value) => updateSearch("branchPriorMean", value)} /></label>
            <label title="Uniform-length gives branches prior mass proportional to length; uniform-edge gives every observed-tree edge equal prior mass."><span>Attachment prior</span><select value={options.search.edgePrior} onChange={(event) => updateSearch("edgePrior", event.target.value as PhyloUcaOptions["search"]["edgePrior"])}><option value="uniform-length">Uniform over tree length</option><option value="uniform-edge">Uniform over edges</option></select></label>
          </>}
        </fieldset>
      </div>
      <div className="phylo-uca-default-reset"><div><strong>Restore calibrated starting values</strong><span>Resets every UCA character, candidate, HMM, and placement option; it does not alter the curated alignment.</span></div><button type="button" onClick={resetOptions}>Reset all UCA settings to defaults</button></div>
    </details>
    {running && <div className="phylo-uca-progress" role="status"><div><strong>{treeStage ? "Observed-only FastTree" : phaseLabel(progress)}</strong><span>{treeStage ? "Removing the germline guide and inferring the fixed tree" : progress?.detail}</span></div><progress max="1" value={progressFraction} /><b>{treeStage ? "tree" : progress?.total ? `${progress.processed.toLocaleString()} / ${progress.total.toLocaleString()}` : "working"}</b></div>}
    {error && <div className="inline-method-error" role="alert"><strong>Phylogenetic UCA inference stopped</strong><span>{error}</span><button type="button" onClick={() => setError("")}>Dismiss</button></div>}
    {result && <div ref={resultRef} className="phylo-uca-results">
      <div className="phylo-uca-result-head"><div><span className="section-kicker">{resultInferenceMode === "maximum-likelihood" ? "Conditional-ML placement + exact HMM posterior" : resultInferenceMode === "grid-marginalization" ? "Grid-marginalized placement + exact HMM posterior" : "Continuous Gibbs/MH placement + sampled HMM posterior"}</span><h4>{result.mapVCall || "V?"} · {result.mapDCalls.length ? result.mapDCalls.join(" → ") : "no D"} · {result.mapJCall || "J?"}</h4><p>{result.characterModel === "nucleotide-gtr4" ? "Four-state nucleotide GTR" : "Five-state A/C/G/T/gap GTR"} · edge {result.bestPlacement.endpointA} ↔ {result.bestPlacement.endpointB} at {(result.bestPlacement.edgeFraction * 100).toFixed(1)}% · UCA branch {result.bestPlacement.ucaBranchLength.toFixed(6)} · {resultInferenceMode === "maximum-likelihood" ? "no placement marginalization" : resultInferenceMode === "gibbs-mh" ? `${result.mcmcDiagnostics?.retainedSamples ?? result.placements.length} retained MCMC draws` : `effective grid points ${result.effectivePlacementCount.toFixed(2)}`}</p></div><div className="result-actions"><button type="button" onClick={() => download(`>phylo_UCA_joint_MAP_aligned\n${result.mapAlignedSequence}\n>phylo_UCA_marginal_consensus_aligned\n${result.posteriorConsensusAligned}\n`, `${base}.fasta`)}>UCA FASTA ↓</button><button type="button" onClick={() => download(posteriorTsv(result), `${base}.nucleotide-posterior.tsv`)}>Nucleotide TSV ↓</button>{Boolean(result.codonPosterior?.length) && <button type="button" onClick={() => download(codonPosteriorTsv(result), `${base}.codon-posterior.tsv`)}>Codon TSV ↓</button>}<button type="button" onClick={() => download(result.placedTreeNewick, `${base}.placed-tree.nwk`)}>Placed tree ↓</button><button type="button" onClick={() => download(JSON.stringify(result, null, 2), `${base}.json`, "application/json")}>Complete JSON ↓</button></div></div>
      <div className="phylo-uca-stats"><article><span>Inference route</span><strong>{resultInferenceMode === "maximum-likelihood" ? "Conditional ML" : resultInferenceMode === "grid-marginalization" ? "Grid marginalization" : "Continuous Gibbs/MH"}</strong></article><article><span>Candidate set</span><strong>{result.candidateReport.v.length} V · {result.candidateReport.d.length} D · {result.candidateReport.j.length} J</strong></article><article><span>Placement log marginal</span><strong>{result.logMarginalLikelihood.toFixed(2)}</strong></article><article><span>Runtime</span><strong>{(result.elapsedMs / 1000).toFixed(2)} s</strong></article></div>
      {result.mcmcDiagnostics && <PhyloUcaMcmcMixing diagnostics={result.mcmcDiagnostics} />}
      {result.dCountPosterior?.length ? <section className="phylo-uca-d-count"><header><div><span className="section-kicker">Retained joint draws</span><h5>Posterior number of D segments</h5></div><p>No additional HMM evaluations: these counts come directly from the retained Gibbs paths.</p></header><div>{result.dCountPosterior.map((point) => <article key={point.dCount}><span><b>{point.dCount}</b> D{point.dCount === 1 ? "" : "s"}</span><div><i style={{ width: `${100 * point.probability}%` }} /></div><strong>{(100 * point.probability).toFixed(1)}%</strong><small>{point.samples.toLocaleString()} draws</small></article>)}</div></section> : null}
      <PhyloUcaPlacementMap newick={result.observedTreeNewick} placements={result.placements} inferenceMode={resultInferenceMode} title={`${lineageLabel} UCA placement likelihood surface`} />
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
