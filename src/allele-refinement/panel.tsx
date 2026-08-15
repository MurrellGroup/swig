import { CommitNumberInput } from "../commit-number-input.tsx";
import { DATASET_SCOPE_LABELS, type DatasetScope } from "../study-design.ts";
import type { AlleleRefinementOptions, AlleleRefinementResult, RefinementSegment } from "./types.ts";
import { adaptiveNeighbourOdds } from "./evidence.ts";

interface Props {
  options: AlleleRefinementOptions;
  onOptionsChange: (options: AlleleRefinementOptions) => void;
  result: AlleleRefinementResult | null;
  applied: boolean;
  applyMinimumPosterior: number;
  onApplyMinimumPosteriorChange: (value: number) => void;
  busy: boolean;
  progress: { processed: number; total: number; phase: string } | null;
  onRun: () => void;
  onApply: () => void;
  onReset: () => void;
  onDownloadModel: () => void;
  onDownloadSidecar: () => void;
  onDownloadAirr: () => void;
}

const SEGMENT_LABELS: Record<RefinementSegment, string> = {
  V: "V alleles",
  D: "D alleles · experimental",
  J: "J alleles",
};

export function AlleleRefinementPanel({
  options, onOptionsChange, result, applied, applyMinimumPosterior,
  onApplyMinimumPosteriorChange, busy, progress, onRun, onApply, onReset,
  onDownloadModel, onDownloadSidecar, onDownloadAirr,
}: Props) {
  const update = <K extends keyof AlleleRefinementOptions>(key: K, value: AlleleRefinementOptions[K]) => onOptionsChange({ ...options, [key]: value });
  const toggleSegment = (segment: RefinementSegment) => {
    const selected = options.segments.includes(segment)
      ? options.segments.filter((value) => value !== segment)
      : [...options.segments, segment];
    update("segments", selected);
  };
  const segmentResults = Object.values(result?.segments ?? {}).filter(Boolean);
  const modeledRows = segmentResults.reduce((sum, segment) => sum + segment!.modeledRows, 0);
  const changedRows = segmentResults.reduce((sum, segment) => sum + segment!.changedMapRows, 0);
  const nonZeros = segmentResults.reduce((sum, segment) => sum + segment!.matrixNonZeros, 0);
  const models = segmentResults.flatMap((segment) => segment!.models);
  const topAlleles = models.flatMap((model) => model.alleles.slice(0, 8).map((allele) => ({ model, allele })))
    .sort((left, right) => left.model.scopeValue.localeCompare(right.model.scopeValue)
      || left.model.locus.localeCompare(right.model.locus)
      || left.model.segment.localeCompare(right.model.segment)
      || right.allele.posteriorMean - left.allele.posteriorMean)
    .slice(0, 250);
  const progressFraction = progress ? Math.max(0, Math.min(1, progress.processed / Math.max(1, progress.total))) : 0;
  const leakageExamples = [0, 0.05, 0.1].map((shm) => {
    const odds = adaptiveNeighbourOdds(shm, options);
    return `${Math.round(shm * 100)}% SHM → ${(100 * odds / (1 + odds)).toFixed(2)}%`;
  }).join(" · ");

  return <>
    <div className="allele-refinement-config">
      <div className="control-grid three">
        <label title="Reference usage is learned independently inside this boundary. Donor/subject is the conservative default for longitudinal or compartmental studies."><span>Pooling boundary</span><select value={options.scope} onChange={(event) => update("scope", event.target.value as DatasetScope)}>{Object.entries(DATASET_SCOPE_LABELS).map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select></label>
        <fieldset className="allele-segment-picker"><legend title="V and J are enabled by default. D inference is less identifiable because D segments are short and heavily trimmed.">Reference segments</legend>{(["V", "J", "D"] as RefinementSegment[]).map((segment) => <label key={segment}><input type="checkbox" checked={options.segments.includes(segment)} onChange={() => toggleSegment(segment)} /><span>{SEGMENT_LABELS[segment]}</span></label>)}</fieldset>
        <label title="Unique weighting prevents expanded clones or PCR abundance from acting as independent genotype evidence. Abundance weighting is available for deliberate usage estimation."><span>Record weighting</span><select value={options.weighting} onChange={(event) => update("weighting", event.target.value as AlleleRefinementOptions["weighting"])}><option value="unique">One vote per active record · default</option><option value="abundance">Weight by duplicate_count</option></select></label>
      </div>
      <div className="allele-refinement-action"><div><strong>{result ? `${result.activeRecords.toLocaleString()} active records modeled` : "No repertoire posterior fitted"}</strong><small>Explicit co-optimal calls enter with equal local weight. Only reference-graph neighbours inside the selected edit radius receive leakage support.</small></div><button className="post-primary" type="button" disabled={busy || !options.segments.length} onClick={onRun}>{result ? "Refit repertoire allele model" : "Fit repertoire allele model"}</button></div>
      <details className="post-advanced"><summary>Advanced evidence-kernel and variational settings</summary><div><div className="control-grid four allele-refinement-advanced">
        <label title="Irreducible relative evidence odds for a one-nucleotide reference neighbour when estimated SHM is zero. This is an assignment-model leakage term, not a sequencing-error estimate."><span>Zero-SHM neighbour odds</span><CommitNumberInput min="0" max="0.5" step="0.001" value={options.baselineNeighbourOdds} onCommit={(value) => update("baselineNeighbourOdds", value)} /></label>
        <label title="Multiplier on the mechanistic mu/[3(1-mu)] contribution for a specific nucleotide substitution at a diagnostic allele position. Set to zero for SHM-independent leakage."><span>SHM sensitivity</span><CommitNumberInput min="0" max="10" step="0.1" value={options.shmLeakageSensitivity} onCommit={(value) => update("shmLeakageSensitivity", value)} /></label>
        <label title="Upper bound on one-SNP neighbour evidence odds for highly mutated or poorly matched reads."><span>Maximum neighbour odds</span><CommitNumberInput min="0.001" max="0.99" step="0.01" value={options.maximumNeighbourOdds} onCommit={(value) => update("maximumNeighbourOdds", value)} /></label>
        <label title="The best-reference V mutation fraction is clamped at this value before the SHM leakage term is calculated."><span>SHM estimate cap</span><CommitNumberInput min="0" max="0.9" step="0.01" value={options.maximumShm} onCommit={(value) => update("maximumShm", value)} /></label>
        <label title="Maximum bounded nucleotide edit distance from any reported or retained alternative reference. Candidate rows remain sparse."><span>Neighbour edit radius</span><CommitNumberInput min="0" max="5" step="1" value={options.neighbourRadius} onCommit={(value) => update("neighbourRadius", Math.floor(value))} /></label>
        <label title="Symmetric Dirichlet pseudo-count for every locus-matched reference node. Nodes absent from all sparse candidate rows retain only this prior mass and are not materialized per record."><span>Dirichlet alpha / allele</span><CommitNumberInput min="0.000001" max="100" step="0.01" value={options.alphaPerAllele} onCommit={(value) => update("alphaPerAllele", value)} /></label>
        <label title="Scale converting SwiftIG score differences for retained alternative hits into relative local evidence. Smaller values penalize lower scores more strongly."><span>Alternative score temperature</span><CommitNumberInput min="0.001" max="100" step="0.1" value={options.alternativeScoreTemperature} onCommit={(value) => update("alternativeScoreTemperature", value)} /></label>
        <label title="Fallback local weight for a retained alternative that lacks a numeric alignment score."><span>Unscored alternative weight</span><CommitNumberInput min="0" max="1" step="0.01" value={options.unscoredAlternativeWeight} onCommit={(value) => update("unscoredAlternativeWeight", value)} /></label>
        <label title="Hard memory guard for unusually dense reference neighbourhoods. Candidates are retained in descending local evidence order."><span>Candidate cap / record</span><CommitNumberInput min="1" max="10000" step="1" value={options.maxCandidatesPerRow} onCommit={(value) => update("maxCandidatesPerRow", Math.floor(value))} /></label>
        <label title="Maximum coordinate-ascent updates for each independent pooling group."><span>Maximum iterations</span><CommitNumberInput min="1" max="10000" step="10" value={options.maxIterations} onCommit={(value) => update("maxIterations", Math.floor(value))} /></label>
        <label title="Stop when the maximum relative change in any Dirichlet parameter falls below this value."><span>Convergence tolerance</span><CommitNumberInput min="0.000000001" max="0.1" step="0.000001" value={options.convergenceTolerance} onCommit={(value) => update("convergenceTolerance", value)} /></label>
      </div><div className="algorithm-note"><strong>Implied one-SNP neighbour probability</strong><span>{leakageExamples}. Displayed probabilities are odds/(1+odds); distance-d evidence uses the per-edit odds to the dth power.</span></div></div></details>
    </div>
    {progress && <div className="allele-refinement-progress" role="status"><div><strong>{progress.phase}</strong><span>{progress.processed.toLocaleString()} / {progress.total.toLocaleString()}</span></div><progress max="1" value={progressFraction} /></div>}
    {result && <div className="allele-refinement-results">
      <div className="post-stat-grid"><article><span>Independent models</span><strong>{models.length.toLocaleString()}</strong></article><article><span>Modeled segment-rows</span><strong>{modeledRows.toLocaleString()}</strong></article><article><span>Repertoire MAP differs</span><strong>{changedRows.toLocaleString()}</strong></article><article><span>Sparse nonzeros</span><strong>{nonZeros.toLocaleString()}</strong></article></div>
      <div className="allele-apply-row"><div><strong>{applied ? "Posterior MAP calls are active for downstream lineage analysis" : "Original AIRR calls remain active downstream"}</strong><small>Only MAP calls at or above the selected posterior threshold replace enabled segment calls. Original calls remain in the AIRR source and all refined exports.</small></div><label title="Calls below this posterior probability remain unchanged when refinement is applied."><span>Apply at posterior ≥</span><CommitNumberInput min="0" max="1" step="0.01" value={applyMinimumPosterior} onCommit={onApplyMinimumPosteriorChange} /></label><button type="button" className="post-primary" disabled={busy} onClick={onApply}>Apply to downstream calls</button>{applied && <button type="button" disabled={busy} onClick={onReset}>Restore original calls</button>}</div>
      {result.warnings.map((warning) => <p className="scientific-note warning" key={warning}><span>!</span>{warning}</p>)}
      <div className="result-actions"><button type="button" onClick={onDownloadModel}>Model summary</button><button type="button" onClick={onDownloadSidecar}>Full sparse posterior sidecar</button><button type="button" onClick={onDownloadAirr}>AIRR with thresholded refined calls</button></div>
      <div className="allele-model-table"><table><thead><tr><th>Pool</th><th>Locus</th><th>Segment</th><th>Reference allele / identical class</th><th>Posterior use</th><th>Expected assignments</th><th>Local evidence</th></tr></thead><tbody>{topAlleles.map(({ model, allele }) => <tr key={`${model.key}-${allele.nodeIndex}`}><td>{model.scopeValue}</td><td>{model.locus}</td><td>{model.segment}</td><td><code>{allele.names.join(", ")}</code></td><td>{(allele.posteriorMean * 100).toFixed(3)}%</td><td>{allele.expectedAssignments.toFixed(2)}</td><td>{allele.localEvidenceAssignments.toFixed(2)}</td></tr>)}</tbody></table></div>
      <p className="scientific-note"><span>i</span>Posterior-use denominators include every locus-matched database node. Prior-only nodes absent from all sparse read neighbourhoods remain implicit; their counts are included in the model-summary export.</p>
      {topAlleles.length >= 250 && <p className="scientific-note"><span>i</span>The interactive table is limited to 250 entries. The model-summary export contains every modeled reference node.</p>}
    </div>}
  </>;
}
