import { DEFAULT_ALLELE_REFINEMENT_OPTIONS, type AlleleRefinementOptions } from "./allele-refinement/types.ts";
import { DEFAULT_MISSING_ALLELE_OPTIONS, type MissingAlleleOptions } from "./germline-evidence.ts";
import type { DenoiseOptions } from "./post-analysis-core.ts";
import type { CompiledReferences, ScopeKey } from "./reference-pack.ts";
import { DEFAULT_REPERTOIRE_SELECTION, type RepertoireSelectionOptions } from "./repertoire-selection.ts";
import { DEFAULT_PIPELINE_PLAN, type DatasetManifestEntry, type PipelinePlan, type StudyDesign } from "./study-design.ts";
import type { AssignerStrategy, CallingProfile, DoubleDScreenOptions } from "./swiftig-runtime.ts";

export const SWIG_CLI_CONFIG_SCHEMA = 1 as const;

export type CliInputFormat = "auto" | "fasta" | "fastq" | "airr";

export interface CliInputDataset {
  path: string;
  format?: CliInputFormat;
  datasetId?: string;
  sampleId?: string;
  /** Samples sharing this value are treated as belonging to the same donor. */
  subjectId?: string;
  cohort?: string;
  timepoint?: string;
  compartment?: string;
}

export interface CliReferenceConfig {
  species: string;
  scope: ScopeKey;
  /** Exact references exported by Swig Web. These take precedence over built-ins. */
  inline?: Partial<Record<"V" | "D" | "J" | "C", string>>;
  /** Optional local FASTA replacements, resolved relative to the config file. */
  files?: Partial<Record<"V" | "D" | "J" | "C", string>>;
}

export interface CliAnnotationConfig {
  workers: number;
  batchRecords: number;
  callingProfile: CallingProfile;
  assignerStrategy: AssignerStrategy;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  doubleD: DoubleDScreenOptions;
}

export interface CliPipelineConfig {
  collapse: PipelinePlan["collapse"] & { denoise: DenoiseOptions };
  chimera: PipelinePlan["chimera"] & {
    priorProbability: number;
    baseMutationProbability: number;
    mutationRates: number[];
    mutationSwitchProbability: number;
    minimumDfr: number;
  };
  selection: RepertoireSelectionOptions & { enabled: boolean };
  alleleRefinement: AlleleRefinementOptions & {
    enabled: boolean;
    reassignmentPolicy: "best" | "confidence";
    applyMinimumPosterior: number;
  };
  lineage: PipelinePlan["lineage"] & { maxCandidateComparisons: number };
  shm: PipelinePlan["shm"];
  missingAlleles: MissingAlleleOptions & { enabled: boolean };
}

export interface SwigCliConfig {
  schema: typeof SWIG_CLI_CONFIG_SCHEMA;
  application: "swig-cli";
  studyName: string;
  studyDesign: StudyDesign;
  inputs: CliInputDataset[];
  references: CliReferenceConfig;
  annotation: CliAnnotationConfig;
  pipeline: CliPipelineConfig;
  output: {
    directory: string;
    prefix: string;
    writeAnnotatedAirr: boolean;
    writeLineageStudy: boolean;
  };
}

export type PartialSwigCliConfig = Partial<Omit<SwigCliConfig, "inputs" | "references" | "annotation" | "pipeline" | "output">> & {
  inputs?: CliInputDataset[];
  references?: Partial<CliReferenceConfig>;
  annotation?: Partial<CliAnnotationConfig> & { doubleD?: Partial<DoubleDScreenOptions> };
  pipeline?: Partial<{
    collapse: Partial<Omit<CliPipelineConfig["collapse"],"denoise">> & { denoise?: Partial<DenoiseOptions> };
    chimera: Partial<CliPipelineConfig["chimera"]>;
    selection: Partial<CliPipelineConfig["selection"]>;
    alleleRefinement: Partial<CliPipelineConfig["alleleRefinement"]>;
    lineage: Partial<CliPipelineConfig["lineage"]>;
    shm: Partial<CliPipelineConfig["shm"]>;
    missingAlleles: Partial<CliPipelineConfig["missingAlleles"]>;
  }>;
  output?: Partial<SwigCliConfig["output"]>;
};

const DEFAULT_DENOISE: DenoiseOptions = {
  mode: "conservative",
  errorRate: 0.00473,
  alpha: 0.01,
  callResolution: "allele",
  ambiguity: "strict",
  minimumParentCount: 2,
  ambiguousPolicy: "exclude",
  unresolvedPolicy: "discard",
  fadNeighborThreshold: 1,
  fadMethod: 2,
  expectedZeroErrorFraction: 1,
  maximumHammingDistance: 1,
  maximumEditDistance: 2,
  minimumIndelParentRatio: 2,
  maxCandidatesPerVariant: 50_000,
  scope: "sample",
  respectConstantCall: true,
};

export const DEFAULT_CLI_CONFIG: SwigCliConfig = {
  schema: SWIG_CLI_CONFIG_SCHEMA,
  application: "swig-cli",
  studyName: "swig-study",
  studyDesign: "independent",
  inputs: [],
  references: { species: "Homo sapiens", scope: "BCR" },
  annotation: {
    workers: 0,
    batchRecords: 2_000,
    callingProfile: "truth_optimized",
    assignerStrategy: "aer",
    minimumIdentity: 0.6,
    strand: 0,
    doubleD: { mode: "off", minimumVjSpan: 40, seedLength: 11, pseudoTrim: 5, maximumPseudoMismatches: 3, minimumScoreGain: 8 },
  },
  pipeline: {
    collapse: { ...DEFAULT_PIPELINE_PLAN.collapse, enabled: true, denoise: { ...DEFAULT_DENOISE } },
    chimera: {
      ...DEFAULT_PIPELINE_PLAN.chimera,
      enabled: false,
      priorProbability: 0.05,
      baseMutationProbability: 0.005,
      mutationRates: [0,0.0179,0.0357,0.0536,0.0714,0.0893,0.1071,0.125,0.1429,0.1607,0.1786,0.1964,0.2143,0.2321,0.25],
      mutationSwitchProbability: 0.01,
      minimumDfr: 1,
    },
    selection: { ...DEFAULT_REPERTOIRE_SELECTION, enabled: false },
    alleleRefinement: { ...DEFAULT_ALLELE_REFINEMENT_OPTIONS, segments: [...DEFAULT_ALLELE_REFINEMENT_OPTIONS.segments], enabled: false, reassignmentPolicy: "confidence", applyMinimumPosterior: 0.8 },
    lineage: { ...DEFAULT_PIPELINE_PLAN.lineage, enabled: true, maxCandidateComparisons: 50_000 },
    shm: { ...DEFAULT_PIPELINE_PLAN.shm, enabled: true },
    missingAlleles: { ...DEFAULT_MISSING_ALLELE_OPTIONS, enabled: false },
  },
  output: { directory: "swig-output", prefix: "swig", writeAnnotatedAirr: true, writeLineageStudy: true },
};

function finite(value: unknown, fallback: number): number {
  const parsed=Number(value);
  return Number.isFinite(parsed)?parsed:fallback;
}

function safeIdentifier(value: string, fallback: string): string {
  return value.trim().replace(/[^A-Za-z0-9_.-]+/g,"-").replace(/^-+|-+$/g,"")||fallback;
}

export function normalizeCliConfig(value: PartialSwigCliConfig): SwigCliConfig {
  const collapse=value.pipeline?.collapse??{};
  const chimera=value.pipeline?.chimera??{};
  const selection=value.pipeline?.selection??{};
  const allele=value.pipeline?.alleleRefinement??{};
  const lineage=value.pipeline?.lineage??{};
  const shm=value.pipeline?.shm??{};
  const missing=value.pipeline?.missingAlleles??{};
  const inputs=(value.inputs??[]).map((input,index)=>{
    const stem=input.path.split(/[\\/]/).pop()?.replace(/\.(?:fasta?|fna|fas|fastq|fq|airr\.tsv|tsv)(?:\.gz)?$/i,"")||`dataset-${index+1}`;
    const datasetId=safeIdentifier(input.datasetId??stem,`dataset-${index+1}`);
    const sampleId=safeIdentifier(input.sampleId??datasetId,`sample-${index+1}`);
    return {...input,format:input.format??"auto",datasetId,sampleId,subjectId:safeIdentifier(input.subjectId??sampleId,`subject-${index+1}`),cohort:input.cohort??"",timepoint:input.timepoint??"",compartment:input.compartment??""};
  });
  const annotation={...DEFAULT_CLI_CONFIG.annotation,...value.annotation,doubleD:{...DEFAULT_CLI_CONFIG.annotation.doubleD,...value.annotation?.doubleD}};
  annotation.workers=Math.max(0,Math.floor(finite(annotation.workers,0)));
  annotation.batchRecords=Math.max(1,Math.floor(finite(annotation.batchRecords,2_000)));
  annotation.minimumIdentity=Math.max(0,Math.min(1,finite(annotation.minimumIdentity,0.6)));
  return {
    ...DEFAULT_CLI_CONFIG,
    ...value,
    schema:SWIG_CLI_CONFIG_SCHEMA,
    application:"swig-cli",
    inputs,
    references:{...DEFAULT_CLI_CONFIG.references,...value.references,inline:value.references?.inline?{...value.references.inline}:undefined,files:value.references?.files?{...value.references.files}:undefined},
    annotation,
    pipeline:{
      collapse:{...DEFAULT_CLI_CONFIG.pipeline.collapse,...collapse,denoise:{...DEFAULT_DENOISE,...collapse.denoise,mode:(collapse.mode??collapse.denoise?.mode??DEFAULT_DENOISE.mode) as DenoiseOptions["mode"],scope:collapse.scope??collapse.denoise?.scope??DEFAULT_DENOISE.scope,respectConstantCall:collapse.respectConstantCall??collapse.denoise?.respectConstantCall??true}},
      chimera:{...DEFAULT_CLI_CONFIG.pipeline.chimera,...chimera,mutationRates:[...(chimera.mutationRates??DEFAULT_CLI_CONFIG.pipeline.chimera.mutationRates)]},
      selection:{...DEFAULT_CLI_CONFIG.pipeline.selection,...selection},
      alleleRefinement:{...DEFAULT_CLI_CONFIG.pipeline.alleleRefinement,...allele,segments:[...(allele.segments??DEFAULT_CLI_CONFIG.pipeline.alleleRefinement.segments)]},
      lineage:{...DEFAULT_CLI_CONFIG.pipeline.lineage,...lineage},
      shm:{...DEFAULT_CLI_CONFIG.pipeline.shm,...shm},
      missingAlleles:{...DEFAULT_CLI_CONFIG.pipeline.missingAlleles,...missing},
    },
    output:{...DEFAULT_CLI_CONFIG.output,...value.output},
  };
}

export interface BrowserCliExport {
  studyName: string;
  studyDesign: StudyDesign;
  datasets: DatasetManifestEntry[];
  references: CompiledReferences;
  species: string;
  scope: ScopeKey;
  workers: number;
  callingProfile: CallingProfile;
  assignerStrategy: AssignerStrategy;
  minimumIdentity: number;
  strand: 0|1|2;
  doubleD: DoubleDScreenOptions;
  pipeline: PipelinePlan;
  collapseOptions?: Partial<DenoiseOptions>;
  selectionOptions?: Partial<RepertoireSelectionOptions>;
  alleleOptions?: Partial<AlleleRefinementOptions> & { reassignmentPolicy?: "best"|"confidence"; applyMinimumPosterior?: number };
  missingAlleleOptions?: Partial<MissingAlleleOptions>;
}

/** Build the portable, human-editable configuration emitted by Swig Web. */
export function cliConfigFromBrowser(run: BrowserCliExport): SwigCliConfig {
  const config=normalizeCliConfig({
    studyName:run.studyName,
    studyDesign:run.studyDesign,
    inputs:run.datasets.map((dataset)=>({path:dataset.inputName,format:"auto",datasetId:dataset.datasetId,sampleId:dataset.sampleId,subjectId:dataset.subjectId,cohort:dataset.cohort,timepoint:dataset.timepoint,compartment:dataset.compartment??""})),
    references:{species:run.species,scope:run.scope,inline:{V:run.references.V,D:run.references.D,J:run.references.J,C:run.references.C}},
    annotation:{workers:run.workers,callingProfile:run.callingProfile,assignerStrategy:run.assignerStrategy,minimumIdentity:run.minimumIdentity,strand:run.strand,doubleD:run.doubleD},
    pipeline:{
      collapse:{...run.pipeline.collapse,denoise:run.collapseOptions},
      chimera:{...run.pipeline.chimera},
      selection:{...run.pipeline.selection,...run.selectionOptions},
      alleleRefinement:{...run.pipeline.alleleRefinement,...run.alleleOptions},
      lineage:{...run.pipeline.lineage},
      shm:{...run.pipeline.shm},
      missingAlleles:{...run.pipeline.missingAlleles,...run.missingAlleleOptions},
    },
    output:{prefix:safeIdentifier(run.studyName,"swig")},
  });
  return config;
}
