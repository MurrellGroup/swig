export interface MissingAlleleOptions {
  /** Retained for version-one session compatibility. Screening is always lineage based. */
  unit: "lineage" | "record";
  minimumIndependentUnits: number;
  minimumCoveredUnits: number;
  minimumAlleleFraction: number;
  maximumShmRate: number;
  minimumAlignedBases: number;
  maximumCandidateSnps: number;
  maximumPValue: number;
  minimumDistinctJCalls: number;
  minimumDistinctJunctionLengths: number;
  minimumDistinctCdr3s: number;
  minimumLinkedFraction: number;
  minimumNearGermlineUnits: number;
  maximumOtherAlternateFraction: number;
}

export const DEFAULT_MISSING_ALLELE_OPTIONS: MissingAlleleOptions = {
  unit: "lineage",
  minimumIndependentUnits: 6,
  minimumCoveredUnits: 20,
  minimumAlleleFraction: 0.2,
  maximumShmRate: 0.08,
  minimumAlignedBases: 180,
  maximumCandidateSnps: 6,
  maximumPValue: 1e-6,
  minimumDistinctJCalls: 3,
  minimumDistinctJunctionLengths: 3,
  minimumDistinctCdr3s: 6,
  minimumLinkedFraction: 0.9,
  minimumNearGermlineUnits: 3,
  maximumOtherAlternateFraction: 0.02,
};

export interface CandidateSubstitution {
  position: number;
  reference: string;
  alternate: string;
  support: number;
  coverage: number;
  fraction: number;
  pValue: number;
  hotspot: boolean;
}

export interface MissingAlleleCandidate {
  id: string;
  vCall: string;
  parentAllele: string;
  substitutions: CandidateSubstitution[];
  independentUnits: number;
  coveredUnits: number;
  alleleFraction: number;
  pValue: number;
  distinctJCalls: number;
  distinctJunctionLengths: number;
  distinctCdr3s: number;
  distinctSubjects: number;
  nearGermlineUnits: number;
  referencePresentUnits: number;
  conflictingUnits: number;
  otherAlternateFraction: number;
  meanBackgroundShm: number;
  sequence: string;
  supportingUnits: string[];
  supportingUnitsTruncated: boolean;
  caution: string;
}

export interface MissingAlleleDashboard {
  mode: "lineage";
  validationPasses: 2;
  inputRecords: number;
  eligibleRecords: number;
  independentUnits: number;
  vAllelesTested: number;
  candidatePatternsTested: number;
  referenceVetoedLineagePatterns: number;
  conflictingLineagePatterns: number;
  proposalTruncations: number;
  candidates: MissingAlleleCandidate[];
  warnings: string[];
}

interface Mutation {
  position: number;
  reference: string;
  alternate: string;
}

interface AlignedBase {
  reference: string;
  query: string;
}

interface Observation {
  ordinal: number;
  unit: number;
  subject: string;
  vCall: string;
  jCall: string;
  cdr3HashA: number;
  cdr3HashB: number;
  junctionLength: number;
  aligned: number;
  rate: number;
  mutations: Mutation[];
  coverageStart: number;
  coverageEnd: number;
  /** Flat start/end pairs exist only when internal query gaps split coverage. */
  coverageBreaks?: number[];
  baseByPosition?: Map<number, AlignedBase>;
  validationIndex?: number;
}

interface EventEvidence {
  position: number;
  reference: string;
  alternate: string;
  support: number;
  coverage: number;
  fraction: number;
  pValue: number;
  hotspot: boolean;
}

interface DiscoveryEvent extends EventEvidence {
  supportingUnits: Set<number>;
}

interface CandidateProposal {
  vCall: string;
  sequence: string;
  events: EventEvidence[];
  preliminarySupport: Set<number>;
  meanBackgroundShm: number;
  nullProbability: number;
  pValue: number;
}

interface ProposalGroup {
  proposals: CandidateProposal[];
  observations: Observation[];
  flags: Uint8Array;
  mixedVCall: Uint8Array;
}

const FLAG_JOINT_COVERAGE = 1;
const FLAG_JOINT_ALTERNATE = 2;
const FLAG_REFERENCE_PRESENT = 4;
const FLAG_CONFLICT_PRESENT = 8;
const MAX_PROPOSALS_PER_V = 64;
const SUPPORTING_UNIT_EXPORT_LIMIT = 1_000;
const EMPTY_MUTATIONS:Mutation[]=[];

function topCall(value: string): string { return value.split(",")[0]?.trim() ?? ""; }
function cleanAlignment(value: string): string { return value.toUpperCase().replace(/\./g, "-").replace(/[^ACGTN-]/g, "N"); }
function cleanCdr3(value: string): string { return value.toUpperCase().replace(/[^ACGTN]/g, ""); }
function cleanReference(value: string): string { return value.toUpperCase().replace(/[^ACGTN]/g, ""); }
function lineageUnit(lineageId: number): number { return Number.isFinite(lineageId) && lineageId > 0 ? Math.floor(lineageId) : 0; }
function displayUnit(item:Observation): string { return `${item.subject} · lineage:${item.unit}`; }
function cdr3Fingerprint(value:string):[number,number]{let first=2166136261,second=2246822519;for(let index=0;index<value.length;index+=1){const code=value.charCodeAt(index);first=Math.imul(first^code,16777619);second=Math.imul(second^code,3266489917);}return[first>>>0,second>>>0];}

function parseFasta(value: string): Map<string, string> {
  const result = new Map<string, string>(); let name = "", sequence = "";
  const commit = () => { if (name && sequence) result.set(name, cleanReference(sequence)); };
  for (const line of value.split(/\r?\n/)) {
    if (line.startsWith(">")) { commit(); name = line.slice(1).trim().split(/\s+/)[0]; sequence = ""; }
    else sequence += line.trim();
  }
  commit(); return result;
}

function logGamma(value: number): number {
  const coefficients = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, .001208650973866179, -.000005395239384953];
  let x = value, y = value, tmp = x + 5.5; tmp -= (x + .5) * Math.log(tmp); let series = 1.000000000190015;
  for (const coefficient of coefficients) { y += 1; series += coefficient / y; }
  return -tmp + Math.log(2.5066282746310005 * series / x);
}

function betaFraction(a: number, b: number, x: number): number {
  const maximum = 200, epsilon = 3e-12, tiny = 1e-300; const qab = a + b, qap = a + 1, qam = a - 1; let c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < tiny) d = tiny; d = 1 / d; let h = d;
  for (let m = 1; m <= maximum; m += 1) {
    const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny; c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny; d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny; c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny; d = 1 / d; const delta = d * c; h *= delta; if (Math.abs(delta - 1) < epsilon) break;
  }
  return h;
}

function regularizedBeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0; if (x >= 1) return 1;
  const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? bt * betaFraction(a, b, x) / a : 1 - bt * betaFraction(b, a, 1 - x) / b;
}

export function binomialSurvival(successes: number, trials: number, probability: number): number {
  if (successes <= 0) return 1; if (successes > trials) return 0;
  const p = Math.min(1 - 1e-12, Math.max(1e-12, probability));
  return Math.min(1, Math.max(0, regularizedBeta(p, successes, trials - successes + 1)));
}

/** AID WRCY/RGYW context at the substituted base. */
function hotspotAt(sequence: string, zeroBasedPosition: number): boolean {
  const base = sequence[zeroBasedPosition]; const isW = (value = "") => value === "A" || value === "T"; const isR = (value = "") => value === "A" || value === "G"; const isY = (value = "") => value === "C" || value === "T";
  if (base === "C") return isW(sequence[zeroBasedPosition - 2]) && isR(sequence[zeroBasedPosition - 1]) && isY(sequence[zeroBasedPosition + 1]);
  if (base === "G") return isR(sequence[zeroBasedPosition - 1]) && isY(sequence[zeroBasedPosition + 1]) && isW(sequence[zeroBasedPosition + 2]);
  return false;
}

function parseObservation(row: Record<string, string>, ordinal: number, unit: number, retainBases = false): Observation | null {
  const vCall = topCall(row.v_call); const query = cleanAlignment(row.v_sequence_alignment || ""); const germline = cleanAlignment(row.v_germline_alignment || "");
  if (!unit || !vCall || !query || !germline) return null;
  const length = Math.min(query.length, germline.length); let mutations: Mutation[] = EMPTY_MUTATIONS; const covered:number[]=[]; const baseByPosition = retainBases ? new Map<number, AlignedBase>() : undefined;
  const reportedStart = Math.floor(Number(row.v_germline_start)); let position = Number.isFinite(reportedStart) && reportedStart > 0 ? reportedStart - 1 : 0; let aligned = 0, intervalStart = 0, lastCovered = 0;
  for (let column = 0; column < length; column += 1) {
    const q = query[column], g = germline[column];
    if (g !== "-") position += 1;
    if (/^[ACGT]$/.test(q) && /^[ACGT]$/.test(g)) {
      aligned += 1;if(!intervalStart||position!==lastCovered+1){if(intervalStart)covered.push(intervalStart,lastCovered);intervalStart=position;}lastCovered=position;baseByPosition?.set(position, { reference: g, query: q });
      if (q !== g) {if(mutations===EMPTY_MUTATIONS)mutations=[];mutations.push({ position, reference: g, alternate: q });}
    }
  }
  if(intervalStart)covered.push(intervalStart,lastCovered);
  const cdr3 = cleanCdr3(row.cdr3 || row.junction || ""),[cdr3HashA,cdr3HashB]=cdr3Fingerprint(cdr3);
  return { ordinal, unit, subject: (row.subject_id || "unassigned-subject").trim() || "unassigned-subject", vCall, jCall: topCall(row.j_call), cdr3HashA,cdr3HashB,junctionLength: cdr3.length, aligned, rate: aligned ? mutations.length / aligned : 1, mutations, coverageStart:covered[0]??0,coverageEnd:covered[covered.length-1]??0,coverageBreaks:covered.length>2?covered:undefined,baseByPosition };
}

function eventKey(event: Pick<EventEvidence, "position" | "reference" | "alternate">): string { return `${event.position}:${event.reference}>${event.alternate}`; }
function eventAt(item: Observation, event: Pick<EventEvidence, "position" | "reference" | "alternate">): boolean {
  const mutation = item.mutations.find(value=>value.position===event.position); return mutation?.reference === event.reference && mutation.alternate === event.alternate;
}
function coversEvent(item: Observation, event: Pick<EventEvidence, "position" | "reference">): boolean { const intervals=item.coverageBreaks;if(!intervals)return event.position>=item.coverageStart&&event.position<=item.coverageEnd;for(let index=0;index<intervals.length;index+=2)if(event.position>=intervals[index]&&event.position<=intervals[index+1])return true;return false; }
function setIntersectionSize(left: Set<number>, right: Set<number>): number { let result = 0; for (const value of left) if (right.has(value)) result += 1; return result; }
function isSubset(left: CandidateProposal | MissingAlleleCandidate, right: CandidateProposal | MissingAlleleCandidate): boolean {
  const rightEvents = "substitutions" in right ? right.substitutions : right.events; const leftEvents = "substitutions" in left ? left.substitutions : left.events;
  const rightKeys = new Set(rightEvents.map(eventKey)); return leftEvents.every((event) => rightKeys.has(eventKey(event)));
}

function mutationNullProbability(rate: number, hotspot: boolean): number {
  return Math.min(0.25, Math.max(1e-5, rate / 3) * (hotspot ? 5 : 1));
}

function diversity(items: Observation[]) {
  return {
    jCalls: new Set(items.map((item) => item.jCall).filter(Boolean)),
    junctionLengths: new Set(items.map((item) => item.junctionLength).filter((value) => value > 0)),
    cdr3s: new Set(items.filter(item=>item.junctionLength>0).map((item) => `${item.cdr3HashA.toString(36)}:${item.cdr3HashB.toString(36)}`)),
    subjects: new Set(items.map((item) => item.subject).filter(Boolean)),
  };
}

export class MissingAlleleAccumulator {
  private readonly options: MissingAlleleOptions;
  private readonly best = new Map<number, Observation>();
  private inputRecords = 0;
  private eligibleRecords = 0;

  constructor(options: Partial<MissingAlleleOptions> = {}) {
    this.options = { ...DEFAULT_MISSING_ALLELE_OPTIONS, ...options, unit: "lineage" };
  }

  add(row: Record<string, string>, ordinal: number, lineageId = 0) {
    this.inputRecords += 1; const unit = lineageUnit(lineageId); const item = parseObservation(row, ordinal, unit); if (!item) return;
    if (item.aligned < this.options.minimumAlignedBases || item.rate > this.options.maximumShmRate) return;
    this.eligibleRecords += 1;const existing = this.best.get(unit);
    if (!existing || item.rate < existing.rate || (item.rate === existing.rate && item.aligned > existing.aligned) || (item.rate === existing.rate && item.aligned === existing.aligned && item.ordinal < existing.ordinal)) this.best.set(unit, item);
  }

  prepareValidation(referenceFasta = ""): MissingAlleleValidator {
    const references = parseFasta(referenceFasta); const allKnown = new Set([...references.values()].map(cleanReference)); const groups = new Map<string, Observation[]>();
    for (const item of this.best.values()) { const group = groups.get(item.vCall) ?? []; group.push(item); groups.set(item.vCall, group); }
    const proposals: CandidateProposal[] = []; let proposalTruncations = 0;
    for (const [vCall, items] of groups) {
      if (items.length < this.options.minimumCoveredUnits) continue;
      const parent = references.get(vCall) || ""; if (!parent) continue;
      const coverageDelta=new Int32Array(parent.length+2);const changes = new Map<string, { position:number;reference: string; alternate: string; supporting: Observation[] }>();let summedRates=0;
      for(const item of items){summedRates+=item.rate;if(item.coverageBreaks){for(let interval=0;interval<item.coverageBreaks.length;interval+=2){const start=Math.max(1,item.coverageBreaks[interval]),end=Math.min(parent.length,item.coverageBreaks[interval+1]);if(start<=end){coverageDelta[start]+=1;coverageDelta[end+1]-=1;}}}else{const start=Math.max(1,item.coverageStart),end=Math.min(parent.length,item.coverageEnd);if(start<=end){coverageDelta[start]+=1;coverageDelta[end+1]-=1;}}for(const mutation of item.mutations){if(mutation.position<1||mutation.position>parent.length||parent[mutation.position-1]!==mutation.reference)continue;const key=eventKey(mutation);const found=changes.get(key)??{...mutation,supporting:[]};found.supporting.push(item);changes.set(key,found);}}
      const coverageByPosition=new Uint32Array(parent.length+1);let runningCoverage=0;for(let position=1;position<=parent.length;position+=1){runningCoverage+=coverageDelta[position];coverageByPosition[position]=runningCoverage;}
      const passing: DiscoveryEvent[] = [];
      for (const change of changes.values()) {
        const coverage=coverageByPosition[change.position];const support=change.supporting.length;const fraction=support/Math.max(1,coverage);const background=Math.max(0,(summedRates-change.supporting.reduce((sum,item)=>sum+1/Math.max(1,item.aligned),0))/items.length);
        const hotspot=hotspotAt(parent,change.position-1);const pValue=binomialSurvival(support,coverage,mutationNullProbability(background,hotspot));
        if(support>=this.options.minimumIndependentUnits&&coverage>=this.options.minimumCoveredUnits&&fraction>=this.options.minimumAlleleFraction&&pValue<=this.options.maximumPValue)passing.push({position:change.position,reference:change.reference,alternate:change.alternate,support,coverage,fraction,pValue,hotspot,supportingUnits:new Set(change.supporting.map(item=>item.unit))});
      }
      passing.sort((left, right) => left.pValue - right.pValue || right.fraction - left.fraction || right.support - left.support || left.position - right.position);
      const signatures = new Set<string>(); const vProposals: CandidateProposal[] = [];const itemByUnit=new Map(items.map(item=>[item.unit,item]));
      for (const seed of passing) {
        const events:DiscoveryEvent[] = [seed];let currentSupport=new Set(seed.supportingUnits);
        for (const candidate of passing) {
          if (events.length >= this.options.maximumCandidateSnps || events.some((event) => eventKey(event) === eventKey(candidate))) continue;
          let currentEligible=0,candidateEligible=0,both=0;for(const unit of currentSupport){const item=itemByUnit.get(unit)!;if(coversEvent(item,candidate)){currentEligible+=1;if(candidate.supportingUnits.has(unit))both+=1;}}for(const unit of candidate.supportingUnits){const item=itemByUnit.get(unit)!;if(events.every(event=>coversEvent(item,event)))candidateEligible+=1;}
          if (both < this.options.minimumIndependentUnits || both / Math.max(1, currentEligible) < this.options.minimumLinkedFraction || both / Math.max(1, candidateEligible) < this.options.minimumLinkedFraction) continue;
          events.push(candidate);currentSupport=new Set([...currentSupport].filter(unit=>candidate.supportingUnits.has(unit)));
        }
        events.sort((left, right) => left.position - right.position); const signature = events.map(eventKey).join("|"); if (signatures.has(signature)) continue; signatures.add(signature);
        const covered = items.filter((item) => events.every((event) => coversEvent(item, event))); const supporting = [...currentSupport].map(unit=>itemByUnit.get(unit)!); const fraction = supporting.length / Math.max(1, covered.length); const groupsSeen = diversity(supporting);
        if (supporting.length < this.options.minimumIndependentUnits || covered.length < this.options.minimumCoveredUnits || fraction < this.options.minimumAlleleFraction || groupsSeen.jCalls.size < this.options.minimumDistinctJCalls || groupsSeen.junctionLengths.size < this.options.minimumDistinctJunctionLengths || groupsSeen.cdr3s.size < this.options.minimumDistinctCdr3s) continue;
        const meanBackgroundShm = items.reduce((sum, item) => sum + Math.max(0, item.mutations.length - events.filter((event) => eventAt(item, event)).length) / Math.max(1, item.aligned), 0) / items.length;
        const nullProbability = events.reduce((probability, event) => probability * mutationNullProbability(meanBackgroundShm, event.hotspot), 1); const pValue = binomialSurvival(supporting.length, covered.length, nullProbability); if (pValue > this.options.maximumPValue) continue;
        const sequence = [...parent]; for (const event of events) if (sequence[event.position - 1] === event.reference) sequence[event.position - 1] = event.alternate; const candidateSequence = sequence.join(""); if (allKnown.has(candidateSequence)) continue;
        const storedEvents:EventEvidence[]=events.map(event=>({position:event.position,reference:event.reference,alternate:event.alternate,support:event.support,coverage:event.coverage,fraction:event.fraction,pValue:event.pValue,hotspot:event.hotspot}));
        vProposals.push({ vCall, sequence: candidateSequence, events:storedEvents, preliminarySupport: new Set(supporting.map((item) => item.unit)), meanBackgroundShm, nullProbability, pValue });
      }
      vProposals.sort((left, right) => left.pValue - right.pValue || right.preliminarySupport.size - left.preliminarySupport.size || right.events.length - left.events.length);
      const maximal = vProposals.filter((candidate, index) => !vProposals.some((other, otherIndex) => otherIndex !== index && other.events.length > candidate.events.length && isSubset(candidate, other) && setIntersectionSize(candidate.preliminarySupport, other.preliminarySupport) / Math.max(1, candidate.preliminarySupport.size) >= this.options.minimumLinkedFraction));
      if (maximal.length > MAX_PROPOSALS_PER_V) proposalTruncations += maximal.length - MAX_PROPOSALS_PER_V;for(const proposal of maximal.slice(0,MAX_PROPOSALS_PER_V)){proposal.preliminarySupport.clear();proposals.push(proposal);}
    }
    return new MissingAlleleValidator(this.options, this.inputRecords, this.eligibleRecords, this.best, groups.size, proposals, proposalTruncations, allKnown);
  }
}

export class MissingAlleleValidator {
  private readonly groups = new Map<string, ProposalGroup>();
  private readonly options: MissingAlleleOptions;
  private readonly inputRecords: number;
  private readonly eligibleRecords: number;
  private readonly discoveryUnits: Map<number, Observation>;
  private readonly vAllelesTested: number;
  private readonly proposals: CandidateProposal[];
  private readonly proposalTruncations: number;
  private readonly allKnown: Set<string>;

  constructor(
    options: MissingAlleleOptions,
    inputRecords: number,
    eligibleRecords: number,
    discoveryUnits: Map<number, Observation>,
    vAllelesTested: number,
    proposals: CandidateProposal[],
    proposalTruncations: number,
    allKnown: Set<string>,
  ) {
    this.options=options;this.inputRecords=inputRecords;this.eligibleRecords=eligibleRecords;this.discoveryUnits=discoveryUnits;this.vAllelesTested=vAllelesTested;this.proposals=proposals;this.proposalTruncations=proposalTruncations;this.allKnown=allKnown;
    for(const proposal of proposals){const group:ProposalGroup=this.groups.get(proposal.vCall)??{proposals:[],observations:[],flags:new Uint8Array(0),mixedVCall:new Uint8Array(0)};group.proposals.push(proposal);this.groups.set(proposal.vCall,group);}
    for(const item of discoveryUnits.values()){const group=this.groups.get(item.vCall);if(!group)continue;item.validationIndex=group.observations.length;group.observations.push(item);}
    for(const group of this.groups.values()){group.flags=new Uint8Array(group.observations.length*group.proposals.length);group.mixedVCall=new Uint8Array(group.observations.length);}
  }

  add(row: Record<string, string>, ordinal: number, lineageId = 0) {
    const unit = lineageUnit(lineageId);const discovery=this.discoveryUnits.get(unit);if(!discovery)return;const group=this.groups.get(discovery.vCall);const unitIndex=discovery.validationIndex;if(!group||unitIndex===undefined)return;
    const observedVCall = topCall(row.v_call); if (observedVCall && observedVCall !== discovery.vCall) { group.mixedVCall[unitIndex]=1; return; }
    const item = parseObservation(row, ordinal, unit, true); if (!item) return;
    group.proposals.forEach((proposal, proposalIndex) => {
      let allCovered = true, allAlternate = true, referencePresent = false, conflictPresent = false;
      for (const event of proposal.events) {
        const base = item.baseByPosition?.get(event.position); if (!base) { allCovered = false; allAlternate = false; continue; }
        if (base.reference !== event.reference) { allCovered = false; allAlternate = false; conflictPresent = true; continue; }
        if (base.query === event.reference) { referencePresent = true; allAlternate = false; }
        else if (base.query !== event.alternate) { conflictPresent = true; allAlternate = false; }
      }
      const offset=unitIndex*group.proposals.length+proposalIndex;
      if (allCovered) group.flags[offset] |= FLAG_JOINT_COVERAGE;
      if (allCovered && allAlternate) group.flags[offset] |= FLAG_JOINT_ALTERNATE;
      if (referencePresent) group.flags[offset] |= FLAG_REFERENCE_PRESENT;
      if (conflictPresent) group.flags[offset] |= FLAG_CONFLICT_PRESENT;
    });
  }

  finish(): MissingAlleleDashboard {
    const candidates: MissingAlleleCandidate[] = []; const candidateSupportSets = new Map<MissingAlleleCandidate, Set<number>>(); let referenceVetoedLineagePatterns = 0, conflictingLineagePatterns = 0;
    for (const [vCall, group] of this.groups) {
      group.proposals.forEach((proposal, proposalIndex) => {
        const supporting: Observation[] = []; let coveredUnits = 0, referencePresentUnits = 0, conflictingUnits = 0;
        for (let unitIndex=0;unitIndex<group.observations.length;unitIndex+=1) {
          const item=group.observations[unitIndex],flags = group.flags[unitIndex*group.proposals.length+proposalIndex],mixedVCall=Boolean(group.mixedVCall[unitIndex]); if ((flags & FLAG_JOINT_COVERAGE) && !mixedVCall) coveredUnits += 1;
          const hadAlternate = Boolean(flags & FLAG_JOINT_ALTERNATE); const hadReference = Boolean(flags & FLAG_REFERENCE_PRESENT); const hadConflict = Boolean(flags & FLAG_CONFLICT_PRESENT) || mixedVCall;
          if (hadAlternate && hadReference) referencePresentUnits += 1;
          if ((flags & FLAG_CONFLICT_PRESENT) && !mixedVCall) conflictingUnits += 1;
          if (hadAlternate && !hadReference && !hadConflict) supporting.push(item);
        }
        referenceVetoedLineagePatterns += referencePresentUnits; conflictingLineagePatterns += conflictingUnits;
        const alleleFraction = supporting.length / Math.max(1, coveredUnits); const groupsSeen = diversity(supporting); const candidateEventKeys = new Set(proposal.events.map(eventKey));
        const nearGermlineUnits = supporting.filter((item) => item.mutations.filter((mutation) => !candidateEventKeys.has(eventKey(mutation))).length <= 2).length;
        const pValue = binomialSurvival(supporting.length, coveredUnits, proposal.nullProbability);
        const otherAlternateFraction=conflictingUnits/Math.max(1,coveredUnits);
        if (supporting.length < this.options.minimumIndependentUnits || coveredUnits < this.options.minimumCoveredUnits || alleleFraction < this.options.minimumAlleleFraction || otherAlternateFraction > this.options.maximumOtherAlternateFraction || pValue > this.options.maximumPValue || groupsSeen.jCalls.size < this.options.minimumDistinctJCalls || groupsSeen.junctionLengths.size < this.options.minimumDistinctJunctionLengths || groupsSeen.cdr3s.size < this.options.minimumDistinctCdr3s || nearGermlineUnits < this.options.minimumNearGermlineUnits || this.allKnown.has(proposal.sequence)) return;
        const units = supporting.map(displayUnit).sort(); const substitutions:CandidateSubstitution[] = proposal.events.map((event) => ({ position:event.position,reference:event.reference,alternate:event.alternate,hotspot:event.hotspot,support: supporting.length, coverage: coveredUnits, fraction: alleleFraction, pValue }));
        const candidate:MissingAlleleCandidate = { id: "", vCall, parentAllele: vCall, substitutions, independentUnits: supporting.length, coveredUnits, alleleFraction, pValue, distinctJCalls: groupsSeen.jCalls.size, distinctJunctionLengths: groupsSeen.junctionLengths.size, distinctCdr3s: groupsSeen.cdr3s.size, distinctSubjects: groupsSeen.subjects.size, nearGermlineUnits, referencePresentUnits, conflictingUnits, otherAlternateFraction, meanBackgroundShm: proposal.meanBackgroundShm, sequence: proposal.sequence, supportingUnits: units.slice(0, SUPPORTING_UNIT_EXPORT_LIMIT), supportingUnitsTruncated: units.length > SUPPORTING_UNIT_EXPORT_LIMIT, caution: "Diagnostic candidate only. Confirm with genotype-aware germline inference and independent data before adding it to a reference set." };
        candidates.push(candidate);candidateSupportSets.set(candidate,new Set(supporting.map(item=>item.unit)));
      });
    }
    candidates.sort((left, right) => left.pValue - right.pValue || right.independentUnits - left.independentUnits || right.substitutions.length - left.substitutions.length);
    const maximal = candidates.filter((candidate, index) => !candidates.some((other, otherIndex) => otherIndex !== index && other.vCall === candidate.vCall && other.substitutions.length > candidate.substitutions.length && isSubset(candidate, other) && setIntersectionSize(candidateSupportSets.get(candidate)!,candidateSupportSets.get(other)!) / Math.max(1, candidate.independentUnits) >= this.options.minimumLinkedFraction));
    maximal.forEach((candidate, index) => { candidate.id = `${candidate.vCall.replace(/[^A-Za-z0-9_.-]/g, "_")}_candidate_${index + 1}`; });
    const warnings: string[] = [];
    if (referenceVetoedLineagePatterns) warnings.push(`${referenceVetoedLineagePatterns.toLocaleString()} lineage–candidate pattern${referenceVetoedLineagePatterns === 1 ? " was" : "s were"} vetoed because another covered member of the same lineage contained a parent-reference nucleotide at a proposed site.`);
    if (conflictingLineagePatterns) warnings.push(`${conflictingLineagePatterns.toLocaleString()} lineage–candidate pattern${conflictingLineagePatterns === 1 ? " contained" : "s contained"} a third nucleotide state; candidates above the configured cross-lineage fraction were rejected.`);
    if (this.proposalTruncations) warnings.push(`${this.proposalTruncations.toLocaleString()} lower-ranked preliminary pattern${this.proposalTruncations === 1 ? " was" : "s were"} omitted by the 64-pattern-per-V browser memory guard.`);
    if (maximal.length) warnings.push(`${maximal.length} linked low-SHM haplotype${maximal.length === 1 ? "" : "s"} passed the two-pass diagnostic. These are referral candidates, not inferred genotypes.`);
    return { mode: "lineage", validationPasses: 2, inputRecords: this.inputRecords, eligibleRecords: this.eligibleRecords, independentUnits: this.discoveryUnits.size, vAllelesTested: this.vAllelesTested, candidatePatternsTested: this.proposals.length, referenceVetoedLineagePatterns, conflictingLineagePatterns, proposalTruncations: this.proposalTruncations, candidates: maximal, warnings };
  }
}

export function candidateFasta(dashboard: MissingAlleleDashboard): string {
  return dashboard.candidates.map((candidate) => `>${candidate.id} parent=${candidate.parentAllele} lineages=${candidate.independentUnits} cdr3s=${candidate.distinctCdr3s} subjects=${candidate.distinctSubjects} reference_vetoes=${candidate.referencePresentUnits} fraction=${candidate.alleleFraction.toFixed(4)}\n${candidate.sequence}\n`).join("");
}

export function candidateEvidenceRows(dashboard: MissingAlleleDashboard): Array<Record<string, string | number>> {
  return dashboard.candidates.map((candidate) => ({ candidate_id: candidate.id, v_call: candidate.vCall, parent_allele: candidate.parentAllele, linked_substitutions: candidate.substitutions.map((item) => `${item.reference}${item.position}${item.alternate}`).join(";"), independent_lineages: candidate.independentUnits, jointly_covered_lineages: candidate.coveredUnits, allele_fraction: candidate.alleleFraction, screening_tail_probability: candidate.pValue, distinct_j_calls: candidate.distinctJCalls, distinct_cdr3_lengths: candidate.distinctJunctionLengths, distinct_cdr3_sequences: candidate.distinctCdr3s, distinct_subjects: candidate.distinctSubjects, near_germline_lineages: candidate.nearGermlineUnits, reference_state_vetoes: candidate.referencePresentUnits, other_alternate_lineages: candidate.conflictingUnits, other_alternate_fraction: candidate.otherAlternateFraction, mean_background_shm: candidate.meanBackgroundShm, hotspot_substitutions: candidate.substitutions.filter((item) => item.hotspot).length, supporting_lineages: candidate.supportingUnits.join(";"), supporting_lineages_truncated: candidate.supportingUnitsTruncated ? 1 : 0 }));
}
