export interface MissingAlleleOptions {
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
}

export const DEFAULT_MISSING_ALLELE_OPTIONS: MissingAlleleOptions = {
  unit: "lineage", minimumIndependentUnits: 6, minimumCoveredUnits: 20, minimumAlleleFraction: 0.15,
  maximumShmRate: 0.12, minimumAlignedBases: 180, maximumCandidateSnps: 6, maximumPValue: 1e-6,
  minimumDistinctJCalls: 3, minimumDistinctJunctionLengths: 3,
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
  meanBackgroundShm: number;
  sequence: string;
  caution: string;
}

export interface MissingAlleleDashboard {
  mode: "lineage" | "record";
  inputRecords: number;
  eligibleRecords: number;
  independentUnits: number;
  vAllelesTested: number;
  candidates: MissingAlleleCandidate[];
  warnings: string[];
}

interface Observation {
  ordinal: number;
  unit: string;
  vCall: string;
  jCall: string;
  junctionLength: number;
  aligned: number;
  rate: number;
  referenceSequence: string;
  mutations: Array<[number,string]>;
  covered: Array<[number,number]>;
}

function topCall(value: string): string { return value.split(",")[0]?.trim() ?? ""; }
function cleanSequence(value: string): string { return value.toUpperCase().replace(/[^ACGTN-]/g, "N"); }
function parseFasta(value: string): Map<string,string> {
  const result=new Map<string,string>(); let name="",sequence="";
  const commit=()=>{if(name&&sequence)result.set(name,sequence.replace(/\s/g,"").toUpperCase());};
  for(const line of value.split(/\r?\n/)){if(line.startsWith(">")){commit();name=line.slice(1).trim().split(/\s+/)[0];sequence="";}else sequence+=line.trim();} commit(); return result;
}

function logGamma(value:number):number {
  const coefficients=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,.001208650973866179,-.000005395239384953];
  let x=value,y=value,tmp=x+5.5; tmp-=(x+.5)*Math.log(tmp); let series=1.000000000190015;
  for(const coefficient of coefficients){y+=1;series+=coefficient/y;} return -tmp+Math.log(2.5066282746310005*series/x);
}

function betaFraction(a:number,b:number,x:number):number {
  const max=200,epsilon=3e-12,tiny=1e-300; let qab=a+b,qap=a+1,qam=a-1,c=1,d=1-qab*x/qap;
  if(Math.abs(d)<tiny)d=tiny; d=1/d; let h=d;
  for(let m=1;m<=max;m++){const m2=2*m; let aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=1+aa*d;if(Math.abs(d)<tiny)d=tiny;c=1+aa/c;if(Math.abs(c)<tiny)c=tiny;d=1/d;h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));d=1+aa*d;if(Math.abs(d)<tiny)d=tiny;c=1+aa/c;if(Math.abs(c)<tiny)c=tiny;d=1/d;const delta=d*c;h*=delta;if(Math.abs(delta-1)<epsilon)break;}
  return h;
}

function regularizedBeta(x:number,a:number,b:number):number {
  if(x<=0)return 0;if(x>=1)return 1;
  const bt=Math.exp(logGamma(a+b)-logGamma(a)-logGamma(b)+a*Math.log(x)+b*Math.log(1-x));
  return x<(a+1)/(a+b+2)?bt*betaFraction(a,b,x)/a:1-bt*betaFraction(b,a,1-x)/b;
}

export function binomialSurvival(successes:number,trials:number,probability:number):number {
  if(successes<=0)return 1;if(successes>trials)return 0;const p=Math.min(1-1e-12,Math.max(1e-12,probability));
  return Math.min(1,Math.max(0,regularizedBeta(p,successes,trials-successes+1)));
}

function hotspotAt(sequence:string,position:number):boolean {
  const context=sequence.slice(Math.max(0,position-2),Math.min(sequence.length,position+3));
  return /[AT]A[CT]|[AG]T[AT]|[AT]GC[CT]|[AG]GC[AT]/.test(context);
}

function observation(row:Record<string,string>,ordinal:number,unit:string):Observation|null {
  const vCall=topCall(row.v_call);const query=cleanSequence(row.v_sequence_alignment||"");const germline=cleanSequence(row.v_germline_alignment||"");
  if(!vCall||!query||!germline)return null;const length=Math.min(query.length,germline.length);const mutations:Array<[number,string]>=[];const covered:Array<[number,number]>=[];let position=0,aligned=0,intervalStart=0,lastCovered=0;
  let referenceSequence="";
  for(let column=0;column<length;column++){const q=query[column],g=germline[column];if(g!=="-"){position+=1;referenceSequence+=g;}if(/^[ACGT]$/.test(q)&&/^[ACGT]$/.test(g)){aligned+=1;if(!intervalStart||position!==lastCovered+1){if(intervalStart)covered.push([intervalStart,lastCovered]);intervalStart=position;}lastCovered=position;if(q!==g)mutations.push([position,`${g}>${q}`]);}}
  if(intervalStart)covered.push([intervalStart,lastCovered]);
  return {ordinal,unit,vCall,jCall:topCall(row.j_call),junctionLength:(row.junction||row.cdr3||"").replace(/[-.\s]/g,"").length,aligned,rate:aligned?mutations.length/aligned:1,referenceSequence,mutations,covered};
}

function covers(item:Observation,position:number):boolean{return item.covered.some(([start,end])=>position>=start&&position<=end);}
function mutationAt(item:Observation,position:number):string|undefined{return item.mutations.find(([site])=>site===position)?.[1];}

export class MissingAlleleAccumulator {
  private readonly options:MissingAlleleOptions;
  private readonly best=new Map<string,Observation>();
  private inputRecords=0;
  private eligibleRecords=0;

  constructor(options:Partial<MissingAlleleOptions>={}){this.options={...DEFAULT_MISSING_ALLELE_OPTIONS,...options};}

  add(row:Record<string,string>,ordinal:number,lineageId=0){
    this.inputRecords+=1;const unit=this.options.unit==="lineage"?`lineage:${lineageId||`record:${ordinal}`}`:`record:${ordinal}`;const item=observation(row,ordinal,unit);if(!item)return;
    if(item.aligned<this.options.minimumAlignedBases||item.rate>this.options.maximumShmRate)return;this.eligibleRecords+=1;
    const key=`${unit}\u0000${item.vCall}`;const existing=this.best.get(key);if(!existing||item.rate<existing.rate||(item.rate===existing.rate&&item.aligned>existing.aligned))this.best.set(key,item);
  }

  finish(referenceFasta=""):MissingAlleleDashboard {
    const references=parseFasta(referenceFasta);const allKnown=new Set([...references.values()].map((sequence)=>sequence.replace(/[^ACGTN]/g,"")));
    const groups=new Map<string,Observation[]>();for(const item of this.best.values()){const group=groups.get(item.vCall)??[];group.push(item);groups.set(item.vCall,group);}
    const candidates:MissingAlleleCandidate[]=[];
    for(const [vCall,items] of groups){
      if(items.length<this.options.minimumCoveredUnits)continue;const eventCounts=new Map<string,{position:number;change:string;support:number;coverage:number;supporting:Observation[]}>();
      const positions=new Set<number>();items.forEach((item)=>item.mutations.forEach(([position])=>positions.add(position)));
      for(const position of positions){const coverage=items.filter((item)=>covers(item,position)).length;const changes=new Map<string,Observation[]>();for(const item of items){const change=mutationAt(item,position);if(change){const list=changes.get(change)??[];list.push(item);changes.set(change,list);}}
        for(const [change,supporting] of changes)eventCounts.set(`${position}:${change}`,{position,change,support:supporting.length,coverage,supporting});}
      const meanBackground=items.reduce((sum,item)=>sum+item.rate,0)/items.length;
      const passing=[...eventCounts.values()].filter((event)=>event.support>=this.options.minimumIndependentUnits&&event.coverage>=this.options.minimumCoveredUnits&&event.support/event.coverage>=this.options.minimumAlleleFraction&&binomialSurvival(event.support,event.coverage,Math.max(1e-4,meanBackground/3))<=this.options.maximumPValue)
        .sort((a,b)=>b.support/a.coverage-a.support/b.coverage||a.position-b.position);
      const consumed=new Set<string>();
      for(const seed of passing){const seedKey=`${seed.position}:${seed.change}`;if(consumed.has(seedKey))continue;const linked=passing.filter((event)=>event.position!==seed.position&&seed.supporting.filter((item)=>mutationAt(item,event.position)===event.change).length/seed.support>=.8).slice(0,this.options.maximumCandidateSnps-1);
        const events=[seed,...linked].sort((a,b)=>a.position-b.position);const supporting=items.filter((item)=>events.every((event)=>mutationAt(item,event.position)===event.change));const covered=items.filter((item)=>events.every((event)=>covers(item,event.position)));
        if(supporting.length<this.options.minimumIndependentUnits||supporting.length/Math.max(1,covered.length)<this.options.minimumAlleleFraction)continue;
        const jCalls=new Set(supporting.map((item)=>item.jCall).filter(Boolean));const lengths=new Set(supporting.map((item)=>item.junctionLength).filter(Boolean));
        if(jCalls.size<this.options.minimumDistinctJCalls||lengths.size<this.options.minimumDistinctJunctionLengths)continue;
        const parent=references.get(vCall)||references.get(vCall.replace(/\*.*$/, ""))||items[0].referenceSequence;const sequence=[...parent];const substitutions:CandidateSubstitution[]=events.map((event)=>{const [reference,alternate]=event.change.split(">");const support=supporting.filter((item)=>mutationAt(item,event.position)===event.change).length;const coverage=covered.filter((item)=>covers(item,event.position)).length;sequence[event.position-1]=alternate;return {position:event.position,reference,alternate,support,coverage,fraction:support/Math.max(1,coverage),pValue:binomialSurvival(support,coverage,Math.max(1e-4,meanBackground/3)),hotspot:hotspotAt(parent,event.position-1)};});
        const candidateSequence=sequence.join("");if(allKnown.has(candidateSequence))continue;events.forEach((event)=>consumed.add(`${event.position}:${event.change}`));
        const pValue=Math.max(...substitutions.map((item)=>item.pValue));const id=`${vCall.replace(/[^A-Za-z0-9_.-]/g,"_")}_candidate_${candidates.length+1}`;
        candidates.push({id,vCall,parentAllele:vCall,substitutions,independentUnits:supporting.length,coveredUnits:covered.length,alleleFraction:supporting.length/Math.max(1,covered.length),pValue,distinctJCalls:jCalls.size,distinctJunctionLengths:lengths.size,meanBackgroundShm:meanBackground,sequence:candidateSequence,
          caution:"Diagnostic candidate only. Confirm with genotype-aware germline inference and independent data before adding it to a reference set."});
      }
    }
    candidates.sort((a,b)=>a.pValue-b.pValue||b.independentUnits-a.independentUnits);
    const warnings:string[]=[];
    if(this.options.unit==="record")warnings.push("Record-level mode is exploratory: clonally expanded mutations can look like alleles. Lineage-level mode is the default.");
    if(candidates.length)warnings.push(`${candidates.length} recurrent low-SHM pattern${candidates.length===1?"":"s"} may be compatible with missing germline alleles. These are warnings, not inferred genotypes.`);
    return {mode:this.options.unit,inputRecords:this.inputRecords,eligibleRecords:this.eligibleRecords,independentUnits:this.best.size,vAllelesTested:groups.size,candidates,warnings};
  }
}

export function candidateFasta(dashboard:MissingAlleleDashboard):string{return dashboard.candidates.map((candidate)=>`>${candidate.id} parent=${candidate.parentAllele} units=${candidate.independentUnits} fraction=${candidate.alleleFraction.toFixed(4)}\n${candidate.sequence}\n`).join("");}

export function candidateEvidenceRows(dashboard:MissingAlleleDashboard):Array<Record<string,string|number>>{
  return dashboard.candidates.map((candidate)=>({candidate_id:candidate.id,v_call:candidate.vCall,parent_allele:candidate.parentAllele,substitutions:candidate.substitutions.map((item)=>`${item.reference}${item.position}${item.alternate}`).join(";"),independent_units:candidate.independentUnits,covered_units:candidate.coveredUnits,allele_fraction:candidate.alleleFraction,p_value:candidate.pValue,distinct_j_calls:candidate.distinctJCalls,distinct_junction_lengths:candidate.distinctJunctionLengths,mean_background_shm:candidate.meanBackgroundShm,hotspot_substitutions:candidate.substitutions.filter((item)=>item.hotspot).length}));
}
