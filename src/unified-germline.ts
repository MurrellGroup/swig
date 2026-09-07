import {PersonalizedGermlineAccumulator,DEFAULT_PERSONALIZED_GERMLINE_OPTIONS,proposeNovelCandidates,type Observation} from './personalized-germline.ts';
import {AlignedRateCalibration} from './shm-model/discovery.ts';
import {substitutionCompatible} from './shm-model/alignment.ts';
import {boundaryLogLikelihood} from './shm-model/boundary.ts';
import {integrate,splitEvidence,HAZARDS,type Row} from './shm-model/unified-mixture.ts';
import {transition} from './shm-model/unified-kernel.ts';
export interface UnifiedOptions {maximumCandidates:number;maximumIterations:number;splitSeed:number}
export const DEFAULT_UNIFIED_OPTIONS:UnifiedOptions={maximumCandidates:64,maximumIterations:300,splitSeed:0};
export type UnifiedGeneResult=Omit<ReturnType<typeof splitEvidence>,'evidence'> & {key:string;trainingLineages:number;testLineages:number;prior:{tau:number;naive:number};times:number[];weights:number[];proposalTruncated:boolean;seconds:number;evidence:Array<ReturnType<typeof splitEvidence>['evidence'][number]&{id:string;sequence:string;known:boolean}>};
export interface UnifiedDashboard {evidenceRule:{alpha:number;testsBySubject:Record<string,number>;description:string};version:1;mode:'joint-inherited-somatic-split';options:UnifiedOptions;scope:string;results:UnifiedGeneResult[];warnings:string[]}
export function inferUnifiedGermline(snapshot:ReturnType<PersonalizedGermlineAccumulator['researchSnapshot']>,options:UnifiedOptions=DEFAULT_UNIFIED_OPTIONS,onProgress?:(done:number,total:number)=>void):UnifiedDashboard {
for(const field of ['maximumCandidates','maximumIterations'] as const)if(!Number.isSafeInteger(options[field])||options[field]<1)throw new Error(field+' must be a positive integer');
const {observations,references}=snapshot;
// Split by lineage ID, before ANY discovery/calibration. Stable across input order.
const isTraining=(o:Observation)=>{let x=(o.lineageId^options.splitSeed)>>>0;x=Math.imul(x^(x>>>16),0x45d9f3b);x=Math.imul(x^(x>>>16),0x45d9f3b);return ((x^(x>>>16))>>>0)%2===0;};
const training=observations.filter(isTraining),calibration=new AlignedRateCalibration(training),groups=new Map<string,Observation[]>();
// Reference-defined compatibility components allow cross-gene competition.
const parents=references,roots=new Map(parents.map(p=>[p.index,p.index]));
const root=(n:number):number=>{while(roots.get(n)!==n)n=roots.get(n)!;return n;};
for(let i=0;i<parents.length;i++)for(let j=0;j<i;j++)if(parents[i].locus===parents[j].locus&&substitutionCompatible(parents[i].sequence,parents[j].sequence,8))roots.set(root(parents[i].index),root(parents[j].index));
for(const o of observations){const key=`${o.subjectId}|${o.locus}|${root(o.parent.index)}|${o.parent.sequence.length}`;const g=groups.get(key)??[];g.push(o);groups.set(key,g);}
const results:UnifiedGeneResult[]=[];let done=0;
for(const [key,items] of groups){
 onProgress?.(done++,groups.size);
 const train=items.filter(isTraining),test=items.filter(o=>!isTraining(o));if(train.length<4||test.length<4)continue;
 const known=references.filter(p=>root(p.index)===root(items[0].parent.index));
 const proposed=proposeNovelCandidates(train,{...DEFAULT_PERSONALIZED_GERMLINE_OPTIONS,maximumNovelCandidatesPerGene:options.maximumCandidates},calibration);
 const candidates=[...new Map([...known.map(p=>({id:p.names.join(','),sequence:p.sequence,known:true})),...proposed.proposals.map(p=>({id:p.candidate.id,sequence:p.candidate.sequence,known:false}))].map(c=>[c.sequence,c])).values()];
 if(candidates.every(c=>c.known))continue;
 // Coarsen only the final two template bases; raw boundary is marginalized
 // identically in inherited and somatic likelihood components.
 const fullLength=items[0].parent.sequence.length;const L=fullLength-2;
 const projected=[...new Map([...candidates].reverse().map(c=>[c.sequence.slice(0,L),c])).values()].map(c=>({...c,sequence:c.sequence.slice(0,L)}));
 // Restore known status when a novel terminal extension projects onto a reference.
 for(const c of projected)if(known.some(k=>k.sequence.slice(0,L)===c.sequence))c.known=true;
 if(projected.every(c=>c.known))continue;
 const prior=calibration.componentModel(train[0],new Set(known.map(p=>p.gene))),times=[0,...Array.from({length:16},(_,j)=>-Math.log(1-(j+.5)/16)*prior.tau)],weights=[prior.naive,...new Array(16).fill((1-prior.naive)/16)];
 // Cache CTMC probabilities by frozen five-base context and exposure state.
 const cache=new Map<string,Float64Array[]>();
 const tables=projected.map(c=>Array.from(c.sequence,(_,p)=>{const context=('NN'+c.sequence+'NN').slice(p,p+5);let table=cache.get(context);if(!table){table=times.map(t=>transition(context,t));cache.set(context,table);}return table;}));
 const logTables=projected.map((c,ci)=>Array.from(c.sequence,(base,p)=>Float64Array.from({length:times.length*4},(_,i)=>Math.log('ACGT'.includes(base)?tables[ci][p][Math.floor(i/4)]['ACGT'.indexOf(base)*4+i%4]:[0,1,2,3].reduce((sum,a)=>sum+tables[ci][p][Math.floor(i/4)][a*4+i%4]/4,0)))));
 const differences=projected.map((c,ci)=>Array.from({length:fullLength-12},(_,p)=>p).filter(p=>c.sequence[p]!==projected[0].sequence[p]||tables[ci][p]!==tables[0][p]));
 const boundaryKeys=projected.map(c=>c.sequence.slice(fullLength-14));
 const build=(obs:Observation[])=>{
  const rows=HAZARDS.map(()=>[] as Row[]);
  for(const o of obs){
   const logs=new Float64Array(projected.length*times.length),bases=new Int8Array(fullLength).fill(-1),baseline=new Float64Array(times.length);
   for(let i=0;i<o.positions.length;i++){const p=o.positions[i]-1,b='ACGT'.indexOf(o.query[i]);if(p>=fullLength-12||b<0)continue;bases[p]=b;for(let k=0;k<times.length;k++)baseline[k]+=logTables[0][p][k*4+b];}
   const tails=new Map<string,Float64Array>();
   for(let c=0;c<projected.length;c++){
    let tail=tails.get(boundaryKeys[c]);
    if(!tail){tail=new Float64Array(times.length);if(o.terminal)for(let k=0;k<times.length;k++)tail[k]=boundaryLogLikelihood(projected[c].sequence+'NN',o.terminal.query,o.terminal.start,(position,base)=>{const p=position-1,b='ACGT'.indexOf(base);return p>=L||b<0?Math.log(.25):logTables[c][p][k*4+b];},prior.trimming,prior.junction);tails.set(boundaryKeys[c],tail);}
    for(let k=0;k<times.length;k++){let value=baseline[k]+tail[k];for(const p of differences[c]){const b=bases[p];if(b>=0)value+=logTables[c][p][k*4+b]-logTables[0][p][k*4+b];}logs[c*times.length+k]=value;}
   }
   const integrated=integrate(logs,times,weights);for(let h=0;h<rows.length;h++)rows[h].push(integrated[h]);
  }return rows;
 };
 const start=performance.now();
 const fit=splitEvidence(build(train),build(test),projected.map(c=>c.known),options.maximumIterations);
 results.push({key,trainingLineages:train.length,testLineages:test.length,prior,times,weights,proposalTruncated:proposed.truncated,seconds:(performance.now()-start)/1000,...fit,evidence:fit.evidence.map(e=>({...projected[e.candidate],...e}))});

}

const testsBySubject:Record<string,number>={};for(const r of results){const subject=r.key.split('|')[0];testsBySubject[subject]=(testsBySubject[subject]??0)+r.evidence.length;}
return {evidenceRule:{alpha:.05,testsBySubject,description:'log evidence >= log(M/0.05), conditional on the fitted observation model; not an established real-data error guarantee'},version:1,mode:"joint-inherited-somatic-split",options,scope:"V sequence excluding final two bases; fixed training-calibrated rearrangement nuisance",results,warnings:["Experimental model: scores are conditional on the mutation, exposure and rearrangement assumptions; real-data FDR is not established.","Candidate generation and calibration use the training half only. Small/rare signals may be unidentifiable in the held-out half.","Final two V bases are unresolved and are not inferred.","One lowest-current-SHM representative per lineage is used; this sampling choice and imperfect lineage calls remain limitations.",...(results.some(r=>r.proposalTruncated)?["Some training proposal queues reached the candidate budget; absence from the tested set is not evidence of absence."]:[])]};
}