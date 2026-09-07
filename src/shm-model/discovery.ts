import { trimWeights,logSum,profileBoundaryAllele,terminalChangeIdentifiable } from "./boundary.ts";
import { hs5fRate } from './hs5f.ts';
import { globalCoordinateMap } from './alignment.ts';

export interface EvidenceObservation {
 ordinal:number; subjectId:string; locus:string; gene:string;
 parent:{index:number;sequence:string;gene:string}; positions:Uint16Array; query:string; terminal?:{start:number;query:string}; junction?:string;
}
export interface Change {position:number;reference:string;alternate:string}
export interface RatePrior {shape:number;rate:number;genes:number}
const clamp=(x:number)=>Math.max(1e-12,Math.min(1-1e-12,x));
export function baseAt(o:EvidenceObservation,p:number):string|undefined {
 let low=0,high=o.positions.length-1;
 while(low<=high){const mid=(low+high)>>>1;if(o.positions[mid]===p)return o.query[mid];if(o.positions[mid]<p)low=mid+1;else high=mid-1;}
 return undefined;
}
export function otherBurden(o:EvidenceObservation,excluded:ReadonlySet<number>,error=.001):number {
 let n=0,k=0,w=0;
 for(let i=0;i<o.positions.length;i++){
  const p=o.positions[i];if(excluded.has(p)||p>o.parent.sequence.length-12)continue;
  n++;k+=Number(o.query[i]!==o.parent.sequence[p-1]);w+=hs5fRate(o.parent.sequence,p-1);
 }
 return n?Math.min(1,(Math.max(0,k-n*error)+.5)/Math.max(1,w)):0.02;
}
const family=(gene:string)=>gene.match(/^[A-Z]+\d+/)?.[0]??gene;
/** Per-subject, explicit family-anchor homology; no tested physical gene calibrates itself. */
export class AlignedRateCalibration {
 private burdenTotals=new Map<string,Map<string,{n:number;zeros:number;k:number;w:number;bases:number}>>();
 private burdenCache=new Map<string,{naive:number;beta:number}>();
 private error:number;
 private junctionCounts=new Map<string,Map<string,Float64Array>>();
 private transitionCounts=new Map<string,Map<string,Float64Array>>();
 private trimCounts=new Map<string,Map<string,Float64Array>>();
 private boundaryCache=new Map<string,number[]>();
 private anchors=new Map<string,string>();
 private cache=new Map<string,RatePrior>();
 private maps=new Map<string,Int32Array>();
 private events=new Map<string,Map<string,{k:number,e:number}>>();
 constructor(observations:readonly EvidenceObservation[],error=.001){
  this.error=error;
  for(const o of observations){
   let k=0,w=0,bases=0;
   for(let i=0;i<o.positions.length;i++){const p=o.positions[i];if(p>o.parent.sequence.length-12)continue;k+=Number(o.query[i]!==o.parent.sequence[p-1]);w+=hs5fRate(o.parent.sequence,p-1);bases++;}
   const key=`${o.subjectId}|${o.locus}`,genes=this.burdenTotals.get(key)??new Map();
   const v=genes.get(o.gene)??{n:0,zeros:0,k:0,w:0,bases:0};v.n++;v.zeros+=Number(k===0);v.k+=k;v.w+=w;v.bases+=bases;genes.set(o.gene,v);this.burdenTotals.set(key,genes);
  }
  for(const o of observations){
   const key=`${o.subjectId}|${o.locus}`,genes=this.junctionCounts.get(key)??new Map<string,Float64Array>();
   const counts=genes.get(o.gene)??new Float64Array(4);for(const b of o.junction??""){const i="ACGT".indexOf(b);if(i>=0)counts[i]++;}genes.set(o.gene,counts);this.junctionCounts.set(key,genes);
  }
  for(const o of observations){
   const key=`${o.subjectId}|${o.locus}`,genes=this.transitionCounts.get(key)??new Map<string,Float64Array>();
   const counts=genes.get(o.gene)??new Float64Array(16),q=o.junction??'';
   for(let i=1;i<q.length;i++){const a='ACGT'.indexOf(q[i-1]),b='ACGT'.indexOf(q[i]);if(a>=0&&b>=0)counts[a*4+b]++;}
   genes.set(o.gene,counts);this.transitionCounts.set(key,genes);
  }
  // One fixed-prior E step per row; target-gene sufficient statistics are then
  // subtracted, so neither its boundary nor its SHM signal trains its own null.
  for(const o of observations){if(!o.terminal)continue;const q=this.junction(o),tau=otherBurden(o,new Set(),error);
   const logp=trimWeights.map((w,d)=>{
    let score=Math.log(w);for(let i=0;i<o.terminal!.query.length;i++){
     const p=o.terminal!.start+i,b=o.terminal!.query[i];if(!'ACGT'.includes(b))continue;
     if(p>o.parent.sequence.length-d)score+=Math.log(q['ACGT'.indexOf(b)]);
     else {const rate=1-Math.exp(-tau*hs5fRate(o.parent.sequence,p-1));score+=Math.log(b===o.parent.sequence[p-1]?Math.max(1e-9,(1-error)*(1-rate)):Math.max(1e-9,error/3+(1-error)*rate/3));}
    }return score;
   });const norm=logSum(logp),key=`${o.subjectId}|${o.locus}`,genes=this.trimCounts.get(key)??new Map<string,Float64Array>(),counts=genes.get(o.gene)??new Float64Array(13);
   logp.forEach((v,i)=>counts[i]+=Math.exp(v-norm));genes.set(o.gene,counts);this.trimCounts.set(key,genes);
  }
  const families=new Map<string,EvidenceObservation[]>();
  for(const o of observations){const key=`${o.subjectId}|${o.locus}|${family(o.gene)}`;const list=families.get(key)??[];list.push(o);families.set(key,list);}
  for(const [key,items] of families){
   const parents=[...new Map(items.map(o=>[o.parent.index,o.parent])).values()].sort((a,b)=>b.sequence.length-a.sequence.length||a.gene.localeCompare(b.gene));
   const anchor=parents[0];this.anchors.set(key,anchor.sequence);
   for(const p of parents)this.maps.set(`${key}|${p.index}`,globalCoordinateMap(p.sequence,anchor.sequence).map);
   for(const o of items){const map=this.maps.get(`${key}|${o.parent.index}`)!;const exposure=this.mutationExposure(o,new Set());
    for(let i=0;i<o.positions.length;i++){
     const p=o.positions[i]-1,homology=map[p];if(homology<0||p>=o.parent.sequence.length-12)continue;
     const ref=o.parent.sequence[p];
     for(const alt of 'ACGT'){if(alt===ref)continue;
      const event=`${key}|${homology}|${ref}>${alt}`;
      const genes=this.events.get(event)??new Map();const v=genes.get(o.gene)??{k:0,e:0};
      v.k+=Number(o.query[i]===alt);v.e+=error/3+(1-exposure.naive)*(1-error/3)*(-Math.expm1(-exposure.tau*hs5fRate(o.parent.sequence,p,alt)));
      genes.set(o.gene,v);this.events.set(event,genes);
     }
    }
   }
  }
 }
 /** Training-only calibration for a complete competitive component. This new
  * reader leaves all existing per-gene methods unchanged. */
 componentModel(o:EvidenceObservation,excluded:ReadonlySet<string>):{tau:number;naive:number;trimming:number[];junction:number[]}{
  const key=`${o.subjectId}|${o.locus}`,total={n:0,zeros:0,k:0,w:0,bases:0};
  for(const [gene,v] of this.burdenTotals.get(key)??[])if(!excluded.has(gene))for(const field of ['n','zeros','k','w','bases'] as const)total[field]+=v[field];
  let naive=0,beta=1;
  if(total.n>=20){const errorCount=this.error*total.bases/total.n,mean=Math.max(0,total.k/total.n-errorCount),zeros=Math.min(.999999,total.zeros/total.n/Math.exp(-errorCount)),activeMean=Math.max(.01,mean/Math.max(1e-6,1-zeros)-1);naive=Math.max(0,Math.min(.999,1-mean/activeMean));beta=(total.w/total.n)/activeMean;}
  const pool=(source:Map<string,Map<string,Float64Array>>,fallback:readonly number[])=>{const counts=fallback.map(p=>p*4);for(const [gene,v] of source.get(key)??[])if(!excluded.has(gene))v.forEach((n,i)=>counts[i]+=n);const sum=counts.reduce((a,b)=>a+b,0);return counts.map(n=>n/sum);};
  return {tau:1/beta,naive,trimming:pool(this.trimCounts,trimWeights),junction:pool(this.junctionCounts,[.25,.25,.25,.25])};
 }
 /** Empirical zero-inflated Gamma-Poisson exposure, learned outside the tested gene.
  * The active component has exponential intensity; the atom represents truly
  * unmutated sequences. Sequencing errors remain possible in that component.
  */
 mutationExposure(o:EvidenceObservation,excluded:ReadonlySet<number>):{tau:number;naive:number}{
  let k=0,w=0,n=0;
  for(let i=0;i<o.positions.length;i++){const p=o.positions[i];if(excluded.has(p)||p>o.parent.sequence.length-12)continue;k+=Number(o.query[i]!==o.parent.sequence[p-1]);w+=hs5fRate(o.parent.sequence,p-1);n++;}
  return this.mutationExposureFromCounts(o,k,w,n);
 }
 mutationExposureFromCounts(o:EvidenceObservation,k:number,w:number,n:number):{tau:number;naive:number}{
  const key=`${o.subjectId}|${o.locus}`,ck=key+'|'+o.gene;
  let prior=this.burdenCache.get(ck);
  if(!prior){
   const total={n:0,zeros:0,k:0,w:0,bases:0};
   for(const [gene,v] of this.burdenTotals.get(key)??[])if(gene!==o.gene)for(const field of ['n','zeros','k','w','bases'] as const)total[field]+=v[field];
   if(total.n<20)prior={naive:0,beta:1};
   else {
    const errorCount=this.error*total.bases/total.n,mean=Math.max(0,total.k/total.n-errorCount);
    const zeros=Math.min(.999999,total.zeros/total.n/Math.exp(-errorCount));
    const activeMean=Math.max(.01,mean/Math.max(1e-6,1-zeros)-1);
    prior={naive:Math.max(0,Math.min(.999,1-mean/activeMean)),beta:(total.w/total.n)/activeMean};
   }
   this.burdenCache.set(ck,prior);
  }
  const errorCount=this.error*n,q=w/(prior.beta+w),success=1-q;
  // Poisson sequencing error convolved with geometric SHM count.
  let poisson=Math.exp(-errorCount),active=poisson*success*q**k,activeShape=(k+1)*active;
  for(let errors=1;errors<=k;errors++){poisson*=errorCount/errors;const term=poisson*success*q**(k-errors);active+=term;activeShape+=(k-errors+1)*term;}
  const naiveLikelihood=poisson;
  const posterior=prior.naive*naiveLikelihood/Math.max(1e-300,prior.naive*naiveLikelihood+(1-prior.naive)*active);
  return {tau:activeShape/Math.max(1e-300,active)/(prior.beta+w),naive:posterior};
 }
 private pooled(o:EvidenceObservation,source:Map<string,Map<string,Float64Array>>,fallback:readonly number[],kind:string):number[]{
  const key=`${o.subjectId}|${o.locus}`,ck=key+'|'+o.gene+'|'+kind,cached=this.boundaryCache.get(ck);if(cached)return cached;
  const counts=fallback.map(p=>p*4);for(const [gene,v] of source.get(key)??[])if(gene!==o.gene)v.forEach((n,i)=>counts[i]+=n);
  const total=counts.reduce((a,b)=>a+b,0),result=counts.map(n=>n/total);this.boundaryCache.set(ck,result);return result;
 }
 junction(o:EvidenceObservation):number[]{return this.pooled(o,this.junctionCounts,[.25,.25,.25,.25],'junction');}
 transitions(o:EvidenceObservation):number[]{
  const counts=this.pooled(o,this.transitionCounts,new Array(16).fill(1/16),'transition');
  return counts.map((v,i)=>v/Math.max(1e-12,counts.slice(Math.floor(i/4)*4,Math.floor(i/4)*4+4).reduce((a,b)=>a+b,0)));
 }
 trimming(o:EvidenceObservation):number[]{return this.pooled(o,this.trimCounts,trimWeights,'trim');}
 prior(o:EvidenceObservation,c:Change):RatePrior{
  const fk=`${o.subjectId}|${o.locus}|${family(o.gene)}`,mk=`${fk}|${o.parent.index}`;
  if(!this.maps.has(mk)&&this.anchors.has(fk))this.maps.set(mk,globalCoordinateMap(o.parent.sequence,this.anchors.get(fk)!).map);
  const pos=this.maps.get(mk)?.[c.position-1];
  const event=`${o.subjectId}|${o.locus}|${family(o.gene)}|${pos}|${c.reference}>${c.alternate}`;
  const cacheKey=event+"|"+o.gene;const cached=this.cache.get(cacheKey);if(cached)return cached;
  const values=[...(this.events.get(event)?.entries()??[])].filter(([g,v])=>g!==o.gene&&v.e>=.5).map(([,v])=>v);
  if(values.length<3)return {shape:1,rate:1,genes:values.length};
  // Equal gene votes: millions of bases cannot collapse between-gene dispersion.
  const logs=values.map(v=>Math.log((v.k+.5)/(v.e+.5))).sort((a,b)=>a-b);
  const center=logs[Math.floor(logs.length/2)];
  const deviations=logs.map(x=>(x-center)**2);
  const variance=Math.max(.25,deviations.reduce((a,b)=>a+b,0)/logs.length-values.reduce((a,v)=>a+1/(v.k+.5),0)/values.length);
  const shape=Math.max(.1,Math.min(4,1/Math.expm1(variance)));
  const mean=Math.exp(center+Math.min(variance,4)/2);
  const result={shape,rate:shape/mean,genes:values.length};this.cache.set(cacheKey,result);return result;
 }
}
export function discoveryCost(length:number,changes:number,parents:number):number {
 let cost=Math.log(Math.max(1,parents));for(let i=1;i<=changes;i++)cost+=Math.log((length-i+1)/i)+Math.log(3);return cost;
}
interface Group {tau:number;naive:number;states:number[];n:number}
function maximize(fn:(eta:number)=>number):number {
 // The heavy-tailed prior need not be log-concave: bracket every coarse local
 // maximum, rather than relying on one Newton trajectory or one global bracket.
 const grid=Array.from({length:25},(_,i)=>i-12),scores=grid.map(fn);
 let best=scores.indexOf(Math.max(...scores)),result=grid[best],value=scores[best];
 for(let k=1;k<grid.length-1;k++){
  if(scores[k]<scores[k-1]||scores[k]<scores[k+1])continue;
  let a=grid[k-1],b=grid[k+1];const r=(Math.sqrt(5)-1)/2;
  let x=b-r*(b-a),y=a+r*(b-a),fx=fn(x),fy=fn(y);
  for(let i=0;i<30;i++){if(fx>fy){b=y;y=x;fy=fx;x=b-r*(b-a);fx=fn(x);}else{a=x;x=y;fx=fy;y=a+r*(b-a);fy=fn(y);}}
  const eta=(a+b)/2,score=fn(eta);if(score>value){result=eta;value=score;}
 }
 return result;
}
/** Linked lineage mixture against alternate-specific SHM hazards, with log-rate random effects. */
export function testHaplotype(observations:readonly EvidenceObservation[],changes:readonly Change[],allChanges:readonly Change[],calibration:AlignedRateCalibration,error:number,parents:number,skipBoundaryProfile=false):{gain:number;frequency:number;reason?:string} {
 if(!observations.length||!changes.length)return {gain:-Infinity,frequency:0};
 const parentSequence=observations[0].parent.sequence;
 if(!skipBoundaryProfile&&changes.every(c=>c.position>parentSequence.length-12)){
  // These sequences are also generated by shortening V and changing the freely
  // fitted initial N word. SHM/context approximations must not break that alias.
  if(!terminalChangeIdentifiable(parentSequence.length,changes.map(c=>c.position)))return {gain:-Infinity,frequency:0,reason:"Terminal change is equivalent to V trimming plus the fitted initial junction word"};
  const alternate=[...parentSequence];for(const c of allChanges)alternate[c.position-1]=c.alternate;
  const nullSequence=[...alternate];for(const c of changes)nullSequence[c.position-1]=c.reference;
  const omitted=new Set(allChanges.map(c=>c.position));
  const junction=calibration.junction(observations[0]),transitions=calibration.transitions(observations[0]);
  const rows=observations.filter(o=>o.terminal).map(o=>{
   const tau=otherBurden(o,omitted,error);
   const emit=(seq:string)=>trimWeights.map((_,d)=>{
    let value=0;
    for(let i=0;i<o.terminal!.query.length;i++){
     const p=o.terminal!.start+i,b=o.terminal!.query[i];if(!'ACGT'.includes(b))continue;
     if(p>seq.length-d){
      if(p<=seq.length-d+2)continue; // Initial inserted dinucleotide is a fitted gene-specific nuisance.
      const previous='ACGT'.indexOf(o.terminal!.query[i-1]??'N');
      value+=Math.log(previous<0?junction['ACGT'.indexOf(b)]:transitions[previous*4+'ACGT'.indexOf(b)]);
     }
     else {const mutation=-Math.expm1(-tau*hs5fRate(seq,p-1));value+=Math.log(b===seq[p-1]?(1-error)*(1-mutation):error/3+(1-error)*mutation*hs5fRate(seq,p-1,b)/Math.max(1e-12,hs5fRate(seq,p-1)));}
    }return value;
   });return {nullLog:emit(nullSequence.join('')),alleleLog:emit(alternate.join('')),first:trimWeights.map((_,d)=>{
    const offset=parentSequence.length-d+1-o.terminal!.start,first='ACGT'.indexOf(o.terminal!.query[offset]??'N'),second='ACGT'.indexOf(o.terminal!.query[offset+1]??'N');
    return first<0?-1:second<0?16+first:first*4+second;
   })};
  });
  const fit=profileBoundaryAllele(rows);
  const rearrangementGain=fit.gain-discoveryCost(parentSequence.length,changes.length,parents);
  if(!(rearrangementGain>0))return {...fit,gain:rearrangementGain};
  const shm=testHaplotype(observations,changes,allChanges,calibration,error,parents,true);
  return {...fit,gain:Math.min(rearrangementGain,shm.gain),reason:shm.gain<=0?'Terminal change is explained by the SHM null':undefined};
 }
 const omitted=new Set(allChanges.map(c=>c.position));const grouped=new Map<string,Group>();
 for(const o of observations){const exposure=calibration.mutationExposure(o,omitted),tau=Math.round(exposure.tau*2000)/2000,naive=Math.round(exposure.naive*1000)/1000;
  const states=changes.map(c=>{const b=baseAt(o,c.position);return b===undefined?-1:Number(b===c.alternate);});
  if(states.every(s=>s<0))continue;
  const key=tau+'|'+naive+'|'+states.join('');const g=grouped.get(key)??{tau,naive,states,n:0};g.n++;grouped.set(key,g);
 }
 const groups=[...grouped.values()],n=groups.reduce((s,g)=>s+g.n,0),parent=observations[0].parent.sequence;
 const priors=changes.map(c=>calibration.prior(observations[0],c));
 const rates=changes.map(c=>hs5fRate(parent,c.position-1,c.alternate));
 const probability=(tau:number,j:number,eta:number)=>clamp(error/3+(1-error/3)*(-Math.expm1(-tau*rates[j]*Math.exp(eta))));
 const trim=calibration.trimming(observations[0]),junction=calibration.junction(observations[0]);
 const terminal=changes.some(c=>c.position>parent.length-12);
 // A heavy-tailed log-rate random effect permits multiplicative SHM/selection
 // deviations. It cannot mutate the explicit naive component to invent an intercept.
 const centers=priors.map(p=>Math.log(p.shape/p.rate)-.5*Math.log1p(1/p.shape));
 const scales=priors.map(p=>Math.max(.25,Math.log1p(1/p.shape)/2));
 const units:number[][]=[];
 for(const j of changes.map((_,j)=>j).sort((a,b)=>changes[a].position-changes[b].position)){
  const last=units[units.length-1];
  if(last&&changes[j].position-changes[last[last.length-1]].position<=4)last.push(j);else units.push([j]);
 }
 const clusters=units.filter(unit=>unit.length>1);
 const clusterParameters=new Map(clusters.map((unit,i)=>[unit,changes.length+i]));
 const burstRates=clusters.map(unit=>Math.max(...unit.map(j=>rates[j])));
 for(const unit of clusters){centers.push(Math.max(...unit.map(j=>centers[j])));scales.push(Math.max(...unit.map(j=>scales[j])));}
 const parameterIndices=centers.map((_,i)=>i);
 const logPrior=(j:number,eta:number)=>-2.5*Math.log1p((eta-centers[j])**2/(4*scales[j]));
 const bernoulli=(y:number,p:number)=>y?Math.log(p):Math.log1p(-p);
 const trimStates=new Map<string,{weight:number;retained:boolean[]}>();
 for(let d=0;d<trim.length;d++){
  const retained=changes.map(c=>c.position<=parent.length-d),key=retained.map(Number).join('');
  const previous=trimStates.get(key);if(previous)previous.weight+=trim[d];else trimStates.set(key,{weight:trim[d],retained});
 }
 const states=[...trimStates.values()];
 const jointLog=(g:Group,eta:readonly number[],signal:boolean)=> {
  const p=changes.map((c,j)=>signal?clamp((1-error)*Math.exp(-g.tau*hs5fRate(parent,c.position-1))):probability(g.tau,j,eta[j]));
  const likelihood=(retained?:readonly boolean[])=>units.reduce((total,unit)=>{
   const independent=unit.reduce((sum,j)=>g.states[j]<0?sum:sum+bernoulli(g.states[j],retained&&!retained[j]?junction["ACGT".indexOf(changes[j].alternate)]:p[j]),0);
   const parameter=clusterParameters.get(unit);
   if(signal||parameter===undefined)return total+independent;
   // Nearby substitutions can arise in one repair/context-interaction event.
   // The event is restricted to the mutated component and still vanishes at tau=0.
   const burst=clamp(-Math.expm1(-g.tau*burstRates[parameter-changes.length]*Math.exp(eta[parameter])));
   const linked=unit.reduce((sum,j)=>g.states[j]<0?sum:sum+bernoulli(g.states[j],retained&&!retained[j]?junction['ACGT'.indexOf(changes[j].alternate)]:clamp((1-error)*Math.exp(-g.tau*hs5fRate(parent,changes[j].position-1)))),0);
   return total+logSum([Math.log1p(-burst)+independent,Math.log(burst)+linked]);
  },0);
  const active=terminal?logSum(states.map(s=>Math.log(s.weight)+likelihood(s.retained))):likelihood();
  if(g.naive===0)return active;
  const naiveLog=(retained?:readonly boolean[])=>changes.reduce((sum,c,j)=>g.states[j]<0?sum:sum+bernoulli(g.states[j],retained&&!retained[j]?junction['ACGT'.indexOf(c.alternate)]:signal?1-error:error/3),0);
  const unmutated=terminal?logSum(states.map(s=>Math.log(s.weight)+naiveLog(s.retained))):naiveLog();
  return logSum([Math.log1p(-g.naive)+active,Math.log(g.naive)+unmutated]);
 };
 const optimize=(j:number,weights:Float64Array,etaAll:readonly number[]=parameterIndices.map(j=>centers[j]))=>maximize(value=>{
  if(changes.length===1&&!terminal)return logPrior(j,value)+groups.reduce((sum,g,i)=>{
   if(g.states[j]<0)return sum;
   const p=clamp((1-g.naive)*probability(g.tau,j,value)+g.naive*error/3);
   return sum+g.n*weights[i]*bernoulli(g.states[j],p);
  },0);
  const trial=[...etaAll];trial[j]=value;
  return logPrior(j,value)+groups.reduce((sum,g,i)=>sum+g.n*weights[i]*jointLog(g,trial,false),0);
 });
 const weights=new Float64Array(groups.length).fill(1),nullEta=parameterIndices.map(j=>optimize(j,weights));
 if(changes.length>1)for(let round=0;round<4;round++)for(let j=0;j<parameterIndices.length;j++)nullEta[j]=optimize(j,weights,nullEta);
 const logNull=groups.map(g=>jointLog(g,nullEta,false));
 // A linked allele retains each diagnostic nucleotide, with subsequent SHM allowed.
 const logAllele=groups.map(g=>jointLog(g,nullEta,true));
 const nullObjective=groups.reduce((s,g,i)=>s+g.n*logNull[i],0)+nullEta.reduce((s,e,j)=>s+logPrior(j,e),0);
 let best=-Infinity,bestF=0;
 for(const start of [.1,.9]){
  let f=start,eta=[...nullEta],objective=-Infinity;
  for(let it=0;it<30;it++){
   let count=0,next=0;
   for(let i=0;i<groups.length;i++){const g=groups[i],ln=jointLog(g,eta,false);
    const delta=Math.max(-700,Math.min(700,logAllele[i]-ln));const z=f*Math.exp(delta)/(1-f+f*Math.exp(delta));weights[i]=1-z;count+=g.n*z;
    next+=g.n*(ln+Math.log(1-f+f*Math.exp(delta)));
   }
   next+=eta.reduce((s,e,j)=>s+logPrior(j,e),0);
   f=Math.max(1e-8,Math.min(1-1e-8,count/Math.max(1,n)));for(const j of parameterIndices)eta[j]=optimize(j,weights,eta);
   if(Math.abs(next-objective)<1e-5){objective=next;break;}objective=next;
  }
  // Re-evaluate at returned parameters (not the previous EM iterate).
  objective=eta.reduce((s,e,j)=>s+logPrior(j,e),0);
  for(let i=0;i<groups.length;i++){const g=groups[i],ln=jointLog(g,eta,false);const d=Math.max(-700,Math.min(700,logAllele[i]-ln));objective+=g.n*(ln+Math.log(1-f+f*Math.exp(d)));}
  if(objective>best){best=objective;bestF=f;}
 }
 return {gain:best-nullObjective-.5*Math.log(Math.max(2,n))-discoveryCost(parent.length,changes.length,parents),frequency:bestF};
}
