import { trimWeights,logSum } from "./boundary.ts";
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
 private junctionCounts=new Map<string,Map<string,Float64Array>>();
 private trimCounts=new Map<string,Map<string,Float64Array>>();
 private boundaryCache=new Map<string,number[]>();
 private anchors=new Map<string,string>();
 private cache=new Map<string,RatePrior>();
 private maps=new Map<string,Int32Array>();
 private events=new Map<string,Map<string,{k:number,e:number}>>();
 constructor(observations:readonly EvidenceObservation[],error=.001){
  for(const o of observations){
   const key=`${o.subjectId}|${o.locus}`,genes=this.junctionCounts.get(key)??new Map<string,Float64Array>();
   const counts=genes.get(o.gene)??new Float64Array(4);for(const b of o.junction??""){const i="ACGT".indexOf(b);if(i>=0)counts[i]++;}genes.set(o.gene,counts);this.junctionCounts.set(key,genes);
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
   for(const o of items){const map=this.maps.get(`${key}|${o.parent.index}`)!;const tau=otherBurden(o,new Set(),error);
    for(let i=0;i<o.positions.length;i++){
     const p=o.positions[i]-1,homology=map[p];if(homology<0||p>=o.parent.sequence.length-12)continue;
     const ref=o.parent.sequence[p];
     for(const alt of 'ACGT'){if(alt===ref)continue;
      const event=`${key}|${homology}|${ref}>${alt}`;
      const genes=this.events.get(event)??new Map();const v=genes.get(o.gene)??{k:0,e:0};
      v.k+=Number(o.query[i]===alt);v.e+=error/3+tau*hs5fRate(o.parent.sequence,p,alt);
      genes.set(o.gene,v);this.events.set(event,genes);
     }
    }
   }
  }
 }
 private pooled(o:EvidenceObservation,source:Map<string,Map<string,Float64Array>>,fallback:readonly number[],kind:string):number[]{
  const key=`${o.subjectId}|${o.locus}`,ck=key+'|'+o.gene+'|'+kind,cached=this.boundaryCache.get(ck);if(cached)return cached;
  const counts=fallback.map(p=>p*4);for(const [gene,v] of source.get(key)??[])if(gene!==o.gene)v.forEach((n,i)=>counts[i]+=n);
  const total=counts.reduce((a,b)=>a+b,0),result=counts.map(n=>n/total);this.boundaryCache.set(ck,result);return result;
 }
 junction(o:EvidenceObservation):number[]{return this.pooled(o,this.junctionCounts,[.25,.25,.25,.25],'junction');}
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
interface Group {tau:number;states:number[];n:number}
function maximize(fn:(eta:number)=>number):number {
 let a=-12,b=12;const r=(Math.sqrt(5)-1)/2;let x=b-r*(b-a),y=a+r*(b-a),fx=fn(x),fy=fn(y);
 for(let i=0;i<45;i++){if(fx>fy){b=y;y=x;fy=fx;x=b-r*(b-a);fx=fn(x);}else{a=x;x=y;fx=fy;y=a+r*(b-a);fy=fn(y);}}
 return (a+b)/2;
}
/** Linked lineage mixture against alternate-specific SHM hazards, with log-rate random effects. */
export function testHaplotype(observations:readonly EvidenceObservation[],changes:readonly Change[],allChanges:readonly Change[],calibration:AlignedRateCalibration,error:number,parents:number):{gain:number;frequency:number} {
 if(!observations.length||!changes.length)return {gain:-Infinity,frequency:0};
 const omitted=new Set(allChanges.map(c=>c.position));const grouped=new Map<string,Group>();
 for(const o of observations){const tau=Math.round(otherBurden(o,omitted,error)*2000)/2000;
  const states=changes.map(c=>{const b=baseAt(o,c.position);return b===undefined?-1:Number(b===c.alternate);});
  if(states.every(s=>s<0))continue;
  const key=tau+'|'+states.join('');const g=grouped.get(key)??{tau,states,n:0};g.n++;grouped.set(key,g);
 }
 const groups=[...grouped.values()],n=groups.reduce((s,g)=>s+g.n,0),parent=observations[0].parent.sequence;
 const priors=changes.map(c=>calibration.prior(observations[0],c));
 const rates=changes.map(c=>hs5fRate(parent,c.position-1,c.alternate));
 const probability=(tau:number,j:number,eta:number)=>clamp(error/3+(1-error/3)*(-Math.expm1(-tau*rates[j]*Math.exp(eta))));
 const trim=calibration.trimming(observations[0]),junction=calibration.junction(observations[0]);
 const terminal=changes.some(c=>c.position>parent.length-12);
 const logPrior=(j:number,eta:number)=>priors[j].shape*eta-priors[j].rate*Math.exp(eta); // Gamma density transformed to log lambda, including Jacobian.
 const bernoulli=(y:number,p:number)=>y?Math.log(p):Math.log1p(-p);
 const trimStates=new Map<string,{weight:number;retained:boolean[]}>();
 for(let d=0;d<trim.length;d++){
  const retained=changes.map(c=>c.position<=parent.length-d),key=retained.map(Number).join('');
  const previous=trimStates.get(key);if(previous)previous.weight+=trim[d];else trimStates.set(key,{weight:trim[d],retained});
 }
 const states=[...trimStates.values()];
 const retainProbability=changes.map((_,j)=>states.reduce((sum,s)=>sum+(s.retained[j]?s.weight:0),0));
 const jointLog=(g:Group,eta:readonly number[],signal:boolean)=> {
  const p=changes.map((c,j)=>signal?clamp((1-error)*Math.exp(-g.tau*hs5fRate(parent,c.position-1))):probability(g.tau,j,eta[j]));
  const likelihood=(retained?:readonly boolean[])=>changes.reduce((sum,c,j)=>g.states[j]<0?sum:sum+bernoulli(g.states[j],retained&&!retained[j]?junction["ACGT".indexOf(c.alternate)]:p[j]),0);
  return terminal?logSum(states.map(s=>Math.log(s.weight)+likelihood(s.retained))):likelihood();
 };
 const optimize=(j:number,weights:Float64Array,etaAll:readonly number[]=changes.map(()=>0))=>{
  if(terminal&&changes.length>1)return maximize(value=>{
   const trial=[...etaAll];trial[j]=value;
   return logPrior(j,value)+groups.reduce((sum,g,i)=>sum+g.n*weights[i]*jointLog(g,trial,false),0);
  });
  let eta=Math.log(priors[j].shape/priors[j].rate);
  for(let it=0;it<16;it++){
   const lambda=Math.exp(eta);let gradient=priors[j].shape-priors[j].rate*lambda,curvature=-priors[j].rate*lambda;
   for(let i=0;i<groups.length;i++){
    const g=groups[i],y=g.states[j];if(y<0)continue;
    const w=g.n*weights[i],t=g.tau*rates[j]*lambda;
    const retained=terminal?retainProbability[j]:1;
    const p=clamp(retained*probability(g.tau,j,eta)+(1-retained)*junction["ACGT".indexOf(changes[j].alternate)]);
    const numerator=retained*(1-error/3)*Math.exp(-t)*t;
    if(y===0){const d=numerator/(1-p);gradient-=w*d;curvature+=w*(-d*(1-t)-d*d);}
    else {const d=numerator/p;gradient+=w*d;curvature+=w*(d*(1-t)-d*d);}
   }
   const step=Math.max(-2,Math.min(2,gradient/Math.min(-1e-9,curvature)));eta=Math.max(-12,Math.min(12,eta-step));if(Math.abs(step)<1e-5)break;
  }
  return eta;
 };
 const weights=new Float64Array(groups.length).fill(1),nullEta=changes.map((_,j)=>optimize(j,weights));
 if(terminal)for(let round=0;round<4;round++)for(let j=0;j<changes.length;j++)nullEta[j]=optimize(j,weights,nullEta);
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
   f=Math.max(1e-8,Math.min(1-1e-8,count/Math.max(1,n)));eta=changes.map((_,j)=>optimize(j,weights,eta));
   if(Math.abs(next-objective)<1e-5){objective=next;break;}objective=next;
  }
  // Re-evaluate at returned parameters (not the previous EM iterate).
  objective=eta.reduce((s,e,j)=>s+logPrior(j,e),0);
  for(let i=0;i<groups.length;i++){const g=groups[i],ln=jointLog(g,eta,false);const d=Math.max(-700,Math.min(700,logAllele[i]-ln));objective+=g.n*(ln+Math.log(1-f+f*Math.exp(d)));}
  if(objective>best){best=objective;bestF=f;}
 }
 return {gain:best-nullObjective-.5*Math.log(Math.max(2,n))-discoveryCost(parent.length,changes.length,parents),frequency:bestF};
}
