/** Broad geometric V deletion prior, with the residual mass assigned to >=12 nt. */
export const TERMINAL_WINDOW=12;
/** A gene-specific initial N dinucleotide is allowed under the rearrangement null. */
export const INITIAL_JUNCTION_BASES=2;
export function terminalChangeIdentifiable(length:number,positions:readonly number[]):boolean {
 return positions.some(p=>p<=length-INITIAL_JUNCTION_BASES);
}
export const trimWeights=Array.from({length:13},(_,d)=>d===12?0.8**12:0.2*0.8**d);
export function logSum(values:readonly number[]):number { const m=Math.max(...values);return m+Math.log(values.reduce((s,v)=>s+Math.exp(v-m),0)); }
export function boundaryLogLikelihood(sequence:string,query:string,start:number,baseLog:(position:number,base:string)=>number, weights:readonly number[]=trimWeights, junction:readonly number[]=[.25,.25,.25,.25]):number {
 return logSum(weights.map((w,d)=>{
  let score=Math.log(w);for(let i=0;i<query.length;i++){const position=start+i;if(!'ACGT'.includes(query[i]))continue;score+=position<=sequence.length-d?baseLog(position,query[i]):Math.log(junction["ACGT".indexOf(query[i])]);}return score;
 }));
}
/** AIRR D/J starts are one-based; the first templated base is excluded. */
export function downstreamJunction(raw:string,start:number,dStart:number,jStart:number):string {
 const end=dStart>0?dStart-1:jStart>0?jStart-1:raw.length;
 return raw.slice(start,Math.min(start+12,end));
}

/** Profile gene-specific deletion probabilities from the complete raw boundary.
 * Both hypotheses estimate the same deletion nuisance; the alternative adds only
 * an allele fraction. A changed-site-only marginal cannot identify this nuisance.
 */
export function profileBoundaryAllele(rows:readonly {nullLog:number[];alleleLog:number[];first?:number[]}[]):{gain:number;frequency:number} {
 if(!rows.length)return {gain:-Infinity,frequency:0};
 const grouped=new Map<string,{nullP:number[];alleleP:number[];n:number;first:number[]}>();
 for(const row of rows){const maximum=Math.max(...row.nullLog,...row.alleleLog);
  const key=row.nullLog.join(',')+'|'+row.alleleLog.join(',')+'|'+(row.first??[]).join(',');const old=grouped.get(key);
  if(old)old.n++;else grouped.set(key,{nullP:row.nullLog.map(x=>Math.exp(x-maximum)),alleleP:row.alleleLog.map(x=>Math.exp(x-maximum)),n:1,first:row.first??new Array(trimWeights.length).fill(-1)});
 }
 const values=[...grouped.values()];
 const startProbability=(q:readonly number[],index:number)=>index<0?1:index<16?q[index]:q.slice((index-16)*4,(index-16)*4+4).reduce((a,b)=>a+b,0);
 const fit=(alternative:boolean,start:number)=>{
  let junction=new Array(16).fill(1/16),weights=[...trimWeights],f=alternative?start:0,previous=-Infinity,value=-Infinity;
  const counts=new Float64Array(weights.length);
  for(let iteration=0;iteration<200;iteration++){
   counts.fill(0);const baseCounts=new Float64Array(16);let alleleCount=0;value=0;
   for(const row of values){
    let normalizer=0,allele=0;
    for(let d=0;d<weights.length;d++){const q=startProbability(junction,row.first[d]);normalizer+=weights[d]*q*((1-f)*row.nullP[d]+f*row.alleleP[d]);allele+=weights[d]*q*f*row.alleleP[d];}
    normalizer=Math.max(1e-300,normalizer);value+=row.n*Math.log(normalizer);alleleCount+=row.n*allele/normalizer;
    for(let d=0;d<weights.length;d++){
     const q=startProbability(junction,row.first[d]),z=row.n*weights[d]*q*((1-f)*row.nullP[d]+f*row.alleleP[d])/normalizer;
     counts[d]+=z;const first=row.first[d];
     if(first>=0&&first<16)baseCounts[first]+=z;
     else if(first>=16)for(let b=0;b<4;b++){const index=(first-16)*4+b;baseCounts[index]+=z*junction[index]/Math.max(1e-300,q);}
    }
   }
   if(Math.abs(value-previous)<1e-7)break;previous=value;
   weights=Array.from(counts,c=>Math.max(1e-12,c/rows.length));const total=weights.reduce((a,b)=>a+b,0);weights=weights.map(w=>w/total);
   const baseTotal=baseCounts.reduce((a,b)=>a+b,0);if(baseTotal>0){junction=Array.from(baseCounts,c=>Math.max(1e-12,c/baseTotal));const sum=junction.reduce((a,b)=>a+b,0);junction=junction.map(q=>q/sum);}
   if(alternative)f=Math.max(1e-9,Math.min(1-1e-9,alleleCount/rows.length));
  }
  return {value,frequency:f};
 };
 const nullFit=fit(false,0),a=fit(true,.1),b=fit(true,.9),best=a.value>b.value?a:b;
 return {gain:best.value-nullFit.value-.5*Math.log(Math.max(2,rows.length)),frequency:best.frequency};
}
