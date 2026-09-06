/** Broad geometric V deletion prior, with the residual mass assigned to >=12 nt. */
export const TERMINAL_WINDOW=12;
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
