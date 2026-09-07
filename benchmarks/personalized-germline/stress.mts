import {writeFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import {PersonalizedGermlineAccumulator} from '../../src/personalized-germline.ts';
/** Deliberately different from HS5F: uniform background SHM, exceptional hotspots,
 * adjacent correlated repair events, sequencing error, and known allele truth. */
export function simulateDiscovery(frequency:number,seedStart:number,sites=2,regime='mixed',targetLineages=5000){
 let seed=seedStart;
const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/2**32;};const bases='ACGT';
const refs=Array.from({length:5},()=>Array.from({length:296},()=>bases[Math.floor(random()*4)]).join(''));
const altered=(b:string)=>bases[(bases.indexOf(b)+1)%4];const truth=[...refs[0]];truth[99]=altered(truth[99]);if(sites!==1)truth[199]=altered(truth[199]);
const acc=new PersonalizedGermlineAccumulator(refs.map((s,i)=>`>IGHV1-${i+1}*01\n${s}\n`).join(''));
let ordinal=0,carriers=0;
for(let gene=0;gene<5;gene++)for(let i=0;i<(gene===0?targetLineages:500);i++){
 const isAllele=gene===0&&random()<frequency;if(isAllele)carriers++;
 const source=isAllele?truth.join(''):refs[gene],query=[...source];const memory=regime==='memory',tau=random()<(memory?.05:.35)?0:-Math.log(Math.max(1e-12,random()))*(memory?.08:.035);
 for(let p=0;p<284;p++)if(random()<1-Math.exp(-tau)||random()<.001)query[p]=bases[(bases.indexOf(query[p])+1+Math.floor(random()*3))%4];
 // A strong, gene-specific hotspot and a linked repair-patch null.
 if(random()<1-Math.exp(-tau*(gene===0?40:5)))query[149]=altered(source[149]);
 if(random()<1-Math.exp(-tau*8)){query[69]=altered(source[69]);query[70]=altered(source[70]);}
 acc.add({subject_id:'animal',v_call:`IGHV1-${gene+1}*01`,v_germline_start:'1',v_sequence_alignment:query.join(''),v_germline_alignment:refs[gene]},ordinal,++ordinal);
}
const t=performance.now(),d=acc.finish(),novel=d.pools.flatMap(p=>p.genes.flatMap(g=>g.activeAlleles.filter(a=>!a.known)));
const result={seed:seedStart,regime:regime,sites:sites,frequency,carriers,recovered:novel.some(a=>a.sequence===truth.join('')),extras:novel.filter(a=>a.sequence!==truth.join('')).map(a=>a.id),finite:d.pools.every(p=>p.genes.every(g=>Number.isFinite(g.relativeLogLikelihood)&&g.activeAlleles.every(a=>Number.isFinite(a.frequency)))),seconds:(performance.now()-t)/1000};
 return {result,dashboard:d};
}
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){
 const run=simulateDiscovery(Number(process.argv[2]??.01),Number(process.argv[3]??101),process.argv[4]==='single'?1:2,process.argv[5]??'mixed',Number(process.argv[7]??5000));
 if(process.argv[6])writeFileSync(process.argv[6],JSON.stringify(run,null,2));
 console.log(JSON.stringify(run.result));
}
