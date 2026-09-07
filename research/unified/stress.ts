import {writeFileSync} from 'node:fs';
import {PersonalizedGermlineAccumulator} from '../../src/personalized-germline.ts';
import {inferUnifiedGermline} from '../../src/unified-germline.ts';
const frequency=Number(process.argv[2]??0),seedStart=Number(process.argv[3]??201),sites=Number(process.argv[4]??2),regime=process.argv[5]??'mixed',targetLineages=Number(process.argv[6]??3000),out=process.argv[7];
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

const start=performance.now(),dashboard=inferUnifiedGermline(acc.researchSnapshot());const evidence=dashboard.results.flatMap(g=>g.evidence),target=truth.join('').slice(0,-2),positive=evidence.filter(e=>e.logEvidence>0);
const result={frequency,seedStart,sites,regime,targetLineages,carriers,targetProposed:evidence.some(e=>e.sequence===target),targetLogEvidence:evidence.find(e=>e.sequence===target)?.logEvidence??null,positiveExtras:positive.filter(e=>e.sequence!==target).map(e=>({id:e.id,logEvidence:e.logEvidence})),seconds:(performance.now()-start)/1000};
if(out)writeFileSync(out,JSON.stringify({result,dashboard},null,2));console.log(JSON.stringify(result));
