import test from 'node:test';
import assert from 'node:assert/strict';
import { globalCoordinateMap, substitutionCompatible } from '../src/shm-model/alignment.ts';
import { discoveryCost,AlignedRateCalibration,testHaplotype,type EvidenceObservation } from '../src/shm-model/discovery.ts';

test('equal length does not authorize compensating-indel coordinate swaps',()=>{
 const left='G'.repeat(30)+'ACGTAC'+'T'.repeat(30),right='G'.repeat(30)+'CGTACA'+'T'.repeat(30);
 assert.equal(left.length,right.length);assert.equal(substitutionCompatible(left,right),false);
 const point=left.slice(0,20)+'A'+left.slice(21);assert.equal(substitutionCompatible(left,point),true);
 const mapped=globalCoordinateMap(left,right);assert.ok([...mapped.map].some(p=>p<0));
});
test('discrete search cost includes positions, nucleotides and parent identity',()=>{
 assert.ok(Math.abs(discoveryCost(300,2,100)-Math.log(300*299/2*9*100))<1e-10);
});
test('aligned leave-gene-out SHM rejects accumulation but retains an intercept',()=>{
 let seed=17;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/2**32;};
 const seq=Array.from({length:240},()=> 'ACGT'[Math.floor(random()*4)]).join('');
 const position=80,ref=seq[position-1],alt=ref==='A'?'C':'A';
 const observations:EvidenceObservation[]=[];
 for(let gene=0;gene<4;gene++)for(let i=0;i<400;i++){
  const mutations=i%20,tau=mutations/240;const query=[...seq];const chosen=new Set<number>();
  while(chosen.size<mutations){const p=Math.floor(random()*210);if(p===position-1)continue;chosen.add(p);}
  for(const p of chosen)query[p]=query[p]==='A'?'C':'A';
  if(random()<(gene===3?.35:1-Math.exp(-10*tau)))query[position-1]=alt;
  observations.push({ordinal:gene*400+i,subjectId:'animal',locus:'IGH',gene:`IGHV1-${gene}`,parent:{index:gene,sequence:seq,gene:`IGHV1-${gene}`},positions:Uint16Array.from({length:240},(_,p)=>p+1),query:query.join('')});
 }
 const cal=new AlignedRateCalibration(observations),changes=[{position,reference:ref,alternate:alt}];
 const hotspot=testHaplotype(observations.slice(0,400),changes,changes,cal,.001,4);
 const allele=testHaplotype(observations.slice(1200),changes,changes,cal,.001,4);
 assert.ok(hotspot.gain<0,`hotspot gain ${hotspot.gain}`);assert.ok(allele.gain>0,`allele gain ${allele.gain}`);
});

test('boundary likelihood treats downstream junction bases symmetrically',async()=>{
 const {boundaryLogLikelihood,trimWeights}=await import('../src/shm-model/boundary.ts');
 const sequence='A'.repeat(30),junction=[.1,.2,.6,.1];
 const emit=(p:number,b:string)=>Math.log(b==='A'?.999:.001/3);
 const a=boundaryLogLikelihood(sequence,'A',30,emit,trimWeights,junction);
 const g=boundaryLogLikelihood(sequence,'G',30,emit,trimWeights,junction);
 assert.ok(Math.abs(Math.exp(a)-(.2*.999+.8*.1))<1e-12);
 assert.ok(Math.abs(Math.exp(g)-(.2*.001/3+.8*.6))<1e-12);
});

test('a novel allele competes against alternatives derived from different parents',async()=>{
 const {PersonalizedGermlineAccumulator}=await import('../src/personalized-germline.ts');
 const truth='ACGT'.repeat(60),first=[...truth],second=[...truth];first[80]='C';second[120]='C';
 const a=first.join(''),b=second.join('');
 const acc=new PersonalizedGermlineAccumulator(`>IGHV1-1*01\n${a}\n>IGHV1-2*01\n${b}\n`);
 for(let i=0;i<120;i++)acc.add({subject_id:'animal',v_call:i<100?'IGHV1-1*01':'IGHV1-2*01',v_germline_start:'1',v_sequence_alignment:truth,v_germline_alignment:i<100?a:b},i,i+1);
 const result=acc.finish();
 assert.ok(result.pools.flatMap(p=>p.genes).flatMap(g=>g.activeAlleles).some(a=>a.sequence===truth&&!a.known));
});

test('adding a second subject cannot change the first subject inference',async()=>{
 const {PersonalizedGermlineAccumulator}=await import('../src/personalized-germline.ts');
 const ref='ACGT'.repeat(60),other=[...ref],novel=[...ref];other[99]='A';novel[39]='A';novel[79]='A';
 const fasta=`>IGHV1-1*01\n${ref}\n>IGHV1-1*02\n${other.join('')}\n`;
 const single=new PersonalizedGermlineAccumulator(fasta),pooled=new PersonalizedGermlineAccumulator(fasta);
 for(let i=0;i<24;i++){
  const row={subject_id:'first',v_call:'IGHV1-1*01',v_germline_start:'1',v_sequence_alignment:novel.join(''),v_germline_alignment:ref};single.add(row,i,i+1);pooled.add(row,i,i+1);
  pooled.add({subject_id:'second',v_call:'IGHV1-1*02',v_germline_start:'1',v_sequence_alignment:other.join(''),v_germline_alignment:other.join('')},100+i,i+1);
 }
 assert.deepEqual(pooled.finish().pools.find(p=>p.subjectId==='first')?.genes,single.finish().pools[0].genes);
});

test('junction composition excludes the first templated D/J nucleotide',async()=>{
 const {downstreamJunction}=await import('../src/shm-model/boundary.ts');
 assert.equal(downstreamJunction('AAAATGCCC',4,7,9),'TG');
 assert.equal(downstreamJunction('AAAATGCCC',4,0,7),'TG');
 assert.equal(downstreamJunction('AAAATGCCC',4,3,7),'');
});
