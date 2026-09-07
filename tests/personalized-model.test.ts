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

test('terminal identity requires evidence beyond the free initial junction word',()=>{
 const sequence='ACGT'.repeat(60),make=(position:number)=>{
  const reference=sequence[position-1],alternate=reference==='A'?'C':'A';
  const observations:EvidenceObservation[]=Array.from({length:300},(_,i)=>{
   const query=[...sequence];if(i<120)query[position-1]=alternate;
   return {ordinal:i,subjectId:'animal',locus:'IGH',gene:'IGHV1-1',parent:{index:0,sequence,gene:'IGHV1-1'},positions:Uint16Array.from({length:240},(_,p)=>p+1),query:query.join(''),terminal:{start:229,query:query.slice(228).join('')}};
  });const changes=[{position,reference,alternate}];
  const controls=observations.map((o,i)=>({...o,ordinal:300+i,gene:'IGHV1-2',parent:{index:1,sequence,gene:'IGHV1-2'},query:sequence,terminal:{start:229,query:sequence.slice(228)}}));
  return testHaplotype(observations,changes,changes,new AlignedRateCalibration([...observations,...controls]),.001,2);
 };
 assert.equal(make(240).gain,-Infinity);
 assert.equal(make(239).gain,-Infinity);
 assert.ok(make(237).gain>0,'an identifiable variant with downstream templated sequence must retain power');
});

test('an unmutated component is learned without using the tested gene',()=>{
 const sequence='ACGT'.repeat(60),observations:EvidenceObservation[]=[];
 for(let gene=0;gene<3;gene++)for(let i=0;i<200;i++){
  const query=[...sequence],k=i<100?0:1+i%20;
  for(let p=0;p<k;p++)query[p]=query[p]==='A'?'C':'A';
  observations.push({ordinal:gene*200+i,subjectId:'animal',locus:'IGH',gene:`IGHV1-${gene}`,parent:{index:gene,sequence,gene:`IGHV1-${gene}`},positions:Uint16Array.from({length:240},(_,p)=>p+1),query:query.join('')});
 }
 const calibration=new AlignedRateCalibration(observations);
 assert.ok(calibration.mutationExposure(observations[0],new Set()).naive>.5);
 assert.ok(calibration.mutationExposure(observations[110],new Set()).naive<.01);
});

test('raw V-end evidence stops before an aligned D segment',async()=>{
 const {PersonalizedGermlineAccumulator}=await import('../src/personalized-germline.ts');
 const sequence='ACGT'.repeat(60),acc=new PersonalizedGermlineAccumulator(`>IGHV1-1*01\n${sequence}\n`);
 acc.add({v_call:'IGHV1-1*01',v_sequence_start:'1',v_germline_start:'1',v_sequence_alignment:sequence,v_germline_alignment:sequence,sequence,d_sequence_start:'235'},0,1);
 const observation=[...(acc as unknown as {selected:Map<string,EvidenceObservation>}).selected.values()][0];
 assert.equal(observation.terminal?.query.length,6);
 assert.equal(observation.positions[observation.positions.length-1],234);
});

test('correlated SHM is rejected and high-SHM rare-allele fits stay finite',async()=>{
 const {simulateDiscovery}=await import('../benchmarks/personalized-germline/stress.mts');
 const nullRun=simulateDiscovery(0,104).result;
 assert.equal(nullRun.finite,true);
 assert.deepEqual(nullRun.extras,[],'an adjacent repair event must not become an allele');
 const rare=simulateDiscovery(.01,105,2,'memory').result;
 assert.equal(rare.finite,true,'survival probabilities must stay finite at large hazards');
 assert.equal(rare.recovered,true,'linked support with extra SHM must reach and survive the allele-set fit');
 assert.deepEqual(rare.extras,[]);
});
