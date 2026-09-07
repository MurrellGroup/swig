import test from 'node:test';import assert from 'node:assert/strict';
import {fit,integrate,splitEvidence,HAZARDS,likelihood} from '../src/shm-model/unified-mixture.ts';
import {transition} from '../src/shm-model/unified-kernel.ts';
test('CTMC and sequencing error preserve probability at zero and high exposure',()=>{for(const t of [0,.001,.1,1]){const p=transition('AACGT',t);for(let a=0;a<4;a++)assert.ok(Math.abs(p.slice(a*4,a*4+4).reduce((s,v)=>s+v,0)-1)<1e-12);}assert.equal(transition('AACGT',0)[0],.999);});
test('concave optimum bound encloses exhaustive simplex grid',()=>{const rows=[{a:Float64Array.of(.3,.1),b:Float64Array.of(.1,.2),count:7},{a:Float64Array.of(.2,.4),b:Float64Array.of(.4,.3),count:13}];const f=fit(rows,[0,1]);let maximum=-Infinity;for(let a=0;a<=100;a++)for(let b=0;b<=100;b++)maximum=Math.max(maximum,likelihood(rows,[a/100,1-a/100],[b/100,1-b/100]));assert.ok(f.upper>=maximum-1e-9);assert.ok(f.ll>=maximum-1e-5);});
function dataset(germline:boolean){const rows=HAZARDS.map(()=>[] as ReturnType<typeof integrate>);for(const t of [0,.03])for(const x of [0,1]){const n=t===0?(x===1?(germline?30:0):(germline?270:300)):(x===1?100:600);const r=integrate(Float64Array.of(Math.log(x===0?.999:.001),Math.log(x===1?.999:.001)),[t],[1]);for(let h=0;h<r.length;h++)rows[h].push({...r[h],count:n});}return rows;}
test('same haplotype can be inherited or somatic in the same likelihood',()=>{const nullRows=dataset(false),altRows=dataset(true);const n=splitEvidence(nullRows,nullRows,[true,false]);const a=splitEvidence(altRows,altRows,[true,false]);assert.ok(n.evidence[0].logEvidence<=0);assert.ok(a.evidence[0].logEvidence>50);});
test('no unmutated observations implies inherited/somatic nonidentifiability',()=>{const rows=dataset(true).map(r=>r.slice(2));const a=splitEvidence(rows,rows,[true,false]);assert.ok(a.evidence[0].logEvidence<=1e-4);});

import {readFileSync} from 'node:fs';
test('the joint box keeps its state, runtime and session separate from the existing method',()=>{
 const ui=readFileSync(new URL('../src/post-analysis.tsx',import.meta.url),'utf8');
 assert.match(ui,/Infer expressed V allele set/);assert.match(ui,/Joint inherited \/ SHM model/);
 const start=ui.indexOf('async function runUnifiedGermlineAnalysis()'),end=ui.indexOf('async function downloadDeduplicated()',start),run=ui.slice(start,end);
 assert.match(run,/setUnifiedGermline\(result\)/);assert.doesNotMatch(run,/setPersonalizedGermline/);
 assert.match(ui,/unifiedGermlineRuntime\.cancel\(\)/);assert.match(ui,/initialSession\.unifiedGermline/);
});
test('hazard integration gives the same total probability in every regime',()=>{
 const sums=HAZARDS.map(()=>({a:[0,0],b:[0,0]}));
 // Every candidate/exposure kernel sums to one over the two observations.
 for(const x of [0,1]){const r=integrate(Float64Array.from([x===0?.9:.1,x===0?.7:.3,x===1?.9:.1,x===1?.7:.3],Math.log),[0,.1],[.3,.7]);
 // integrate row-scales internally, so restore this example's common maximum .9.
 for(let h=0;h<r.length;h++)for(let c=0;c<2;c++){sums[h].a[c]+=.9*r[h].a[c];sums[h].b[c]+=.9*r[h].b[c];}}
 for(const s of sums)for(let c=0;c<2;c++)assert.ok(Math.abs(s.a[c]+s.b[c]-1)<1e-12);
});
import {AlignedRateCalibration,type EvidenceObservation} from '../src/shm-model/discovery.ts';
test('component calibration excludes all tested genes and is row-order invariant',()=>{
 const make=(gene:string,ordinal:number,mutations:number):EvidenceObservation=>{const sequence='ACGT'.repeat(60),q=[...sequence];for(let j=0;j<mutations;j++)q[j]=q[j]==='A'?'C':'A';return {ordinal,subjectId:'s',locus:'IGH',gene,parent:{index:Number(gene.slice(-1)),gene,sequence},positions:Uint16Array.from({length:240},(_,i)=>i+1),query:q.join('')};};
 const external=Array.from({length:40},(_,i)=>make('IGHV1-3',i,i%3===0?0:5));
 const a=[make('IGHV1-1',100,0),make('IGHV1-2',101,30),...external],b=[...external].reverse().concat([make('IGHV1-2',101,0),make('IGHV1-1',100,30)]),excluded=new Set(['IGHV1-1','IGHV1-2']);
 const p=new AlignedRateCalibration(a).componentModel(a[0],excluded),q=new AlignedRateCalibration(b).componentModel(b[0],excluded);
 assert.ok(Math.abs(p.tau-q.tau)<1e-12);assert.ok(Math.abs(p.naive-q.naive)<1e-12);assert.deepEqual(p.junction,q.junction);
});
