import { performance } from "node:perf_hooks";
import { getHeapStatistics } from "node:v8";
import { MissingAlleleAccumulator } from "../src/germline-evidence.ts";

const records=50_000,reference="ACGT".repeat(55),alternatePosition=50;
function encoded(index:number){return index.toString(4).split("").map(value=>"ACGT"[Number(value)]).join("");}
function row(index:number){const query=[...reference];if(index<12_000)query[alternatePosition-1]=query[alternatePosition-1]==="A"?"G":"A";const token=encoded(index);const length=30+index%4;return {subject_id:"benchmark-subject",v_call:"IGHV1-2*01",j_call:`IGHJ${index%4+1}*01`,cdr3:`TGT${token}${"A".repeat(Math.max(0,length-token.length-6))}TGG`,v_germline_start:"1",v_sequence_alignment:query.join(""),v_germline_alignment:reference};}

const collect=()=>{(globalThis as typeof globalThis&{gc?:()=>void}).gc?.();};const heapUsed=()=>getHeapStatistics().used_heap_size;collect();const baseline=heapUsed(),start=performance.now(),discovery=new MissingAlleleAccumulator();
for(let index=0;index<records;index+=1)discovery.add(row(index),index,index+1);
collect();const afterDiscovery=heapUsed(),validator=discovery.prepareValidation(`>IGHV1-2*01\n${reference}\n`);
for(let index=0;index<records;index+=1)validator.add(row(index),index,index+1);
const dashboard=validator.finish(),end=performance.now();collect();const peak=heapUsed();
console.log(JSON.stringify({records,candidates:dashboard.candidates.length,support:dashboard.candidates[0]?.independentUnits??0,seconds:Number(((end-start)/1000).toFixed(3)),heapDeltaMiB:Number(((peak-baseline)/1048576).toFixed(1)),discoveryHeapDeltaMiB:Number(((afterDiscovery-baseline)/1048576).toFixed(1))}));
