import assert from "node:assert/strict";
import test from "node:test";

import "fake-indexeddb/auto";
import { alignmentText, tableHeader, tableRow, treeNexus } from "../src/export-formats.ts";
import { MissingAlleleAccumulator } from "../src/germline-evidence.ts";
import { DEFAULT_REPERTOIRE_SELECTION, selectRepertoire } from "../src/repertoire-selection.ts";
import { AirrResultStore } from "../src/result-store.ts";
import { decodeSession, encodeSession, packSessionVector, SWIG_SESSION_SCHEMA, unpackSessionVector, type SwigSession } from "../src/session-state.ts";
import { computeShmMetric, ShmAccumulator } from "../src/shm-analysis.ts";

test("tabular, alignment, and tree serializers cover standard interchange formats", () => {
  const fields=["sequence_id","note"];
  assert.equal(tableHeader(fields,"csv"),"sequence_id,note\n");
  assert.equal(tableRow(fields,{sequence_id:"read,1",note:'a "quote"'},"csv"),'"read,1","a ""quote"""\n');
  assert.deepEqual(JSON.parse(tableRow(fields,{sequence_id:"r1",note:"ok"},"jsonl")),{sequence_id:"r1",note:"ok"});
  const records=[{name:"a one",sequence:"AC-G"},{name:"b",sequence:"ACTG"}];
  assert.match(alignmentText(records,"fasta"),/>a one/);
  assert.match(alignmentText(records,"clustal"),/^CLUSTAL W/);
  assert.match(alignmentText(records,"phylip"),/^2 4/);
  assert.match(alignmentText(records,"stockholm"),/^# STOCKHOLM 1.0/);
  assert.match(alignmentText(records,"nexus"),/dimensions ntax=2 nchar=4/);
  assert.match(treeNexus("(a:0.1,b:0.2);"),/\[&R\] \(a:0.1,b:0.2\);/);
});

test("SHM metrics distinguish synonymous and replacement changes and retain multiplicity", () => {
  const row={sequence_id:"r1",v_call:"IGHV1-2*01",j_call:"IGHJ4*02",duplicate_count:"3",v_sequence_alignment:"GCCGATGCT",v_germline_alignment:"GCTGCTGCT",v_sequence_start:"1",sequence_frame:"1",cdr1_start:"1",cdr1_end:"3",fwr1_start:"4",fwr1_end:"9"};
  const metric=computeShmMetric(row,0,7);
  assert.ok(metric);assert.equal(metric.vNtMutations,2);assert.equal(metric.synonymous,1);assert.equal(metric.vAaReplacements,1);assert.equal(metric.cdrNtMutations,1);assert.equal(metric.frameworkNtMutations,1);assert.equal(metric.duplicateCount,3);
  const accumulator=new ShmAccumulator({metric:"vNtRate"});accumulator.add(row,0,7);const dashboard=accumulator.finish();assert.equal(dashboard.analyzedAbundance,3);assert.equal(dashboard.lineages[0].label,"Lineage 7");
});

test("lineage-aware missing-allele warning uses independent units and emits a candidate sequence", () => {
  const reference="ACGT".repeat(55);const accumulator=new MissingAlleleAccumulator();
  for(let index=0;index<24;index+=1){const germline=reference;const query=[...reference];if(index<9)query[49]=query[49]==="A"?"G":"A";accumulator.add({v_call:"IGHV1-2*01",j_call:`IGHJ${index%4+1}*01`,junction:"A".repeat(30+index%4),v_sequence_alignment:query.join(""),v_germline_alignment:germline},index,index+1);}
  const result=accumulator.finish(`>IGHV1-2*01\n${reference}\n`);assert.equal(result.mode,"lineage");assert.equal(result.independentUnits,24);assert.ok(result.candidates.length>=1);assert.equal(result.candidates[0].independentUnits,9);assert.notEqual(result.candidates[0].sequence,reference);
});

test("double-D positive repertoire selection uses sparse evidence and composes call/CDR3 filters", async () => {
  const store=new AirrResultStore();const header=["sequence_id","sequence","locus","v_call","d_call","j_call","productive","cdr3","cdr3_aa","v_identity","j_identity"].join("\t");
  const body=["r0\tACGTACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ4*01\tT\tTGTAAA\tCAK\t0.98\t0.97","r1\tACGTACGT\tIGH\tIGHV3*01\tIGHD2*01\tIGHJ6*01\tT\tTGTCCC\tCCP\t0.96\t0.95"].join("\n")+"\n";
  const ddHeader=["swig_batch_record_index","sequence_id","d_call","d2_call","swig_double_d_score_gain"].join("\t");const ddBody="1\tr1\tIGHD2*01\tIGHD3*01\t12\n";
  await store.appendBatch(header,body,{header:ddHeader,body:ddBody});await store.finalize();
  const selection=await selectRepertoire(store,{...DEFAULT_REPERTOIRE_SELECTION,doubleD:"positive",vCall:"IGHV3",cdr3Nt:"CCC"});assert.equal(selection.inputRecords,2);assert.equal(selection.retainedRecords,1);assert.equal(selection.mask[1],1);assert.equal(selection.mask[0],0);
  const page=await store.page({sequenceId:"",cdr3:"",locus:"",productive:"",vCall:"",dCall:"",jCall:"",cCall:"",isotype:"",minVIdentity:0,minDIdentity:0,minJIdentity:0,minCIdentity:0,minCdr3AaLength:0,maxCdr3AaLength:0,vjInFrame:"",stopCodon:"",completeVdj:"",revComp:"",hasD:false,hasDoubleD:true,hasCdr3:false},0,10);assert.equal(page.rows.length,1);assert.equal(page.rows[0].d2Call,"IGHD3*01");
  await store.clear();
});

test("AIRR session fingerprint is invariant to batch boundaries and changes with content", async()=>{
  const header="sequence_id\tsequence\tv_call\tj_call";const rows=["a\tACGT\tV1\tJ1\n","b\tTGCA\tV2\tJ2\n"];
  const first=new AirrResultStore();await first.appendBatch(header,rows.join(""));await first.finalize();
  const second=new AirrResultStore();await second.appendBatch(header,rows[0]);await second.appendBatch(header,rows[1]);await second.finalize();
  const changed=new AirrResultStore();await changed.appendBatch(header,`${rows[0]}b\tTGCC\tV2\tJ2\n`);await changed.finalize();
  assert.equal(first.fingerprint,second.fingerprint);assert.notEqual(first.fingerprint,changed.fingerprint);
  await first.clear();await second.clear();await changed.clear();
});

test("linked sessions gzip round-trip typed analysis vectors without embedding AIRR rows", async () => {
  const mask=new Uint8Array([1,0,1,1]);
  const session:SwigSession={schema:SWIG_SESSION_SCHEMA,application:"Swig",applicationVersion:"0.13.2",savedAt:new Date(0).toISOString(),linkedAirr:{name:"x.airr.tsv",size:100,lastModified:0,records:4,headers:["sequence_id","sequence"],fingerprint:"abcd"},analysis:{inputName:"x.fasta",species:"Homo sapiens",scope:"IGH",workers:2,callingProfile:"igblast_compatible",minimumIdentity:.6,strand:0,references:{V:">v\nACG\n",D:"",J:">j\nACG\n",C:"",counts:{V:1,D:0,J:1,C:0},annotation:{V:{annotated:1,total:1},J:{annotated:1,total:1}},loci:["IGH"]}},doubleD:[],postAnalysis:{workingStages:[],activeMask:packSessionVector(mask)}};
  const restored=await decodeSession(await encodeSession(session));assert.equal(restored.linkedAirr.records,4);assert.equal(restored.analysis.callingProfile,"igblast_compatible");assert.equal(JSON.stringify(restored).includes("ACGTACGT"),false);assert.deepEqual([...unpackSessionVector(restored.postAnalysis.activeMask!)],[1,0,1,1]);
});
