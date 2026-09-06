import assert from 'node:assert/strict';
import test from 'node:test';
import { gzipSync } from 'node:zlib';
import { assignmentCompletion, overallAssignmentProgress } from '../src/assignment-progress.ts';
import { estimateSequenceCharacters, streamSequenceBatches } from '../src/sequence-stream.ts';

test('equal-sized remaining files retain their share of progress',()=>{
 assert.ok(overallAssignmentProgress(100,100,200,.14)<.6);
 assert.ok(overallAssignmentProgress(100,100,200,0)<.55);
 assert.equal(overallAssignmentProgress(200,0,200,1),.99);
});
test('buffered or dispatched records are not committed work',()=>{
 const state={acknowledged:0,parsed:2000,inputDone:false,inputRecords:2000,eligibleRecords:2000,countHint:10000};
 assert.equal(assignmentCompletion(state),0);
 assert.equal(assignmentCompletion({...state,acknowledged:1000}),.1);
 assert.equal(assignmentCompletion({...state,countHint:null,inputFraction:.2,acknowledged:1000}),.1);
 assert.equal(assignmentCompletion({...state,inputDone:true,acknowledged:1000}),.5);
});
test('subsampling tracks scan plus committed selected reads',()=>{
 const state={acknowledged:0,parsed:0,inputDone:false,inputRecords:10000,eligibleRecords:10000,countHint:10000,subsampleSize:1000};
 assert.equal(assignmentCompletion(state),.1);
 assert.equal(assignmentCompletion({...state,parsed:1000,acknowledged:500}),.55);
 assert.equal(assignmentCompletion({...state,parsed:1000,acknowledged:1000,inputDone:true}),1);
});
test('gzip workload and streaming fraction use expanded characters, not read-ahead bytes',async()=>{
 const input=Array.from({length:1000},(_,i)=>`>r${i}\n${'ACGT'.repeat(100)}\n`).join('');
 const file=new File([gzipSync(input)],'reads.fasta.gz');
 assert.equal(await estimateSequenceCharacters(file),input.length);
 let state:any;
 const generator=streamSequenceBatches({source:file,format:1,batchSize:100,onProgress:p=>{state=p;}});
 const first=await generator.next();assert.equal(first.value.count,100);
 assert.ok(state.charactersRead/state.totalCharacters<.12);
 await generator.return(undefined);
});
