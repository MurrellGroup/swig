import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { join, resolve } from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { gzipSync, gunzipSync } from "node:zlib";

function cliInvocation(root,configPath){
  const standalone=process.env.SWIG_CLI_EXECUTABLE;
  return {command:standalone||process.execPath,arguments_:standalone?["run","--config",configPath]:[join(root,"cli/swig-cli.mjs"),"run","--config",configPath]};
}

function runCli(root,configPath){
  const {command,arguments_}=cliInvocation(root,configPath);
  return spawnSync(command,arguments_,{encoding:"utf8",timeout:120_000});
}

test("the CLI pipeline runs annotation through lazy lineage-study export",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-test-"));
  try{
    const output=join(temporary,"out");
    const configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({
      inputs:[{path:join(root,"tests/fixtures/cli-smoke.fasta"),sampleId:"sample-A",subjectId:"donor-A"}],
      annotation:{workers:1},
      pipeline:{lineage:{productiveOnly:false}},
      output:{directory:output,prefix:"smoke"},
    }));
    const result=runCli(root,configPath);
    assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);
    assert.equal(summary.inputRecords,2);
    assert.equal(summary.annotatedRecords,2);
    assert.equal(summary.retainedRecords,1);
    assert.equal(summary.lineages,1);
    const airr=await readFile(join(output,"smoke.lineages.airr.tsv"));
    const manifest=JSON.parse(gunzipSync(await readFile(join(output,"smoke.swig-lineage-study.json.gz"))).toString("utf8"));
    assert.equal(manifest.linkedAirr.size,airr.byteLength);
    assert.equal(manifest.linkedAirr.sha256,createHash("sha256").update(airr).digest("hex"));
    assert.deepEqual(manifest.summaries[0].sampleCounts,[{sampleId:"sample-A",uniqueMembers:1,abundance:2}]);
    const range=manifest.ranges[0];
    const lineageBytes=airr.subarray(range.start,range.end);
    assert.equal(lineageBytes.toString("utf8").trimEnd().split("\n").length,1);
    assert.equal(range.sha256,createHash("sha256").update(lineageBytes).digest("hex"));
  }finally{
    await rm(temporary,{recursive:true,force:true});
  }
});

test("CLI FASTQ preprocessing matches browser QC, wrapped parsing, and per-dataset reservoir sampling",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-fastq-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const sequence=fasta.split(/\r?\n/).slice(1).join("").split(">read2",1)[0];
    const wrapped=(name,quality)=>`@${name}\n${sequence.slice(0,187)}\n${sequence.slice(187)}\n+\n${quality.slice(0,120)}\n${quality.slice(120)}\n`;
    const fastq=[wrapped("good-1","I".repeat(sequence.length)),wrapped("bad","!".repeat(sequence.length)),wrapped("good-2","I".repeat(sequence.length)),wrapped("good-3","I".repeat(sequence.length))].join("");
    const input=join(temporary,"reads.fastq");const output=join(temporary,"out");const configPath=join(temporary,"config.json");
    await writeFile(input,fastq);
    await writeFile(configPath,JSON.stringify({inputs:[{path:input,sampleId:"sample-QC",subjectId:"donor-QC"}],preprocessing:{fastqFilter:{enabled:true,maximumExpectedErrors:0.1,phredOffset:33,trim3Prime:{enabled:false}},subsample:{enabled:true,size:2,seed:7}},annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"qc"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);
    assert.equal(summary.inputRecords,4);
    assert.equal(summary.eligibleRecords,3);
    assert.equal(summary.annotatedRecords,2);
    assert.equal(summary.preprocessing.fastqFilter.stats.recordsEvaluated,4);
    assert.equal(summary.preprocessing.fastqFilter.stats.recordsRejectedExpectedErrors,1);
    assert.deepEqual(summary.preprocessing.subsample,{enabled:true,size:2,seed:7});
    const resolved=JSON.parse(await readFile(join(output,"qc.resolved-config.json"),"utf8"));
    assert.equal(resolved.preprocessing.fastqFilter.maximumExpectedErrors,0.1);
    assert.equal(resolved.preprocessing.subsample.size,2);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI can address separate samples inside one concatenated gzip and preserves member names",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-members-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const records=fasta.trim().split(/(?=>)/).filter(Boolean);
    const first=gzipSync(`${records[0]}\n`),second=gzipSync(`${records[1]}\n`);
    const input=join(temporary,"joined.fasta.gz");const output=join(temporary,"out");const configPath=join(temporary,"config.json");
    await writeFile(input,Buffer.concat([first,second]));
    await writeFile(configPath,JSON.stringify({studyName:"members",inputs:[{path:input,inputName:"day-0.fasta.gz",format:"fasta",gzipRange:{start:0,end:first.length},datasetId:"d0",sampleId:"day-0",subjectId:"donor-A"},{path:input,inputName:"day-30.fasta.gz",format:"fasta",gzipRange:{start:first.length,end:first.length+second.length},datasetId:"d30",sampleId:"day-30",subjectId:"donor-A"}],annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"members"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);assert.equal(summary.inputRecords,2);assert.equal(summary.annotatedRecords,2);assert.equal(summary.retainedRecords,2);assert.equal(summary.lineages,1);
    const manifest=JSON.parse(gunzipSync(await readFile(join(output,"members.swig-lineage-study.json.gz"))).toString("utf8"));
    assert.deepEqual(manifest.analysis.datasets.map((dataset)=>dataset.inputName),["day-0.fasta.gz","day-30.fasta.gz"]);
    const airr=await readFile(join(output,"members.annotated.airr.tsv"),"utf8");
    assert.match(airr,/d0::read1/);assert.match(airr,/d30::read2/);assert.match(airr,/swig_source_sequence_id/);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI AIRR mode can either preserve calls or reannotate sequence, including browser-accepted CSV input",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-airr-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const sequence=fasta.split(/\r?\n/)[1];
    const input=join(temporary,"reads.csv");await writeFile(input,`sequence_id,sequence,v_call,j_call,productive\noriginal,${sequence},FAKEV*01,FAKEJ*01,T\n`);
    for(const airrMode of ["preserve","reannotate"]){
      const output=join(temporary,airrMode);const configPath=join(temporary,`${airrMode}.json`);
      await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"airr",datasetId:"csv",sampleId:"sample",subjectId:"donor"}],annotation:{workers:1,airrMode},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:airrMode,writeLineageStudy:false}}));
      const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
      const table=await readFile(join(output,`${airrMode}.annotated.airr.tsv`),"utf8");
      assert.match(table,/csv::original/);assert.match(table,/swig_source_sequence_id/);
      if(airrMode==="preserve")assert.match(table,/FAKEV\*01/);else assert.doesNotMatch(table,/FAKEV\*01/);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI executes a browser-exported pasted input embedded in the JSON config",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-inline-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const output=join(temporary,"out"),configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({inputs:[{path:"pasted-sequences.txt",inputName:"pasted-sequences.txt",inline:fasta,format:"fasta",sampleId:"paste",subjectId:"donor"}],annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"inline"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);assert.equal(summary.inputRecords,2);assert.equal(summary.annotatedRecords,2);
  }finally{await rm(temporary,{recursive:true,force:true});}
});
