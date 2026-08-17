import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { mkdtemp, readFile, readdir, rm, writeFile } from "node:fs/promises";
import { join, resolve } from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { gzipSync, gunzipSync } from "node:zlib";

function cliInvocation(root,configPath){
  const standalone=process.env.SWIG_CLI_EXECUTABLE;
  return {command:standalone||process.execPath,arguments_:standalone?["run","--config",configPath]:[join(root,"cli/swig-cli.mjs"),"run","--config",configPath]};
}

function runCli(root,configPath,extra=[]){
  const {command,arguments_}=cliInvocation(root,configPath);
  return spawnSync(command,[...arguments_,...extra],{encoding:"utf8",timeout:120_000});
}

function runRawCli(root,args){
  const standalone=process.env.SWIG_CLI_EXECUTABLE;
  return spawnSync(standalone||process.execPath,standalone?args:[join(root,"cli/swig-cli.mjs"),...args],{encoding:"utf8",timeout:120_000});
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
    const headers=airr.toString("utf8").split("\n",1)[0].split("\t");
    for(const field of ["v_support","d_support","j_support","c_support"])assert.ok(headers.includes(field));
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

test("CLI prepares local V/J metadata by default and honors the config opt-out",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-reference-prep-"));
  try{
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const v=human.loci.IGH.V.find((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0));
    const d=human.loci.IGH.D[0];
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);
    assert.ok(v&&d&&j);
    await writeFile(join(temporary,"custom-v.fasta"),`>${v[0]}\n${v[1]}\n`);
    await writeFile(join(temporary,"custom-j.fasta"),`>${j[0]}\n${j[1]}\n`);
    const input=join(temporary,"read.fasta");
    await writeFile(input,`>custom-read\n${v[1]}${d[1]}${j[1]}\n`);
    for(const prepareMetadata of [true,false]){
      const label=prepareMetadata?"prepared":"raw";
      const output=join(temporary,label),configPath=join(temporary,`${label}.json`);
      const references={species:"Homo sapiens",scope:"IGH",files:{V:"custom-v.fasta",J:"custom-j.fasta"}};
      if(!prepareMetadata)references.prepareMetadata=false;
      await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"fasta",sampleId:"sample",subjectId:"donor"}],references,annotation:{workers:1},pipeline:{collapse:{enabled:false},lineage:{productiveOnly:false},shm:{enabled:false}},output:{directory:output,prefix:label}}));
      const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
      const summary=JSON.parse(result.stdout);
      assert.equal(summary.references.metadataPreparation.enabled,prepareMetadata);
      const manifest=JSON.parse(gunzipSync(await readFile(join(output,`${label}.swig-lineage-study.json.gz`))).toString("utf8"));
      if(prepareMetadata){
        assert.equal(summary.references.metadataPreparation.segments.V.annotated,1);
        assert.equal(summary.references.metadataPreparation.segments.J.annotated,1);
        assert.match(manifest.analysis.references.V,/SWIGMETA=/);
        assert.match(manifest.analysis.references.J,/SWIGMETA=/);
        assert.match(result.stderr,/Preparing custom V reference metadata/);
      }else{
        assert.deepEqual(summary.references.metadataPreparation.segments,{});
        assert.doesNotMatch(manifest.analysis.references.V,/SWIGMETA=/);
        assert.doesNotMatch(manifest.analysis.references.J,/SWIGMETA=/);
        assert.doesNotMatch(result.stderr,/Preparing custom V reference metadata/);
      }
      const resolved=JSON.parse(await readFile(join(output,`${label}.resolved-config.json`),"utf8"));
      assert.equal(resolved.references.prepareMetadata,prepareMetadata);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("pipeline CLI requires an explicit output directory and command-line workers override config",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-output-"));
  try{
    const input=join(root,"tests/fixtures/cli-smoke.fasta");
    const missingPath=join(temporary,"missing-output.json");
    await writeFile(missingPath,JSON.stringify({inputs:[{path:input}],annotation:{workers:1}}));
    const missing=runCli(root,missingPath);
    assert.notEqual(missing.status,0);
    assert.match(missing.stderr,/requires an explicit output directory/);

    const output=join(temporary,"out"),configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({inputs:[{path:input}],annotation:{workers:1},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:"workers"}}));
    const result=runCli(root,configPath,["--workers","2"]);
    assert.equal(result.status,0,result.stderr);
    const resolved=JSON.parse(await readFile(join(output,"workers.resolved-config.json"),"utf8"));
    assert.equal(resolved.annotation.workers,2);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj streams assignment-only, IgBLAST-data, and Swig-annotation modes without downstream state",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-vdj-"));
  try{
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const v=human.loci.IGH.V.find((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0));
    const d=human.loci.IGH.D[0];
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);
    assert.ok(v&&d&&j);
    const vPath=join(temporary,"V.fasta"),dPath=join(temporary,"D.fasta"),jPath=join(temporary,"J.fasta"),queryPath=join(temporary,"query.fasta");
    await writeFile(vPath,`>IMGT|${v[0]}|Homo sapiens\n${v[1].slice(0,30)}...${v[1].slice(30)}\n`);
    await writeFile(dPath,`>${d[0]}\n${d[1]}\n`);
    await writeFile(jPath,`>${j[0]}\n${j[1]}\n`);
    await writeFile(queryPath,`>read-1\n${v[1]}AACCGG${d[1]}TTG${j[1]}\n`);
    const bounds=v[2].slice(2,12);
    const internalPath=join(temporary,"human.ndm.imgt"),auxPath=join(temporary,"human_gl.aux");
    await writeFile(internalPath,`# FWR/CDR positions are 1-based; frame is 0-based\n${v[0]} ${bounds.map((value,index)=>index%2===0?value+1:value).join(" ")} VH ${v[2][0]}\n`);
    await writeFile(auxPath,`# name, frame, chain, CDR3 stop, trailing non-coding bases; all but intervals are 0-based\n${j[0]} ${j[2][0]} JH ${j[2][1]} 1\n`);
    const common=["--vdj","-query",queryPath,"-germline_db_V",vPath,"-germline_db_D",dPath,"-germline_db_J",jPath,"-outfmt","19","-num_threads","3","--workers","1","--batch-records","1"];
    const readRow=async(path)=>{
      const [header,line]=String(await readFile(path,"utf8")).trimEnd().split("\n");
      const names=header.split("\t"),values=line.split("\t");return Object.fromEntries(names.map((name,index)=>[name,values[index]??""]));
    };

    const plainPath=join(temporary,"plain.airr.tsv");
    const plain=runRawCli(root,[...common,"-out",plainPath]);
    assert.equal(plain.status,0,plain.stderr);assert.match(plain.stderr,/assignments-only; 1 worker/);
    const plainRow=await readRow(plainPath);
    assert.equal(plainRow.v_call,v[0]);assert.equal(plainRow.j_call,j[0]);assert.equal(plainRow.cdr1,"");assert.equal(plainRow.cdr3,"");assert.equal(plainRow.region_definition,"");
    assert.ok(plainRow.v_support&&Number.isFinite(Number(plainRow.v_support))&&Number(plainRow.v_support)>=0);
    assert.ok(plainRow.d_support&&Number.isFinite(Number(plainRow.d_support))&&Number(plainRow.d_support)>=0);
    assert.ok(plainRow.j_support&&Number.isFinite(Number(plainRow.j_support))&&Number(plainRow.j_support)>=0);
    assert.equal(plainRow.c_support,"");

    const igblastPath=join(temporary,"igblast-data.airr.tsv");
    const igblast=runRawCli(root,[...common,"-custom_internal_data",internalPath,"-auxiliary_data",auxPath,"-out",igblastPath]);
    assert.equal(igblast.status,0,igblast.stderr);assert.match(igblast.stderr,/IgBLAST internal data/);assert.match(igblast.stderr,/IgBLAST auxiliary data/);
    const igblastRow=await readRow(igblastPath);
    assert.ok(igblastRow.cdr1);assert.ok(igblastRow.cdr3);assert.equal(igblastRow.region_definition,"IMGT");
    assert.equal(Number(igblastRow.j_sequence_end)-Number(igblastRow.fwr4_end),1);

    const swigPath=join(temporary,"swig.airr.tsv");
    const swig=runRawCli(root,[...common,"--swigannots","-organism","human","-ig_seqtype","Ig","-out",swigPath]);
    assert.equal(swig.status,0,swig.stderr);assert.match(swig.stderr,/Swig metadata preparation/);
    const swigRow=await readRow(swigPath);assert.ok(swigRow.cdr1);assert.ok(swigRow.cdr3);
    assert.deepEqual((await readdir(temporary)).filter((name)=>/summary|resolved-config|processed|lineage/i.test(name)),[]);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj refuses to accumulate output when no destination is supplied",()=>{
  const root=resolve(import.meta.dirname,"..");
  const result=runRawCli(root,["--vdj","-germline_db_V","V.fasta","-germline_db_J","J.fasta"]);
  assert.notEqual(result.status,0);assert.match(result.stderr,/requires -out/);
});

test("CLI exposes the embedded IMGT data notice",()=>{
  const root=resolve(import.meta.dirname,"..");
  const result=runRawCli(root,["notices"]);
  assert.equal(result.status,0,result.stderr);
  assert.match(result.stdout,/IMGT\/GENE-DB reference data/);
  assert.match(result.stdout,/CC BY 4\.0/);
  assert.match(result.stdout,/Swig modifies the source data/);
  assert.match(result.stdout,/gki010/);
});
