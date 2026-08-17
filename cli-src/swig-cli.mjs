import { createReadStream, createWriteStream, readFileSync, writeFileSync } from "node:fs";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { basename, dirname, extname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { fileURLToPath } from "node:url";
import { Worker as NodeWorker } from "node:worker_threads";
import { createGunzip, gzipSync, gunzipSync } from "node:zlib";
import { once } from "node:events";
import { availableParallelism } from "node:os";
import { Readable } from "node:stream";

import { SparseEvidenceAccumulator } from "../src/allele-refinement/evidence.ts";
import { refinedCall } from "../src/allele-refinement/export.ts";
import { toRefinementInputRow } from "../src/allele-refinement/input.ts";
import { fitSparseAlleleModel } from "../src/allele-refinement/model.ts";
import { buildReferenceAlleleGraph } from "../src/allele-refinement/reference-graph.ts";
import { DEFAULT_CLI_CONFIG, normalizeCliConfig } from "../src/pipeline-config.ts";
import { MissingAlleleAccumulator } from "../src/germline-evidence.ts";
import {
  DenoiseAccumulator,
  assignLineages,
  chmmairraDistanceFromReference,
  deduplicate,
  prepareReferenceMsa,
  runChmm,
  threadSequenceToMsa,
} from "../src/post-analysis-core.ts";
import { airrRowToPostAnalysisRecord, denoiseVdjSequence } from "../src/post-analysis-record.ts";
import { compileReferences } from "../src/reference-pack.ts";
import { repertoireRowMatches, validateRepertoireSelection } from "../src/repertoire-selection.ts";
import { ShmAccumulator } from "../src/shm-analysis.ts";
import {
  addFastqQualityFilterStats,
  airrRecords,
  canonicalFastq,
  emptyFastqQualityFilterStats,
  fastaRecords,
  fastqRecords,
  filterFastqRecord,
  seededRandom,
  validateFastqQualityFilter,
} from "../src/sequence-stream.ts";
import { annotateAirrBatch, annotateDoubleDBatch, stableDatasetSeed } from "../src/study-design.ts";

const VERSION="0.30.0";
const CLI_DIRECTORY=dirname(fileURLToPath(import.meta.url));

export function defaultCliAssets(){
  const directory=join(CLI_DIRECTORY,"assets");
  return {wasmPath:join(directory,"swiftig.wasm"),referencePackPath:join(directory,"imgt-reference-pack.json.gz")};
}

function usage(){
  return `swig-cli ${VERSION}\n\n`+
    `Run a complete non-phylogenetic Swig pipeline:\n`+
    `  swig-cli run reads.fastq.gz --out swig-output\n`+
    `  swig-cli run --config swig.config.json\n\n`+
    `Create an editable config:\n`+
    `  swig-cli init swig.config.json\n\n`+
    `Single-input metadata options:\n`+
    `  --sample SAMPLE_ID  --donor SUBJECT_ID  --dataset DATASET_ID\n\n`+
    `Samples with the same subjectId/--donor are treated as the same donor.\n`+
    `Lineage phylogenetics is intentionally not run by swig-cli.`;
}

function argumentValue(args,name){const index=args.indexOf(name);return index>=0?args[index+1]:undefined;}
function hasFlag(args,name){return args.includes(name);}
function positional(args){const value=[];for(let i=0;i<args.length;i+=1){if(args[i].startsWith("--")){if(!["--help","--version"].includes(args[i]))i+=1;continue;}value.push(args[i]);}return value;}
function cleanCell(value){return String(value??"").replace(/[\t\r\n]+/g," ");}

function resolveFrom(base,value){return isAbsolute(value)?value:resolve(base,value);}

function detectFormat(path,requested="auto"){
  if(requested!=="auto")return requested;
  const plain=path.replace(/\.gz$/i,"");
  const extension=extname(plain).toLowerCase();
  if([".fa",".fasta",".fna",".fas"].includes(extension))return "fasta";
  if([".fq",".fastq"].includes(extension))return "fastq";
  if([".tsv",".airr"].includes(extension)||plain.toLowerCase().endsWith(".airr.tsv"))return "airr";
  throw new Error(`Cannot infer the input format for ${path}; set format to fasta, fastq, or airr in the config.`);
}

function inputLines(input){
  const range=input.gzipRange;
  if(typeof input.inline==="string"){
    if(range)throw new Error(`${input.path} cannot combine inline input with gzipRange.`);
    return createInterface({input:Readable.from([input.inline]),crlfDelay:Infinity});
  }
  if(range&&(!Number.isSafeInteger(range.start)||!Number.isSafeInteger(range.end)||range.start<0||range.end<=range.start))throw new Error(`${input.path} has an invalid gzipRange.`);
  const raw=createReadStream(input.path,range?{start:range.start,end:range.end-1}:undefined);
  const stream=/\.gz$/i.test(input.path)||range?raw.pipe(createGunzip()):raw;
  return createInterface({input:stream,crlfDelay:Infinity});
}

async function* sequenceBatches(input,batchRecords,preprocessing,airrMode,datasetIndex,state){
  const format=detectFormat(input.path,input.format);
  const lines=inputLines(input);
  const filter=preprocessing.fastqFilter;
  validateFastqQualityFilter(filter);
  state.fastqFilter=emptyFastqQualityFilterStats(filter.enabled,filter.enabled&&format==="fastq");
  const parsed=async function*(){
    if(format==="airr"){
      for await(const record of airrRecords(lines)){
        state.inputRecords+=1;state.eligibleRecords+=1;
        if(filter.enabled){state.fastqFilter.recordsRetained+=1;state.fastqFilter.recordsPassedThrough+=1;}
        const delimiter=record.header.includes("\t")?"\t":",";
        const header=delimiter==="\t"?record.header:record.header.split(delimiter).join("\t");
        const row=delimiter==="\t"?record.row:record.row.split(delimiter).join("\t");
        yield {ordinal:state.inputRecords-1,text:`${row}\n`,header};
      }
      return;
    }
    if(format==="fasta"){
      for await(const text of fastaRecords(lines)){
        state.inputRecords+=1;state.eligibleRecords+=1;
        if(filter.enabled){state.fastqFilter.recordsRetained+=1;state.fastqFilter.recordsPassedThrough+=1;}
        yield {ordinal:state.inputRecords-1,text};
      }
      return;
    }
    for await(const record of fastqRecords(lines)){
      state.inputRecords+=1;
      const text=filter.enabled?filterFastqRecord(record,filter,state.fastqFilter):canonicalFastq(record);
      if(text!==null){state.eligibleRecords+=1;yield {ordinal:state.inputRecords-1,text};}
    }
  };
  let selected=parsed();
  if(preprocessing.subsample.enabled){
    const size=preprocessing.subsample.size;
    const random=seededRandom(stableDatasetSeed(preprocessing.subsample.seed,datasetIndex));
    const reservoir=[];
    for await(const record of selected){
      if(reservoir.length<size)reservoir.push(record);
      else{const replacement=Math.floor(random()*state.eligibleRecords);if(replacement<size)reservoir[replacement]=record;}
    }
    reservoir.sort((left,right)=>left.ordinal-right.ordinal);
    selected=(async function*(){yield* reservoir;})();
  }
  let batch=[];
  const emit=()=>{
    const header=format==="airr"?batch[0]?.header:undefined;
    const body=batch.map((record)=>record.text).join("");
    const result=format==="airr"&&airrMode==="preserve"
      ? {format,count:batch.length,header,body}
      : {format,count:batch.length,text:header?`${header}\n${body}`:body};
    batch=[];
    return result;
  };
  for await(const record of selected){batch.push(record);state.selectedRecords+=1;if(batch.length>=batchRecords)yield emit();}
  if(batch.length)yield emit();
  if(!state.inputRecords)throw new Error(`No sequence records were found in ${input.path}.`);
  if(!state.eligibleRecords&&state.fastqFilter.applicable)throw new Error(`The FASTQ quality filter rejected all ${state.inputRecords.toLocaleString()} reads in ${input.path}.`);
}

class WasmPool {
  constructor(size,init){this.size=size;this.init=init;this.workers=[];this.pending=new Map();this.nextId=1;this.nextWorker=0;}
  async start(){
    for(let index=0;index<this.size;index+=1){
      const webWorker=Boolean(process.versions.bun&&globalThis.Worker);
      const worker=webWorker?new globalThis.Worker(new URL("./swig-worker.js",import.meta.url)):new NodeWorker(new URL("./swig-worker.mjs",import.meta.url));
      const receive=(message)=>{const pending=this.pending.get(message.id);if(!pending)return;this.pending.delete(message.id);if(message.error)pending.reject(new Error(message.error));else pending.resolve(message.result);};
      const fail=(error)=>{for(const [id,pending] of this.pending){if(pending.worker===worker){this.pending.delete(id);pending.reject(error instanceof Error?error:new Error(error?.message??String(error)));}}};
      if(webWorker){worker.addEventListener("message",(event)=>receive(event.data));worker.addEventListener("error",fail);}
      else{worker.on("message",receive);worker.on("error",fail);}
      this.workers.push(worker);
    }
    await Promise.all(this.workers.map((worker)=>this.request(worker,{type:"init",...this.init})));
  }
  request(worker,message){
    const id=this.nextId++;
    return new Promise((resolvePromise,reject)=>{this.pending.set(id,{resolve:resolvePromise,reject,worker});worker.postMessage({id,...message});});
  }
  run(message){const worker=this.workers[this.nextWorker++%this.workers.length];return this.request(worker,{type:"annotate",...message});}
  async close(){for(const worker of this.workers)await worker.terminate();this.workers=[];}
}

function parseTable(headerText,bodyText){
  const headers=headerText.replace(/\r$/,"").split("\t");
  const rows=[];
  for(const line of bodyText.split(/\r?\n/)){if(!line)continue;const fields=line.split("\t");const values={};headers.forEach((header,index)=>{values[header]=fields[index]??"";});rows.push(values);}
  return {headers,rows};
}

function serializeRows(headers,rows){return `${headers.join("\t")}\n${rows.map((row)=>headers.map((header)=>cleanCell(row[header])).join("\t")).join("\n")}${rows.length?"\n":""}`;}
function addHeaders(headers,names){for(const name of names)if(!headers.includes(name))headers.push(name);}

async function loadReferences(config,base,referencePackPath){
  const packed=JSON.parse(gunzipSync(readFileSync(referencePackPath)).toString("utf8"));
  const species=packed.species.find((entry)=>entry.name===config.references.species);
  if(!species)throw new Error(`The bundled reference pack has no exact species named ${config.references.species}.`);
  const overrides={...config.references.inline};
  for(const segment of ["V","D","J","C"]){const path=config.references.files?.[segment];if(path)overrides[segment]=(await readFile(resolveFrom(base,path))).toString("utf8");}
  const references=compileReferences(species,config.references.scope,overrides);
  if(!references.V.trim()||!references.J.trim())throw new Error("The selected reference composition must contain V and J records.");
  return references;
}

function compactAlleleResult(result){
  return {version:result.version,options:result.options,totalRecords:result.totalRecords,activeRecords:result.activeRecords,runAt:result.runAt,warnings:result.warnings,segments:Object.fromEntries(Object.entries(result.segments).map(([segment,value])=>[segment,{segment:value.segment,models:value.models,modeledRows:value.modeledRows,changedMapRows:value.changedMapRows,skippedRows:value.skippedRows,matrixNonZeros:value.matrixNonZeros,truncatedRows:value.truncatedRows,exactDuplicateLabels:value.exactDuplicateLabels}]))};
}

function fitAlleles(rows,references,options){
  const segments={};
  for(const segment of options.segments){
    const fasta=references[segment];if(!fasta?.trim())continue;
    const graph=buildReferenceAlleleGraph(fasta,segment,options.neighbourRadius);
    const accumulator=new SparseEvidenceAccumulator(rows.length,graph,options);
    rows.forEach((values,ordinal)=>accumulator.add(toRefinementInputRow({ordinal,values},segment)));
    segments[segment]=fitSparseAlleleModel(accumulator.finish(),graph,options,rows.length);
  }
  const warnings=[];
  if(options.model==="active-set")warnings.push("The hurdle model estimates repertoire-active usage, not literal genomic presence.");
  if(options.weighting==="abundance")warnings.push("Abundance weighting allows clonal expansion to influence the fitted mixture.");
  return {version:1,options,totalRecords:rows.length,activeRecords:rows.length,segments,runAt:new Date().toISOString(),warnings};
}

function applyAlleles(rows,result,policy,threshold,headers){
  for(const segment of ["V","D","J"]){
    const prefix=segment.toLowerCase();
    if(!result.segments[segment])continue;
    addHeaders(headers,[`swig_original_${prefix}_call`,`swig_repertoire_${prefix}_call`]);
    rows.forEach((row,ordinal)=>{const call=refinedCall(result,segment,ordinal,policy,threshold);if(!call)return;row[`swig_original_${prefix}_call`]=row[`${prefix}_call`]??"";row[`swig_repertoire_${prefix}_call`]=call;row[`${prefix}_call`]=call;});
  }
}

function runChimera(rows,activeMask,config,scope,headers,references){
  const selectedSegment=config.segment.toUpperCase();
  const msaText=config.uploadedMsa?.trim()||(config.msaSource==="selected"?references[selectedSegment]?.trim():"");
  if(!msaText)throw new Error("CLI chimera filtering requires either pipeline.chimera.uploadedMsa or an aligned selected-segment reference.");
  let msa;
  try{msa=prepareReferenceMsa(msaText);}
  catch(error){
    if(config.msaSource==="selected"&&!config.uploadedMsa?.trim())throw new Error(`The selected ${selectedSegment} references are not already aligned. Export this pipeline from Swig Web to embed its Kalign MSA, or set pipeline.chimera.uploadedMsa explicitly. ${error instanceof Error?error.message:String(error)}`);
    throw error;
  }
  const segment=config.segment.toLowerCase();
  const method=config.model==="auto"?(scope==="TCR"||String(scope).startsWith("TR")?"DB":"BW"):config.model;
  const options={method,priorProbability:config.priorProbability,baseMutationProbability:config.baseMutationProbability,mutationRates:config.mutationRates,mutationSwitchProbability:config.mutationSwitchProbability,detailed:config.detailed};
  const cache=new Map();let evaluated=0,flagged=0,unevaluated=0;
  addHeaders(headers,["swig_chimera_probability","swig_chimera_status"]);
  if(config.detailed)addHeaders(headers,["swig_chimera_starting_reference","swig_chimera_recombinations"]);
  for(let ordinal=0;ordinal<rows.length;ordinal+=1){
    if(!activeMask[ordinal])continue;
    const row=rows[ordinal],call=row[`${segment}_call`]??"",sequence=row[`${segment}_sequence_alignment`]??"",germline=row[`${segment}_germline_alignment`]??"";
    if(!call||!sequence||!germline){row.swig_chimera_status="missing_alignment";unevaluated+=1;if(!config.retainUnevaluated)activeMask[ordinal]=0;continue;}
    const dfr=chmmairraDistanceFromReference(sequence,germline);
    if(dfr<config.minimumDfr){row.swig_chimera_status="low_dfr";unevaluated+=1;if(!config.retainUnevaluated)activeMask[ordinal]=0;continue;}
    const key=`${call}\0${sequence}\0${germline}`;
    try{
      let result=cache.get(key);
      if(result===undefined){result=runChmm(msa,threadSequenceToMsa(sequence,germline,call,msa),sequence,germline,options);cache.set(key,result);}
      row.swig_chimera_probability=String(result.probability);row.swig_chimera_status="evaluated";evaluated+=1;
      if(config.detailed){row.swig_chimera_starting_reference=result.startingReference;row.swig_chimera_recombinations=result.recombinations.map((event)=>`${event.left}->${event.right}@${event.position}`).join(";");}
      if(result.probability>=config.posteriorThreshold){activeMask[ordinal]=0;flagged+=1;}
    }catch(error){row.swig_chimera_status="error";unevaluated+=1;if(!config.retainUnevaluated)activeMask[ordinal]=0;}
  }
  return {evaluated,flagged,unevaluated,threshold:config.posteriorThreshold};
}

function writeChunk(stream,chunk){if(stream.write(chunk))return Promise.resolve();return once(stream,"drain");}

async function writeLineageStudy(path,manifestPath,headers,rows,activeMask,lineages,config,references,shm){
  const stream=createWriteStream(path);const hash=createHash("sha256");let offset=0;
  const append=async(text)=>{const chunk=Buffer.from(text,"utf8");hash.update(chunk);offset+=chunk.byteLength;await writeChunk(stream,chunk);};
  await append(`${headers.join("\t")}\n`);
  const retainedIds=new Set(lineages.summaries.map((summary)=>summary.id));
  const indexed=rows.map((row,ordinal)=>({row,ordinal,lineageId:lineages.assignments[ordinal]??0})).filter((item)=>activeMask[item.ordinal]&&retainedIds.has(item.lineageId)).sort((left,right)=>left.lineageId-right.lineageId||left.ordinal-right.ordinal);
  const ranges=[];let current=0,start=offset,count=0,rangeHash=createHash("sha256");
  const finishRange=()=>{if(current>0)ranges.push({lineageId:current,start,end:offset,records:count,sha256:rangeHash.digest("hex")});};
  for(const item of indexed){
    if(item.lineageId!==current){finishRange();current=item.lineageId;start=offset;count=0;rangeHash=createHash("sha256");}
    const text=`${headers.map((header)=>cleanCell(item.row[header])).join("\t")}\n`;
    rangeHash.update(text,"utf8");await append(text);count+=1;
  }
  finishRange();stream.end();await once(stream,"finish");
  const shmSummaries=shm?.lineages.flatMap((group)=>{const match=/^Lineage\s+(\d+)$/.exec(group.label);return match?[{lineageId:Number(match[1]),mean:group.mean,p95:group.p95??0}]:[]})??[];
  const manifest={schema:1,application:"Swig lineage study",applicationVersion:VERSION,createdAt:new Date().toISOString(),linkedAirr:{name:basename(path),size:offset,records:indexed.length,headers,sha256:hash.digest("hex")},analysis:{inputName:config.studyName,species:config.references.species,scope:config.references.scope,references,datasets:config.inputs.map((input)=>({datasetId:input.datasetId,inputName:input.inputName||basename(input.path),sampleId:input.sampleId,subjectId:input.subjectId,cohort:input.cohort,timepoint:input.timepoint,compartment:input.compartment})),callingProfile:config.annotation.callingProfile,assignerStrategy:config.annotation.assignerStrategy,minimumIdentity:config.annotation.minimumIdentity,strand:config.annotation.strand,lineage:{scope:config.pipeline.lineage.scope,identity:config.pipeline.lineage.identity,resolution:config.pipeline.lineage.resolution,ambiguity:config.pipeline.lineage.ambiguity,productiveOnly:config.pipeline.lineage.productiveOnly,maxCandidateComparisons:config.pipeline.lineage.maxCandidateComparisons}},summaries:lineages.summaries,ranges,shm:shm?{metric:shm.metric,summaries:shmSummaries}:undefined};
  await writeFile(manifestPath,gzipSync(JSON.stringify(manifest)));
  return manifest;
}

async function runPipeline(config,base,assets){
  if(!config.inputs.length)throw new Error("No input datasets were specified.");
  const outputDirectory=resolveFrom(base,config.output.directory);await mkdir(outputDirectory,{recursive:true});
  const references=await loadReferences(config,base,assets.referencePackPath);
  if(config.annotation.airrMode==="preserve"&&config.annotation.doubleD.mode!=="off"&&config.inputs.some((input)=>detectFormat(input.path,input.format)==="airr"))throw new Error("Double-D screening of AIRR input requires annotation.airrMode = \"reannotate\".");
  const needsAnnotation=config.inputs.some((input)=>detectFormat(input.path,input.format)!=="airr"||config.annotation.airrMode==="reannotate");
  const pool=needsAnnotation?new WasmPool(config.annotation.workers,{wasmPath:assets.wasmPath,references,callingProfile:config.annotation.callingProfile,assignerStrategy:config.annotation.assignerStrategy}):null;
  if(pool)await pool.start();
  const rows=[];const headers=[];let annotatedRecords=0;let inputRecords=0;let eligibleRecords=0;
  let fastqFilterStats=emptyFastqQualityFilterStats(config.preprocessing.fastqFilter.enabled,false);
  try{
    for(let datasetIndex=0;datasetIndex<config.inputs.length;datasetIndex+=1){
      const input=config.inputs[datasetIndex];
      const format=detectFormat(input.path,input.format);
      const preserveAirr=format==="airr"&&config.annotation.airrMode==="preserve";
      process.stderr.write(`${preserveAirr?"Loading existing AIRR calls from":"Annotating"} ${input.inputName||basename(input.path)} (${input.sampleId}; donor ${input.subjectId})…\n`);
      const preprocessingState={inputRecords:0,eligibleRecords:0,selectedRecords:0,fastqFilter:emptyFastqQualityFilterStats(false,false)};
      const pending=[];
      const consume=async(item)=>{
        const result=await item.promise;
        const manifest={datasetId:input.datasetId,inputName:input.inputName||basename(input.path),sampleId:input.sampleId,subjectId:input.subjectId,cohort:input.cohort,timepoint:input.timepoint,compartment:input.compartment};
        const annotated=annotateAirrBatch(result.header,result.body,manifest);
        const table=parseTable(annotated.header,annotated.body);
        const doubleDAnnotated=result.doubleDHeader?annotateDoubleDBatch(result.doubleDHeader,result.doubleDBody,manifest):null;
        const doubleD=doubleDAnnotated?parseTable(doubleDAnnotated.header,doubleDAnnotated.body):null;
        const ddById=new Map((doubleD?.rows??[]).map((row)=>[row.sequence_id,row]));
        addHeaders(headers,table.headers);
        for(const row of table.rows){
          const dd=ddById.get(row.sequence_id);if(dd){Object.assign(row,dd);addHeaders(headers,Object.keys(dd));}
          rows.push(row);
        }
        annotatedRecords+=table.rows.length;
      };
      for await(const batch of sequenceBatches(input,config.annotation.batchRecords,config.preprocessing,config.annotation.airrMode,datasetIndex,preprocessingState)){
        const promise=batch.format==="airr"&&config.annotation.airrMode==="preserve"
          ? Promise.resolve({direct:true,header:batch.header,body:batch.body,doubleDHeader:"",doubleDBody:""})
          : pool.run({text:batch.text,count:batch.count,format:batch.format==="fasta"?1:batch.format==="fastq"?2:3,minimumIdentity:config.annotation.minimumIdentity,strand:config.annotation.strand,doubleD:config.annotation.doubleD});
        pending.push({promise,count:batch.count});
        if(pending.length>=Math.max(2,config.annotation.workers*2))await consume(pending.shift());
      }
      while(pending.length)await consume(pending.shift());
      inputRecords+=preprocessingState.inputRecords;eligibleRecords+=preprocessingState.eligibleRecords;
      fastqFilterStats=addFastqQualityFilterStats(fastqFilterStats,preprocessingState.fastqFilter);
    }
  }finally{await pool?.close();}
  rows.forEach((row,ordinal)=>{if(!row.sequence_id)row.sequence_id=`swig_${ordinal+1}`;});
  const prefix=config.output.prefix;
  if(config.output.writeAnnotatedAirr)await writeFile(join(outputDirectory,`${prefix}.annotated.airr.tsv`),serializeRows(headers,rows));

  let alleleResult=null;
  if(config.pipeline.alleleRefinement.enabled){
    process.stderr.write("Fitting repertoire-level allele model…\n");
    alleleResult=fitAlleles(rows,references,config.pipeline.alleleRefinement);
    applyAlleles(rows,alleleResult,config.pipeline.alleleRefinement.reassignmentPolicy,config.pipeline.alleleRefinement.applyMinimumPosterior,headers);
    await writeFile(join(outputDirectory,`${prefix}.allele-models.json`),JSON.stringify(compactAlleleResult(alleleResult),null,2));
  }

  const records=rows.map((values,ordinal)=>airrRowToPostAnalysisRecord({ordinal,values}));
  let collapseResult=null;let activeMask=new Uint8Array(rows.length);activeMask.fill(1);
  if(config.pipeline.collapse.enabled){
    process.stderr.write(`${config.pipeline.collapse.mode==="exact"?"Collapsing exact duplicates":"Denoising reads"}…\n`);
    if(config.pipeline.collapse.mode==="exact")collapseResult=deduplicate(records,config.pipeline.collapse.key,config.pipeline.collapse.unresolvedPolicy,config.pipeline.collapse.scope,config.pipeline.collapse.respectConstantCall);
    else{
      const options={...config.pipeline.collapse.denoise,mode:config.pipeline.collapse.mode,scope:config.pipeline.collapse.scope,respectConstantCall:config.pipeline.collapse.respectConstantCall,unresolvedPolicy:config.pipeline.collapse.unresolvedPolicy};
      const accumulator=new DenoiseAccumulator(records,options);rows.forEach((values,ordinal)=>accumulator.add(ordinal,denoiseVdjSequence({ordinal,values})));collapseResult=accumulator.finish();
    }
    activeMask=Uint8Array.from(collapseResult.counts,(count)=>count>0?1:0);addHeaders(headers,["duplicate_count"]);
    rows.forEach((row,ordinal)=>{if(collapseResult.counts[ordinal])row.duplicate_count=String(collapseResult.counts[ordinal]);});
  }

  let chimeraSummary=null;
  if(config.pipeline.chimera.enabled){process.stderr.write("Filtering candidate chimeras…\n");chimeraSummary=runChimera(rows,activeMask,config.pipeline.chimera,config.references.scope,headers,references);}

  if(config.pipeline.selection.enabled){
    const errors=validateRepertoireSelection(config.pipeline.selection);if(errors.length)throw new Error(errors.join(" "));
    process.stderr.write("Applying explicit repertoire selection…\n");
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal]&&!repertoireRowMatches(row,config.pipeline.selection))activeMask[ordinal]=0;});
  }

  let lineages=null;
  if(config.pipeline.lineage.enabled){
    process.stderr.write("Assigning lineages…\n");
    const doubleDMask=Uint8Array.from(rows,(row)=>row.d2_call?1:0);
    lineages=assignLineages(records,{identity:config.pipeline.lineage.identity,callResolution:config.pipeline.lineage.resolution,ambiguity:config.pipeline.lineage.ambiguity,productiveOnly:config.pipeline.lineage.productiveOnly,requireSameLocus:true,maxCandidateComparisons:config.pipeline.lineage.maxCandidateComparisons,scope:config.pipeline.lineage.scope},collapseResult??undefined,activeMask,doubleDMask);
    addHeaders(headers,["clone_id"]);rows.forEach((row,ordinal)=>{row.clone_id=lineages.assignments[ordinal]>0?String(lineages.assignments[ordinal]):"";});
  }

  let shm=null;
  if(config.pipeline.shm.enabled){
    process.stderr.write("Calculating SHM summaries…\n");
    const accumulator=new ShmAccumulator({metric:config.pipeline.shm.metric,maxSamplesPerLineage:2000});
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])accumulator.add(row,ordinal,lineages?.assignments[ordinal]??0,"All selected");});shm=accumulator.finish();
    await writeFile(join(outputDirectory,`${prefix}.shm.json`),JSON.stringify(shm,null,2));
  }

  let missingAlleles=null;
  if(config.pipeline.missingAlleles.enabled){
    if(!lineages)throw new Error("Missing-allele screening requires lineage assignment.");
    process.stderr.write("Screening for possible missing V alleles…\n");
    const accumulator=new MissingAlleleAccumulator(config.pipeline.missingAlleles);
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])accumulator.add(row,ordinal,lineages.assignments[ordinal]??0);});
    const validator=accumulator.prepareValidation(references.V);
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])validator.add(row,ordinal,lineages.assignments[ordinal]??0);});missingAlleles=validator.finish();
    await writeFile(join(outputDirectory,`${prefix}.missing-v-alleles.json`),JSON.stringify(missingAlleles,null,2));
  }

  addHeaders(headers,["swig_retained"]);rows.forEach((row,ordinal)=>{row.swig_retained=activeMask[ordinal]?"T":"F";});
  const processed=rows.filter((_,ordinal)=>activeMask[ordinal]);
  await writeFile(join(outputDirectory,`${prefix}.processed.airr.tsv`),serializeRows(headers,processed));

  let lineageStudy=null;
  if(config.output.writeLineageStudy&&lineages){
    process.stderr.write("Writing lazy lineage-study bundle…\n");
    lineageStudy=await writeLineageStudy(join(outputDirectory,`${prefix}.lineages.airr.tsv`),join(outputDirectory,`${prefix}.swig-lineage-study.json.gz`),headers,rows,activeMask,lineages,config,references,shm);
  }

  const retained=activeMask.reduce((sum,value)=>sum+(value?1:0),0);
  const summary={application:"swig-cli",version:VERSION,completedAt:new Date().toISOString(),inputRecords,eligibleRecords,annotatedRecords,retainedRecords:retained,lineages:lineages?.lineageCount??0,preprocessing:{subsample:config.preprocessing.subsample,fastqFilter:{options:config.preprocessing.fastqFilter,stats:fastqFilterStats}},collapse:collapseResult?{mode:collapseResult.mode,inputRecords:collapseResult.inputRecords,inputAbundance:collapseResult.inputAbundance,uniqueRecords:collapseResult.uniqueRecords,collapsedRecords:collapseResult.collapsedRecords,warnings:collapseResult.warnings}:null,chimera:chimeraSummary,alleleRefinement:alleleResult?compactAlleleResult(alleleResult):null,shm:shm?{analyzedRecords:shm.analyzedRecords,analyzedAbundance:shm.analyzedAbundance,skippedRecords:shm.skippedRecords,metric:shm.metric}:null,missingAlleles:missingAlleles?{candidates:missingAlleles.candidates.length,warnings:missingAlleles.warnings}:null,lineageStudy:lineageStudy?{manifest:`${prefix}.swig-lineage-study.json.gz`,airr:`${prefix}.lineages.airr.tsv`,records:lineageStudy.linkedAirr.records}:null};
  await writeFile(join(outputDirectory,`${prefix}.summary.json`),JSON.stringify(summary,null,2));
  await writeFile(join(outputDirectory,`${prefix}.resolved-config.json`),JSON.stringify(config,null,2));
  process.stdout.write(`${JSON.stringify(summary,null,2)}\n`);
}

export async function runCli(assets=defaultCliAssets()){
  const args=process.argv.slice(2);const command=args[0]&&!args[0].startsWith("-")?args[0]:"run";const rest=command===args[0]?args.slice(1):args;
  if(hasFlag(args,"--help")||command==="help"){process.stdout.write(`${usage()}\n`);return;}
  if(hasFlag(args,"--version")||command==="version"){process.stdout.write(`${VERSION}\n`);return;}
  if(command==="init"){
    const target=positional(rest)[0]??"swig.config.json";
    const example=normalizeCliConfig({...DEFAULT_CLI_CONFIG,inputs:[{path:"reads.fastq.gz",sampleId:"sample-1",subjectId:"donor-1"}]});
    writeFileSync(target,`${JSON.stringify(example,null,2)}\n`,{flag:"wx"});process.stdout.write(`Created ${target}\n`);return;
  }
  if(command!=="run")throw new Error(`Unknown command ${command}.\n\n${usage()}`);
  const configPath=argumentValue(rest,"--config");let base=process.cwd();let raw={};
  if(configPath){const absolute=resolve(configPath);base=dirname(absolute);raw=JSON.parse(await readFile(absolute,"utf8"));}
  const inputs=positional(rest);
  if(inputs.length){raw.inputs=inputs.map((path,index)=>({path,datasetId:index?undefined:argumentValue(rest,"--dataset"),sampleId:index?undefined:argumentValue(rest,"--sample"),subjectId:index?undefined:argumentValue(rest,"--donor")}));}
  if(argumentValue(rest,"--out"))raw.output={...(raw.output??{}),directory:argumentValue(rest,"--out")};
  const config=normalizeCliConfig(raw);
  if(config.annotation.workers===0)config.annotation.workers=Math.max(1,Math.min(8,availableParallelism()));
  config.inputs=config.inputs.map((input)=>({...input,path:typeof input.inline==="string"?input.path:resolveFrom(base,input.path)}));
  await runPipeline(config,base,assets);
}
