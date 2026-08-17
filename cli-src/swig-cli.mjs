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
  applyIgblastAuxiliaryData,
  applyIgblastDFrameData,
  applyIgblastInternalData,
  prepareIgblastStyleGermlineFasta,
  preprocessGermlineFastaAcrossTiers,
} from "../src/germline-preprocess.ts";
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
import { compileReferences, germlineTemplateTiers, lociForScope } from "../src/reference-pack.ts";
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

const VERSION="0.35.0";
const CLI_DIRECTORY=dirname(fileURLToPath(import.meta.url));

export function defaultCliAssets(){
  const directory=join(CLI_DIRECTORY,"assets");
  return {wasmPath:join(directory,"swiftig.wasm"),referencePackPath:join(directory,"imgt-reference-pack.json.gz")};
}

function usage(){
  return `swig-cli ${VERSION}\n\n`+
    `Run a complete non-phylogenetic Swig pipeline:\n`+
    `  swig-cli run reads.fastq.gz --out swig-output\n`+
    `  swig-cli run --config swig.config.json [--out DIRECTORY] [--workers N]\n\n`+
    `Run only streaming V(D)J assignment (AIRR outfmt 19):\n`+
    `  swig-cli --vdj -query reads.fasta -germline_db_V V.fasta -germline_db_D D.fasta \\\n`+
    `    -germline_db_J J.fasta -out calls.airr.tsv\n\n`+
    `Display bundled-data attribution and license:\n`+
    `  swig-cli notices\n\n`+
    `Create an editable config:\n`+
    `  swig-cli init swig.config.json\n\n`+
    `Single-input metadata options:\n`+
    `  --sample SAMPLE_ID  --donor SUBJECT_ID  --dataset DATASET_ID\n\n`+
    `Samples with the same subjectId/--donor are treated as the same donor.\n`+
    `Lineage phylogenetics is intentionally not run by swig-cli.`;
}

function vdjUsage(){
  return `swig-cli ${VERSION} --vdj\n\n`+
    `Low-overhead, streaming SwiftIG V(D)J assignment with IgBLAST-style option names.\n\n`+
    `Required:\n`+
    `  -germline_db_V FASTA  -germline_db_J FASTA  -out AIRR_TSV\n\n`+
    `Input and optional references:\n`+
    `  -query FASTA            Query FASTA or '-' for stdin (default '-')\n`+
    `  -germline_db_D FASTA    D germline FASTA\n`+
    `  -c_region_db FASTA      Constant-region FASTA\n\n`+
    `Annotation modes (default: assignments only; CDR/FWR fields remain empty):\n`+
    `  -custom_internal_data FILE  IgBLAST V .ndm.imgt data (1-based inclusive intervals)\n`+
    `  -auxiliary_data FILE        IgBLAST J .aux data (0-based frame/CDR3 stop)\n`+
    `  -d_frame_data FILE          IgBLAST D frame-one starts\n`+
    `  --swigannots                Infer/validate metadata as in Swig Web\n\n`+
    `Execution:\n`+
    `  -num_threads N          Worker count; --workers N overrides it\n`+
    `  --workers N             Exact workers, or 0 for automatic\n`+
    `  --batch-records N       Records per bounded WASM batch (default 2000)\n`+
    `  -strand both|plus|minus -outfmt 19 -organism NAME -ig_seqtype Ig|TCR\n\n`+
    `The germline options take FASTA files, not makeblastdb binary prefixes. Output is SwiftIG AIRR,\n`+
    `not IgBLAST pairwise/tabular formatting. The output path is mandatory and is written incrementally.`;
}

function thirdPartyNotices(){
  return `Bundled IMGT/GENE-DB reference data\n\n`+
    `Source: IMGT/GENE-DB release 202632-7, retrieved 2026-08-08.\n`+
    `Copyright © 1995-2026 IMGT®, the international ImMunoGeneTics information system®.\n`+
    `Attribution: IMGT®, the international ImMunoGeneTics information system®, https://www.imgt.org/, Institute of Human Genetics, Université de Montpellier and CNRS.\n`+
    `License: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/\n`+
    `Terms: https://www.imgt.org/about/termsofuse.php\n`+
    `Citation: Giudicelli V, Chaume D, Lefranc M-P. Nucleic Acids Research. 2005;33:D593-D597. https://doi.org/10.1093/nar/gki010\n\n`+
    `Swig modifies the source data by selecting and reorganizing IG/TR V/D/J/C records, normalizing and ungapping nucleotide sequences, deriving compact coordinate metadata, selecting one source sequence per allele identifier, and joining selected coding IGH/TR constant exons. Membrane-only and untranslated constant exons are omitted. IMGT, Université de Montpellier, and CNRS do not endorse Swig or warrant the modified pack or its use.`;
}

function argumentValue(args,name){const index=args.indexOf(name);if(index>=0)return args[index+1];const inline=args.find((value)=>value.startsWith(`${name}=`));return inline?.slice(name.length+1);}
function hasFlag(args,name){return args.includes(name);}
function positional(args){const value=[];for(let i=0;i<args.length;i+=1){if(args[i].startsWith("--")){if(!args[i].includes("=")&&!["--help","--version"].includes(args[i]))i+=1;continue;}value.push(args[i]);}return value;}
function cleanCell(value){return String(value??"").replace(/[\t\r\n]+/g," ");}

function parseIntegerOption(value,label,{minimum=0,allowZero=true}={}){
  if(value===undefined||!/^\d+$/.test(value))throw new Error(`${label} requires an integer value.`);
  const parsed=Number(value);
  if(!Number.isSafeInteger(parsed)||parsed<minimum||(!allowZero&&parsed===0))throw new Error(`${label} has an invalid value: ${value}.`);
  return parsed;
}

function parseFiniteOption(value,label){
  const parsed=Number(value);if(value===undefined||!Number.isFinite(parsed))throw new Error(`${label} requires a numeric value.`);return parsed;
}

function parseVdjArguments(rawArgs){
  const aliases=new Map([
    ["--query","-query"],["--output","-out"],["--out","-out"],
    ["--germline-db-v","-germline_db_V"],["--germline-db-d","-germline_db_D"],
    ["--germline-db-j","-germline_db_J"],["--c-region-db","-c_region_db"],
  ]);
  const valued=new Set([
    "-query","-out","-germline_db_V","-germline_db_D","-germline_db_J","-c_region_db",
    "-custom_internal_data","-auxiliary_data","-d_frame_data","-organism","-domain_system",
    "-ig_seqtype","-strand","-outfmt","-num_threads","--workers","--batch-records",
    "--minimum-identity","--assigner","--calling-profile","-min_D_match","-min_J_length",
    "-num_alignments_D","-num_alignments_J","-D_penalty","-J_penalty",
  ]);
  const flags=new Set(["--swigannots","-show_translation"]);
  const options={};
  for(let index=0;index<rawArgs.length;index+=1){
    let token=rawArgs[index];
    if(token==="--vdj")continue;
    if(["-h","-help","--help"].includes(token)){options.help=true;continue;}
    if(["-version","--version"].includes(token)){options.version=true;continue;}
    const equals=token.indexOf("=");let inline;
    if(equals>0){inline=token.slice(equals+1);token=token.slice(0,equals);}
    token=aliases.get(token)??token;
    if(flags.has(token)){
      if(inline!==undefined)throw new Error(`${token} does not take a value.`);
      options[token]=true;continue;
    }
    if(!valued.has(token))throw new Error(`Unsupported --vdj option ${token}.\n\n${vdjUsage()}`);
    const value=inline!==undefined?inline:rawArgs[++index];
    if(value===undefined)throw new Error(`${token} requires a value.`);
    options[token]=value;
  }
  return options;
}

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
  if(input.path==="-"&&range)throw new Error("Standard input cannot be combined with gzipRange.");
  if(range&&(!Number.isSafeInteger(range.start)||!Number.isSafeInteger(range.end)||range.start<0||range.end<=range.start))throw new Error(`${input.path} has an invalid gzipRange.`);
  const raw=input.path==="-"?process.stdin:createReadStream(input.path,range?{start:range.start,end:range.end-1}:undefined);
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

function trimAuxiliaryFwr4(headerText,bodyText,endOffsets,jLengths){
  if(!Object.keys(endOffsets).length)return bodyText;
  const headers=headerText.split("\t");const at=new Map(headers.map((header,index)=>[header,index]));
  const required=["j_call","j_sequence_start","j_germline_start","j_sequence_alignment","j_germline_alignment","fwr4","fwr4_aa","fwr4_end"];
  if(required.some((header)=>!at.has(header)))return bodyText;
  const lines=[];
  for(const line of bodyText.split(/\r?\n/)){
    if(!line)continue;
    const values=line.split("\t");const call=(values[at.get("j_call")]??"").split(",",1)[0];
    const offset=endOffsets[call],referenceLength=jLengths.get(call);
    if(offset===undefined||!referenceLength||offset===0){lines.push(line);continue;}
    const queryStart=Number(values[at.get("j_sequence_start")]);
    const referenceStart=Number(values[at.get("j_germline_start")]);
    const queryAlignment=values[at.get("j_sequence_alignment")]??"";
    const germlineAlignment=values[at.get("j_germline_alignment")]??"";
    const currentEnd=Number(values[at.get("fwr4_end")]);
    const targetReferenceEnd=referenceLength-offset;
    if(!Number.isFinite(queryStart)||!Number.isFinite(referenceStart)||!Number.isFinite(currentEnd)||targetReferenceEnd<referenceStart-1){lines.push(line);continue;}
    let queryPosition=queryStart-1,referencePosition=referenceStart-1;
    const columns=Math.min(queryAlignment.length,germlineAlignment.length);
    for(let column=0;column<columns;column+=1){
      if(referencePosition>=targetReferenceEnd)break;
      if(queryAlignment[column]!=="-")queryPosition+=1;
      if(germlineAlignment[column]!=="-")referencePosition+=1;
    }
    if(referencePosition<targetReferenceEnd||queryPosition>=currentEnd){lines.push(line);continue;}
    const trim=currentEnd-queryPosition;
    const nucleotide=values[at.get("fwr4")]??"";
    const kept=Math.max(0,nucleotide.length-trim);
    values[at.get("fwr4")]=nucleotide.slice(0,kept);
    values[at.get("fwr4_aa")]=(values[at.get("fwr4_aa")]??"").slice(0,Math.floor(kept/3));
    values[at.get("fwr4_end")]=String(queryPosition);
    lines.push(values.join("\t"));
  }
  return lines.length?`${lines.join("\n")}\n`:"";
}

function serializeRowBody(headers,rows){return `${rows.map((row)=>headers.map((header)=>cleanCell(row[header])).join("\t")).join("\n")}${rows.length?"\n":""}`;}
function addHeaders(headers,names){for(const name of names)if(!headers.includes(name))headers.push(name);}

async function loadReferences(config,base,referencePackPath){
  const packed=JSON.parse(gunzipSync(readFileSync(referencePackPath)).toString("utf8"));
  const species=packed.species.find((entry)=>entry.name===config.references.species);
  if(!species)throw new Error(`The bundled reference pack has no exact species named ${config.references.species}.`);
  const overrides={...config.references.inline};
  const preparation={};
  const allowedLoci=lociForScope(species,config.references.scope);
  if(!allowedLoci.length)throw new Error(`The bundled reference pack has no ${config.references.scope} loci for ${config.references.species}.`);
  for(const segment of ["V","D","J","C"]){
    const path=config.references.files?.[segment];if(!path)continue;
    const raw=(await readFile(resolveFrom(base,path))).toString("utf8");
    if(config.references.prepareMetadata){
      process.stderr.write(`Preparing custom ${segment} reference metadata from ${basename(path)}…\n`);
      const report=preprocessGermlineFastaAcrossTiers(raw,segment,germlineTemplateTiers(packed,species,config.references.scope,segment),allowedLoci);
      overrides[segment]=report.fasta;
      const {fasta,...summary}=report;preparation[segment]={file:path,...summary};
      const detail=segment==="V"?`${report.annotated.toLocaleString()} with validated FWR/CDR metadata`:segment==="J"?`${report.annotated.toLocaleString()} with validated frame/CDR3-anchor metadata`:"validated and normalized";
      process.stderr.write(`Prepared ${report.count.toLocaleString()} ${segment} record${report.count===1?"":"s"}; ${detail}.\n`);
      for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
    }else overrides[segment]=raw;
  }
  const references=compileReferences(species,config.references.scope,overrides);
  if(!references.V.trim()||!references.J.trim())throw new Error("The selected reference composition must contain V and J records.");
  return {references,preparation};
}

async function readMaybeCompressedText(path,label){
  if(path==="-")throw new Error(`${label} must be a FASTA/data file; standard input is reserved for -query.`);
  let bytes;
  try{bytes=await readFile(resolve(path));}
  catch(error){
    if(error?.code==="ENOENT")throw new Error(`${label} was not found at ${path}. The IgBLAST-style germline options require source FASTA files, not makeblastdb prefixes.`);
    throw error;
  }
  if(bytes[0]===0x1f&&bytes[1]===0x8b)bytes=gunzipSync(bytes);
  return bytes.toString("utf8");
}

async function readVdjFasta(path,label){
  const text=await readMaybeCompressedText(path,label);
  if(!text.trimStart().startsWith(">"))throw new Error(`${label} must point to a nucleotide FASTA file. Binary BLAST database files/prefixes are not decoded by swig-cli.`);
  return text;
}

function vdjScope(options){
  const value=String(options["-ig_seqtype"]??"Ig").toUpperCase();
  if(value==="IG")return "BCR";
  if(value==="TCR")return "TCR";
  throw new Error("-ig_seqtype must be Ig or TCR.");
}

function vdjSpecies(pack,requested){
  const aliases={
    human:"Homo sapiens",
    mouse:"Mus musculus_C57BL/6",
    rat:"Rattus norvegicus_BN; Sprague-Dawley",
    rabbit:"Oryctolagus cuniculus",
    rhesus_monkey:"Macaca mulatta_AG07107",
  };
  const value=String(requested??"human");
  const target=aliases[value.toLowerCase()]??value;
  const species=pack.species.find((entry)=>entry.name.toLowerCase()===target.toLowerCase());
  if(!species)throw new Error(`The embedded Swig reference pack has no species matching -organism ${value}. Use one of human, mouse, rat, rabbit, rhesus_monkey, or an exact pack species name.`);
  return species;
}

function fastaLengths(fasta){
  const result=new Map();let name="",sequence="";
  const finish=()=>{if(name)result.set(name,sequence.replace(/\s/g,"").length);};
  for(const line of fasta.split(/\r?\n/)){
    if(line.startsWith(">")){finish();name=line.slice(1).trim().split(/\s+/,1)[0];sequence="";}
    else if(name)sequence+=line.trim();
  }
  finish();return result;
}

function unmatchedMessage(kind,application){
  const preview=application.unmatched.slice(0,8).join(", ");
  return `${kind} did not match ${application.unmatched.length.toLocaleString()} selected germline identifier${application.unmatched.length===1?"":"s"}${preview?`: ${preview}${application.unmatched.length>8?", …":""}`:""}.`;
}

async function prepareVdjReferences(options,assets){
  const paths={V:options["-germline_db_V"],D:options["-germline_db_D"],J:options["-germline_db_J"],C:options["-c_region_db"]};
  if(!paths.V||!paths.J)throw new Error("--vdj requires -germline_db_V and -germline_db_J source FASTA files.");
  const raw={
    V:await readVdjFasta(paths.V,"-germline_db_V"),
    D:paths.D?await readVdjFasta(paths.D,"-germline_db_D"):"",
    J:await readVdjFasta(paths.J,"-germline_db_J"),
    C:paths.C?await readVdjFasta(paths.C,"-c_region_db"):"",
  };
  const useSwig=Boolean(options["--swigannots"]);
  const internalPath=options["-custom_internal_data"];
  const auxiliaryPath=options["-auxiliary_data"];
  if(useSwig&&(internalPath||auxiliaryPath))throw new Error("--swigannots cannot be combined with -custom_internal_data or -auxiliary_data; choose one annotation source.");
  const scope=vdjScope(options);
  let pack;let species;let allowedLoci;
  const needPack=useSwig||Boolean(auxiliaryPath&&!internalPath);
  if(needPack){
    pack=JSON.parse(gunzipSync(readFileSync(assets.referencePackPath)).toString("utf8"));
    species=vdjSpecies(pack,options["-organism"]);
    allowedLoci=lociForScope(species,scope);
    if(!allowedLoci.length)throw new Error(`The embedded reference pack has no ${scope} loci for ${species.name}.`);
  }

  const references={V:"",D:"",J:"",C:""};
  let mode="assignments-only";
  if(useSwig){
    mode="swig-metadata";
    for(const segment of ["V","D","J","C"]){
      if(!raw[segment])continue;
      const report=preprocessGermlineFastaAcrossTiers(raw[segment],segment,germlineTemplateTiers(pack,species,scope,segment),allowedLoci);
      references[segment]=report.fasta;
      process.stderr.write(`Swig metadata preparation: ${segment} ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} annotated.\n`);
      for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
    }
  }else{
    for(const segment of ["V","D","J","C"]){
      if(raw[segment])references[segment]=prepareIgblastStyleGermlineFasta(raw[segment],segment).fasta;
    }
  }

  if(internalPath){
    mode="igblast-data";
    const application=applyIgblastInternalData(references.V,await readMaybeCompressedText(internalPath,"-custom_internal_data"));
    if(application.matched!==application.total)throw new Error(`${unmatchedMessage("-custom_internal_data",application)} IgBLAST requires custom internal data for every selected V sequence.`);
    references.V=application.fasta;
    process.stderr.write(`IgBLAST internal data: ${application.matched.toLocaleString()}/${application.total.toLocaleString()} V records annotated.\n`);
  }else if(auxiliaryPath){
    // IgBLAST normally obtains V domains from its organism-specific internal
    // directory when only -auxiliary_data is supplied. The standalone Swig
    // equivalent uses the fixed embedded pack and the web metadata transfer.
    mode="igblast-data";
    const report=preprocessGermlineFastaAcrossTiers(raw.V,"V",germlineTemplateTiers(pack,species,scope,"V"),allowedLoci);
    references.V=report.fasta;
    process.stderr.write(`Embedded ${species.name} V metadata: ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} records annotated.\n`);
    for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
  }

  let fwr4EndOffsets={};
  if(auxiliaryPath){
    mode="igblast-data";
    const application=applyIgblastAuxiliaryData(references.J,await readMaybeCompressedText(auxiliaryPath,"-auxiliary_data"));
    if(!application.matched)throw new Error("-auxiliary_data has no exact identifiers matching the selected J FASTA.");
    references.J=application.fasta;fwr4EndOffsets=application.fwr4EndOffsets??{};
    process.stderr.write(`IgBLAST auxiliary data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} J records have CDR3-stop annotations.\n`);
    if(application.unmatched.length)process.stderr.write(`Reference warning: ${unmatchedMessage("-auxiliary_data",application)}\n`);
  }

  const dFramePath=options["-d_frame_data"];
  if(dFramePath){
    if(!references.D)throw new Error("-d_frame_data requires -germline_db_D.");
    const application=applyIgblastDFrameData(references.D,await readMaybeCompressedText(dFramePath,"-d_frame_data"));
    if(!application.matched)throw new Error("-d_frame_data has no exact identifiers matching the selected D FASTA.");
    references.D=application.fasta;
    process.stderr.write(`IgBLAST D-frame data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} D records annotated.\n`);
    if(application.unmatched.length)process.stderr.write(`Reference warning: ${unmatchedMessage("-d_frame_data",application)}\n`);
  }
  return {references,mode,fwr4EndOffsets,jLengths:fastaLengths(references.J)};
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

async function finishWritable(stream){stream.end();await once(stream,"finish");}

async function writeRowsFile(path,headers,rows,include=()=>true){
  const stream=createWriteStream(path);
  try{
    await once(stream,"open");
    await writeChunk(stream,`${headers.join("\t")}\n`);
    let batch=[];
    for(let ordinal=0;ordinal<rows.length;ordinal+=1){
      if(!include(rows[ordinal],ordinal))continue;
      batch.push(rows[ordinal]);
      if(batch.length>=2_000){await writeChunk(stream,serializeRowBody(headers,batch));batch=[];}
    }
    if(batch.length)await writeChunk(stream,serializeRowBody(headers,batch));
    await finishWritable(stream);
  }catch(error){stream.destroy();throw error;}
}

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
  const prefix=config.output.prefix;
  const loadedReferences=await loadReferences(config,base,assets.referencePackPath);
  const references=loadedReferences.references;
  if(config.annotation.airrMode==="preserve"&&config.annotation.doubleD.mode!=="off"&&config.inputs.some((input)=>detectFormat(input.path,input.format)==="airr"))throw new Error("Double-D screening of AIRR input requires annotation.airrMode = \"reannotate\".");
  const needsAnnotation=config.inputs.some((input)=>detectFormat(input.path,input.format)!=="airr"||config.annotation.airrMode==="reannotate");
  const pool=needsAnnotation?new WasmPool(config.annotation.workers,{wasmPath:assets.wasmPath,references,callingProfile:config.annotation.callingProfile,assignerStrategy:config.annotation.assignerStrategy}):null;
  if(pool)await pool.start();
  const annotatedStream=config.output.writeAnnotatedAirr?createWriteStream(join(outputDirectory,`${prefix}.annotated.airr.tsv`)):null;
  let annotatedOutputHeaders=null;
  const rows=[];const headers=[];let annotatedRecords=0;let inputRecords=0;let eligibleRecords=0;
  let fastqFilterStats=emptyFastqQualityFilterStats(config.preprocessing.fastqFilter.enabled,false);
  try{
    if(annotatedStream)await once(annotatedStream,"open");
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
        const batchHeaders=[...table.headers];
        if(doubleD)addHeaders(batchHeaders,doubleD.headers);
        for(const row of table.rows){
          const dd=ddById.get(row.sequence_id);if(dd)Object.assign(row,dd);
          rows.push(row);
        }
        addHeaders(headers,batchHeaders);
        if(annotatedStream){
          if(!annotatedOutputHeaders){annotatedOutputHeaders=batchHeaders;await writeChunk(annotatedStream,`${annotatedOutputHeaders.join("\t")}\n`);}
          else{
            const newHeaders=batchHeaders.filter((header)=>!annotatedOutputHeaders.includes(header));
            if(newHeaders.length)throw new Error(`Streaming annotated AIRR output cannot add columns after its header was written (${newHeaders.join(", ")}). Use matching AIRR schemas for all preserved inputs or reannotate them.`);
          }
          await writeChunk(annotatedStream,serializeRowBody(annotatedOutputHeaders,table.rows));
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
  }catch(error){annotatedStream?.destroy();throw error;}
  finally{await pool?.close();}
  if(annotatedStream)await finishWritable(annotatedStream);
  rows.forEach((row,ordinal)=>{if(!row.sequence_id)row.sequence_id=`swig_${ordinal+1}`;});

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
  await writeRowsFile(join(outputDirectory,`${prefix}.processed.airr.tsv`),headers,rows,(_,ordinal)=>Boolean(activeMask[ordinal]));

  let lineageStudy=null;
  if(config.output.writeLineageStudy&&lineages){
    process.stderr.write("Writing lazy lineage-study bundle…\n");
    lineageStudy=await writeLineageStudy(join(outputDirectory,`${prefix}.lineages.airr.tsv`),join(outputDirectory,`${prefix}.swig-lineage-study.json.gz`),headers,rows,activeMask,lineages,config,references,shm);
  }

  const retained=activeMask.reduce((sum,value)=>sum+(value?1:0),0);
  const summary={application:"swig-cli",version:VERSION,completedAt:new Date().toISOString(),inputRecords,eligibleRecords,annotatedRecords,retainedRecords:retained,lineages:lineages?.lineageCount??0,references:{metadataPreparation:{enabled:config.references.prepareMetadata,segments:loadedReferences.preparation}},preprocessing:{subsample:config.preprocessing.subsample,fastqFilter:{options:config.preprocessing.fastqFilter,stats:fastqFilterStats}},collapse:collapseResult?{mode:collapseResult.mode,inputRecords:collapseResult.inputRecords,inputAbundance:collapseResult.inputAbundance,uniqueRecords:collapseResult.uniqueRecords,collapsedRecords:collapseResult.collapsedRecords,warnings:collapseResult.warnings}:null,chimera:chimeraSummary,alleleRefinement:alleleResult?compactAlleleResult(alleleResult):null,shm:shm?{analyzedRecords:shm.analyzedRecords,analyzedAbundance:shm.analyzedAbundance,skippedRecords:shm.skippedRecords,metric:shm.metric}:null,missingAlleles:missingAlleles?{candidates:missingAlleles.candidates.length,warnings:missingAlleles.warnings}:null,lineageStudy:lineageStudy?{manifest:`${prefix}.swig-lineage-study.json.gz`,airr:`${prefix}.lineages.airr.tsv`,records:lineageStudy.linkedAirr.records}:null};
  await writeFile(join(outputDirectory,`${prefix}.summary.json`),JSON.stringify(summary,null,2));
  await writeFile(join(outputDirectory,`${prefix}.resolved-config.json`),JSON.stringify(config,null,2));
  process.stdout.write(`${JSON.stringify(summary,null,2)}\n`);
}

function vdjTuning(options,callingProfile){
  const keys=["-min_D_match","-min_J_length","-num_alignments_D","-num_alignments_J","-D_penalty","-J_penalty"];
  if(!keys.some((key)=>options[key]!==undefined))return undefined;
  const agreement=callingProfile!=="truth_optimized";
  const minD=options["-min_D_match"]===undefined?(agreement?5:6):parseIntegerOption(options["-min_D_match"],"-min_D_match",{minimum:5});
  const minJ=options["-min_J_length"]===undefined?10:parseIntegerOption(options["-min_J_length"],"-min_J_length",{minimum:0});
  const topD=options["-num_alignments_D"]===undefined?(agreement?3:2):parseIntegerOption(options["-num_alignments_D"],"-num_alignments_D",{minimum:1,allowZero:false});
  const topJ=options["-num_alignments_J"]===undefined?2:parseIntegerOption(options["-num_alignments_J"],"-num_alignments_J",{minimum:1,allowZero:false});
  const dMismatch=options["-D_penalty"]===undefined?(agreement?-4:-3):parseFiniteOption(options["-D_penalty"],"-D_penalty");
  const jMismatch=options["-J_penalty"]===undefined?(agreement?-4:-3):parseFiniteOption(options["-J_penalty"],"-J_penalty");
  if(options["-D_penalty"]!==undefined&&(!Number.isInteger(dMismatch)||dMismatch<=-5||dMismatch>=0))throw new Error("-D_penalty must be an integer greater than -5 and less than 0.");
  if(options["-J_penalty"]!==undefined&&(!Number.isInteger(jMismatch)||jMismatch<=-4||jMismatch>=0))throw new Error("-J_penalty must be an integer greater than -4 and less than 0.");
  return {dMatch:2,dMismatch,dGapOpen:agreement?-11:-13,dGapExtend:-1,topD,minDMatch:minD,jMatch:2,jMismatch,jGapOpen:agreement?-13:-17,jGapExtend:agreement?-1:-2,topJ,minJLength:Math.max(1,minJ)};
}

async function runVdj(rawArgs,assets){
  const options=parseVdjArguments(rawArgs);
  if(options.help){process.stdout.write(`${vdjUsage()}\n`);return;}
  if(options.version){process.stdout.write(`${VERSION}\n`);return;}
  const outputValue=options["-out"];
  if(!outputValue)throw new Error("--vdj requires -out (or --out) so AIRR rows can be streamed to an explicit destination.");
  const outfmt=String(options["-outfmt"]??"19").trim();
  if(outfmt!=="19")throw new Error("swig-cli --vdj currently emits only AIRR rearrangement format; use -outfmt 19.");
  const domain=String(options["-domain_system"]??"imgt").toLowerCase();
  if(domain!=="imgt")throw new Error("swig-cli --vdj supports only -domain_system imgt; Kabat coordinates must not be labeled as IMGT/AIRR annotations.");
  const strandValue=String(options["-strand"]??"both").toLowerCase();
  const strand={both:0,plus:1,minus:2}[strandValue];
  if(strand===undefined)throw new Error("-strand must be both, plus, or minus.");
  const minimumIdentity=options["--minimum-identity"]===undefined?0.6:parseFiniteOption(options["--minimum-identity"],"--minimum-identity");
  if(minimumIdentity<0||minimumIdentity>1)throw new Error("--minimum-identity must be between 0 and 1.");
  const batchRecords=options["--batch-records"]===undefined?2_000:parseIntegerOption(options["--batch-records"],"--batch-records",{minimum:1,allowZero:false});
  const threadValue=options["--workers"]??options["-num_threads"];
  let workers=threadValue===undefined?Math.max(1,Math.min(4,availableParallelism())):parseIntegerOption(threadValue,options["--workers"]!==undefined?"--workers":"-num_threads",{minimum:0});
  if(options["--workers"]===undefined&&threadValue!==undefined&&workers===0)throw new Error("-num_threads must be at least 1; use --workers 0 for automatic selection.");
  if(workers===0)workers=Math.max(1,Math.min(8,availableParallelism()));
  const assigner=String(options["--assigner"]??"aer");
  if(!["standard","riat_mp","aer"].includes(assigner))throw new Error("--assigner must be standard, riat_mp, or aer.");
  const callingProfile=String(options["--calling-profile"]??"truth_optimized");
  if(!["truth_optimized","igblast_compatible","igblast_balanced"].includes(callingProfile))throw new Error("--calling-profile must be truth_optimized, igblast_compatible, or igblast_balanced.");
  const queryValue=String(options["-query"]??"-");
  const queryPath=queryValue==="-"?"-":resolve(queryValue);
  const outputPath=String(outputValue)==="-"?"-":resolve(String(outputValue));
  if(outputPath!=="-")await mkdir(dirname(outputPath),{recursive:true});
  const prepared=await prepareVdjReferences(options,assets);
  const tuning=vdjTuning(options,callingProfile);
  const pool=new WasmPool(workers,{wasmPath:assets.wasmPath,references:prepared.references,callingProfile,assignerStrategy:assigner,tuning});
  await pool.start();
  const output=outputPath==="-"?process.stdout:createWriteStream(outputPath);
  let outputHeader=null;let records=0;let completed=false;
  try{
    if(outputPath!=="-")await once(output,"open");
    process.stderr.write(`Streaming SwiftIG V(D)J assignments (${prepared.mode}; ${workers} worker${workers===1?"":"s"}) to ${outputPath}.\n`);
    const state={inputRecords:0,eligibleRecords:0,selectedRecords:0,fastqFilter:emptyFastqQualityFilterStats(false,false)};
    const pending=[];
    const consume=async(item)=>{
      const result=await item.promise;
      if(outputHeader===null){outputHeader=result.header;await writeChunk(output,`${outputHeader}\n`);}
      else if(result.header!==outputHeader)throw new Error("SwiftIG changed the AIRR schema between V(D)J batches.");
      const body=trimAuxiliaryFwr4(result.header,result.body,prepared.fwr4EndOffsets,prepared.jLengths);
      await writeChunk(output,body);records+=result.count;
    };
    const preprocessing={fastqFilter:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter,trim3Prime:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter.trim3Prime}},subsample:{...DEFAULT_CLI_CONFIG.preprocessing.subsample,enabled:false}};
    const input={path:queryPath,format:"fasta"};
    for await(const batch of sequenceBatches(input,batchRecords,preprocessing,"reannotate",0,state)){
      pending.push({promise:pool.run({text:batch.text,count:batch.count,format:1,minimumIdentity,strand,doubleD:{mode:"off"}})});
      if(pending.length>=Math.max(2,workers*2))await consume(pending.shift());
    }
    while(pending.length)await consume(pending.shift());
    if(outputPath!=="-")await finishWritable(output);
    completed=true;
  }finally{
    if(!completed&&outputPath!=="-")output.destroy();
    await pool.close();
  }
  process.stderr.write(`Completed ${records.toLocaleString()} streaming V(D)J assignment${records===1?"":"s"}.\n`);
}

export async function runCli(assets=defaultCliAssets()){
  const args=process.argv.slice(2);
  if(args.includes("--vdj")){await runVdj(args,assets);return;}
  const command=args[0]&&!args[0].startsWith("-")?args[0]:"run";const rest=command===args[0]?args.slice(1):args;
  if(hasFlag(args,"--help")||command==="help"){process.stdout.write(`${usage()}\n`);return;}
  if(hasFlag(args,"--version")||command==="version"){process.stdout.write(`${VERSION}\n`);return;}
  if(hasFlag(args,"--notices")||command==="notices"){process.stdout.write(`${thirdPartyNotices()}\n`);return;}
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
  const commandOutput=argumentValue(rest,"--out");
  const configuredOutput=typeof raw.output?.directory==="string"&&raw.output.directory.trim();
  if(!commandOutput&&!configuredOutput)throw new Error("The pipeline CLI requires an explicit output directory in output.directory or --out so results can be streamed to disk.");
  if(commandOutput)raw.output={...(raw.output??{}),directory:commandOutput};
  const commandWorkers=argumentValue(rest,"--workers");
  if(commandWorkers!==undefined)raw.annotation={...(raw.annotation??{}),workers:parseIntegerOption(commandWorkers,"--workers",{minimum:0})};
  const config=normalizeCliConfig(raw);
  if(config.annotation.workers===0)config.annotation.workers=Math.max(1,Math.min(8,availableParallelism()));
  config.inputs=config.inputs.map((input)=>({...input,path:typeof input.inline==="string"?input.path:resolveFrom(base,input.path)}));
  await runPipeline(config,base,assets);
}
