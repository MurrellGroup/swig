import { readFileSync } from "node:fs";
import { parentPort } from "node:worker_threads";
import { WASI } from "@bjorn3/browser_wasi_shim";
import { applyBalancedDFilter, reconcileBalancedDoubleD } from "../src/balanced-calling-profile.ts";

const encoder=new TextEncoder();
const decoder=new TextDecoder();
let runtime=null;
let callingProfile="truth_optimized";

function put(bytes){
  const pointer=runtime.swig_alloc(bytes.byteLength);
  if(!pointer&&bytes.byteLength)throw new Error("SwiftIG ran out of WebAssembly memory.");
  new Uint8Array(runtime.memory.buffer,pointer,bytes.byteLength).set(bytes);
  return [pointer,bytes.byteLength];
}

function read(pointer,length){return decoder.decode(new Uint8Array(runtime.memory.buffer,pointer,length));}
function errorText(){return read(runtime.swig_error_ptr(),runtime.swig_error_len())||"SwiftIG could not complete the annotation.";}

async function initialize(message){
  if(message.callingProfile==="r_optimized"&&message.assignerStrategy!=="aer_robust")throw new Error("The r-optimized calling profile requires AER-R.");
  const wasi=new WASI([],[],[]);
  const module=await WebAssembly.compile(readFileSync(message.wasmPath));
  const instance=await WebAssembly.instantiate(module,{wasi_snapshot_preview1:wasi.wasiImport});
  wasi.initialize(instance);
  runtime=instance.exports;
  const strategy=message.assignerStrategy==="standard"?0:message.assignerStrategy==="aer"?2:message.assignerStrategy==="aer_robust"?3:1;
  if(runtime.swig_set_assigner_strategy(strategy)!==0)throw new Error("SwiftIG rejected the assignment strategy.");
  // The browser's balanced profile intentionally starts from compatible calls
  // and then applies the shared post-filter below.
  const profile=message.callingProfile==="truth_optimized"?0:message.callingProfile==="r_optimized"?2:1;
  if(runtime.swig_set_calling_profile(profile)!==0)throw new Error("SwiftIG rejected the calling profile.");
  callingProfile=message.callingProfile;
  if(message.hasTuning){
    if(typeof runtime.swig_set_tuning_options!=="function")throw new Error("This SwiftIG build does not expose the requested D/J compatibility controls.");
    const accepted=runtime.swig_set_tuning_options(
      message.tuningDMatch,message.tuningDMismatch,message.tuningDGapOpen,message.tuningDGapExtend,message.tuningTopD,message.tuningMinDMatch,
      message.tuningJMatch,message.tuningJMismatch,message.tuningJGapOpen,message.tuningJGapExtend,message.tuningTopJ,message.tuningMinJLength,
    );
    if(accepted!==0)throw new Error("SwiftIG rejected the requested D/J compatibility controls.");
    if(message.rOptimized){
      if(typeof runtime.swig_set_v_tuning_options!=="function"||typeof runtime.swig_set_aer_r_decision_tuning!=="function")throw new Error("This SwiftIG build does not expose r-optimized tuning controls.");
      if(runtime.swig_set_v_tuning_options(2,-4,-13,-1,1)!==0||runtime.swig_set_aer_r_decision_tuning(10)!==0)throw new Error("SwiftIG rejected the r-optimized tuning controls.");
    }
  }
  const allocations=[message.referenceV,message.referenceD,message.referenceJ,message.referenceC].map((value)=>put(encoder.encode(value||"")));
  try{
    const genes=runtime.swig_init_database(...allocations.flat());
    if(genes<0)throw new Error(errorText());
    return {genes};
  }finally{allocations.forEach(([pointer])=>runtime.swig_free(pointer));}
}

function annotate(message){
  const bytes=encoder.encode(message.text);
  const [pointer,size]=put(bytes);
  let count;
  try{
    count=message.doubleDMode!=="off"
      ? runtime.swig_annotate_double_d(pointer,size,message.format,Math.round(message.minimumIdentity*1000),message.strand,message.doubleDMode==="all"?1:2,Math.round(message.doubleDMinimumVjSpan),Math.round(message.doubleDSeedLength),Math.round(message.doubleDPseudoTrim),Math.round(message.doubleDMaximumPseudoMismatches),Math.round(message.doubleDMinimumScoreGain))
      : runtime.swig_annotate(pointer,size,message.format,Math.round(message.minimumIdentity*1000),message.strand);
  }finally{runtime.swig_free(pointer);}
  if(count<0)throw new Error(errorText());
  if(count!==message.count)throw new Error(`SwiftIG returned ${count} rows for a ${message.count}-record batch.`);
  const result=new Uint8Array(runtime.memory.buffer,runtime.swig_result_ptr(),runtime.swig_result_len());
  const newline=result.indexOf(10);
  if(newline<0)throw new Error("SwiftIG returned an invalid AIRR table.");
  const header=decoder.decode(result.subarray(0,newline)).replace(/\r$/,"");
  let body=result.subarray(newline+1);
  const balanced=callingProfile==="igblast_balanced"?applyBalancedDFilter(header,body):null;
  if(balanced)body=balanced.body;
  const response={header,body:decoder.decode(body),count,doubleDHeader:"",doubleDBody:""};
  if(message.doubleDMode!=="off"){
    const dd=new Uint8Array(runtime.memory.buffer,runtime.swig_double_d_result_ptr(),runtime.swig_double_d_result_len());
    const ddNewline=dd.indexOf(10);
    if(ddNewline<0)throw new Error("SwiftIG returned invalid double-D evidence.");
    response.doubleDHeader=decoder.decode(dd.subarray(0,ddNewline)).replace(/\r$/,"");
    let ddBody=dd.slice(ddNewline+1);
    if(balanced)ddBody=reconcileBalancedDoubleD(response.doubleDHeader,ddBody,balanced.suppressedSequenceIds);
    response.doubleDBody=decoder.decode(ddBody);
  }
  return response;
}

const send=(message)=>parentPort?parentPort.postMessage(message):globalThis.postMessage(message);
const receive=async(message)=>{
  try{
    const result=message.type==="init"?await initialize(message):annotate(message);
    send({id:message.id,...result});
  }catch(error){send({id:message.id,error:error instanceof Error?error.message:String(error)});}
};

if(parentPort)parentPort.on("message",receive);
else globalThis.addEventListener("message",(event)=>void receive(event.data));
