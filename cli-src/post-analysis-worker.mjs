import { parentPort } from "node:worker_threads";

import { runDenoisePartitionJob, runExactDedupJob } from "../src/post-analysis-core.ts";

const send=(message,transfer)=>parentPort?parentPort.postMessage(message,transfer):globalThis.postMessage(message,transfer);
const receive=(message)=>{
  try{
    if(message.kind==="denoise"){
      const result=runDenoisePartitionJob(message.job);
      send({id:message.id,result},[result.targets.buffer]);
    }else if(message.kind==="exact"){
      const result=runExactDedupJob(message.job);
      send({id:message.id,result},[result.representatives.buffer,result.representativeOrdinals.buffer,result.representativeCounts.buffer]);
    }else throw new Error(`Unknown post-analysis worker task: ${String(message.kind)}`);
  }catch(error){send({id:message.id,error:error instanceof Error?error.message:String(error)});}
};

if(parentPort)parentPort.on("message",receive);
else globalThis.addEventListener("message",(event)=>receive(event.data));
