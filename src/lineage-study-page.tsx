import { useEffect, useMemo, useRef, useState } from "react";

import { PostAnalysisWorkbench } from "./post-analysis";
import type { LineageSummary } from "./post-analysis-core";
import { readLineageAirrSlice, type LineageStudyManifest } from "./lineage-study";
import { AirrResultStore, type AirrDetailRow } from "./result-store";
import { createSampleColorMap } from "./sample-colors";
import { packSessionVector, type PostAnalysisSessionSnapshot } from "./session-state";

interface LoadedLineage {
  id:number;
  store:AirrResultStore;
  session:PostAnalysisSessionSnapshot;
  summary:LineageSummary;
}

function shmText(value:number|undefined,metric:string):string {
  if(value===undefined)return "—";
  return metric.toLowerCase().includes("rate")?`${(value*100).toFixed(1)}%`:value.toFixed(value<10?2:0);
}

export function LineageStudyPage({manifest,airrFile,onClose}:{manifest:LineageStudyManifest;airrFile:File;onClose:()=>void}){
  const [search,setSearch]=useState("");
  const [sample,setSample]=useState("");
  const [minimumAbundance,setMinimumAbundance]=useState(0);
  const [loaded,setLoaded]=useState<LoadedLineage|null>(null);
  const [loading,setLoading]=useState("");
  const [error,setError]=useState("");
  const [detail,setDetail]=useState<AirrDetailRow|null>(null);
  const activeStoreRef=useRef<AirrResultStore|null>(null);
  useEffect(()=>()=>{void activeStoreRef.current?.clear();},[]);
  const samples=useMemo(()=>[...new Set(manifest.summaries.flatMap((summary)=>summary.sampleIds??[]))].sort((left,right)=>left.localeCompare(right,undefined,{numeric:true})),[manifest]);
  const shm=useMemo(()=>new Map(manifest.shm?.summaries.map((value)=>[value.lineageId,value])??[]),[manifest]);
  const summaries=useMemo(()=>{
    const needle=search.trim().toLowerCase();
    return manifest.summaries.filter((summary)=>{
      if(sample&&!(summary.sampleIds??[]).includes(sample))return false;
      if(minimumAbundance&&summary.abundance<minimumAbundance)return false;
      return !needle||`${summary.id} ${summary.studyGroup} ${summary.locus} ${summary.vCalls.join(" ")} ${summary.jCalls.join(" ")} ${(summary.sampleIds??[]).join(" ")}`.toLowerCase().includes(needle);
    });
  },[manifest,minimumAbundance,sample,search]);

  async function openLineage(summary:LineageSummary){
    setLoading(`Loading lineage ${summary.id} from its AIRR byte range`);setError("");setDetail(null);
    try{
      const slice=await readLineageAirrSlice(airrFile,manifest,summary.id);
      const store=new AirrResultStore();store.configureStudyMetadataForImport(manifest.analysis.datasets);await store.appendBatch(slice.header,slice.body);await store.finalize();
      const previous=activeStoreRef.current;activeStoreRef.current=store;if(previous)await previous.clear();
      const assignments=new Int32Array(store.count);assignments.fill(summary.id);
      const adjusted={...summary,representativeOrdinal:0};
      const dashboard={summaries:[adjusted],lineageCount:1,sizeHistogram:[{label:"Selected lineage",count:1}],vUsage:adjusted.vCalls.map((call)=>({call,lineages:1,abundance:adjusted.abundance})),jUsage:adjusted.jCalls.map((call)=>({call,lineages:1,abundance:adjusted.abundance})),assignedRecords:adjusted.abundance,unassignedRecords:0,candidateComparisons:0,truncatedCandidates:0};
      const session:PostAnalysisSessionSnapshot={skippedModules:["alleles","dedup","chimera","selection","lineage","diagnostics","query"],workingStages:[],lineage:{options:{identity:manifest.analysis.lineage.identity,resolution:manifest.analysis.lineage.resolution,ambiguity:manifest.analysis.lineage.ambiguity,productiveOnly:manifest.analysis.lineage.productiveOnly,lineageScope:manifest.analysis.lineage.scope},assignments:packSessionVector(assignments),dashboard},selectedLineageIds:[summary.id]};
      setLoaded({id:summary.id,store,session,summary:adjusted});
    }catch(problem){setError(problem instanceof Error?problem.message:String(problem));}
    finally{setLoading("");}
  }

  async function inspect(ordinal:number){if(!loaded)return;const values=await loaded.store.detailMany([ordinal]);setDetail(values[0]??null);}
  const sampleColors=useMemo(()=>createSampleColorMap(manifest.analysis.datasets),[manifest]);

  return <main className="lineage-study-page">
    <header className="lineage-study-heading"><div><span className="section-kicker">Lazy lineage-study session</span><h1>{manifest.analysis.inputName}</h1><p>{manifest.summaries.length.toLocaleString()} indexed lineage summaries · linked AIRR rows are read only after a lineage is selected.</p></div><div><strong>{airrFile.name}</strong><small>{airrFile.size.toLocaleString()} bytes · no whole-file indexing</small><button type="button" onClick={onClose}>Close lineage study</button></div></header>
    <section className="lineage-study-browser">
      <header><div><span className="section-kicker">Choose a lineage</span><h2>{summaries.length.toLocaleString()} matching lineages</h2></div><div className="lineage-study-filters"><label><span>Search</span><input type="search" value={search} onChange={(event)=>setSearch(event.target.value)} placeholder="lineage, V, J, sample…"/></label><label><span>Contains sample</span><select value={sample} onChange={(event)=>setSample(event.target.value)}><option value="">Any sample</option>{samples.map((value)=><option key={value}>{value}</option>)}</select></label><label><span>Minimum reads</span><input type="number" min="0" value={minimumAbundance||""} placeholder="Any" onChange={(event)=>setMinimumAbundance(Math.max(0,Number(event.target.value)||0))}/></label></div></header>
      {loading&&<div className="post-progress"><div><span>{loading}</span><strong>working</strong></div><progress/></div>}
      {error&&<p className="run-error" role="alert">{error}</p>}
      <div className="lineage-table-wrap"><table><thead><tr><th>Lineage</th><th>Reads</th><th>Unique</th><th>Samples and read counts</th><th>Locus</th><th>V calls</th><th>J calls</th>{manifest.shm&&<><th>Mean SHM</th><th>SHM upper q95</th></>}<th/></tr></thead><tbody>{summaries.slice(0,1000).map((summary)=><tr key={summary.id} className={loaded?.id===summary.id?"selected":""} onClick={()=>void openLineage(summary)}><td><strong>{summary.id}</strong><small>{summary.studyGroup}</small></td><td>{summary.abundance.toLocaleString()}</td><td>{summary.uniqueMembers.toLocaleString()}</td><td className="lineage-sample-counts">{(summary.sampleCounts??[]).map((value)=><small key={value.sampleId}>{value.sampleId}: {value.abundance.toLocaleString()} reads</small>)}</td><td>{summary.locus}</td><td>{summary.vCalls.join(", ")}</td><td>{summary.jCalls.join(", ")}</td>{manifest.shm&&<><td>{shmText(shm.get(summary.id)?.mean,manifest.shm.metric)}</td><td>{shmText(shm.get(summary.id)?.p95,manifest.shm.metric)}</td></>}<td><button type="button">Load lineage →</button></td></tr>)}</tbody></table></div>
    </section>
    {loaded&&<section className="lineage-study-workbench"><PostAnalysisWorkbench key={loaded.id} store={loaded.store} references={manifest.analysis.references} scope={manifest.analysis.scope} loci={loaded.store.facets().loci} resultFacets={loaded.store.facets()} inputName={manifest.analysis.inputName} workers={Math.max(1,Math.min(4,navigator.hardwareConcurrency||1))} callingProfile={manifest.analysis.callingProfile} assignerStrategy={manifest.analysis.assignerStrategy} minimumIdentity={manifest.analysis.minimumIdentity} strand={manifest.analysis.strand} datasets={manifest.analysis.datasets} sampleColors={sampleColors} defaultCollapseScope="sample" defaultLineageScope={manifest.analysis.lineage.scope} onInspect={(ordinal)=>void inspect(ordinal)} initialSession={loaded.session}/></section>}
    {detail&&<div className="output-modal-backdrop" role="presentation" onMouseDown={(event)=>{if(event.target===event.currentTarget)setDetail(null);}}><section className="output-modal lineage-study-detail" role="dialog" aria-modal="true"><button className="output-modal-close" type="button" onClick={()=>setDetail(null)}>×</button><span className="section-kicker">Selected AIRR record</span><h2>{detail.values.sequence_id}</h2><p>{detail.values.sample_id} · {detail.values.v_call} · {detail.values.d_call} · {detail.values.j_call}</p><code>{detail.values.cdr3_aa||detail.values.cdr3||"CDR3 unavailable"}</code></section></div>}
  </main>;
}
