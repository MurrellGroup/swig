import assert from "node:assert/strict";
import test from "node:test";

import { cliConfigFromBrowser, normalizeCliConfig } from "../src/pipeline-config.ts";
import type { CompiledReferences } from "../src/reference-pack.ts";
import { DEFAULT_PIPELINE_PLAN } from "../src/study-design.ts";

const references:CompiledReferences={
  V:">IGHV1-2*01\nACGT\n",
  D:">IGHD1-1*01\nAC\n",
  J:">IGHJ4*02\nTGCA\n",
  C:"",
  counts:{V:1,D:1,J:1,C:0},
  annotation:{V:{annotated:0,total:1},J:{annotated:0,total:1}},
  loci:["IGH"],
};

test("small CLI configs receive analysis defaults and explicit donor grouping",()=>{
  const config=normalizeCliConfig({inputs:[
    {path:"day 0.fastq.gz",sampleId:"day_0",subjectId:"donor A"},
    {path:"day 30.fastq.gz",sampleId:"day_30",subjectId:"donor A"},
  ]});
  assert.equal(config.inputs[0].subjectId,"donor-A");
  assert.equal(config.inputs[1].subjectId,"donor-A");
  assert.equal(config.pipeline.collapse.enabled,true);
  assert.equal(config.pipeline.selection.enabled,false);
  assert.equal(config.pipeline.lineage.enabled,true);
  assert.equal(config.pipeline.lineage.scope,"subject");
  assert.equal(config.pipeline.shm.enabled,true);
  assert.equal(config.annotation.workers,0,"zero means choose the host default at CLI startup");
});

test("browser CLI export preserves exact references, methods, and sample-to-donor map",()=>{
  const config=cliConfigFromBrowser({
    studyName:"vaccine study",
    studyDesign:"longitudinal",
    datasets:[
      {datasetId:"d0",inputName:"day0.fastq.gz",sampleId:"day_0",subjectId:"person_1",cohort:"vaccinated",timepoint:"0",compartment:"blood"},
      {datasetId:"d30",inputName:"day30.fastq.gz",sampleId:"day_30",subjectId:"person_1",cohort:"vaccinated",timepoint:"30",compartment:"lymph"},
    ],
    references,
    species:"Homo sapiens",
    scope:"IGH",
    workers:3,
    callingProfile:"truth_optimized",
    assignerStrategy:"aer",
    minimumIdentity:0.7,
    strand:0,
    doubleD:{mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8},
    pipeline:{...DEFAULT_PIPELINE_PLAN,enabled:true,selection:{...DEFAULT_PIPELINE_PLAN.selection,enabled:true,sampleId:"day_30"}},
  });
  assert.equal(config.output.prefix,"vaccine-study");
  assert.equal(config.inputs[0].subjectId,config.inputs[1].subjectId);
  assert.equal(config.references.inline?.V,references.V);
  assert.equal(config.references.inline?.J,references.J);
  assert.equal(config.pipeline.selection.enabled,true);
  assert.equal(config.pipeline.selection.sampleId,"day_30");
  assert.equal(config.annotation.workers,3);
});
