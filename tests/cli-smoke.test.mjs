import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { join, resolve } from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { gunzipSync } from "node:zlib";

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
    const standalone=process.env.SWIG_CLI_EXECUTABLE;
    const command=standalone||process.execPath;
    const arguments_=standalone?["run","--config",configPath]:[join(root,"cli/swig-cli.mjs"),"run","--config",configPath];
    const result=spawnSync(command,arguments_,{encoding:"utf8",timeout:120_000});
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
