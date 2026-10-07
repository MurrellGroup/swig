import assert from "node:assert/strict";
import fs from "node:fs";
import zlib from "node:zlib";
import readline from "node:readline";
import { WASI } from "@bjorn3/browser_wasi_shim";
import { generateVdjDataset, recordsToFasta, seededRandom } from "./vdj-simulator.mjs";
const before = process.argv[2];
if (!before) throw new Error("Usage: node tests/benchmark-junction-anchor.mjs BEFORE.wasm");
const after = new URL("../public/swiftig.wasm", import.meta.url);
const pack = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz", import.meta.url))));
let referenceText;
const fasta = (records) => records.map(([id,seq,meta]) => `>${id}${meta ? ` SWIGMETA=${meta.join(",")}` : ""}\n${seq}\n`).join("");
function parseAirr(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const headers = lines.shift().split("\t");
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
  });
}

async function makeRuntime(binary, strategy) {
  const bytes = fs.readFileSync(binary);
  const wasi = new WASI([], [], []);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
  wasi.initialize(instance);
  const runtime = instance.exports;
  const encoder = new TextEncoder();
  const decoder = new TextDecoder();
  const put = (value) => {
    const bytesValue = encoder.encode(value);
    const pointer = runtime.swig_alloc(bytesValue.byteLength);
    if (!pointer && bytesValue.byteLength) throw new Error("WASM allocation failed");
    new Uint8Array(runtime.memory.buffer, pointer, bytesValue.byteLength).set(bytesValue);
    return [pointer, bytesValue.byteLength];
  };
  const error = () => decoder.decode(new Uint8Array(
    runtime.memory.buffer, runtime.swig_error_ptr(), runtime.swig_error_len(),
  ));
  if (runtime.swig_set_assigner_strategy(strategy) !== 0) {
    throw new Error(`WASM rejected strategy ${strategy}`);
  }
  if (runtime.swig_set_calling_profile(2) !== 0) throw new Error("WASM rejected the R-optimized profile");
  const allocations = [referenceText.V, referenceText.D, referenceText.J, ""].map(put);
  const genes = runtime.swig_init_database(...allocations.flat());
  allocations.forEach(([pointer]) => runtime.swig_free(pointer));
  if (genes < 0) throw new Error(error());
  return {
    genes,
    annotate(fasta) {
      const [pointer, length] = put(fasta);
      const returned = runtime.swig_annotate(pointer, length, 1, 600, 0);
      runtime.swig_free(pointer);
      if (returned < 0) throw new Error(error());
      const result = decoder.decode(new Uint8Array(
        runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len(),
      ));
      return parseAirr(result);
    },
  };
}


const derived = new Set(["productive", "vj_in_frame", "stop_codon", "v_frameshift", "j_frameshift", "complete_vdj", "sequence_frame", "sequence_aa", "junction", "junction_aa", "junction_length", "junction_aa_length", "cdr3", "cdr3_aa", "cdr3_start", "cdr3_end", "fwr4", "fwr4_aa", "fwr4_start", "fwr4_end", "d_frame", "sequence_alignment_aa", "germline_alignment_aa"]);
async function compare(name, records, validate) {
  const original = await makeRuntime(before, 3);
  const changed = await makeRuntime(after, 3);
  const stats = { name, count: records.length, rescued: 0, productiveRescues: 0, falseAnchor: 0, junctionRestoredCodons: 0, productiveJunctionRestoredCodons: 0, existingAnnotationChanges: 0, assignmentChanges: 0, assignmentChangedFields: {} };
  for (let offset = 0; offset < records.length; offset += 100) {
    const batch = records.slice(offset, offset + 100);
    const input = recordsToFasta(batch);
    const oldRows = original.annotate(input), newRows = changed.annotate(input);
    assert.equal(newRows.length, oldRows.length);
    for (let i = 0; i < newRows.length; ++i) {
      const old = oldRows[i], row = newRows[i];
      for (const field of Object.keys(old)) if (!derived.has(field) && old[field] !== row[field]) { stats.assignmentChanges++; stats.assignmentChangedFields[field] = (stats.assignmentChangedFields[field] || 0) + 1; }
      if (old.productive && JSON.stringify(old) !== JSON.stringify(row)) stats.existingAnnotationChanges++;
      if (!old.cdr3 && row.cdr3) {
        stats.rescued++;
        if (row.productive === "T") stats.productiveRescues++;
        if (validate) {
          const evidence = validate(batch[i], row);
          if (!evidence.coordinateCorrect) stats.falseAnchor++;
          if (!evidence.fullyVEncoded) {
            stats.junctionRestoredCodons++;
            if (row.productive === "T") stats.productiveJunctionRestoredCodons++;
          }
        }
      }
    }
    if ((offset + batch.length) % 1000 === 0) process.stderr.write(`${name}: ${offset + batch.length}/${records.length}\n`);
  }
  console.log(JSON.stringify(stats));
  assert.equal(stats.assignmentChanges, 0);
  assert.equal(stats.existingAnnotationChanges, 0);
  return stats;
}
for (const species of ["Vicugna pacos", "Homo sapiens"]) {
  const locus = pack.species.find((entry) => entry.name === species).loci.IGH;
  referenceText = { V: fasta(locus.V), D: fasta(locus.D), J: fasta(locus.J), C: "" };
  const refs = Object.fromEntries(["V", "D", "J"].map((seg) => [seg, locus[seg].map(([name,sequence]) => ({name,sequence}))]));
  const records = generateVdjDataset({ ...refs, count: 3000, seed: 765127, indelRate: 0, sequencingErrorRate: 0, ambiguousRate: 0, reverseRate: 0, doubleDRate: 0 });
  const stats = await compare(`${species} simulated`, records, (record, row) => {
    const allele = locus.V.find((a) => a[0] === record.truth.vCall);
    const anchor = allele[2]?.[11] - 3;
    if (!(anchor >= 0) || allele[1].slice(anchor,anchor+3).match(/^TG[TC]$/) === null) return {coordinateCorrect: false, fullyVEncoded: false};
    const expected = record.truth.spans.V.start + anchor - record.truth.trims.read5;
    return {coordinateCorrect: Number(row.cdr3_start) - 4 === expected, fullyVEncoded: anchor + 3 <= record.truth.retained.V};
  });
  assert.equal(stats.falseAnchor, 0, `${species}: incorrect anchor coordinate`);
  assert.equal(stats.productiveJunctionRestoredCodons, 0, `${species}: productive rescue without a retained V codon`);
}
const kimdb = new URL("../benchmarks/personalized-germline/v0388/", import.meta.url);
referenceText = Object.fromEntries(["V", "D", "J"].map((s) => [s, fs.readFileSync(new URL(`kimdb.${s}.fasta`, kimdb), "utf8")]));
referenceText.C = "";
for (const animal of ["ERR4238110", "ERR4238104"]) {
  const input = fs.createReadStream(new URL(`inputs/${animal}.full.airr.tsv.gz`, kimdb)).pipe(zlib.createGunzip());
  const lines = readline.createInterface({input, crlfDelay: Infinity});
  const records = [], random = seededRandom(75612);
  let header, count = 0;
  for await (const line of lines) {
    if (!header) { header = line.split("\t"); continue; }
    const values = line.split("\t");
    const record = { id: values[header.indexOf("sequence_id")], sequence: values[header.indexOf("sequence")] };
    ++count;
    if (records.length < 2000) records.push(record);
    else { const slot = Math.floor(random() * count); if (slot < 2000) records[slot] = record; }
  }
  await compare(`${animal} real paired`, records);
}
