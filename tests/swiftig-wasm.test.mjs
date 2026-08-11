import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import { WASI } from "@bjorn3/browser_wasi_shim";

const pack = JSON.parse(
  zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7.json.gz", import.meta.url))),
);
const wasmBytes = fs.readFileSync(new URL("../public/swiftig.wasm", import.meta.url));
const encoder = new TextEncoder();
const decoder = new TextDecoder();

function asFasta(alleles) {
  return alleles.map(([name, sequence, metadata]) => `>${name}${metadata ? ` SWIGMETA=${metadata.join(",")}` : ""}\n${sequence}\n`).join("");
}

function reverseComplement(sequence) {
  const complements = { A: "T", C: "G", G: "C", T: "A", N: "N" };
  return [...sequence].reverse().map((base) => complements[base] ?? "N").join("");
}

async function makeRuntime() {
  const wasi = new WASI([], [], []);
  const wasmModule = await WebAssembly.compile(wasmBytes);
  const instance = await WebAssembly.instantiate(wasmModule, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  wasi.initialize(instance);
  const exports = instance.exports;

  function put(value) {
    const bytes = encoder.encode(value);
    const pointer = exports.swig_alloc(bytes.length);
    assert.ok(pointer || bytes.length === 0, "WASM allocation failed");
    new Uint8Array(exports.memory.buffer, pointer, bytes.length).set(bytes);
    return [pointer, bytes.length];
  }

  function read(pointer, length) {
    return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
  }

  function initialize(references) {
    const allocations = [references.V, references.D ?? "", references.J, references.C ?? ""].map(put);
    const count = exports.swig_init_database(...allocations.flat());
    allocations.forEach(([pointer]) => exports.swig_free(pointer));
    if (count < 0) {
      throw new Error(read(exports.swig_error_ptr(), exports.swig_error_len()));
    }
    return count;
  }

  function annotate(query, format, strand = 0) {
    const [pointer, size] = put(query);
    const count = exports.swig_annotate(pointer, size, format, 600, strand);
    exports.swig_free(pointer);
    if (count < 0) {
      throw new Error(read(exports.swig_error_ptr(), exports.swig_error_len()));
    }
    const tsv = read(exports.swig_result_ptr(), exports.swig_result_len());
    const lines = tsv.trimEnd().split("\n");
    const headers = lines.shift().split("\t");
    const rows = lines.map((line) => {
      const values = line.split("\t");
      return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
    });
    return { count, headers, rows, tsv };
  }

  return { initialize, annotate };
}

function referenceFor(locus, customJName) {
  const v = locus.V.find((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0)) ?? locus.V[0];
  const d = locus.D?.[0];
  const originalJ = locus.J.find((allele) => allele[2]?.[0] >= 0 && allele[2]?.[1] >= 0) ?? locus.J[0];
  const j = customJName ? [customJName, originalJ[1], originalJ[2]] : originalJ;
  return {
    references: {
      V: asFasta([v]),
      D: d ? asFasta([d]) : "",
      J: asFasta([j]),
    },
    sequence: `${v[1]}AACCGG${d?.[1] ?? ""}TTG${j[1]}`,
    names: { V: v[0], D: d?.[0] ?? "", J: j[0] },
  };
}

test("reference pack covers complete IG and TR loci", () => {
  assert.ok(pack.species.length >= 60);
  const observed = new Set();
  for (const species of pack.species) {
    for (const [locusName, locus] of Object.entries(species.loci)) {
      observed.add(locusName);
      assert.ok(locus.V.length, `${species.name} ${locusName} has no V records`);
      assert.ok(locus.J.length, `${species.name} ${locusName} has no J records`);
    }
  }
  assert.deepEqual([...observed].sort(), ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]);
});

test("WASM annotates FASTA, FASTQ, and AIRR; handles heavy, light, TCR, strand, and J-only swaps", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human, "human reference set is missing");
  const runtime = await makeRuntime();

  const customJName = "IGHJ_SWIGTEST*01";
  const heavy = referenceFor(human.loci.IGH, customJName);
  assert.equal(runtime.initialize(heavy.references), heavy.names.D ? 3 : 2);

  const fasta = runtime.annotate(`>fasta_case\n${heavy.sequence}\n`, 1);
  assert.equal(fasta.count, 1);
  assert.equal(fasta.rows[0].sequence_id, "fasta_case");
  assert.equal(fasta.rows[0].locus, "IGH");
  assert.equal(fasta.rows[0].j_call, customJName);
  assert.ok(fasta.headers.includes("junction_aa"));
  assert.ok(fasta.headers.includes("germline_alignment"));
  assert.equal(fasta.rows[0].region_definition, "IMGT");
  assert.equal(fasta.rows[0].v_annotation_source, "IMGT-gapped");
  assert.equal(fasta.rows[0].j_annotation_source, "validated-J-motif");
  assert.ok(Number(fasta.rows[0].sequence_frame) >= 1 && Number(fasta.rows[0].sequence_frame) <= 3);
  for (const region of ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]) {
    assert.ok(fasta.rows[0][region], `${region.toUpperCase()} nucleotide sequence is empty`);
    assert.ok(fasta.rows[0][`${region}_aa`], `${region.toUpperCase()} amino-acid sequence is empty`);
    assert.ok(Number(fasta.rows[0][`${region}_start`]) > 0);
    assert.ok(Number(fasta.rows[0][`${region}_end`]) >= Number(fasta.rows[0][`${region}_start`]));
  }
  assert.ok(fasta.rows[0].cdr3);
  assert.ok(fasta.rows[0].cdr3_aa);
  assert.equal(fasta.rows[0].sequence_alignment_aa.length, fasta.rows[0].germline_alignment_aa.length);
  for (const segment of ["v", "j"]) {
    assert.ok(fasta.rows[0][`${segment}_sequence_alignment`], `${segment.toUpperCase()} query alignment is empty`);
    assert.ok(fasta.rows[0][`${segment}_germline_alignment`], `${segment.toUpperCase()} germline alignment is empty`);
    assert.equal(
      fasta.rows[0][`${segment}_sequence_alignment`].length,
      fasta.rows[0][`${segment}_germline_alignment`].length,
      `${segment.toUpperCase()} alignment strings differ in length`,
    );
  }
  if (heavy.names.D) {
    assert.ok(fasta.rows[0].d_sequence_alignment);
    assert.ok(fasta.rows[0].d_germline_alignment);
  }

  const fastq = runtime.annotate(
    `@fastq_case\n${heavy.sequence}\n+\n${"I".repeat(heavy.sequence.length)}\n`,
    2,
  );
  assert.equal(fastq.rows[0].sequence_id, "fastq_case");
  assert.equal(fastq.rows[0].quality.length, heavy.sequence.length);

  const airr = runtime.annotate(`sequence_id\tsequence\nairr_case\t${heavy.sequence}\n`, 3);
  assert.equal(airr.rows[0].sequence_id, "airr_case");
  assert.equal(airr.rows[0].j_call, customJName);

  const reverse = runtime.annotate(`>reverse_case\n${reverseComplement(heavy.sequence)}\n`, 1, 2);
  assert.equal(reverse.rows[0].rev_comp, "T");
  assert.equal(reverse.rows[0].j_call, customJName);

  const light = referenceFor(human.loci.IGK);
  assert.equal(runtime.initialize(light.references), 2);
  const lightResult = runtime.annotate(`>light_case\n${light.sequence}\n`, 1);
  assert.equal(lightResult.rows[0].locus, "IGK");
  assert.equal(lightResult.rows[0].d_call, "");
  assert.ok(lightResult.rows[0].v_call);
  assert.ok(lightResult.rows[0].j_call);

  const tcr = referenceFor(human.loci.TRB);
  assert.equal(runtime.initialize(tcr.references), tcr.names.D ? 3 : 2);
  const tcrResult = runtime.annotate(`sequence_id\tsequence\ntrb_case\t${tcr.sequence}\n`, 3);
  assert.equal(tcrResult.rows[0].locus, "TRB");
  assert.ok(tcrResult.rows[0].v_call);
  assert.ok(tcrResult.rows[0].j_call);
  if (tcr.names.D) assert.ok(tcrResult.rows[0].d_call);
});
