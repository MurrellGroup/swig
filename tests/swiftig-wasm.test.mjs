import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import { WASI } from "@bjorn3/browser_wasi_shim";
import { preprocessGermlineFasta } from "../src/germline-preprocess.ts";
import { inferLineageGermline } from "../src/lineage-alignment.ts";

const pack = JSON.parse(
  zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz", import.meta.url))),
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

  function setCallingProfile(profile) {
    return exports.swig_set_calling_profile(profile === "igblast_compatible" ? 1 : profile === "igblast_balanced" ? 2 : profile === "truth_optimized" ? 0 : Number(profile));
  }

  function setAssignerStrategy(strategy) {
    const value = strategy === "standard" ? 0 : strategy === "riat_mp" ? 1 : strategy === "aer" ? 2 : Number(strategy);
    return exports.swig_set_assigner_strategy(value);
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

  function annotateDoubleD(query, format, options, strand = 0) {
    const [pointer, size] = put(query);
    const count = exports.swig_annotate_double_d(
      pointer,
      size,
      format,
      600,
      strand,
      options.mode,
      options.minimumVjSpan,
      options.seedLength,
      options.pseudoTrim,
      options.maximumPseudoMismatches,
      options.minimumScoreGain,
    );
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
    const doubleDTsv = read(exports.swig_double_d_result_ptr(), exports.swig_double_d_result_len());
    const doubleDLines = doubleDTsv.trimEnd().split("\n");
    const doubleDHeaders = doubleDLines.shift().split("\t");
    const doubleDRows = doubleDLines.filter(Boolean).map((line) => {
      const values = line.split("\t");
      return Object.fromEntries(doubleDHeaders.map((header, index) => [header, values[index] ?? ""]));
    });
    return {
      count,
      headers,
      rows,
      tsv,
      doubleDTsv,
      doubleDHeaders,
      doubleDRows,
      doubleDCount: exports.swig_double_d_count(),
    };
  }

  return { initialize, annotate, annotateDoubleD, setCallingProfile, setAssignerStrategy };
}

function referenceFor(locus, customJName) {
  const v = locus.V.find((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0)) ?? locus.V[0];
  const d = locus.D?.[0];
  const originalJ = locus.J.find((allele) => allele[2]?.[0] >= 0 && allele[2]?.[1] >= 0) ?? locus.J[0];
  const j = customJName ? [customJName, originalJ[1], originalJ[2]] : originalJ;
  const c = locus.C?.[0];
  return {
    references: {
      V: asFasta([v]),
      D: d ? asFasta([d]) : "",
      J: asFasta([j]),
      C: c ? asFasta([c]) : "",
    },
    sequence: `${v[1]}AACCGG${d?.[1] ?? ""}TTG${j[1]}${c?.[1] ?? ""}`,
    names: { V: v[0], D: d?.[0] ?? "", J: j[0], C: c?.[0] ?? "" },
  };
}

function referenceCount(reference) {
  return Object.values(reference).reduce((count, fasta) => count + (fasta.match(/^>/gm) ?? []).length, 0);
}

function fastaRecords(text) {
  const records = [];
  let name = "";
  let sequence = "";
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (name) records.push({ name, sequence });
      name = line.slice(1).split(/\s+/, 1)[0];
      sequence = "";
    } else sequence += line.trim();
  }
  if (name) records.push({ name, sequence });
  return records;
}

function imgtTiers(speciesName, segment) {
  const selected = pack.species.find((entry) => entry.name === speciesName);
  const baseTaxon = speciesName.split("_", 1)[0];
  const genus = baseTaxon.split(" ", 1)[0];
  return [
    [selected],
    pack.species.filter((entry) => entry.name !== speciesName && entry.name.split("_", 1)[0] === baseTaxon),
    pack.species.filter((entry) => entry.name.split("_", 1)[0] !== baseTaxon && entry.name.startsWith(`${genus} `)),
    pack.species.filter((entry) => !entry.name.startsWith(`${genus} `)),
  ].map((group) => group.flatMap((entry) => entry?.loci.IGH?.[segment] ?? [])).filter((group) => group.length);
}

function preprocessTiered(text, speciesName, segment) {
  let report;
  for (const templates of imgtTiers(speciesName, segment)) {
    report = preprocessGermlineFasta(report?.fasta ?? text, segment, templates, ["IGH"]);
    if (report.annotated === report.count) break;
  }
  assert.ok(report);
  assert.equal(report.annotated, report.count, `${speciesName} ${segment}`);
  return report.fasta;
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

test("calling profiles are explicit, switchable, and reject unknown profile identifiers", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V[0];
  const j = human.loci.IGH.J[0];
  const d = ["IGHD_PROFILE_TEST*01", "ACGTA"];
  const runtime = await makeRuntime();
  runtime.initialize({ V: asFasta([v]), D: asFasta([d]), J: asFasta([j]), C: "" });
  const query = `>profile_case\n${v[1]}AAA${d[1]}TTT${j[1]}\n`;

  assert.equal(runtime.setCallingProfile("truth_optimized"), 0);
  const defaultResult = runtime.annotate(query, 1);
  assert.equal(runtime.setCallingProfile("igblast_compatible"), 0);
  const compatibilityResult = runtime.annotate(query, 1);
  assert.equal(runtime.setCallingProfile(2), -1);
  assert.equal(runtime.setCallingProfile("truth_optimized"), 0);
  assert.equal(runtime.annotate(query, 1).tsv, defaultResult.tsv, "resetting the profile did not restore default calls");
  assert.equal(defaultResult.rows[0].d_call, "");
  assert.equal(compatibilityResult.rows[0].d_call, d[0]);
});

test("terminal SHM does not masquerade as missing V 5-prime or J 3-prime sequence", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V.find((allele) => allele[1].length >= 120) ?? human.loci.IGH.V[0];
  const j = human.loci.IGH.J.find((allele) => allele[1].length >= 35) ?? human.loci.IGH.J[0];
  const mutate = (base) => ({ A: "C", C: "A", G: "T", T: "G" })[base] ?? "A";
  const terminalRun = 9;
  const mutatedV = [...v[1]].map((base, index) => index < terminalRun ? mutate(base) : base).join("");
  const mutatedJ = [...j[1]].map((base, index, values) => index >= values.length - terminalRun ? mutate(base) : base).join("");
  const query = `${mutatedV}AACCGG${mutatedJ}`;
  const runtime = await makeRuntime();
  runtime.initialize({ V: asFasta([v]), D: "", J: asFasta([j]), C: "" });
  const row = runtime.annotate(`>terminal_shm\n${query}\n`, 1).rows[0];
  assert.equal(row.v_sequence_start, "1");
  assert.equal(row.v_germline_start, "1");
  assert.equal(Number(row.j_sequence_end), query.length);
  assert.equal(Number(row.j_germline_end), j[1].length);
});

test("standard, RIAT-MP, and AER assignment strategies initialize explicitly", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V.slice(0, 8);
  const j = human.loci.IGH.J.slice(0, 3);
  const query = `>strategy_case\n${v[0][1]}AACCGG${j[0][1]}\n`;
  const observed = new Map();
  for (const strategy of ["standard", "riat_mp", "aer"]) {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize({ V: asFasta(v), D: "", J: asFasta(j), C: "" });
    const result = runtime.annotate(query, 1);
    assert.equal(result.count, 1);
    assert.ok(result.rows[0].v_call, `${strategy} did not assign V`);
    assert.ok(result.rows[0].j_call, `${strategy} did not assign J`);
    observed.set(strategy, [result.rows[0].v_call, result.rows[0].j_call]);
  }
  assert.equal((await makeRuntime()).setAssignerStrategy(3), -1);
  assert.deepEqual(observed.get("aer"), observed.get("standard"));
});

test("WASM annotates FASTA, FASTQ, and AIRR; handles heavy, light, TCR, strand, and J-only swaps", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human, "human reference set is missing");
  const runtime = await makeRuntime();

  const customJName = "IGHJ_SWIGTEST*01";
  const heavy = referenceFor(human.loci.IGH, customJName);
  assert.equal(runtime.initialize(heavy.references), referenceCount(heavy.references));

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
  assert.equal(fasta.rows[0].c_call, heavy.names.C);
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
  assert.equal(runtime.initialize(light.references), referenceCount(light.references));
  const lightResult = runtime.annotate(`>light_case\n${light.sequence}\n`, 1);
  assert.equal(lightResult.rows[0].locus, "IGK");
  assert.equal(lightResult.rows[0].d_call, "");
  assert.ok(lightResult.rows[0].v_call);
  assert.ok(lightResult.rows[0].j_call);
  assert.equal(lightResult.rows[0].c_call, light.names.C);

  const tcr = referenceFor(human.loci.TRB);
  assert.equal(runtime.initialize(tcr.references), referenceCount(tcr.references));
  const tcrResult = runtime.annotate(`sequence_id\tsequence\ntrb_case\t${tcr.sequence}\n`, 3);
  assert.equal(tcrResult.rows[0].locus, "TRB");
  assert.ok(tcrResult.rows[0].v_call);
  assert.ok(tcrResult.rows[0].j_call);
  assert.equal(tcrResult.rows[0].c_call, tcr.names.C);
  if (tcr.names.D) assert.ok(tcrResult.rows[0].d_call);
});

test("bundled macaque KIMDB references produce complete IMGT region and CDR calls end to end", { timeout: 45_000 }, async () => {
  const runtime = await makeRuntime();
  for (const [speciesName, slug] of [["Macaca mulatta_AG07107", "Macaca_mulatta"], ["Macaca fascicularis", "Macaca_fascicularis"]]) {
    const species = pack.species.find((entry) => entry.name === speciesName);
    assert.ok(species?.loci.IGH);
    const root = new URL(`../public/references/kimdb-1.1/${slug}/IGH/`, import.meta.url);
    const V = preprocessTiered(fs.readFileSync(new URL("V.fasta", root), "utf8"), speciesName, "V");
    const J = preprocessTiered(fs.readFileSync(new URL("J.fasta", root), "utf8"), speciesName, "J");
    const D = preprocessGermlineFasta(fs.readFileSync(new URL("D.fasta", root), "utf8"), "D", [], ["IGH"]).fasta;
    const C = asFasta(species.loci.IGH.C ?? []);
    const v = fastaRecords(V)[0];
    const d = fastaRecords(D)[0];
    const j = fastaRecords(J)[0];
    const c = fastaRecords(C)[0];
    assert.ok(v && d && j && c);
    const references = { V, D, J, C };
    assert.equal(runtime.initialize(references), referenceCount(references));
    const result = runtime.annotate(`>${slug}_kimdb\n${v.sequence}AACCGG${d.sequence}TTG${j.sequence}${c.sequence}\n`, 1);
    const row = result.rows[0];
    assert.equal(row.region_definition, "IMGT");
    assert.equal(row.v_annotation_source, "IMGT-boundary-transfer");
    assert.equal(row.j_annotation_source, "J-anchor-transfer");
    for (const region of ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3"]) {
      assert.ok(row[region], `${speciesName} ${region} is empty`);
      assert.ok(row[`${region}_aa`], `${speciesName} ${region}_aa is empty`);
    }
    assert.ok(row.c_call, `${speciesName} constant call is empty`);
  }
});

test("opt-in double-D screening is sparse and leaves the standard AIRR result byte-for-byte unchanged", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V.find((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0));
  const j = human.loci.IGH.J.find((allele) => allele[2]?.[0] >= 0 && allele[2]?.[1] >= 0);
  assert.ok(v && j);
  const d1 = ["IGHDTEST1*01", "ATGCCGTACGTTAGC"];
  const d2 = ["IGHDTEST2*01", "CGATTCGGAACCTGA"];
  const references = {
    V: asFasta([v]),
    D: asFasta([d1, d2]),
    J: asFasta([j]),
    C: "",
  };
  const sequence = `${v[1]}AAAA${d1[1]}GGGGGG${d2[1]}TTTT${j[1]}`;
  const query = `>vddj_case\n${sequence}\n`;
  const runtime = await makeRuntime();
  runtime.initialize(references);
  const baseline = runtime.annotate(query, 1);
  const options = {
    mode: 1,
    minimumVjSpan: 0,
    seedLength: 11,
    pseudoTrim: 5,
    maximumPseudoMismatches: 3,
    minimumScoreGain: 8,
  };
  const screened = runtime.annotateDoubleD(query, 1, options);
  assert.equal(screened.tsv, baseline.tsv, "the opt-in sidecar changed the ordinary AIRR table");
  assert.equal(screened.doubleDCount, 1);
  assert.equal(screened.doubleDRows.length, 1);
  assert.equal(screened.doubleDRows[0].swig_batch_record_index, "0");
  assert.equal(screened.doubleDRows[0].sequence_id, "vddj_case");
  assert.equal(screened.doubleDRows[0].d_call, d1[0]);
  assert.equal(screened.doubleDRows[0].d2_call, d2[0]);
  assert.ok(Number(screened.doubleDRows[0].d_sequence_end) < Number(screened.doubleDRows[0].d2_sequence_start));
  assert.equal(screened.doubleDRows[0].np2, "GGGGGG");
  assert.ok(Number(screened.doubleDRows[0].swig_double_d_score_gain) >= options.minimumScoreGain);

  const values = { ...screened.rows[0], ...screened.doubleDRows[0] };
  const inferred = inferLineageGermline([{
    record: {
      ordinal: 0,
      sequenceId: values.sequence_id,
      locus: values.locus,
      vIdentity: Number(values.v_identity),
      jIdentity: Number(values.j_identity),
    },
    values,
  }]);
  assert.equal(inferred.doubleDTemplate, true);
  assert.equal(inferred.selectedDCall, d1[0]);
  assert.equal(inferred.selectedD2Call, d2[0]);

  const leftPadding = Math.max(0, Number(values.v_germline_start || 1) - 1);
  const coordinateToColumn = new Map();
  let coordinate = Number(values.v_sequence_start);
  for (let column = 0; column < values.sequence_alignment.length; column += 1) {
    if (values.sequence_alignment[column] === "-") continue;
    coordinateToColumn.set(coordinate, leftPadding + column);
    coordinate += 1;
  }
  for (const prefix of ["d", "d2"]) {
    let queryCoordinate = Number(values[`${prefix}_sequence_start`]);
    for (let column = 0; column < values[`${prefix}_sequence_alignment`].length; column += 1) {
      if (values[`${prefix}_sequence_alignment`][column] === "-") continue;
      assert.equal(
        inferred.template[coordinateToColumn.get(queryCoordinate)],
        values[`${prefix}_germline_alignment`][column],
        `${prefix.toUpperCase()} germline base was not retained in the lineage root`,
      );
      queryCoordinate += 1;
    }
  }
  for (let position = Number(values.d_sequence_end) + 1; position < Number(values.d2_sequence_start); position += 1) {
    assert.equal(inferred.template[coordinateToColumn.get(position)], "N", "NP2 must remain unknown in the tree root");
  }
  for (let position = Number(values.d2_sequence_end) + 1; position < Number(values.j_sequence_start); position += 1) {
    assert.equal(inferred.template[coordinateToColumn.get(position)], "N", "NP3 must remain unknown in the tree root");
  }

  const gated = runtime.annotateDoubleD(query, 1, { ...options, mode: 2, minimumVjSpan: 10_000 });
  assert.equal(gated.tsv, baseline.tsv);
  assert.equal(gated.doubleDCount, 0);
  assert.equal(gated.doubleDRows.length, 0);
});
