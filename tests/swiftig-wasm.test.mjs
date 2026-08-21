import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import { WASI } from "@bjorn3/browser_wasi_shim";
import { preprocessGermlineFasta } from "../src/germline-preprocess.ts";
import { inferLineageGermline } from "../src/lineage-alignment.ts";
import {
  generateVdjDataset,
  parseFasta as parseSimulatorFasta,
  recordsToFasta,
} from "./vdj-simulator.mjs";

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
    return exports.swig_set_calling_profile(profile === "igblast_compatible" || profile === "igblast_balanced" ? 1 : profile === "r_optimized" ? 2 : profile === "truth_optimized" ? 0 : Number(profile));
  }

  function setAssignerStrategy(strategy) {
    const value = strategy === "standard" ? 0 : strategy === "riat_mp" ? 1 :
      strategy === "aer" ? 2 : strategy === "aer_robust" ? 3 : Number(strategy);
    return exports.swig_set_assigner_strategy(value);
  }

  function setOptimizedKernels(enabled) {
    return exports.swig_set_optimized_kernels(enabled ? 1 : 0);
  }

  function setOptimizedOutput(enabled) {
    return exports.swig_set_optimized_output(enabled ? 1 : 0);
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

  return {
    initialize,
    annotate,
    annotateDoubleD,
    setCallingProfile,
    setAssignerStrategy,
    setOptimizedKernels,
    setOptimizedOutput,
  };
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
  assert.equal(pack.license?.id, "CC-BY-4.0");
  assert.equal(pack.license?.modified, true);
  assert.match(pack.license?.attribution ?? "", /IMGT®.*CNRS/);
  assert.match(pack.license?.modifications ?? "", /ungapped nucleotide sequences/);
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
  assert.equal(runtime.setCallingProfile("r_optimized"), 0);
  assert.equal(runtime.setCallingProfile(3), -1);
  assert.equal(runtime.setCallingProfile("truth_optimized"), 0);
  assert.equal(runtime.annotate(query, 1).tsv, defaultResult.tsv, "resetting the profile did not restore default calls");
  assert.equal(defaultResult.rows[0].d_call, "");
  assert.equal(compatibilityResult.rows[0].d_call, d[0]);
});

test("joint boundary scoring recovers a strong KIMDB D tract swallowed by a gappy V tail", async () => {
  const referenceRoot = new URL("../public/references/kimdb-1.1/Macaca_mulatta/IGH/", import.meta.url);
  const references = {
    V: fs.readFileSync(new URL("V.fasta", referenceRoot), "utf8"),
    D: fs.readFileSync(new URL("D.fasta", referenceRoot), "utf8"),
    J: fs.readFileSync(new URL("J.fasta", referenceRoot), "utf8"),
    C: "",
  };
  const query = fs.readFileSync(
    new URL("fixtures/kimdb-perfect-d-regression.fasta", import.meta.url),
    "utf8",
  );

  for (const strategy of ["aer", "aer_robust", "riat_mp", "standard"]) {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize(references);
    for (const profile of ["truth_optimized", "igblast_compatible"]) {
      assert.equal(runtime.setCallingProfile(profile), 0);
      const row = runtime.annotate(query, 1, 1).rows[0];
      assert.ok(row.d_call, `${strategy}/${profile} omitted the strong D call`);
      if (profile === "igblast_compatible") {
        assert.match(row.d_call, /IGHD5-32\*01_S0263/, `${strategy} omitted the simulated D allele`);
      }
      assert.ok(
        row.d_sequence_alignment.includes("ATACAGTGGGTACAGTT"),
        `${strategy}/${profile} did not retain the complete perfect D tract`,
      );
      assert.equal(row.d_sequence_alignment, row.d_germline_alignment);
      assert.ok(Number(row.v_sequence_end) < Number(row.d_sequence_start));
      assert.ok(Number(row.d_sequence_start) <= 421 && Number(row.d_sequence_end) >= 437);
    }
  }
  const optimized = await makeRuntime();
  assert.equal(optimized.setAssignerStrategy("aer_robust"), 0);
  assert.equal(optimized.setCallingProfile("r_optimized"), 0);
  optimized.initialize(references);
  const optimizedRow = optimized.annotate(query, 1, 1).rows[0];
  assert.match(optimizedRow.d_call, /IGHD5-32\*01_S0263/);
  assert.ok(optimizedRow.d_sequence_alignment.includes("ATACAGTGGGTACAGTT"));
});

test("R-optimized corrects representative junction-facing V over-extension regressions", async () => {
  const roots = ["Macaca_mulatta", "Macaca_fascicularis"].map((species) =>
    new URL(`../public/references/kimdb-1.1/${species}/IGH/`, import.meta.url));
  const combined = (segment) => {
    const records = new Map();
    for (const root of roots) {
      for (const record of parseSimulatorFasta(fs.readFileSync(new URL(`${segment}.fasta`, root), "utf8"))) {
        records.set(record.name, record.sequence);
      }
    }
    return [...records].map(([name, sequence]) => `>${name}\n${sequence}\n`).join("");
  };
  const references = { V: combined("V"), D: combined("D"), J: combined("J"), C: "" };
  const records = parseSimulatorFasta(fs.readFileSync(
    new URL("fixtures/r-optimized-v-boundary-regressions.fasta", import.meta.url), "utf8"));
  const query = records.map((record) => `>${record.name}\n${record.sequence}\n`).join("");
  const annotate = async (profile) => {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy("aer_robust"), 0);
    assert.equal(runtime.setCallingProfile(profile), 0);
    runtime.initialize(references);
    return runtime.annotate(query, 1, 1).rows;
  };
  const [legacy, optimized] = await Promise.all([
    annotate("truth_optimized"), annotate("r_optimized"),
  ]);
  const goldEnds = [330, 428, 349];
  const legacyLoss = legacy.reduce((sum, row, index) =>
    sum + Math.abs(Number(row.v_sequence_end) - goldEnds[index]), 0);
  const optimizedLoss = optimized.reduce((sum, row, index) =>
    sum + Math.abs(Number(row.v_sequence_end) - goldEnds[index]), 0);
  assert.ok(optimizedLoss < legacyLoss);
  assert.deepEqual(optimized.map((row) => Number(row.v_sequence_end)), [330, 428, 352]);
  for (let index = 0; index < optimized.length; index += 1) {
    assert.ok(optimized[index].v_call.split(",").includes([
      "IGHV4-NL_36*01_S9644", "IGHV3-122*01", "IGHV3-172*01_S0577",
    ][index]));
    assert.ok(Math.abs(Number(optimized[index].v_sequence_end) - goldEnds[index]) <= 3);
  }
});

test("R-optimized conditions D presence on score and independent template support", async () => {
  const roots = ["Macaca_mulatta", "Macaca_fascicularis"].map((species) =>
    new URL(`../public/references/kimdb-1.1/${species}/IGH/`, import.meta.url));
  const combined = (segment) => {
    const records = new Map();
    for (const root of roots) {
      for (const record of parseSimulatorFasta(fs.readFileSync(new URL(`${segment}.fasta`, root), "utf8"))) {
        records.set(record.name, record.sequence);
      }
    }
    return [...records].map(([name, sequence]) => `>${name}\n${sequence}\n`).join("");
  };
  const runtime = await makeRuntime();
  assert.equal(runtime.setAssignerStrategy("aer_robust"), 0);
  assert.equal(runtime.setCallingProfile("r_optimized"), 0);
  runtime.initialize({ V: combined("V"), D: combined("D"), J: combined("J"), C: "" });
  const query = fs.readFileSync(
    new URL("fixtures/r-optimized-d-decision-regressions.fasta", import.meta.url), "utf8");
  const [supported, zeroD] = runtime.annotate(query, 1, 1).rows;

  assert.match(supported.d_call, /IGHD1-5\*01/);
  assert.equal(supported.d_sequence_alignment, supported.d_germline_alignment);
  assert.equal(zeroD.d_call, "");
  assert.match(zeroD.j_call, /IGHJ3-2\*01/);
});

test("AER-R accepts strong distributed D evidence instead of a short exact-seed decoy", async () => {
  const roots = ["Macaca_mulatta", "Macaca_fascicularis"].map((species) =>
    new URL(`../public/references/kimdb-1.1/${species}/IGH/`, import.meta.url));
  const combined = (segment) => {
    const records = new Map();
    for (const root of roots) {
      for (const record of parseSimulatorFasta(fs.readFileSync(new URL(`${segment}.fasta`, root), "utf8"))) {
        assert.ok(!records.has(record.name) || records.get(record.name) === record.sequence);
        records.set(record.name, record.sequence);
      }
    }
    return [...records].map(([name, sequence]) => `>${name}\n${sequence}\n`).join("");
  };
  const references = { V: combined("V"), D: combined("D"), J: combined("J"), C: "" };
  const sequence = "GGAGCCCGGAACACTGAGGGTCGTGGCTAAAGTCGTTTGGTGCAGGGTGCCCTGCTTCTCAGGGATGCTGTACACGTTCCTAACAGGTGCAGCTTCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCTGTCACCTGCGCTGTCTCTGGTGGCTCTATCAGCAGTAGCTACTGGAATTGGATCCGCCAGCCCCCAGCGAAGGGACGGGAGTGTATTGGCTACATCTATGGTAGTGGTGGGAGCACCAGCGACAACCCCTCCCTCCTGAGTCGAGTCACCCTGTCAGTAGATACGTCTGAGCACCAGCTCTCCTTGAAGCTGAGCCCTGTGTCCGCCGCGGACCCGGTCGTGTATTACTGTGCGATAGGATTAACGGTCTAGAAGCATCTACATCCTCAATCTGATATCTGGGGCCCTGGCACCCCAATCACCATGTCCTCAGGCTCCTATGTTGGGCAGGCTGGTCCGAGGAGCCTCTCATCAAAACTCAGGATTGGCACCTTGCGCGCCGCCAGAGGTAATAACCACACTTTCAGCAAGGTTAATCGACCGCGTCTGATTACATTT";
  const query = `>distributed_d\n${sequence}\n`;

  for (const profile of ["truth_optimized", "igblast_compatible"]) {
    const ordinary = await makeRuntime();
    assert.equal(ordinary.setAssignerStrategy("aer"), 0);
    assert.equal(ordinary.setCallingProfile(profile), 0);
    ordinary.initialize(references);
    const ordinaryRow = ordinary.annotate(query, 1).rows[0];

    const robust = await makeRuntime();
    assert.equal(robust.setAssignerStrategy("aer_robust"), 0);
    assert.equal(robust.setCallingProfile(profile), 0);
    robust.initialize(references);
    const row = robust.annotate(query, 1).rows[0];

    assert.doesNotMatch(ordinaryRow.d_call, /IGHD6-39\*01/);
    assert.equal(row.v_call, "IGHV4-NL_1*02_S8056");
    assert.equal(row.d_call, "IGHD6-39*01");
    assert.equal(row.j_call, "IGHJ2*01_S5087");
    assert.equal(row.d_sequence_start, "386");
    assert.equal(row.d_sequence_end, "402");
    assert.equal(row.d_germline_start, "2");
    assert.equal(row.d_germline_end, "18");
    assert.equal(row.d_sequence_alignment, "GGTCTAGAAGCATCTAC");
    assert.equal(row.d_germline_alignment, "GGTATAGCAGCAGCTAC");
    assert.equal(row.np1, "GATTAAC");
    assert.equal(row.np2, "ATCCTCAATCT");
  }
  const optimized = await makeRuntime();
  assert.equal(optimized.setAssignerStrategy("aer_robust"), 0);
  assert.equal(optimized.setCallingProfile("r_optimized"), 0);
  optimized.initialize(references);
  assert.match(optimized.annotate(query, 1).rows[0].d_call, /IGHD6-39\*01/);
});

test("AER-R short D calls cannot crowd a stronger distributed match out of traceback", async () => {
  const roots = ["Macaca_mulatta", "Macaca_fascicularis"].map((species) =>
    new URL(`../public/references/kimdb-1.1/${species}/IGH/`, import.meta.url));
  const combined = (segment) => {
    const records = new Map();
    for (const root of roots) {
      for (const record of parseSimulatorFasta(fs.readFileSync(new URL(`${segment}.fasta`, root), "utf8"))) {
        records.set(record.name, record.sequence);
      }
    }
    return [...records].map(([name, sequence]) => `>${name}\n${sequence}\n`).join("");
  };
  const sequence = "GGAGCCCGGAACACTGAGGGTCGTGGCTAAAGTCGTTTGGTGCAGGGTGCCCTGCTTCTCAGGGATGCTGTACACGTTCCTAACAGGTGCAGCTTCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCTGTCACCTGCGCTGTCTCTGGTGGCTCTATCAGCAGTAGCTACTGGAATTGGATCCGCCAGCCCCCAGCGAAGGGACGGGAGTGTATTGGCTACATCTATGGTAGTGGTGGGAGCACCAGCGACAACCCCTCCCTCCTGAGTCGAGTCACCCTGTCAGTAGATACGTCTGAGCACCAGCTCTCCTTGAAGCTGAGCCCTGTGTCCGCCGCGGACCCGGTCGTGTATTACTGTGCGATAGGATTAACGGTCTAGAAGCATCTACATCCTCAATCTGATATCTGGGGCCCTGGCACCCCAATCACCATGTCCTCAGGCTCCTATGTTGGGCAGGCTGGTCCGAGGAGCCTCTCATCAAAACTCAGGATTGGCACCTTGCGCGCCGCCAGAGGTAATAACCACACTTTCAGCAAGGTTAATCGACCGCGTCTGATTACATTT";
  // Fifty exact 7-mers have stronger seed rankings than IGHD6-39's mutated
  // 17-nt tract but lower affine scores. They deliberately overflow the
  // ordinary bounded D candidate/traceback set while remaining inferior joint
  // explanations. A binary "no D" retry would never run because these decoys
  // all produce short calls.
  const decoys = Array.from({ length: 50 }, (_, index) => {
    const start = 376 + (index % 35);
    return `>DECOY_D_${String(index).padStart(2, "0")}\n${sequence.slice(start, start + 7)}\n`;
  }).join("");
  for (const profile of ["truth_optimized", "igblast_compatible"]) {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy("aer_robust"), 0);
    assert.equal(runtime.setCallingProfile(profile), 0);
    runtime.initialize({
      V: combined("V"),
      D: decoys + combined("D"),
      J: combined("J"),
      C: "",
    });
    const row = runtime.annotate(`>crowded_distributed_d\n${sequence}\n`, 1).rows[0];
    assert.equal(row.d_call, "IGHD6-39*01");
    assert.equal(row.d_sequence_alignment, "GGTCTAGAAGCATCTAC");
    assert.equal(row.d_germline_alignment, "GGTATAGCAGCAGCTAC");
  }
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

test("standard, RIAT-MP, AER, and experimental AER-R strategies initialize explicitly", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V.slice(0, 8);
  const j = human.loci.IGH.J.slice(0, 3);
  const query = `>strategy_case\n${v[0][1]}AACCGG${j[0][1]}\n`;
  const observed = new Map();
  for (const strategy of ["standard", "riat_mp", "aer", "aer_robust"]) {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize({ V: asFasta(v), D: "", J: asFasta(j), C: "" });
    const result = runtime.annotate(query, 1);
    assert.equal(result.count, 1);
    assert.ok(result.rows[0].v_call, `${strategy} did not assign V`);
    assert.ok(result.rows[0].j_call, `${strategy} did not assign J`);
    observed.set(strategy, [result.rows[0].v_call, result.rows[0].j_call]);
  }
  assert.equal((await makeRuntime()).setAssignerStrategy(4), -1);
  assert.deepEqual(observed.get("aer"), observed.get("standard"));
});

test("development simulator finds net AER-R boundary gains without changing ordinary AER", { timeout: 30_000 }, async () => {
  const referenceRoot = new URL("../public/references/kimdb-1.1/Macaca_mulatta/IGH/", import.meta.url);
  const references = {
    V: fs.readFileSync(new URL("V.fasta", referenceRoot), "utf8"),
    D: fs.readFileSync(new URL("D.fasta", referenceRoot), "utf8"),
    J: fs.readFileSync(new URL("J.fasta", referenceRoot), "utf8"),
    C: "",
  };
  const simulated = generateVdjDataset({
    V: parseSimulatorFasta(references.V),
    D: parseSimulatorFasta(references.D),
    J: parseSimulatorFasta(references.J),
    count: 400,
    seed: 913_771,
    doubleDRate: 0,
    reverseRate: 0,
  });
  const query = recordsToFasta(simulated);
  const annotateWith = async (strategy) => {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize(references);
    return runtime.annotate(query, 1).rows;
  };
  const [aer, robust] = await Promise.all([annotateWith("aer"), annotateWith("aer_robust")]);
  const atoms = (value) => String(value ?? "").split(/[,/]/).filter(Boolean);
  let aerStrongD = 0;
  let robustStrongD = 0;
  let recovered = 0;
  let aerVj = 0;
  let robustVj = 0;
  for (let index = 0; index < simulated.length; index += 1) {
    const truth = simulated[index].truth;
    const a = aer[index];
    const r = robust[index];
    const strongD = (truth.spans.D1?.bases ?? 0) >= 10;
    const aD = atoms(a.d_call).includes(truth.dCalls[0]);
    const rD = atoms(r.d_call).includes(truth.dCalls[0]);
    if (strongD) {
      aerStrongD += Number(aD);
      robustStrongD += Number(rD);
      recovered += Number(!aD && rD);
    }
    aerVj += Number(atoms(a.v_call).includes(truth.vCall) && atoms(a.j_call).includes(truth.jCall));
    robustVj += Number(atoms(r.v_call).includes(truth.vCall) && atoms(r.j_call).includes(truth.jCall));
  }
  assert.ok(recovered > 0, "the stress cohort did not exercise a recoverable AER boundary failure");
  assert.ok(robustStrongD >= aerStrongD, "AER-R reduced strong-D truth recall on the fixed stress cohort");
  assert.ok(robustVj >= aerVj, "AER-R reduced joint V/J truth recall on the fixed stress cohort");
});

test("AER-R's gated complete-J retry recovers strong-V seed-pruning failures", { timeout: 30_000 }, async () => {
  const referenceRoot = new URL("../public/references/kimdb-1.1/Macaca_mulatta/IGH/", import.meta.url);
  const references = {
    V: fs.readFileSync(new URL("V.fasta", referenceRoot), "utf8"),
    D: fs.readFileSync(new URL("D.fasta", referenceRoot), "utf8"),
    J: fs.readFileSync(new URL("J.fasta", referenceRoot), "utf8"),
    C: "",
  };
  const simulated = generateVdjDataset({
    V: parseSimulatorFasta(references.V),
    D: parseSimulatorFasta(references.D),
    J: parseSimulatorFasta(references.J),
    count: 320,
    seed: 913_771,
  });
  const selected = simulated.filter((record) =>
    record.id === "sim_0000292" || record.id === "sim_0000319");
  assert.equal(selected.length, 2);
  const query = recordsToFasta(selected);
  const annotateWith = async (strategy, profile = "truth_optimized") => {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    assert.equal(runtime.setCallingProfile(profile), 0);
    runtime.initialize(references);
    return runtime.annotate(query, 1).rows;
  };
  const [aer, robust, optimized] = await Promise.all([
    annotateWith("aer"), annotateWith("aer_robust"), annotateWith("aer_robust", "r_optimized"),
  ]);
  const atoms = (value) => String(value ?? "").split(/[,/]/).filter(Boolean);
  let ordinaryMisses = 0;
  for (let index = 0; index < selected.length; index += 1) {
    const truth = selected[index].truth;
    const a = aer[index];
    const r = robust[index];
    ordinaryMisses += Number(
      !atoms(a.j_call).includes(truth.jCall) || !atoms(a.d_call).includes(truth.dCalls[0]));
    assert.ok(atoms(r.v_call).includes(truth.vCall));
    assert.ok(atoms(r.j_call).includes(truth.jCall));
    assert.ok(atoms(r.d_call).includes(truth.dCalls[0]));
    assert.ok(atoms(optimized[index].v_call).includes(truth.vCall));
    assert.ok(atoms(optimized[index].j_call).includes(truth.jCall));
  }
  assert.ok(ordinaryMisses > 0, "the fixed cohort no longer exercises the gated rescue");
});

test("AER-R double-D rescue improves mutation-tolerant tandem detection only when selected", { timeout: 30_000 }, async () => {
  const referenceRoot = new URL("../public/references/kimdb-1.1/Macaca_mulatta/IGH/", import.meta.url);
  const references = {
    V: fs.readFileSync(new URL("V.fasta", referenceRoot), "utf8"),
    D: fs.readFileSync(new URL("D.fasta", referenceRoot), "utf8"),
    J: fs.readFileSync(new URL("J.fasta", referenceRoot), "utf8"),
    C: "",
  };
  const simulated = generateVdjDataset({
    V: parseSimulatorFasta(references.V),
    D: parseSimulatorFasta(references.D),
    J: parseSimulatorFasta(references.J),
    count: 150,
    seed: 404,
    doubleDRate: 0.20,
    reverseRate: 0,
  });
  const query = recordsToFasta(simulated);
  const options = {
    mode: 1, minimumVjSpan: 40, seedLength: 11,
    pseudoTrim: 5, maximumPseudoMismatches: 3, minimumScoreGain: 8,
  };
  const screenWith = async (strategy) => {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize(references);
    return new Set(runtime.annotateDoubleD(query, 1, options).doubleDRows.map((row) => row.sequence_id));
  };
  const [aer, robust] = await Promise.all([screenWith("aer"), screenWith("aer_robust")]);
  const tandem = simulated.filter((record) => record.truth.tandem);
  const baselineDetected = tandem.filter((record) => aer.has(record.id)).length;
  const robustDetected = tandem.filter((record) => robust.has(record.id)).length;
  const single = simulated.filter((record) => !record.truth.tandem);
  const addedFalsePositives = single.filter((record) => robust.has(record.id) && !aer.has(record.id)).length;
  assert.ok(robustDetected > baselineDetected, "the robust seed tier did not rescue tandem-D simulations");
  assert.ok(addedFalsePositives <= Math.ceil(single.length * 0.01), "the robust seed tier admitted too many single-D screens");
});

test("optimized RIAT-MP/AER kernels and AIRR writer are byte-identical to retained reference implementations", { timeout: 45_000 }, async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const references = {
    V: asFasta(human.loci.IGH.V.slice(0, 48)),
    D: asFasta(human.loci.IGH.D),
    J: asFasta(human.loci.IGH.J),
    C: asFasta(human.loci.IGH.C.slice(0, 4)),
  };
  const mutate = (base) => ({ A: "C", C: "G", G: "T", T: "A" })[base] ?? "N";
  const records = [];
  for (let index = 0; index < 64; index += 1) {
    const v = human.loci.IGH.V[index % 48][1];
    const d = human.loci.IGH.D[index % human.loci.IGH.D.length][1];
    const j = human.loci.IGH.J[index % human.loci.IGH.J.length][1];
    const c = human.loci.IGH.C[index % 4][1];
    let sequence = `${v}${"ACGTN".slice(0, index % 5)}${d}${"TGCA".slice(0, index % 4)}${j}${c}`;
    if (index % 6 === 1) {
      const at = Math.min(31, sequence.length - 1);
      sequence = `${sequence.slice(0, at)}${mutate(sequence[at])}${sequence.slice(at + 1)}`;
    } else if (index % 6 === 2) {
      sequence = `${sequence.slice(0, 27)}NNN${sequence.slice(30)}`;
    } else if (index % 6 === 3) {
      sequence = `${sequence.slice(0, 43)}ACT${sequence.slice(43)}`;
    } else if (index % 6 === 4) {
      sequence = `${sequence.slice(0, 55)}${sequence.slice(57)}`;
    } else if (index % 6 === 5) {
      sequence = reverseComplement(sequence);
    }
    records.push(`>equivalence_${index}\n${sequence}\n`);
  }
  const query = records.join("");

  for (const strategy of ["riat_mp", "aer", "aer_robust"]) {
    const runtime = await makeRuntime();
    assert.equal(runtime.setAssignerStrategy(strategy), 0);
    runtime.initialize(references);

    assert.equal(runtime.setOptimizedKernels(false), 0);
    assert.equal(runtime.setOptimizedOutput(false), 0);
    const reference = runtime.annotate(query, 1).tsv;

    assert.equal(runtime.setOptimizedKernels(true), 0);
    const optimizedKernels = runtime.annotate(query, 1).tsv;
    assert.equal(
      optimizedKernels,
      reference,
      `${strategy} optimized compute kernels changed AIRR bytes`,
    );

    assert.equal(runtime.setOptimizedOutput(true), 0);
    const optimizedOutput = runtime.annotate(query, 1).tsv;
    assert.equal(
      optimizedOutput,
      reference,
      `${strategy} direct AIRR writer changed AIRR bytes`,
    );
  }
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
  for (const segment of ["v", "d", "j", "c"]) {
    assert.ok(fasta.headers.includes(`${segment}_support`));
    if (fasta.rows[0][`${segment}_call`]) {
      assert.notEqual(fasta.rows[0][`${segment}_support`], "");
      const support = Number(fasta.rows[0][`${segment}_support`]);
      assert.ok(Number.isFinite(support) && support >= 0, `${segment.toUpperCase()} support is not a finite E-value`);
    } else {
      assert.equal(fasta.rows[0][`${segment}_support`], "");
    }
  }
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

test("calibrated AIRR support uses the supplied segment database search space", async () => {
  const human = pack.species.find((entry) => entry.name === "Homo sapiens");
  assert.ok(human?.loci.IGH);
  const v = human.loci.IGH.V[0];
  const j = human.loci.IGH.J[0];
  const query = `>search_space\n${v[1]}AACCGG${j[1]}\n`;
  const runtime = await makeRuntime();
  runtime.initialize({ V: asFasta([v]), D: "", J: asFasta([j]), C: "" });
  const one = runtime.annotate(query, 1).rows[0];
  const duplicate = [`${v[0]}_DECOY`, v[1], v[2]];
  runtime.initialize({ V: asFasta([v, duplicate]), D: "", J: asFasta([j]), C: "" });
  const two = runtime.annotate(query, 1).rows[0];
  assert.equal(one.v_score, two.v_score);
  assert.ok(Number(one.v_support) > 0);
  assert.ok(Number(two.v_support) > Number(one.v_support), "a doubled V search space did not increase the E-value");
  assert.equal(one.d_support, "");
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
