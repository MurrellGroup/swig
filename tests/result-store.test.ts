import assert from "node:assert/strict";
import test from "node:test";

import "fake-indexeddb/auto";
import { AirrResultStore, EMPTY_FILTERS } from "../src/result-store.ts";

const encoder = new TextEncoder();
const header = [
  "sequence_id", "sequence", "locus", "v_call", "d_call", "j_call", "productive",
  "cdr3", "cdr3_aa", "junction_aa", "v_identity", "d_identity", "j_identity",
  "vj_in_frame", "stop_codon", "complete_vdj", "rev_comp",
].join("\t");

function makeBody(start: number, count: number): Uint8Array {
  let body = "";
  for (let offset = 0; offset < count; offset += 1) {
    const ordinal = start + offset;
    const locus = ordinal % 4 < 2 ? "IGH" : ordinal % 4 === 2 ? "IGK" : "IGL";
    body += [
      `read_${ordinal}`,
      "ACGT".repeat(75),
      locus,
      `IGHV${ordinal % 8 + 1}*01`,
      locus === "IGH" ? `IGHD${ordinal % 4 + 1}*01` : "",
      `${locus}J${ordinal % 5 + 1}*01`,
      ordinal % 3 ? "T" : "F",
      "TGTGCCAGAGATCGT",
      "CARD R".replace(" ", ""),
      "CARDRWG",
      "0.98",
      locus === "IGH" ? "1" : "",
      "0.97",
      "T",
      "F",
      "T",
      "F",
    ].join("\t") + "\n";
  }
  return encoder.encode(body);
}

test("50k AIRR records are committed, filtered, retrieved, and exported batchwise", async () => {
  const store = new AirrResultStore();
  const batchSize = 1000;
  for (let start = 0; start < 50_000; start += batchSize) {
    await store.appendBatch(header, makeBody(start, batchSize));
  }
  await store.finalize();

  assert.equal(store.count, 50_000);
  assert.equal(store.summary.assigned, 50_000);
  assert.equal(store.summary.withCdr3, 50_000);
  assert.equal(store.facets().loci.reduce((sum, value) => sum + value.count, 0), 50_000);
  const repertoire = store.repertoire({ locus: "IGH", ambiguity: "fractional" });
  assert.equal(repertoire.records, 25_000);
  assert.equal(repertoire.vCalls.reduce((sum, value) => sum + value.count, 0), 25_000);
  assert.ok(repertoire.vjPairs.length > 0);

  const page = await store.page({ ...EMPTY_FILTERS, locus: "IGK" }, 0, 25);
  assert.equal(page.rows.length, 25);
  assert.ok(page.rows.every((row) => row.locus === "IGK"));
  const detail = await store.detail(page.rows[0]);
  assert.equal(detail.sequence_id, page.rows[0].sequenceId);
  assert.equal(detail.sequence.length, 300);

  let exportedBytes = 0;
  await store.writeAirr(async (part) => {
    exportedBytes += typeof part === "string" ? encoder.encode(part).byteLength : part instanceof Blob ? part.size : part.byteLength;
  });
  assert.equal(exportedBytes, store.outputBytes);
  await store.clear();
});

test("constant-region evidence yields isotype facets and repertoire summaries", async () => {
  const store = new AirrResultStore();
  const constantHeader = [
    "sequence_id", "sequence", "locus", "v_call", "d_call", "j_call", "c_call",
    "productive", "cdr3_aa", "v_identity", "d_identity", "j_identity", "c_identity",
    "c_sequence_alignment", "vj_in_frame", "stop_codon", "complete_vdj", "rev_comp",
  ].join("\t");
  const constantBody = [
    "ig_m\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHM*01\tT\tCARDR\t.98\t1\t.97\t.99\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_g_short\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHG1*01\tT\tCARDR\t.98\t1\t.97\t.99\tAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_a\tACGT\tIGH\tIGHV2*01\tIGHD1*01\tIGHJ2*01\tIGHA2*01\tF\tCARDRR\t.95\t1\t.96\t.97\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
  ].join("\n") + "\n";
  await store.appendBatch(constantHeader, encoder.encode(constantBody));
  await store.finalize();
  assert.deepEqual(store.facets().isotypes.map((item) => item.value).sort(), ["IgA2", "IgM"]);
  assert.equal(store.repertoire().isotypes.reduce((sum, item) => sum + item.count, 0), 2);
  const page = await store.page({ ...EMPTY_FILTERS, isotype: "IgM" }, 0, 10);
  assert.equal(page.rows.length, 1);
  assert.equal(page.rows[0].cCall, "IGHM*01");
  await store.clear();
});

test("lineage export adds AIRR clone_id values and leaves excluded records blank", async () => {
  const store = new AirrResultStore();
  await store.appendBatch(header, makeBody(0, 3));
  await store.finalize();
  let output = "";
  await store.writeLineageAirr(Int32Array.from([7, 0, 12]), async (part) => {
    output += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part);
  });
  const lines = output.trimEnd().split("\n");
  assert.equal(lines[0].split("\t").at(-1), "clone_id");
  assert.equal(lines[1].split("\t").at(-1), "swig_lineage_7");
  assert.equal(lines[2].split("\t").at(-1), "");
  assert.equal(lines[3].split("\t").at(-1), "swig_lineage_12");
  await store.clear();
});

test("masked AIRR scans yield only active ordinals and report filtered progress totals", async () => {
  const store = new AirrResultStore();
  await store.appendBatch(header, makeBody(0, 8));
  await store.finalize();
  const mask = Uint8Array.from([0, 1, 0, 0, 1, 0, 0, 1]);
  const observed: number[] = [];
  const progress: Array<[number, number]> = [];
  await store.scanAirrRows(["sequence_id"], (rows) => {
    observed.push(...rows.map((row) => row.ordinal));
  }, {
    includeMask: mask,
    batchSize: 100,
    onProgress: (processed, total) => progress.push([processed, total]),
  });
  assert.deepEqual(observed, [1, 4, 7]);
  assert.deepEqual(progress.at(-1), [3, 3]);
  await store.clear();
});

test("double-D evidence is indexed and exported separately without rewriting the AIRR table", async () => {
  const store = new AirrResultStore();
  const doubleDHeader = [
    "swig_batch_record_index", "sequence_id", "standard_d_call", "d_call", "d2_call",
    "d_identity", "d2_identity", "d_sequence_start", "d_sequence_end",
    "d2_sequence_start", "d2_sequence_end", "np2", "np3", "swig_double_d_mode",
  ].join("\t");
  const doubleDBody = [
    "1", "read_1", "IGHD2*01", "IGHDA*01", "IGHDB*01", "0.99", "0.98",
    "110", "122", "129", "141", "AACCGG", "TT", "all",
  ].join("\t") + "\n";
  await store.appendBatch(header, makeBody(0, 3), { header: doubleDHeader, body: doubleDBody });
  const secondDoubleDBody = [
    "0", "read_3", "IGHD4*01", "IGHDC*01", "IGHDD*01", "0.97", "0.96",
    "112", "124", "130", "143", "AACCC", "T", "long_span",
  ].join("\t") + "\n";
  await store.appendBatch(header, makeBody(3, 2), { header: doubleDHeader, body: secondDoubleDBody });
  await store.finalize();

  assert.equal(store.doubleDCount, 2);
  const doubleDPage = await store.page({ ...EMPTY_FILTERS, hasDoubleD: true }, 0, 10);
  assert.equal(doubleDPage.rows.length, 2);
  assert.equal(doubleDPage.rows[0].ordinal, 1);
  assert.equal(doubleDPage.rows[1].ordinal, 3, "batch-local evidence indexes were not mapped to global AIRR ordinals");
  assert.equal(doubleDPage.rows[0].dCall, "IGHDA*01");
  assert.equal(doubleDPage.rows[0].d2Call, "IGHDB*01");
  assert.equal(doubleDPage.rows[0].dIdentity, 0.99);
  const detail = await store.detail(doubleDPage.rows[0]);
  assert.equal(detail.standard_d_call, "IGHD2*01");
  assert.equal(detail.d_call, "IGHDA*01");
  assert.equal(detail.d2_call, "IGHDB*01");
  assert.equal(detail.swig_batch_record_index, undefined);

  const main = await store.airrBlob();
  const mainText = await main.text();
  assert.ok(!mainText.split("\n", 1)[0].includes("d2_call"));
  const mainHeader = mainText.split("\n", 1)[0].split("\t");
  const mainRows = mainText.trimEnd().split("\n").slice(1).map((line) => line.split("\t"));
  assert.equal(mainRows[1][mainHeader.indexOf("d_call")], "IGHD2*01");

  let evidence = "";
  await store.writeDoubleD(async (part) => {
    evidence += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part);
  });
  const evidenceLines = evidence.trimEnd().split("\n");
  assert.ok(evidenceLines[0].startsWith("swig_airr_ordinal\tsequence_id"));
  assert.ok(!evidenceLines[0].includes("swig_batch_record_index"));
  assert.equal(evidenceLines[1].split("\t", 1)[0], "2");
  assert.ok(evidenceLines[1].includes("IGHDB*01"));
  assert.equal(evidenceLines[2].split("\t", 1)[0], "4");
  await store.clear();
});

test("direct output keeps byte offsets usable for on-demand detail", async () => {
  const parts: BlobPart[] = [];
  let file = new File([], "direct.airr.tsv");
  const store = new AirrResultStore({
    handle: { getFile: async () => file },
    writable: {
      write: async (part) => { parts.push(part instanceof Uint8Array ? part.slice().buffer : part); },
      close: async () => { file = new File(parts, "direct.airr.tsv", { type: "text/tab-separated-values" }); },
      abort: async () => { parts.length = 0; },
    },
  });
  await store.appendBatch(header, makeBody(0, 2));
  await store.finalize();
  const page = await store.page({ ...EMPTY_FILTERS }, 0, 2);
  const detail = await store.detail(page.rows[1]);
  assert.equal(detail.sequence_id, "read_1");
  assert.equal(file.size, store.outputBytes);
  assert.equal((await store.airrBlob()).size, file.size);
  await store.clear();
});
