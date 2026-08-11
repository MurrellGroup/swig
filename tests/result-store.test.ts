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
