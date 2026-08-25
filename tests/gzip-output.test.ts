import assert from "node:assert/strict";
import test from "node:test";
import { gunzipSync } from "node:zlib";

import { compressedExtension, compressedMime, compressedName, outputWriter } from "../src/gzip-output.ts";

function collectingSink(parts: Uint8Array[]) {
  return {
    write: async (value: string | Blob | Uint8Array) => {
      if (typeof value === "string") parts.push(new TextEncoder().encode(value));
      else if (value instanceof Blob) parts.push(new Uint8Array(await value.arrayBuffer()));
      else parts.push(value.slice());
    },
  };
}

test("bounded gzip writer round-trips strings, bytes, and streamed Blobs", async () => {
  const parts: Uint8Array[] = [];
  const writer = outputWriter(collectingSink(parts), "gzip");
  await writer.write("sequence_id\tsequence\n");
  await writer.write(new Uint8Array([114, 49, 9]));
  await writer.write(new Blob(["ACGT\n"]));
  await writer.close();
  const compressed = Buffer.concat(parts.map((part) => Buffer.from(part)));
  assert.equal(compressed[0], 0x1f);
  assert.equal(compressed[1], 0x8b);
  assert.equal(gunzipSync(compressed).toString("utf8"), "sequence_id\tsequence\nr1\tACGT\n");
});

test("uncompressed writer preserves bytes and output naming is explicit", async () => {
  const parts: Uint8Array[] = [];
  const writer = outputWriter(collectingSink(parts), "none");
  await writer.write("x\n");
  await writer.close();
  assert.equal(Buffer.concat(parts.map((part) => Buffer.from(part))).toString(), "x\n");
  assert.equal(compressedName("calls.airr.tsv", "gzip"), "calls.airr.tsv.gz");
  assert.equal(compressedName("calls.airr.tsv.gz", "gzip"), "calls.airr.tsv.gz");
  assert.equal(compressedExtension(".tsv", "gzip"), ".gz");
  assert.equal(compressedMime("text/plain", "gzip"), "application/gzip");
});
