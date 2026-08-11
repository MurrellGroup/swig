import assert from "node:assert/strict";
import test from "node:test";
import { gzipSync } from "node:zlib";

import {
  streamSequenceBatches,
  type SequenceStreamProgress,
} from "../src/sequence-stream.ts";

async function collect(source: string | File, format: 1 | 2 | 3, batchSize = 2) {
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format,
    batchSize,
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  return { batches, progress: progress as SequenceStreamProgress | null };
}

test("incremental FASTA parser normalizes wrapped records into bounded batches", async () => {
  const result = await collect(">one notes\nAC GT\nTG\n>two\nAAAA\n>three\nCCCC\n", 1);
  assert.deepEqual(result.batches.map((batch) => batch.count), [2, 1]);
  assert.equal(result.batches[0].text, ">one notes\nACGTTG\n>two\nAAAA\n");
  assert.equal(result.progress?.recordsRead, 3);
});

test("incremental FASTQ parser supports wrapped sequence and quality lines", async () => {
  const result = await collect(
    "@one notes\nACGT\nTG\n+one\nIIII\nII\n@two\nAAAA\n+\n@@@@\n",
    2,
    10,
  );
  assert.equal(result.batches[0].count, 2);
  assert.equal(result.batches[0].text, "@one notes\nACGTTG\n+one\nIIIIII\n@two\nAAAA\n+\n@@@@\n");
});

test("gzip FASTQ is decompressed and parsed without whole-file text conversion", async () => {
  const source = "@gzip_one\nACGTACGT\n+\nIIIIIIII\n@gzip_two\nTGCATGCA\n+\nHHHHHHHH\n";
  const file = new File([gzipSync(source)], "reads.fastq.gz", { type: "application/gzip" });
  const result = await collect(file, 2, 1);
  assert.deepEqual(result.batches.map((batch) => batch.count), [1, 1]);
  assert.equal(result.progress?.bytesRead, file.size);
  assert.ok(result.batches[1].text.startsWith("@gzip_two\n"));
});

test("AIRR input retains its header in every compute batch", async () => {
  const source = "sequence_id\tsequence\tlocus\na\tACGT\tIGH\nb\tTGCA\tTRB\nc\t\tIGH\nd\tAAAA\tIGK\n";
  const result = await collect(source, 3, 2);
  assert.deepEqual(result.batches.map((batch) => batch.count), [2, 1]);
  assert.ok(result.batches.every((batch) => batch.text.startsWith("sequence_id\tsequence\tlocus\n")));
  assert.ok(!result.batches.some((batch) => batch.text.includes("\nc\t")));
});

test("50k gzip FASTA stress path stays within one 1000-record batch", async () => {
  const sequence = "ACGT".repeat(75);
  const records = Array.from({ length: 50_000 }, (_, index) => `>read_${index}\n${sequence}\n`).join("");
  const compressed = gzipSync(records, { level: 1 });
  const file = new File([compressed], "stress.fasta.gz", { type: "application/gzip" });
  let count = 0;
  let batches = 0;
  let finalProgress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source: file,
    format: 1,
    batchSize: 1000,
    onProgress: (next) => { finalProgress = next; },
  })) {
    count += batch.count;
    batches += 1;
    assert.ok(batch.count <= 1000);
    assert.ok(batch.text.length < 400_000);
  }
  assert.equal(count, 50_000);
  assert.equal(batches, 50);
  assert.equal(finalProgress?.bytesRead, file.size);
  assert.equal(finalProgress?.recordsRead, 50_000);
  assert.ok((finalProgress?.maxBatchCharacters ?? Infinity) < 400_000);
  assert.ok((finalProgress?.maxCarryCharacters ?? Infinity) < 10_000);
});
