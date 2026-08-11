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

async function sampledIds(source: string, size: number, seed: number) {
  const ids: string[] = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format: 1,
    batchSize: 3,
    subsample: { size, seed },
    onProgress: (next) => { progress = next; },
  })) {
    ids.push(...[...batch.text.matchAll(/^>([^\n]+)/gm)].map((match) => match[1]));
  }
  return { ids, progress: progress as SequenceStreamProgress | null };
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

test("seeded reservoir subsampling is exact, reproducible, and input-order stable", async () => {
  const source = Array.from({ length: 100 }, (_, index) => `>read_${index}\nACGT${index % 10}\n`).join("");
  const first = await sampledIds(source, 10, 73);
  const again = await sampledIds(source, 10, 73);
  const other = await sampledIds(source, 10, 74);
  assert.deepEqual(first.ids, again.ids);
  assert.notDeepEqual(first.ids, other.ids);
  assert.equal(first.ids.length, 10);
  assert.deepEqual(first.ids, [...first.ids].sort((a, b) => Number(a.slice(5)) - Number(b.slice(5))));
  assert.equal(first.progress?.recordsRead, 100);
  assert.equal(first.progress?.recordsSelected, 10);
});

test("AIRR reservoir sampling keeps the header and scans all source rows", async () => {
  const source = `sequence_id\tsequence\n${Array.from({ length: 25 }, (_, index) => `r${index}\tACGT`).join("\n")}\n`;
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format: 3,
    batchSize: 4,
    subsample: { size: 7, seed: 9 },
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  assert.deepEqual(batches.map((batch) => batch.count), [4, 3]);
  assert.ok(batches.every((batch) => batch.text.startsWith("sequence_id\tsequence\n")));
  assert.equal(progress?.recordsRead, 25);
  assert.equal(progress?.recordsSelected, 7);
});

test("50k gzip FASTA stress path and 10k reservoir sample remain batch-bounded", async () => {
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

  let sampledCount = 0;
  let sampledBatches = 0;
  let sampledProgress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source: file,
    format: 1,
    batchSize: 1000,
    subsample: { size: 10_000, seed: 2024 },
    onProgress: (next) => { sampledProgress = next; },
  })) {
    sampledCount += batch.count;
    sampledBatches += 1;
    assert.ok(batch.count <= 1000);
  }
  assert.equal(sampledCount, 10_000);
  assert.equal(sampledBatches, 10);
  assert.equal(sampledProgress?.recordsRead, 50_000);
  assert.equal(sampledProgress?.recordsSelected, 10_000);
  assert.ok((sampledProgress?.maxBatchCharacters ?? Infinity) < 400_000);
});
