import assert from "node:assert/strict";
import test from "node:test";

import { donorForFlatRoot, inferDirectoryDonors } from "../src/directory-input.ts";
import { activeProjectRun, attachProjectDirectory, loadActiveProjectFiles, prepareProjectRun, saveProjectCheckpoint, writeProjectDatasetManifest } from "../src/project-directory.ts";
import type { SwigSession } from "../src/session-state.ts";

async function asBytes(value: unknown): Promise<Uint8Array> {
  if (typeof value === "string") return new TextEncoder().encode(value);
  if (value instanceof Blob) return new Uint8Array(await value.arrayBuffer());
  if (value instanceof ArrayBuffer) return new Uint8Array(value);
  if (ArrayBuffer.isView(value)) return new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
  throw new Error("Unsupported memory write");
}

class MemoryFileHandle {
  readonly kind = "file";
  bytes = new Uint8Array();
  readonly name: string;
  constructor(name: string) { this.name = name; }
  async getFile() { return new File([this.bytes], this.name); }
  async createWritable(options?: { keepExistingData?: boolean }) {
    let draft = options?.keepExistingData ? this.bytes.slice() : new Uint8Array();
    let position = 0;
    return {
      write: async (value: unknown) => {
        const bytes = await asBytes(value); const end = position + bytes.length;
        if (end > draft.length) { const expanded = new Uint8Array(end); expanded.set(draft); draft = expanded; }
        draft.set(bytes, position); position = end;
      },
      seek: async (next: number) => { position = next; },
      close: async () => { this.bytes = draft; },
      abort: async () => undefined,
    };
  }
}

class MemoryDirectoryHandle {
  readonly kind = "directory";
  private directories = new Map<string, MemoryDirectoryHandle>();
  private files = new Map<string, MemoryFileHandle>();
  readonly name: string;
  constructor(name: string) { this.name = name; }
  async getDirectoryHandle(name: string, options?: { create?: boolean }) {
    let result = this.directories.get(name);
    if (!result && options?.create) { result = new MemoryDirectoryHandle(name); this.directories.set(name, result); }
    if (!result) throw new DOMException("missing", "NotFoundError");
    return result;
  }
  async getFileHandle(name: string, options?: { create?: boolean }) {
    let result = this.files.get(name);
    if (!result && options?.create) { result = new MemoryFileHandle(name); this.files.set(name, result); }
    if (!result) throw new DOMException("missing", "NotFoundError");
    return result;
  }
}

test("nested directory names supply donor suggestions", () => {
  const plan = inferDirectoryDonors([
    { relativePath: "study/donor_A/day0.fasta", rootName: "study", fromDirectory: true },
    { relativePath: "study/donor_A/day14.fastq.gz", rootName: "study", fromDirectory: true },
    { relativePath: "study/donor_B/blood/sample.fasta", rootName: "study", fromDirectory: true },
  ]);
  assert.deepEqual(plan.flatRoots, []);
  assert.deepEqual(plan.suggestions.map((entry) => entry.donor), ["donor_A", "donor_A", "donor_B"]);
});

test("flat selected directories require a donor decision", () => {
  const plan = inferDirectoryDonors([
    { relativePath: "donor-run/a.fasta", rootName: "donor-run", fromDirectory: true },
    { relativePath: "donor-run/b.fastq", rootName: "donor-run", fromDirectory: true },
  ]);
  assert.deepEqual(plan.flatRoots, ["donor-run"]);
  assert.deepEqual(plan.suggestions.map((entry) => entry.donor), [null, null]);
  assert.equal(donorForFlatRoot("Donor 17 files"), "Donor_17_files");
});

test("ordinary independently selected files do not trigger directory grouping", () => {
  const plan = inferDirectoryDonors([{ relativePath: "a.fasta", rootName: "a.fasta", fromDirectory: false }]);
  assert.deepEqual(plan.flatRoots, []);
  assert.equal(plan.suggestions[0].donor, null);
});

test("project directories write a numbered run, latest state, checkpoint, manifest, and log", async () => {
  const root = new MemoryDirectoryHandle("study-project") as unknown as FileSystemDirectoryHandle;
  const attached = await attachProjectDirectory(root, "0.18.0");
  assert.equal(attached.existing, false);
  const prepared = await prepareProjectRun(attached.workspace, "longitudinal study", "0.18.0");
  await writeProjectDatasetManifest(attached.workspace,prepared.run,[{datasetId:"dataset_1",inputName:"donor_A/day0.fasta",sampleId:"day0",subjectId:"donor_A",cohort:"study",timepoint:"0",compartment:"blood",records:1}]);
  await prepared.writable.write("sequence_id\tsequence\nread_1\tACGT\n");
  await prepared.writable.close();
  const session = {
    schema: 1, application: "Swig", applicationVersion: "0.18.0", savedAt: new Date(0).toISOString(),
    linkedAirr: { name: "longitudinal-study.airr.tsv", size: 34, lastModified: 0, records: 1, headers: ["sequence_id", "sequence"], fingerprint: "abcd" },
    analysis: { inputName: "longitudinal study", species: "Homo sapiens", scope: "IGH", workers: 1, minimumIdentity: .6, strand: 0, references: { V: "", D: "", J: "", C: "", counts: { V: 0, D: 0, J: 0, C: 0 }, annotation: { V: { annotated: 0, total: 0 }, J: { annotated: 0, total: 0 } }, loci: ["IGH"] } },
    doubleD: [], postAnalysis: { workingStages: [] },
  } satisfies SwigSession;
  await saveProjectCheckpoint(attached.workspace, session, new Blob([JSON.stringify(session)]), "test");
  assert.equal(activeProjectRun(attached.workspace)?.id, "001-longitudinal-study");
  assert.equal(activeProjectRun(attached.workspace)?.checkpointCount, 1);
  const reopened=await attachProjectDirectory(root,"0.18.0");
  assert.equal(reopened.existing,true);
  const restored = await loadActiveProjectFiles(reopened.workspace);
  assert.match(await restored.sessionFile.text(), /longitudinal study/);
  assert.match(await restored.airrFile.text(), /read_1/);
  const runDirectory=await root.getDirectoryHandle("runs");
  const activeDirectory=await runDirectory.getDirectoryHandle("001-longitudinal-study");
  const inputDirectory=await activeDirectory.getDirectoryHandle("inputs");
  assert.match(await (await inputDirectory.getFileHandle("datasets.tsv")).getFile().then((file)=>file.text()),/donor_A\/day0\.fasta/);
  const manifest = await (await root.getFileHandle("swig-project.json")).getFile();
  assert.match(await manifest.text(), /"checkpointCount": 1/);
});
