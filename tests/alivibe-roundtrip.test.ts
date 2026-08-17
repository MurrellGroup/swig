import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

import {
  ALIVIBE_BRIDGE_VERSION,
  ALIVIBE_SOURCE_REVISION,
  assertAlivibeInitialLoad,
  assertAlivibeRoundTripTarget,
  getAlivibeBridge,
  loadAlivibeNucleotideFasta,
  readAlivibeNucleotideFasta,
  type AlivibeEditorWindow,
} from "../src/alivibe-roundtrip.ts";

function editorWithSnapshot(snapshot: Record<string, unknown>): AlivibeEditorWindow {
  return {
    swigAlivibeBridge: {
      version: ALIVIBE_BRIDGE_VERSION,
      sourceRevision: ALIVIBE_SOURCE_REVISION,
      loadNucleotideFasta: () => snapshot,
      snapshotNucleotide: () => snapshot,
      installMsaRunner: () => undefined,
      createMsaJob: () => null,
    },
  } as unknown as AlivibeEditorWindow;
}

test("Alivibe returns the exact complete ordered NT rows, including manually placed gaps", () => {
  const fasta = ">germline_N_masked\nACGT---NN\n>member_1\nA-GTT--NN\n";
  const snapshot = {
    version: ALIVIBE_BRIDGE_VERSION,
    sourceRevision: ALIVIBE_SOURCE_REVISION,
    alphabet: "NT",
    mode: "NT",
    frameOffset: 2,
    records: [
      { name: "germline_N_masked", sequence: "ACGT---NN" },
      { name: "member_1", sequence: "A-GTT--NN" },
    ],
    fasta,
  };
  const editor = editorWithSnapshot(snapshot);
  assert.equal(getAlivibeBridge(editor)?.sourceRevision, ALIVIBE_SOURCE_REVISION);
  assert.deepEqual(loadAlivibeNucleotideFasta(editor, fasta, 2), {
    fasta,
    frameOffset: 2,
    records: snapshot.records,
    sourceRevision: ALIVIBE_SOURCE_REVISION,
  });
  assert.deepEqual(readAlivibeNucleotideFasta(editor).records, snapshot.records);
});

test("Alivibe rejects stale, AA, empty, or internally inconsistent returns", () => {
  const base = {
    version: ALIVIBE_BRIDGE_VERSION,
    sourceRevision: ALIVIBE_SOURCE_REVISION,
    alphabet: "NT",
    mode: "NT",
    frameOffset: 0,
    records: [{ name: "a", sequence: "AC-G" }],
    fasta: ">a\nAC-G\n",
  };
  assert.throws(() => readAlivibeNucleotideFasta(editorWithSnapshot({ ...base, mode: "AA" })), /nucleotide view/i);
  assert.throws(
    () => readAlivibeNucleotideFasta(editorWithSnapshot({ ...base, sourceRevision: "stale" })),
    /does not match this Swig release/i,
  );
  assert.throws(
    () => readAlivibeNucleotideFasta(editorWithSnapshot({ ...base, fasta: ">a\nACGT\n" })),
    /records and nucleotide FASTA disagree/i,
  );
  assert.throws(
    () => readAlivibeNucleotideFasta(editorWithSnapshot({ ...base, records: [] })),
    /empty nucleotide alignment/i,
  );
});

test("a nucleotide return cannot cross lineage or alignment boundaries", () => {
  const transfer = readAlivibeNucleotideFasta(editorWithSnapshot({
    version: ALIVIBE_BRIDGE_VERSION,
    sourceRevision: ALIVIBE_SOURCE_REVISION,
    alphabet: "NT",
    mode: "NT",
    frameOffset: 0,
    records: [{ name: "lineage_7", sequence: "AC-G" }],
    fasta: ">lineage_7\nAC-G\n",
  }));
  assert.doesNotThrow(() => assertAlivibeInitialLoad(">lineage_7\nAC-G\n", transfer));
  assert.throws(() => assertAlivibeInitialLoad(">lineage_9\nTT-G\n", transfer), /did not load the exact/i);

  const origin = { groupKey: "7", alignmentFingerprint: "fasta-7" };
  assert.doesNotThrow(() => assertAlivibeRoundTripTarget(origin, { ...origin }));
  assert.throws(
    () => assertAlivibeRoundTripTarget(origin, { groupKey: "9", alignmentFingerprint: "fasta-9" }),
    /selected lineage changed/i,
  );
  assert.throws(
    () => assertAlivibeRoundTripTarget(origin, { groupKey: "7", alignmentFingerprint: "new-fasta-7" }),
    /alignment changed/i,
  );
});

test("the bundled Alivibe bridge snapshots the same NT state used by its viewer/export and initializes controls first", () => {
  const source = fs.readFileSync(new URL("../public/tools/alivibe.html", import.meta.url), "utf8");
  assert.match(source, new RegExp(ALIVIBE_SOURCE_REVISION));
  assert.match(source, new RegExp(`version:\\s*${ALIVIBE_BRIDGE_VERSION}`));
  assert.match(source, /window\.swigAlivibeBridge\s*=\s*Object\.freeze/);
  assert.match(source, /installMsaRunner\(runner\)/);
  assert.match(source, /createMsaJob\(sequences\)/);
  assert.match(source, /swigMsaRunner\(sequences\.map\(sequence => String\(sequence\)\), state\.mode === 'NT' \? 'nucleotide' : 'amino-acid'\)/);
  assert.match(source, /Alivibe-compatible WebAssembly worker/);
  assert.match(source, /activeMsaJob\.cancel\(\)/);
  assert.match(source, /setMode\('NT'\);\s*const records = state\.viewSequences\.map/);
  assert.match(source, /while\(entry\.seq\.length < width\)[\s\S]*entry\.seq\.push\('-'\)/);
  assert.match(source, /buildFastaContent\(rows, 0, null\)/);
  assert.ok(source.indexOf("dom.file.addEventListener") < source.indexOf("new Aioli"));
  const modeFunction = source.match(/function setMode\(m\) \{[\s\S]*?\n\}\n\nfunction recalc/)?.[0] ?? "";
  assert.match(modeFunction, /grpFrame\.style\.pointerEvents = 'auto'/);
  assert.doesNotMatch(modeFunction, /grpFrame\.style\.pointerEvents = 'none'/);
});

test("MurrellGroup components do not carry duplicate third-party license notices", () => {
  assert.equal(fs.existsSync(new URL("../public/tools/WebWidgets-LICENSE.txt", import.meta.url)), false);
  assert.equal(fs.existsSync(new URL("../wasm/LICENSE.swiftig", import.meta.url)), false);
  const toolsReadme = fs.readFileSync(new URL("../public/tools/README.md", import.meta.url), "utf8");
  const frameClean = fs.readFileSync(new URL("../public/tools/frameclean.js", import.meta.url), "utf8");
  assert.doesNotMatch(toolsReadme, /WebWidgets-LICENSE|license is retained/i);
  assert.doesNotMatch(frameClean, /License:\s*MIT/i);
});
