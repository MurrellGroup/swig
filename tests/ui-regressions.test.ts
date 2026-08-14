import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

test("clickable individual-sequence views retain CDR3 context", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");

  // Paged AIRR records show AA with nucleotide fallback and retain nucleotide
  // context when both forms are present.
  assert.match(app, /row\.cdr3Aa \|\| row\.cdr3 \|\| "—"/);
  assert.match(app, /row\.cdr3Aa && row\.cdr3/);

  // The sparse VDDJ explorer hydrates CDR3 from the main index rather than
  // expecting it to be duplicated in the double-D sidecar.
  assert.match(app, /cdr3:record\.values\.cdr3\|\|index\?\.cdr3/);
  assert.match(app, /<th>CDR3<\/th>/);

  // Chimera candidates and lineage-member click targets also expose CDR3.
  assert.match(post, /indexed\?\.cdr3Aa\|\|indexed\?\.cdr3\|\|"CDR3 —"/);
  assert.match(post, /row\.values\.cdr3_aa\|\|row\.values\.cdr3\|\|"CDR3 —"/);
});

test("FASTQ quality control is an initially collapsed pre-assignment step", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");

  // No `open` property means the native details box is initially collapsed.
  assert.match(app, /<details className=\{`fastq-quality-control/);
  assert.doesNotMatch(app, /<details[^>]*fastq-quality-control[^>]*\sopen(?:=|\s|>)/);
  assert.match(app, /3′ trim → minimum retained length → expected-error filter → optional random subsample → V\(D\)J assignment/);
  assert.match(app, /FASTA and AIRR records will pass through this step unchanged/);
});

test("analysis and result tools use contextual single-panel workspaces", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const repertoire = fs.readFileSync(new URL("../src/repertoire-charts.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const css = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");

  assert.match(app, /aria-label="Analysis setup sections"/);
  assert.match(app, /analysisWorkspace!=="data"/);
  assert.match(app, /aria-label="Sequence result panels"/);
  assert.match(app, /setSequenceWorkspace\("detail"\)/);
  assert.match(app, /aria-label="Double-D panels"/);

  assert.match(repertoire, /aria-label="Repertoire panels"/);
  assert.match(repertoire, /panel==="usage"/);

  assert.match(post, /type PostWorkspaceId = "overview" \| PostModuleId/);
  assert.match(post, /aria-label="Post-analysis sections"/);
  assert.match(post, /new Set<PostModuleId>\(activeWorkspace === "overview" \? \[\] : \[activeWorkspace\]\)/);
  assert.match(css, /\.post-context-main > \.post-module\.is-collapsed \{ display: none !important; \}/);
});
