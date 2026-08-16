import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

import { logoGlyphRun } from "../src/logo-glyphs.ts";
import { AMINO_ACID_LOGO_SYMBOLS, NUCLEOTIDE_LOGO_SYMBOLS } from "../src/probability-logo.ts";

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

test("concatenated gzip uploads require an explicit merged-versus-separate sample choice", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const gzip = fs.readFileSync(new URL("../src/gzip-members.ts", import.meta.url), "utf8");

  assert.match(app, /CONCATENATED GZIP · SAMPLE STRUCTURE/);
  assert.match(app, /Import as \{pendingGzipImport\.members\.length\} separate samples/);
  assert.match(app, /Merge into one sample/);
  assert.match(app, /Independent sample, donor, cohort, timepoint, and tissue fields/);
  assert.match(app, /gzipMemberSource\(candidate\.file,members\)/);
  assert.match(gzip, /startsAtRecordBoundary/);
  assert.match(gzip, /multi-member candidates are each decompressed once/);
});

test("session restoration applies saved metadata during the first index pass", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const restore = app.slice(app.indexOf("async function restoreSavedSession"), app.indexOf("async function restoreLinkedAirr"));
  const configure = restore.indexOf("configureStudyMetadataForImport(restoredDatasets)");
  const append = restore.indexOf("store.appendBatch(header,body)");
  assert.ok(configure >= 0 && configure < append, "saved metadata must be installed before AIRR indexing starts");
  assert.doesNotMatch(restore, /updateStudyMetadata/);
  assert.doesNotMatch(restore, /Applying saved study metadata to local indexes/);
});

test("assignment is one progressive action page while independent result tools remain contextual", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const repertoire = fs.readFileSync(new URL("../src/repertoire-charts.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const css = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");

  assert.match(app, /analysis-layout single-action-layout/);
  assert.doesNotMatch(app, /AnalysisWorkspace|analysisWorkspace|Analysis setup sections/);
  assert.doesNotMatch(app, /settings-strip/);
  assert.match(app, /<b>Advanced options<\/b>/);
  assert.match(app, /Customize individual loci, V\/D\/J\/C sources, or allele inclusion/);
  assert.match(app, /aria-label="Sequence result panels"/);
  assert.match(app, /setSequenceWorkspace\("detail"\)/);
  assert.match(app, /Records \+ filters/);
  assert.match(app, /sequence-records-stack/);
  assert.doesNotMatch(app, /useState<"filters"\|"records"/);
  assert.match(app, /context-rail-secondary/);
  assert.match(app, /<span>Double-D explorer<small>/);
  assert.doesNotMatch(app, /results-tab-double-d|results-panel-double-d/);
  assert.match(app, /results-hero" hidden aria-hidden="true"/);
  assert.match(app, /className="results-rail-tools"/);
  assert.match(app, /sidebarTools=\{resultsSidebarTools\(\)\}/);
  assert.match(app, /overview=\{repertoireOverview\}/);
  assert.match(app, /className="results-tool-drawer"/);
  assert.match(app, /<strong>Study design<\/strong>/);
  assert.match(app, /<strong>Sample colors<\/strong>/);
  assert.match(app, /<b>01<\/b><span>Repertoire/);
  assert.match(app, /FacetPicker label=\{`\$\{segment\} gene or allele`\}/);
  assert.match(app, /Include multi-call assignments containing this/);
  assert.match(app, /addFieldHelp\(root\)/);
  assert.match(app, /<select value=\{species\?\.name\?\?""\}/);
  assert.doesNotMatch(app, /speciesDraft/);
  assert.match(app, /setPendingDatabaseId\(nextId\)/);
  assert.match(app, /aria-busy=\{busy\} value=\{value\} onChange=/);

  assert.match(repertoire, /aria-label="Repertoire panels"/);
  assert.match(repertoire, /panel==="usage"/);
  assert.doesNotMatch(repertoire, /panel==="controls"|setPanel\("controls"\)/);
  assert.match(repertoire, /repertoire-global-settings/);
  assert.match(repertoire, /\{sidebarTools\}/);
  assert.match(repertoire, /\{overview\}/);

  assert.match(post, /type PostWorkspaceId = "overview" \| PostModuleId/);
  assert.match(post, /aria-label="Post-analysis sections"/);
  assert.match(post, /new Set<PostModuleId>\(activeWorkspace === "overview" \? \[\] : \[activeWorkspace\]\)/);
  assert.match(post, /Advanced \{collapseMode === "fad"/);
  assert.match(post, /Advanced call matching and performance settings/);
  assert.match(post, /Advanced missing-allele evidence thresholds/);
  assert.match(post, /Study metadata filters/);
  assert.match(post, /Additional D, constant-gene, isotype, and rare-event filters/);
  assert.match(post, /FacetPicker label="Sample"/);
  assert.match(post, /FacetPicker label="V gene or allele"/);
  assert.match(post, /vCallIncludeAmbiguous/);
  assert.match(post, /<span>Allele pooling<small>/);
  assert.match(post, /Resolve ambiguous germline calls by pooling repertoire evidence/);
  assert.match(post, /<AlleleRefinementPanel references=\{references\} options=/);
  assert.match(post, /\{sidebarTools\}/);
  assert.match(fs.readFileSync(new URL("../src/allele-refinement/panel.tsx", import.meta.url), "utf8"), /Advanced evidence-kernel and model settings/);
  assert.match(css, /\.post-context-main > \.post-module\.is-collapsed \{ display: none !important; \}/);
  assert.match(css, /\.pipeline-stage-grid article:not\(\.enabled\) \.pipeline-fields/);
  assert.match(css, /\.post-module > label\.constant-collapse-policy/);
  assert.match(css, /\.results-table col\.column-cdr3/);
  assert.match(css, /\.sequence-context-main \.detail-shell[\s\S]*overflow-x: auto/);
  assert.match(css, /results are an application workspace, not a report masthead/);
  assert.match(css, /\.results-application-page \.results-view-tabs[\s\S]*background: var\(--ink\)/);
  assert.match(css, /\.post-analysis-heading \{ display: none !important; \}/);
});

test("allele pooling exposes parameter-responsive reference and hard assignment diagnostics", () => {
  const panel = fs.readFileSync(new URL("../src/allele-refinement/panel.tsx", import.meta.url), "utf8");
  const views = fs.readFileSync(new URL("../src/allele-refinement/diagnostic-views.tsx", import.meta.url), "utf8");
  const diagnostics = fs.readFileSync(new URL("../src/allele-refinement/diagnostics.ts", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  assert.match(panel, /ReferenceKernelInspector references=\{references\} options=\{options\}/);
  assert.match(panel, /AlleleAssignmentShiftChart results=\{segmentResults\}/);
  assert.match(panel, /Best posterior if confidence passes/);
  assert.match(panel, /Best posterior for every modeled record/);
  assert.match(panel, /Fast hurdle active set/);
  assert.match(panel, /Continuous Dirichlet mixture/);
  assert.match(panel, /Inclusion posterior threshold/);
  assert.doesNotMatch(panel, /allele-model-table/);
  assert.doesNotMatch(panel, /model\.alleles\.slice\(0,8\)/);
  assert.match(views, /label="Reference allele"/);
  assert.match(views, />Differences<\/button>/);
  assert.match(views, /row\.primary \? <small>bar omitted<\/small>/);
  assert.match(views, /Best-match allele counts before and after reassignment/);
  assert.match(views, /Vanished alleles only/);
  assert.match(views, /local_best_count/);
  assert.match(views, /SVG ↓/);
  assert.match(views, /Download surviving allele reference/);
  assert.match(views, /minimumReferenceReads/);
  assert.match(diagnostics, /unstripped\.some\(\(row\) => row\.sequence\[column\] !== "-"\)/);
  assert.match(diagnostics, /export function hardAssignmentShiftData/);
  assert.match(diagnostics, /export function survivingAlleleReference/);
  assert.match(diagnostics, /vanishes: beforeCount > 0 && afterCount === 0/);
  assert.match(panel, /Cross-donor override/);
  assert.match(panel, /Evidence never crosses participants/);
  assert.match(post, /alleleRuntime\.run\(store, references, alleleOptions, null,/);
  const automaticAlleles = post.indexOf("if (autoPipeline.alleleRefinement.enabled)");
  const automaticCollapse = post.indexOf("if (autoPipeline.collapse.enabled)");
  assert.ok(automaticAlleles >= 0 && automaticAlleles < automaticCollapse, "automatic pooling must run before collapse");
  const manualCollapse = post.slice(post.indexOf("async function runDedup()"), post.indexOf("async function applyDedupFilter()"));
  assert.doesNotMatch(manualCollapse, /setRepertoireCallOverrides\(null\)/);
  assert.match(app, /1 · Repertoire allele pooling/);
  assert.match(app, /Uses policy-selected V\/D\/J calls/);
});

test("phylogenetic UCA posterior uses contour-bounded embedded glyphs and aligned HMM tracks", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const component = fs.readFileSync(new URL("../src/probability-logo.tsx", import.meta.url), "utf8");
  const glyphs = fs.readFileSync(new URL("../src/logo-glyphs.ts", import.meta.url), "utf8");
  const panel = fs.readFileSync(new URL("../src/phylo-uca/panel.tsx", import.meta.url), "utf8");
  const annotation = fs.readFileSync(new URL("../src/phylo-uca/hmm-annotation.tsx", import.meta.url), "utf8");
  const annotationModel = fs.readFileSync(new URL("../src/phylo-uca/hmm-annotation-model.ts", import.meta.url), "utf8");
  const hmm = fs.readFileSync(new URL("../src/phylo-uca/hmm.ts", import.meta.url), "utf8");
  const styles = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");
  assert.match(component, /LOGO_MONOSPACE_FONT/);
  assert.match(component, /logoGlyphRun/);
  assert.match(component, /<path/);
  assert.match(component, /run\.yMax \* scaleY/);
  assert.match(component, /\$\{-scaleY\}/);
  assert.doesNotMatch(component, /getBBox|useLogoGlyphBounds|data-logo-measurement/);
  assert.match(glyphs, /Literal glyph contours from DejaVu Sans Mono Bold/);
  assert.match(glyphs, /"Q": \{[^\n]+yMin: -281/);
  assert.match(component, /width=\{columnWidth\}/);
  assert.doesNotMatch(component, /textLength=/);
  assert.match(component, /letter height is marginal frequency and is not scaled by entropy/i);
  assert.match(panel, /UCA posterior frequency logo/);
  assert.match(panel, /Logo SVG ↓/);
  assert.match(panel, /Full tracks \+ logo SVG ↓/);
  assert.match(panel, /Visible tracks \+ logo SVG ↓/);
  assert.match(panel, /serializePhyloUcaTrackLogoSvg/);
  assert.match(panel, /bottomAnnotations=\{logoBottomAnnotations\}/);
  assert.match(panel, /PhyloUcaPlacementMap/);
  assert.match(panel, /V\/J nucleotide mixture · default/);
  assert.match(panel, /"Full-HMM edges"[\s\S]{0,120}?\(0 = all\)/);
  assert.match(panel, /Conditional ML · fastest/);
  assert.match(panel, /Explicit grid marginalization/);
  assert.match(panel, /Gibbs\/MH · default · continuous placement/);
  assert.match(panel, /Exact pendant-length grid/);
  assert.match(panel, /does not use the grid settings/);
  assert.match(panel, /Gibbs\/MH mixing/);
  assert.match(panel, />Best path</);
  assert.match(panel, />Marginalized</);
  assert.match(panel, /PhyloUcaHmmAnnotationTracks/);
  assert.match(panel, /collapseAndOrderHmmAnnotationTracks/);
  assert.match(panel, /labelOffset=\{annotationScrollLeft\}/);
  assert.match(panel, /Reset all UCA settings to defaults/);
  assert.match(panel, /Additional-D probability[\s\S]{0,500}?max="1"/);
  assert.doesNotMatch(panel, /Additional-D probability[\s\S]{0,500}?max="0\.5"/);
  assert.match(panel, /adoptedResultSnapshotRef\.current === incomingGeneratedAt/);
  assert.match(panel, /adoptedResultSnapshotRef\.current = inference\.generatedAt/);
  assert.doesNotMatch(panel, /\[matchingInitial, result\?\.generatedAt\]/);
  assert.match(app, /addEventListener\("beforeunload", warnBeforeLeaving\)/);
  assert.match(app, /addEventListener\("popstate", interceptHistoryDeparture\)/);
  assert.match(app, /setLeavePrompt\(true\)/);
  assert.match(app, /Stay on this page/);
  assert.match(app, /Leave anyway/);
  assert.match(styles, /overscroll-behavior-x:\s*none/);
  assert.match(panel, /position \* logoColumnWidth \/ 3/);
  assert.match(panel, />Codon</);
  assert.match(panel, /Amino acid/);
  assert.match(panel, /does not multiply nucleotide marginals/);
  assert.match(annotation, /total glyph height is the posterior track occupancy/i);
  assert.match(annotation, /FittedLogoGlyph/);
  assert.match(annotation, /nucleotide mass/);
  assert.match(annotation, /phylo-uca-track-label-panel/);
  assert.match(annotation, /exact scrolled crop/i);
  assert.match(annotation, /viewBox.*\$\{x\}.*\$\{y\}.*\$\{width\}.*\$\{height\}/);
  assert.match(annotation, /root\.append\(logo\)/);
  assert.match(annotation, /trackHeight \+ logoSize\.height/);
  assert.match(annotationModel, /combinedTrack\(`display\|\$\{key\}`/);
  assert.match(annotationModel, /weightedCenter/);
  assert.match(annotationModel, /mode === "marginalized"/);
  assert.match(annotationModel, /NT · all non-template mass/);
  assert.match(styles, /\.phylo-uca-track-legend \.v/);
  assert.match(styles, /\.phylo-uca-track-label-panel/);
  assert.match(hmm, /alignment site minus D-reference position/);
  assert.match(hmm, /Exact P\(x_i,x_\{i\+1\},x_\{i\+2\}/);

  for (const symbol of new Set([...NUCLEOTIDE_LOGO_SYMBOLS, ...AMINO_ACID_LOGO_SYMBOLS])) {
    const run = logoGlyphRun(symbol);
    assert.ok(run, `missing embedded outline for ${symbol}`);
    assert.equal(run.width, run.xMax - run.xMin);
    assert.equal(run.height, run.yMax - run.yMin);
    assert.ok(run.width > 0 && run.height > 0, `degenerate embedded outline for ${symbol}`);
  }
  const codon = logoGlyphRun("A-Q");
  assert.ok(codon && codon.paths.length === 3, "multi-character codon glyph runs must retain all contours");
});
