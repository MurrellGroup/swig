import assert from "node:assert/strict";
import test from "node:test";
import { inferLineageGermline, lineageInputFasta, lineageTemplateChoices, projectLineageTemplate, restoreLineageOuterFlanks, retainLineageTemplate, GERMLINE_OUTGROUP } from "../src/lineage-alignment.ts";
import { parseFasta } from "../src/post-analysis-core.ts";
import { defaultPhyloUcaOptions } from "../src/phylo-uca/defaults.ts";
import { inferPhyloUca } from "../src/phylo-uca/inference.ts";
import { mappedTemplateBoundaries, prepareObservedOnlyAlignment, preparePhyloUcaReferences } from "../src/phylo-uca/references.ts";
import type { AirrDetailRow } from "../src/result-store.ts";

const references = { V: ">IGHV1*01 metadata=preserved\nACGTCCTA\n", D: "", J: ">IGHJ1*01\nGGTTAC\n" };

function member(ordinal: number, partial = false, junction = "AAA"): AirrDetailRow {
  const vQuery = partial ? "GTCC" : "TCGTCC"; // All complete tips have an SHM at the V 5′ end.
  const vGermline = partial ? "GTCC" : "ACGTCC";
  const jQuery = partial ? "TT" : "TTAC";
  return {
    record: { ordinal, sequenceId: `read_${ordinal}`, locus: "IGH" },
    values: {
      sequence_id: `read_${ordinal}`, locus: "IGH", productive: "T",
      sequence_alignment: vQuery + junction + jQuery,
      germline_alignment: vGermline + "N".repeat(junction.length) + jQuery,
      v_call: "IGHV1*01", j_call: "IGHJ1*01",
      v_germline_start: partial ? "3" : "1", v_germline_end: "6",
      j_germline_start: "3", j_germline_end: partial ? "4" : "6",
      v_sequence_start: "1", v_sequence_end: String(vQuery.length),
      j_sequence_start: String(vQuery.length + junction.length + 1), j_sequence_end: String(vQuery.length + junction.length + jQuery.length),
      v_sequence_alignment: vQuery, v_germline_alignment: vGermline,
      j_sequence_alignment: jQuery, j_germline_alignment: jQuery,
      sequence_frame: partial ? "2" : "1",
    },
  } as AirrDetailRow;
}

test("complete outer V/J alignments outrank a cleaner truncated member; explicit selection overrides", () => {
  const rows = [member(0, true), member(1)];
  const choices = lineageTemplateChoices(rows, { references });
  assert.deepEqual(choices.map((choice) => [choice.ordinal, choice.complete]), [[1, true], [0, false]]);
  assert.equal(inferLineageGermline(rows, "closest", { references }).selectedOrdinal, 1);
  const manual = inferLineageGermline(rows, "closest", { references, templateOrdinal: 0 });
  assert.equal(manual.selectedOrdinal, 0);
  assert.equal(manual.template, "ACGTCCNNNTTAC");
  assert.equal(manual.restoredVPrefix, 2);
  assert.equal(manual.restoredJSuffix, 2);
  // V's junction-facing TA and J's junction-facing GG are NOT restored.
  assert.equal(manual.uca, "ACGTCCAAATTAC");
  assert.throws(() => inferLineageGermline(rows, "closest", { references, templateOrdinal: 99 }), /not in the current lineage/);
});

test("missing outer read bases stay gaps while the complete guide reaches both reference ends", () => {
  const rows = [member(0, true), member(1, true, "ATA"), member(2, true, "ACA")];
  const input = lineageInputFasta(rows, "closest", { references });
  assert.match(input.fasta, />read_0__1\n--GTCCAAATT--/);
  assert.match(input.fasta, />__germline_N_masked__\nACGTCCNNNTTAC/);
  assert.equal(input.frames.at(-1), 0); // Restored prefix adjusts the guide's codon phase.
  const observed = prepareObservedOnlyAlignment(input.fasta, GERMLINE_OUTGROUP);
  const metadata = rows.map((row) => ({ ordinal: row.record.ordinal, sequenceId: row.record.sequenceId, locus: "IGH", values: row.values }));
  assert.deepEqual(mappedTemplateBoundaries(observed.posteriorFasta, metadata, 0), { vEnd: 5, jStart: 9 });
  assert.equal(observed.columns, 9); // Tree-only columns with no tip evidence are removed.
  assert.equal(parseFasta(observed.posteriorFasta, true)[0].sequence.length, 13);
  assert.equal(observed.posteriorColumns.length, 13); // HMM still integrates all reference columns.
  assert.equal(inferLineageGermline(rows, "consensus", { references }).template, input.germline);
});

test("flank completion validates exact calls and coordinates, and retains uncertainty between tied alleles", () => {
  const row = member(0, true);
  row.values.v_call = "IGHV1*01,IGHV1*02";
  const ambiguous = { ...references, V: references.V + ">IGHV1*02\nTCGTCCTA\n" };
  const inferred = inferLineageGermline([row], "closest", { references: ambiguous });
  assert.equal(inferred.template, "NCGTCCNNNTTAC");
  const prepared = preparePhyloUcaReferences(inferred.template, [row], ambiguous, "IGH", defaultPhyloUcaOptions().candidates);
  assert.equal(prepared.v.find((candidate) => candidate.name === "IGHV1*01")?.projection.slice(0, 6), "ACGTCC");
  assert.equal(prepared.v.find((candidate) => candidate.name === "IGHV1*02")?.projection.slice(0, 6), "TCGTCC");
  row.values.v_call = "missing_allele";
  assert.equal(inferLineageGermline([row], "closest", { references }).template.slice(0, 2), "--");
  row.values.v_call = "IGHV1*01";
  row.values.v_germline_alignment = "AAAA";
  assert.equal(inferLineageGermline([row], "closest", { references }).template.slice(0, 2), "--");
});

test("manual template controls junction trimming and survives an MSA display limit", () => {
  const a = member(0), b = member(1);
  b.values.v_germline_end = "5";
  b.values.v_sequence_alignment = "TCGTC";
  b.values.v_germline_alignment = "ACGTC";
  b.values.sequence_alignment = "TCGTCAAATTAC";
  b.values.germline_alignment = "ACGTCNNNTTAC";
  const sample = retainLineageTemplate([a, b], [a], 1);
  assert.deepEqual(sample.map((row) => row.record.ordinal), [1]);
  assert.equal(lineageInputFasta(sample, "closest", { references, templateOrdinal: 1 }).germline, "ACGTCNNNTTAC");
  const guide = "ACGTCCNNNTTAC";
  // The guide row's selected member supplies the anchors, rather than the
  // median from the longer V alignments in other members.
  const opts = defaultPhyloUcaOptions().candidates;
  const shortBoundary = preparePhyloUcaReferences(guide, [a, { ...a, record: { ...a.record, ordinal: 2 } }, b].map((row) => ({ ordinal: row.record.ordinal, sequenceId: row.record.sequenceId, locus: "IGH", values: row.values })), references, "IGH", opts, 1);
  assert.equal(shortBoundary.vEndColumn, 4);
});

test("changing the template reprojects the guide through a truncated query and preserves curated tip gaps", () => {
  const rows = [member(0, true), member(1), member(2)];
  const old = ">read_0__1\n--GTCCAAATT--\n>read_1__2\nTCGTCCAAATTAC\n>read_2__3\nTCGTCCAAATTAC\n>__germline_N_masked__\n--GTCCNNNTT--\n";
  const projected = projectLineageTemplate(old, rows, { references, templateOrdinal: 0 });
  assert.equal(projected.addedLeftColumns, 0);
  const records = parseFasta(projected.fasta, true);
  assert.equal(records.find((record) => record.name === GERMLINE_OUTGROUP)?.sequence, "ACGTCCNNNTTAC");
  assert.deepEqual(records.filter((record) => record.name !== GERMLINE_OUTGROUP).map((record) => record.sequence), parseFasta(old, true).slice(0, 3).map((record) => record.sequence));
  const allPartial = ">read_0__1\nGTCCAAATT\n>read_3__4\nGTCCATATT\n>__germline_N_masked__\nGTCCNNNTT\n";
  const padded = projectLineageTemplate(allPartial, [member(0, true), member(3, true, "ATA")], { references, templateOrdinal: 0 });
  assert.equal(padded.addedLeftColumns, 2);
  assert.match(padded.fasta, />read_0__1\n--GTCCAAATT--/);
  assert.throws(() => projectLineageTemplate(old.replace("--GTCCAAATT--", "---TCCAAATT--"), rows, { references, templateOrdinal: 0 }), /edited or deleted/);
  const curated = old.replace("--GTCCNNNTT--", "--GTCNNNNTT--");
  const repaired = restoreLineageOuterFlanks(curated, rows, { references });
  assert.equal(repaired.changed, true);
  assert.match(repaired.fasta, />__germline_N_masked__\nACGTCNNNNTTAC/);
  assert.equal(restoreLineageOuterFlanks(repaired.fasta, rows, { references }).changed, false);
});

test("phyloHMM constrains truncated-template flanks to reference despite unanimous mutated tips", async () => {
  const rows = [member(0, true), member(1), member(2, false, "ATA")];
  const input = lineageInputFasta(rows, "closest", { references, templateOrdinal: 0 });
  const observed = prepareObservedOnlyAlignment(input.fasta, GERMLINE_OUTGROUP);
  const options = defaultPhyloUcaOptions();
  options.search = { ...options.search, inferenceMode: "maximum-likelihood", fullHmmEdges: 0, edgeGridPoints: 2, branchGridPoints: 2, localRefinementRounds: 0, maximumUcaBranchLength: 0.1 };
  options.hmm.maximumDSegments = 0;
  const result = await inferPhyloUca({
    curatedAlignmentFasta: input.fasta, observedAlignmentFasta: observed.posteriorFasta,
    retainedColumns: observed.posteriorColumns, germlineGuideName: GERMLINE_OUTGROUP,
    templateOrdinal: 0,
    observedTreeNewick: "((read_0__1:0.01,read_1__2:0.01):0.01,read_2__3:0.02);",
    lineageRows: rows.map((row) => ({ ordinal: row.record.ordinal, sequenceId: row.record.sequenceId, locus: "IGH", values: row.values })),
    references, locus: "IGH", lineageLabel: "truncated", alignmentFingerprint: "truncated", frameOffset: 0, options,
  });
  assert.equal(result.mapAlignedSequence.slice(0, 2), "AC");
  assert.equal(result.templateOrdinal, 0);
  assert.equal(result.posterior[0].probabilities[0], 1);
  assert.equal(result.mapAlignedSequence.slice(-2), "AC");
  assert.match(result.observedAlignmentFasta, /TCGTCC/); // Observed SHM is preserved, not rewritten.
});
