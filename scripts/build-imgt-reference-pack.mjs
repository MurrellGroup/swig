#!/usr/bin/env node

import fs from "node:fs";
import path from "node:path";
import zlib from "node:zlib";

const [, , inputPath, release, retrieved, outputPath] = process.argv;

if (!inputPath || !release || !retrieved || !outputPath) {
  console.error(
    "Usage: node scripts/build-imgt-reference-pack.mjs " +
      "<IMGT-all-species.fasta> <release> <YYYY-MM-DD> <output.json.gz>",
  );
  process.exit(2);
}

if (!/^\d{4}-\d{2}-\d{2}$/.test(retrieved)) {
  console.error("The retrieval date must use YYYY-MM-DD.");
  process.exit(2);
}

const LOCI = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
const SEGMENTS = ["V", "D", "J", "C"];
const IMGT_V_GAPPED_ENDS = [78, 114, 165, 195, 312];

// Compact metadata layout shared with src/reference-pack.ts and the WASM
// adapter: frame, J CDR3 stop, five start/end pairs, provenance code.
// Provenance 1 is an exact delineation from an IMGT-gapped V-REGION;
// provenance 4 is a frame-validated J-motif annotation.
function metadata(frame = -1, cdr3Stop = -1, bounds = Array(10).fill(-1), source = 0) {
  return [frame, cdr3Stop, ...bounds, source];
}

function codonFrame(fields) {
  const value = Number.parseInt(fields[7] ?? "", 10);
  return value >= 1 && value <= 3 ? value - 1 : -1;
}

function ungappedLengthBefore(sequence, end) {
  return sequence.slice(0, end).replace(/[.\-]/g, "").length;
}

function nearestFrameCysEnd(sequence, predictedEnd, frame) {
  const candidates = [];
  const start = Math.max(0, predictedEnd - 60);
  const stop = Math.min(sequence.length - 2, predictedEnd + 60);
  for (let index = start; index <= stop; index += 1) {
    const codon = sequence.slice(index, index + 3);
    if (codon !== "TGT" && codon !== "TGC") continue;
    if (frame >= 0 && ((index - frame) % 3 + 3) % 3 !== 0) continue;
    candidates.push(index + 3);
  }
  candidates.sort((left, right) => Math.abs(left - predictedEnd) - Math.abs(right - predictedEnd) || right - left);
  return candidates[0];
}

function imgtVMetadata(rawSequence, fields) {
  if (rawSequence.length < IMGT_V_GAPPED_ENDS.at(-1)) return undefined;
  const ends = IMGT_V_GAPPED_ENDS.map((end) => ungappedLengthBefore(rawSequence, end));
  const sequence = rawSequence.replace(/[.\-]/g, "");
  const anchorEnd = nearestFrameCysEnd(sequence, ends[4], codonFrame(fields));
  if (anchorEnd && anchorEnd > ends[3]) ends[4] = anchorEnd;
  const bounds = [0, ends[0], ends[0], ends[1], ends[1], ends[2], ends[2], ends[3], ends[3], ends[4]];
  if (bounds.some((value, index) => index && value < bounds[index - 1])) return undefined;
  return metadata(codonFrame(fields), -1, bounds, 1);
}

function isAnchorCodon(codon) {
  return codon === "TGG" || codon === "TTT" || codon === "TTC";
}

function jMetadata(sequence, fields) {
  const frame = codonFrame(fields);
  if (frame < 0) return undefined;
  const candidates = [];
  for (let index = frame; index + 5 < sequence.length; index += 3) {
    if (isAnchorCodon(sequence.slice(index, index + 3)) && /^GG[ACGT]$/.test(sequence.slice(index + 3, index + 6))) {
      candidates.push(index);
    }
  }
  if (!candidates.length) return metadata(frame, -1, Array(10).fill(-1), 0);
  return metadata(frame, candidates[0] - 1, Array(10).fill(-1), 4);
}

function compareText(a, b) {
  return a < b ? -1 : a > b ? 1 : 0;
}

function parseFasta(text) {
  const records = [];
  let header = "";
  let sequence = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (header) records.push([header, sequence.join("")]);
      header = line.slice(1);
      sequence = [];
    } else if (header && line.trim()) {
      sequence.push(line.trim());
    }
  }
  if (header) records.push([header, sequence.join("")]);
  return records;
}

function normalizedSequence(rawSequence) {
  const gapped = rawSequence
    .toUpperCase()
    .replaceAll("U", "T")
    .replace(/\s/g, "");
  return { gapped, ungapped: gapped.replace(/[.\-]/g, "") };
}

// IMGT/GENE-DB exports light-chain constants as one C-REGION record, but
// heavy-chain and TCR constants as individual exon records.  SwiftIG expects
// one reference per allele, so assemble the coding exons in the biological
// order used by the IMGT export.  Membrane/UTR exons are intentionally omitted:
// the secreted coding path is long enough for isotype calling and does not mix
// mutually exclusive transcript tails.
function constantExonKind(locus, region) {
  if (region === "C-REGION") return "direct";
  if (locus === "IGH" && /^(?:CH(?:\d+(?:D\d*)?|X)|H\d*|CHS)(?:-(?:CH(?:\d+(?:D\d*)?|X)|H\d*|CHS))*$/.test(region)) {
    return "exon";
  }
  if (locus.startsWith("TR") && /^EX(?:[1-3](?:[A-Z]+)?|4)$/.test(region)) return "exon";
  return undefined;
}

function ensureLocusEntry(species, locus) {
  let speciesEntry = speciesMap.get(species);
  if (!speciesEntry) {
    speciesEntry = new Map();
    speciesMap.set(species, speciesEntry);
  }
  let locusEntry = speciesEntry.get(locus);
  if (!locusEntry) {
    locusEntry = { V: new Map(), D: new Map(), J: new Map(), C: new Map() };
    speciesEntry.set(locus, locusEntry);
  }
  return locusEntry;
}

function retainAllele(species, locus, segment, gene, sequence, annotation, score = sequence.length) {
  const target = ensureLocusEntry(species, locus)[segment];
  const previous = target.get(gene);
  if (!previous || score > previous[2] || (score === previous[2] && sequence.length > previous[0].length)) {
    target.set(gene, [sequence, annotation, score]);
  }
}

const input = fs.readFileSync(inputPath, "utf8");
const speciesMap = new Map();
const constantGroups = new Map();

for (const [order, [header, rawSequence]] of parseFasta(input).entries()) {
  const fields = header.split("|");
  const gene = fields[1]?.trim();
  const species = fields[2]?.trim();
  const region = fields[4]?.trim();
  const segmentMatch = /^([VDJC])-REGION$/.exec(region ?? "");
  const locus = LOCI.find((candidate) => gene?.toUpperCase().startsWith(candidate));
  if (!gene || !species || !locus || !region) continue;

  const { gapped: gappedSequence, ungapped: sequence } = normalizedSequence(rawSequence);
  if (!sequence) continue;

  const constantKind = constantExonKind(locus, region);
  if (!segmentMatch && constantKind === "exon") {
    const key = `${species}\u0000${locus}\u0000${gene}`;
    let group = constantGroups.get(key);
    if (!group) {
      group = { species, locus, gene, parts: new Map() };
      constantGroups.set(key, group);
    }
    const previous = group.parts.get(region);
    if (!previous || sequence.length > previous.sequence.length) {
      group.parts.set(region, { order, sequence, fields });
    }
    continue;
  }
  if (!segmentMatch) continue;

  let annotation;
  if (segmentMatch[1] === "V") annotation = imgtVMetadata(gappedSequence, fields);
  else if (segmentMatch[1] === "J") annotation = jMetadata(sequence, fields);
  else {
    const frame = codonFrame(fields);
    if (frame >= 0) annotation = metadata(frame, -1, Array(10).fill(-1), 1);
  }

  const segment = segmentMatch[1];
  retainAllele(species, locus, segment, gene, sequence, annotation);
}

for (const group of constantGroups.values()) {
  const parts = [...group.parts.values()].sort((left, right) => left.order - right.order);
  if (!parts.length) continue;
  const sequence = parts.map((part) => part.sequence).join("");
  const frame = parts.map((part) => codonFrame(part.fields)).find((value) => value >= 0) ?? -1;
  const annotation = frame >= 0 ? metadata(frame, -1, Array(10).fill(-1), 1) : undefined;
  const hasFirstDomain = group.locus === "IGH"
    ? group.parts.has("CH1")
    : group.parts.has("EX1");
  // Prefer a complete multi-exon path over a longer partial source when the
  // same allele appears more than once in the IMGT export.
  const score = (hasFirstDomain ? 1_000_000_000 : 0) + parts.length * 1_000_000 + sequence.length;
  retainAllele(group.species, group.locus, "C", group.gene, sequence, annotation, score);
}

const species = [];
for (const speciesName of [...speciesMap.keys()].sort(compareText)) {
  const sourceLoci = speciesMap.get(speciesName);
  const loci = {};
  for (const locus of LOCI) {
    const source = sourceLoci.get(locus);
    if (!source || !source.V.size || !source.J.size) continue;
    const reference = {};
    for (const segment of SEGMENTS) {
      const alleles = [...source[segment].entries()]
        .sort(([a], [b]) => compareText(a, b))
        .map(([name, [sequence, annotation]]) => annotation ? [name, sequence, annotation] : [name, sequence]);
      if (alleles.length) reference[segment] = alleles;
    }
    loci[locus] = reference;
  }
  if (Object.keys(loci).length) species.push({ name: speciesName, loci });
}

const pack = {
  schemaVersion: 2,
  source: "IMGT/GENE-DB",
  release,
  retrieved,
  species,
};
const json = `${JSON.stringify(pack)}\n`;
const compressed = zlib.gzipSync(Buffer.from(json), { level: 9, mtime: 0 });
fs.mkdirSync(path.dirname(outputPath), { recursive: true });
fs.writeFileSync(outputPath, compressed);

const alleleCount = species.reduce(
  (total, entry) => total + Object.values(entry.loci).reduce(
    (locusTotal, locus) => locusTotal + SEGMENTS.reduce(
      (segmentTotal, segment) => segmentTotal + (locus[segment]?.length ?? 0),
      0,
    ),
    0,
  ),
  0,
);
console.log(`Wrote ${species.length} species/strain sets and ${alleleCount} alleles to ${outputPath}`);
