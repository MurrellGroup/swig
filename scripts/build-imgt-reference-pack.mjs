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

function imgtVMetadata(rawSequence, fields) {
  if (rawSequence.length < IMGT_V_GAPPED_ENDS.at(-1)) return undefined;
  const ends = IMGT_V_GAPPED_ENDS.map((end) => ungappedLengthBefore(rawSequence, end));
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

const input = fs.readFileSync(inputPath, "utf8");
const speciesMap = new Map();

for (const [header, rawSequence] of parseFasta(input)) {
  const fields = header.split("|");
  const gene = fields[1]?.trim();
  const species = fields[2]?.trim();
  const region = fields[4]?.trim();
  const segmentMatch = /^([VDJC])-REGION$/.exec(region ?? "");
  const locus = LOCI.find((candidate) => gene?.toUpperCase().startsWith(candidate));
  if (!gene || !species || !segmentMatch || !locus) continue;

  const gappedSequence = rawSequence
    .toUpperCase()
    .replaceAll("U", "T")
    .replace(/\s/g, "");
  const sequence = gappedSequence.replace(/[.\-]/g, "");
  if (!sequence) continue;

  let annotation;
  if (segmentMatch[1] === "V") annotation = imgtVMetadata(gappedSequence, fields);
  else if (segmentMatch[1] === "J") annotation = jMetadata(sequence, fields);
  else {
    const frame = codonFrame(fields);
    if (frame >= 0) annotation = metadata(frame, -1, Array(10).fill(-1), 1);
  }

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
  const segment = segmentMatch[1];
  const previous = locusEntry[segment].get(gene);
  if (!previous || sequence.length > previous[0].length) {
    locusEntry[segment].set(gene, [sequence, annotation]);
  }
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
