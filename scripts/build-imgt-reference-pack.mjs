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

  const sequence = rawSequence
    .toUpperCase()
    .replaceAll("U", "T")
    .replace(/[.\-\s]/g, "");
  if (!sequence) continue;

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
  if (!previous || sequence.length > previous.length) {
    locusEntry[segment].set(gene, sequence);
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
      const alleles = [...source[segment].entries()].sort(([a], [b]) => compareText(a, b));
      if (alleles.length) reference[segment] = alleles;
    }
    loci[locus] = reference;
  }
  if (Object.keys(loci).length) species.push({ name: speciesName, loci });
}

const pack = {
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
