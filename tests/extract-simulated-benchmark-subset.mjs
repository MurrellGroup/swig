import fs from "node:fs";

const [queriesPath, callsPath, modulusText, residueText, outputPrefix] = process.argv.slice(2);
if (!outputPrefix) throw new Error("Usage: node extract-simulated-benchmark-subset.mjs QUERIES CALLS MODULUS RESIDUE OUTPUT_PREFIX");
const modulus = Number(modulusText);
const residue = Number(residueText);

function fnv1a(value) {
  let hash = 0x811c9dc5;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193) >>> 0;
  }
  return hash;
}

const selected = (id) => fnv1a(id) % modulus === residue;
const fasta = fs.readFileSync(queriesPath, "utf8");
const fastaOutput = [];
let header = "";
let sequence = "";
const commit = () => {
  if (!header) return;
  const id = header.slice(1).trim().split(/\s+/, 1)[0];
  if (selected(id)) fastaOutput.push(`${header}\n${sequence}\n`);
};
for (const line of fasta.split(/\r?\n/)) {
  if (line.startsWith(">")) {
    commit();
    header = line;
    sequence = "";
  } else sequence += line.trim();
}
commit();

const callLines = fs.readFileSync(callsPath, "utf8").trimEnd().split(/\r?\n/);
const callHeader = callLines.shift();
const selectedCalls = callLines.filter((line) => selected(line.slice(0, line.indexOf("\t"))));
fs.writeFileSync(`${outputPrefix}.fasta`, fastaOutput.join(""));
fs.writeFileSync(`${outputPrefix}.tsv`, `${callHeader}\n${selectedCalls.join("\n")}\n`);
console.log(JSON.stringify({ records: selectedCalls.length, fasta: `${outputPrefix}.fasta`, calls: `${outputPrefix}.tsv` }));
