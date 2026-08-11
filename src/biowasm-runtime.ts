import Aioli from "@biowasm/aioli";

import { projectCodonAlignment } from "./alignment-model";
import { prepareFastTreeInput } from "./fasttree-input";
import { parseFasta, type FastaRecord } from "./post-analysis-core";
import { extractNewick } from "./phylogeny";

interface AioliRuntime {
  write(options: { path: string; buffer: Uint8Array }): Promise<void>;
  exec(command: string): Promise<string>;
}

export interface FastTreeRun {
  newick: string;
  inputFasta: string;
  alignmentFasta: string;
  command: string;
  rows: number;
  columns: number;
  fingerprint: string;
}

let runtime: Promise<AioliRuntime> | null = null;

async function tools(): Promise<AioliRuntime> {
  runtime ??= Promise.resolve(new Aioli(["kalign/3.3.1", "fasttree/2.1.11"]) as unknown as Promise<AioliRuntime>);
  return runtime;
}

function safeFasta(records: FastaRecord[], preserveGaps = false): { fasta: string; names: string[] } {
  return {
    fasta: records.map((record, index) => `>${index}\n${preserveGaps ? record.sequence : record.sequence.replaceAll("-", "")}`).join("\n") + "\n",
    names: records.map((record) => record.name),
  };
}

function restoreNames(fasta: string, names: string[]): string {
  const records = parseFasta(fasta, true);
  return records.map((record) => {
    const index = Number(record.name);
    return `>${Number.isInteger(index) && names[index] ? names[index] : record.name}\n${record.sequence}`;
  }).join("\n") + "\n";
}

async function writeInput(cli: AioliRuntime, path: string, value: string): Promise<string> {
  // Aioli.mount() retains every Blob it has ever mounted. Direct w+ writes
  // truncate this browser-local MEMFS path on every invocation, so a corrected
  // alignment can never resolve to bytes from an earlier run.
  await cli.write({ path, buffer: new TextEncoder().encode(value) });
  return path;
}

const CODONS: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S",
  TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
  CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

function translate(sequence: string, frame: number): string {
  const clean = sequence.replaceAll("-", "").toUpperCase().replaceAll("U", "T");
  let aa = "";
  for (let index = frame; index + 2 < clean.length; index += 3) {
    const codon = clean.slice(index, index + 3);
    aa += CODONS[codon] ?? "X";
  }
  return aa;
}

export async function runKalign(fasta: string): Promise<string> {
  const records = parseFasta(fasta, true);
  if (records.length < 2) throw new Error("Kalign needs at least two sequences.");
  const input = safeFasta(records);
  const cli = await tools();
  const path = await writeInput(cli, "/shared/data/swig_kalign_nt.fa", input.fasta);
  const output = await cli.exec(`kalign ${path} -f fasta`);
  if (!output.includes(">")) throw new Error("Kalign did not return a FASTA alignment.");
  return restoreNames(output, input.names);
}

export async function runCodonAwareKalign(fasta: string, frames?: number[]): Promise<string> {
  const records = parseFasta(fasta, true);
  if (records.length < 2) throw new Error("Codon-aware Kalign needs at least two sequences.");
  const normalizedFrames = records.map((_, index) => Math.max(0, Math.min(2, frames?.[index] ?? 0)));
  const aminoAcids = records.map((record, index) => ({ name: String(index), sequence: translate(record.sequence, normalizedFrames[index]) }));
  const cli = await tools();
  const path = await writeInput(cli, "/shared/data/swig_kalign_codon_aa.fa", aminoAcids.map((record) => `>${record.name}\n${record.sequence}`).join("\n") + "\n");
  const output = await cli.exec(`kalign ${path} -f fasta`);
  const aligned = parseFasta(output, true);
  const byIndex = new Map(aligned.map((record) => [Number(record.name), record.sequence]));
  const projected = records.map((record, index) => ({
    name: record.name,
    sequence: projectCodonAlignment(record.sequence, byIndex.get(index) ?? aminoAcids[index].sequence, normalizedFrames[index]),
  }));
  const maximum = Math.max(...projected.map((record) => record.sequence.length));
  return projected.map((record) => `>${record.name}\n${record.sequence.padEnd(maximum, "-")}`).join("\n") + "\n";
}

export async function runFastTree(alignedFasta: string, model: "gtr" | "jc" = "gtr", fast = false): Promise<FastTreeRun> {
  const input = prepareFastTreeInput(alignedFasta, model, fast);
  const cli = await tools();
  const path = await writeInput(cli, "/shared/data/swig_fasttree_input.fa", input.inputFasta);
  let output = extractNewick(await cli.exec(`fasttree ${input.flags} ${path}`));
  input.names.forEach((name, index) => {
    const safeName = name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") || `tip_${index + 1}`;
    output = output.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), (_match, prefix: string) => `${prefix}${safeName}`);
  });
  return {
    newick: output.trim(),
    inputFasta: input.inputFasta,
    alignmentFasta: input.alignmentFasta,
    command: input.command,
    rows: input.rows,
    columns: input.columns,
    fingerprint: input.fingerprint,
  };
}

export async function warmBiowasm(): Promise<void> {
  await tools();
}
