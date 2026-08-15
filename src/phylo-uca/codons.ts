/**
 * Codon-state utilities for the phylogenetic UCA module.
 *
 * A codon state is an ordered triple over A/C/G/T/gap. The 125-state order is
 * base-5 lexicographic order, so state (a,b,c) has index 25a + 5b + c.
 */

export const PHYLO_UCA_CODON_CHARACTERS = ["A", "C", "G", "T", "-"] as const;
export const PHYLO_UCA_CODON_STATE_COUNT = 125;

const GENETIC_CODE: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S", TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
  CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P", CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T", AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A", GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

export const PHYLO_UCA_CODON_SYMBOLS: readonly string[] = Array.from(
  { length: PHYLO_UCA_CODON_STATE_COUNT },
  (_, index) => {
    const first = Math.floor(index / 25);
    const second = Math.floor(index / 5) % 5;
    const third = index % 5;
    return `${PHYLO_UCA_CODON_CHARACTERS[first]}${PHYLO_UCA_CODON_CHARACTERS[second]}${PHYLO_UCA_CODON_CHARACTERS[third]}`;
  },
);

export function phyloUcaCodonStateIndex(first: number, second: number, third: number): number {
  return first * 25 + second * 5 + third;
}

/** Complete gap codons are alignment gaps; partly gapped codons are unresolved. */
export function translatePhyloUcaCodon(codon: string): string {
  const normalized = codon.toUpperCase().replaceAll("U", "T").replaceAll(".", "-");
  if (normalized === "---") return "-";
  if (normalized.length !== 3 || normalized.includes("-") || /[^ACGT]/.test(normalized)) return "X";
  return GENETIC_CODE[normalized] ?? "X";
}

export function translatePhyloUcaCodonState(index: number): string {
  return translatePhyloUcaCodon(PHYLO_UCA_CODON_SYMBOLS[index] ?? "NNN");
}
