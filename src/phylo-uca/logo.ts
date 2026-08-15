import {
  AMINO_ACID_LOGO_SYMBOLS,
  NUCLEOTIDE_LOGO_SYMBOLS,
  normalizeProbabilityVector,
  probabilityLogoColor,
  type ProbabilityLogoColumn,
} from "../probability-logo.ts";
import {
  PHYLO_UCA_CODON_STATE_COUNT,
  PHYLO_UCA_CODON_SYMBOLS,
  translatePhyloUcaCodonState,
} from "./codons.ts";
import type { PhyloUcaCodonPosterior, PhyloUcaSitePosterior } from "./types.ts";

function summary(symbols: readonly string[], probabilities: readonly number[], maximum = 16): string {
  return symbols.map((symbol, index) => ({ symbol, probability: probabilities[index] ?? 0 }))
    .filter((entry) => entry.probability > 0)
    .sort((left, right) => right.probability - left.probability || left.symbol.localeCompare(right.symbol))
    .slice(0, maximum)
    .map((entry) => `${entry.symbol} ${(entry.probability * 100).toFixed(3)}%`)
    .join(" · ");
}

export function nucleotideUcaLogoColumns(posterior: readonly PhyloUcaSitePosterior[]): ProbabilityLogoColumn[] {
  return posterior.map((site) => {
    const probabilities = normalizeProbabilityVector(site.probabilities);
    return {
      label: String(site.alignmentColumn),
      title: `Alignment column ${site.alignmentColumn} · ${summary(NUCLEOTIDE_LOGO_SYMBOLS, probabilities)}`,
      entries: NUCLEOTIDE_LOGO_SYMBOLS.map((symbol, index) => ({ symbol, probability: probabilities[index] })),
    };
  });
}

export function codonUcaLogoColumns(posterior: readonly PhyloUcaCodonPosterior[]): ProbabilityLogoColumn[] {
  return posterior.map((codon) => {
    const probabilities = normalizeProbabilityVector(codon.probabilities);
    return {
      label: String(codon.codonIndex),
      title: `Codon ${codon.codonIndex} · alignment columns ${codon.alignmentColumns.join(", ")} · ${summary(PHYLO_UCA_CODON_SYMBOLS, probabilities)}`,
      entries: PHYLO_UCA_CODON_SYMBOLS.map((symbol, index) => ({
        symbol,
        probability: probabilities[index],
        color: probabilityLogoColor(translatePhyloUcaCodonState(index), "amino-acid"),
      })),
    };
  });
}

/** Sum the exact 125-state codon posterior over synonymous codon states. */
export function aminoAcidProbabilitiesForCodon(codon: PhyloUcaCodonPosterior): number[] {
  const codonProbabilities = normalizeProbabilityVector(codon.probabilities);
  const masses = new Map<string, number>(AMINO_ACID_LOGO_SYMBOLS.map((symbol) => [symbol, 0]));
  for (let state = 0; state < PHYLO_UCA_CODON_STATE_COUNT; state += 1) {
    const aminoAcid = translatePhyloUcaCodonState(state);
    masses.set(aminoAcid, (masses.get(aminoAcid) ?? 0) + codonProbabilities[state]);
  }
  return normalizeProbabilityVector(AMINO_ACID_LOGO_SYMBOLS.map((symbol) => masses.get(symbol) ?? 0));
}

/** Exact AA marginals obtained by summing, never multiplying, codon states. */
export function aminoAcidUcaLogoColumns(posterior: readonly PhyloUcaCodonPosterior[]): ProbabilityLogoColumn[] {
  return posterior.map((codon) => {
    const probabilities = aminoAcidProbabilitiesForCodon(codon);
    return {
      label: String(codon.codonIndex),
      title: `Codon ${codon.codonIndex} · alignment columns ${codon.alignmentColumns.join(", ")} · ${summary(AMINO_ACID_LOGO_SYMBOLS, probabilities)}`,
      entries: AMINO_ACID_LOGO_SYMBOLS.map((symbol, index) => ({ symbol, probability: probabilities[index] })),
    };
  });
}
