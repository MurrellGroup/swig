export type ProbabilityLogoAlphabet = "nucleotide" | "amino-acid" | "custom";

export interface ProbabilityLogoEntry {
  symbol: string;
  probability: number;
  color?: string;
}

export interface ProbabilityLogoColumn {
  label: string;
  entries: ProbabilityLogoEntry[];
  title?: string;
}

export const NUCLEOTIDE_LOGO_SYMBOLS = ["A", "C", "G", "T", "-"] as const;
export const AMINO_ACID_LOGO_SYMBOLS = [..."ACDEFGHIKLMNPQRSTVWY", "*", "X", "-"] as const;

/**
 * Normalize non-negative finite masses and repair the final floating-point
 * residual on the largest entry. Logo heights and exported posteriors can
 * therefore share one normalization rule.
 */
export function normalizeProbabilityVector(values: readonly number[]): number[] {
  if (!values.length) return [];
  const safe = values.map((value) => Number.isFinite(value) && value > 0 ? value : 0);
  const total = safe.reduce((sum, value) => sum + value, 0);
  if (!(total > 0)) {
    const uniform = 1 / safe.length;
    const result = safe.map(() => uniform);
    result[result.length - 1] = 1 - result.slice(0, -1).reduce((sum, value) => sum + value, 0);
    return result;
  }
  const result = safe.map((value) => value / total);
  let largest = 0;
  for (let index = 1; index < result.length; index += 1) if (result[index] > result[largest]) largest = index;
  const residual = 1 - result.reduce((sum, value) => sum + value, 0);
  result[largest] = Math.max(0, result[largest] + residual);
  return result;
}

export function normalizedLogoEntries(entries: readonly ProbabilityLogoEntry[]): ProbabilityLogoEntry[] {
  const probabilities = normalizeProbabilityVector(entries.map((entry) => entry.probability));
  return entries.map((entry, index) => ({ ...entry, probability: probabilities[index] }));
}

export function probabilityLogoColor(symbol: string, alphabet: ProbabilityLogoAlphabet): string {
  const value = symbol.toUpperCase();
  if (alphabet === "nucleotide") return ({
    A: "#2d7d46", C: "#2e63a6", G: "#c58416", T: "#b84236", "-": "#737d79",
  } as Record<string, string>)[value] ?? "#5c6763";
  if (alphabet === "amino-acid") {
    if ("AVLIMFWY".includes(value)) return "#386fa4";
    if ("STNQC".includes(value)) return "#2f8b62";
    if ("KRH".includes(value)) return "#bd433e";
    if ("DE".includes(value)) return "#94539a";
    if ("GP".includes(value)) return "#c07822";
    if (value === "*") return "#222a27";
    if (value === "X" || value === "-") return "#737d79";
  }
  return "#385d53";
}
