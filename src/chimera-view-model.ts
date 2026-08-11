export type ChimeraMatchCategory = "parent_a" | "parent_b" | "neutral";

export function chimeraVisiblePositions(query: string, parentA: string, parentB: string): number[] {
  return Array.from({ length: Math.max(query.length, parentA.length, parentB.length) }, (_, index) => index)
    .filter((index) => !(query[index] === "-" && parentA[index] === "-" && parentB[index] === "-"));
}

export function classifyChimeraQuerySite(query: string, parentA: string, parentB: string): ChimeraMatchCategory {
  const value = query.toUpperCase();
  if (!/[ACGT]/.test(value)) return "neutral";
  const matchesA = value === parentA.toUpperCase();
  const matchesB = value === parentB.toUpperCase();
  if (matchesA && !matchesB) return "parent_a";
  if (matchesB && !matchesA) return "parent_b";
  return "neutral";
}
