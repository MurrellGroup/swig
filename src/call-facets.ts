import type { FacetValue } from "./result-store";

/**
 * Expand complete AIRR call strings into discoverable gene and allele tokens.
 * Ambiguous comma-separated assignments contribute to each member, but are
 * never exposed as a single unwieldy picker option.
 */
export function callFacetItems(values: FacetValue[]): FacetValue[] {
  const counts = new Map<string, number>();
  for (const facet of values) {
    for (const rawCall of facet.value.split(",")) {
      const allele = rawCall.trim();
      if (!allele) continue;
      counts.set(allele, (counts.get(allele) ?? 0) + facet.count);
      const gene = allele.replace(/\*.*$/, "");
      if (gene && gene !== allele) counts.set(gene, (counts.get(gene) ?? 0) + facet.count);
    }
  }
  return [...counts.entries()]
    .map(([value, count]) => ({ value, count }))
    .sort((left, right) => left.value.localeCompare(right.value, undefined, { numeric: true }));
}
