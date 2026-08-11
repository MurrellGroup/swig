export type SegmentKey = "V" | "D" | "J" | "C";
export type LocusKey = "IGH" | "IGK" | "IGL" | "TRA" | "TRB" | "TRD" | "TRG";
export type ScopeKey = "BCR" | "TCR" | LocusKey;

export type Allele = [name: string, sequence: string];

export interface ReferenceLocus {
  V: Allele[];
  D?: Allele[];
  J: Allele[];
}

export interface ReferenceSpecies {
  name: string;
  loci: Partial<Record<LocusKey, ReferenceLocus>>;
}

export interface ReferencePack {
  source: string;
  release: string;
  retrieved: string;
  species: ReferenceSpecies[];
}

export interface CompiledReferences {
  V: string;
  D: string;
  J: string;
  C: string;
  counts: Record<SegmentKey, number>;
  loci: LocusKey[];
}

const BCR_LOCI: LocusKey[] = ["IGH", "IGK", "IGL"];
const TCR_LOCI: LocusKey[] = ["TRA", "TRB", "TRD", "TRG"];

export async function loadReferencePack(): Promise<ReferencePack> {
  const response = await fetch(`${import.meta.env.BASE_URL}references/imgt-202632-7.json.gz`);
  if (!response.ok) {
    throw new Error("The built-in reference library could not be loaded.");
  }
  const bytes = await response.arrayBuffer();
  const signature = new Uint8Array(bytes, 0, Math.min(2, bytes.byteLength));
  let json: string;
  if (signature[0] === 0x1f && signature[1] === 0x8b) {
    if (!("DecompressionStream" in globalThis)) {
      throw new Error("This browser cannot open the compressed reference library.");
    }
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
    json = await new Response(stream).text();
  } else {
    // Some static hosts transparently decode .gz responses via Content-Encoding.
    json = new TextDecoder().decode(bytes);
  }
  return JSON.parse(json) as ReferencePack;
}

export function availableScopes(species: ReferenceSpecies): ScopeKey[] {
  const loci = Object.keys(species.loci) as LocusKey[];
  const scopes: ScopeKey[] = [];
  if (loci.some((locus) => BCR_LOCI.includes(locus))) scopes.push("BCR");
  if (loci.some((locus) => TCR_LOCI.includes(locus))) scopes.push("TCR");
  return [...scopes, ...loci];
}

export function lociForScope(species: ReferenceSpecies, scope: ScopeKey): LocusKey[] {
  const requested = scope === "BCR" ? BCR_LOCI : scope === "TCR" ? TCR_LOCI : [scope];
  return requested.filter((locus) => species.loci[locus]);
}

export function allelesToFasta(alleles: Allele[]): string {
  return alleles.map(([name, sequence]) => `>${name}\n${sequence}\n`).join("");
}

export function compileReferences(
  species: ReferenceSpecies,
  scope: ScopeKey,
  overrides: Partial<Record<SegmentKey, string>> = {},
): CompiledReferences {
  const loci = lociForScope(species, scope);
  const segments: Record<SegmentKey, Allele[]> = { V: [], D: [], J: [], C: [] };
  for (const locus of loci) {
    const reference = species.loci[locus];
    if (!reference) continue;
    segments.V.push(...reference.V);
    segments.D.push(...(reference.D ?? []));
    segments.J.push(...reference.J);
  }
  const compiled: CompiledReferences = {
    V: overrides.V ?? allelesToFasta(segments.V),
    D: overrides.D ?? allelesToFasta(segments.D),
    J: overrides.J ?? allelesToFasta(segments.J),
    C: overrides.C ?? "",
    counts: {
      V: overrides.V ? countFastaRecords(overrides.V) : segments.V.length,
      D: overrides.D ? countFastaRecords(overrides.D) : segments.D.length,
      J: overrides.J ? countFastaRecords(overrides.J) : segments.J.length,
      C: overrides.C ? countFastaRecords(overrides.C) : 0,
    },
    loci,
  };
  return compiled;
}

export function countFastaRecords(value: string): number {
  return value.split(/\r?\n/).filter((line) => line.startsWith(">")).length;
}

export function makeDemoFasta(references: ReferenceSpecies, scope: ScopeKey): string {
  const loci = lociForScope(references, scope).slice(0, 3);
  const records: string[] = [];
  for (const locus of loci) {
    const reference = references.loci[locus];
    if (!reference) continue;
    const v = reference.V[0]?.[1];
    const d = reference.D?.[0]?.[1] ?? "";
    const j = reference.J[0]?.[1];
    if (!v || !j) continue;
    const sequence = `${v}AACCGG${d}TTG${j}`;
    records.push(`>demo_${locus}\n${sequence}`);
  }
  return `${records.join("\n")}\n`;
}
