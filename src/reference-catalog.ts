import type { LocusKey, SegmentKey } from "./reference-pack";

export interface ReferenceCollection {
  id: string;
  name: string;
  provider: string;
  version: string;
  speciesPrefixes: string[];
  locus: LocusKey;
  summary: string;
  sourceUrl: string;
  citationUrl: string;
  terms?: { label: string; url: string };
  segments: Partial<Record<SegmentKey, ReferenceSource>>;
}

interface ReferenceSource {
  url: string;
  headers?: Record<string, string>;
  local?: boolean;
}

const KI_TCR_ROOT = "https://gkhlab.gitlab.io/tcr/sequences";
const KIARVA_API = "https://kiarva.scilifelab.se/api/fasta/genomic";
const KIARVA_HEADERS = { "X-api-key": "kiarvafrontend" };

function kiTcrCollection(locus: LocusKey, segments: SegmentKey[]): ReferenceCollection {
  return {
    id: `ki-human-${locus.toLowerCase()}`,
    name: `KI human ${locus} germline set`,
    provider: "Karlsson Hedestam laboratory, Karolinska Institutet",
    version: "Corcoran et al. 2023",
    speciesPrefixes: ["Homo sapiens"],
    locus,
    summary: `Human ${locus} V(D)J reference sequences from the KI human TCR database.`,
    sourceUrl: "https://gkhlab.gitlab.io/tcr/sequences/",
    citationUrl: "https://doi.org/10.1016/j.immuni.2023.02.007",
    segments: Object.fromEntries(segments.map((segment) => [segment, {
      url: `${KI_TCR_ROOT}/${locus}${segment}.fasta`,
    }])) as Partial<Record<SegmentKey, ReferenceSource>>,
  };
}

function kimdbCollection(species: "Macaca mulatta" | "Macaca fascicularis"): ReferenceCollection {
  const slug = species.replace(" ", "_");
  return {
    id: `kimdb-${slug.toLowerCase()}`,
    name: `KIMDB 1.1 · ${species} IGH`,
    provider: "Karlsson Hedestam laboratory, Karolinska Institutet",
    version: "1.1",
    speciesPrefixes: [species],
    locus: "IGH",
    summary: `Curated ${species} heavy-chain V, D, and J germline sequences from KIMDB.`,
    sourceUrl: "http://kimdb.gkhlab.se/",
    citationUrl: "https://pubmed.ncbi.nlm.nih.gov/33484642/",
    segments: Object.fromEntries((["V", "D", "J"] as SegmentKey[]).map((segment) => [segment, {
      url: `references/kimdb-1.1/${slug}/IGH/${segment}.fasta`,
      local: true,
    }])) as Partial<Record<SegmentKey, ReferenceSource>>,
  };
}

export const REFERENCE_COLLECTIONS: ReferenceCollection[] = [
  {
    id: "kiarva-human-igh",
    name: "KIARVA · human IGH",
    provider: "Karlsson Hedestam laboratory, Karolinska Institutet",
    version: "current genomic release",
    speciesPrefixes: ["Homo sapiens"],
    locus: "IGH",
    summary: "Human population IGHV, IGHD, and IGHJ alleles inferred from 2,486 1000 Genomes Project cases.",
    sourceUrl: "https://kiarva.scilifelab.se/about",
    citationUrl: "https://kiarva.scilifelab.se/about",
    terms: { label: "CC BY-NC 4.0", url: "https://kiarva.scilifelab.se/license" },
    segments: {
      V: { url: `${KIARVA_API}?file_name=IGHV`, headers: KIARVA_HEADERS },
      D: { url: `${KIARVA_API}?file_name=IGHD`, headers: KIARVA_HEADERS },
      J: { url: `${KIARVA_API}?file_name=IGHJ`, headers: KIARVA_HEADERS },
    },
  },
  kiTcrCollection("TRA", ["V", "J"]),
  kiTcrCollection("TRB", ["V", "D", "J"]),
  kiTcrCollection("TRD", ["V", "D", "J"]),
  kiTcrCollection("TRG", ["V", "J"]),
  kimdbCollection("Macaca mulatta"),
  kimdbCollection("Macaca fascicularis"),
];

export function collectionsFor(species: string, locus: LocusKey | null): ReferenceCollection[] {
  if (!locus) return [];
  return REFERENCE_COLLECTIONS.filter((collection) => collection.locus === locus
    && collection.speciesPrefixes.some((prefix) => species === prefix || species.startsWith(`${prefix}_`)));
}

export async function loadCollectionSegment(
  collection: ReferenceCollection,
  segment: SegmentKey,
): Promise<string> {
  const source = collection.segments[segment];
  if (!source) throw new Error(`${collection.name} does not contain a ${segment} set.`);
  const url = source.local ? `${import.meta.env.BASE_URL}${source.url}` : source.url;
  const response = await fetch(url, { headers: source.headers });
  if (!response.ok) throw new Error(`${collection.name} ${segment} download failed (HTTP ${response.status}).`);
  const text = await response.text();
  if (!text.trimStart().startsWith(">")) throw new Error(`${collection.name} returned data that are not FASTA.`);
  return text;
}
