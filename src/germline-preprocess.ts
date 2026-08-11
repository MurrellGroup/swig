export type GermlineSegment = "V" | "D" | "J" | "C";
export type GermlineLocus = "IGH" | "IGK" | "IGL" | "TRA" | "TRB" | "TRD" | "TRG";

// frame, J CDR3 stop, FWR1/CDR1/FWR2/CDR2/FWR3 start/end pairs, source.
export type CompactMetadata = [
  frame: number,
  cdr3Stop: number,
  fwr1Start: number,
  fwr1End: number,
  cdr1Start: number,
  cdr1End: number,
  fwr2Start: number,
  fwr2End: number,
  cdr2Start: number,
  cdr2End: number,
  fwr3Start: number,
  fwr3End: number,
  source: number,
];

export type MetadataAllele = [name: string, sequence: string, metadata?: CompactMetadata];

export interface GermlinePreprocessReport {
  fasta: string;
  count: number;
  annotated: number;
  unannotated: number;
  exactImgt: number;
  transferred: number;
  motifValidated: number;
  ambiguousBases: number;
  loci: GermlineLocus[];
  warnings: string[];
}

interface ParsedRecord {
  header: string;
  name: string;
  rawSequence: string;
  sequence: string;
  locus: GermlineLocus;
  metadata?: CompactMetadata;
}

interface Alignment {
  query: string;
  reference: string;
  identity: number;
}

const LOCI: GermlineLocus[] = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
const IMGT_V_GAPPED_ENDS = [78, 114, 165, 195, 312] as const;
const EMPTY_BOUNDS = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1] as const;
const SOURCE_IMGT_GAPPED = 1;
const SOURCE_AIRR_C = 2;
const SOURCE_TRANSFERRED_IMGT = 3;
const SOURCE_VALIDATED_J_MOTIF = 4;
const SOURCE_PROVIDED = 5;
const SOURCE_TRANSFERRED_J = 6;

export const METADATA_SOURCE_LABELS: Record<number, string> = {
  [SOURCE_IMGT_GAPPED]: "IMGT-gapped delineation",
  [SOURCE_AIRR_C]: "AIRR-C annotation",
  [SOURCE_TRANSFERRED_IMGT]: "validated IMGT-boundary transfer",
  [SOURCE_VALIDATED_J_MOTIF]: "frame-validated J motif",
  [SOURCE_PROVIDED]: "provided annotation",
  [SOURCE_TRANSFERRED_J]: "validated J-anchor transfer",
};

function positiveModulo(value: number, modulus: number): number {
  return ((value % modulus) + modulus) % modulus;
}

function compactMetadata(
  frame = -1,
  cdr3Stop = -1,
  bounds: readonly number[] = EMPTY_BOUNDS,
  source = 0,
): CompactMetadata {
  return [frame, cdr3Stop, ...bounds, source] as CompactMetadata;
}

function inferLocus(value: string): GermlineLocus | null {
  const upper = value.toUpperCase();
  return LOCI.find((locus) => upper.includes(locus)) ?? null;
}

function inferSegment(name: string): GermlineSegment | null {
  const match = /^(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)([VDJC])/.exec(name.toUpperCase());
  return match ? match[1] as GermlineSegment : null;
}

function germlineName(header: string): string {
  const identifier = header.trim().split(/\s+/, 1)[0] ?? "";
  const fields = identifier.split("|");
  return fields.find((field) => inferLocus(field) && field.includes("*"))
    ?? fields.find((field) => inferLocus(field))
    ?? identifier;
}

function normalizeSequence(raw: string): { sequence: string; ambiguous: number } {
  const compact = raw.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
  if (/[^ACGTNRYKMSWBDHV.\-]/.test(compact)) {
    throw new Error("Germline FASTA contains characters outside the IUPAC nucleotide alphabet.");
  }
  const withoutGaps = compact.replace(/[.\-]/g, "");
  const ambiguous = (withoutGaps.match(/[^ACGT]/g) ?? []).length;
  return { sequence: withoutGaps.replace(/[^ACGT]/g, "N"), ambiguous };
}

function parseFasta(text: string): Array<{ header: string; rawSequence: string }> {
  const records: Array<{ header: string; rawSequence: string }> = [];
  let header = "";
  let sequence: string[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (header) records.push({ header, rawSequence: sequence.join("") });
      header = line.slice(1).trim();
      sequence = [];
    } else if (header && line.trim()) {
      sequence.push(line.trim());
    } else if (!header && line.trim()) {
      throw new Error("Expected a FASTA header beginning with '>'.");
    }
  }
  if (header) records.push({ header, rawSequence: sequence.join("") });
  if (!records.length) throw new Error("No germline FASTA records were found.");
  return records;
}

function parseImgtFields(header: string): string[] | null {
  const identifier = header.split(/\s+/, 1)[0] ?? "";
  const fields = identifier.split("|");
  return fields.length >= 8 && fields[1] && fields[4] ? fields : null;
}

function imgtFrame(fields: string[] | null): number {
  const value = Number.parseInt(fields?.[7] ?? "", 10);
  return value >= 1 && value <= 3 ? value - 1 : -1;
}

function metadataFromHeader(header: string): CompactMetadata | undefined {
  const match = /(?:^|\s)SWIGMETA=([\d,\-]+)/i.exec(header);
  if (!match) return undefined;
  const values = match[1].split(",").map(Number);
  if (values.length !== 13 || values.some((value) => !Number.isInteger(value))) return undefined;
  return values as CompactMetadata;
}

function metadataHeader(metadata?: CompactMetadata): string {
  return metadata ? ` SWIGMETA=${metadata.join(",")}` : "";
}

function hasRegionMetadata(metadata?: CompactMetadata): boolean {
  return Boolean(metadata && metadata.slice(2, 12).every((value) => value >= 0));
}

function hasJMetadata(metadata?: CompactMetadata): boolean {
  return Boolean(metadata && metadata[0] >= 0 && metadata[1] >= 0);
}

function validateMetadata(metadata: CompactMetadata, sequenceLength: number, segment: GermlineSegment): boolean {
  if (metadata[0] < -1 || metadata[0] > 2 || metadata[1] < -1 || metadata[1] >= sequenceLength) return false;
  if (segment === "V" && hasRegionMetadata(metadata)) {
    const bounds = metadata.slice(2, 12);
    if (bounds.some((value, index) => value < 0 || value > sequenceLength || (index && value < bounds[index - 1]))) return false;
    for (let index = 0; index < bounds.length; index += 2) {
      if (bounds[index + 1] <= bounds[index]) return false;
    }
  }
  return true;
}

function imgtVMetadata(rawSequence: string, fields: string[] | null): CompactMetadata | undefined {
  if (fields?.[4]?.toUpperCase() !== "V-REGION") return undefined;
  const gapped = rawSequence.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
  if (gapped.length < IMGT_V_GAPPED_ENDS.at(-1)!) return undefined;
  const ends = IMGT_V_GAPPED_ENDS.map((end) => gapped.slice(0, end).replace(/[.\-]/g, "").length);
  const bounds = [0, ends[0], ends[0], ends[1], ends[1], ends[2], ends[2], ends[3], ends[3], ends[4]];
  const result = compactMetadata(imgtFrame(fields), -1, bounds, SOURCE_IMGT_GAPPED);
  return validateMetadata(result, gapped.replace(/[.\-]/g, "").length, "V") ? result : undefined;
}

function isJAnchor(sequence: string, index: number): boolean {
  const aromatic = sequence.slice(index, index + 3);
  return (aromatic === "TGG" || aromatic === "TTT" || aromatic === "TTC")
    && /^GG[ACGT]$/.test(sequence.slice(index + 3, index + 6));
}

function jMetadata(sequence: string, fields: string[] | null): CompactMetadata | undefined {
  const providedFrame = imgtFrame(fields);
  const candidates: number[] = [];
  const limit = Math.min(sequence.length - 6, 54);
  for (let index = Math.max(0, providedFrame); index <= limit; index += providedFrame >= 0 ? 3 : 1) {
    if (isJAnchor(sequence, index)) candidates.push(index);
  }
  if (!candidates.length || (providedFrame < 0 && candidates.length !== 1)) return undefined;
  const anchor = candidates[0];
  const frame = providedFrame >= 0 ? providedFrame : anchor % 3;
  if (anchor % 3 !== frame) return undefined;
  return compactMetadata(frame, anchor - 1, EMPTY_BOUNDS, SOURCE_VALIDATED_J_MOTIF);
}

function globalAlignment(query: string, reference: string): Alignment {
  const columns = reference.length + 1;
  const rows = query.length + 1;
  const scores = new Int16Array(columns);
  const next = new Int16Array(columns);
  const trace = new Uint8Array(rows * columns);
  for (let column = 1; column < columns; column += 1) {
    scores[column] = -3 * column;
    trace[column] = 2;
  }
  for (let row = 1; row < rows; row += 1) {
    next[0] = -3 * row;
    trace[row * columns] = 1;
    for (let column = 1; column < columns; column += 1) {
      const diagonal = scores[column - 1] + (query[row - 1] === reference[column - 1] ? 2 : -2);
      const up = scores[column] - 3;
      const left = next[column - 1] - 3;
      const best = Math.max(diagonal, up, left);
      next[column] = best;
      trace[row * columns + column] = best === diagonal ? 0 : best === up ? 1 : 2;
    }
    scores.set(next);
  }
  let row = query.length;
  let column = reference.length;
  const alignedQuery: string[] = [];
  const alignedReference: string[] = [];
  let matches = 0;
  while (row || column) {
    const direction = trace[row * columns + column];
    if (row && column && direction === 0) {
      const queryBase = query[--row];
      const referenceBase = reference[--column];
      alignedQuery.push(queryBase);
      alignedReference.push(referenceBase);
      if (queryBase === referenceBase) matches += 1;
    } else if (row && (direction === 1 || !column)) {
      alignedQuery.push(query[--row]);
      alignedReference.push("-");
    } else {
      alignedQuery.push("-");
      alignedReference.push(reference[--column]);
    }
  }
  alignedQuery.reverse();
  alignedReference.reverse();
  return {
    query: alignedQuery.join(""),
    reference: alignedReference.join(""),
    identity: matches / Math.max(1, alignedQuery.length),
  };
}

function canonicalAllele(name: string): string {
  return name.replace(/_S\d+$/i, "").toUpperCase();
}

function canonicalGene(name: string): string {
  return canonicalAllele(name).split("*", 1)[0];
}

function kmerSet(sequence: string, size = 9): Set<string> {
  const output = new Set<string>();
  for (let index = 0; index + size <= sequence.length; index += 1) {
    const kmer = sequence.slice(index, index + size);
    if (!kmer.includes("N")) output.add(kmer);
  }
  return output;
}

function templateCandidates(
  queryName: string,
  query: string,
  templates: MetadataAllele[],
  eligible: (metadata?: CompactMetadata) => boolean = hasRegionMetadata,
): MetadataAllele[] {
  // Functional/pseudogene labels are deliberately irrelevant here. Any locus-matched
  // template with a complete IMGT delineation may transfer coordinates after the
  // sequence-level validation below.
  const delineated = templates.filter((template) => eligible(template[2]));
  const alleleName = canonicalAllele(queryName);
  const geneName = canonicalGene(queryName);
  const exactAllele = delineated.filter(([name]) => canonicalAllele(name) === alleleName);
  if (exactAllele.length) return exactAllele;
  const exactGene = delineated.filter(([name]) => canonicalGene(name) === geneName);
  if (exactGene.length) return exactGene;
  const queryKmers = kmerSet(query);
  return delineated
    .map((template) => {
      let shared = 0;
      for (const kmer of kmerSet(template[1])) if (queryKmers.has(kmer)) shared += 1;
      return { template, shared };
    })
    .sort((a, b) => b.shared - a.shared)
    .slice(0, 6)
    .map(({ template }) => template);
}

function mapReferenceCoordinates(alignment: Alignment, referenceLength: number): number[] {
  const mapped = new Array<number>(referenceLength + 1).fill(-1);
  let queryPosition = 0;
  let referencePosition = 0;
  mapped[0] = 0;
  for (let column = 0; column < alignment.query.length; column += 1) {
    const queryBase = alignment.query[column];
    const referenceBase = alignment.reference[column];
    if (referenceBase === "-") {
      if (queryBase !== "-") queryPosition += 1;
      mapped[referencePosition] = queryPosition;
      continue;
    }
    if (mapped[referencePosition] < 0) mapped[referencePosition] = queryPosition;
    if (queryBase !== "-") queryPosition += 1;
    referencePosition += 1;
    mapped[referencePosition] = queryPosition;
  }
  return mapped;
}

function transferMetadata(name: string, sequence: string, templates: MetadataAllele[]): CompactMetadata | undefined {
  const candidates = templateCandidates(name, sequence, templates);
  let selected: { template: MetadataAllele; alignment: Alignment } | undefined;
  for (const template of candidates) {
    const alignment = globalAlignment(sequence, template[1]);
    if (!selected || alignment.identity > selected.alignment.identity) selected = { template, alignment };
  }
  if (!selected) return undefined;
  const named = canonicalGene(selected.template[0]) === canonicalGene(name);
  if (selected.alignment.identity < (named ? 0.80 : 0.72)) return undefined;

  const templateMetadata = selected.template[2]!;
  const templateBounds = templateMetadata.slice(2, 12);
  const mapped = mapReferenceCoordinates(selected.alignment, selected.template[1].length);
  const bounds = templateBounds.map((boundary) => mapped[boundary]);
  if (bounds.some((value, index) => value < 0 || value > sequence.length || (index && value < bounds[index - 1]))) return undefined;
  for (let index = 0; index < bounds.length; index += 2) {
    if (bounds[index + 1] <= bounds[index]) return undefined;
  }
  const templateFrame = templateMetadata[0] >= 0 ? templateMetadata[0] : templateBounds[0] % 3;
  const frame = positiveModulo(bounds[0] + templateFrame - templateBounds[0], 3);
  const result = compactMetadata(frame, -1, bounds, SOURCE_TRANSFERRED_IMGT);
  return validateMetadata(result, sequence.length, "V") ? result : undefined;
}

function transferJMetadata(name: string, sequence: string, templates: MetadataAllele[]): CompactMetadata | undefined {
  const candidates = templateCandidates(name, sequence, templates, hasJMetadata);
  let selected: { template: MetadataAllele; alignment: Alignment } | undefined;
  for (const template of candidates) {
    const alignment = globalAlignment(sequence, template[1]);
    if (!selected || alignment.identity > selected.alignment.identity) selected = { template, alignment };
  }
  if (!selected) return undefined;
  const named = canonicalGene(selected.template[0]) === canonicalGene(name);
  if (selected.alignment.identity < (named ? 0.75 : 0.68)) return undefined;
  const metadata = selected.template[2]!;
  const referenceAnchor = metadata[1] + 1;
  if (referenceAnchor < 0 || referenceAnchor + 6 > selected.template[1].length) return undefined;
  const mapped = mapReferenceCoordinates(selected.alignment, selected.template[1].length);
  const anchor = mapped[referenceAnchor];
  const anchorEnd = mapped[referenceAnchor + 6];
  if (anchor < 0 || anchorEnd - anchor !== 6 || !isJAnchor(sequence, anchor)) return undefined;
  const result = compactMetadata(anchor % 3, anchor - 1, EMPTY_BOUNDS, SOURCE_TRANSFERRED_J);
  return validateMetadata(result, sequence.length, "J") ? result : undefined;
}

function annotationPresent(segment: GermlineSegment, metadata?: CompactMetadata): boolean {
  if (segment === "V") return hasRegionMetadata(metadata);
  if (segment === "J") return hasJMetadata(metadata);
  return true;
}

export function preprocessGermlineFasta(
  text: string,
  segment: GermlineSegment,
  templates: MetadataAllele[] = [],
  allowedLoci: readonly GermlineLocus[] = LOCI,
): GermlinePreprocessReport {
  const inputRecords = parseFasta(text);
  const records: ParsedRecord[] = [];
  const seen = new Set<string>();
  const loci = new Set<GermlineLocus>();
  const warnings: string[] = [];
  let exactImgt = 0;
  let transferred = 0;
  let motifValidated = 0;
  let ambiguousBases = 0;

  for (const input of inputRecords) {
    const name = germlineName(input.header);
    if (!name) throw new Error(`A ${segment} germline record has an empty identifier.`);
    if (seen.has(name)) throw new Error(`Duplicate ${segment} germline identifier: ${name}.`);
    seen.add(name);
    const locus = inferLocus(`${name} ${input.header}`);
    if (!locus) throw new Error(`${name} does not identify one of the supported IG/TR loci.`);
    if (!allowedLoci.includes(locus)) {
      throw new Error(`${name} belongs to ${locus}, outside the selected ${allowedLoci.join("/")} search space.`);
    }
    const detectedSegment = inferSegment(name);
    if (detectedSegment && detectedSegment !== segment) {
      throw new Error(`${name} looks like a ${detectedSegment} gene, not a ${segment} gene.`);
    }
    const normalized = normalizeSequence(input.rawSequence);
    if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
    ambiguousBases += normalized.ambiguous;
    const fields = parseImgtFields(input.header);
    let metadata = metadataFromHeader(input.header);
    if (metadata && !validateMetadata(metadata, normalized.sequence.length, segment)) {
      throw new Error(`${name} contains invalid SWIGMETA coordinates.`);
    }
    if (metadata) metadata = [...metadata.slice(0, 12), metadata[12] || SOURCE_PROVIDED] as CompactMetadata;
    if (!metadata && segment === "V") metadata = imgtVMetadata(input.rawSequence, fields);
    if (!metadata && segment === "V" && templates.length) metadata = transferMetadata(name, normalized.sequence, templates);
    if (!metadata && segment === "J" && templates.length) metadata = transferJMetadata(name, normalized.sequence, templates);
    if (!metadata && segment === "J") metadata = jMetadata(normalized.sequence, fields);
    if (!metadata && segment === "D") {
      const frame = imgtFrame(fields);
      if (frame >= 0) metadata = compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_IMGT_GAPPED);
    }
    if (metadata?.[12] === SOURCE_IMGT_GAPPED) exactImgt += 1;
    if (metadata?.[12] === SOURCE_TRANSFERRED_IMGT || metadata?.[12] === SOURCE_TRANSFERRED_J) transferred += 1;
    if (metadata?.[12] === SOURCE_VALIDATED_J_MOTIF) motifValidated += 1;
    if (!annotationPresent(segment, metadata) && warnings.length < 8) {
      warnings.push(`${name}: ${segment === "V" ? "no validated IMGT FWR/CDR delineation" : segment === "J" ? "no unique frame-consistent F/W-G J motif" : "no segment metadata"}.`);
    }
    loci.add(locus);
    records.push({ header: input.header, name, rawSequence: input.rawSequence, sequence: normalized.sequence, locus, metadata });
  }

  const annotated = records.filter((record) => annotationPresent(segment, record.metadata)).length;
  const fasta = records.map((record) => `>${record.name}${metadataHeader(record.metadata)}\n${record.sequence}\n`).join("");
  return {
    fasta,
    count: records.length,
    annotated,
    unannotated: records.length - annotated,
    exactImgt,
    transferred,
    motifValidated,
    ambiguousBases,
    loci: [...loci].sort(),
    warnings,
  };
}

export function annotationCoverage(fasta: string, segment: GermlineSegment): { annotated: number; total: number } {
  const records = parseFasta(fasta);
  let annotated = 0;
  for (const record of records) {
    const normalized = normalizeSequence(record.rawSequence);
    const metadata = metadataFromHeader(record.header);
    if (metadata && validateMetadata(metadata, normalized.sequence.length, segment) && annotationPresent(segment, metadata)) annotated += 1;
  }
  return { annotated, total: records.length };
}

export function alleleMetadataHeader(allele: MetadataAllele): string {
  return `>${allele[0]}${metadataHeader(allele[2])}\n${allele[1]}\n`;
}

export function metadataSource(metadata?: CompactMetadata): string {
  return metadata ? METADATA_SOURCE_LABELS[metadata[12]] ?? "unclassified annotation" : "none";
}
