import { useMemo } from "react";

type AirrRow = Record<string, string>;
type AlignmentMode = "nt" | "aa";

interface AlignmentResult {
  query: string;
  reference: string;
  score: number;
  queryFrame: number;
  referenceFrame: number;
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
  const clean = sequence.toUpperCase().replace(/[^ACGTN]/g, "");
  let protein = "";
  for (let index = frame; index + 2 < clean.length; index += 3) {
    protein += CODONS[clean.slice(index, index + 3)] ?? "X";
  }
  return protein;
}

function proteinAlign(query: string, reference: string): Omit<AlignmentResult, "queryFrame" | "referenceFrame"> {
  const columns = reference.length + 1;
  const rows = query.length + 1;
  const scores = new Int16Array(rows * columns);
  const trace = new Uint8Array(rows * columns);
  for (let row = 1; row < rows; row += 1) {
    scores[row * columns] = -2 * row;
    trace[row * columns] = 1;
  }
  for (let column = 1; column < columns; column += 1) {
    scores[column] = -2 * column;
    trace[column] = 2;
  }
  for (let row = 1; row < rows; row += 1) {
    for (let column = 1; column < columns; column += 1) {
      const diagonal = scores[(row - 1) * columns + column - 1] +
        (query[row - 1] === reference[column - 1] ? 3 : -1);
      const up = scores[(row - 1) * columns + column] - 2;
      const left = scores[row * columns + column - 1] - 2;
      const best = Math.max(diagonal, up, left);
      scores[row * columns + column] = best;
      trace[row * columns + column] = best === diagonal ? 0 : best === up ? 1 : 2;
    }
  }
  let row = query.length;
  let column = reference.length;
  let alignedQuery = "";
  let alignedReference = "";
  while (row || column) {
    const direction = trace[row * columns + column];
    if (row && column && direction === 0) {
      alignedQuery = query[--row] + alignedQuery;
      alignedReference = reference[--column] + alignedReference;
    } else if (row && (direction === 1 || !column)) {
      alignedQuery = query[--row] + alignedQuery;
      alignedReference = `-${alignedReference}`;
    } else {
      alignedQuery = `-${alignedQuery}`;
      alignedReference = reference[--column] + alignedReference;
    }
  }
  return { query: alignedQuery, reference: alignedReference, score: scores[scores.length - 1] };
}

function bestProteinAlignment(queryNucleotides: string, referenceNucleotides: string): AlignmentResult | null {
  const queryClean = queryNucleotides.replace(/[-.]/g, "");
  const referenceClean = referenceNucleotides.replace(/[-.]/g, "");
  if (queryClean.length < 3 || referenceClean.length < 3) return null;
  let best: AlignmentResult | null = null;
  for (let queryFrame = 0; queryFrame < 3; queryFrame += 1) {
    for (let referenceFrame = 0; referenceFrame < 3; referenceFrame += 1) {
      const aligned = proteinAlign(translate(queryClean, queryFrame), translate(referenceClean, referenceFrame));
      const candidate = { ...aligned, queryFrame, referenceFrame };
      if (!best || candidate.score > best.score) best = candidate;
    }
  }
  return best;
}

function matchLine(query: string, reference: string): string {
  return [...query].map((value, index) => {
    const other = reference[index];
    if (value === "-" || other === "-") return " ";
    return value === other ? "│" : "·";
  }).join("");
}

function blocks(query: string, reference: string, width = 72) {
  const output: Array<{ query: string; reference: string; match: string }> = [];
  for (let index = 0; index < Math.max(query.length, reference.length); index += width) {
    const queryBlock = query.slice(index, index + width);
    const referenceBlock = reference.slice(index, index + width);
    output.push({ query: queryBlock, reference: referenceBlock, match: matchLine(queryBlock, referenceBlock) });
  }
  return output;
}

function percent(value: string): string {
  const parsed = Number(value);
  return Number.isFinite(parsed) && value !== "" ? `${(parsed * 100).toFixed(1)}%` : "—";
}

const SEGMENTS = [
  { key: "v", label: "V", color: "var(--segment-v)" },
  { key: "d", label: "D", color: "var(--segment-d)" },
  { key: "j", label: "J", color: "var(--segment-j)" },
  { key: "c", label: "C", color: "var(--segment-c)" },
] as const;

function SegmentAlignment({ row, segment, mode }: {
  row: AirrRow;
  segment: typeof SEGMENTS[number];
  mode: AlignmentMode;
}) {
  const call = row[`${segment.key}_call`];
  const nucleotideQuery = row[`${segment.key}_sequence_alignment`] ?? "";
  const nucleotideReference = row[`${segment.key}_germline_alignment`] ?? "";
  const protein = useMemo(
    () => mode === "aa" ? bestProteinAlignment(nucleotideQuery, nucleotideReference) : null,
    [mode, nucleotideQuery, nucleotideReference],
  );
  if (!call) return null;
  const query = mode === "aa" ? protein?.query ?? "" : nucleotideQuery;
  const reference = mode === "aa" ? protein?.reference ?? "" : nucleotideReference;
  const renderedBlocks = blocks(query, reference);

  return (
    <article className="alignment-card" style={{ "--segment-color": segment.color } as React.CSSProperties}>
      <header>
        <span className="segment-chip">{segment.label}</span>
        <div><strong>{call}</strong><small>{percent(row[`${segment.key}_identity`])} identity · {row[`${segment.key}_cigar`] || "no CIGAR"}</small></div>
        <dl>
          <div><dt>Query</dt><dd>{row[`${segment.key}_sequence_start`] || "—"}–{row[`${segment.key}_sequence_end`] || "—"}</dd></div>
          <div><dt>Germline</dt><dd>{row[`${segment.key}_germline_start`] || "—"}–{row[`${segment.key}_germline_end`] || "—"}</dd></div>
          <div><dt>Score</dt><dd>{row[`${segment.key}_score`] || "—"}</dd></div>
        </dl>
      </header>
      {mode === "aa" && protein && (
        <p className="translation-note">On-demand translation · query frame +{protein.queryFrame + 1}, germline frame +{protein.referenceFrame + 1}</p>
      )}
      {renderedBlocks.length ? renderedBlocks.map((block, index) => (
        <div className="alignment-block" key={index}>
          <div><span>query</span><code>{block.query}</code></div>
          <div className="match-row"><span /><code>{block.match}</code></div>
          <div><span>{segment.label.toLowerCase()} ref</span><code>{block.reference}</code></div>
        </div>
      )) : <p className="empty-alignment">No aligned sequence was reported for this segment.</p>}
    </article>
  );
}

export function AlignmentViewer({ row, mode, onMode }: {
  row: AirrRow;
  mode: AlignmentMode;
  onMode: (mode: AlignmentMode) => void;
}) {
  return (
    <section className="alignment-viewer">
      <div className="alignment-toolbar">
        <div>
          <span className="section-kicker">Per-segment evidence</span>
          <h3>Query ↔ germline alignments</h3>
        </div>
        <div className="mode-toggle" role="group" aria-label="Alignment alphabet">
          <button className={mode === "nt" ? "active" : ""} type="button" onClick={() => onMode("nt")}>Nucleotide</button>
          <button className={mode === "aa" ? "active" : ""} type="button" onClick={() => onMode("aa")}>Amino acid</button>
        </div>
      </div>
      {mode === "aa" && (
        <p className="alignment-explainer">Protein views are calculated only for the open record by translating all frame pairs and retaining the highest-scoring alignment.</p>
      )}
      <div className="alignment-stack">
        {SEGMENTS.map((segment) => <SegmentAlignment key={segment.key} row={row} segment={segment} mode={mode} />)}
      </div>
    </section>
  );
}

