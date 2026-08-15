export interface AlignmentGapSemantics {
  firstObserved: number;
  lastObserved: number;
  internalGaps: number;
  terminalMissingGaps: number;
}

export function alignmentGapSemantics(rawSequence: string): AlignmentGapSemantics {
  const sequence = rawSequence.toUpperCase().replaceAll(".", "-").replaceAll("U", "T");
  let firstObserved = -1;
  let lastObserved = -1;
  for (let site = 0; site < sequence.length; site += 1) if (sequence[site] !== "-" && sequence[site] !== "?") {
    if (firstObserved < 0) firstObserved = site;
    lastObserved = site;
  }
  let internalGaps = 0;
  let terminalMissingGaps = 0;
  for (let site = 0; site < sequence.length; site += 1) if (sequence[site] === "-") {
    if (isInternalAlignmentGap(site, { firstObserved, lastObserved })) internalGaps += 1;
    else terminalMissingGaps += 1;
  }
  return { firstObserved, lastObserved, internalGaps, terminalMissingGaps };
}

export function isInternalAlignmentGap(
  site: number,
  semantics: Pick<AlignmentGapSemantics, "firstObserved" | "lastObserved">,
): boolean {
  return semantics.firstObserved >= 0 && site > semantics.firstObserved && site < semantics.lastObserved;
}

export function isTerminalAlignmentGap(
  site: number,
  semantics: Pick<AlignmentGapSemantics, "firstObserved" | "lastObserved">,
): boolean {
  return !isInternalAlignmentGap(site, semantics);
}
