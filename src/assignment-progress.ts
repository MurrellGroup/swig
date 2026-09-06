/** Work fraction based on completed AIRR commits, never bytes buffered by gzip. */
export function assignmentCompletion(state: {
  acknowledged: number; parsed: number; inputDone: boolean;
  inputRecords: number; eligibleRecords: number; countHint: number | null;
  subsampleSize?: number; inputFraction?: number;
}): number {
  const { acknowledged, parsed, inputDone, inputRecords, eligibleRecords, countHint, subsampleSize } = state;
  const clamp = (value: number) => Math.max(0, Math.min(1, value));
  if (subsampleSize !== undefined) {
    const scan = countHint ? clamp(inputRecords / countHint) : (parsed || inputDone ? 1 : clamp(state.inputFraction ?? 0));
    const selected = inputDone ? parsed : Math.min(subsampleSize, countHint ?? subsampleSize);
    return 0.1 * scan + 0.9 * (selected ? clamp(acknowledged / selected) : (inputDone ? 1 : 0));
  }
  if (inputDone) return parsed ? clamp(acknowledged / parsed) : 1;
  // Unknown totals must not be inferred from compressed-byte read-ahead.
  if (!countHint) return Math.min(0.95, clamp(state.inputFraction ?? 0) * (inputRecords ? clamp((acknowledged + Math.max(0, inputRecords - eligibleRecords)) / inputRecords) : 0));
  const rejected = Math.max(0, inputRecords - eligibleRecords);
  return clamp((acknowledged + rejected) / countHint);
}

export function overallAssignmentProgress(completed: number, weight: number, total: number, fraction: number): number {
  return Math.min(0.99, 0.04 + 0.95 * (completed + Math.max(0, Math.min(1, fraction)) * weight) / Math.max(1, total));
}
