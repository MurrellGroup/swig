export const ASSIGNMENT_TELEMETRY_INTERVAL_MS = 500;
export const ASSIGNMENT_RATE_TIME_CONSTANT_MS = 2_000;

export type AssignmentPhase = "initializing" | "streaming" | "draining" | "finalizing";

export interface AssignmentTelemetry {
  phase: AssignmentPhase;
  activeWorkers: number;
  totalWorkers: number;
  queriesPerSecond: number | null;
  recordsParsed: number;
  recordsCommitted: number;
  recordsOutstanding: number;
}

/**
 * Time-aware EMA for the end-to-end rate at which annotated records clear the
 * browser result-store acknowledgement barrier. A zero delta deliberately
 * decays an established rate, while the display remains unset until the first
 * record is committed.
 */
export function updateQueriesPerSecondEma(
  previous: number | null,
  committedDelta: number,
  elapsedMilliseconds: number,
): number | null {
  if (!Number.isFinite(elapsedMilliseconds) || elapsedMilliseconds <= 0) return previous;
  const delta = Math.max(0, Number.isFinite(committedDelta) ? committedDelta : 0);
  const sample = delta * 1_000 / elapsedMilliseconds;
  if (previous === null) return delta > 0 ? sample : null;
  const alpha = 1 - Math.exp(-elapsedMilliseconds / ASSIGNMENT_RATE_TIME_CONSTANT_MS);
  return Math.max(0, previous + alpha * (sample - previous));
}

export function initializingAssignmentTelemetry(totalWorkers: number): AssignmentTelemetry {
  return {
    phase: "initializing",
    activeWorkers: 0,
    totalWorkers: Math.max(1, Math.floor(totalWorkers || 1)),
    queriesPerSecond: null,
    recordsParsed: 0,
    recordsCommitted: 0,
    recordsOutstanding: 0,
  };
}
