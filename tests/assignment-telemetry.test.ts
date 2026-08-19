import assert from "node:assert/strict";
import test from "node:test";

import {
  ASSIGNMENT_TELEMETRY_INTERVAL_MS,
  initializingAssignmentTelemetry,
  updateQueriesPerSecondEma,
} from "../src/assignment-telemetry.ts";

test("assignment throughput is a time-aware EMA of acknowledged AIRR records", () => {
  assert.equal(ASSIGNMENT_TELEMETRY_INTERVAL_MS, 500, "the visible rate must update at most twice per second");
  assert.equal(updateQueriesPerSecondEma(null, 0, 500), null, "do not show a false zero before the first commit");
  const initial = updateQueriesPerSecondEma(null, 1_000, 500);
  assert.equal(initial, 2_000);
  const decayed = updateQueriesPerSecondEma(initial, 0, 500);
  assert.ok(decayed !== null && decayed < initial && decayed > 0, "a stalled commit stream should decay smoothly toward zero");
  const recovered = updateQueriesPerSecondEma(decayed, 1_000, 500);
  assert.ok(recovered !== null && recovered > decayed && recovered < initial, "one fast interval should not make the displayed rate jump immediately");
});

test("initial assignment telemetry has stable, valid worker and counter fields", () => {
  assert.deepEqual(initializingAssignmentTelemetry(0), {
    phase: "initializing",
    activeWorkers: 0,
    totalWorkers: 1,
    queriesPerSecond: null,
    recordsParsed: 0,
    recordsCommitted: 0,
    recordsOutstanding: 0,
  });
  assert.equal(initializingAssignmentTelemetry(8.9).totalWorkers, 8);
});
