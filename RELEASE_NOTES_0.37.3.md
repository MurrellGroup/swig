# Swig 0.37.3

## Stable browser-assignment progress and live scheduler telemetry

- The **Calling rearrangements** card now reserves stable space for its stage text, percentage, telemetry, progress phases, and controls. Changing record counts or phase labels no longer shifts the card's surrounding layout.
- A committed queries-per-second reporter updates at most twice per second and uses a time-aware two-second exponential moving average. It measures AIRR rows only after their batch clears the browser result-store acknowledgement barrier, so it reflects the throughput users actually experience.
- The card reports occupied assignment workers against the effective pool size and labels initialization, streaming, final-batch drain, and local-index finalization separately.
- Per-dataset WASM loading and concurrent germline-index construction are now visible as initialization. This explains the normal brief all-core burst at a dataset boundary without conflating it with sustained query assignment.

Assignment algorithms, AIRR fields and ordering, browser worker limits, batch sizes, result-store behavior, and CLI behavior are unchanged.
