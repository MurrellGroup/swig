# Swig 0.12.0

## Correct cumulative post-analysis scope

CHMMAIRRa now reads the authoritative cumulative mask from the post-analysis worker immediately before starting. The AIRR store skips upstream-excluded ordinals before parsing or dispatching rows, so chimera inference, lineage assignment, sequence querying, and expansion all consume the same explicit working set. The CHMMAIRRa result reports its working-set input and upstream-excluded count before a posterior threshold is applied downstream.

Exact collapse and denoising now expose an ineligible-record policy. Records without a usable selected exact key, or without both a trimmed sequence and V/J calls for a denoising partition, are discarded from the downstream representative set by default. “Retain unchanged” is available beside the method controls. The preview, warning, applied-stage summary, multiplicity vector, and downstream mask use the selected policy.

## Opt-in double-D / VDDJ evidence

The ordinary caller is unchanged by default. “Off” calls the legacy WebAssembly annotation ABI. Two optional modes screen either every eligible IGH/TRB/TRD junction or only baseline V-end/J-start intervals at or above a user-set length.

The screen uses an indexed exact D seed (11 nt by default), anchored independent extensions, ordered non-overlap, a configurable two-D score gain, and an IgScout-style single-D pseudo-tandem Δ-distance rejection. Expensive pseudo-tandem calculations are evaluated in evidence order and cached by outer span. Supported calls are sparse sidecar records; they never rewrite the normal AIRR TSV. The result index, filter, detail map, D1/D2 nucleotide or amino-acid alignments, and a separate AIRR-2-compatible evidence TSV all understand the second D segment. A regression test asserts that the baseline AIRR output is byte-for-byte identical when the sidecar screen is enabled.

## Background-tab execution

Long analysis and post-analysis operations acquire a semantic same-origin Web Lock. In Chrome this excludes the browsing-context group from Energy Saver freezing while the lock is held; the progress UI reports the state. Cooperative AIRR scan yielding now uses `scheduler.yield()` or `MessageChannel` rather than a timer that hidden-tab throttling can stretch to one second per batch. A navigation warning protects an active main run. Operating-system or browser emergency page discard remains outside a site's control; Chrome's site performance exception is still appropriate for unattended critical runs.

## Verification

- optimized SIMD/LTO WebAssembly rebuilt from the updated C++ core;
- 50,000-record IndexedDB/storage/filter/export regression passes;
- cumulative-mask, exact-collapse, denoising, lineage, query, and expansion tests pass;
- double-D C++ call, sparse store, viewer fields, separate export, long-span gate, and unchanged-baseline invariants pass;
- complete TypeScript check and production build pass.
