# Swig 0.24.7

## Legible HMM-source tracks

- Made the HMM track-label rail remain visible during horizontal scrolling.
- Added distinct V, NT, D, and J row backgrounds and a compact legend.
- Collapsed route/register duplicates of the same allele into one visualization row. Conflicting registrations are rendered as nucleotide mixtures; the underlying HMM and serialized route-specific tracks are unchanged.
- Grouped non-template and unresolved boundary mass into spatial NT rows, then ordered tracks by posterior-weighted V → NT → D-block → NT → J position.
- Added detailed hover text for full track identity, route/register provenance, occupancy, and numeric nucleotide composition.

## UCA controls

- Expanded Additional-D probability to its full 0–1 input domain without changing HMM routing mathematics.
- Added a one-click reset for every advanced UCA option.

## Verification

- Added deterministic aggregation, mixture, biological-ordering, control-bound, floating-label, and color regression coverage.
- Verified the complete automated test suite, TypeScript check, and production build.
