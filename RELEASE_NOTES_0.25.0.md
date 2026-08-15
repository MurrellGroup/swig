# Swig 0.25.0

## Fast exact-zero allele pooling

- Keeps the existing continuous Dirichlet allele-pooling model intact and selected by default.
- Adds a selectable fast hurdle active-set model with exact-zero excluded alleles, a configurable per-allele activity prior, a long-tailed positive-use slab, an inclusion-posterior threshold, and bounded one-dimensional quadrature.
- Reuses the existing compressed sparse evidence matrix and fits each donor/study-boundary, locus, and segment pool independently in the worker.
- Applies the same best-posterior and confidence-gated reassignment policies after either model.
- Recomputes hard before/after bars, vanished/appeared flags, CSV data, surviving-reference FASTA, AIRR overlay, and posterior sidecar from the selected model and current policy.
- Adds active/inactive state and approximate inclusion probability to model and sidecar exports.
- Exposes model selection in both the manual allele-pooling workspace and automatic pipeline.

## UCA SVG composition

- Full and visible-window HMM-track SVG exports now include the aligned UCA posterior frequency logo directly below the tracks.
- The visible export preserves the current horizontal window and vertical track scroll, keeps the floating track-label rail, and uses the same horizontal crop for the logo, numbering, and CDR bands.
- The standalone logo SVG remains available.
