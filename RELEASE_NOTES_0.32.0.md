# Swig 0.32.0 release notes

## CLI custom-reference metadata parity

- Local V/D/J/C FASTA supplied through `references.files` is now prepared by default before SwiftIG builds its per-worker indexes.
- The CLI and browser now call the same multi-tier germline preprocessing driver and the same species-first taxonomic template construction. V FWR/CDR boundaries and J frame/CDR3 anchors are projected and validated by the existing shared kernel; D/C files receive the same parsing, normalization, locus/segment, and embedded-metadata checks.
- Preparation uses the fixed reference pack already embedded in every standalone CLI executable and requires no network access.
- `references.prepareMetadata: false` is the explicit pass-through opt-out. Existing browser-exported `references.inline` FASTA remains exact and is not needlessly reprocessed.
- Per-segment preparation counts and warnings are reported on stderr and in the machine-readable run summary. The resolved config records the effective flag.

Regression coverage verifies default V/J `SWIGMETA` transfer, the opt-out's literal raw-reference behavior, config normalization, browser/CLI use of the shared taxonomic tiers, and the complete standalone CLI pipeline.
