# Swig 0.24.6

## Explicit repertoire-allele reassignment

- Replaced the soft responsibility-frequency figure with weighted hard assignment counts: local-evidence argmax before pooling versus the policy-selected posterior MAP afterward.
- Added All, Changed, and Vanished filters plus complete CSV fields for count changes, appearance, and disappearance.
- Removed the truncated fitted-allele table from the page; the complete downloadable model summary remains available.
- Added two explicit application policies: posterior MAP for every modeled record, or posterior MAP only above a confidence gate (80% by default).
- Propagated the selected policy through live downstream overlays, pipeline mode, sessions, detail views, sidecar export, and refined AIRR export.
- Preserved the fitted sparse variational Dirichlet method unchanged; only the saved hard-projection metadata, reassignment choice, and visualization changed.

## Verification

- Added deterministic tests for hard MAP counts, confidence gating, vanished alleles, session round trips, and policy-aware posterior exports.
- Verified the complete automated test suite and production build.
