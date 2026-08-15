# Swig 0.24.0

## Repertoire-level allele pooling

- Added optional post-selection, pre-lineage germline-call refinement for overcomplete databases.
- Explicit co-optimal calls enter equally; retained scored alternatives and bounded reference-graph neighbours form sparse evidence rows.
- Added SHM-adaptive substitution-neighbour leakage using `baseline + sensitivity × mu/[3(1-mu)]`, with controls for the zero-SHM floor, sensitivity, cap, and SHM-estimate cap. Indel neighbours receive baseline leakage only.
- Added independent sparse Dirichlet mixtures by study boundary, locus, and segment. Donor/subject pooling, V+J modeling, unique-record weighting, and alpha 0.1 are defaults; D is experimental. The prior normalization includes every locus-matched database node, while prior-only nodes remain implicit rather than being allocated per record or donor.
- The worker stores no dense reads-by-database matrix and collapses repeated evidence patterns during updates.
- Applying posterior MAP calls is explicit and thresholded. Original AIRR calls remain recoverable; changing the overlay invalidates lineage-dependent results.
- Added model-summary, full sparse per-record posterior-sidecar, and provenance-preserving refined AIRR exports in TSV, CSV, and JSON Lines. Full responsibilities are reconstructed and streamed on demand rather than retained as another large in-memory matrix.
- Saved sessions retain posterior vectors, options, apply threshold, and active/reset state.
- Pipeline mode can fit and apply the same model between repertoire selection and lineage assignment.
- Added detailed algorithm documentation in `public/REPERTOIRE_ALLELE_REFINEMENT.md`.

## Validation

- Added tests for identical-reference equivalence classes, bounded neighbours, literal equal-weight ties, SHM-adaptive leakage, donor-isolated updates, repertoire resolution of ambiguous reads, and call apply/reset behavior.
