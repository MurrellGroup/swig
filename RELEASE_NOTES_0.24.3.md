# Swig 0.24.3

## Error-model parameter inspection

- Added a reference-neighbour evidence inspector inside the allele-pooling advanced controls.
- V, D, or J alleles can be located with a searchable database selector; sequence-identical labels remain one explicit unresolved node.
- An adjustable assumed-SHM control recomputes the exact sparse evidence kernel used by the fitting worker as model parameters are committed.
- Non-primary candidates are shown as horizontal bars with absolute percentage labels. The selected primary sequence remains in the alignment but has no bar, so small alternatives remain legible.
- Candidate sequences share one selected-reference coordinate. Columns that are gaps in every displayed sequence are removed.
- The alignment switches between standard nucleotide colors and a difference highlighter for matches, substitutions, and indel columns.
- The configured reference edit radius now consistently supports the full UI range through five edits.

## Before/after assignment frequencies

- Added a fitted-pool selector for donor/study-boundary, locus, and segment models.
- Each modeled allele has paired bars for normalized local evidence before repertoire pooling and normalized variational assignment responsibility after pooling.
- Rows are sorted by post-pooling frequency. Dirichlet prior-only mass is excluded from both plotted series so the figure isolates reassignment of observed records.
- The displayed allele limit is adjustable; CSV always includes the complete fitted pool.
- The paired-bar figure exports as a self-contained SVG with the selected pool and both percentage series.

## Verification

- Added regression tests for SHM-responsive alternative probability, exact row normalization, omission of the primary from alternative candidates, rectangular/all-gap-free reference alignments, and post-frequency sorting.
- Added UI wiring checks for both diagnostics and their SVG/data exports.
