# Swig 0.16.0

This release refines multi-dataset lineage analysis, guided post-analysis, alignment editing, and tree annotation while retaining the static GitHub Pages deployment model.

## Lineage roots and related-lineage search

- The default lineage germline/UCA template is the single loaded member with the best equal-weight V/J identity across its matched germline regions. Combined identity, informative coverage, and AIRR order provide deterministic tie-breaks.
- Equal-weight AIRR-anchored member consensus remains available as an explicit alternative. The selected method is saved in sessions and used consistently for alignment, phylogeny, and exact lineage-neighbour rescoring.
- Targeted repertoire search has a lineage-return mode. It ranks each lineage by its best match to any member sequence and opens a selected hit directly in the lineage workbench.

## Correctable study metadata

- Sample, donor/subject, cohort, and timepoint labels can be edited after V(D)J assignment.
- Applying an edit updates the compact IndexedDB metadata index, filters, facets, display colors, and AIRR exports without rerunning annotation.
- Existing downstream state is invalidated deliberately; rerunning collapse, chimera filtering, and lineage assignment uses the corrected biological boundaries.

## Guided post-analysis interface

- Every post-analysis step card is collapsible. Applying a step closes it and opens the next step.
- The next scientifically valid default action is highlighted so a first-time user can follow the standard workflow one button at a time.
- The lineage table sorts from any column header and starts with abundance from largest to smallest.
- Targeted sequence/lineage querying is placed at the end of the workflow.

## Alignment and phylogeny details

- Alivibe round trips may deliberately remove bad biological rows, alignment columns, or nucleotide characters. Added/renamed rows, substituted bases, reordered content, and removal of the N-masked germline remain rejected.
- Manually edited alignments, including deletions, remain session-persistent; reproducible unedited alignments are not duplicated unless a saved tree depends on their exact bytes.
- Mutation labels expose independent controls for the number shown before `+N` and their font size. Both are preserved in the composite SVG.

## Verification

- Added regression coverage for closest-member versus consensus germline construction, equal weighting of V and J identity, metadata reindex/export overlays, deletion-preserving Alivibe imports, dataset-wide best-member lineage query aggregation, and the existing 50,000-record IndexedDB path.
- Type checking, the complete test suite, and a production Vite build are required before packaging.
