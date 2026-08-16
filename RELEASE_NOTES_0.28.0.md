# Swig 0.28.0 release notes

## Cancellation is now transactional

Every analysis-scale running state now has an explicit Cancel route. This covers input and concatenated-gzip inspection, reference preparation, initial V(D)J annotation, study-metadata re-indexing, the automatic pipeline, interactive repertoire-scale analyses, local-index searches, Kalign, FastTree, phylogenetic UCA inference, linked-session restoration, large table exports, and portable-session compression.

Cancellation no longer means merely hiding a spinner:

- dedicated analysis workers are terminated and recreated before the next run;
- streaming scans receive an `AbortSignal` and stop between bounded batches;
- post-analysis mutations run inside a checkpoint spanning the active mask, call overlay, collapse state, and lineage state;
- a cancelled mutation re-indexes the replacement worker and restores the complete preceding checkpoint before controls unlock;
- the automatic pipeline has one outer transaction, so no intermediate stage becomes authoritative;
- metadata correction uses a single abortable IndexedDB transaction and preserves the previous index on abort;
- alignment, tree, CHMMAIRRa, query, and UCA results are installed only after successful completion, leaving a previous result visible after cancellation; and
- streamed file writes abort their writable where supported, while Blob/download fallbacks are not handed to the browser until complete.

The progress control says **Restoring…** when rollback itself requires a bounded AIRR-index scan. A regression test interrupts metadata re-indexing after cursor updates have begun and verifies both indexed fields and facets remain unchanged.

## Complete methods audit

[`public/METHODS_INDEX.md`](public/METHODS_INDEX.md) maps every runnable interface block to a specification:

1. input, FASTQ QC, subsampling, and study structure;
2. reference composition and preparation;
3. V(D)J assignment and AIRR annotation;
4. Double-D evidence screening;
5. local storage, dashboard, filtering, and selection;
6. repertoire allele refinement;
7. collapse and all three denoising routes;
8. CHMMAIRRa;
9. lineage assignment;
10. SHM and possible-missing-allele diagnostics;
11. targeted retrieval and lineage neighbours;
12. lineage root construction, alignment, Alivibe correction, FastTree, and rooting;
13. phylogenetic UCA inference; and
14. execution order, persistence, export, cancellation, and rollback.

The documents give implemented equations, exact defaults, scope boundaries, output semantics, complexity, and limitations. Citations are labeled as direct implementation, browser port, compatible rule, or related precedent so a similarity is not presented as an implementation claim. Direct method links are now available from the post-analysis blocks and a complete Methods link is in the application footer.

No biological decision rule, prior, or phylogenetic likelihood was changed for this documentation/cancellation release.
