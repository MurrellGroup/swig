# Swig 0.29.0 release notes

## Post-analysis flow and lineage selection

- Repertoire selection (step 04) is skipped by default.
- Every post-analysis card has **Skip step** / **Include step** controls, retained in sessions and respected by guided navigation.
- The lineage explorer can require membership from one selected sample without excluding lineages that also contain other samples.
- Each lineage row reports exact `duplicate_count`-weighted reads for every represented sample ID.
- Once SHM is available, the lineage selector reports weighted mean SHM and upper 95% quantile.
- “Diagnostics” is now **Post-lineage analysis**. SHM and possible-missing-V evidence are separate full-width sections stacked vertically.
- One possible missing allele supported by more than 50 independent lineages, or more than five lower-support candidates, triggers a prominent recommendation to use a more complete database or personalized germline discovery with IgDiscover.

## Standalone `swig-cli`

`swig-cli` runs V(D)J assignment and the non-phylogenetic downstream pipeline from the same WASM and TypeScript analysis kernels as Swig Web. It accepts FASTA, FASTQ, AIRR TSV, and gzip-compressed inputs; supports explicit sample-to-donor metadata; and can run allele reassignment, collapse/denoising, CHMMAIRRa, selection, lineage assignment, SHM, and possible-missing-V evidence.

End users do not install Node, npm, Bun, Python, a container engine, or a WASM runtime. Release automation embeds everything into one native executable and smoke-tests Linux x86-64/ARM64, macOS Intel/Apple-silicon, and Windows x86-64 artifacts. Node remains only a maintainer dependency for building the web repository and its inexpensive development CLI bundle.

The CLI consumes a compact JSON config with defaults. Swig Web can export the exact references, methods, settings, dataset labels, and donor map from a browser analysis as a reusable CLI config. Phylogenetics is intentionally excluded from CLI pipeline execution because the lineage MSA requires inspection and possible curation.

## Lazy lineage-study loading

CLI lineage runs emit a lineage-sorted AIRR TSV plus a compressed manifest. The manifest contains exact references/options, sample-level lineage counts, optional SHM mean/q95, and SHA-256-protected byte ranges. **Load lineage study** opens the summaries without indexing the whole AIRR file. Selecting a lineage reads and verifies only its byte range, then opens the normal alignment, tree, and UCA workbench on that lineage.

## Reproducibility and tests

The release adds regression coverage for sample/donor config export, per-sample lineage abundance, incomplete-reference escalation, manifest validation, lazy byte-range integrity, and an end-to-end CLI run through lineage-study export. The standalone release workflow repeats the end-to-end test against the actual native binary on every supported platform.
