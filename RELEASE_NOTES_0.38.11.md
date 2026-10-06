# Swig 0.38.11

Assignment reference controls now download the active database as one FASTA, by V/D/J/C segment, or by locus/segment cell. Exports retain prepared metadata and allele exclusions.

Lineage guides now retain the V 5′ and J 3′ reference flanks when their selected template read is truncated. Completion requires exact called records and coordinate-consistent AIRR segment alignments. It extends the germline guide while keeping missing tip coverage as gaps, and preserves V 3′ / J 5′ recombination trimming. Tied-reference differences remain uncertain rather than choosing an arbitrary alias. PhyloHMM projections retain ambiguous guide coordinates and candidate-specific reference constraints.

Automatic template selection favors complete outer V/J alignments before mutation identity, within the existing safely projected VDDJ preference. A manual trimming-template selector overrides the automatic choice, saves it per lineage group/productivity view, and retains it under MSA row limits. An intact template query can reproject its guide onto an existing curated MSA without changing observed nucleotide content. Unsafe query/gap mappings require a fresh MSA. Older guides missing outer flanks are repaired when that projection is safe; their interior is retained and old tree/UCA results are invalidated. PhyloHMM boundary anchors use the chosen member when retained in the MSA.

Validation includes an AER-R/R-optimized WASM-to-AIRR-to-guide regression with both read ends truncated, complete-versus-partial template ranking, manual trimming and row-limit retention, tied allele flanks, and an end-to-end phyloHMM case where a truncated template correctly constrains the germline prefix despite unanimous mutated full-length tips. The underlying AIRR/MSA from the supplied screenshots was not attached, so that exact biological lineage has not been rerun.

All 286 tests and the production build pass. The final template-coordinate changes also pass the focused 40-test lineage/UCA suite and the TypeScript check.

Assignment algorithms, AER-R/R-optimized defaults, the 0.38.10 speed optimizations, and both repertoire germline-discovery methods are retained.
