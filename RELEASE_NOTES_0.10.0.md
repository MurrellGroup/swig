# Swig 0.10.0

## Tree and alignment display

- Corrects phylogram normalization for ordinary fractional substitution/site branch lengths; the maximum root-to-tip distance now spans the complete tree viewport.
- Lets tree width and vertical tip spacing reach zero, with proportional padding rather than a hidden minimum geometry.
- Places parsimony mutation labels at branch midpoints with a light text halo instead of node-adjacent boxes.
- Adds a substitution/site scale bar, cleaner branch and node styling, larger labels, and area-proportional multiplicity bubbles.
- Adds Newick-order, large-clade-first, and small-clade-first tip ordering; the displayed ordering is retained when downloading Newick.
- Removes vertical padding between rows in the standalone lineage alignment preview.

## Germline rooting and Kabat numbering

- Roots exactly at the N-masked germline: the germline root edge is zero and the ingroup receives the full original connecting-edge length. Pairwise patristic distances are unchanged.
- Continues Kabat numbering when a terminal partial codon or stop codon is detected and shows a local reliability warning instead of rejecting the alignment.

## Collapse and denoising

- Adds three explicit post-assignment modes: exact deduplication, FAD-compatible denoising, and an experimental conservative exact-neighbor error model.
- FAD mode implements the published corrected 6-mer distance, method 1/2 template rule, Poisson tail, sequence-length correction, and final nearest accepted centroid.
- Replaces dense all-pairs FAD candidate scans with a complete threshold index and exact VP-tree nearest-neighbor search after locus/V/J partitioning.
- Adds a conservative alternative using a `d+1` block index, exact Hamming verification, and a sequence-specific Poisson error model; isolated singletons remain independent.
- Streams trimmed VDJ sequences into two-bit chunked storage and releases temporary neighbor profiles one V/J partition at a time.
- Preserves and sums existing `duplicate_count` values through every collapse mode, lineage abundance, AIRR export, and phylogeny bubble.
- Adds a reproducible 50,000-record denoising benchmark.
