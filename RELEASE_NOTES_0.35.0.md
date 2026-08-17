# Swig 0.35.0 release notes

## Three web workflows

- Replaced the scrolling overview with three explicit entry routes: focused VDJ annotation/visualization, the unchanged end-to-end repertoire workflow, and direct single-lineage analysis.
- Focused VDJ mode accepts one uploaded or pasted dataset and exposes references, advanced mapping, and Double-D controls without study metadata, filtering, clustering, or post-analysis stages.
- Direct single-lineage mode annotates one input, treats it as one predefined lineage, bypasses clustering and selection, and opens the MSA/tree/UCA workbench directly.
- Concatenated gzip member detection remains explicit. A focused input can be merged, or separate members can move into the end-to-end workflow with editable sample metadata.

## Tree and UCA separation

- The lineage workbench can validate and display an uploaded Newick tree for the exact current alignment, including an optional synthetic germline guide for ordinary visualization.
- Phylogenetic UCA inference always excludes the synthetic germline guide. An eligible observed-only uploaded tree can be selected explicitly, while a fresh observed-only FastTree inference remains available independently.
- A display tree containing the guide is marked display-only and is never passed to UCA inference.

## Defaults and persistence

- RIAT-MP is now the default assignment strategy for all web workflows. CLI defaults and CLI pipeline behavior are unchanged.
- Saved sessions retain the selected web workflow and restore into the corresponding focused or end-to-end result surface.

## Verification

- Added deterministic uploaded-tree validation tests for exact/original/numeric identifiers, branch lengths, duplicate or missing tips, and germline-guide eligibility.
- Added UI regression checks for the three-route landing page, direct-lineage bypass, RIAT-MP web default, and the observed-only UCA tree boundary.
