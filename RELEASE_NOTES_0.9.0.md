# Swig 0.9.0

## Lineage phylogeny workbench

- Displays the exact named nucleotide alignment consumed by FastTree alongside the corresponding tree.
- Locks alignment rows to leaf vertical positions and draws dotted connectors from shorter/basal tips.
- Sizes tip bubbles by retained `duplicate_count` / deduplicated multiplicity.
- Selects full alignment, variable domain, all CDRs, individual FWR/CDR regions, custom coordinate ranges, or motif-matched columns.
- Uses standard nucleotide or amino-acid colors, or assigns a distinct color to each user-supplied motif.
- Controls branch-length/topology mode, tree width, row height, residue width, and full-screen display.
- Exports the tree, aligned leaf sequences, region track, mutation annotations, labels, and bubbles in a single SVG.

## Ancestral reconstruction

- Maps nucleotide substitutions and indels to branches with equal-cost Sankoff/Fitch-equivalent parsimony.
- Constrains the UCA at known reconstructed-germline bases.
- Treats germline `N` as an unknown A/C/G/T state, infers it from descendants, and does not emit an artificial `N→base` event.
- Retains gaps as a fifth state so insertion/deletion changes remain visible.
- Toggles branch labels and limits them to the currently displayed alignment columns.

## Kabat coordinates

- Adds lazy-loaded in-browser Kabat numbering for IGH, IGK, and IGL with Immunum 1.2.0 WASM.
- Projects numbering through codon-alignment gaps and takes a per-column vote across up to 24 aligned lineage members.
- Accepts Kabat insertion labels and ranges in the custom column selector.
- Leaves TCR views on alignment coordinates because Kabat does not define TCR numbering.

## Verification

- Adds tests for germline-constrained UCA reconstruction, unknown germline bases, indel states, motif matching, range parsing, and Kabat projection through codon gaps.
- Retains the existing 50,000-record streaming, storage, germline, post-analysis, and WebAssembly tests.
