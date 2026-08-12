# Swig 0.15.0

## Lineage alignment and germline inference

- The default in-Swig lineage alignment is now **Ref-anchored quick view**. Kalign nucleotide and codon-aware modes remain available.
- The lineage root template is no longer inherited from the first or closest member. Equal-weight member votes are placed on the shared AIRR V-reference origin; known V/D/J bases require 80% agreement, unresolved junction sites remain `N` for parsimony, and endpoints require 20% member support. Large lineages use a deterministic bounded reservoir so the estimate is not simply the first file/sample block.
- Manually corrected/imported alignments are stored per original-lineage set in saved sessions. Reproducible unedited alignments are not duplicated unless an inferred tree needs its exact input.

## Study-wide sample palette

- A central Results-page palette assigns one editable color to every biological sample.
- The same association is used in sample-composition figures, result labels, sample-stratified SHM summaries, and phylogeny tip circles.
- Tree tips can be colored by sample, original lineage, or uniformly. Palette changes propagate without rerunning analysis and persist in sessions.

## Lineage-neighbour explorer and explicit merges

- Search already-separated lineages by exact V/J-aware CDR3 identity below the assignment cutoff, inferred-lineage-germline identity, or either criterion.
- Germline retrieval uses a bounded one-sketch-per-lineage screen, followed by exact multi-member germline reconstruction and banded edit verification for the shortlist. The sketch is never the merge decision.
- Select multiple neighbours and view, align, and infer a tree on them as one temporary group; combined alignment input is sampled round-robin across original lineages.
- **Merge + view together** creates a derived merged-lineage assignment. Original `clone_id` values remain unchanged; exports add `swig_merged_lineage_id`, and saved sessions retain the merge register.
- SHM and lineage-aware missing-allele summaries use the derived merged assignments after an explicit merge while preserving the originals for audit/export.
