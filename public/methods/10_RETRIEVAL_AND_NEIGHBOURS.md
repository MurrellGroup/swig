# Targeted retrieval, expansion, and lineage-neighbour search

## Query parsing and constraints

The query box accepts one or more raw sequences or FASTA records. Targets are CDR3 nucleotide, CDR3 amino acid, or V–J-trimmed nucleotide sequence. Manual locus/V/J constraints are checked with the selected call resolution/ambiguity policy.

For complete nucleotide rearrangements, `infer constraints` runs the same SwiftIG references, V strategy, calling profile, identity floor, and strand setting as the main analysis independently for each query. The inferred V/J calls constrain that seed only; non-empty manual overrides replace the corresponding inferred field. This is a fresh query annotation, not an inference from the repertoire allele-pooling model.

## Similarity modes

- **exact:** normalized strings must be identical;
- **substring:** either string may contain the other; score is shorter/longer length;
- **Hamming:** equal lengths only; exact distance must satisfy the identity-derived bound;
- **edit:** banded Levenshtein distance with bound `floor((1−identity) × max length)`;
- **VDJ sketch:** eight independently seeded MinHash values over 7-mers; score is the fraction of equal minima.

Sequence mode returns the best-scored ordinals. Lineage mode scores every assigned lineage by its best member/query pair and additionally reports matched-member and matched-query counts. The result cap is explicit.

The VDJ sketch is a **custom compact retrieval application of MinHash**, not an alignment likelihood. MinHash's resemblance estimator originates with Broder: [On the resemblance and containment of documents](https://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/broder97resemblance.pdf).

## Fixed-point CDR3 expansion

Expansion starts from current sequence hits. For the selected CDR3 identity it builds the same complete equal-length `d+1` block index used by lineage assignment, then performs a breadth-first traversal of exact Hamming-compatible V/J/study-scope edges. It stops only at a fixed point or the displayed result cap. This is read-only; it does not change lineage assignments.

## CDR3 lineage neighbours

For one or more selected source lineages, Swig indexes active representatives by study scope, locus, compatible V/J token, CDR3 length, and `d+1` blocks. Every candidate is exact-Hamming verified. A candidate lineage is ranked by its best source-member/candidate-member identity. This search is designed to review groups separated by the original clustering threshold.

## Inferred-germline neighbours

To avoid loading every lineage, Swig makes one eight-word provisional-ancestor sketch per lineage in one AIRR scan, using its least-mutated/highest-coverage representative only for screening. The top sketch candidates are loaded and re-inferred with the currently selected exact lineage-root method. Final similarity is banded Levenshtein identity between trimmed inferred roots, with `N` treated as unknown/wildcard. The sketch never supplies the final score.

Combining lineages is always explicit. The merge register adds a derived merged ID and preserves every original lineage ID.

## Literature relationship and limitations

These retrieval and neighbour algorithms are **custom Swig utilities**. They are not partis clonal-family inference, BLAST, MMseqs2, or a phylogenetic placement method.

- Hamming and edit identities are user-selected heuristics.
- MinHash estimates k-mer-set similarity and can miss positional/indel structure.
- A provisional root sketch can screen out a biologically relevant candidate when the representative is poor; broaden candidate count/threshold for sensitivity checks.
- Explicit merging is a human adjudication layer, not statistical evidence that the groups are one clone.
