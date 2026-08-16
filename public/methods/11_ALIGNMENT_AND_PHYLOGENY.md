# Lineage rooting, alignment, and ordinary phylogeny

**Evidence label:** custom lineage-root reconstruction; direct use of Kalign 3.3.1 and FastTree 2.1.11 through browser WebAssembly.

This document specifies the ordinary lineage workbench. It is distinct from Swig's phylogenetic-UCA model: the ordinary tree contains a reconstructed, N-masked germline guide as a rooting anchor, whereas phylogenetic UCA inference removes that guide, estimates a tree from observed sequences only, and places an unobserved ancestor against the fixed tree.

## Rows entering the workbench

The workbench consumes the current cumulative working set and one selected lineage, or an explicit set of lineages after a user merge. `duplicate_count` affects displayed abundance but does not create repeated alignment rows. If the selection exceeds the user-set maximum, Swig groups rows by their original lineage ID and takes one row from each group in round-robin order until the cap is reached. Selection is deterministic for a fixed AIRR order. This preserves representation from every merged component but is not a random or posterior sample.

Every row is projected to a common reference origin using the one-based AIRR `v_germline_start`. The caller's aligned sequence, germline alignment, gaps, and insertion structure are retained; shorter rows receive terminal gap padding. A supported VDDJ row is reconstructed from its D1 and D2 AIRR projections when both can be mapped safely. Otherwise the original single-D germline composite is retained and the incomplete Double-D status is reported.

## Root/guide alternatives

### Closest member (default)

This is a **custom heuristic**, not an implementation of a published ancestral-reconstruction method.

For each loaded representative, Swig calculates nucleotide identity independently on informative V and J aligned query/germline columns. If both segments are informative, their identities receive equal weight, regardless of segment length. Candidates are ranked lexicographically by:

1. number of informative end segments (V and J);
2. equal-weight mean V/J identity;
3. identity over the combined informative V/J columns;
4. number of compared columns; and
5. original AIRR ordinal.

Older AIRR imports without segment-specific alignments fall back to identity on informative columns of the combined query/germline alignment. If any safely reconstructed VDDJ members exist, ranking is restricted to them so a baseline single-D composite cannot erase D2.

The selected member's aligned germline is the **tree guide**. Known reference bases remain fixed. At guide positions marked `N`, the selected member's observed nucleotide may fill a separate comparison-UCA sequence, but the guide used to root the ordinary tree deliberately remains `N`. Thus observed junction bases do not masquerade as known germline in the tree.

### Equal-weight member consensus

Each loaded unique representative contributes one vote per covered column; `duplicate_count` is not a vote multiplier. If safely reconstructed VDDJ members exist, only those members vote across the complete root so an unresolved single-D composite cannot erase D2.

- A reference nucleotide is accepted into the tree guide only when one untied base has at least 80% of informative germline votes.
- Columns with insufficient or conflicting reference support remain `N` when covered and `-` when uncovered.
- A separate comparison UCA may fill an unresolved column when one untied observed query base has at least 60% of query votes.
- The retained interval begins and ends where at least `ceil(0.20 n)` voting members have query or germline coverage.

These 80%, 60%, and 20% values are fixed implementation thresholds, not literature-derived universal constants.

## Alignment routes

### Reference-anchored quick view (default)

Swig keeps the caller's pairwise alignment columns, left-pads by `v_germline_start - 1`, and right-pads every row to a common width. It does not optimize a de-novo multiple alignment. The route is fast and preserves assignment coordinates, but unrelated insertions can occupy nominally corresponding columns; inspect the junction before phylogenetic use.

### Nucleotide Kalign

Swig removes input gaps, replaces identifiers with numeric IDs, runs **Kalign 3.3.1** in nucleotide mode, and restores the original names. This is direct use of Kalign's executable, not a reimplementation of its alignment algorithm.

### Codon-aware Kalign

For each row, Swig removes gaps, translates from the AIRR-derived frame offset, aligns the amino-acid strings with Kalign 3.3.1, and projects amino-acid gaps back as complete nucleotide codons. Leading nucleotides before the first complete codon and terminal incomplete codons are retained by the projection. This preserves codon-sized gaps, but it is **not** a native codon-substitution alignment and does not optimize a codon evolutionary likelihood.

The displayed amino-acid frame is an alignment property with allowed offsets 0, 1, or 2. FASTA does not encode it, so a manually imported alignment must be checked against the adjacent frame control.

## Manual correction through Alivibe

Swig opens a pinned Alivibe build and exchanges the exact ordered nucleotide records through a versioned browser bridge. The return is accepted only for the originating lineage selection and alignment fingerprint. A corrected alignment may:

- move or add gap characters;
- delete complete alignment columns;
- delete nucleotide characters while retaining their order; and
- delete biological rows.

It may not add, rename, or reorder nucleotide content, insert a new nucleotide, or substitute one nucleotide for another. The `__germline_N_masked__` row must remain because ordinary rooting depends on it. Swig verifies equal aligned width, unique names, the nucleotide alphabet, and that every returned ungapped row is a subsequence of its original. The manually selected amino-acid frame is stored separately.

## FastTree execution

Before every run Swig validates the exact current nucleotide MSA, writes it afresh to the browser-local WebAssembly filesystem, and substitutes numeric tip IDs. It runs the double-precision FastTree 2.1.11 executable with:

- `-nt -gtr` for the default nucleotide GTR approximation;
- `-nt` for the selectable Jukes-Cantor route; and
- optional `-fastest` when the user enables the fastest heuristic.

Original names are restored in the returned Newick. The executed named FASTA, exact numeric FASTA, command, row/column counts, and a deterministic alignment fingerprint are exposed. Swig does not recompute or silently repair the MSA during tree inference.

This is direct execution of FastTree, but the surrounding lineage selection, germline guide, rooting, and display are Swig behavior. It is not IgPhyML, HLP17/HLP19, or another antibody-specific codon phylogeny.

## Rooting and stable display

FastTree first estimates an unrooted tree containing the N-masked guide. Swig then reroots **exactly at the guide tip**: the guide-to-root edge becomes zero and the full original connecting-edge length is assigned to the ingroup side. No midpoint split is introduced, so pairwise path lengths are preserved.

Three outputs remain available:

- **Raw FastTree:** untouched returned topology and lengths;
- **Rooted, resolved:** rooted at the N-masked guide with all inferred resolutions retained; and
- **Rooted, floor-collapsed:** internal edges of length at most `1e-8` are contracted for a stable polytomy display. Leaves are never contracted.

Child order is canonicalized for deterministic rendering and can be ladderized for display without changing topology or branch lengths. A phylogram uses non-negative FastTree branch lengths exactly, including zero; cladogram mode uses unit depth.

## Exports and invalidation

The current MSA can be exported as aligned FASTA, Clustal, relaxed PHYLIP, Stockholm, or NEXUS. Tree exports include raw/rooted/floor-collapsed Newick, NEXUS, and the current SVG. A manually corrected MSA and its frame are retained in a saved session. Changing lineage membership, the germline method, the alignment, or its frame invalidates dependent tree and phylogenetic-UCA results rather than reusing stale coordinates.

## Limitations

- Neither guide method integrates over recombination or ancestral uncertainty.
- The quick view is coordinate preserving, not a de-novo MSA.
- A bounded lineage is a deterministic subset when it exceeds the selected maximum.
- FastTree is a practical approximate maximum-likelihood estimator and does not model SHM context, selection, or codon structure.
- The ordinary guide-rooted tree must not be interpreted as the posterior UCA analysis documented separately.

## References and exact relationship

- Lassmann T. [Kalign 3: multiple sequence alignment of large datasets](https://doi.org/10.1093/bioinformatics/btz795). *Bioinformatics* (2020). Swig executes Kalign 3.3.1; the reference-anchored quick view and codon projection are Swig wrappers, not methods from this paper.
- Price MN, Dehal PS, Arkin AP. [FastTree 2—approximately maximum-likelihood trees for large alignments](https://doi.org/10.1371/journal.pone.0009490). *PLOS ONE* (2010). Swig executes FastTree 2.1.11; lineage capping, the N-masked guide, exact-tip rooting, and floor-collapse display are Swig additions.
- Hoehn KB, Lunter G, Pybus OG. [A phylogenetic codon substitution model for antibody lineages](https://doi.org/10.1534/genetics.116.196303). *Genetics* (2017). HLP17 motivates caution about ordinary nucleotide models for antibody evolution; Swig's ordinary FastTree route does not implement HLP17.
