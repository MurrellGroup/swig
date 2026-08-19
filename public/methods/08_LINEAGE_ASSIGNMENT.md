# CDR3-based lineage assignment

## Eligibility and hard boundaries

The default fit uses productive records in the current working set, partitions by donor/subject, requires the same locus, resolves V/J at gene level, treats ambiguous calls as compatible when any assignment overlaps, and requires exact CDR3 nucleotide length. Collapse multiplicities do not create additional graph vertices; they contribute abundance to the retained representative.

Study scope is a hard boundary. No candidate pair is generated across different scope values. V/J compatibility policies are:

- **overlap:** any normalized V and any normalized J label overlap;
- **top:** only the first comma-separated call is used;
- **strict:** complete normalized call sets must match.

Gene resolution removes allele suffixes; allele resolution retains them.

## Exact accelerated single linkage

For CDR3 length \(L\) and identity threshold \(q\), the maximum Hamming distance is

\[
d=\left\lfloor(1-q)L\right\rfloor.
\]

Every eligible CDR3 is split into \(d+1\) disjoint blocks. Two equal-length strings within distance \(d\) must share at least one exact block, so indexing by partition/call token/block generates a complete candidate set unless the explicit pathological-bucket cap is reached. Every candidate is then verified with exact bounded Hamming distance. Accepted edges are combined by union–find; connected components are the reported lineages.

This is single linkage: A may link to B and B to C even when A and C are farther apart than the threshold. The default identity is 85%, a starting value rather than a universal clonal boundary.

## Summaries and exports

Each component receives a deterministic numeric ID. Summary abundance is the sum of multiplicities; unique members count active representatives. For every retained summary, Swig separately stores each sample ID's unique-representative count and `duplicate_count`-weighted read abundance. The explorer can therefore require that a lineage contain reads from one selected sample while still allowing that lineage to contain any number of other samples. This is a membership predicate, not an “only this sample” predicate.

The interactive table retains the 10,000 largest summaries with exact sample/donor/timepoint/compartment membership and displays the read count for every represented sample. Each row also shows a representative nucleotide/amino-acid CDR3. Before SHM calculation this is the component's deterministic representative; afterward Swig uses the member with the lowest V-nucleotide SHM rate (then mutation count and input ordinal as tie breaks). That minimum-member selection scans every active row and is not taken from the bounded plot reservoir. Weighted mean SHM and upper 95% quantile are joined by lineage ID and shown in the same selector. The assignment vector still covers every row, including lineages outside the bounded interactive summary table. Export writes `clone_id=swig_lineage_N` and an optional separate `swig_merged_lineage_id`.

Checkboxes can open several original lineages together in one alignment/tree/UCA workbench. This temporary combined view does not rewrite the original assignment vector or declare a biological merge. An explicit neighbour-review merge similarly preserves original lineage IDs; removing it recovers the original groups.

## Literature relationship

The biological partitioning and single-linkage/Hamming choice are **methodologically similar** to widely used distance-based B-cell clone definitions. Gupta et al. found single-linkage Hamming clustering effective under evaluated conditions: [Hierarchical clustering can identify B cell clones with high confidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC5340603/).

Swig is **not an implementation of Change-O, SHazaM threshold inference, SCOPer, or partis**. It uses a user-set threshold, its own ambiguity/study-scope semantics, a complete `d+1` candidate index, and union–find rather than calling those packages. Partis instead uses a likelihood model of rearrangement and clonal family structure: [Ralph and Matsen, 2016](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005086).

## Limitations

- Threshold choice, junction definition, sequencing error, and incomplete germline calls can split or join clones.
- Single linkage can chain through intermediates.
- Heavy-chain-only grouping cannot use paired light-chain evidence.
- A candidate-cap warning means completeness is not guaranteed for affected buckets.
- Cross-sample lineages are meaningful only when scope and donor metadata are correct.
