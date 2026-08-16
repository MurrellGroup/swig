# Optional Double-D / VDDJ evidence screen

## Scope

The screen runs only after an ordinary SwiftIG V(D)J annotation and only for loci with D segments (IGH, TRB, and TRD). It is off by default. A positive result is retained in a sparse sidecar and interactive overlay; it does not silently rewrite the ordinary assignment table. This separation is intentional because tandem-D inference from short, mutable junctions is substantially less certain than ordinary V/J assignment.

## Candidate construction

The search window is the oriented query interval from the end of the selected V alignment to the start of the selected J alignment. In `long span` mode it is skipped unless this interval reaches the configured threshold (40 nt by default); `all` mode screens every eligible interval.

All exact canonical \(k\)-mers from locus-compatible D references are indexed; \(k=11\) by default. An exact seed hit is ungapped-extended left and right with match `+2`, mismatch `−2`, retaining the highest-scoring extension in each direction. Duplicate gene/diagonal seed hits are collapsed. Two hits form a candidate only when they are ordered and non-overlapping on the query.

For a pair \((D_1,D_2)\), the pair score is the sum of its two extension scores. It must exceed the best supported single-D extension by the configured minimum score gain (8 by default). Equal-score/equal-coordinate allele labels remain co-optimal.

## Pseudo-tandem veto

Two seeds from one D allele with a mutation/deletion between them can mimic two D segments. To reject that explanation, Swig takes the complete query span from the start of \(D_1\) through the end of \(D_2\), removes a configurable delta length (5 nt by default), and asks how closely every remaining substring can match any single locus-compatible D substring of the same length. If the minimum Hamming distance is at most the configured maximum (3 mismatches by default), the two-D pair is vetoed.

Candidates are tested in evidence order. Alternatives within two raw score units of the best are retained; lower-scoring combinations are not materialized. Reported fields include D1/D2 calls and coordinates, score gain, pseudo distance, V–J span, and NP1/NP2/NP3 intervals.

## Literature relationship

The 11-nt exact seed starting value and single-D mimic concern are informed by IgScout/tandem-D work, but this is a **custom SwiftIG evidence screen**, not an implementation of IgScout. See Safonova and Pevzner, [V(DD)J recombination is an important and evolutionarily conserved mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC7605257/), and Briney et al., [Frequency and genetic characterization of V(DD)J recombinants](https://pmc.ncbi.nlm.nih.gov/articles/PMC3449247/).

## Explorer and downstream use

The Double-D explorer filters the sparse positive sidecar; opening a row overlays the tandem interpretation on the original AIRR record. Repertoire selection may explicitly require or exclude positives. Where both D alignments can be safely projected, lineage-root construction can use V–D1–D2–J rather than allowing the unchanged baseline single-D composite to erase D2.

## Limitations

- A positive is supported VDDJ evidence, not a definitive genomic recombination call.
- Exact seeds lose sensitivity under mutation or sequencing error; shorter seeds reduce specificity.
- The ungapped seed extension is deliberately fast and is not a full pair-HMM.
- The pseudo-tandem veto is a bounded alternative-explanation test, not a full probability ratio over recombination histories.
- Defaults are starting values and should be checked against negative controls and organism-specific D references.
