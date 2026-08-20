# Swig 0.37.5

## Joint V/D/J boundary adjudication

- Fixes a D false negative in which a gappy terminal V local alignment consumed 13 bases of a perfect simulated D tract, leaving only a five-base V/J-bounded D window.
- Weak-D cases with a junction-facing V or J gap now receive one narrowly bounded rescue. A 10-nt exact D seed is required before any second dynamic-programming refinement is run.
- Candidate V/D/J boundaries are scored jointly. The retained V and J alignments are clipped and rescored, and a rescued D is accepted only if the disjoint V+D+J affine score is strictly better than the original V+(weak D)+J score.
- The supplied KIMDB regression now recovers the complete 17-nt simulated tract with AER, RIAT-MP, and Standard SwiftIG under both truth-optimized and agreement settings. The observed 18-nt perfect span is co-optimal for `IGHD5-27*01`, `IGHD5-32*01`, and `IGHD5-32*01_S0263`; the agreement/balanced profile's three-call retention therefore reports all three, including the simulated allele.

The ordinary successful D path is unchanged and performs no additional alignment. Web and CLI use the same rebuilt SwiftIG WASM binary.
