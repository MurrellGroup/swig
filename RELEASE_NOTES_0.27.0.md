# Swig 0.27.0 release notes

## Phylogenetic UCA HMM audit

- J entry is now possible only on a concrete projected J nucleotide. The old unresolved early-J route that could behave like absorbing NT has been removed.
- V3, D5, D3, and J5 trimming use normalized finite-support deletion distributions. D3 is a terminal-trim prior rather than a repeated per-retained-base penalty.
- D alleles receive equal prior mass before their within-allele 5′ trim distribution.
- First-D use, additional-D use, non-empty-N probability, one-nucleotide N mass, phased N-length tail, N nucleotide composition, junction search support, leakage, and gap behavior are explicit parameters.
- Template leakage defaults to exactly zero. Mutations between a templated UCA base and the tree attachment remain possible through the inferred GTR pendant branch.
- The default human-IGH starting moments were cross-checked against public IGoR/OLGA parameters and the model structure against partis/ham and linearham. These are visible starting values, not universal or per-allele biological truth.
- A deterministic prior-predictive panel reports generated trim lengths, retained D length, N lengths, junction span, and D-count frequencies from the current settings and active D lengths.

## Inference and mixing

- Continuous Gibbs/MH remains the default, with 320 iterations, 80 burn-in, thin 2, and exact-zero leakage.
- Full-HMM initializer scores are reused in a full-support collapsed independence proposal for edge, continuous within-edge position, and continuous pendant length.
- Exact collapsed refreshes run every three iterations by default, chosen from the supplied lineage by ESS per wall-clock second.
- The collapsed move evaluates the HMM-marginal likelihood first and samples a new path/UCA only after acceptance. Ordinary MH updates remain conditional on the current sampled UCA.
- The mixing panel now reports collapsed acceptance, branch and log-target ESS, ESS per second, and mean time for an FFBS draw, collapsed marginal evaluation, and fixed-UCA proposal.
- Runs with branch or log-target ESS below 20 carry a visible warning.
- Gibbs/MH results report the posterior over the number of D segments by counting retained joint paths, with no additional HMM evaluations.
- Conditional ML and grid marginalization remain available. Grid mode shows and stores its exact pendant-length grid; Gibbs/MH does not discretize branch length or attachment position.

## Performance

- Directed all-to-all tree messages are computed once; placement moves use only local cached half-edge messages.
- HMM state catalogs, D entry priors, D endpoint hazards, and fixed D emissions are constructed once and reused.
- D-state computation and backward storage are sparse outside the user-visible junction window.
- Gibbs backward typed arrays are reused rather than rebuilt for every draw.
- Accepted placement surfaces are carried into the next iteration.
- Rejected collapsed proposals no longer pay for a useless conditional HMM/UCA draw.

The supplied lineage-40 test completes the 320-iteration default in about 62 seconds in the benchmark container. See `BENCHMARK_PHYLO_UCA_0.27.0.md` for route timings, ESS-per-second tuning, seed replication, and the zero-leakage D-extension regression.

## Controls and output

- Every recombination-HMM and Gibbs/MH control is editable, including after a completed run; Additional-D accepts the full 0–1 range without being overwritten by the previous result.
- Reset restores the complete UCA option set to version 0.27.0 defaults.
- Complete and currently visible multi-track SVG exports retain the corresponding UCA posterior logo below the tracks.
- Complete JSON/session output uses result schema 6 and includes the parameter set, timing diagnostics, D-count posterior, and sampling trace. Older saved HMM options are migrated through retained legacy aliases.

## Navigation safety

- Browser back gestures and history navigation are intercepted while loaded data, results, pasted/file input, a pending session, or an active run exists.
- A custom dialog defaults to **Stay on this page** and requires an explicit **Leave anyway** action.
- Reload, tab close, external navigation, and browser shutdown use the browser's native unsaved-work warning.
- Horizontal overscroll navigation is disabled in CSS where the browser honors it.
- Starting a replacement analysis while a result is loaded asks for confirmation.

## Scientific caveat

The observed lineage tree is still fixed. Swig marginalizes or samples UCA placement, pendant length, recombination path, and UCA sequence according to the selected route, but it does not integrate over observed-tree topology, the fitted internal branch lengths, or per-allele learned recombination parameters.
