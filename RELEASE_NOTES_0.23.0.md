# Swig 0.23.0

This release adds optional fixed-tree phylogenetic UCA inference to the selected-lineage workbench.

- The exact current nucleotide alignment is used, including accepted Alivibe edits and row/column deletion. The ordinary N-masked germline guide is removed before a new observed-only double-precision FastTree is inferred.
- Directed Felsenstein cavity messages are computed in both directions on every observed-tree edge. Candidate UCA attachment position and pendant length are screened, evaluated with the complete recombination model, and refined locally.
- Automatic character selection uses ordinary nucleotide GTR4 when the retained observed alignment is ungapped. An explicit reversible `A/C/G/T/gap` GTR5 fixed-alignment approximation is activated only when an observed gap is present; `N` remains nucleotide ambiguity. No 5-mer or other context-dependent likelihood is used.
- The default nucleotide exchangeabilities are a documented reversible, context-averaged projection of the human HS5F substitution tables. The context is discarded before likelihood evaluation; gap rates are separate Swig defaults.
- A sparse/factorized V–N–D–N–J HMM reruns a broad V/J candidate screen for every lineage member, unions those hypotheses with all reported co-optimal/near-tied calls, and retains every active D. It marginalizes trim points, N additions, D identity and endpoints, J entry, and zero through three D segments by default.
- Repeated-D transitions use probability hubs rather than a dense all-D-to-all-D matrix. D states are restricted to a conservative junction window. A representative human-IGH stress shape with 48 V, 48 D, 15 J, 350 columns, and up to three D segments takes about 0.14 seconds per marginal placement and 0.45 seconds for a posterior/Viterbi pass in the release test environment.
- The result distinguishes one joint MAP Viterbi UCA/path from the sitewise marginal consensus and per-column `A/C/G/T/gap` probabilities. Nearby placement hypotheses can be locally marginalized, with an effective-placement count reported.
- The placed Newick is rooted at the inferred UCA. Its named sequence carrier has length zero; the complete inferred pendant length is assigned to the observed-tree side rather than split.
- Results export as aligned UCA FASTA, posterior TSV, placed Newick, complete JSON, and the standard publication-oriented tree SVG. Full UCA state and settings are retained by Save session and bound to the exact lineage set and alignment fingerprint.
- [`public/PHYLO_UCA_INFERENCE.md`](public/PHYLO_UCA_INFERENCE.md) documents the full likelihood equations, HMM topology, candidate union, gap semantics, priors, placement search, local marginalization, defaults, complexity, exports, and limitations.

The observed tree is fixed and gap evolution is not modeled with a generative indel process. This release is an empirical-Bayes approximation intended for interactive sensitivity analysis, not a replacement for full phylo-HMM MCMC.
