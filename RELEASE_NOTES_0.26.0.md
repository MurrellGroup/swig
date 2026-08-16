# Swig 0.26.0

> Historical release record. The 0.26.0 terminal-D explanation was superseded by the biological HMM audit in 0.27.0; see `RELEASE_NOTES_0.27.0.md` and `BENCHMARK_PHYLO_UCA_0.27.0.md`.

## Phylogenetic UCA placement routes

- Added a fast conditional-maximum-likelihood route, now the default. It continuously optimizes within-edge attachment position and UCA pendant length under the complete recombination-HMM marginal likelihood; it does not apply placement priors or marginalize these coordinates.
- Replaced the old coarse pendant grid in marginalization mode with an explicit zero-plus-logarithmic grid and quadrature cell weights. The settings show the exact numerical branch-length list before a run, and the result stores it.
- Added an exact Metropolis-within-Gibbs alternative. Each iteration uses forward-filter/backward-sample for a coherent HMM path and UCA sequence, followed by cheap placement MH updates against cached directed tree messages.
- Gibbs/MH pendant length and within-edge attachment position are continuous. No discretized branch lengths, grid lookup, or branch-grid interpolation is used.
- Added user controls for inference route, ML optimization, grid extent/resolution, MCMC length/burn-in/thinning/seed, MH steps and proposal scales, global-edge moves, and placement priors.
- Added sampler trace plots, acceptance rates, edge-switch counts, and branch/log-target ESS summaries.

## Near-zero branch correction

- The V/J-only initializer now optimizes both attachment and pendant length continuously before the complete HMM stage.
- Conditional ML optimizes every edge admitted by the user-selected search breadth, avoiding the earlier failure in which a one-point full-HMM score discarded an edge with a narrow near-zero optimum.
- Template leakage accepts zero in the UI and evaluates it at a `1e-9` numerical floor.
- On the supplied 112-sequence lineage-40 alignment with KIARVA and zero leakage, all three routes recover `IGHD3-3*01` occupancy through the six-base extension that the old coarse grid suppressed. See `BENCHMARK_PHYLO_UCA_0.26.0.md`.

## Verification

- Added a Monte Carlo regression showing that 5,000 exact FFBS samples reproduce the analytic forward-backward nucleotide posterior.
- Added end-to-end tests for the exact reported branch grid and for off-grid continuous Gibbs/MH branch values and mixing diagnostics.
- The observed tree's directed half-edge messages remain a one-time calculation; placement evaluations never rerun full-tree pruning.
