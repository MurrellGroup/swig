# Phylogenetic-UCA audit and benchmark — Swig 0.27.0

## Scope

This is a targeted regression and performance check on the user-supplied `swig-study.lineage-40.fasttree-alignment.fasta`, using the KIARVA germline collections requested for the analysis. It is not a general accuracy benchmark and the wall times are not portable across browsers or machines.

- observed lineage sequences: 112
- aligned nucleotide columns: 433
- alignment gaps: 7
- active KIARVA references: 561 V, 51 D, 12 J
- tree: a fresh double-precision FastTree inferred after removing the germline-guide row
- timing host: Node 24.19.0, Linux x86-64
- template leakage: exactly 0
- default Gibbs/MH output: 320 iterations, 80 burn-in, thin 2, 120 retained joint draws

With the supplied alignment and a local `FastTreeDbl`, the three routes can be rerun with:

```bash
node --experimental-strip-types scripts/benchmark-phylo-uca-lineage40.ts \
  swig-study.lineage-40.fasttree-alignment.fasta /path/to/FastTreeDbl
```

An optional fourth argument restricts the run to `maximum-likelihood`, `grid-marginalization`, or `gibbs-mh`. The script downloads the requested KIARVA V/D/J collections at run time.

## Biological regression

The pre-audit automaton could enter a J state before a concrete projected J nucleotide and let that unresolved J projection emit an effectively uniform base. That created an absorbing, NT-like competitor without the N-duration prior. It also encoded the D 3′ parameter as a repeated per-retained-base exit process rather than as a prior over terminal trimming.

Version 0.27.0 requires J entry on a concrete J nucleotide, uses normalized finite-support deletion distributions at V3, D5, D3, and J5, gives D alleles equal prior mass before within-allele trimming, and represents positive N lengths with a one-base atom plus a short phased tail. Exact-zero template leakage does not prohibit a template/tree disagreement: a substitution is still available through the exact GTR transition on the inferred UCA-to-tree branch.

For the disputed terminal `IGHD3-3*01` extension, all three inference routes now agree at zero leakage:

| Route | Total wall time | D occupancy at column 350 | D occupancy at column 356 | D occupancy at column 357 |
|---|---:|---:|---:|---:|
| Conditional ML | 37.14 s | 99.34% | 85.81% | 0.58% |
| Explicit grid | 74.39 s | 99.19% | 85.71% | 1.11% |
| Gibbs/MH, seed 1729 | 62.10 s | 99.17% | 87.50% | 0.00% |
| Gibbs/MH, seed 99173 | 62.52 s | 100.00% | 90.00% | 0.83% |

The retained Gibbs paths had `P(number of D segments = 1) = 1.00` in both full default runs. This posterior is lineage-specific, not a claim that one D is universal.

## Placement routes

### Conditional ML

- wall time: 37.14 s
- HMM search: 30.34 s
- final posterior: 4.42 s
- best placement: edge 49 at its `internal49` endpoint
- UCA pendant length: 0.002635 substitutions/character
- MAP recombination call: `IGHV3-30*18/IGHV3-30-5*01 · IGHD3-3*01 · IGHJ3*02`

Each scalar objective call uses the full HMM, but the observed-tree pruning likelihood is not recomputed. It uses the directed all-to-all half-edge messages cached once for the fixed tree.

### Explicit grid marginalization

- wall time: 74.39 s
- evaluated points: 234
- posterior points retained for character/HMM output: 12
- quadrature mass retained by those points: 99.069%
- effective grid points: 6.72
- reported pendant length: 0.0070644 substitutions/character

The complete zero-plus-logarithmic pendant grid is shown before inference and stored with the result.

### Continuous Gibbs/MH

Both within-edge attachment position and UCA pendant length are continuous. No attachment or branch-length grid is used.

| Seed | Wall time | Sampling time | Branch ESS / 120 | Log-target ESS / 120 | Branch ESS/s | Log-target ESS/s | Collapsed accepted | Edge switches |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1729 | 62.10 s | 54.79 s | 42.41 | 36.92 | 0.774 | 0.674 | 1 / 106 | 2 |
| 99173 | 62.52 s | about 55 s | 64.99 | 64.65 | 1.182 | 1.176 | 2 / 106 | 2 |

The two seeds settle on the same principal attachment region and give compatible branch lengths and D occupancy. Mixing is adequate for this short default run but not spectacular; the seed-to-seed ESS spread is why the UI reports ESS and warns below 20 instead of presenting the iteration count as a convergence guarantee.

## Collapsed-refresh tuning by ESS per second

The benchmark compared how often to perform the expensive placement proposal that marginalizes the complete HMM. Each ordinary iteration also performs four cheap fixed-UCA placement proposals.

| Exact collapsed refresh | Sampling time, 240 iterations | Branch ESS/s | Log-target ESS/s |
|---|---:|---:|---:|
| every 20 iterations | 42.0 s | 0.246 | 0.083 |
| every 5 iterations | 45.8 s | 0.624 | 0.094 |
| every 2 iterations | 54.1 s | 0.608 | 0.075 |
| every 3 iterations, seed 1729 | 48.2 s | 0.846 | 0.747 |
| every 3 iterations, seed 99173 | 48.8 s | 1.149 | 1.076 |

The default is therefore every three iterations—approximately one collapsed proposal per twelve ordinary fixed-UCA proposals—not the setting with the most proposals per second.

Representative costs on this lineage were:

- exact full-HMM conditional draw: about 133–137 ms
- exact collapsed full-HMM marginal proposal: about 109–119 ms
- fixed-UCA continuous placement proposal: about 0.14 ms

Because a collapsed proposal is roughly 800 times the cost of a fixed-UCA proposal here, marginalizing the HMM at every MH step is a poor use of time unless it produces an enormous ESS improvement; it did not in this test.

## Why the delayed latent draw is valid

Let `theta` be tree edge, continuous within-edge position, and UCA pendant length; let `z` be the HMM path and UCA sequence. The collapsed proposal accepts or rejects `theta'` using the exact marginal target

`p(y | theta') p(theta') = sum_z p(y,z | theta') p(theta')`

and the exact forward/reverse proposal-density ratio. No proposed `z'` appears in that acceptance probability. If `theta'` is accepted, Swig immediately samples

`z' ~ p(z | y, theta')`.

If it is rejected, the current `(theta,z)` is retained. This is a valid partially collapsed kernel because after either outcome the latent state has the correct conditional distribution for the retained placement. Drawing `z'` before acceptance would also be possible with a different augmented-state construction, but would waste nearly every rejected FFBS draw. Accepting a collapsed placement and keeping the old `z` would be invalid; the implementation explicitly does not do that. A regression test also checks that the total conditional-draw count is exactly `iterations + accepted collapsed proposals`.

## Runtime changes

- one immutable V/D/J HMM catalog is reused across ML, grid, and Gibbs/MH evaluations;
- fixed D emissions and finite endpoint hazards are precomputed;
- D entry scans and D backward storage are omitted outside the bounded junction window;
- Gibbs backward storage uses typed sparse rows and is reused between iterations;
- tree likelihood surfaces use cached directed half-edge messages, never a fresh whole-tree pruning pass;
- a collapsed proposal first computes only the rolling marginal likelihood and performs FFBS only after acceptance;
- the current surface is carried across accepted proposals rather than recomputed at the start of the next iteration.

The audited 240-iteration sampling loop is substantially faster than the pre-optimization implementation while retaining the same statistical target. Final browser timing will still depend on JavaScript engine, CPU, active reference breadth, alignment length, D-window width, and number of D states.

## External calibration cross-checks

Defaults are compact moment matches and regularizers, not imported per-allele fits. The public human-IGH IGoR/OLGA parameterization informed trim, N-length, non-empty-N, and N-composition moments. Partis/ham was used as a structural cross-check: it supports concrete germline states, allele-specific deletion distributions, observed N-length distributions, and tiered aggregation when data are sparse. Swig deliberately retains a smaller shared-tail model so every parameter remains visible and inference stays fast. Linearham is a fuller joint phylo-HMM precedent; Swig still conditions on one observed tree.

Relevant sources are listed in `public/PHYLO_UCA_INFERENCE.md`.
