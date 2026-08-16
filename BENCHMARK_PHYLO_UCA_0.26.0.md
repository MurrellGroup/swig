# Phylogenetic UCA placement benchmark (0.26.0)

## Reproduction target

- Input: user-supplied `swig-study.lineage-40.fasttree-alignment.fasta`
- Alignment: 112 observed sequences × 433 nucleotide columns
- Observed tree: native double-precision FastTree 2.1.11, nucleotide GTR
- References: KIARVA genomic coding FASTA fetched on 2026-08-16 (561 IGHV, 51 IGHD, 12 IGHJ records before Swig candidate filtering)
- HMM template leakage setting: exactly `0`; the emission implementation uses a `1e-9` numerical floor
- Search breadth: default six V/J-screen-ranked full-HMM edges; five screen positions per edge
- Runtime: Node 24.19.0, x86-64; reference download and FastTree construction excluded from the three inference timings

Reproduce with:

```bash
node --experimental-strip-types scripts/benchmark-phylo-uca-lineage40.ts \
  swig-study.lineage-40.fasttree-alignment.fasta /path/to/FastTreeDbl
```

An optional fourth argument selects one route: `maximum-likelihood`, `grid-marginalization`, or `gibbs-mh`.

## Results

| Route | Wall time | Best edge / fraction | UCA branch | Best full-HMM log L | Points/draws used |
|---|---:|---|---:|---:|---:|
| Conditional ML | 62.24 s | edge 49 / 0.000 | 0.00559479 | -3773.0573 | 1 |
| Grid marginalization | 120.61 s | edge 49 / 0.500 | 0.00706443 | -3773.7793 | 12 retained of 234 evaluated |
| Continuous Gibbs/MH | 46.56 s | edge 48 / 0.824 | 0.00517141 | -3768.6432 | 60 retained draws |

The ML timing includes continuous full-HMM coordinate optimization on every one of the six admitted edges. Grid mode evaluated `6 × 3 × 13 = 234` complete HMM likelihoods, then computed exact posteriors at the 12 leading quadrature points. Gibbs/MH ran 160 exact HMM Gibbs updates, discarded 40, retained every second draw, and performed four cheap conditional placement moves per Gibbs update.

The best Gibbs/MH draw is on edge 48, whereas the restricted default ML/grid runs report edge 49. This is expected behavior, not an inconsistency: the V/J screen is an explicit speed/coverage approximation for ML/grid, while Gibbs/MH global proposals retain nonzero support on every edge. Setting **Full-HMM edges** to zero searches every edge in ML/grid; intermediate values provide a less expensive larger search.

## The reported D-extension failure

With leakage effectively zero, the corrected near-zero branch treatment extends the `IGHD3-3*01` source through the previously missing six-base region. Marginalized source occupancy was:

| Alignment column | Conditional ML | Grid | Gibbs/MH |
|---:|---:|---:|---:|
| 346 | 0.990 | 0.982 | 0.983 |
| 347 | 0.983 | 0.967 | 0.967 |
| 348 | 0.963 | 0.949 | 0.917 |
| 349 | 0.904 | 0.893 | 0.850 |
| 350 | 0.726 | 0.724 | 0.717 |
| 351 | 0.726 | 0.723 | 0.717 |
| 352 | 0.724 | 0.721 | 0.717 |
| 353 | 0.718 | 0.716 | 0.683 |
| 354 | 0.702 | 0.699 | 0.683 |
| 355 | 0.651 | 0.647 | 0.567 |
| 356 | 0.494 | 0.490 | 0.333 |

All three routes keep at least 10% `IGHD3-3*01` occupancy through column 356. The MCMC column estimates are visibly coarser because they are frequencies from 60 retained joint draws; longer chains reduce that Monte Carlo noise.

The old branch grid (`0`, `0.075`, `0.3`, with a few coarse refinements) could behave as if the UCA-to-tree branch were zero and suppress the mutation-supported D extension. The new conditional-ML optimum is near `0.0056`; the explicit logarithmic grid includes `0.00706443`, close enough to recover the same HMM explanation.

## Gibbs/MH mixing

- Branch random-walk acceptance: 52.8%
- Within-edge position acceptance: 64.9%
- Global-edge acceptance: 4.6%
- Accepted edge switches: 4
- Branch-length ESS: 20.1 of 60 retained draws
- Marginal-log-target ESS: 13.3 of 60 retained draws
- Distinct pendant lengths visited across the 160-iteration trace: 108

The branch values are continuous exact GTR evaluations. They are not members of, interpolations over, or indices into the explicit grid used by grid marginalization.

## Tree-likelihood reuse

Every route constructs the observed tree's two directed cavity messages per edge once. A new attachment or pendant proposal combines the two cached half-edge messages locally and applies the exact GTR transition for its floating-point distances. Conditional ML/grid then run one complete HMM marginal evaluation for that local surface. In Gibbs/MH, only the single Gibbs update per iteration runs the HMM; its subsequent MH placement proposals condition on the sampled UCA and therefore need no HMM pass and no observed-tree pruning traversal.
