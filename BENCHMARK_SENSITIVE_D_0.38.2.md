# AER-R Sensitive-D and boundary audit — Swig 0.38.2

## Scope

Swig 0.38.2 adds the opt-in AER-R `sensitive_d` calling profile and corrects junction-end materialization for both V and J. The default assignment strategy and default truth-optimized calling profile are unchanged. `r_optimized` and `sensitive_d` both require `aer_robust`.

The boundary correction follows one invariant: an examined D candidate cannot mutate an independent V or J alignment. Each complete V-(D)-J alternative owns temporary boundaries, and no boundary is committed until that complete alternative wins. AER-R compares the ordinary allocation (retain D, clip an overlapping V/J) with the symmetric allocation (retain V/J, clip D). The symmetric allocation must improve both the matching same-D allocation and the current global partition strictly. The sole exact-tie exception restores one matching V nucleotide when the independent V alignment still has at least ten unaligned 3′ reference nucleotides; that condition was selected on development records and made no errors in either held-out simulation partition. There is no analogous automatic one-base J preference: it increased exact J starts but worsened overall J-start MAE.

In high-SHM reads, a local alignment can omit a five-nucleotide junction-facing block whose substitution score is exactly zero under `+2/−3` (three matches and two mismatches). For reads with selected V identity at most 0.88, 0.38.2 materializes such a block at V end or J start only when no shorter prefix has positive score and the block lies entirely between the selected segments. The bound comes from the actual selected D call, or from the opposite selected segment when no D is called; provisional D candidates are not consulted. This rule was selected on development records and improved held-out V-end and J-start absolute error without changing the low-SHM endpoint distribution.

## Sensitive-D setting

Sensitive-D starts from the complete 0.38.1 R-optimized tuple and changes one search gate:

| Parameter | R-optimized | Sensitive-D |
|---|---:|---:|
| V scoring | `+2/−3/−9/−1` | same |
| D scoring | `+2/−3/−13/−1` | same |
| Minimum consecutive exact D run | 5 nt | **4 nt** |
| D candidates | 2 | same |
| D-presence cost | 12; relaxed to 10 for the existing strong-score/template-support rule | same |
| D ambiguity window | same-span candidates within 1 raw point | same |
| J scoring | `+2/−3/−17/−2` | same |

The lower floor changes candidate admission, not the final D-presence decision. It is intended for mutation-interrupted D tracts and knowingly accepts a small increase in false D calls on no-D truth records.

## Evaluation

The evaluation used the supplied 20,000-record lower-SHM/IgM-like and 20,000-record IgG-like simulations with the deduplicated union of bundled rhesus and cynomolgus KIMDB 1.1 V/D/J references. Identical full-reference sequences were collapsed to one class. If a report contained `K` distinct predicted classes and the truth class was present, call accuracy received `1/K`; otherwise it received zero. Absence of D was a class. The 19 lower-SHM and 22 IgG-like tandem-D truths were excluded from the single-D primary score.

### Fair V, D, and J call accuracy

| Dataset | Version and method | V | D | J |
|---|---|---:|---:|---:|
| Lower-SHM | 0.37.6 plain AER | 94.6675% | 76.3150% | 91.2600% |
| Lower-SHM | 0.38.2 AER-R R-optimized | **94.7154%** | 76.8958% | **91.4748%** |
| Lower-SHM | 0.38.2 AER-R Sensitive-D | **94.7154%** | **76.8977%** | **91.4748%** |
| IgG-like | 0.37.6 plain AER | 87.3450% | 67.6794% | 85.9750% |
| IgG-like | 0.38.2 AER-R R-optimized | **87.4435%** | 68.7936% | **86.2114%** |
| IgG-like | 0.38.2 AER-R Sensitive-D | **87.4435%** | **68.8658%** | 86.2089% |

### Sensitive-D recovery/specificity trade-off

| Dataset | Profile | Truth-D call set contains truth | Extra correct truth-D sets vs R | No-D false calls | Extra no-D false calls vs R |
|---|---|---:|---:|---:|---:|
| Lower-SHM | R-optimized | 16,375 / 18,942 (86.4481%) | — | 182 / 1,039 | — |
| Lower-SHM | Sensitive-D | **16,383 / 18,942 (86.4903%)** | **+8** | 187 / 1,039 | +5 |
| IgG-like | R-optimized | 14,964 / 18,925 (79.0700%) | — | 214 / 1,053 | — |
| IgG-like | Sensitive-D | **14,994 / 18,925 (79.2285%)** | **+30** | 220 / 1,053 | +6 |

### Full endpoint losses

Endpoint metrics include only records whose reported call set contains the truth sequence class. Exact and within-one are proportions; MAE is nucleotides.

| Dataset | Version/profile | V-end MAE | V exact | V within 1 | J-start MAE | J exact | J within 1 |
|---|---|---:|---:|---:|---:|---:|---:|
| Lower-SHM | 0.38.1 R-optimized | 0.4842 | 71.3409% | 89.9011% | 0.5244 | 70.0183% | 88.8121% |
| Lower-SHM | 0.38.2 R-optimized | **0.4836** | **71.3967%** | **89.9062%** | **0.5244** | **70.0183%** | **88.8121%** |
| Lower-SHM | 0.38.2 Sensitive-D | **0.4835** | **71.4018%** | **89.9062%** | 0.5246 | 70.0081% | 88.8070% |
| IgG-like | 0.38.1 R-optimized | 0.7987 | 61.4784% | 82.7253% | 0.8351 | 60.6852% | 81.8980% |
| IgG-like | 0.38.2 R-optimized | **0.7866** | **61.6537%** | **82.8952%** | **0.8205** | **60.8481%** | **82.1554%** |
| IgG-like | 0.38.2 Sensitive-D | **0.7862** | **61.6537%** | **82.8952%** | **0.8201** | **60.8481%** | **82.1554%** |

The change was selected on the complete loss distribution, not on ±1 alone. In the IgG-like R-optimized output, V errors with absolute magnitude at least 5 fell from 526 to 496 and J errors from 644 to 601. The V −5 count fell from 157 to 129 while +5 rose from 58 to 68; the net MAE and tail counts both improved. J +5 fell from 160 to 120, although J −5 rose from 87 to 96. On lower-SHM R-optimized, J's complete endpoint distribution is byte-for-byte unchanged; V −1 falls from 517 to 507 and exact endpoints rise from 14,072 to 14,083.

## Supplied regressions

The supplied read `sim_igsim_indelp0.003_n20000__assign=assign__sim=all_14597` contains a 15-column D alignment with three substitutions and a longest exact run of four:

```text
query     ATACTGTGGCGACAG
reference ATACAGTGGGTACAG
```

R-optimized correctly continues to omit it because of the five-base candidate floor. Sensitive-D reports `IGHD5-27*01,IGHD5-32*01,IGHD5-32*01_S0263` at query positions 430–444. On the supplied 760-read historical missed-D fixture, fair D accuracy rises from 78.3564% to 78.4881%, truth-in-reported-set rises from 87.3518% to 87.4835%, and no-D specificity remains 75%.

## Timing

One full one-worker pass per configuration used the same current harness, initialized WASM runtime, references, and 20,000-record input. Timing starts after database initialization. CPU time is the more stable comparison because the lower-SHM and IgG-like chains ran concurrently, while configurations within each dataset ran sequentially.

| Dataset | Version and method | CPU seconds | CPU reads/s | Change vs 0.37.6 plain AER | Change vs 0.37.6 AER-R |
|---|---|---:|---:|---:|---:|
| Lower-SHM | 0.37.6 plain AER | 27.075 | 738.7 | — | — |
| Lower-SHM | 0.37.6 AER-R | 28.257 | 707.8 | +4.4% | — |
| Lower-SHM | 0.38.2 R-optimized | 29.990 | 666.9 | +10.8% | **+6.1%** |
| Lower-SHM | 0.38.2 Sensitive-D | 30.519 | 655.3 | +12.7% | **+8.0%** |
| IgG-like | 0.37.6 plain AER | 49.867 | 401.1 | — | — |
| IgG-like | 0.37.6 AER-R | 51.083 | 391.5 | +2.4% | — |
| IgG-like | 0.38.2 R-optimized | 54.115 | 369.6 | +8.5% | **+5.9%** |
| IgG-like | 0.38.2 Sensitive-D | 51.792 | 386.2 | +3.9% | **+1.4%** |

There is therefore a measured slowdown, not a speed-neutral result: about 1–8% versus the same v0.37.6 AER-R strategy in these single passes, or 4–13% versus v0.37.6 plain AER. The apparent IgG-like Sensitive-D speed advantage over R-optimized is treated as run noise, not as a performance claim; repeated medians would be required for a finer comparison.

## Limitations

- Both simulations use the same supplied simulator/reference regime; the held-out partitions prevent direct threshold reuse but are not external biological validation.
- Sensitive-D intentionally trades a small number of no-D false positives for additional truth-D recovery.
- The endpoint repair is conservative and does not eliminate every −1, +1, or −5 feature. It was retained because aggregate exact, within-one, MAE, and tail losses improved without moving the low-SHM J distribution.
- Very short D tracts and severely mutated junctions remain intrinsically ambiguous. Study-critical calls should be checked against suitable controls and an independent workflow.
