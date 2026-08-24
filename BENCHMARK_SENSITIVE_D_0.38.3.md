# AER-R Sensitive-D decision audit — Swig 0.38.3

## Scope

Swig 0.38.3 promotes the accuracy-improving four-base D-candidate floor from the former Sensitive-D profile into AER-R R-optimized, then redefines Sensitive-D as a lower final joint D/no-D threshold. The default assignment strategy and truth-optimized calling profile remain unchanged. Both optimized profiles remain opt-in and require `aer_robust`.

The change separates two questions that had previously been conflated:

1. **Is there enough D-like signal to admit a candidate?** Both profiles require an ordinary exact run of at least 4 nt or AER-R's existing conservative distributed-D alternative.
2. **Does that admitted candidate improve a complete V–D–J partition enough to report D?** R-optimized charges 12 points, relaxed to 10 for raw D score at least 20 or exact support from at least two distinct D template sequences. Sensitive-D applies a fixed 10-point cost.

Sensitive-D therefore does not synthesize a D label without an alignment and does not force D on every read. Under the `+2/−3` tuple, a bare four-base exact seed scores only 8 and cannot clear the fixed 10-point joint cost by itself. Sensitive-D asks for the best D-bearing complete partition at a lower, explicit decision cost after signal admission.

| Parameter | R-optimized | Sensitive-D |
|---|---:|---:|
| V scoring | `+2/−3/−9/−1` | same |
| D scoring | `+2/−3/−13/−1` | same |
| Minimum consecutive exact D run | 4 nt | same |
| Distributed-D fallback | ≥14 matches across 16 aligned bases plus aggregate-score gate | same |
| D candidates | 2 | same |
| Weak-evidence D-presence cost | 12 | **10** |
| Strong/template-supported D-presence cost | 10 | **10** |
| D ambiguity window | same span, within 1 raw point | same |
| J scoring | `+2/−3/−17/−2` | same |

## Evaluation and scoring

The evaluation used the supplied 20,000-record lower-SHM/IgM-like and 20,000-record IgG-like simulations with the deduplicated union of bundled rhesus and cynomolgus KIMDB 1.1 V/D/J references. Identical full-reference sequences were collapsed to one class. If a report contained `K` distinct predicted sequence classes and the truth class was present, accuracy received `1/K`; otherwise it received zero. Absence of D was an explicit class. The 19 IgM-like and 22 IgG-like tandem-D truths were excluded from the primary single-D score.

The same deterministic ID hash used previously was also used to inspect development, validation, and test partitions. Because all partitions were examined while choosing the released threshold, those partitions are stability checks, not an untouched external test.

## Threshold selection

The former Sensitive-D setting—4-nt candidate floor with R-optimized's 12-to-10 evidence-conditioned decision cost—improved fair D accuracy over Swig 0.38.1 R-optimized on both complete simulations. It is therefore the new main R-optimized setting.

The subsequent joint cost was then varied without changing search, alignments, ambiguity reporting, or V/J scoring:

| Dataset | Joint D cost | Fair D accuracy | Truth-D detection | D-positive reports | No-D false calls |
|---|---|---:|---:|---:|---:|
| IgM-like | R: 12→10 | 76.8977% | 96.9644% | 18,554 | 187 |
| IgM-like | fixed 11 | **76.9330%** | 97.8672% | 18,774 | 236 |
| IgM-like | **fixed 10** | 76.7019% | **98.4373%** | **18,932** | 286 |
| IgG-like | R: 12→10 | 68.8658% | 95.8890% | 18,367 | 220 |
| IgG-like | fixed 11 | **68.9229%** | 96.9775% | 18,640 | 287 |
| IgG-like | **fixed 10** | 68.6664% | **97.7701%** | **18,856** | 353 |

Fixed 11 looked slightly better on each complete dataset, but it reduced IgG-like fair D accuracy relative to R on both the validation partition (−0.0716 percentage points) and the test partition (−0.1385 points). It was therefore not treated as a robust main-profile improvement. Fixed 10 was selected for the explicitly sensitivity-oriented profile: it produces a much larger call increase while keeping the complete-data fair D-accuracy loss to about 0.20 percentage points in each regime.

## Full V, D, and J accuracy

| Dataset | Version/profile | V | D | J |
|---|---|---:|---:|---:|
| IgM-like | 0.38.1 R-optimized | 94.7154% | 76.8958% | 91.4748% |
| IgM-like | **0.38.3 R-optimized** | **94.7154%** | **76.8977%** | **91.4748%** |
| IgM-like | **0.38.3 Sensitive-D** | **94.7154%** | 76.7019% | 91.4698% |
| IgG-like | 0.38.1 R-optimized | 87.4435% | 68.7932% | **86.2114%** |
| IgG-like | **0.38.3 R-optimized** | **87.4435%** | **68.8658%** | 86.2089% |
| IgG-like | **0.38.3 Sensitive-D** | 87.4426% | 68.6664% | 86.1914% |

## D recovery and specificity

These counts make the intended trade-off clearer than accuracy alone. Detection asks only whether any D was reported on a truth-D record; it does not claim that the D identity is correct.

| Dataset | Profile | Truth-D detected | Remaining missed | D-positive reports | Correct truth-D sets | No-D false calls | No-D specificity |
|---|---|---:|---:|---:|---:|---:|---:|
| IgM-like | R-optimized | 18,367 / 18,942 (96.9644%) | 575 | 18,554 | 16,383 | 187 / 1,039 | 82.0019% |
| IgM-like | Sensitive-D | **18,646 / 18,942 (98.4373%)** | **296** | **18,932** | 16,382 | 286 / 1,039 | 72.4735% |
| IgG-like | R-optimized | 18,147 / 18,925 (95.8890%) | 778 | 18,367 | 14,994 | 220 / 1,053 | 79.1073% |
| IgG-like | Sensitive-D | **18,503 / 18,925 (97.7701%)** | **422** | **18,856** | 15,002 | 353 / 1,053 | 66.4767% |

Relative to R-optimized, Sensitive-D adds 378 IgM-like and 489 IgG-like D-positive reports. It recovers 279 of 575 (48.5%) and 356 of 778 (45.8%) remaining truth-D misses, respectively. The added detections are deliberately not all correct identities: fair D accuracy falls by 0.1959 and 0.1993 percentage points, and false D calls increase by 99 and 133.

## V/J boundary check

Endpoint metrics include only records whose reported call set contains the truth sequence class. The decision-threshold change was checked over complete distributions, including absolute-error tails and the previously noted −5 feature.

| Dataset | Profile | V-end MAE | V exact | V within 1 | V \|error\|≥5 | V −5 | J-start MAE | J exact | J within 1 | J \|error\|≥5 | J −5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IgM-like | R-optimized | 0.48355 | 71.4018% | 89.9062% | 203 | 18 | 0.52459 | 70.0081% | 88.8070% | 243 | 78 |
| IgM-like | Sensitive-D | 0.48370 | 71.3714% | 89.9062% | 203 | 18 | 0.52487 | 69.9812% | 88.8166% | 244 | 78 |
| IgG-like | R-optimized | 0.78615 | 61.6537% | 82.8952% | 495 | 129 | 0.82008 | 60.8481% | 82.1554% | 601 | 96 |
| IgG-like | Sensitive-D | 0.78680 | 61.6185% | 82.8749% | 494 | 128 | 0.81991 | 60.8072% | 82.1956% | 601 | 95 |

The new D threshold has negligible V/J endpoint effects: the largest MAE movement is 0.00065 nt. Tail counts and −5 counts are essentially unchanged. The 0.38.2 invariant remains in force: provisional D hypotheses cannot mutate V or J boundaries; only the actual winning complete partition can trade occupied query bases between selected segments.

## Supplied regressions

- The supplied mutation-interrupted read `sim_igsim_indelp0.003_n20000__assign=assign__sim=all_14597` has a 15-column D alignment with three substitutions and a longest exact run of four. Because the four-base gate is now promoted, both 0.38.3 profiles recover the expected `IGHD5-32` sequence class.
- A permanent threshold regression uses a 16-column, 75%-identity `IGHD1-7*01` alignment with raw D score 12 and an exact four-base run. R-optimized retains no-D; Sensitive-D reports the D. This confirms that the profiles differ at the post-signal joint decision, not candidate pruning.
- A strong singleton D-like decoy in a simulated no-D junction remains no-D under both profiles. This guards against turning Sensitive-D into an unconditional D emitter.
- On the supplied 760-read historical missed-D fixture, R-optimized detects D on 97.1114% of truth-D records; Sensitive-D reaches 98.7620%. Fair D accuracy changes from 78.4881% to 78.2027%, with no-D false calls increasing from 8/32 to 12/32.

## Timing

Timing is reported from three strictly serial one-worker runs of every configuration on each complete 20,000-record simulation. Configuration order rotates by repetition; medians are reported after initialized WASM and reference setup.

| Dataset | Version/profile | Median CPU s | CPU reads/s | Median wall s | Wall reads/s | CPU change vs 0.38.1 |
|---|---|---:|---:|---:|---:|---:|
| IgM-like | 0.38.1 R-optimized | 29.874 | 669.5 | 29.454 | 679.0 | — |
| IgM-like | 0.38.3 R-optimized | **29.440** | **679.4** | **29.011** | **689.4** | **−1.45%** |
| IgM-like | 0.38.3 Sensitive-D | 29.804 | 671.1 | 29.223 | 684.4 | −0.23% |
| IgG-like | 0.38.1 R-optimized | **52.449** | **381.3** | **51.840** | **385.8** | — |
| IgG-like | 0.38.3 R-optimized | 52.869 | 378.3 | 52.299 | 382.4 | +0.80% |
| IgG-like | 0.38.3 Sensitive-D | 52.810 | 378.7 | 52.149 | 383.5 | +0.69% |

There is no meaningful slowdown in this repeated comparison. The IgM-like medians are slightly faster than 0.38.1; the IgG-like medians are 0.7–0.8% slower. R-optimized and Sensitive-D themselves are effectively speed-equivalent because they share candidate search and alignment work and differ only in the final integer decision cost. One IgM-like Sensitive-D repetition was a 32.592-second outlier; the other two were 28.727 and 29.804 seconds, so the reported median is more representative than that single run.

## Limitations

- Both simulations use one supplied simulator/reference regime. They are useful paired calibration datasets, not external biological validation.
- The fixed-10 threshold intentionally reduces no-D specificity. Users prioritizing overall accuracy or specificity should use R-optimized; users prioritizing visibility of weak but nonzero D evidence can select Sensitive-D.
- A four-base exact tract can occur by chance. It is only an admission gate: the complete D alignment and V–D–J partition must still satisfy the joint score comparison. Very short or heavily mutated D assignments remain intrinsically uncertain.
- Sensitive-D changes actual selected partitions on some reads, so small V/J endpoint movements are expected even though rejected hypotheses cannot alter boundaries.
