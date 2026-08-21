# AER-R R-optimized calibration — Swig 0.38.0

## Scope and frozen profile

R-optimized is an alternative AER-R calling profile, not a replacement for legacy AER-R and not a new assignment strategy. It is selected as `r_optimized` only with `aer_robust`; every prior profile remains available and unchanged.

The frozen settings are:

| Component | R-optimized setting | Runtime consequence |
|---|---|---|
| V affine score | `+2/−4/−13/−1` | Same candidate alignments and recurrence; endpoint ranking changes. |
| D affine score | `+2/−4/−13/−1` | Same SIMD kernel. |
| D evidence floor | 5-nt exact run; 2 leading candidates | Existing AER-R distributed-evidence route remains available. |
| D presence | subtract 10 raw-score points, relaxed to 8 for score ≥18 or ≥3 distinct exact template classes | Cached scan of existing D hits; no new alignment or DP. |
| D reported ambiguity | include same-locus, same-query-span D hits within 1 raw score point | Linear scan of the already-ranked small hit set; boundary is unchanged. |
| J affine score | `+2/−4/−17/−2`; 2 candidates | Same candidate alignments and recurrence. |

The stronger V gap-open cost addresses a specific failure: a local V optimum could cheaply open a junction-facing gap, resume matching, and consume N/D bases. The D-presence cost addresses a separate multiple-comparison effect: after maximizing over D alleles, positions, and lengths, any positive short D score previously received a free advantage over the no-D partition. The ten-point base cost is reduced to eight for a D alignment scoring at least 18 or an exact ungapped tract present in at least three distinct locus-matched D template sequences. These conditions preserved the prior short-D and zero-D/J regression gates while improving validation proper loss; they are not claimed to be a universal germline-recombination prior.

## Data and partitions

The supplied truth table has 20,000 simulated macaque IGH records. Assignment used the deduplicated union of the bundled KIMDB 1.1 rhesus and cynomolgus V/D/J references: 1,276 named records. The split is deterministic from FNV-1a(`sequence_id`) modulo 10:

| Partition | Residues | Records | Single-D score records | Use |
|---|---|---:|---:|---|
| Training | 0–5 | 12,001 | 11,993 | Failure classification and parameter search. |
| Validation | 6–7 | 3,991 | 3,986 | Profile/threshold selection. |
| Test | 8–9 | 4,008 | 4,002 | Opened after freezing parameters. |

Nineteen comma-separated tandem-D truths (8/5/6 by split) are excluded only from the primary single-D call score. They remain annotated and counted in V/J and boundary summaries.

The user-supplied 760-record V-stop cohort and prior 43-record J-miss/93-record D-miss cohorts are derived from the same 20,000 records. Because those lists were inspected as targeted regressions, a second clean test report excludes every record in those cohorts plus the two individually investigated D cases (`1359`, `7769`). That leaves 3,844 test records. This exclusion is the appropriate endpoint generalization check; the unfiltered test result is also reported for completeness.

Input SHA-256 values:

| Input | SHA-256 |
|---|---|
| 20,000-record AIRR truth | `1a11279992c8e51698a30429f75da45c926d3f15d770fce70739fea8e6f10ddd` |
| 760 V-stop cohort | `041c5092b82c69051dd940df616710596b9e105f8b07fe14b59b9636f47562a3` |
| 43 prior J misses | `ea85a31ccf8abf39f1ed8987d246618889881a2c75aac1dd1025afc4cfb20f55` |
| 93 prior D misses | `d2094aa6e7b378926ef0c1f50e19f5d806db542da426b72bad5d4007618f5dde` |
| V reference union | `5081474c5f67e141a5871bff85bf4298e27d51f16dae7d19cd1ff0a01eaad7af` |
| D reference union | `0ee3cc470002350ef2ed416427935e7d838d17cce83556d31725f0c84cbd1145` |
| J reference union | `cc4cdcc09f4248e26a3ad504cf64e9a2a35704f843c09085113fb91e1253e0fa` |

## Proper call score

Named references with identical complete nucleotide sequences are collapsed into one equivalence class. A comma-separated reported call is interpreted as a uniform distribution over its unique equivalence classes. Truth is one-hot, including an explicit no-D class. The normalized multiclass Brier loss is

\[
L = \frac{1}{2}\sum_k (p_k-y_k)^2,
\]

so a correct singleton scores 0, a wrong singleton scores 1, and adding weak alternatives is penalized. The uniform interpretation is an evaluation convention; comma-separated SwiftIG calls are score uncertainty, not posterior probabilities. Endpoint exactness, within-one accuracy, MAE, signed quantiles, and absolute-error tails are scored separately so ambiguity cannot hide trimming error.

An all-reference candidate oracle was also run on 80 prioritized development errors. At roughly 0.42 reads/s it rescued only 2/20 sampled V-call errors, no sampled J-call error, and one sampled D-call error. Most remaining failures were therefore joint scoring or genuine sequence ambiguity rather than candidate pruning; no runtime exhaustive fallback was added.

## Clean test result

The clean 3,844-record test excludes all previously inspected target cases.

| Call metric | Legacy AER-R | R-optimized | Change |
|---|---:|---:|---:|
| V mean Brier | 0.031781 | 0.031391 | −1.23% |
| V truth in reported set | 98.647% | 98.777% | +0.130 pp |
| D mean Brier | 0.179954 | 0.160987 | −10.54% |
| D truth in reported set | 85.070% | 86.868% | +1.798 pp |
| J mean Brier | 0.051791 | 0.049211 | −4.98% |
| J truth in reported set | 98.491% | 98.881% | +0.390 pp |

| Endpoint, all records | Legacy exact | R exact | Legacy within 1 | R within 1 | Legacy MAE | R MAE |
|---|---:|---:|---:|---:|---:|---:|
| V sequence end | 69.875% | 72.399% | 88.840% | 91.571% | 0.539 | 0.427 |
| J sequence start | 68.314% | 70.760% | 86.941% | 89.958% | 0.595 | 0.468 |

D-presence sensitivity changes from 98.769% to 99.344% (+0.574 pp), while no-D specificity changes from 50.276% to 62.431% (+12.155 pp). Thus the conditioned decision improves both sides of the D/no-D decision on the clean test rather than trading one against the other.

For records whose truth allele is reported by both profiles, V-end MAE improves and J-start MAE improves. Paired D endpoints are essentially unchanged: the very small D-start movement is dominated by difficult newly truth-compatible calls rather than a general shift in D trimming.

## Complete test result

On all 4,008 test records, including targeted cases, V/D/J Brier changes are `0.033891 → 0.032560`, `0.180776 → 0.162226`, and `0.051065 → 0.048653`. V-end exactness changes `67.465% → 72.991%`; J-start exactness changes `68.388% → 71.493%`. D-presence sensitivity changes `98.530% → 99.134%`, and no-D specificity changes `53.125% → 64.063%`.

## Targeted regression cohorts

| Cohort | Legacy AER-R | R-optimized | Interpretation |
|---|---|---|---|
| 760 dramatic V-stop cases | 4.21% exact; 10.39% within 1; MAE 4.518; p95 9 nt | 75.66% exact; 93.09% within 1; MAE 0.376; p95 2 nt | Junction-facing V gap failure is corrected without forcing signed error to zero. |
| 43 prior J misses | J truth-in-set 100%; no-D specificity 100% | Identical | Existing zero-D complete-J rescue is retained. |
| Strong 17-nt perfect D case `1359` | Legacy AER-R already rescued after 0.37.5 | Truth allele retained in a three-class D uncertainty set | Ten-point D cost does not remove strong exact D evidence. |
| Distributed 14/17 D case `7769` | AER-R calls `IGHD6-39*01` | Same D retained | Distributed-D rescue survives the new profile. |
| 93 prior D misses | Truth in set 18.28%; Brier 0.800 | Truth in set 22.58%; Brier 0.678 | Deliberately adverse short-D enrichment: 69 have ≤5 retained bases and 17 have 6–9. Seven legacy truth-compatible calls change away from truth while eleven previously missed calls become truth-compatible: a net gain of four, not a claim that every intrinsically ambiguous short tract is preserved. |

The user explicitly requested individually best boundaries rather than a bias-corrected mean. Selection therefore used exact/within-one/MAE/tail loss. R-optimized retains a small positive V-end signed mean (`+0.179 nt` on clean test), rather than shifting all endpoints merely to force mean error to zero.

## Runtime

Four alternating one-worker runs per profile used one initialized WASM instance and the 3,991-record validation partition, so startup could not dominate the comparison. Median elapsed time was 7.133 s for legacy AER-R and 6.369 s for R-optimized (−10.7%). Individual runs were noisy (`6.108–12.467 s` legacy; `6.207–7.429 s` R-optimized), so the apparent speedup should not be overinterpreted; importantly, no slowdown was observed. The conditioned D cost reuses existing hits, stops its exact-template scan after three distinct sequences, and performs no additional alignment or dynamic programming.

## Reproduction

`tests/benchmark-r-optimized.mjs` reads gzipped/plain AIRR truth or the targeted header-annotated FASTA format, performs the deterministic split, applies sequence-equivalence-aware scoring, and optionally writes AIRR predictions/failure JSONL. For example:

```bash
node tests/benchmark-r-optimized.mjs \
  --truth sim_igsim_indelp0.003_n20000__assign_assign__sim_all.airr.tsv.gz \
  --v V.fasta --d D.fasta --j J.fasta \
  --profile r_optimized --split validation
```

The simulator, truth data, and calibration outputs are development inputs, not runtime dependencies. This is internal validation within one simulator and KIMDB macaque reference regime, not evidence of universal superiority on biological repertoires.
