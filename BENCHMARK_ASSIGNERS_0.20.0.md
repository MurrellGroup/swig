# Assignment-strategy benchmark — Swig 0.20.0

## Scope

This report compares three selectable V-assignment strategies while holding the
truth-optimized D and J implementation and parameters fixed. It also records a
separate structured-J experiment. The supplied low-SHM and IgG simulations were
used during development, so these are tuning-set measurements rather than
held-out validation.

- **AER — Adaptive Exact Refinement:** uses the complete 9-mer V index and exact
  affine alignment. A strong-seed query normally refines three candidates. When
  the next candidates remain within the larger of 8 weighted votes or 5% of the
  leading vote count, the exact refinement set grows, up to 16 alleles.
- **RIAT-MP — Root-Indexed Allele Tree, Multipath:** clusters close V alleles,
  indexes only representative roots, aligns at most three selected roots, and
  propagates sparse reference substitutions through each allele tree. If the
  provisional winner contains an indel, it evaluates at most two root traceback
  geometries within four raw-score units, with a 1,024-state cap. Descendant V
  alleles are not independently aligned.
- **Standard SwiftIG:** the existing complete V index with three exact
  strong-seed refinement alignments and its prior weak-seed/seedless safety
  pools.

All three use the same truth-oriented D and J settings: D `+2/-3/-13/-1`, two
reported hits, and a six-base support floor; J `+2/-3/-17/-2`, two reported
hits, and a ten-base aligned-length floor. Here the four scores are
match/mismatch/gap-open/gap-extension.

## Scoring and sample

The benchmark retains records for which `FNV1a(sequence_id) % 10 == 0`. This
selects 9,994 records independently from each 100,000-record simulation. No
outcome-dependent subsampling was used.

- **First** gives one point only when the first reported allele equals the first
  reference allele. Two empty calls also agree.
- **Fair** gives `|prediction ∩ reference| / |prediction|`. Thus extra reported
  alleles dilute the score instead of creating free correct calls. Joint fair
  VDJ is the product of the three per-segment fractions for each record.

Every table cell below is **first / fair, percent**.

## Agreement with simulated truth

### Low-SHM simulation

| Method | V | D | J | Joint VDJ |
| --- | ---: | ---: | ---: | ---: |
| AER | 92.56 / 92.69 | 77.84 / 77.99 | 97.60 / 97.37 | 70.47 / 70.47 |
| RIAT-MP | 92.39 / 92.52 | 77.84 / 77.99 | 97.61 / 97.38 | 70.31 / 70.30 |
| Standard SwiftIG | 90.27 / 90.37 | 77.84 / 77.99 | 97.59 / 97.36 | 68.65 / 68.68 |
| IgBLAST reference calls | 92.39 / 92.52 | 76.50 / 76.76 | 97.87 / 97.61 | 69.37 / 69.44 |

### IgG / higher-SHM simulation

| Method | V | D | J | Joint VDJ |
| --- | ---: | ---: | ---: | ---: |
| AER | 81.37 / 81.29 | 68.47 / 68.92 | 93.15 / 92.97 | 52.77 / 53.05 |
| RIAT-MP | 80.43 / 80.47 | 68.51 / 68.97 | 93.15 / 92.97 | 52.18 / 52.63 |
| Standard SwiftIG | 77.01 / 76.92 | 68.47 / 68.92 | 93.14 / 92.96 | 50.20 / 50.33 |
| IgBLAST reference calls | 81.90 / 81.86 | 67.30 / 67.72 | 92.75 / 92.57 | 52.07 / 52.49 |

AER is the default because it gives the best truth score of the tested Swig
strategies in both mutation regimes while retaining exact allele alignment.
RIAT-MP is a faster experimental alternative with substantially better V
accuracy than the fixed-depth path, but it remains approximate.

## Agreement with IgBLAST calls

### Low-SHM simulation

| Method | V | D | J | Joint VDJ |
| --- | ---: | ---: | ---: | ---: |
| AER | 99.02 / 99.34 | 94.10 / 94.80 | 99.48 / 99.63 | 92.71 / 93.83 |
| RIAT-MP | 99.13 / 99.23 | 94.12 / 94.82 | 99.49 / 99.64 | 92.87 / 93.78 |
| Standard SwiftIG | 96.40 / 96.73 | 94.11 / 94.81 | 99.46 / 99.61 | 90.25 / 91.36 |

### IgG / higher-SHM simulation

| Method | V | D | J | Joint VDJ |
| --- | ---: | ---: | ---: | ---: |
| AER | 96.86 / 98.00 | 89.47 / 90.61 | 98.33 / 99.00 | 85.61 / 88.18 |
| RIAT-MP | 96.69 / 97.20 | 89.49 / 90.65 | 98.34 / 99.01 | 85.47 / 87.52 |
| Standard SwiftIG | 91.40 / 92.80 | 89.47 / 90.61 | 98.32 / 98.99 | 80.68 / 83.43 |

These agreement values describe the truth-optimized calling profile. Choosing
one of Swig's IgBLAST-oriented calling profiles changes D/J parameters and is a
separate control.

## Throughput

The harness uses one Node/WASI instance, one thread, both strands, the complete
supplied human IGH references, AIRR formatting, and batches of 500. Browser
throughput depends on worker count, hardware, browser, selected reference set,
input length, and storage. Rates are therefore useful only as a matched relative
comparison here.

| Method | Low-SHM reads/s | IgG reads/s | Relative to standard, low / IgG |
| --- | ---: | ---: | ---: |
| AER | 332 | 227 | 0.96× / 0.93× |
| RIAT-MP | 504 | 348 | 1.45× / 1.43× |
| Standard SwiftIG | 347 | 244 | 1.00× / 1.00× |

This workload is not the same as SwiftIG's published native benchmark and must
not be interpreted as a WebAssembly/native ratio.

## Structured J experiment

The supplied J reference has 12 alleles in seven close-allele clusters. A
separate prototype aligned all seven roots and propagated descendant scores,
with at most four near-optimal root paths within eight raw-score units and a
2,048-state cap. With the calibrated J scoring, its J truth score was:

| J implementation | Low-SHM truth | IgG truth | Low-SHM IgBLAST agreement | IgG IgBLAST agreement |
| --- | ---: | ---: | ---: | ---: |
| Exact J retained in 0.20.0 RIAT-MP | 97.61 / 97.38 | 93.15 / 92.97 | 99.49 / 99.64 | 98.34 / 99.01 |
| Root/multipath J prototype | 97.59 / 97.36 | 93.06 / 92.91 | 99.43 / 99.62 | 98.22 / 98.94 |

The prototype's combined elapsed time was about 5% lower in this run, but the
direction reversed between the two datasets (slower on low-SHM, faster on IgG).
J sets are small, J alignments are short and variably 5'-truncated, and endpoint
placement determines which allele-specific positions are scored. The measured
speed change was therefore not sufficiently stable to justify the small loss or
the additional risk across the many bundled species/loci. Swig 0.20.0 keeps
exact J alignment for AER, RIAT-MP, and standard assignment.

## Reproducibility

Input and core SHA-256 values:

| File | SHA-256 |
| --- | --- |
| `sim_100k(1).tsv` | `911eaa3027d249e065e264ff3c826ac7619ae7d562efe5eff118c510ded27e8c` |
| `sim_100k_igg(1).tsv` | `b2842c7de8adc748a90f748e2a49aaad175a352b91869be1ea58425e8a75d4af` |
| `KI+1KGP-IGHV-SHORT.fasta` | `2df3dbb5d8c936c55bb28cc9d53d0f188a760a8ec442e2615f3ec77ebd796abe` |
| `KI+1KGP-IGHD-SHORT.fasta` | `70f9ac7b0a4c86e68e018b14c0eaf4bc32db802b877a6173503519b329df70fc` |
| `KI+1KGP-IGHJ-SHORT.fasta` | `4bc809d2d0081adafb8c21fd478f989a133e09b97079dcb5092ff2944f904ca0` |
| `public/swiftig.wasm` | `d6566c0e7683c462f66cd22c124c0b952eb4b28a84910cf707f28b6ff98d3b7c` |

The raw JSON summaries are included under `benchmarks/0.20.0/`.
