# KIMDB D-boundary regression and hot-path timing — Swig 0.37.5

## Reproduction

The regression fixture is `tests/fixtures/kimdb-perfect-d-regression.fasta`, annotated against the bundled `kimdb-1.1/Macaca_mulatta/IGH` V/D/J FASTAs. The simulated header identifies `IGHD5-32*01_S0263`, the 17-nt D sequence `ATACAGTGGGTACAGTT`, and one-based gold span 421–437.

Swig 0.37.4 selected `IGHV2-69*01_S2931` through query base 433 and `IGHJ4-3*01_S9191` from base 439. Its V alignment ended as:

```text
query     ...TGTGCGATACAGTGGGTAC
V allele  ...TGTGC---AC---GGGTAC
```

The corresponding V CIGAR was `126S293M3I2M3I6M122S`. Two inexpensive three-base insertions plus terminal matches made this a valid independent V local-alignment optimum, but it consumed query positions belonging to D. The ordinary D search was subsequently restricted to the five remaining bases between V and J. The agreement profile could form only a five-base provisional D hit, which the balanced rule correctly removed; the long D evidence had never entered the D search.

## Correction

The rescue runs only when the selected ordinary D has fewer than 10 consecutive exact matches and the V 3′ or J 5′ edge contains an affine gap within 24 query bases. A cheap exact 10-mer test against locus-compatible D references must pass before the expanded D dynamic program is invoked. For every refined D candidate, Swig:

1. converts its alignment back to complete-query coordinates;
2. clips any overlapping V suffix and J prefix;
3. recomputes the retained V/J affine scores; and
4. compares disjoint V+D+J against the original V+(weak D)+J score.

Only a strict score improvement is accepted. A failed rescue leaves the original call and endpoints unchanged.

## Regression result

The same result was obtained with AER, RIAT-MP, and Standard SwiftIG under the agreement settings:

| Field | Swig 0.37.4 balanced output | Swig 0.37.5 |
|---|---:|---:|
| `v_sequence_end` | 433 | 419 |
| `d_sequence_start` | empty | 420 |
| `d_sequence_end` | empty | 437 |
| `j_sequence_start` | 439 | 439 |
| D aligned query/reference | empty | `GATACAGTGGGTACAGTT` |
| `d_call` | empty | `IGHD5-27*01,IGHD5-32*01,IGHD5-32*01_S0263` |

The recovered match starts one base before the simulated D boundary because that junction base also matches all three references. It contains the complete gold 17-mer exactly. The three names are therefore genuine co-optimal calls over the observed tract, not uncertainty introduced by the rescue.

The automated test initializes each assignment strategy before building its strategy-specific V index, asserts recovery of the complete gold tract, verifies identical D query/germline alignment, and requires non-overlapping V/D coordinates.

## Hot-path timing

Old and fixed release WASM binaries were instantiated together against the same bundled human IGH references. After warm-up, identical 1,000-record ordinary VDJ batches were sent to the two runtimes in alternating order; AIRR outputs were asserted byte-for-byte equal. The shared host occasionally paused a process for seconds, so runs exceeding three times the minimum for that binary were classified as scheduler stalls before calculating the median.

| Strategy | 0.37.4 median | 0.37.5 median | Difference |
|---|---:|---:|---:|
| AER | 805 ms | 800 ms | −0.6% |
| RIAT-MP | 775 ms | 747 ms | −3.6% |

These differences are ordinary timing noise; no slowdown was detected. The structural reason is stronger than the wall-clock result: the successful-D path already computed the retained exact-run length, and 0.37.5 merely carries that scalar forward and evaluates a branch. It performs no extra seed scan or alignment unless D evidence is weak and a junction-facing gap is present.
