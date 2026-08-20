# Experimental AER-R robustness audit — Swig 0.37.6

## Question and scope

AER-R is a separately selectable experimental caller (`aer_robust`). It tests whether narrowly gated joint V/D/J boundary and candidate fallbacks can recover failures hidden by independent local alignments or fixed top-N cutoffs without slowing the ordinary path. Ordinary AER, RIAT-MP, and standard SwiftIG were not replaced.

The supplied KIMDB regression contains a retained 17-nt `IGHD5-32*01_S0263` tract that independent V alignment had consumed through a junction-facing gap. Both the WebAssembly regression test and AER-R recover a D call for that sequence under the IgBLAST-balanced profile.

## Development simulator

`tests/vdj-simulator.mjs` is development/test code and is not linked into the browser or CLI runtime. It samples from the bundled KIMDB 1.1 rhesus IGH references with:

- skewed allele usage;
- V 3′, D 5′/3′, and J 5′ deletion;
- P and N addition and occasional tandem D;
- hotspot-weighted SHM with transition bias;
- short insertions/deletions, sequencing substitutions, and `N` calls;
- partial 5′/3′ reads, non-reference flanks, and reverse complements; and
- base-level provenance, exact truth spans, and retained-tract-compatible D labels.

The release audit used seed `913771`, 5,000 records, 774 V / 52 D / 14 J alleles, double-D rate 0.025, per-base indel rate 0.0015, sequencing-substitution rate 0.0008, `N` rate 0.00025, and reverse-complement rate 0.08. The realized cohort contained 140 tandem-D records, 417 reverse complements, 2,094 indel-bearing records, and a mean 19.958 simulated SHM substitutions. This deliberately stresses the fallbacks; it is not a fitted generative model of every macaque repertoire and is not held-out validation.

## Single-D assignment results

Both strategies used the same truth-optimized calling profile and exact AIRR evaluation. A retained-tract-compatible D call accepts another allele only when its retained truth tract cannot distinguish it from the sampled allele.

| Metric | AER | AER-R |
|---|---:|---:|
| Exact V call | 95.46% | 95.72% |
| Exact J call | 96.84% | 97.28% |
| D detection, truth retained ≥6 nt | 96.86% | 98.42% |
| Exact D call, truth retained ≥6 nt | 84.86% | 86.97% |
| Retained-tract-compatible D | 88.77% | 89.91% |
| Strong-D detection, truth retained ≥10 nt | 98.17% | 98.88% |
| Strong-D retained-tract compatible | 91.57% | 92.18% |
| Mean fraction of truth D span overlapped | 90.69% | 92.57% |
| D call on a non-tandem truth retaining <6 D nt | 2.36% | 4.32% |

Case classification found 39 strong-D and 29 shorter-D retained-tract recoveries, 16 D regressions, one J regression, and 306 unresolved strong-D cases. Many unresolved labels are intrinsically weak partial-D distinctions; the exact provenance and emitted AIRR rows are written by `tests/benchmark-aer-robustness.mjs --out PREFIX` for inspection.

The increased call rate on records retaining fewer than six true D bases is the principal observed specificity tradeoff. This is why AER-R is experimental and opt-in rather than the default.

## Candidate oracle

`tests/build-aer-robust-oracle-wasm.sh` compiles a validation-only binary that disables production candidate pruning, aligns every supplied V/D/J allele, and disables the nominal-window fast path. It is intentionally too slow for runtime use.

On 30 prioritized difficult records from a 1,000-record replicate, production AER-R had no D recovery available only to the exhaustive oracle. It agreed with the oracle's evaluated V/J/D correctness classification on 29/30 records. The sole disagreement was a V allele recovered by exhaustive search; both production and oracle still failed the same D/J truth classification. The oracle required 28.14 seconds for 30 records, confirming why it is a diagnostic rather than a shipping method.

## Timing and output invariance

All timings are one-worker WebAssembly wall time on the same host. They should be read as local regression checks, not portable throughput promises.

For two repetitions of 5,000 identical clean human-IGH records:

| Strategy | Run 1 | Run 2 | Mean |
|---|---:|---:|---:|
| AER | 4.759 s | 4.811 s | 4.785 s |
| AER-R | 4.795 s | 4.787 s | 4.791 s |

The mean difference was 0.13%, within run-to-run timing noise. Every output had the same SHA-256, `62bd4dcee34401f46fb58bea75beee5341db74da3928ef3dcf239f5415399614`. Ordinary AER also produced that exact hash with the pinned 0.37.5 WebAssembly binary, demonstrating that the new strategy did not alter the old one.

On the difficult 5,000-record simulation, AER took 8.920 s and AER-R took 9.167 s: AER-R retained 97.30% of AER throughput, a 2.7% cost when its fallback work is frequently triggered. Anchor projection is lazy, the maximum D length is cached once per worker, and the complete-J retry is gated to a strongly V-supported orientation with no valid partition.

## Optional double-D screen

With double-D screening explicitly enabled:

| Metric | AER screen | AER-R screen |
|---|---:|---:|
| Tandem-D detection | 42.14% | 55.71% |
| Ordered sequence-equivalent D1/D2 | 36.43% | 41.43% |
| False positive among simulated single-D records | 0.103% | 0.267% |
| Time, 5,000 records | 8.573 s | 9.045 s |

AER-R's shorter seed must extend to the ordinary configured seed length before pair scoring. It improves mutation tolerance, but its measured sensitivity/specificity tradeoff is another reason to keep the method explicit.

## Reproduction

```bash
npm run build:wasm
npm run benchmark:aer-robust -- --count 5000 --seed 913771 --screen-double-d --out tmp/aer-r
bash tests/build-aer-robust-oracle-wasm.sh
node tests/benchmark-aer-robustness.mjs --count 1000 --seed 913771 \
  --oracle-wasm .build/aer-robust-oracle.wasm --oracle-count 30
```

The oracle build requires `WASI_SDK`. Generated simulation FASTA/truth/case reports are development artifacts and are not included in the runtime package.
