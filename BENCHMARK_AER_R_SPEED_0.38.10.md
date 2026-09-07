# AER-R / R-optimized speed audit for Swig 0.38.10

## Invariant

This work changes execution only. AER-R scoring, R-optimized thresholds, candidate and fallback coverage, affine recurrences, tie-breaking, boundary adjudication, AIRR columns, and worker semantics are unchanged. The production kernel is checked against the retained scalar/reference implementations at byte level.

## Profile

A symbolized production-flags WASM profile of 5,000 real ERR4238110 reads against the complete KIMDB set assigned 43.6% of samples to full four-lane affine alignment, 28.0% to exact candidate score screening, and 16.4% to seed collection/voting/reduction. Parsing, JavaScript orchestration, and AIRR writing were individually below 1%.

## Paired real-data benchmark

Inputs were read directly from the supplied gzip FASTA files. Each row used the complete rhesus KIMDB 1.1 reference (774 V, 52 D, 14 J records), `aer_robust`, `r_optimized`, plus strand, one worker, and 1,000 records per batch. Runs used the same first 5,000 records in source order. The host was the restricted Codex Linux runtime; absolute rates are host-specific.

| Input | Swig 0.38.9 wall time | Swig 0.38.10 wall time | Throughput gain | AIRR identity |
|---|---:|---:|---:|---|
| ERR4238110, paired run 1 | 19.729 s | 14.492 s | 1.361x | byte-identical |
| ERR4238110, paired run 2 | 21.013 s | 14.678 s | 1.432x | byte-identical |
| ERR4238104 | 20.287 s | 14.566 s | 1.393x | byte-identical |

The two-run ERR4238110 midpoint improves from 20.371 to 14.585 seconds: **1.397x throughput, or 28.4% less wall time**. ERR4238104 independently shows **1.393x throughput, or 28.2% less wall time**.

A separate 20,000-read ERR4238110 run with four workers improved from 25.674 to 16.254 seconds: **1.579x throughput, or 36.7% less wall time**. Its two AIRR files were byte-identical with SHA-256 `9751f7fd2720e841478f4a59f6f06057a0fcb89c3ea195a3ef4878d6ea82724c`.

The first 10,000 ERR4238110 records were also compared end to end. The optimized AIRR and 0.38.9 AIRR both have SHA-256 `a413026d823bf001fdd4c7770c035fcb7d028359f3aec37841c29646304b9129`. The paired ERR4238104 5,000-record outputs both have SHA-256 `b30c88e49b669f6377257f2552b84cbfd8c6b8ad10d8ff1d1045c9db213ce229`.

## Interpretation

The measured gain is a single-worker kernel gain and therefore composes with the existing worker pool. Exact duplicates were only 1.3% of the first 10,000 ERR4238110 records, so the result is not a duplicate-cache artifact. After optimization, approximately 72% of sampled time remains in the intended exact affine alignment and score recurrences; further large gains would require a new exact DP implementation rather than wrapper or allocation tuning.
