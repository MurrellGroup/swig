# SwiftIG `swig_annotate` kernel audit for 0.37.0

## Scope and invariants

This audit profiles the work below the WebAssembly `swig_annotate` export and optimizes the AER and RIAT-MP production paths. It does **not** change scoring, candidate limits, band selection, tie breaking, acceptance thresholds, AIRR fields, or assignment/annotation semantics. The pre-existing scalar candidate search, scalar affine alignment, and ostream AIRR writer remain in the source as reference implementations. RIAT multipath recurrence, priority, and state limits are unchanged while its partial-trace storage is allocation-reduced. Diagnostic WASM exports select reference or production kernels for equivalence tests; ordinary web and CLI callers use the optimized path.

Two controls are used:

1. A native C++ profiling harness annotates 10,000 deterministic full-length human-IGH queries containing exact, substitution, indel, and `N` variants against all 543 supplied human IGH V/D/J alleles. It makes function-level sampling practical while exercising the same `swig_annotate` entry point.
2. The actual release WASM runs in one Node worker through `@bjorn3/browser_wasi_shim`, with one bounded 10,000-record batch and the complete human IGH reference. This includes query parsing, annotation, AIRR serialization, WASM/JS copying, and worker messaging. It is the end-to-end number reported below.

The timing host was the same restricted nine-logical-CPU allocation used for the 0.36.x CLI audits. Single-worker results avoid scheduling and oversubscription effects. Absolute throughput will vary by browser/runtime and CPU; the before/after ratio is the useful result.

## Where the reference implementation spent time

The original no-LTO native sampling profiles showed that `swig_annotate` was already compute-bound, but much of that compute was avoidable machinery around the intended alignment calculation.

| Reference hot-path group | AER sampled self time | RIAT-MP sampled self time | Cause |
|---|---:|---:|---|
| Scalar affine local alignment | 47.1% | 51.2% | One DP and traceback allocation per candidate; candidates handled serially |
| Seed candidate construction | about 31% | about 24% | Node-based hash lookups for seed hits and vote bins, rehashing, sorting every diagonal bin, then repeatedly rescanning bins by gene |
| RIAT multipath traceback | — | 8.1% | Alternative partial traces copied growing strings/vectors at every branch |
| String/vector allocation and destruction | about 18.5% | about 8.2% | DP scratch buffers, temporary hit containers, and a roughly 100-string AIRR row |

FASTA parsing, CDR/FWR projection, E-value calculation, translation, and AIRR field logic were not individually important. Optimizing them first would not have materially changed throughput.

## Production-kernel changes

### Exact four-candidate affine SIMD

Four independent candidate alignments occupy the four 32-bit lanes of one 128-bit vector. Match/mismatch, affine insertion/deletion recurrences, local-zero resets, direction ties, best-cell ties, band boundaries, and traceback flags use the same comparisons and ordering as the scalar reference. This maps directly to WebAssembly SIMD128 and does not require architecture-specific intrinsics.

Wide safety pools first receive the same exact affine calculation in a score-only SIMD kernel. Candidates are sorted by that exact score; complete DP plus traceback is then performed through the score group that can still enter the requested top-N result. Invalid high-scoring traces cause evaluation to continue, and score ties at the cutoff are retained. This removes trace-matrix writes and alignment-string construction for candidates that provably cannot be reported without introducing an approximate prefilter.

An `AlignmentWorkspace` owns reusable score rows, insertion rows, trace matrices, and CIGAR-operation storage for the lifetime of an annotation batch. Banded and unbanded scalar fallback loops are separately specialized so the inner loop does not repeatedly branch on band mode.

### Exact seed-vote reduction

Canonical DNA k-mers are already integers in `[0, 4^k)`. All production indexes use `k <= 9`, so AER/RIAT use a direct-address offset table with contiguous CSR seed hits instead of hashing each query k-mer. The original unordered-map index remains populated for the reference path.

Candidate votes use a reusable, generation-stamped dense `(gene, diagonal-bin)` table on ordinary reads, with an allocation-free flat-hash fallback for exceptional dimensions. A two-pass per-gene reduction computes exactly the same best bin and strong-bin span as the reference implementation while sorting one candidate per touched gene rather than every occupied bin. Per-gene scratch storage is also generation-stamped, so reset touches only genes observed in the current query.

Forward/reverse orientation screening already computes complete V and J candidate rankings. The selected orientation now reuses those rankings instead of repeating the same seed search.

### RIAT-MP-specific work

RIAT root alignments use the same four-lane full-DP kernel. Allele-tree child links are constructed once with stable first-child/next-sibling arrays rather than rebuilding nested adjacency vectors per query. Near-optimal multipath traceback retains the same priority and state limit but stores partial traces as constant-size parent-linked arena nodes; strings and CIGAR operations are materialized only for terminal paths.

### AIRR serialization

The production writer appends directly to the result string. Integers use `to_chars`; fixed/scientific floating-point formatting preserves the reference formatting; optional and hit fields are emitted without first constructing a vector of roughly 100 temporary strings. Alternative-evidence fields are appended directly as well. The ostream implementation remains callable as the reference writer.

## Final profile

The optimized no-LTO function-level sample is now dominated by the intended exact calculations:

| Final hot-path group | AER | RIAT-MP |
|---|---:|---:|
| Four-lane full DP and traceback | 58.3% | 52.3% |
| Four-lane exact score screening | 14.2% | 11.5% |
| RIAT near-optimal multipath DP/traceback | — | 9.6% |
| Seed scanning, voting, reduction, and candidate sort | 22.7% | 12.2% |
| RIAT fixed-path tree projection/materialization | — | about 6.5% |
| Direct AIRR row writer | 0.4% | 0.1% |

Thus about 72.5% of sampled AER time and 73.4% of sampled RIAT-MP time is now affine DP/traceback. For RIAT-MP, another roughly 6.5% is the algorithm's actual sparse allele-tree projection. Parsing, annotation projection, statistics, and serialization collectively sit in the residual tail. Further large gains would require a substantially different exact Smith-Waterman implementation rather than removal of obvious orchestration, allocation, or hashing overhead.

## End-to-end release-WASM result

| Strategy, 10,000 reads, one worker | 0.36.1 reference WASM | 0.37.0 production WASM | Throughput gain | Peak RSS |
|---|---:|---:|---:|---:|
| AER | 27.608 s; 362 reads/s | median 12.576 s; about 795 reads/s | **2.20x** | 268.1 to about 233.5 MiB |
| RIAT-MP | 18.733 s; 534 reads/s | 10.580 s; 945 reads/s | **1.77x** | 266.5 to 229.4 MiB |

The native diverse-query control produced the same original 64-bit FNV-1a output hashes after every compute and writer optimization:

- AER: 32,171,680 bytes, `7416024721985682163`
- RIAT-MP: 32,182,426 bytes, `1953680844651455677`

The permanent WASM regression test separately compares (a) reference compute plus reference writer, (b) optimized compute plus reference writer, and (c) optimized compute plus optimized writer for both AER and RIAT-MP on exact, mutated, ambiguous, insertion, deletion, and reverse-complement records. It requires byte-for-byte equality, not merely equal calls or scores.

## Interpretation

This is a same-computation optimization. It does not claim that AER or RIAT-MP is biologically superior, change their previously documented calibration scope, or make runtime independent of read length, germline size, mutation/indel burden, or seed ambiguity. Hard seedless reads still deliberately pay for broad exact search. Multiple CLI workers remain independent WASM instances; the per-worker kernel speedup composes with, but does not replace, the separate worker-scaling and I/O improvements documented for 0.36.x.
