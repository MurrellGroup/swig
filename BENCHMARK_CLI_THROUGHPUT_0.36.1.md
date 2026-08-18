# Swig CLI throughput audit for 0.36.1

## Scope

The benchmark compares the packaged 0.36.0 CLI with the optimized 0.36.1 source bundle. It uses 20,000 identical full-length human IGH rearrangements, the complete bundled human IGH V/D/J set, AER assignment, and an eight-worker Node 24 control runtime on an allocation exposing nine logical CPUs. Each output is 46.2 MB without CDR/FWR metadata and 64.4 MB with a prepared annotated reference.

Node is used as the controlled before/after runtime because this sandbox prevents Bun from executing JIT code. The standalone release uses Bun 1.3.14. Bun-specific gains from its flat-message fast path and the modern AVX2 executable therefore are not included in the measured percentages below.

## Results

| Case | 0.36.0 | 0.36.1 | Wall-time reduction | Output comparison |
|---|---:|---:|---:|---|
| Assignment only, explicit 1,000-record batches; mean of two runs | 10.74 s | 9.68 s | 9.9% | byte-identical |
| Prepared CDR/FWR annotations, explicit 1,000-record batches; mean of two runs | 10.54 s | 9.54 s | 9.5% | byte-identical |
| Gzip input, explicit 1,000-record batches | 10.11 s | 9.12 s | 9.9% | byte-identical |
| Assignment only, each version's default batch size | 13.00 s | 8.94 s | 31.3% | byte-identical |
| Prepared CDR/FWR annotations, each version's default batch size | 12.25 s | 10.17 s | 16.9% | byte-identical |

The larger default improvement comes from replacing the fixed 2,000-record batch with worker-aware sizing. On this finite 20,000-read job, 2,000-record batches produce only ten tasks for eight workers and leave a long final tail; 500-record batches produce forty tasks. A one-worker 5,000-read control was 2.3% slower with 500 than with 2,000 records, which is why automatic sizing retains 2,000 for one or two workers.

Writing the 46.2 MB output to `/dev/null` took 8.98 s versus 8.94 s for the workspace file, and parsing/canonicalizing the 7.8 MB FASTA alone took 0.19 s. Local input and disk output therefore are not the principal bottlenecks at eight workers; assignment remains compute-dominant. Larger stream buffers mainly protect high-worker-count and slower-filesystem runs from avoidable callback/backpressure stalls. Per-worker temporary output would add writes and a collation pass and was not adopted.

At 32 requested workers on the same nine-CPU allocation, 20,000 reads took 13.19 s with 2,000-record batches and 10.22 s with 500-record batches. This is an oversubscription/granularity check, not a 32-core scaling measurement.

## Build targets

The previous Linux x64 release used `bun-linux-x64-baseline`. Bun defines that target for pre-2013/Nehalem-compatible CPUs. The release workflow now additionally emits `bun-linux-x64-modern`, which Bun defines for Haswell/AVX2 and describes as faster. Baseline remains available under the existing `swig-cli-linux-x64` name.

`--minify --bytecode` is enabled for standalone builds. Bun documents bytecode compilation as moving JavaScript parsing to build time, so it improves startup rather than the steady-state SwiftIG kernel. The SwiftIG WebAssembly itself was already built with `-O3`, LTO, SIMD128, bulk memory, `wasm-opt -O4 --converge`, and stripped release assertions; no conservative WASM optimization flag was found to remove without changing compatibility or numerical behavior.

Official Bun references:

- <https://bun.sh/docs/bundler/executables>
- <https://bun.sh/docs/runtime/workers>
