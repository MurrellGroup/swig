# SwiftIG worker scaling audit for Swig 0.36.0

## Question

The V(D)J kernel has no cross-read dependency, so why is end-to-end worker scaling not perfectly linear, and can the CLI start 32 workers?

## Reproducible benchmark

`tests/benchmark-wasm-pool.mjs` constructs 50,000 identical full-length human-IGH rearrangements, uses the bundled human IGH V/D/J references (543 alleles per worker), runs the standard SwiftIG assignment path in 1,000-record batches, and transfers the complete AIRR output back to the main Node process. Timing starts after every worker has compiled/instantiated WASM and built its private reference indexes; it therefore measures steady-state annotation and worker messaging, not startup or disk writing.

The test environment reported nine available logical CPUs. Hardware topology and dedicated-versus-shared-core status were not exposed, so the 16- and 32-worker runs are oversubscription checks, not estimates for a dedicated 16/32-core host.

```bash
SWIG_BENCH_WORKERS=N SWIG_BENCH_BATCH=1000 \
  node tests/benchmark-wasm-pool.mjs
```

| Workers | Reads | Reads/s | Speedup | Parallel efficiency | Peak RSS |
|---:|---:|---:|---:|---:|---:|
| 1 | 50,000 | 399 | 1.00× | 100% | 176.0 MiB |
| 2 | 50,000 | 743 | 1.86× | 93% | 220.5 MiB |
| 4 | 50,000 | 1,353 | 3.39× | 85% | 319.4 MiB |
| 8 | 50,000 | 2,263 | 5.67× | 71% | 512.5 MiB |
| 16 | 50,000 | 2,376 | 5.95× | 37% | 905.5 MiB |

The 1–8 worker curve corresponds to an effective serial/contended fraction of approximately 6% under Amdahl's law. That is enough to turn an ideal 8× speedup into approximately 5.7×. The 16-worker plateau is expected because the allocation exposes only nine CPUs.

A separate 32-worker run used 32,000 reads in 500-record batches so every worker received work. It successfully initialized all 32 workers and processed 2,223 reads/s while reaching 1,259 MiB peak RSS. This verifies that an explicit worker count is not capped by the CLI; it says nothing favorable or unfavorable about scaling on a genuine 32-core machine because this run oversubscribed the nine-CPU allocation.

## Sources of sublinearity

There is no global alignment lock and no cross-read biological calculation. The non-worker-local work is:

- one main thread parses and, when relevant, decompresses input; constructs and UTF-8-encodes batches; dispatches and receives worker messages; preserves batch order; and writes one output stream;
- every worker has a private WASM linear memory, germline database, and k-mer index, so memory use grows approximately linearly and high worker counts contend for last-level cache and memory bandwidth;
- the production CLI decodes AIRR bytes to JavaScript strings in workers and structured-clones those strings to the parent, whereas this benchmark uses transferable `ArrayBuffer` output and is therefore optimistic for messaging;
- assignment-only scheduling submits work round-robin rather than always selecting the next free worker, and ordered output can expose a slow-batch head-of-line stall;
- one gzip stream and one AIRR TSV writer ultimately cap throughput; and
- single-core turbo frequency normally exceeds all-core frequency, while logical/shared CPUs do not behave like independent physical cores.

Large batches amortize message overhead but require enough reads to occupy every worker. With the CLI default of 2,000 reads per batch, fewer than 64,000 reads cannot give 32 workers even one complete batch each.

## Worker-count rules

- `swig-cli --vdj` with no thread option defaults to `min(4, availableParallelism)`.
- `--workers 0` selects `min(8, availableParallelism)`.
- A positive `--workers N` is exact and has no CLI maximum; `--workers 32` starts 32 workers.
- The automatic ceiling is deliberate because every additional worker duplicates reference-index memory. It is not a scientific or WASM limit.

Near-linear 32-core scaling would require reducing the main-thread fraction and per-worker memory/cache pressure—most immediately by transferring AIRR byte buffers directly to the writer, using availability-driven scheduling, and avoiding redundant WASM compilation/index startup. These are execution changes and need not alter assignment or annotation semantics.
