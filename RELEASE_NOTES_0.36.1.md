# Swig 0.36.1

## Faster streaming CLI

- Hot worker requests and responses are flat primitive/string messages. This activates Bun's optimized `postMessage` path instead of the general structured-clone path; no AIRR fields or values change.
- Workers take the next available batch instead of receiving round-robin queues. This improves utilization when read lengths or batch runtimes differ while ordered output remains unchanged.
- Removed an unconditional copy of every AIRR result body.
- Increased bounded input, gzip, and output buffering to reduce JavaScript callbacks and allow compute and filesystem writes to overlap. Output remains one ordered AIRR stream; no temporary collation files are introduced.
- Automatic batch sizing now uses 2,000 records for one or two workers, 1,000 for three or four, and 500 for five or more. `--batch-records N` and `annotation.batchRecords` still override this; zero selects automatic sizing.
- Release builds now use Bun minification and bytecode compilation to reduce startup parsing. Linux releases retain the compatible baseline x64 binary and add `swig-cli-linux-x64-modern`, built for Haswell/AVX2 CPUs.

The optimized and 0.36.0 paths produced byte-identical AIRR output in assignment-only, prepared-annotation, plain-FASTA, and gzip-input benchmarks. See `BENCHMARK_CLI_THROUGHPUT_0.36.1.md`.
