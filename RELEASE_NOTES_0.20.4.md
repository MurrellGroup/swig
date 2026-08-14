# Swig 0.20.4

## Streaming FASTQ quality control

- Analysis parameters now contain a collapsed **FASTQ quality filter** stage.
  It is optional and disabled initially; its default maximum expected errors is
  `0.01` when enabled.
- Expected errors are evaluated per read as `Σ 10^(-Q/10)`. Reads are retained
  when the sum is less than or equal to the selected threshold. Standard
  Phred+33 and explicit legacy Phred+64 encodings are supported; invalid quality
  characters stop the run with the read and base position identified.
- Optional 3′ terminal-window trimming runs before filtering. Window width,
  minimum terminal mean Phred, and minimum retained length are independently
  configurable. Expected errors are evaluated only over retained bases.
- FASTA and AIRR inputs pass through unchanged when the stage is enabled. The
  controls and run manifest state this explicitly for non-FASTQ and mixed-format
  analyses.

## Throughput and provenance

- Filtering is integrated into the gzip-capable streaming parser before random
  subsampling and before WASM batch allocation. Phred error probabilities are
  precomputed, rejected records never receive a canonical output string, and
  trimmed-tail contributions are subtracted from the full-read sum without a
  second retained-read scan.
- When quality filtering and random subsampling are both enabled, the exact
  reservoir sample is drawn only from quality-retained reads.
- Progress reports raw reads scanned, reads retained, and reads rejected while
  filtering. Final results, portable sessions, project checkpoints, and project
  logs retain the settings and aggregate counts for expected-error rejection,
  minimum-length rejection, trimmed reads/bases, and non-FASTQ passthrough.
- A 50,000-read × 120-nt regression benchmark covers the hot path and verifies
  bounded 1,000-read batches. Correctness tests cover inclusive thresholds,
  trim-before-filter order, minimum retained length, mixed-format passthrough,
  and quality-retained reservoir sampling.
