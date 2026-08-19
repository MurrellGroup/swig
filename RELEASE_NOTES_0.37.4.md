# Swig 0.37.4

## Parallel post-analysis

- Browser and CLI exact collapse use deterministic key shards for large runs.
- FAD-compatible, conservative, and indel-aware denoising dynamically schedule complete independent study/C/locus/V/J partitions across the configured worker pool.
- CLI CHMMAIRRa now distributes bounded row batches across the same configured worker count.
- Serial and parallel routes share the same pure kernels, preserve ordinal tie breaks, and produce identical representatives and counts regardless of task completion order.
- A single large denoising partition remains indivisible because splitting it would alter global parent/centroid decisions.

## Lineage exploration

- Lineage selectors show a representative nucleotide/amino-acid CDR3. After SHM calculation, this is taken from the exact lowest-V-SHM member rather than the bounded plot sample.
- Checkboxes and **Open together** load several original lineages into one alignment/tree/UCA workbench without rewriting their assignments. Lazy CLI lineage studies fetch and verify only the selected indexed ranges.
- CDR3 neighbour search now offers equal-length Hamming or indel-aware Levenshtein distance and independently configurable V/J compatibility (`both`, `either`, or `ignore`). Results expose edit count, length delta, and V/J agreement.

## Assignment telemetry carried forward

- The stable browser progress card reports committed queries/s, active workers, initialization, streaming, drain, and finalization phases without layout jumping.

Scientific scoring, denoising decisions, tie breaks, and AIRR output ordering are unchanged by the parallel execution paths.
