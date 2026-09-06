# Swig 0.38.7

Reconstructed from the supplied 0.38.5 source release and the surviving 0.38.6 conversation. The exact lost patch was unavailable. This is a tested reconstruction with newly selected, seeded dropout panels, not a byte-identical restoration.

- AER-R and R-optimized are the Web and CLI assignment defaults. Explicit legacy settings remain supported.
- Multi-file assignment progress uses consistent file-size estimates when counts are unknown, discounts read-ahead by actual AIRR commits, and reserves final completion for indexing. No counting pass is required.
- Personalized V inference adds HS5F targeting/substitution tables, aligned leave-gene-out SHM calibration, linked-haplotype evidence and sequence-search penalties, compatible cross-gene competition, conditional reference-backbone checks, and raw-sequence V trim/N uncertainty.
- The browser worker and new `swig-cli personalized-germline` command use the same engine. Evidence and personalized V FASTA exports are included.
- The active-allele setting remains a computational guard, not a diploid/copy-number prior.

## Empirical validation

Each animal used an independent seeded 300,000-read sample, full rhesus KIMDB versus six dropped V records, AER-R/R-optimized, eight workers, indel-aware collapse and lineage assignment. Full and dropout arms used identical raw reads.

| Animal | Recovery | Unrelated dropout-induced calls | Shared non-target candidates |
|---|---|---:|---:|
| ERR4238110 | 5/6 exact; sixth differs only at terminal base 296 | 0 | 10 |
| ERR4238104 | 6/6 exact | 0 | 8 |

The shared background candidates are **not independently validated alleles**. Zero dropout-induced extras is not a claim of zero overall false positives. The exact panels, settings, source hashes and full outputs are in `benchmarks/personalized-germline/results`. See `BENCHMARK_PERSONALIZED_GERMLINE_0.38.7.md` for interpretation and reproduction commands.

## Verification and packaging

All **267 tests** pass, including new progress, homology, cross-parent competition, donor-isolation, terminal-coordinate and gzip CLI recovery checks. TypeScript checking, generated CLI bundling and the production Web build pass. Vite retains its pre-existing large-main-chunk advisory; the build succeeds.

The ZIP contains complete source, generated Node CLI, bundled WASM assets, references, tests, methods and the GitHub Actions release workflow. It excludes dependencies, build output and large raw/intermediate repertoire files. The assignment WASM is unchanged from the supplied 0.38.5 package. Platform executables have not been compiled or published here; the included workflow builds them when tagged on GitHub. Million-lineage browser memory/performance has not been validated.
