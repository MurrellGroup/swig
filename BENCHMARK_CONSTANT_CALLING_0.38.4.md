# Swig 0.38.4 constant-calling and output audit

## Scope

This audit compares the released 0.38.3 WebAssembly binary with 0.38.4 on the two supplied 20,000-record KIMDB macaque simulations. Assignment used AER-R with the R-optimized profile, the deduplicated rhesus/cynomolgus V/D/J union, one Node process, 500-record batches, and the same inputs and benchmark harness for each paired run. Constant-enabled runs additionally used all 8 bundled `Macaca mulatta_AG07107` IGH C records (1,284 total references including V/D/J). Neither supplied simulation has a true `c_call`, so those runs measure false-positive control, not biological C sensitivity.

## Mechanistic change

0.38.3 aligned every long C reference against the complete oriented read. A local match anywhere upstream could therefore become a false C call, and AER-R could use such a candidate to clip the already-selected J endpoint. The full-read dynamic programs also dominated runtime when a C database was loaded.

0.38.4 treats C as a post-J annotation. It searches from the permitted one-base J/C overlap (`selected J query_end - 1`), never clips J, and never adds C score to V/D/J orientation or call selection. All C hits require at least 90% local identity. A hit then passes one of two calibrated evidence routes:

- ordinary/unanchored local hit: at least 30 aligned query nucleotides and `c_support <= 1e-5`, using the actual post-J search length and complete supplied C database; or
- short 5′-anchored hit: gapless, query offset ≤8 nt, C-reference offset ≤3 nt, and anchored `c_support <= 1e-4`, where the opportunity count is the number of locus-matched C references times the permitted query/reference offsets.

Twelve aligned bases is only a candidate-generation floor, not an acceptance threshold. The short-route score determines the effective minimum and becomes stricter as the C database grows. The alignment remains local, so a strong contiguous C tract is retained even when a downstream primer or unrelated tail mismatches the reference. An optional covered-prefix identity gate (off by default) includes the prefix bases that the local HSP skipped; failing it discards only C.

## Constant-negative controls

The fixed first 1,000 records from each no-C simulation were run with the C database loaded.

| Dataset | Version | False C calls | C false-positive rate | CPU seconds | Reads/s (wall) |
| --- | ---: | ---: | ---: | ---: | ---: |
| IgM-like | 0.38.3 | 29 / 1,000 | 2.9% | 33.750 | 29.91 |
| IgM-like | 0.38.4 | 0 / 1,000 | 0.0% | 5.114 | 209.43 |
| IgG-like | 0.38.3 | 40 / 1,000 | 4.0% | 34.095 | 29.62 |
| IgG-like | 0.38.4 | 0 / 1,000 | 0.0% | 7.152 | 147.95 |

The new C-enabled path used 84.8% less CPU on IgM-like (6.60-fold) and 79.0% less on IgG-like (4.77-fold). In a wider final-archive adaptive-route check, it made 0 false C calls among the first 5,000 records of either dataset. Those 5,000-record runs took 24.712 CPU-seconds (IgM-like; 24.162 wall seconds) and 31.282 CPU-seconds (IgG-like; 30.786 wall seconds).

The old false C calls clipped the J endpoint in 2/1,000 IgM-like and 5/1,000 IgG-like records. The post-J implementation removes that cross-region boundary effect.

## V/D/J invariance and ordinary speed

With no C database loaded—the normal setup used in the previous V/D/J audit—the complete 20,000-record 0.38.4 prediction files were byte-identical to 0.38.3 after accounting for the benchmark harness's newly added empty `c_call` column. Fair call accuracy therefore remains:

| Dataset | V accuracy | D accuracy | J accuracy |
| --- | ---: | ---: | ---: |
| IgM-like | 94.7154% | 76.8977% | 91.4748% |
| IgG-like | 87.4435% | 68.8658% | 86.2089% |

Accuracy gives `1/K` credit when one of `K` unique predicted sequence classes matches truth. Tandem-D truth records are excluded only from the primary single-D score (19,981 IgM-like and 19,978 IgG-like D-scored records).

Paired full-dataset timings show no ordinary-path slowdown:

| Dataset | Version | Wall seconds | CPU seconds | Reads/s |
| --- | ---: | ---: | ---: | ---: |
| IgM-like | 0.38.3 | 34.842 | 35.898 | 574.02 |
| IgM-like | 0.38.4 | 34.749 | 35.641 | 575.56 |
| IgG-like | 0.38.3 | 61.540 | 62.558 | 324.99 |
| IgG-like | 0.38.4 | 61.705 | 62.669 | 324.12 |

CPU changed by -0.72% and +0.18%, respectively, consistent with run-to-run noise around an unchanged no-C search path. The final full runs retained the same fair V/D/J accuracies shown above.

## Positive and delimiter regressions

Deterministic WebAssembly tests establish that:

- with the complete 111-record human IGH C database, an exact 15-nt 5′ C prefix passes while an exact 13-nt prefix remains below the database-size-adjusted evidence gate;
- a 72-nt C tract remains called while a downstream primer-like mismatch tail is clipped from the local C alignment;
- an 80%-identity post-J tract is rejected;
- a 99% optional covered-prefix threshold rejects a C tract with two leading mismatches that its local HSP omitted, while leaving V/D/J fields unchanged;
- supplying a C database cannot change V/D/J calls, orientation ranking, scores, alignments, or coordinates for the positive fixture; and
- optimized AIRR output replaces embedded tab/CR/LF characters in user-controlled fields so each record retains the header's field count.

The direct CLI now recognizes a `.gz` output suffix and writes a real gzip stream. Pipeline JSON/CLI can gzip annotated and processed AIRR tables, and Web Results/post-analysis exports offer bounded gzip streams for TSV/CSV/JSONL. Smoke tests verify gzip magic bytes, decompression, a constant field count, and that AIRR field 13 is `j_call`, so pipelines such as `zcat file.tsv.gz | cut -f13 | sort | uniq -c` work directly.

## Limits

The no-C simulations cannot estimate sensitivity or subclass accuracy on real constant-region reads. The regression fixtures test evidence-boundary behavior, not repertoire prevalence. Isoform resolution remains limited by the C paths present in the selected FASTA and by whether the read covers discriminating sequence; co-optimal calls remain explicit rather than being resolved arbitrarily.
