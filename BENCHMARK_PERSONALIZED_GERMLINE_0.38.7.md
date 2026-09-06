# Personalized-germline reconstruction validation — Swig 0.38.7

This release reconstructs the documented work after 0.38.5. The lost 0.38.6 patch, its exact nine-record panel and its sampling seed were unavailable. This is a new implementation and a new reproducible dropout experiment, not a byte-identical restoration or an exact replay of the old reported numbers.

## Experimental controls

- Each animal was processed independently; no reads, lineages, mutation-rate estimates or fitted allele frequencies were pooled between animals.
- Inputs contained 666,005 reads (ERR4238110) and 659,889 reads (ERR4238104). Each used a 300,000-read Algorithm R reservoir with Python `random.Random(1)`, restored to original ordinal order. Full and dropout runs used the identical sampled FASTA.
- The complete rhesus KIMDB 1.1 set contains 774 V, 52 D and 14 J records. Metadata preparation annotated every record. The CLI template label resolved to `Macaca mulatta_AG07107`; the tested V/D/J sequences are the uploaded-package **Macaca mulatta KIMDB** references.
- Both arms used AER-R (`aer_robust`), `r_optimized`, eight workers, the indel-aware collapse mode, and lineage assignment. There was no FASTQ filtering or additional subsampling of the prepared reservoirs. The exact resolved configurations are included.
- Six sequence-unique V records per animal were chosen from full-reference assignments before discovery: >=240 aligned V bases, >=20 distinct CDR3s with <=1% V mismatches, and a surviving same-gene same-length neighbor within 1–6 SNPs. Selection maximized support within distance strata, one allele per gene. Animal two had no qualifying 4–6 SNP targets under this rule.
- Inference consumed the **complete processed AIRR table**, not the viewer export limited to the top 10,000 lineage summaries. Exactly one lowest-current-V-SHM representative contributed per lineage.

Sample FASTA SHA256:

| Animal | SHA256 |
|---|---|
| ERR4238110 | `6a87b2a2caf75c5bfbeb5620210afab0b2a9a781d0c07a32b9e0c370f297af8a` |
| ERR4238104 | `287e74867e84f5238a8f9d348a454ff772ae815b038f88e16fa9eba1780ada62` |

## Paired results

| Animal | Exact target recovery | Additional terminal-only recovery | Full/dropout lineages | Unrelated dropout-induced calls | Shared non-target candidates | Dropout fit time |
|---|---:|---:|---|---:|---:|---:|
| ERR4238110 | 5/6 | 1 | 15,101 / 15,104 | 0 | 10 | 31.62 s |
| ERR4238104 | 6/6 | 0 | 13,395 / 13,390 | 0 | 8 | 23.93 s |

**Interpretation:** zero unrelated dropout-induced calls does not mean zero false positives. The ten and eight non-target dropout candidates also occurred with the full reference, but are not independently validated germline alleles. Several involve terminal bases and should be interpreted cautiously. Full-reference fits proposed 12 and 10 novel sequences respectively; candidate sets can lose background sequences when assignment/lineage structure changes.

ERR4238110 recovered IGHV5-157*01_S0308 only up to its interior haplotype; the inferred sequence differs at terminal nucleotide **296**. This is counted as **5/6 exact**, not 6/6 exact. The terminal model did not establish that last base. Its interior-equivalent sequence was also a background candidate in the full-reference control.

Fit times are wall-clock observations from this validation, not standardized throughput benchmarks. Some control fitting overlapped assignment work for the same animal. No million-lineage scalability claim is made.

## ERR4238110 dropout panel

| Target | Nearest surviving SNP distance | Near-exact baseline CDR3s | Recovery | Differing positions |
|---|---:|---:|---|---|
| IGHV1-NL_2*01_S3568 | 6 | 25 | exact | — |
| IGHV3-37*01_S4420 | 5 | 78 | exact | — |
| IGHV4-NL_34*01_S4637 | 4 | 173 | exact | — |
| IGHV2-118*01_S4527 | 3 | 381 | exact | — |
| IGHV5-157*01_S0308 | 2 | 586 | terminal-only | 296 |
| IGHV3-50*01 | 1 | 1170 | exact | — |

## ERR4238104 dropout panel

| Target | Nearest surviving SNP distance | Near-exact baseline CDR3s | Recovery | Differing positions |
|---|---:|---:|---|---|
| IGHV2-NL_1*01_S7874 | 3 | 319 | exact | — |
| IGHV3-NL_1*01_S9854 | 2 | 490 | exact | — |
| IGHV3-NL_11*01_S3736 | 1 | 984 | exact | — |
| IGHV3-NL_17*01_S4736 | 1 | 569 | exact | — |
| IGHV4-149*01_S1166 | 1 | 553 | exact | — |
| IGHV3-94*01_S3427 | 1 | 471 | exact | — |

## Reconstruction and diagnosis

The restored behavior includes explicit allele-alignment eligibility; within-subject competitive components across gene labels; HS5F targeting/substitution data; aligned leave-gene-out rate calibration with between-gene dispersion; candidate-specific linked-haplotype tests and discrete sequence-search cost; conditional checks against better-supported haplotypes on other reference backbones; raw-sequence terminal trim/N modeling; and the same engine in the Web worker and CLI.

The reconstruction initially missed a high-support one-SNP target despite strong proposal evidence. Inspection identified an incoherent likelihood-ratio correction: a ratio under a calibrated null was multiplied into an uncorrected parent likelihood. The final implementation instead uses a common calibrated categorical emission model for all sequences, with a separate linked proposal-validity test. This two-stage formulation is documented explicitly; it is not claimed to reproduce the lost implementation exactly.

Another extra was a projection on a different inherited reference backbone. A sequence-level conditional test now asks whether those inherited distinguishing bases are supported relative to the better-supported recovered allele. No diploid allele cap or case-specific threshold changes were used.

Boundary artifacts motivated learning junction composition and trimming sufficient statistics from other genes in the same animal. Equivalent trimming states are collapsed exactly for each tested site set, and single-site rates use analytic derivatives. These optimizations reduced the first-animal dropout fit to about 30 seconds.

## Reproduce

Build the CLI with `npm ci` and `npm run build:cli`, or use the bundled generated CLI with Node 24. Run animals sequentially:

```bash
python3 benchmarks/personalized-germline/run.py --animal ERR4238110 --fasta "ERR4238110.fasta(1).gz" --workdir validation-8110
python3 benchmarks/personalized-germline/score.py validation-8110 ERR4238110
python3 benchmarks/personalized-germline/run.py --animal ERR4238104 --fasta "ERR4238104.fasta(1).gz" --workdir validation-8104
python3 benchmarks/personalized-germline/score.py validation-8104 ERR4238104
```

The results directory includes prepared references, exact dropouts, resolved configurations, per-target scores, full model/evidence outputs and source hashes. Recorded configurations contain historical absolute paths; the reproduction script regenerates paths for its requested work directory. Raw reads and large AIRR intermediates are not duplicated in the source ZIP.

## Assignment defaults and progress

AER-R/R-optimized are now the Web and CLI defaults. Explicit saved settings remain explicit. Choosing a legacy assigner without a profile uses its compatible truth-optimized profile; explicitly incompatible combinations still fail.

Progress no longer mixes known read counts with a weight of one for unknown files. It uses consistent expanded-size estimates when needed, reads gzip footer sizes without a counting pass, and discounts parsed progress by acknowledged AIRR commits. Final indexing reserves completion until results are ready.

## Software verification

All **267 tests** and the production build pass; details are recorded in RELEASE_NOTES_0.38.7.md. The Web and CLI interfaces are checked, including a gzip AIRR end-to-end personalized-allele recovery test. The WASM assignment binary is unchanged from the supplied 0.38.5 package.
