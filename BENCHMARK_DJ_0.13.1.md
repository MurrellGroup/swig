# Swig 0.13.1 simulated V(D)J accuracy benchmark

## Scope

This report covers the 100,000 simulated human IGH queries, truth/call table, and exact V/D/J references supplied for tuning. As requested, all supplied records were available for tuning; this is therefore a tuning-set result, not a held-out estimate. An unseen set should be used to assess generalization.

Parameter exploration used deterministic 5,004- and 20,017-record subsets before one final 100,000-record evaluation. The explored surface included D/J candidate depths, match/mismatch and gap penalties, D exact-run evidence thresholds, minimum J alignment lengths, and a J 3-prime reference-end penalty. Alternatives that merely increased the number of reported hits were rejected.

## Scoring rules

Comma- and slash-delimited calls are treated as atomic allele alternatives.

- **First prediction:** 1 when the first predicted atomic allele belongs to the truth set, otherwise 0.
- **Fair score:** `number of unique predicted alleles in the truth set / number of predicted alleles`. This gives partial credit but penalizes extra alternatives, so emitting a longer list cannot improve the expected score without adding correct alleles.
- For a truth record with no D, an empty D prediction scores 1 and a non-empty prediction scores 0 under both rules.

## Full 100,000-record result

All values below are percentages. The Swig 0.13.0 baseline is the output of the attached production WASM itself, not the precomputed SwiftIG columns in the supplied table.

| Method | V first | V fair | D first | D fair | J first | J fair |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Swig 0.13.0 WASM | 92.901 | 92.845 | 63.750 | 63.844 | 96.021 | 95.838 |
| **Swig 0.13.1 tuned WASM** | **92.911** | **92.855** | **70.825** | **70.882** | **98.053** | **97.712** |
| Supplied IgBLAST calls | 93.937 | 93.929 | 65.625 | 65.701 | 98.623 | 98.279 |

Relative to Swig 0.13.0, the tuned defaults improve D by 7.075 percentage points under first-call scoring and 7.038 points under fair scoring. J improves by 2.032 and 1.874 points, respectively. Relative to the supplied IgBLAST calls, tuned Swig is 5.200 points higher on overall D first-call accuracy and remains 0.570 points lower on J.

The small V change is indirect: V scoring itself is unchanged, but J scoring participates in orientation and compatible V/J pair selection.

### D behavior and ambiguity

| Method | D exact set | D call rate | D multi-call rate | Mean D calls/read | D-positive first | D-positive fair | D-negative specificity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Swig 0.13.0 WASM | 56.384 | 89.159 | 27.881 | 1.3507 | 64.759 | 64.870 | 58.268 |
| **Swig 0.13.1 tuned WASM** | **62.982** | 77.330 | **23.407** | **1.1292** | 66.787 | 66.854 | **92.768** |
| Supplied IgBLAST calls | 56.645 | 91.297 | 35.523 | 1.4631 | **68.080** | **68.170** | 52.284 |

There are 84,458 truth-positive and 15,542 truth-negative D records. The overall D gain is not hidden multi-hit inflation: Swig's multi-call rate falls by 4.474 points. It comes from a 2.028-point improvement among D-positive reads plus a 34.500-point improvement in D-negative specificity. IgBLAST remains 1.293 points higher than tuned Swig on D-positive first-call accuracy, while tuned Swig is much less likely to force a D call into a truly D-negative junction.

### J behavior and ambiguity

| Method | J exact set | J call rate | J multi-call rate | Mean J calls/read |
| --- | ---: | ---: | ---: | ---: |
| Swig 0.13.0 WASM | 95.121 | 98.886 | 1.456 | 1.0036 |
| **Swig 0.13.1 tuned WASM** | **96.762** | 99.803 | 1.914 | 1.0174 |
| Supplied IgBLAST calls | 97.317 | 100.000 | 1.938 | 1.0196 |

The tuned J ambiguity profile is almost identical to IgBLAST's rather than being expanded beyond it.

## Selected defaults

| Segment setting | 0.13.0 | 0.13.1 |
| --- | ---: | ---: |
| D match | 2 | 2 |
| D mismatch | -2 | **-5** |
| D gap open / extend | -5 / -1 | -5 / -1 |
| Minimum exact D run | 5 nt | **7 nt** |
| D candidate depth | 6 | 6 |
| J match | 2 | 2 |
| J mismatch | -2 | **-5** |
| J gap open / extend | -5 / -1 | **-11 / -1** |
| Minimum J aligned length | 10 nt | 10 nt |
| J candidate depth | 3 | 3 |

The seven-nucleotide D evidence requirement is the main specificity control. The stronger J mismatch and gap-open penalties suppress short accidental local HSPs in junction sequence while retaining long supported J alignments. A tested J 3-prime penalty had no measurable effect and was not shipped.

## Throughput check

In the same single-runtime Node/WASM harness over the deterministic 5,004-record subset, the attached 0.13.0 WASI-SDK production binary processed 207.4 reads/s and the 0.13.1 WASI-SDK production binary processed 240.1 reads/s. A 20,017-record tuned run processed 221.1 reads/s. These are individual timing observations rather than a formal repeated benchmark, but they show no throughput regression from the selected defaults. Browser worker-pool throughput depends on hardware, worker count, cross-origin isolation, input lengths, and active references.

The packaged production binary and the separately source-built tuning binary produced identical metrics on all 20,017 deterministic confirmation records. The full test suite also verifies the ordinary and opt-in double-D paths separately.

## Input and binary identities

| File | SHA-256 |
| --- | --- |
| `seen_queries.fasta` | `ab8ce0d3d5961626d131df55fc8f2233546491090325ccfd9a46059c3c767ce4` |
| `seen_calls.tsv` | `4e1c8642a90866089128330ec937e0e5a927a88f9537995c66449314f45b681f` |
| `V 1.fasta` (558 records) | `4614a3052669c442ca678b2d4c6be89f148354d19ef11813c90e874c7daa7a4e` |
| `D 1.fasta` (51 records) | `70f9ac7b0a4c86e68e018b14c0eaf4bc32db802b877a6173503519b329df70fc` |
| `J 1.fasta` (12 records) | `4bc809d2d0081adafb8c21fd478f989a133e09b97079dcb5092ff2944f904ca0` |
| Swig 0.13.0 production WASM | `5b5dbdf0bb5e3ae38d154dcb4179fadfb2a6c0f24f84cc09880c5f919fd61cc7` |
| Swig 0.13.1 production WASM | `78eb32bf1172b839830066aaef157a7ba4093abf177070c3b4a6cfbd6cd5d161` |

The supplied table's precomputed SwiftIG D column differs from the attached 0.13.0 WASM in 241 rows, exclusively in alternative-call tails in this comparison; the first D call is unchanged. For a reproducible baseline, this report uses freshly generated output from the attached WASM.

## Reproduction utilities

- `tests/benchmark-simulated-accuracy.mjs` runs the WASM directly against supplied references and implements the scoring rules above.
- `tests/score-simulated-accuracy.mjs` scores an existing call table or one or more partitioned prediction tables.
- `tests/extract-simulated-benchmark-subset.mjs` creates deterministic subsets by sequence-ID hash.

Run `npm test` and `npm run build` after rebuilding the core. The GitHub Pages workflow rebuilds `public/swiftig.wasm` with pinned WASI SDK 25, SIMD, LTO, and Binaryen before deployment.
