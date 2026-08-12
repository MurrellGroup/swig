# Swig 0.13.2 IgBLAST-agreement calling profile

## Purpose

Swig 0.13.2 adds an optional `IgBLAST-agreement` calling profile. It was tuned only to reproduce the supplied IgBLAST calls; simulation truth columns were not read during parameter selection. The existing ground-truth-optimized profile remains the default and its output is unchanged.

This profile is a SwiftIG scoring configuration. It is not an IgBLAST implementation and it does not reproduce IgBLAST's search procedure, germline preprocessing, or reporting rules.

## Data and selection protocol

- Input: 100,000 supplied simulated human IGH rearrangements, their supplied IgBLAST calls, and the exact supplied V/D/J FASTA files.
- Search: a deterministic 5,004-record subset was used for the initial parameter grid and a deterministic 20,017-record subset for confirmation.
- Objective: improve first-call agreement, ambiguity-penalized fair agreement, exact-set agreement, and called/uncalled agreement without selecting parameters merely by adding more alternative calls.
- Final report: the locked profile was run once across all 100,000 records in four deterministic partitions. This is a tuning-set report, not held-out validation.

The fair score for one record is

`|Swig calls ∩ IgBLAST calls| / |Swig calls|`.

Thus a Swig call of `A,B` receives one-half credit against an IgBLAST call containing only `A`; an empty prediction receives zero when IgBLAST made a call; and two empty calls receive full agreement. Reciprocal coverage and Jaccard agreement were also inspected so reduced candidate depth could not appear favorable solely by suppressing alternatives.

## Locked settings

| Segment | Default truth-optimized profile | Optional IgBLAST-agreement profile |
| --- | --- | --- |
| V | `+2/−3/−5/−1`; 3 candidates | unchanged |
| D | `+2/−5/−5/−1`; 7-nt exact-run floor; 6 candidates | `+2/−4/−11/−1`; 5-nt exact-run floor; 3 candidates |
| J | `+2/−5/−11/−1`; 3 candidates | `+2/−4/−13/−1`; 2 candidates |

Scoring tuples are match, mismatch, gap-open, and gap-extension scores. All other alignment, ambiguity, AIRR, strand, reference, and double-D behavior is shared.

## Complete 100,000-record IgBLAST agreement

| Segment | Profile | First call | Fair | Exact set | Call-status agreement | Swig call rate | IgBLAST call rate | Mean calls/record |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| V | default | 97.814% | 98.075% | 96.306% | 100.000% | 100.000% | 100.000% | 1.0981 |
| V | **IgBLAST-agreement** | **97.815%** | **98.076%** | **96.307%** | **100.000%** | 100.000% | 100.000% | 1.0981 |
| D | default | 83.339% | 83.328% | 79.447% | 85.857% | 77.330% | 91.297% | 1.1292 |
| D | **IgBLAST-agreement** | **98.825%** | **98.902%** | **93.569%** | **99.279%** | **91.016%** | 91.297% | 1.3739 |
| J | default | 99.298% | 99.320% | 99.276% | **99.803%** | **99.803%** | 100.000% | 1.0174 |
| J | **IgBLAST-agreement** | **99.361%** | **99.382%** | **99.324%** | 99.583% | 99.583% | 100.000% | 1.0148 |

For the 91,297 records on which IgBLAST called D, the compatibility profile called D in 99.451% and achieved 98.954% first-call, 99.038% fair, and 93.197% exact-set agreement. For the 8,703 IgBLAST D-negative records, it also omitted D in 97.472%. The profile therefore improves agreement principally by recovering IgBLAST-supported D calls, not by indiscriminately reporting D on every record.

Relative to the default, D sets changed on 18,625 records. A fair score improved on 17,059 and worsened on 197. The compatibility profile recovered 13,558 IgBLAST-positive D calls that the default omitted, while adding a D to 175 IgBLAST-negative records. J sets changed on 353 records; fair agreement improved on 87 and worsened on 17.

## Default-path regression check

The pre-switch production WASM and the dual-profile WASM were run with the default profile on the same deterministic 5,004 records. Their exported `sequence_id`, `v_call`, `d_call`, and `j_call` tables were byte-for-byte identical (SHA-256 `5776209067edc4697fc6769b30c78e38a5f12f46904f494ed2c9d50e31c2d264`). The optional profile must be selected explicitly in **Analysis parameters → Calling profile**.

## Reproduction and hashes

The browser-facing production binary was evaluated through `tests/benchmark-simulated-accuracy.mjs --target igblast --calling-profile igblast_compatible`. Full default-versus-profile comparisons use `tests/compare-call-agreement.mjs`.

| File | SHA-256 |
| --- | --- |
| Swig 0.13.2 production WASM | `8a3d82c6832d83a8a8274f5bd29ad8b6035c1ff9f5b0c68b591fc5baf0eb4e26` |
| `seen_queries.fasta` | `ab8ce0d3d5961626d131df55fc8f2233546491090325ccfd9a46059c3c767ce4` |
| `seen_calls.tsv` | `4e1c8642a90866089128330ec937e0e5a927a88f9537995c66449314f45b681f` |
| V reference | `4614a3052669c442ca678b2d4c6be89f148354d19ef11813c90e874c7daa7a4e` |
| D reference | `70f9ac7b0a4c86e68e018b14c0eaf4bc32db802b877a6173503519b329df70fc` |
| J reference | `4bc809d2d0081adafb8c21fd478f989a133e09b97079dcb5092ff2944f904ca0` |

## Limitations

The same simulated human IGH data were used for tuning and the complete report. Agreement is not biological accuracy, and IgBLAST itself is not ground truth. The observed settings may not transfer to light chains, TCRs, other species, different read lengths, different SHM/error regimes, or different reference databases. Swig therefore keeps the truth-optimized profile as the default, records the selected profile in saved sessions and result manifests, and applies that same recorded profile when post-analysis query sequences are assigned.
