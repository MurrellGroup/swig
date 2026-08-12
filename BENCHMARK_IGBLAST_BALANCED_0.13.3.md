# Swig 0.13.3 IgBLAST-balanced profile benchmark

## Objective

The optional `IgBLAST-balanced` profile was selected by constrained optimization on the supplied 100,000 simulated human IGH records:

1. maximize agreement with the supplied IgBLAST calls;
2. require mean V/D/J first-call truth accuracy to exceed IgBLAST's mean;
3. independently require mean V/D/J ambiguity-aware truth accuracy to exceed IgBLAST's mean.

This is a tuning-set result. No held-out data were supplied or used, and the profile is not the default.

## Scoring

- **Truth first call:** one point when the first reported allele is in the simulated truth set; an empty call receives one point only for an empty truth set.
- **Truth fair:** the number of distinct reported alleles in the truth set divided by the number of reported alleles. An empty call receives one point only for an empty truth set.
- **IgBLAST first agreement:** one point for identical first calls, including two empty calls.
- **IgBLAST fair agreement:** overlap with the IgBLAST call set divided by the number of Swig calls, including one point for two empty calls.

Both truth constraints were enforced. Increasing or removing alternative calls cannot satisfy both scores merely by changing ambiguity reporting.

## Selected rule

The underlying D/J configuration is the 0.13.2 IgBLAST-agreement profile:

- D: match `+2`, mismatch `-4`, gap open `-11`, gap extend `-1`, three candidates, minimum five-nucleotide exact run;
- J: match `+2`, mismatch `-4`, gap open `-13`, gap extend `-1`, two candidates.

For a selected D hit whose longest exact run is exactly five nucleotides, the balanced profile removes the D call when:

```text
j_sequence_start - v_sequence_end <= 11
```

In AIRR coordinates this corresponds to at most ten unassigned query bases between the V and J alignments. The rule leaves every D hit with a six-base-or-longer exact run unchanged and retains a five-base hit in a longer junction.

Among the evaluated exact-five evidence rules, this simple span boundary had the highest combined first/fair IgBLAST agreement while satisfying both truth constraints. A less restrictive span of 10 failed the constraints; a span of 12 was more accurate but agreed less often. Candidate-count-weighted rules were also dominated. A J-score sweep over mismatch penalties `-3` through `-6`, gap-open penalties `-9` through `-13`, one through six candidates, and minimum lengths five through ten did not create a higher-agreement joint J/D solution.

## Full tuning-set truth accuracy

Percentages are over all 100,000 records.

| Method/profile | V first | D first | J first | Mean first | V fair | D fair | J fair | Mean fair |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| IgBLAST | 93.9370 | 65.6250 | 98.6230 | 86.0617 | 93.9286 | 65.7007 | 98.2790 | 85.9694 |
| Swig agreement-only | 92.9110 | 65.6450 | 98.0100 | 85.5220 | 92.8541 | 65.7078 | 97.6750 | 85.4123 |
| **Swig IgBLAST-balanced** | **92.9110** | **67.3280** | **98.0100** | **86.0830** | **92.8541** | **67.3930** | **97.6750** | **85.9740** |

Relative to IgBLAST, the selected profile is higher by **0.0213 percentage points** for mean first-call accuracy and **0.0046 percentage points** for mean fair accuracy. The margin is deliberately small because agreement, not truth margin, is the optimized quantity.

The rule removed 2,168 D calls. Of those records, 1,732 had an empty simulated D truth set and 80 had an empty IgBLAST D call. D-negative specificity increased from 52.7603% in the agreement-only profile to 63.9043%; D-positive first-call accuracy changed from 68.0161% to 67.9580%.

## Agreement with IgBLAST

| Segment | First-to-first | Fair directional agreement |
|---|---:|---:|
| V | 97.8150% | 98.0758% |
| D | 96.8720% | 96.9440% |
| J | 99.3610% | 99.3815% |
| **Mean** | **98.0160%** | **98.1338%** |

The agreement-only profile remains available when reproducing IgBLAST calls is the sole objective; its mean first/fair agreement was 98.6670%/98.7863% on these records.

## Throughput

On the deterministic 5,004-record subset, the optimized single-worker WebAssembly path processed 262.6 reads/s with the agreement-only profile and 243.8 reads/s with the complete balanced AIRR consistency transformation, a 7.2% reduction. Browser runs still use the same bounded multi-worker queue and streaming output architecture. This local Node/WASI number is a relative implementation check, not a claim about browser or native SwiftIG throughput.

## Reproducibility

The selection and scoring utilities are:

- `tests/benchmark-simulated-accuracy.mjs`
- `tests/score-simulated-accuracy.mjs`
- `tests/compare-call-agreement.mjs`
- `tests/analyze-hybrid-d-evidence.mjs`
- `tests/select-balanced-d-rule.mjs`

Input SHA-256 values:

| File | SHA-256 |
|---|---|
| `seen_queries.fasta` | `ab8ce0d3d5961626d131df55fc8f2233546491090325ccfd9a46059c3c767ce4` |
| `seen_calls.tsv` | `4e1c8642a90866089128330ec937e0e5a927a88f9537995c66449314f45b681f` |
| `V 1(1).fasta` | `4614a3052669c442ca678b2d4c6be89f148354d19ef11813c90e874c7daa7a4e` |
| `D 1(1).fasta` | `70f9ac7b0a4c86e68e018b14c0eaf4bc32db802b877a6173503519b329df70fc` |
| `J 1(1).fasta` | `4bc809d2d0081adafb8c21fd478f989a133e09b97079dcb5092ff2944f904ca0` |

## Limitations

- The rule was calibrated on one simulated human IGH dataset and these exact references. It is not evidence of equivalent behavior for other species, loci, sequencing technologies or mutation regimes.
- The score advantage over IgBLAST is a tuning-set advantage, not independent validation.
- “More accurate” here refers to the predeclared mean across V, D and J under both scoring rules. V and J segment accuracy remain below IgBLAST on this dataset; the D improvement is what raises the combined mean.
- An empty D call means that the available sequence did not meet this profile's evidence threshold. It does not assert that the biological rearrangement lacked a D segment.
