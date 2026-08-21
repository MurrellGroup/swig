# AER-R distributed-D audit — Swig 0.37.7

## Failure and correction

The 0.37.6 all-reference candidate oracle still applied the production D evidence gate. It therefore aligned every supplied D allele but discarded any alignment lacking the calling profile's consecutive exact run. Calling it simply “exhaustive” obscured that distinction.

The motivating complete-KIMDB macaque query contains this best D alignment:

```text
query       GGTCTAGAAGCATCTAC
            |||.|||.||||.||||
IGHD6-39    GGTATAGCAGCAGCTAC
```

It has 14 matches across 17 bases and a truth-profile D score of 19, but its longest exact run is four. The old truth profile emitted no D; the agreement profile's all-reference build could instead retain an unrelated five-base exact match. Both outcomes arose before joint V–D–J comparison.

AER-R now retains an already-scored D alignment that fails the ordinary exact-run floor only when all of these hold:

- at least 16 aligned substitution columns;
- at least 14 nucleotide matches;
- at least 80% identity across substitution columns; and
- score at least 19 under the fixed D evidence tuple `+2/−3/−13/−1`.

The fixed evidence tuple makes this rescue boundary independent of the selected calling profile. The alignment then enters the unchanged non-overlapping joint partition comparison; it is not accepted merely because it passed the aggregate gate. Ordinary AER, RIAT-MP, standard assignment, their exact-run floors, and their output are unchanged.

For the supplied query against the deduplicated union of the bundled rhesus and cynomolgus KIMDB 1.1 V/D/J collections (1,276 unique records), both AER-R calling profiles now emit:

| Segment | Call | AIRR query coordinates |
|---|---|---:|
| V | `IGHV4-NL_1*02_S8056` | 84–378 |
| D | `IGHD6-39*01` | 386–402 |
| J | `IGHJ2*01_S5087` | 414–453 |

The retained D is KIMDB positions 2–18. `np1` is `GATTAAC`; `np2` is `ATCCTCAATCT`. A deterministic WebAssembly regression runs both AER-R profiles and verifies those calls, coordinates, and alignments.

## Candidate-pruning audit

No additional dynamic programming was needed for the motivating read. Production AER-R had already admitted `IGHD6-39*01` through the weak-signal 3-mer candidate tier and scored it exactly; only the exact-run veto removed it.

Production D search uses 5-mer votes, conditionally adds 3-mer votes when the primary signal is weak, retains a 32-candidate safety pool, and computes exact affine scores in four SIMD lanes. The final AER-R rescue is based on candidate quality rather than the binary presence of a call:

1. Trace the ordinary leading strong-seed candidates.
2. If none passes the long distributed-evidence rule—even if a short D was called—score the deferred safety pool through selected-profile score 16, a conservative lower bound for the aggregate rule.
3. If the bounded pool still has no convincing long D, force the 3-mer tier over the complete locus-matched D set, score only new gene/diagonal hypotheses, and traceback only candidates capable of passing the bound.
4. Retain every alignment that actually passes the aggregate rule for the unchanged joint V–D–J comparison, even below the generic top-N cutoff.

An adversarial regression adds fifty perfect 7-nt D decoys. They all yield short calls, outrank the mutated `IGHD6-39*01` tract in seed voting, and overflow the ordinary candidate budget, but their affine scores are lower under both calling profiles. A no-D-only retry cannot solve that construction; the quality-triggered path recovers the 17-nt KIMDB alignment in both profiles.

Two broader alternatives were rejected. Unconditionally scanning the safety pool slowed clean reads by about 44%. Widening weak/no-D tracebacks from two to eight recovered no additional truth-compatible D call on the primary cohort, changed 74 weak D results, and raised the cohort fraction called on non-tandem truths retaining fewer than six D bases from 4.32% to 4.66%.

The compile-time all-reference candidate oracle was rebuilt with the new evidence rule. It removes V/D/J candidate-count pruning but deliberately retains the production evidence gates. Among 60 prioritized difficult records from a 1,000-read replicate, it produced no truth-compatible D rescue over production AER-R. Two classifications differed only in V; D remained unresolved in both production and oracle. This supports retaining the production D candidate search, but it is not held-out validation.

## Primary simulated cohort

The primary cohort is the same deterministic KIMDB rhesus IGH stress simulation used for 0.37.6: seed `913771`, 5,000 reads, hotspot-weighted SHM, V/D/J trimming, P/N additions, occasional tandem D, indels, sequencing substitutions, ambiguous bases, partial reads, flanks, and reverse complements. Metrics below compare the pinned 0.37.6 AER-R WebAssembly with 0.37.7 AER-R. A compatible D call includes alleles indistinguishable from the sampled allele over its retained truth tract.

| Metric | 0.37.6 AER-R | 0.37.7 AER-R |
|---|---:|---:|
| Exact V call | 95.72% | 95.72% |
| Exact J call | 97.28% | 97.28% |
| D detected, truth retained ≥6 nt | 98.42% | 98.55% |
| Exact D, truth retained ≥6 nt | 86.97% | 87.16% |
| Retained-tract-compatible D | 89.91% | 90.11% |
| Strong-D detected, truth retained ≥10 nt | 98.88% | 99.02% |
| Strong-D retained-tract-compatible | 92.18% | 92.40% |
| Mean fraction of truth D span overlapped | 92.57% | 92.75% |
| D-call rate among 310 non-tandem truths retaining <6 D nt | 69.68% | 69.68% |
| Same subthreshold calls as fraction of all 5,000 records | 4.32% | 4.32% |

Twelve AIRR D rows changed: nine previously wrong or absent calls became retained-tract compatible; two wrong calls changed to another wrong or absent result; and one already-compatible call retained its allele but gained a better boundary. No previously compatible call became incompatible. V and J outputs were unchanged.

Three independent 3,000-read seeds added 17 retained-tract-compatible corrections, one wrong-to-different-wrong change, zero compatible-to-incompatible changes, and no new subthreshold calls. Across those seeds and the primary cohort—14,000 simulated reads in total—the 0.37.7 change therefore produced 26 measured compatible corrections, three wrong-to-different-wrong/absent changes, one compatible boundary improvement, and zero measured losses of a previously compatible D call. These are tuning/development simulations, not a held-out accuracy estimate.

## Fully trimmed-D control

The simulator now supports a control in which a D allele is sampled but completely deleted before N addition. With tandem D disabled, all 5,000 records in seed `882001` contain no retained D nucleotide. P/N addition, V/J SHM, sequencing noise, indels, partial reads, and orientation sampling remain active.

Pinned 0.37.6 AER-R called a D in 53.50% of these controls; final 0.37.7 called 53.52%, one additional record (+0.02 percentage points). The high baseline rate is produced chiefly by the pre-existing six-base exact-match rule in random junction sequence. The additional control call is a chance 14/17 match that is sequence-evidence-equivalent to the motivating true D. Rejecting that pattern categorically would also reject the supplied case. The complete-set anti-crowding step introduced no further null call.

## Timing

One worker processed the same 5,000-read batches in six alternating, post-warm-up repetitions with independently initialized pinned-0.37.6 and 0.37.7 WASM runtimes.

| Cohort | 0.37.6 median | 0.37.7 median | Change |
|---|---:|---:|---:|
| Clean exact/near-exact reads | 3.769 s | 3.846 s | +2.0% |
| Difficult SHM/indel reads | 9.501 s | 9.629 s | +1.3% |

Clean output hashes were identical. The difficult output hash changed because the intended D corrections alter AIRR rows. The aggregate-evidence check itself is a short scan over already-traced alignments. Only a junction whose leading hypotheses contain no convincing long D pays for deferred score screening or, if still unresolved, the forced complete-D 3-mer collection. Even there, four-lane score-only affine evaluation precedes bounded traceback.

## Reproduction

```bash
npm run build:wasm
node tests/benchmark-aer-robustness.mjs \
  --count 5000 --seed 913771 --out tmp/aer-r-0.37.7

node tests/benchmark-aer-robustness.mjs \
  --count 5000 --seed 882001 --double-d-rate 0 \
  --fully-trimmed-d-rate 1 --out tmp/aer-r-null-0.37.7

bash tests/build-aer-robust-oracle-wasm.sh
node tests/benchmark-aer-robustness.mjs \
  --count 1000 --seed 913771 \
  --oracle-wasm .build/aer-robust-oracle.wasm --oracle-count 60
```

The oracle build requires `WASI_SDK`. Generated FASTA, truth, AIRR-row, and case reports are development artifacts and are not runtime dependencies.
