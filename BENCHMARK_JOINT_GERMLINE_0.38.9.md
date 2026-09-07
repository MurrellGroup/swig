# Joint inherited / SHM model — 0.38.9 validation

## Outcome

The new method is an independent experimental box, not a replacement for the existing personalized-reference method. It has a consistent likelihood and a conservative held-out comparison, but **does not meet the intended recovery target**. All 24 held-out targets were nominated; 11 pass the model-conditional multiple-testing rule. The existing method recovered all 24 identifiable target haplotypes (22 exact full sequences and two with terminal ambiguity).

| Animal | Panel | Nominated | Positive held-out evidence | Pass evidence rule | Non-target proposals passing |
|---|---|---:|---:|---:|---:|
| ERR4238110 | Original | 6/6 | 4/6 | 3/6 | 0 |
| ERR4238110 | Lower support | 6/6 | 2/6 | 0/6 | 0 |
| ERR4238104 | Original | 6/6 | 5/6 | 5/6 | 0 |
| ERR4238104 | Lower support | 6/6 | 4/6 | 3/6 | 0 |

The passing count is based on log evidence >= log(M/0.05), with M training-nominated hypotheses per subject. This is a model-conditional rule, **not an empirically established 5% false-positive rate**. The full-reference controls have no non-reference proposal passing this rule. Every dropout arm has one non-target proposal with merely positive evidence; none passes the rule. These observations do not establish improved overall specificity: the new method also rejects many genuine alleles.

All target comparisons use the identifiable sequence prefix, excluding the final two V bases. The new method does not claim exact full-length germline recovery at those terminal bases.

## Controlled inputs

Both animals use their original deterministic 300,000-read samples and the preserved AER-R/R-optimized, indel-aware collapse and lineage-assignment outputs. The full, original-dropout and lower-support-dropout arms are paired within each animal. Animals are fitted independently; no repertoire evidence is pooled between them. The new discovery model additionally splits the selected lineage representatives into training/test halves.

The lower-support panels contain different physical genes from the original panel and were selected from full-reference support in the previous validation, before these new fit outcomes. They contain 20–26 near-exact distinct CDR3s in ERR4238110 and 23–31 in ERR4238104. Low absolute support is distinct from a low within-gene allele fraction; the synthetic checks below address that distinction. These are development/regression data, not a blinded external validation.

No assignment, denoising or lineage algorithm changed in 0.38.9. Reusing their exact outputs isolates the downstream discovery change. The previous 0.38.8 report, prepared references, fixed panels, raw-read replay script and all six compact AIRR replay inputs remain in the package.

## Per-target results

Null/empty log evidence means no supported numerical evidence value was reported; it is not a probability of genomic absence. The optimization gap is the maximum certificate gap among fitted hazard models, not a posterior uncertainty interval.

| Animal / panel | Target | Near-exact CDR3s | Log evidence | Required | Null optimization gap | Pass |
|---|---|---:|---:|---:|---:|---|
| ERR4238110 / dropout | IGHV1-NL_2*01_S3568 | 25 | 4.877 | 9.739 | 0.000 | No |
| ERR4238110 / dropout | IGHV3-37*01_S4420 | 78 | -21.424 | 9.739 | 0.002 | No |
| ERR4238110 / dropout | IGHV4-NL_34*01_S4637 | 173 | 41.105 | 9.739 | 1.028 | Yes |
| ERR4238110 / dropout | IGHV2-118*01_S4527 | 381 | -49.195 | 9.739 | 0.643 | No |
| ERR4238110 / dropout | IGHV5-157*01_S0308 | 586 | 18.456 | 9.739 | 0.374 | Yes |
| ERR4238110 / dropout | IGHV3-50*01 | 1170 | 237.554 | 9.739 | 0.019 | Yes |
| ERR4238110 / rare | IGHV4-67*02_S8671 | 20 | 3.889 | 9.793 | 0.000 | No |
| ERR4238110 / rare | IGHV3-36*01_S7141 | 21 | -12.441 | 9.793 | 0.000 | No |
| ERR4238110 / rare | IGHV3-28*01_S5703 | 23 | -37.317 | 9.793 | 0.003 | No |
| ERR4238110 / rare | IGHV1-84*01_S8557 | 24 | 1.761 | 9.793 | 0.000 | No |
| ERR4238110 / rare | IGHV3-9*01 | 26 | -19.040 | 9.793 | 0.000 | No |
| ERR4238110 / rare | IGHV3-NL_19*01_S4000 | 26 | -5.779 | 9.793 | 0.000 | No |
| ERR4238104 / dropout | IGHV2-NL_1*01_S7874 | 319 | -32.883 | 9.505 | 1.612 | No |
| ERR4238104 / dropout | IGHV3-NL_1*01_S9854 | 490 | 44.598 | 9.505 | 0.360 | Yes |
| ERR4238104 / dropout | IGHV3-NL_11*01_S3736 | 984 | 205.867 | 9.505 | 0.785 | Yes |
| ERR4238104 / dropout | IGHV3-NL_17*01_S4736 | 569 | 55.374 | 9.505 | 0.006 | Yes |
| ERR4238104 / dropout | IGHV4-149*01_S1166 | 553 | 75.566 | 9.505 | 20.681 | Yes |
| ERR4238104 / dropout | IGHV3-94*01_S3427 | 471 | 66.782 | 9.505 | 0.004 | Yes |
| ERR4238104 / rare | IGHV3-73*02_S7884 | 23 | 16.568 | 9.507 | 0.000 | Yes |
| ERR4238104 / rare | IGHV5-15*01 | 23 | -62.712 | 9.507 | 0.022 | No |
| ERR4238104 / rare | IGHV3-88*01 | 24 | -6.654 | 9.507 | 0.000 | No |
| ERR4238104 / rare | IGHV1-170*01 | 30 | 9.602 | 9.507 | 0.000 | Yes |
| ERR4238104 / rare | IGHV1-138*01_S2593 | 31 | 11.577 | 9.507 | 0.000 | Yes |
| ERR4238104 / rare | IGHV3-172*01_S8824 | 31 | 3.089 | 9.507 | 0.000 | No |

## Synthetic checks

The generator uses known allele truth, uniform background SHM (deliberately different from the fitted HS5F kernel), exceptional hotspots, adjacent correlated repair events, sequencing error and independent lineage IDs. Mixed repertoires have 35% unmutated lineages and mean active exposure 0.035; the memory-like setting has 5% and 0.08. Each target gene has 3,000 lineages, with four additional calibration genes of 500 each. These few deterministic stress cases diagnose failure modes; they do not estimate general sensitivity/FDR.

| Seed / setting | True allele | Carriers | Nominated | Log evidence | Required | Pass | Other positive proposals |
|---|---|---:|---|---:|---:|---|---:|
| 203 / memory | 1%, 2 site(s) | 27 | Yes | -53.803 | 7.783 | No | 0 |
| 201 / mixed | 0%, 2 site(s) | 0 | No | — | 8.132 | No | 0 |
| 202 / mixed | 1%, 1 site(s) | 36 | Yes | -45.646 | 8.211 | No | 0 |
| 204 / mixed | 5%, 2 site(s) | 154 | Yes | 4.454 | 8.172 | No | 0 |

The recurrent-SHM null produces no positive proposal. However, the 1% single-site and 1% linked memory-like alleles fail, and even the 5% linked allele has only weak positive evidence below the rule. This is a substantive sensitivity failure, not evidence that all rejected sequences were artifacts.

## Interpretation

The redesign separates an inherited sequence distribution from an SHM-dependent somatic endpoint distribution under the same normalized CTMC/rearrangement likelihood. Candidates are generated on training data and evaluated on held-out data. Null optimization uses a concavity certificate, and no later BIC/set-fit objective overrides the evidence. A read-order-dependent calibration choice was corrected: the complete tested compatibility component is excluded from nuisance calibration.

The remaining limitation is biological/statistical, not a reason to loosen the evidence threshold. The somatic endpoint distribution can assign mass to whole linked sequences without a detailed mutation-path cost. This deliberately broad alternative can absorb real inherited signal. Splitting data and estimating a large nuisance mixture also costs substantial predictive power. The failed synthetic cases show that the resulting robustness is too expensive in sensitivity. A useful next model must constrain somatic paths using justified mutation/selection structure and calibrate its uncertainty, rather than treating arbitrary endpoint sequences as equally available somatic explanations. The current implementation is a coherent conservative baseline for that work, not a claim that the problem is solved.

## Software verification and use

All 278 tests pass (including 180 post-analysis tests), and TypeScript/production build pass. The CLI smoke test runs both methods on the same input, checks the new result schema/options, and verifies that the existing result is unchanged. Regression tests cover stochastic kernels, the concavity bound, inherited-versus-somatic discrimination, nonidentifiability without unmutated information, component calibration invariance, and independent UI/session wiring.

The Web boxes have independent workers and saved results. The new box exports evidence JSON, not a confident personalized FASTA. AER-R/R-optimized defaults and the multi-file assignment progress correction are retained. See [the full model specification](public/methods/10_JOINT_INHERITED_SHM_MODEL.md).

Run the new method with:

```
node cli/swig-cli.mjs joint-germline --airr processed.airr.tsv.gz --v-reference V.fasta --out joint-result
python benchmarks/personalized-germline/v0389/replay.py --animal ERR4238110 --out replay-8110
python benchmarks/personalized-germline/v0389/replay.py --animal ERR4238104 --out replay-8104
```

Run the animals separately. The replay script evaluates all three reference arms using the exact representative inputs from v0388. Synthetic reproduction is documented in `research/unified/README.md`.

Fits shared the machine with same-animal control fits and verification work; elapsed times in JSON are not isolated performance benchmarks. The implementation remains TypeScript and has not been validated for a million lineages.

## Statistical reference

The split-likelihood construction follows [Wasserman, Ramdas and Balakrishnan, Universal inference (2020)](https://www.pnas.org/doi/10.1073/pnas.1922664117). Its formal guarantees require a correct null family and appropriate independent observations. Training-estimated nuisance distributions, alignment-conditioned data and lineage representative selection mean that those assumptions cannot simply be presumed here.
