# Personalized germline validation — 0.38.8

## Outcome and limits

| Animal | Panel | Exact recovery | Interior-only recovery | Other candidates | Other candidates absent from full-reference control |
|---|---|---:|---:|---:|---:|
| ERR4238110 | Original six | 4/6 | 2/6 | 5 | 0 |
| ERR4238110 | Six lower-support alleles | 6/6 | 0 | 5 | 0 |
| ERR4238104 | Original six | 6/6 | 0 | 3 | 0 |
| ERR4238104 | Six lower-support alleles | 6/6 | 0 | 5 | 0 |

The full-reference fits each propose five additional sequences. These are **unvalidated**, not established KIMDB omissions. The original dropout runs previously proposed 10/8 other sequences in 0.38.7; they now propose 5/3. Full-reference recurrence demonstrates stability under the paired database perturbation, not biological truth or a measured overall false-discovery rate. These repertoires were used during model development, so this is regression validation rather than a blinded external validation.

The two incomplete recoveries are IGHV4-NL_34*01_S4637 and IGHV5-157*01_S0308: only their final V nucleotides differ. Under flexible trimming and initial junction composition, these bases are not identifiable. The FASTA retains the parent's unresolved terminal base; it must not be interpreted as a confident terminal genotype.

## Biological changes

A germline allele should explain a linked sequence inherited by independent rearrangements. SHM instead depends on mutation exposure, context, repair and selection. A low observed mutation count is not proof that the true mutation exposure was exactly zero. The new model learns an unmutated/active mixture outside the tested gene and integrates that uncertainty when testing candidate haplotypes. Aligned leave-gene-out event calibration remains; a heavy-tailed residual allows exceptional gene-specific mutation rates. Adjacent changes with overlapping five-base contexts also have a coordinated SHM-event null, so their linkage alone is insufficient evidence of an allele.

Terminal proposals must defeat both a flexible SHM explanation and a rearrangement explanation that profiles deletion lengths and initial inserted dinucleotides. Raw evidence stops before assigned D/J sequence. Changes confined to the last two bases remain unresolved. Partial haplotypes compete against accepted more complete haplotypes, rather than being promoted by missing coverage.

Recovery improvements are mechanistic: linked carriers may have additional SHM, proposal priority uses expected unmutated support, and the former 5% minimum within-gene novel fraction is removed. Four supporting representatives are still required. The candidate budget remains 64 per gene and is now exposed as `--max-candidates` and an advanced Web control. Budget warnings remain in these fits; completeness beyond the tested panels is not guaranteed. Eight active components per gene is a computational guard, not a diploid genotype prior.

A separate numerical defect was fixed: high-hazard match probabilities rounded to zero and made final mixture frequencies NaN. Match likelihoods now use direct log survival. Discovery scores and final set-selection gains are exported separately. The final set fit remains a simpler emission model; neither score is a calibrated allele-truth probability or FDR.

## Experimental setup

Each animal was processed independently with a deterministic 300,000-read reservoir (Python random.Random(1), Algorithm R, restored input order). All three reference arms use exactly the same reads. Assignment uses **AER-R / R-optimized**, eight workers, prepared rhesus KIMDB 1.1 (774 V, 52 D, 14 J), indel-aware collapse, then lineage assignment. No cross-animal evidence is pooled. The lower-support panels remove different physical genes from the original panel; selection uses full-reference support before inspecting revised discovery outcomes, not threshold adjustment on failed targets.

The lower-support targets have 20–26 distinct near-exact CDR3s in ERR4238110 and 23–31 in ERR4238104. This tests low absolute expression/support; it is not equivalent to a 0.5% minor allele within a highly expressed gene. Synthetic cases separately probe that situation.

Selected lineages (full/original/rare): ERR4238110 15,101 / 15,104 / 15,102; ERR4238104 13,395 / 13,390 / 13,396. Recorded fits took approximately 120–149 seconds while same-animal arms shared the machine; these are not isolated performance benchmarks. The engine remains TypeScript with a Web Worker and matching CLI, and million-lineage browser performance is not established.

Sample SHA256:

- ERR4238110: `6a87b2a2caf75c5bfbeb5620210afab0b2a9a781d0c07a32b9e0c370f297af8a`
- ERR4238104: `287e74867e84f5238a8f9d348a454ff772ae815b038f88e16fa9eba1780ada62`

## Remaining evidence

See [candidate-level support and mutation strata](benchmarks/personalized-germline/v0388/candidate-evidence.md). They are heterogeneous: some have near-complete linked support with several zero-other-mutation representatives; others have small support and weak discovery margins. In particular the six-change NL_11 candidate has strong model score but sparse linked support, and clustered/selected SHM beyond the short-range null remains a possible explanation. Strong model evidence does not exclude a misspecified biological null. No candidate was manually removed because it was outside the dropout panel.

## Reproduction

`benchmarks/personalized-germline/v0388/` contains exact selected AIRR representatives for all six arms, prepared references, panels, configurations, final dashboards/evidence/FASTA, scores and hashes. From the source root:

```
python benchmarks/personalized-germline/v0388/replay.py --animal ERR4238110 --out /tmp/swig-replay-8110
python benchmarks/personalized-germline/v0388/replay.py --animal ERR4238104 --out /tmp/swig-replay-8104
```

Replay checks active novel sequences against recorded results. These compact inputs reproduce discovery, not raw assignment or denoising. For a complete raw-read rerun, use `python benchmarks/personalized-germline/v0388/run-recorded-panels.py --animal ERR4238110 --fasta INPUT.fasta.gz --out WORKDIR` (then the second animal separately). This verifies the sample hash, resolves reference/configuration paths and runs all three complete pipelines using the fixed recorded panels. Full original FASTA data are not duplicated in the package.

Synthetic generator: `benchmarks/personalized-germline/stress.mts`; results in `v0388/synthetic`. It deliberately uses uniform background mutation rather than the fitted HS5F model, exceptional hotspots, adjacent coordinated mutations, sequencing error and known rare alleles. The tested null has no false novel call, and both a 0.5% single-site mixed-repertoire allele and a 1% linked high-SHM allele are recovered. A 1% single-site high-SHM allele is **missed**: limited low-exposure evidence cannot reliably distinguish it from recurrent mutation. These few deterministic cases diagnose mechanisms, not population sensitivity or specificity.

## Verification and biological references

All 271 tests pass, including synthetic correlated-SHM rejection, finite high-SHM rare-allele inference, terminal ambiguity, subject isolation and CLI output. TypeScript and production build pass. Existing AER-R/R-optimized defaults and multi-file assignment progress corrections are retained.

- [IGoR: rearrangement inference](https://www.nature.com/articles/s41467-018-02832-w): gene-dependent deletions and junction sequence modeling.
- [Tandem substitutions in SHM](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.807015/full): adjacent changes need not be independent mutation events.
- [Correlations among SHM sites](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.618409/full): context overlap and mutation burden can generate correlated substitutions.
