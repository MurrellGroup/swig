# Joint inherited / SHM model (experimental, 0.38.9)

This is a separate analysis box. It does not replace the 0.38.8 personalized-reference method, its results, or its FASTA export. Its output is evidence for sequence hypotheses, not an automatically certified genotype.

## Biological question and identifiable information

A recurrent sequence can be inherited from a germline allele, generated repeatedly by SHM, enriched by selection, or misrepresented by assignment/rearrangement uncertainty. Linked substitutions are not sufficient evidence of inheritance: somatic outcomes can also be linked. The distinction requires information about their occurrence across mutation exposure, especially genuinely unmutated lineages, and competition with other reference/haplotype explanations.

Without constraints on mutation and selection, this problem is not identifiable. For example, if all observations are mutated, a somatic endpoint distribution may reproduce the same sequences and frequencies as an additional expressed allele. No score can resolve that ambiguity from these observations alone. The new model explicitly includes this competing explanation and can abstain.

## A normalized model

Condition on the reference-defined competitive component and the observed coverage structure. Let c index candidate V sequences, t index latent mutation exposure, theta be inherited sequence frequencies, and q be a nuisance distribution over somatic endpoint sequences. Both theta and q are probability vectors. With r_lambda(t) = exp(-lambda t),

```
p(x | theta, q, lambda) = sum_t w_t sum_c [r_lambda(t) theta_c + (1-r_lambda(t)) q_c] K(x | c,t).
```

At t=0 the somatic branch has zero mass. It can nevertheless explain an apparently near-germline observation through a small positive latent exposure. Exposure is integrated, not set to zero merely because no other mutation was observed.

The somatic endpoint distribution q is deliberately broad: a proposed multi-site sequence can occur under the null too. This is a phenomenological nuisance model for recurrence and selection, not a claim that arbitrary multi-site replacement is a literal biochemical event. Its breadth improves robustness at a real cost in power. It is not a detailed lineage-tree or repair-tract simulator.

K is the same observation kernel in every comparison. At each V position it uses a four-state CTMC with HS5F rates, allowing repeated hits and reversions, followed by a separate symmetric sequencing-error channel. The focal nucleotide changes state; flanking context is frozen to the candidate sequence. This approximation does not model evolving motif interactions. Raw terminal sequence is marginalized over V deletion lengths and junction bases. Rearrangement distributions and the zero-plus-exponential exposure distribution are calibrated on training observations; all physical genes in the tested competitive component are excluded from its calibration lookup. The final two template bases are unresolved and marginalized as unknown rather than inferred.

The exposure integral uses an atom at zero and 16 equal-probability exponential-quantile states. Lambda belongs to the declared finite family {0, 1, 10, 100, infinity}; infinity still leaves r(0)=1. Numerical certificates apply to this finite family, not to an unsearched continuum of hazard functions.

For fixed lambda, integrating K gives arrays A and B and

```
log L = sum_i log(sum_c theta_c A_ic + q_c B_ic).
```

This is concave on a product of simplexes. EM fits the weights. The maximum gradient minus its current weighted average in each simplex gives a tangent-plane upper bound on improvement. Consequently an unfinished null fit cannot manufacture a large evidence score by failing to find its optimum: its upper bound, rather than its achieved likelihood, is used in the denominator.

## Discovery and evaluation are separated

Lineage IDs are deterministically hashed into equal-probability training/test halves, before proposal generation or parameter calibration. Existing lineage representatives and upstream assignments are reused. Proposals are generated only from training lineages. Reference-defined, explicitly alignment-compatible connected components permit nearby genes to compete. Search remains bounded (64 novel proposals by default), so failing to nominate an allele is not evidence against it.

For each candidate c, the null sets theta_c=0 but leaves q_c free. Other candidate sequences remain available. The alternative predictor is trained without the test half. It averages two normalized predictors with equal fixed weights: the broad fitted candidate model and a focused model permitting known inherited sequences plus c. Averaging likelihoods is different from choosing whichever test score is larger; the mixture weight is included. This reduces some training overfit without changing the test null.

The reported log evidence is

```
log(training-fitted predictive likelihood on test data)
    - upper_bound(maximum test-data null log likelihood).
```

This follows the split-likelihood principle of [Universal inference](https://www.pnas.org/doi/10.1073/pnas.1922664117). A candidate provably unable to achieve positive evidence against even a feasible null receives zero evidence (log evidence minus infinity); exports encode this as null. No later whole-sequence selection score overrides this comparison. No BIC approximation or copy-number cap is used in the new method.

The displayed model-conditional evidence rule is E >= M/0.05, where M is the number of training-nominated hypotheses tested for that subject. If the null family contains the true conditional distribution and test observations are independent, the split likelihood and union bound provide a familywise statement. **Those assumptions have not been established for real repertoire data. This is not a demonstrated 5% real-world false-positive rate or FDR.** Estimated mutation/rearrangement nuisance distributions, assignment-dependent coverage, clone errors and representative selection can violate them.

## What remains unresolved

The lowest-current-SHM representative is an order statistic, not a randomly sampled lineage tip; this can distort its exposure distribution. Somatic endpoint weights are flexible but not a complete selection model. Training-calibrated nuisance parameters are fixed rather than integrated with full uncertainty. References and covered coordinates are conditioned upon, and errors in these inputs remain relevant. Terminal junction composition is pooled, not freely profiled per gene in this model. The whole likelihood is consistent, but its biological assumptions remain approximations requiring validation.

In particular, removing spurious calls by also losing rare alleles is not an acceptable demonstration of improved discovery. The release benchmark reports recovery and failures separately. The experimental box is provided for comparison, not promoted as superior to the existing box.

## Running either or both

The two Web boxes have independent results and saved-session fields; running one does not clear the other. Both are lazy workers and add no idle inference work. The new box exports its complete evidence JSON; it does not export a confident full-length personalized FASTA.

```
swig-cli personalized-germline --airr processed.airr.tsv.gz --v-reference V.fasta --out existing-method
swig-cli joint-germline --airr processed.airr.tsv.gz --v-reference V.fasta --out joint-method
```

The new CLI accepts `--max-candidates N` and `--iterations N`. Larger iteration budgets tighten the numerical null bound; they do not alter the evidence rule. Sequence and V/D/J boundary fields should be retained in AIRR input. Subjects are analyzed independently.
