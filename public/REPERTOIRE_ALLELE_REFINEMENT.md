# Repertoire-level germline-call refinement

## Purpose

V(D)J assignment is normally performed one sequence at a time. With a deliberately overcomplete reference database, a read may not distinguish several related alleles even when the repertoire collectively supports only a subset. Swig's optional allele-pooling step combines:

1. local alignment evidence retained for each read; and
2. the allele-usage distribution inferred from other active records in the same repertoire pool.

The method refines assignments. It does **not** certify a diploid genotype, infer genomic copy number, discover an allele absent from the reference set, or overwrite the original AIRR result.

## Independent pools

Models are fitted independently by user-selected study boundary (donor/subject by default), receptor locus, and segment. V and J are enabled by default. D is available as an experimental option because short templated spans, trimming, and N addition make D emissions much less identifiable.

No posterior information crosses these boundaries. Donor pooling is usually appropriate for longitudinal or compartmental data; sample pooling is available when samples must remain independent.

## Reference graph and local evidence

Reference records with identical nucleotide sequences are collapsed into one unresolved node carrying all equivalent labels. The remaining nodes are connected within a small bounded nucleotide edit radius. This graph constructs sparse candidate support; Swig never assigns every read nonzero evidence against the complete database.

Let \(E_{ra}\) denote local evidence for read \(r\) and reference node \(a\).

- Every literal co-optimal AIRR call starts at weight 1. A two-way explicit tie is therefore 50:50 before repertoire pooling.
- A retained alternative with alignment score \(S_a\) receives

  \[
  \min\left(1,\exp\left[\frac{S_a-S_{\max}}{T}\right]\right),
  \]

  where \(T\) is an adjustable score temperature.
- An alternative without a numeric score receives a separate fallback weight.

These values form an evidence kernel. SwiftIG alignment scores are not claimed to be fully calibrated read likelihoods.

### SHM-adaptive neighbour leakage

An unreported nearby allele should not receive an artificial exact zero merely because it fell outside the caller's retained hypothesis set. For a same-length one-substitution neighbour, Swig uses per-edit evidence odds

\[
o_1(\mu)=\min\left(o_{\max},\;o_0+s\frac{\mu}{3(1-\mu)}\right).
\]

- \(o_0\): zero-SHM neighbour odds (default 0.01). This is irreducible assignment-model leakage, **not** a sequencing-error estimate.
- \(\mu\): read-level SHM estimate, normally one minus the best retained V-reference identity, clamped at a configurable maximum.
- \(1/3\): equal-substitution probability of landing on the particular nucleotide distinguishing the neighbouring allele.
- \(s\): SHM-sensitivity multiplier (default 1; zero disables SHM adaptation).
- \(o_{\max}\): conservative cap for highly mutated or poorly matched reads.

The term \(\mu/[3(1-\mu)]\) is the likelihood ratio of a particular diagnostic substitution to retention under a context-averaged equal-substitution approximation. It is intentionally simple and inspectable, not a five-mer or other context-dependent model.

A substitution-only neighbour at distance \(d\) receives relative leakage \(o_1(\mu)^d\). Reference neighbours involving an indel receive the non-SHM baseline only because the substitution derivation does not apply. Direct, alternative, and neighbour weights are combined, capped to a maximum sparse row width, and normalized.

### Interactive evidence-kernel inspection

The advanced settings include a reference-sequence inspector so parameter changes can be interpreted before fitting the repertoire model. The user chooses a V, D, or J database allele and an assumed read SHM fraction. Swig then evaluates the exact sparse evidence-row constructor used by the worker with that allele as the sole literal primary call.

- The selected reference sequence remains visible as the alignment anchor, but its dominant probability is deliberately not drawn as a bar.
- Every non-primary candidate receives a horizontal bar and its absolute probability after normalization over the complete row, including the primary.
- Alternative bar lengths are scaled relative to the largest non-primary candidate. The printed percentage, rather than bar length, is the absolute probability.
- Candidate sequences are globally aligned to the selected reference through a common anchor coordinate. Columns that are gaps in every displayed sequence are removed.
- Nucleotide coloring shows ordinary bases. Difference-highlighter coloring makes matches neutral, substitutions red, and insertion/deletion columns amber.
- Sequence-identical database labels remain one unresolved reference node and are displayed together.

The diagnostic responds to the zero-SHM floor, SHM sensitivity and cap, assumed SHM, edit radius, maximum neighbour odds, and candidate cap. Dirichlet alpha and variational stopping controls do not alter this local pre-repertoire kernel. Numeric settings commit on Enter or when focus leaves the field, matching the rest of Swig's large-data controls; the inexpensive diagnostic then updates immediately.

## Dirichlet mixture

For one independent pool, let \(\boldsymbol\theta\) be expressed reference usage:

\[
\boldsymbol\theta \sim \operatorname{Dirichlet}(\boldsymbol\alpha),
\qquad
p(z_r=a\mid\boldsymbol\theta,E_r)\propto\theta_aE_{ra}.
\]

Mean-field coordinate ascent uses

\[
q(z_r=a) \propto E_{ra}\exp\left(\psi(\gamma_a)-\psi\left(\sum_b\gamma_b\right)\right)
\]

and

\[
\gamma_a=\alpha_a+\sum_rw_rq(z_r=a).
\]

The default \(w_r=1\) gives every active record one vote. Optional abundance weighting uses `duplicate_count`, but allows PCR abundance and clonal expansion to influence the inferred mixture.

The symmetric prior covers every locus-matched reference node in the selected segment, including an allele absent from every read's sparse candidate neighbourhood. Such an inactive node has zero row likelihood, receives no assignment responsibility, and retains \(\gamma_a=\alpha_a\). Swig accounts for the combined prior mass of these inactive nodes in the Dirichlet normalization without materializing them once per donor or once per read. Nodes with unknown locus metadata are included only when they occur in an observed candidate row, because there is no defensible locus pool to which an otherwise inactive custom record can be assigned.

This resembles sparse variational machinery used for LDA, but the biological model is a finite mixture with one latent germline node per read, not a document containing multiple topics.

### Before/after assignment-frequency figure

After fitting, Swig displays paired horizontal bars for every modeled allele in a selected donor/study-boundary, locus, and segment pool. Rows are sorted by post-pooling assignment frequency.

- **Before** is the normalized sum of local evidence responsibilities, \(\sum_r w_r E_{ra}\), over materialized candidates.
- **After** is the normalized sum of fitted variational responsibilities, \(\sum_r w_r q(z_r=a)\).

Both denominators describe assignment mass over observed records. Dirichlet prior-only mass is intentionally excluded, so the figure isolates redistribution of assignments rather than conflating it with the usage prior. The model table still reports the full posterior usage mean, which includes the complete reference-set prior denominator. The interactive display has a configurable row limit; CSV export always includes every modeled allele in the selected pool, and the complete paired-bar figure is exportable as SVG.

## Sparse browser implementation

- Candidate evidence uses compressed sparse row storage.
- Repeated sparse evidence patterns are collapsed during coordinate updates and multiplied by their accumulated weight.
- Donor/locus/segment pools are fitted separately.
- Expected-log weights use a stable digamma approximation and responsibilities are normalized in log space.
- The model runs in a dedicated worker. V, D, and J matrices are built and released sequentially to bound peak memory.

For \(N\) modeled records and average sparse width \(K\), iterative work is proportional to \(NK\), rather than \(N\) times the complete reference count.

## Applying and exporting results

Fitting never changes downstream calls. The user must explicitly choose **Apply to downstream calls** and a minimum posterior threshold.

- A V or J MAP call at or above the threshold replaces that call for subsequent lineage assignment, query constraints, lineage alignment construction, SHM stratification, and missing-allele screening.
- A call below threshold remains as originally reported.
- Reset restores immutable original calls.
- Existing lineage-dependent results are invalidated when the overlay changes.

Exports include a complete model summary, a long-form sparse per-record posterior sidecar, a refined AIRR table, paired-frequency chart data as CSV, and the paired-frequency figure as SVG. The sidecar reconstructs every nonzero candidate responsibility from the saved mixture parameters while streaming the AIRR input, so the complete reads-by-candidate posterior is available without retaining a second large matrix in interactive memory. The refined table places threshold-passing calls in the ordinary call columns while `swig_original_*` and `swig_repertoire_*` columns retain provenance. Options, compact MAP/entropy vectors, mixture summaries, threshold, and apply/reset state are included in Swig sessions.

## Interpretation and limitations

- Posterior usage is expressed rearrangement usage, not genomic genotype frequency or copy number.
- Expanded lineages violate independent-read assumptions. Unique weighting is the conservative default.
- A missing true allele cannot be recovered unless represented by a modeled reference node. Use Swig's missing-allele screen and a dedicated workflow such as IgDiscover, TIgGER, partis, or related methods.
- D pooling is experimental.
- The SHM term is a context-averaged approximation. Its sensitivity control exists for explicit sensitivity analysis, not to imply universal calibration.
- A symmetric Dirichlet prior does not force exact zero usage: every known-locus database node retains its prior mass. Local assignment probability is nevertheless exactly zero outside the configured sparse evidence neighbourhood; broadening that neighbourhood is an evidence-model decision, not something the repertoire prior should silently do.

## Related methods

This narrower assignment-refinement model is informed by the broader subject-specific germline/genotype problem:

- [TIgGER genotype inference](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.00129/full)
- [partis per-sample germline inference](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007133)
- [PIgLET allele sequence clusters](https://academic.oup.com/nar/article/51/16/e86/7238142)
