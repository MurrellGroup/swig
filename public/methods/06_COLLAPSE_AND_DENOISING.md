# Exact collapse and read-error denoising

All four methods first respect a hard study boundary (`dataset`, `sample`, `subject`, `cohort`, or complete study). Different normalized top C-gene/isotype calls are separate partitions by default. Existing `duplicate_count` values are treated as multiplicity and summed into the retained representative.

If allele refinement has been applied, any V/J-dependent key uses its policy-selected overlay. The immutable AIRR calls remain recoverable.

## Method A: exact collapse

Method A groups records with byte-identical selected keys. Choices are full input sequence, V–J-trimmed nucleotide sequence (default), locus + CDR3 nucleotide, or locus + normalized V/J calls + CDR3. No mismatches are permitted and no sequencing-error model is invoked. The earliest ordinal is the deterministic representative. Records missing a required key are excluded from the representative set by default or retained unchanged by explicit policy.

## Shared denoising preparation

Methods B–D require a normalized V–J-trimmed nucleotide sequence and a locus/V/J partition at the selected call resolution and ambiguity policy. Sequences are exact-dereplicated into a two-bit packed arena. Only one V/J/C partition's temporary candidate index exists at a time. Records containing ambiguous bases follow the selected exclude/retain policy; records without usable trim or V/J follow the separate unresolved policy.

## Method B: FAD-compatible

For each unique sequence, Swig counts all overlapping 6-mers and uses squared Euclidean distance between 6-mer count vectors. A configured FAD radius \(t\) is represented by `distance² ≤ floor(12t)`. Variants are considered in descending abundance, and only variants meeting the minimum parent abundance enter template selection.

Method 1 uses the published abundance/neighbor rule. Method 2 additionally computes

\[
\lambda=\frac{n_{parent}}{f_0}\epsilon,
\]

where \(\epsilon\) is the configured error rate and \(f_0\) the expected zero-error fraction. The strict Poisson upper tail for the child count is multiplied by sequence length and compared with \(\alpha\), matching the source implementation's decision direction. After template selection, every non-template variant is assigned to its globally nearest accepted corrected-6-mer centroid within the V/J partition, even when outside the template-selection radius. An exact VP tree accelerates this final nearest-centroid search without changing its result.

This is labeled **FAD-compatible** because it reproduces FAD's corrected-6-mer distance, template decision, and forced nearest-centroid semantics while replacing dense scans with exact indexes and adding immunoreceptor V/J/C/study partitions. It is not the Julia executable and is not RAD. Source: Kumar et al., [Long-read amplicon denoising](https://academic.oup.com/nar/article/47/18/e104/5550323).

## Method C: conservative exact-neighbour model

Method C is **custom/experimental**. A sequence split into \(d+1\) exact blocks must share at least one block with any equal-length sequence within Hamming distance \(d\). The block index generates a complete bounded-radius candidate set, followed by exact Hamming verification.

For a candidate parent with count \(n_p\), child length \(L\), and \(d\) substitutions,

\[
p_{exact}=(\epsilon/3)^d(1-\epsilon)^{L-d},\qquad
\lambda=n_pp_{exact}.
\]

The strict Poisson upper tail for the observed child count is multiplied by \({L\choose d}3^d\). The child collapses only when that adjusted value is at least \(\alpha\), i.e. its abundance remains compatible with the sequencing-error explanation. Isolated singletons remain templates and no distant-centroid assignment occurs.

## Method D: indel-aware bounded-edit model

Method D is **custom/experimental**. A complete length-aware \(d+1\)-segment join generates every candidate within configured Levenshtein radius 1 or 2. An allocation-bounded banded dynamic program returns exact distance and substitution/insertion/deletion counts; equal-cost paths prefer fewer indels.

Substitution-only paths use Method C's Poisson rule. A path containing an indel instead collapses when the parent:child abundance ratio reaches the configured minimum (2 by default). The best eligible parent is chosen by distance, fewer substitutions, greater abundance, then ordinal. No distant-centroid assignment is performed.

## Defaults and diagnostics

The shipped error-rate starting value is 0.00473, inherited from the linked MiSeq workflow and not a universal rate. Other defaults are \(\alpha=0.01\), minimum parent count 2, FAD radius 1/method 2, conservative Hamming radius 1, edit radius 2, indel parent:child ratio 2, and candidate cap 50,000 per variant.

The result reports unique representatives, collapsed rows, preserved abundance, partitions, exact candidate comparisons, substitution/indel merges, unresolved/ambiguous exclusions, and cap hits. Hitting the cap makes candidate generation incomplete and produces a warning; increase it and rerun before treating the result as final.

## Parallel execution and reproducibility

Exact collapse assigns complete comparison keys to deterministic hash shards. Identical keys therefore always enter the same shard, and each shard preserves input ordinal order. Swig merges shard results by ordinal, so the representative, abundance, and output are identical regardless of worker completion order. Browser and CLI enable these shards for at least 10,000 records, where worker startup is normally amortized.

Methods B–D run independent study/C/locus/V/J partitions on a dynamically scheduled worker pool once the run contains at least 500 unique variants. A partition is never subdivided: template/centroid decisions can depend on every variant within it, particularly FAD's globally nearest accepted centroid. Consequently a repertoire containing one overwhelmingly large partition cannot usefully occupy every configured worker without changing the method. The serial and parallel routes call the same pure partition kernel and use identical sort keys and ordinal tie breaks. Cancellation terminates the outer worker and all of its active partition workers; results are committed only after every partition succeeds.

## Limitations

- FAD was developed for amplicon denoising; its forced nearest-centroid assignment can merge a genuinely distinct variant.
- Methods C/D use a simple independent, equal-substitution error model rather than per-base qualities or learned context.
- Method D's indel abundance rule can merge true biological length variants.
- Denoising expanded clones or PCR-amplified counts confounds molecule abundance with biological abundance unless preprocessing supplies suitable multiplicities.
- Parallel speedup is bounded by the largest independent V/J partition and worker startup/serialization costs.
