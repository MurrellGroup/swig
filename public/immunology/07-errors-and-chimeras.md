# 7. Errors, duplicates, and chimeras

Two records can resemble one another because they came from the same molecule, from different cells in the same clone, or from independently assembled receptors. Laboratory errors add further possibilities. Preparing a repertoire for biological interpretation requires keeping these explanations separate.

## Exact sequence collapse

**Exact collapse**, also called [deduplication](https://en.wikipedia.org/wiki/Data_deduplication) or dereplication in sequence workflows, replaces records with the same selected key by a representative plus a count.

The key matters. Full-read equality treats differences in constant tails or residual flanking sequence as distinct. Equality over the V–J-aligned region focuses on the rearranged variable region. Equality of CDR3 alone deliberately combines records that may differ elsewhere.

Consider two B-cell sequences with identical CDR3s but different mutations in CDR1. Full-variable-region collapse preserves both variants. CDR3-only collapse groups them under one representative, removing that distinction from the representative sequence set. Which choice is useful depends on the question.

Swig offers full-sequence, V–J-region, CDR3, and V/J-plus-CDR3 keys. Collapse sums multiplicities and respects the selected study boundary. Distinct called constant genes or isotypes are kept separate by default, preserving an otherwise easily lost distinction between class-switched members. ([S06](references.md#s06))

## Denoising proposes an error explanation

A sequencing error can turn an abundant sequence into a nearby low-count variant. If a sequence has thousands of observations, occasional one-base errors can generate many apparent neighbours. **Denoising** uses sequence differences and an error model to decide which neighbours can plausibly arise this way.

An illustrative dataset might contain 2,000 copies of sequence A and three copies of sequence B differing at one position. An error explanation is plausible. The same observations could also contain a genuine rare B-cell variant. Sequence similarity and abundance provide evidence, while the assay and independent replication determine how persuasive that evidence is.

[Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) counts mismatching positions in equal-length sequences. [Edit distance](https://en.wikipedia.org/wiki/Levenshtein_distance) also allows insertions and deletions. These distances describe proximity; a denoising model supplies the rule for turning proximity into a merge.

Swig's choices include a FAD-compatible abundance-based method, a conservative bounded-substitution method, and an indel-aware method. The latter two are experimental. FAD-compatible assignment can place a variant with its nearest accepted template even outside the local template-selection radius; the conservative choices restrict merging to qualifying nearby candidates. Each choice can affect rare biological variants. ([S06](references.md#s06), [R25](references.md#r25))

For SHM studies, the main trade-off is direct. Aggressive correction can erase genuine mutation and short indel variants. Minimal correction can retain technical errors as extra mutations, tiny clones, or terminal branches. Results of particular interest should be checked under more than one defensible correction setting, using controls or upstream UMI consensus where available.

## PCR can join parts of different templates

During PCR, an incompletely extended product can anneal to a related template and continue being copied. The resulting molecule contains sequence from more than one source. Such a product is called a [PCR chimera](https://en.wikipedia.org/wiki/Chimera_(molecular_biology)).

For example, the first part of an apparent V region might resemble allele A while the later part resembles allele B. Treating the entire sequence as one ordinary receptor can create a false pattern of mutations or a misleading germline assignment. A chimera can also bridge otherwise separate sequence groups.

The real V–D–J assembly described in Chapter 2 happens in a developing lymphocyte. A PCR chimera arises during library amplification. Its evidence is an unexpected change of template within the amplified material, such as within a V segment. ([R26](references.md#r26))

## Reading a chimera result

Swig's **CHMMAIRRa** analysis compares an observed V or J sequence with an alignment of references. A [hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) represents the possibility that different positions were copied from different references. Its hidden state describes the current reference source; the observed bases provide evidence for or against that source.

The model produces a posterior probability of a reference switch and can reconstruct a likely path and breakpoint. The inference depends on the available references and on how ordinary mutations are modelled. Missing parental alleles or inherited variants can change the result. The detailed alignment helps distinguish a well-supported switch from a less decisive pattern. ([S07](references.md#s07))

A high-scoring candidate is worth reviewing in the context of its abundance, quality, reference matches, and related sequences. Independent libraries can help assess whether the pattern is reproducible. Some records remain unevaluated because the required reference or alignment information is unavailable; that status deserves its own count.

## Filters define the population being studied

A repertoire may be restricted to productive heavy chains, one isotype, a range of CDR3 lengths, or a particular sample. Such choices change the population to which every downstream summary applies.

Filtering on mutation burden, for example, changes the subsequent mutation distribution by construction. Excluding long CDR3s changes the available evidence for long-junction or Double-D analyses. A useful report records the number of records and total multiplicity retained at each stage, together with the reasons for exclusion.

Swig maintains a cumulative **working set** for downstream analysis while retaining the original assigned rows. Applying collapse, chimera exclusion, or repertoire selection changes that working set. Saving its settings makes a result interpretable and allows the effect of a particular decision to be examined later. ([S05](references.md#s05), [S12](references.md#s12))

---

[Previous: V(D)J assignment](06-vdj-assignment.md) · [Contents](index.md) · [Next: Inherited variation](08-germline-variation.md)
