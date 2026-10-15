# 8. Inherited variation and germline references

Two people can inherit different versions of an immunoglobulin or T-cell receptor segment. Even within one person, the two copies of a chromosome can carry different alleles. Receptor loci also contain closely related genes and can vary in gene copy number. A public reference catalogue brings together sequences found in many individuals; each person's expressed repertoire draws on a subset of that variation. ([R16](references.md#r16))

This matters whenever observed receptors are compared with a reference. Suppose a donor carries a V allele with `G` at a particular position, while the chosen reference has `A`. Receptors assembled from that allele will tend to share the `G`. Counting every such difference as somatic hypermutation would give those receptors an inherited “mutation.”

## Inherited differences and acquired mutations leave different patterns

An inherited V allele can be used in many independently rearranged B cells. Its characteristic bases can therefore occur beside many different CDR3s and J segments. A mutation acquired by one activated B cell is passed to its descendants and initially belongs to that one clonal family.

The distinction becomes harder when SHM repeatedly targets the same position. Mutation hotspots and selection can produce similar changes in independently rearranged families. Conversely, a rare inherited allele may appear in only a few sampled families. The interpretation depends on the combination of sequence pattern, mutation burden, independent rearrangements, and reference coverage. ([R11](references.md#r11), [R16](references.md#r16))

Consider two observations. In the first, one large family contributes 5,000 reads carrying the same change. In the second, 30 families with diverse junctions carry that change, including several otherwise close to their assigned germline sequence. The second observation provides more independent evidence for an inherited difference. The first can be explained by a single mutation followed by expansion.

## Combining ambiguous allele assignments across a repertoire

An individual read may cover too few distinguishing positions to choose between alleles A and B. Other reads from the same donor can help: some may cover those positions clearly, and many may jointly favour one expressed allele over the other.

Swig's **repertoire allele pooling** combines per-record candidate evidence to estimate expressed reference usage. Its Dirichlet-mixture option fits regularized usage proportions. Its hurdle option separately models whether an allele is detectably active and how frequently an active allele is used. The latter can assign zero fitted usage to an inactive candidate. These estimates can refine ambiguous calls, with confidence-gated application available. ([S14](references.md#s14))

A [mixture model](https://en.wikipedia.org/wiki/Mixture_model) treats the observations as coming from several possible sources. Here the sources are reference alleles. An ambiguous sequence contributes evidence to several sources; the overall fit uses the other observations to help distribute that evidence. The hurdle adds a source-presence decision before estimating usage.

The biological quantity being fitted is **expression in the sampled repertoire**. It depends on inherited availability, rearrangement preferences, selection, cell abundance, and the assay. An allele absent from the sample could be genetically absent or simply poorly represented. Genomic sequencing provides more direct evidence about inherited sequence and copy number.

Pooling normally gives each active record one vote; duplicate-count weighting changes that balance. Neither unit automatically equals an independent rearrangement. A heavily sampled B-cell family can contribute many distinct, mutated records. Pooling samples from the same donor can improve coverage, whereas pooling donors mixes potentially different inherited sets. The current implementation primarily targets V and J calls; D pooling is experimental. ([S14](references.md#s14))

## Looking for a missing V reference

Sometimes the best explanation is an allele absent from the supplied catalogue. A recurring set of linked differences is especially informative: the same bases occur together within sequences from several independent families.

Swig's **missing-V warning** starts with a low-SHM representative from each donor–lineage group and looks for these recurring patterns. It then checks the other members for contradictions, such as the original reference base at a proposed inherited site or incompatible assignments. Junction diversity, coverage, and hotspot context help assess the evidence. ([S09](references.md#s09))

For example, changes at positions 40 and 75 might occur together in many families. Finding a well-covered, low-error family member with the original base at position 40 weakens the claim that this family's starting allele carried both changes. Reversion remains a biological possibility, but the simple inherited explanation has become less complete.

A warning is a starting point for reviewing candidate references. Useful checks include raw sequence quality, coverage of every proposed position, support from independent low-mutation families, and consistency across samples from the donor. Independent genomic evidence is especially valuable for establishing a new inherited allele.

## Inferring a personalized expressed V set

The **personalized expressed V-set** analysis takes the next step: it compares candidate sets of known and proposed V sequences that could explain a donor's observed repertoire.

The model gives one usable low-SHM observation per lineage a role in the comparison, reducing the influence of clone expansion. It considers recurrent linked substitution candidates and models expected SHM. Its current search concerns same-length, substitution-compatible sequences; suitably aligned references can compete across gene labels. Candidate complexity is penalized so that adding more alleles must improve the explanation sufficiently. ([S09](references.md#s09))

The principle resembles [model selection](https://en.wikipedia.org/wiki/Model_selection). One explanation might use a single inherited V sequence plus acquired mutations. Another might use two inherited sequences, requiring fewer acquired changes but introducing an extra inherited candidate. The comparison balances fit against complexity.

The resulting set concerns the V sequences expressed and distinguishable in the available data. Identical covered sequence can leave gene labels unresolved. Uncovered or weakly informative positions can retain reference assumptions. Low-SHM representative selection also depends on the current reference assignment, which is one reason to review the result after reannotation.

Swig can export a personalized reference set from this analysis. Applying it changes the starting assumptions for assignment and SHM counting, so the appropriate next step is to rerun annotation and examine the consequences. ([S09](references.md#s09))

## A competing somatic explanation

A separate experimental **joint inherited/SHM model** asks how much evidence for inherited variation remains when recurring somatic outcomes are given a flexible explanation. Mutation and selection can repeatedly produce some sequence states, particularly in highly mutated repertoires.

Swig separates candidate discovery from evaluation using lineage-level training and held-out groups. The held-out observations test whether proposed inherited candidates explain additional data better than the competing somatic model. Its reported evidence is conditional on those models and the available observations. The analysis has a separate JSON result and currently does not provide a confident full-length personalized FASTA export. ([S15](references.md#s15))

A [held-out set](https://en.wikipedia.org/wiki/Training,_validation,_and_test_data_sets) contains observations excluded from the relevant fitting or discovery step. It helps reveal candidates that mainly describe peculiarities of the discovery data. Splitting by lineage is important because relatives share mutations; distributing close relatives between training and evaluation would make the test less independent.

Some datasets contain little information that can separate inherited from repeatedly acquired sequence. A highly mutated repertoire without suitable independent near-germline observations is an example. In that situation, preserving competing explanations is a useful result in itself.

## Consequences for downstream biology

Reference changes can alter V labels, measured SHM, candidate lineage boundaries, and reconstructed ancestors. A correction from allele A to allele B may remove several differences previously counted as mutations. It can also change the plausible pre-SHM sequence at a binding-site residue.

For a study that uses reconstructed antibodies experimentally, retain the original reference set, the personalized set, the evidence for proposed changes, and the results under both defensible choices. This connects a reference decision to its actual biological consequences.

---

[Previous: Errors, duplicates, and chimeras](07-errors-and-chimeras.md) · [Contents](index.md) · [Next: Repertoires and clonal families](09-repertoires-and-clones.md)
