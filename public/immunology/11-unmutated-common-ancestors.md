# 11. Estimating an unmutated common ancestor

An **unmutated common ancestor**, or [UCA](https://en.wikipedia.org/wiki/Ancestral_sequence_reconstruction), is the inferred rearranged receptor sequence from which a B-cell family developed before somatic hypermutation. Reconstructing it asks what the family's receptor looked like at the beginning of its mutation history.

The target already contains the consequences of V(D)J recombination: chosen segments, trimmed ends, and junctional additions. “Unmutated” refers to the subsequent SHM process. It does not remove the diversity created when the receptor was assembled. ([R29](references.md#r29))

An antibody has both a heavy and a light chain. Estimating the ancestral heavy-chain variable region reconstructs one component. A functional ancestral antibody requires an appropriate ancestral light chain and a chosen experimental expression format.

## Three different ancestral ideas

An **inherited V allele** is one segment present before receptor rearrangement. It supplies much of the variable domain but lacks the complete rearranged junction.

The **most recent common ancestor**, or [MRCA](https://en.wikipedia.org/wiki/Most_recent_common_ancestor), of the sampled sequences is their most recent shared ancestor. It can already carry mutations acquired before the sampled branches diverged.

The **UCA** is the earlier, pre-SHM rearranged sequence. If a mutation spread through the lineage before the ancestors of all sampled cells separated, that mutation belongs to the sampled MRCA but differs from the UCA. Figure 3 in [Chapter 10](10-shm-alignments-and-trees.md#what-a-phylogeny-represents) shows this situation.

Sampling can leave a substantial gap between these ancestral stages. A dataset collected late in a response may contain many descendants of an already-mutated sublineage and no surviving record of the earlier states.

## Why a consensus can retain acquired mutations

A [consensus sequence](https://en.wikipedia.org/wiki/Consensus_sequence) selects a representative base at each aligned position, often the most frequent. If every observed family member carries a shared acquired mutation, the consensus will carry it too.

Choosing the least-mutated observed member has a similar difficulty. That member can still contain shared early mutations and its own private changes. Replacing its V and J differences with reference bases helps in confidently assigned templated regions but leaves uncertain junction sequence and reference alternatives to resolve.

The useful question is therefore: which complete starting sequences plausibly arose by recombination and could have produced the related descendants? Statistical ancestor reconstruction combines these two sources of information. ([R18](references.md#r18), [R29](references.md#r29))

## Recombination supplies information about plausible starting sequences

A candidate UCA usually resembles a V reference over a long tract and a J reference near the end. Between them, the sequence may include retained D sequence and bases added during joining. Several trimming boundaries and D assignments can explain the same short stretch.

A **recombination model** assigns probabilities to these alternatives. It describes preferences for segment use, end trimming, and added sequence. Such probabilities are called a [prior](https://en.wikipedia.org/wiki/Prior_probability) when they are used before incorporating the particular descendant family being analysed. The prior can favour biologically plausible arrangements while allowing uncertainty about which one occurred. ([R31](references.md#r31))

For example, a junction tract might be explained by a short retained D segment or by added nucleotides. The alternatives can agree on the actual ancestral bases while disagreeing about their source. Conversely, different source assignments can imply different ancestral bases after allowing for later mutations.

Swig represents recombination with a [hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model). Hidden states track V, D, J, or non-templated contributions along the alignment, including uncertain boundaries and optional additional D contributions. An N-addition state generates concrete nucleotide possibilities; the literal sequence symbol `N` records unresolved nucleotide identity. ([S13](references.md#s13))

A zero-retained-D explanation can be relevant even to a chain normally assembled through D-containing recombination: extensive trimming can leave little recognizable D sequence. The number of supported D tracts in a reconstruction and the number of developmental DNA-joining events are related but distinct quantities.

## The tree supplies information about shared mutations

For each proposed starting sequence, an evolutionary model calculates how compatible the descendants are with that start. This quantity is the [likelihood](https://en.wikipedia.org/wiki/Likelihood_function).

The tree accounts for shared history. Ten sequences carrying the same mutation can descend from one cell in which that mutation occurred. Treating them as ten independently mutated descendants would misrepresent the evidence. A tree-based likelihood distributes changes along ancestral and descendant branches. ([R27](references.md#r27), [R18](references.md#r18))

Swig conditions on a fixed observed-sequence tree, using at least three observed rows and excluding the synthetic germline guide. It estimates where the UCA connects to that tree and the length of an ancestral branch that can accommodate shared early changes. Its nucleotide likelihood is site-independent and context-averaged; internal alignment gaps can be treated as a fifth character, while terminal padding represents missing coverage. ([S13](references.md#s13))

“Fixed” means that this analysis treats the supplied topology and its existing branch lengths as given. The uncertainty explored within that analysis concerns the ancestral reconstruction and its placement, rather than every possible tree. A fifth-character gap model conditions on a chosen alignment; a full model of insertion and deletion histories would ask a broader question.

## Combining the two sources of evidence

In [Bayesian inference](https://en.wikipedia.org/wiki/Bayesian_inference), evidence updates prior probabilities to produce **posterior probabilities**. Schematically:

```text
posterior support for an ancestral explanation
    ∝ prior support from recombination and other assumptions
      × likelihood of the observed descendants
```

The proportionality sign means that the resulting values must be normalized so the probabilities of the considered alternatives sum to one. An ancestral explanation can include a sequence, a recombination path, and a placement on the descendant tree.

Strong sequence evidence can distinguish alternatives despite a moderate prior preference. When the observations contain little information, the prior contributes more strongly. A junction base missing from every observed read illustrates the latter situation: its reconstruction relies on the surrounding sequence, candidate segment paths, and model assumptions.

![Two sources of information, plausible recombination histories and related descendant sequences, combine in ancestor inference. A synthetic codon example shows that uncertainty about a whole sequence induces correlated uncertainty at its individual positions.](assets/04-uca-inference.svg)

*Figure 4. The information combined in UCA inference, followed by a deliberately simplified probability example. If the only ancestral codon possibilities are AAA with probability 0.6 and GGC with probability 0.4, the corresponding amino-acid probabilities are lysine 0.6 and glycine 0.4. Multiplying individual base probabilities would invent additional codons. All sequences and probabilities in this figure are illustrative.*

## Reading the inference controls

Swig offers conditional maximum-likelihood placement, grid averaging, and continuous Gibbs/Metropolis–Hastings sampling. It reports sequence and source uncertainty, coherent joint-path and marginal-consensus reconstructions, nucleotide/codon/amino-acid summaries, and sampling diagnostics where applicable. ([S13](references.md#s13))

The controls refer to familiar statistical operations that answer slightly different questions.

**Conditional optimization** selects a placement and ancestral-branch length with a high likelihood, then reports the sequence inference under that choice. Uncertainty in the chosen placement remains outside that conditional result.

**Grid averaging** evaluates a finite collection of placements and lengths and combines them with probability weights. This is an approximation to [marginalization](https://en.wikipedia.org/wiki/Marginal_distribution): summing or integrating over alternatives instead of choosing only one. The grid's range and resolution determine which alternatives receive consideration.

**Markov chain Monte Carlo**, or [MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo), visits possible explanations through a sequence of dependent random draws. After a suitable initial period, a well-behaved chain spends time in regions in proportion to their posterior probability. [Gibbs sampling](https://en.wikipedia.org/wiki/Gibbs_sampling) and [Metropolis–Hastings](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm) are ways to construct these updates.

A trace plot shows how sampled quantities change over iterations. A chain that stays near one value may reflect either concentrated evidence or poor exploration; multiple runs and movement among alternatives help distinguish these cases. **Effective sample size** estimates how much independent information remains after accounting for correlation between draws. More iterations can improve sampling precision, while model misspecification requires a different remedy.

## Reading the sequence outputs

A **marginal posterior** at one position gives the probabilities of its possible bases after combining the alternatives included in the inference. A position with `A: 0.98, G: 0.02` is more resolved under the model than one with `A: 0.55, G: 0.45`.

In a probability or frequency [sequence logo](https://en.wikipedia.org/wiki/Sequence_logo), each column has total height one, and each letter's height gives its probability. This convention differs from an information-content logo, where the total stack height also varies with conservation.

Sequence states can be correlated because neighbouring positions share a candidate allele or recombination path. A **joint** reconstruction chooses a coherent combination across positions. A **marginal consensus** chooses the most probable base separately at each position. The resulting consensus need not be the most probable whole sequence.

The same issue matters when translating uncertain sequence. In Figure 4, the only possible codons are `AAA` and `GGC`. Multiplying the three separate base distributions would incorrectly introduce possibilities such as `AGA`. Amino-acid uncertainty should instead be obtained from the joint codon distribution, summing the probabilities of codons that encode the same amino acid.

Source uncertainty and base uncertainty also differ. Several D placements may support the same nucleotide, producing uncertain origin but confident base identity. Inspecting both types of result is useful around junction boundaries.

## How much confidence to place in an ancestor

A reported posterior probability is conditional on the data and assumptions supplied. Important dependencies include lineage membership, the alignment, available germline alleles, recombination preferences, the substitution model, and the tree.

A practical sensitivity analysis asks which conclusions survive plausible alternatives. For an uncertain junction, vary defensible trimming or addition assumptions. For an ambiguous inherited allele, compare reconstructions under the supported references. For weak tree resolution, repeat the analysis with suitable alternative trees. For an indel-rich region, review the alignment and the consequences of its treatment.

Default parameter values also have a biological origin. Parameters derived from human heavy-chain data may be a poor description of another species, a light-chain locus, or TCR rearrangement. The relevant assumptions should be matched to that setting. The ordinary lack of TCR SHM further changes the biological question, often making the observed rearranged sequence itself the main object of interest. ([S13](references.md#s13), [R04](references.md#r04))

## From a sequence estimate to an experimental antibody

When an ancestor will be expressed, retain complete supported sequence alternatives at consequential uncertain sites, especially in CDR3. Reconstruct the paired chain with its own evidence and uncertainty. Expressing several plausible ancestors can reveal whether a biological conclusion depends on one uncertain reconstruction.

Binding, neutralization, expression, and stability measurements then test properties of those molecules. A reconstructed sequence provides a hypothesis about an ancestral receptor; experiments establish what that candidate molecule can do under the tested conditions. This connection is particularly useful in vaccine research.

---

[Previous: Mutations, alignments, and trees](10-shm-alignments-and-trees.md) · [Contents](index.md) · [Next: Vaccine design](12-vaccine-design.md)
