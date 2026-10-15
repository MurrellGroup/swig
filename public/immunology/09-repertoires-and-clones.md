# 9. Repertoires and clonal families

A repertoire can be examined at several levels: the distribution of gene use across a sample, groups of related receptors, or the individual variants within one group. These levels answer different questions about diversity, expansion, and history.

## Reading a repertoire summary

A V-gene usage plot counts the observations assigned to each V gene. A CDR3-length distribution counts receptors with junctions of different lengths. A V–J matrix counts combinations of segments used within a chain. An isotype plot summarizes supported constant-region calls.

For each display, the first question is the counting unit. A frequency of 20% could refer to reads, deduplicated sequence records, measured cells, or inferred clones. A single expanded family can dominate the first measures while contributing only one inferred clone. Missing or ambiguous calls also affect the denominator.

Swig's dashboard and selection tools expose gene-use, junction-length, identity, and related summaries. Ambiguous calls can be represented by a leading label or by splitting a record's contribution across labels. Equal fractional allocation is a counting convention; calibrated assignment uncertainty requires the corresponding probabilistic model. A V–J plot describes segment combinations within a chain. Heavy–light or alpha–beta pairing requires paired-chain observations. ([S05](references.md#s05))

These plots can reveal sample composition, technical differences, and candidates for closer inspection. A rise in one V gene might reflect expansion of one family, many independent families using that gene, a change in tissue composition, or an assay preference. Examining those alternatives connects a summary to a biological explanation.

## Richness, evenness, and sampling

Two repertoires can contain the same number of observed clones but distribute their cells very differently. One might have similar numbers in every clone; the other might be dominated by a few expanded clones. [Richness](https://en.wikipedia.org/wiki/Species_richness) describes how many types are present, and [evenness](https://en.wikipedia.org/wiki/Species_evenness) describes how evenly their abundances are distributed. These ecological concepts are also useful for receptor repertoires.

Observed richness grows with sampling depth because rare clones are easier to detect in larger samples. A comparison between 1,000 cells and 100,000 cells therefore needs to account for sampling. Read depth adds another layer because amplification and transcript abundance intervene between cells and reads. The calculation in [Chapter 5](05-samples-sequences-and-counts.md#subsampling-and-detection) gives the simplest detection model.

A descriptive repertoire display provides the observations for such comparisons. Formal diversity estimation or tests of enrichment require an explicitly chosen statistical analysis, including biological replicates and a defined counting unit.

## What makes a B-cell clone?

A biological [clone](https://en.wikipedia.org/wiki/Clone_(cell_biology)) comprises descendants of a common cell. In a B-cell receptor study, a clonal family usually means descendants of a cell with a particular successful receptor rearrangement. Those descendants inherit the rearrangement and may acquire different SHM patterns during the response.

The full receptor includes a heavy and a light chain. Heavy-chain-only studies approximate families using the available chain. A distinctive junction is valuable evidence because it contains information about segment choice, trimming, and added nucleotides. Related V and J assignments provide additional constraints. ([R12](references.md#r12), [R17](references.md#r17))

As mutations accumulate, members can differ in both the V region and CDR3. A grouping rule must tolerate some differences to recover the family. Tolerating too many can join independently rearranged receptors; tolerating too few can split a real family.

## Turning sequence similarity into proposed families

Swig's default lineage grouping uses productive records within a donor and locus, compatible V- and J-gene calls, equal CDR3 nucleotide lengths, and a CDR3 identity threshold. The documented default threshold is 85%. These choices define an operational grouping rule that can be changed for the study. ([S08](references.md#s08))

Identity is the fraction of compared positions that match. For two equal-length 20-base sequences, three differences give 17/20, or 85%, identity. Real junctions are longer or shorter, so the same percentage permits different numbers of differences.

Swig uses [single-linkage clustering](https://en.wikipedia.org/wiki/Single-linkage_clustering): sequences are joined when a path of qualifying pairwise matches connects them. Suppose A differs from B at three positions and B differs from C at three other positions. A and C can differ at six positions while still belonging to one connected group. This chaining behaviour is useful for recovering a diversifying family but can also connect separate families through intermediate sequences. ([S08](references.md#s08))

Confidence in a proposed family increases when several observations agree: distinctive related junctions, compatible inherited segments, shared derived mutations, supported chain pairing, and a plausible distribution across the donor's samples. A very large heterogeneous group deserves inspection even when every connecting pair passes the threshold.

Indels introduce another complication. An SHM-associated insertion or deletion can change CDR3 length and split relatives under an equal-length rule. Missing or incorrect V/J calls can do the same. Conversely, common short junctions and broadly compatible segment labels can make independently generated receptors hard to separate.

## T-cell groups require a different interpretation

Ordinary human and mouse T cells generally preserve their rearranged TCR sequence during an immune response. Expansion therefore produces many cells with the same receptor rather than a characteristic SHM-diversified family. Exact junction nucleotide sequence with compatible segment calls is a common starting point for TCR clonotype definitions; paired alpha–beta data improve specificity. ([R04](references.md#r04), [R22](references.md#r22))

Different T cells can recognize the same peptide–MHC complex using similar or quite different receptors. Some independently generated TCRs are shared across people and are called **public clonotypes**. Similarity can support a specificity hypothesis when combined with appropriate models and experiments. Cross-donor similarity reflects recurrent generation and selection, rather than descent from one lymphocyte in those donors. ([R22](references.md#r22), [R23](references.md#r23))

A permissive CDR3-similarity group is therefore best named according to its construction and intended interpretation. A TCR similarity cluster and a B-cell mutation lineage carry different biological claims.

## Following families across samples

A family observed at several time points provides evidence of persistence and can reveal changing abundance, isotype composition, or sequence diversity. A family found in two tissues can motivate questions about migration or shared responses.

The denominator still matters. A clone's fraction can rise because it expanded, because other cells contracted, or because the sample captured a different mixture. Absolute cell counts, consistent sampling, and biological replicates help distinguish these explanations.

Within-family mutation patterns provide another perspective. A later sample may contain variants nested within a branch represented earlier. That pattern is compatible with continued diversification, while unsampled cells and tissue reservoirs leave additional possible histories. Sampling date constrains when a sequence was observed; its origin can precede detection.

Swig supports dataset and sample labels alongside donor-level lineage groups. Summaries of group size and sample membership are useful starting points for choosing families to inspect in detail. ([S00](references.md#s00), [S08](references.md#s08))

## Searching for a sequence or its neighbours

Sequence retrieval begins with a chosen representation. A nucleotide query retains synonymous differences. An amino-acid query focuses on protein sequence and combines synonymous nucleotide variants because the [genetic code](https://en.wikipedia.org/wiki/Genetic_code) maps several codons to the same amino acid.

An exact search finds matching sequence. A substring search can find a motif within a longer sequence. Hamming distance measures substitutions between equal-length strings; edit distance also permits insertions and deletions. Each measure imposes a particular meaning of “nearby.”

Swig offers CDR3 or broader variable-region searches in nucleotide and amino-acid space, alongside approximate k-mer-based retrieval. A [k-mer](https://en.wikipedia.org/wiki/K-mer) is a short sequence word of length *k*. Shared words provide a fast similarity estimate; [MinHash](https://en.wikipedia.org/wiki/MinHash) compresses sets of words for comparison. Approximate retrieval is useful for finding candidates, followed by inspection with a more direct sequence comparison. ([S10](references.md#s10))

Similarity over a long V region can be dominated by the common inherited segment, whereas similarity within CDR3 gives more weight to junctional history. A short motif can recur in unrelated receptors. The most appropriate search region depends on whether the question concerns ancestry, a known receptor, or a possible binding-site motif.

## Reviewing neighbouring groups

Neighbour searches can reveal families that were split by an initial threshold or assignment difference. Swig can compare junction distances, vary V/J compatibility requirements, and expand from seed sequences through qualifying similarity links. Successive links can move far from the original seed, so the connecting path is informative. Donor and locus remain important boundaries. ([S10](references.md#s10))

Swig also searches for neighbours by comparing inferred starting sequences. A provisional sequence sketch screens candidate families, followed by reconstruction and sequence comparison for the retained candidates. This asks whether the proposed starting receptors resemble one another; shared inherited segments and reconstruction uncertainty affect that resemblance. ([S10](references.md#s10))

Opening several groups in one workbench allows direct comparison. A manual merge is a stronger decision: it changes the family interpretation used downstream. Swig preserves the original clone identifiers while recording a separate merged-lineage identifier. This makes the revised grouping traceable. ([S08](references.md#s08))

Before merging, inspect the junction alignment, shared mutation patterns, reference ambiguity, and paired-chain evidence where available. A tree constructed from the combined sequences will always impose a relationship on its input. Establishing a plausible family before tree inference gives that relationship its biological meaning.

---

[Previous: Inherited variation](08-germline-variation.md) · [Contents](index.md) · [Next: Mutations, alignments, and trees](10-shm-alignments-and-trees.md)
