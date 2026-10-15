# 13. Finding the biology behind an analysis

The chapters follow receptor biology from development to immune responses and sequence-based reconstruction. This page provides a second route, starting from an analysis or result in Swig.

## Samples, assignment, and record preparation

| Analysis or decision | Biological question | Background |
|---|---|---|
| Dataset, sample, donor, time-point, and compartment labels | Which observations belong to the same person, specimen, or sampling occasion? | [Samples and counts](05-samples-sequences-and-counts.md) |
| FASTA, FASTQ, or AIRR input; quality filtering, trimming, subsampling | What was measured, how accurately, and which receptors could the assay detect? | [Samples and counts](05-samples-sequences-and-counts.md) |
| Species/locus reference selection; custom-reference metadata preparation | Which inherited segments and receptor landmarks are relevant? | [Recombination](02-vdj-recombination.md), [Assignment](06-vdj-assignment.md), [Inherited variation](08-germline-variation.md) |
| V(D)J assignment, calling profiles, productivity and constant-region calls | How was this receptor assembled, and what coding and isotype evidence is covered? | [Assignment](06-vdj-assignment.md) |
| Double-D screening | Does the junction support two ordered D contributions? | [Recombination](02-vdj-recombination.md), [Assignment](06-vdj-assignment.md#evidence-for-two-d-segments) |
| Exact collapse and denoising | Which differences should be retained as distinct observations, and which could arise from repeated measurement or error? | [Errors and chimeras](07-errors-and-chimeras.md) |
| CHMMAIRRa | Could amplification have combined parts of different templates? | [Errors and chimeras](07-errors-and-chimeras.md#reading-a-chimera-result) |
| Repertoire selection and cumulative working set | Which biological population and counting unit do the retained records represent? | [Errors and chimeras](07-errors-and-chimeras.md#filters-define-the-population-being-studied) |

The corresponding implementation details are in Swig's input, reference, assignment, Double-D, collapse, chimera, and storage specifications. ([S01](references.md#s01)–[S07](references.md#s07))

## Repertoire-level inference and diagnostics

| Analysis or result | Biological question | Background |
|---|---|---|
| Repertoire allele pooling: Dirichlet mixture or hurdle model | Can evidence across a donor's repertoire resolve ambiguous expressed-allele usage? | [Inherited variation](08-germline-variation.md#combining-ambiguous-allele-assignments-across-a-repertoire) |
| Gene-use, V–J, CDR3-length, identity, and isotype summaries | How is the sampled repertoire composed? | [Repertoires](09-repertoires-and-clones.md#reading-a-repertoire-summary) |
| Lineage assignment and sample membership | Which receptors plausibly descend from the same rearranged B cell, or meet the chosen TCR grouping definition? | [Clonal families](09-repertoires-and-clones.md) |
| V-SHM summaries, regional and codon breakdowns, weighted distributions | How much covered sequence differs from the assigned inherited baseline, and where? | [SHM summaries](10-shm-alignments-and-trees.md#measuring-sequence-differences-from-an-inherited-reference) |
| Missing-V warnings | Do independent families share a pattern suggesting an absent inherited reference? | [Missing V references](08-germline-variation.md#looking-for-a-missing-v-reference) |
| Personalized expressed V-set inference | Which supported inherited candidates best explain the donor's expressed V sequences? | [Personalized V sets](08-germline-variation.md#inferring-a-personalized-expressed-v-set) |
| Joint inherited/SHM model | How much inherited-sequence evidence remains under a competing recurrent-somatic explanation? | [Competing somatic explanations](08-germline-variation.md#a-competing-somatic-explanation) |
| Exact, substring, distance, and approximate sequence search | Where does a sequence or selected kind of similarity occur? | [Sequence retrieval](09-repertoires-and-clones.md#searching-for-a-sequence-or-its-neighbours) |
| CDR3 and inferred-germline neighbours, seed expansion, and manual merge | Was a family split, or does a similarity connection need a different interpretation? | [Neighbouring groups](09-repertoires-and-clones.md#reviewing-neighbouring-groups) |

The relevant method specifications cover lineage grouping, SHM and germline diagnostics, retrieval, repertoire allele refinement, and the separate joint inherited/SHM analysis. ([S08](references.md#s08)–[S10](references.md#s10), [S14](references.md#s14), [S15](references.md#s15))

## The selected-lineage workbench

| Analysis or result | Biological question | Background |
|---|---|---|
| Single-lineage input or opening selected groups | Is the supplied set a plausible family for ancestral interpretation? | [Clonal families](09-repertoires-and-clones.md#what-makes-a-b-cell-clone) |
| Reference-projected, de novo, and codon-aware alignments; manual editing | Which positions are homologous, and where is the evidence incomplete? | [Alignment](10-shm-alignments-and-trees.md#aligning-members-of-a-family) |
| Germline guide, consensus, and heuristic comparison sequence | Which parts of a proposed starting sequence come from references, observations, or unresolved assumptions? | [Rooting and guides](10-shm-alignments-and-trees.md#rooting-and-germline-guides), [UCA](11-unmutated-common-ancestors.md) |
| FastTree inference, uploaded Newick, rooting, branch display | Which branching histories are compatible with the aligned sequences? | [Phylogenies](10-shm-alignments-and-trees.md#what-a-phylogeny-represents) |
| Ancestral nucleotide and amino-acid mutation mapping | Where can changes be placed under the displayed reconstruction? | [Mutation mapping](10-shm-alignments-and-trees.md#mapping-changes-onto-a-tree) |
| Phylogenetic UCA: conditional optimization, grid, or MCMC | Which pre-SHM rearranged sequences and placements explain the sampled family? | [UCA inference](11-unmutated-common-ancestors.md) |
| Source tracks, posterior logos, codon/amino-acid probabilities, and alternative sequences | Which parts of the ancestor are resolved, and which uncertainties affect the molecule? | [UCA outputs](11-unmutated-common-ancestors.md#reading-the-sequence-outputs) |
| Alignment, tree, table, sequence, SVG, and result exports | Which observations and assumptions are needed to reproduce a conclusion or choose an experimental construct? | [UCA confidence](11-unmutated-common-ancestors.md#how-much-confidence-to-place-in-an-ancestor), [Vaccine design](12-vaccine-design.md) |

The single-lineage entry route assumes that its input belongs together. Assessing that assumption remains part of the biological analysis. The standard tree and mutation display and the phylogenetic UCA calculation are separate methods. ([S00](references.md#s00), [S11](references.md#s11), [S13](references.md#s13))

## A worked study: one donor sampled before and after immunization

Consider an invented study with blood samples taken before immunization and on two later occasions. The investigators have paired, UMI-supported B-cell sequences, cell identifiers, sample labels, and binding measurements for several antibodies. The questions are whether particular families expanded, how their variable regions changed, and which starting receptors could have given rise to them.

**Establish the observations.** Preserve the distinction between cells, molecules, and sequence records. Keep heavy and light chains linked through their measured cell identities. Record the receptor coverage and the procedure used to construct consensus sequences. The same donor label connects samples for appropriate donor-level analyses, while time-point and sample labels remain separate.

**Establish the reference comparisons.** Assign each chain using a suitable reference set. Review incomplete calls, productivity, and ambiguous alleles. Repertoire-level evidence can refine expressed-allele usage, and independent low-SHM families can provide evidence about a missing V reference. Reannotation after a reference change makes the downstream comparison internally consistent.

**Define the retained data.** Choose collapse keys that preserve the mutation variants and chain associations needed for the study. Review candidate errors and chimeras. Keep an account of how filters alter both record counts and measured cell counts.

**Propose families and compare samples.** Group heavy-chain sequences with a justified rule, then examine paired-light-chain consistency. Inspect neighbouring groups when shared mutations or plausible indels suggest a split family. Compare abundance on a cell basis where the assay supports it, using each sample's relevant cell denominator. A family absent from the baseline sample may have been present below its detection limit.

**Examine a selected history.** Review the family's nucleotide and codon alignment. Build or supply a tree and inspect mutation patterns alongside sampling times and binding measurements. A late, abundant branch with several replacements becomes a candidate for functional study; the combined observations still leave causality to be tested.

**Reconstruct and test earlier receptors.** Estimate heavy- and light-chain ancestors separately while preserving the biological pairing evidence. Review junction uncertainty, plausible references, tree sensitivity, and alternative complete sequences. Express the selected ancestral and intermediate antibodies and measure the activities relevant to the study. Test whether the conclusion changes across supported ancestral alternatives.

This example uses many analysis blocks because it asks a historical and experimental question. A study of TCR clonotype persistence might stop after careful assignment, counting, grouping, and longitudinal comparison.

## Saving the context of a result

A result is most useful when its inputs and choices remain recoverable. Preserve the exact reference sequences, study metadata, retained-set definition, grouping settings, edited alignment, tree, and any ancestor-model settings used for the reported result.

Swig's portable session records derived state and links it to an AIRR table that must also be retained. Its pipeline configuration export supports repeatable upstream processing; selected-lineage alignment and ancestor inference remain inspection-dependent workbench analyses. ([S12](references.md#s12), [S16](references.md#s16))

These records make it possible to distinguish a biological change between samples from a change in reference choice, filtering, or inference settings.

---

[Previous: Vaccine design](12-vaccine-design.md) · [Contents](index.md) · [Glossary](glossary.md)
