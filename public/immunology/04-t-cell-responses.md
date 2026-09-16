# 4. T-cell responses

T and B cells share the use of rearranged antigen receptors, clonal expansion, and long-lived memory. Their receptors recognize targets in different ways, and their sequence histories usually reflect different processes.

## Presenting peptides to T cells

Most alpha–beta T cells recognize peptides bound to MHC molecules. Two presentation pathways are especially important.

[MHC class I](https://en.wikipedia.org/wiki/MHC_class_I) molecules are expressed on most nucleated cells. They commonly display peptides produced from proteins inside the cell. A virus-infected cell can therefore display peptides derived from viral proteins. Conventional [CD8 T cells](https://en.wikipedia.org/wiki/Cytotoxic_T_cell) inspect peptide–MHC class I complexes.

[MHC class II](https://en.wikipedia.org/wiki/MHC_class_II) molecules are especially associated with professional [antigen-presenting cells](https://en.wikipedia.org/wiki/Antigen-presenting_cell), including dendritic cells, B cells, and macrophages. They commonly display peptides from material processed in intracellular vesicles. Conventional [CD4 T cells](https://en.wikipedia.org/wiki/T_helper_cell) inspect peptide–MHC class II complexes. There are connections and exceptions to these routes, including [cross-presentation](https://en.wikipedia.org/wiki/Cross-presentation), in which externally acquired material is presented on MHC class I. ([R06](references.md#r06))

The CD4 and CD8 proteins act as co-receptors in these interactions. The TCR recognizes a particular peptide in the context of a particular MHC molecule. Because HLA genes vary extensively between people, individuals differ in which peptides they can present effectively. Receptor sequence, peptide, and HLA context all matter when studying a T-cell response.

## Selection in the thymus

T-cell precursors arise from bone-marrow-derived cells and develop in the [thymus](https://en.wikipedia.org/wiki/Thymus). They rearrange receptor genes and undergo [thymic selection](https://en.wikipedia.org/wiki/Thymocyte).

For conventional alpha–beta T cells, **positive selection** favours cells whose receptors can interact sufficiently with self peptide–MHC complexes. **Negative selection** removes many cells whose receptors respond too strongly to self material. Some self-reactive cells enter specialized regulatory pathways. Additional tolerance mechanisms act after cells leave the thymus.

These developmental steps reshape the repertoire produced by recombination. The receptors seen in mature T cells reflect both the generation process and selection for survival and function. Comparing their gene use or CDR3 composition requires both processes to be considered. ([R04](references.md#r04), [R07](references.md#r07), [R22](references.md#r22))

## Activation and cell function

Naive T-cell activation generally involves recognition of peptide–MHC together with [co-stimulation](https://en.wikipedia.org/wiki/Co-stimulation) and an appropriate cytokine environment. Dendritic cells are particularly important for initiating these responses. Signalling through the TCR-associated [CD3 complex](https://en.wikipedia.org/wiki/CD3_(immunology)) helps initiate changes in gene expression, cell division, and differentiation.

Activated CD8 T cells can become cytotoxic cells. They can kill target cells through mechanisms involving [perforin](https://en.wikipedia.org/wiki/Perforin) and [granzymes](https://en.wikipedia.org/wiki/Granzyme), which help trigger programmed cell death. They also produce cytokines.

CD4 T cells can develop several kinds of helper activity. Tfh cells support B-cell responses; other helper populations activate or recruit immune cells in different tissue and pathogen settings. [Regulatory T cells](https://en.wikipedia.org/wiki/Regulatory_T_cell) help restrain immune responses and maintain tolerance.

After expansion, many responding cells die during **contraction** of the response. Surviving [memory T cells](https://en.wikipedia.org/wiki/Memory_T_cell) include populations that circulate and populations retained in tissues. Their location affects which cells appear in a blood sample. ([R05](references.md#r05), [R07](references.md#r07))

## What remains stable in a T-cell clone

During ordinary human and mouse T-cell responses, the rearranged TCR sequence is generally maintained as the cell divides. T cells do not normally diversify their receptors through the AID-dependent affinity-maturation programme described for B cells. The response can still change substantially as different clones expand, contract, migrate, or acquire different functional states. ([R04](references.md#r04), [R06](references.md#r06))

This distinction changes the meaning of sequence comparisons. A set of similar, mutated BCRs can represent descendants of one rearrangement. A set of similar TCRs can instead contain receptors assembled independently in different cells. Sequence similarity may motivate a shared-specificity hypothesis, while evidence for shared cellular ancestry requires a stricter interpretation.

For TCR clone tracking, exact productive junction sequence together with V/J assignments is a useful operational starting point. Paired alpha–beta sequences provide stronger resolution. A matching beta chain alone leaves uncertainty about the alpha chain and can occasionally group independently generated receptors. Different nucleotide sequences can also encode the same amino-acid sequence.

Specialized methods use recurring TCR motifs or structural sequence similarities to identify possible antigen-specificity groups. Those groups answer a different biological question from clonal ancestry. Functional testing, peptide–MHC binding data, and HLA information strengthen their interpretation. ([R23](references.md#r23))

In Swig, a permissive CDR3-similarity group containing TCRs should therefore be treated as a sequence group until the intended clone definition has been checked. A tree of such a group mainly summarizes sequence differences; a B-cell-style maturation history would require additional biological evidence.

## Gamma–delta and other unconventional T cells

[Gamma–delta T cells](https://en.wikipedia.org/wiki/Gamma_delta_T_cell) use rearranged gamma and delta chains. Their recognition systems are diverse and often involve ligands outside conventional peptide–MHC presentation. Their distribution and receptor usage also vary across tissues and species.

Some alpha–beta populations have unconventional restrictions as well. [Natural killer T cells](https://en.wikipedia.org/wiki/Natural_killer_T_cell) can recognize lipid antigens presented by CD1d, while [mucosal-associated invariant T cells](https://en.wikipedia.org/wiki/Mucosal-associated_invariant_T_cell) recognize particular metabolites presented by MR1. These examples explain why a receptor's locus and biological source should remain attached to its sequence record.

---

[Previous: B-cell responses](03-b-cell-responses.md) · [Contents](index.md) · [Next: Samples, sequences, and counts](05-samples-sequences-and-counts.md)
