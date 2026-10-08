# References

The chapter links use `R` identifiers for biological background, research, and standards, and `S` identifiers for Swig specifications. Wikipedia links within the chapters provide additional routes into individual concepts.

## Biological background, research, and standards

### R01

Janeway CA Jr, Travers P, Walport M, Shlomchik MJ. **Antigen Recognition by B-cell and T-cell Receptors.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK10770/).

### R02

Janeway CA Jr and colleagues. **B-cell activation by armed helper T cells.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK27142/). Background on linked recognition and T-cell-dependent B-cell responses.

### R03

Janeway CA Jr and colleagues. **The generation of diversity in immunoglobulins.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK27140/).

### R04

Janeway CA Jr and colleagues. **The rearrangement of antigen-receptor gene segments controls lymphocyte development.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK27113/).

### R05

Janeway CA Jr and colleagues. **Immunological memory.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK27158/).

### R06

Janeway CA Jr and colleagues. **Antigen recognition by T cells.** In *Immunobiology*, 5th edition (2001). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK27098/).

### R07

Alberts B, Johnson A, Lewis J, Raff M, Roberts K, Walter P. **Lymphocytes and the Cellular Basis of Adaptive Immunity.** In *Molecular Biology of the Cell*, 4th edition (2002). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK26921/).

### R08

Alberts B and colleagues. **The Generation of Antibody Diversity.** In *Molecular Biology of the Cell*, 4th edition (2002). [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK26860/). Foundational receptor assembly and expression; the molecular account of AID is supplemented by R09 and R30.

### R09

Muramatsu M et al. (2000). **Class switch recombination and hypermutation require activation-induced cytidine deaminase (AID), a potential RNA editing enzyme.** *Cell* 102:553–563. [DOI: 10.1016/S0092-8674(00)00078-7](https://doi.org/10.1016/S0092-8674(00)00078-7). Experimental evidence for AID's requirement in both processes; the historical title predates the evidence for DNA deamination.

### R10

Gitlin AD, Shulman Z, Nussenzweig MC (2014). **Clonal selection in the germinal centre by regulated proliferation and hypermutation.** *Nature* 509:637–640. [DOI: 10.1038/nature13300](https://doi.org/10.1038/nature13300).

### R11

Yaari G et al. (2013). **Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data.** *Frontiers in Immunology* 4:358. [Article](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2013.00358/full). Source of the human five-nucleotide-context mutation model discussed in the guide.

### R12

Gupta NT et al. (2017). **Hierarchical clustering can identify B cell clones with high confidence in Ig repertoire sequencing data.** *Journal of Immunology* 198:2489–2499. [DOI: 10.4049/jimmunol.1601850](https://doi.org/10.4049/jimmunol.1601850). Research on sequence-based B-cell clone identification.

### R13

AIRR Community. **Rearrangement Schema.** AIRR Standards documentation. [Specification](https://docs.airr-community.org/en/stable/datarep/rearrangements.html). Field meanings, junction conventions, coordinates, and annotated-repertoire representation.

### R14

Vander Heiden JA et al. (2018). **AIRR Community standardized representations for annotated immune repertoires.** *Frontiers in Immunology* 9:2206. [Article](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2018.02206/full).

### R15

Safonova Y, Pevzner PA (2020). **V(DD)J recombination is an important and evolutionarily conserved mechanism for generating antibodies with unusually long CDR3s.** *Genome Research* 30:1547–1558. [DOI: 10.1101/gr.259598.119](https://doi.org/10.1101/gr.259598.119).

### R16

Corcoran MM et al. (2016). **Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.** *Nature Communications* 7:13642. [DOI: 10.1038/ncomms13642](https://doi.org/10.1038/ncomms13642). Primary research introducing IgDiscover and demonstrating the value of individual V-reference inference.

### R17

Ralph DK, Matsen FA IV (2016). **Likelihood-based inference of B cell clonal families.** *PLOS Computational Biology* 12:e1005086. [Article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005086). A probabilistic approach to clonal-family inference, distinct from Swig's threshold-based grouping.

### R18

Dhar A, Ralph DK, Minin VN, Matsen FA IV (2020). **A Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis.** *PLOS Computational Biology* 16:e1008030. [Article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008030). The linearham formulation jointly addresses recombination and phylogeny; Swig's fixed-observed-tree UCA model is a separate implementation with a narrower uncertainty calculation.

### R19

Leggat DJ et al. (2022). **Vaccination induces HIV broadly neutralizing antibody precursors in humans.** *Science* 378:eadd6502. [DOI: 10.1126/science.add6502](https://doi.org/10.1126/science.add6502), [PubMed](https://pubmed.ncbi.nlm.nih.gov/36454825/). Human germline-targeting trial used as a dated example in Chapter 12.

### R20

Liao HX et al. (2013). **Co-evolution of a broadly neutralizing HIV-1 antibody and founder virus.** *Nature* 496:469–476. [DOI: 10.1038/nature12053](https://doi.org/10.1038/nature12053). Longitudinal antibody and viral analysis of the CH103 lineage.

### R21

Shugay M et al. (2014). **Towards error-free profiling of immune repertoires.** *Nature Methods* 11:653–655. [DOI: 10.1038/nmeth.2960](https://doi.org/10.1038/nmeth.2960). Molecular barcoding and error correction in repertoire sequencing.

### R22

Elhanati Y, Murugan A, Callan CG Jr, Mora T, Walczak AM (2014). **Quantifying selection in immune receptor repertoires.** *PNAS* 111:9875–9880. [DOI: 10.1073/pnas.1409572111](https://doi.org/10.1073/pnas.1409572111), [author manuscript](https://arxiv.org/abs/1404.4956). Recombination, selection, and shared TCR sequences.

### R23

Glanville J et al. (2017). **Identifying specificity groups in the T cell receptor repertoire.** *Nature* 547:94–98. [DOI: 10.1038/nature22976](https://doi.org/10.1038/nature22976), [PubMed](https://pubmed.ncbi.nlm.nih.gov/28636589/). TCR similarity and experimentally supported specificity grouping.

### R24

Edgar RC, Flyvbjerg H (2015). **Error filtering, pair assembly and error correction for next-generation sequencing reads.** *Bioinformatics* 31:3476–3482. [DOI: 10.1093/bioinformatics/btv401](https://doi.org/10.1093/bioinformatics/btv401).

### R25

Kumar V et al. (2019). **Long-read amplicon denoising.** *Nucleic Acids Research* 47:e104. [DOI: 10.1093/nar/gkz657](https://doi.org/10.1093/nar/gkz657). Methodological background for the FAD-compatible denoising route; Swig's other denoising variants are documented separately.

### R26

Chernyshev M et al. (2025). **Detection of PCR chimeras in adaptive immune receptor repertoire sequences.** *Bioinformatics* 41:btaf576. [DOI: 10.1093/bioinformatics/btaf576](https://doi.org/10.1093/bioinformatics/btaf576). The CHMMAIRRa method.

### R27

Felsenstein J (1981). **Evolutionary trees from DNA sequences: a maximum likelihood approach.** *Journal of Molecular Evolution* 17:368–376. [DOI: 10.1007/BF01734359](https://doi.org/10.1007/BF01734359). Foundational phylogenetic likelihood calculation.

### R28

Price MN, Dehal PS, Arkin AP (2010). **FastTree 2—approximately maximum-likelihood trees for large alignments.** *PLOS ONE* 5:e9490. [DOI: 10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490).

### R29

Kepler TB (2013). **Reconstructing a B-cell clonal lineage. I. Statistical inference of unobserved ancestors.** *F1000Research* 2:103. [DOI: 10.12688/f1000research.2-103.v1](https://doi.org/10.12688/f1000research.2-103.v1), [author manuscript](https://arxiv.org/abs/1303.0424). Statistical ancestor inference and its uncertainty.

### R30

Petersen-Mahrt SK, Harris RS, Neuberger MS (2002). **AID mutates E. coli suggesting a DNA deamination mechanism for antibody diversification.** *Nature* 418:99–103. [DOI: 10.1038/nature00862](https://doi.org/10.1038/nature00862).

### R31

Marcou Q, Mora T, Walczak AM (2018). **High-throughput immune repertoire analysis with IGoR.** *Nature Communications* 9:561. [DOI: 10.1038/s41467-018-02832-w](https://doi.org/10.1038/s41467-018-02832-w). Probabilistic modelling of rearrangement. Swig uses compact rearrangement defaults informed by public models, rather than running IGoR inference.

## Swig implementation specifications

These references describe the public repository reviewed on **12 September 2026**, whose package and methods index identified **version 0.38.10**. The links follow the repository's `main` branch. Interface details and method defaults can change with subsequent releases.

### S00

MurrellGroup. **Swig repository and methods index.** [Repository](https://github.com/MurrellGroup/swig), [methods index](https://github.com/MurrellGroup/swig/blob/main/public/METHODS_INDEX.md), [package version](https://github.com/MurrellGroup/swig/blob/main/package.json).

### S01

**Input, quality control, and study metadata.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/01_INPUT_QC_AND_STUDY.md).

### S02

**Reference preparation.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/02_REFERENCE_PREPARATION.md).

### S03

**V(D)J assignment.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/03_VDJ_ASSIGNMENT.md).

### S04

**Double-D screening.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/04_DOUBLE_D_SCREEN.md).

### S05

**Storage, dashboard, and selection.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/05_STORAGE_DASHBOARD_SELECTION.md).

### S06

**Collapse and denoising.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/06_COLLAPSE_AND_DENOISING.md).

### S07

**Chimera inference.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/07_CHIMERA_INFERENCE.md).

### S08

**Lineage assignment.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/08_LINEAGE_ASSIGNMENT.md).

### S09

**SHM and germline diagnostics.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/09_SHM_AND_GERMLINE_DIAGNOSTICS.md). Includes missing-V warnings and personalized expressed V-set inference.

### S10

**Retrieval and neighbours.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/10_RETRIEVAL_AND_NEIGHBOURS.md).

### S11

**Alignment and phylogeny.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/11_ALIGNMENT_AND_PHYLOGENY.md).

### S12

**Execution, persistence, exports, and cancellation.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/12_EXECUTION_PERSISTENCE_CANCELLATION.md).

### S13

**Phylogenetic UCA inference.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/PHYLO_UCA_INFERENCE.md). Fixed-tree placement, recombination HMM, character models, inference routes, posterior summaries, and exports.

### S14

**Repertoire allele refinement.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/REPERTOIRE_ALLELE_REFINEMENT.md). Expressed-allele pooling, Dirichlet mixture, hurdle model, and reassignment policies.

### S15

**Joint inherited/SHM model.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/10_JOINT_INHERITED_SHM_MODEL.md). Separate experimental inherited-versus-somatic analysis.

### S16

**Command-line pipeline and lineage-study interchange.** [Specification](https://github.com/MurrellGroup/swig/blob/main/public/methods/13_CLI_AND_LINEAGE_STUDY.md).

---

[Contents](index.md) · [Analysis guide](13-analysis-guide.md) · [Glossary](glossary.md)
