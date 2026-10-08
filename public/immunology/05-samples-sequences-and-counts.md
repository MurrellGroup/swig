# 5. Samples, sequences, and counts

A receptor sequence reaches an analysis through a chain of experimental steps. Each step determines which molecules can be observed and what a count means. These details affect even simple questions such as whether one clone is more abundant than another.

## From cells to reads

An experiment begins with a biological sample: perhaps blood, a lymph-node aspirate, or a tissue biopsy. Researchers may analyse all lymphocytes or enrich a subset, such as memory B cells or CD8 T cells. [Flow cytometry](https://en.wikipedia.org/wiki/Flow_cytometry) can identify cells by surface markers, and fluorescent antigen probes can help collect antigen-binding B cells.

Receptor sequences can be obtained from genomic DNA or from RNA. With RNA, [reverse transcription](https://en.wikipedia.org/wiki/Reverse_transcriptase) first makes a complementary DNA copy. [Polymerase chain reaction](https://en.wikipedia.org/wiki/Polymerase_chain_reaction), or PCR, then amplifies selected regions before sequencing. The instrument produces **reads**: observations of nucleotide sequences from the library.

Targeted assays differ in which chains, isotypes, or sequence spans they recover. Primers can favour some genes over others or contribute sequence that needs trimming. Single-cell assays can retain a cell barcode, allowing receptor records to be associated with the same cell. These assay details belong with the results. ([R14](references.md#r14), [R21](references.md#r21))

## Reads, molecules, cells, and clones

Four levels recur throughout repertoire analysis.

A **read** is an instrument observation. A **molecule** is an original DNA or RNA template. A **cell** can contribute one or more receptor molecules. A **clone** contains descendants of an ancestral lymphocyte and can include many cells and, for B cells, many related receptor variants.

Consider an illustrative experiment containing 20 cells from one clone. Suppose each cell contributes 50 receptor RNA molecules and amplification produces 10 reads per molecule. That clone could yield 10,000 reads. A different clone with the same cell count could yield a different read count because its cells express less receptor RNA or its templates amplify less efficiently.

Genomic-DNA sequencing avoids variation in receptor-transcript abundance, while still depending on template recovery, amplification, rearrangement state, and counting conventions. RNA-based assays can be particularly influenced by antibody-secreting cells, which produce large quantities of immunoglobulin transcript.

A [unique molecular identifier](https://en.wikipedia.org/wiki/Unique_molecular_identifier), or **UMI**, is a short tag attached early enough in library preparation to mark an original molecule. Reads sharing a UMI and appropriate sequence context can be combined into a consensus. Cell barcodes label cells; UMIs label molecules. Their interpretation depends on how the library was constructed and on error correction. ([R14](references.md#r14), [R21](references.md#r21))

Swig propagates multiplicities through `duplicate_count`. The biological meaning comes from the input: raw-read multiplicity, a count supplied by upstream processing, or another explicitly defined unit. Any abundance plot should name that unit. ([S05](references.md#s05), [S06](references.md#s06))

## Paired chains

An antibody binding site combines heavy and light chains. An alpha–beta TCR combines alpha and beta chains. Bulk sequencing of separate chains usually loses their cellular pairing. Single-cell methods can preserve it through shared cell identifiers.

Paired-end sequencing has a different meaning: two reads are obtained from opposite ends of one library fragment. It can improve sequence reconstruction without identifying the naturally paired receptor chain. For example, two reads covering one heavy-chain fragment still provide information about a heavy chain.

When sequences are exported for further analysis, retain available cell and pairing metadata. Swig's chain-level annotation and CDR3 grouping can be used alongside paired data, but pairing evidence has to enter the analysis explicitly. The lineage grouping described here uses individual-chain sequence features. ([R14](references.md#r14), [S08](references.md#s08))

## FASTA, FASTQ, and AIRR tables

[FASTA](https://en.wikipedia.org/wiki/FASTA_format) stores sequence identifiers and nucleotide strings. [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) adds a quality score for each base. An **AIRR rearrangement table** stores receptor sequences with structured fields such as locus, V/D/J calls, junction sequence, productivity, and sample identifiers. The [AIRR Community](https://docs.airr-community.org/) maintains the standard. ([R13](references.md#r13), [R14](references.md#r14))

Swig accepts these input types, including compressed files. In the web workflow, imported AIRR sequences undergo a new annotation using the selected references and settings. Preserve an original table when comparing annotation methods. Reads requiring adapter removal, paired-read assembly, or UMI consensus should undergo those operations in appropriate upstream processing. ([S01](references.md#s01), [S00](references.md#s00))

## Base quality and trimming

A [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score) expresses the estimated probability of a base-calling error. Q20 corresponds to an error probability of 0.01; Q30 corresponds to 0.001. The relationship is `error probability = 10^(-Q/10)`.

Adding these probabilities across a read gives its **expected number of base-call errors**. For example, 300 bases each at Q30 have 0.3 expected errors in total. This is an expectation across comparable reads; any particular read can contain zero, one, or several errors.

Removing a poor-quality terminal region can improve the retained sequence. Trimming also removes biological information. Losing the end of J or the start of C can change whether a complete junction or an isotype can be assigned. Swig offers optional 3′ trimming and expected-error filtering for FASTQ input before assignment. Their thresholds should reflect the assay and the sequence coverage required. ([R24](references.md#r24), [S01](references.md#s01))

## Sample structure is part of the analysis

A **dataset** is an input library or file. A **sample** represents a biological collection. A **donor** or **subject** identifies the individual. A **timepoint** and **compartment** describe when and where it was collected. A **cohort** groups individuals for a study comparison.

Two sequencing files can be technical replicates of one sample. Blood samples taken before and after vaccination can share a donor while remaining distinct samples. Different donors can generate similar receptors independently.

These relationships motivate different processing boundaries. Combining technical replicates may be appropriate for molecule or sequence counts. Following a clone across time requires linking samples from the same donor. Pooling different donors under one ancestry label can create biologically misleading families. Swig represents these boundaries explicitly. ([S01](references.md#s01))

## Subsampling and detection

Randomly retaining a subset of eligible reads can make an exploratory analysis smaller or give datasets comparable sequencing depth. Swig's seeded subsampling provides a repeatable random selection within each dataset. ([S01](references.md#s01))

Subsampling reduces the chance of observing rare variants. In a simple model with independent observations, a sequence present at frequency `p` has probability `(1 - p)^n` of being missed in `n` observations. Real libraries depart from this model through amplification and unequal sampling, but the calculation illustrates why absence in one sample is weak evidence for biological absence.

---

[Previous: T-cell responses](04-t-cell-responses.md) · [Contents](index.md) · [Next: V(D)J assignment](06-vdj-assignment.md)
