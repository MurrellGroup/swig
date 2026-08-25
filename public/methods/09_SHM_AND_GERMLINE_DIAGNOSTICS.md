# Post-lineage analysis: SHM, missing-V evidence, and personalized V inference

## SHM metrics

Swig compares `v_sequence_alignment` with `v_germline_alignment` (falling back to complete aligned fields). Only A/C/G/T aligned pairs enter nucleotide denominators. Query coordinates are advanced through query bases and AIRR FWR1/2/3 and CDR1/2 intervals classify positions.

Reported choices are:

- V nucleotide mismatch count and rate;
- translated codon replacement count and rate;
- synonymous changed-codon count;
- CDR1/2 nucleotide mismatch rate;
- framework 1/2/3 nucleotide mismatch rate.

Codons use the AIRR `sequence_frame`/`v_frame`; an indel clears the partial aligned codon and resumes at the next frame boundary. A changed codon is counted once as synonymous or replacement, not once per changed nucleotide. Missing/ambiguous codons are excluded.

`duplicate_count` weights abundance summaries, means, medians, and 95th percentiles. Scalar analyzed/skipped counts cover every active row. Plot/session records are memory bounded: at most the configured number per lineage × sample and a deterministic global reservoir of 100,000. Thus a plotted point sample can be bounded even while aggregate record/abundance totals remain complete.

These are **direct descriptive mismatch summaries**, not a context-dependent SHM targeting model, selection test, or phylogenetic mutation reconstruction. Yaari et al.'s S5F work is relevant biological context but is not implemented here: [Models of somatic hypermutation targeting and substitution](https://pmc.ncbi.nlm.nih.gov/articles/PMC3828525/).

## Possible missing V allele screen

This diagnostic requires lineage assignments and treats each donor × lineage as one independent unit. It is a two-pass warning system, not genotype inference.

### Pass 1: discovery representatives

For each eligible low-SHM V alignment, Swig maps observed mismatches to germline coordinates. Per reported V allele it identifies recurrent linked substitution patterns subject to minimum aligned bases, unit/coverage support, allele fraction, binomial-tail probability, maximum SNP count, and diversity in J calls, junction lengths, and CDR3 fingerprints. AID WRCY/RGYW context is annotated but does not waive statistical thresholds.

Default discovery requirements include at least 6 independent supporting lineages, 20 covered units, allele fraction 0.2, SHM rate at most 0.08, at least 180 aligned bases, at most 6 candidate SNPs, and binomial survival probability at most \(10^{-6}\).

### Pass 2: all-member validation

For each proposed linked pattern, Swig scans every retained member of the relevant lineages. It requires joint coverage and linked alternate support, counts units retaining the reference or a conflicting base, requires at least 3 near-germline units, minimum linked fraction 0.9, and maximum other-alternate fraction 0.02. It also requires default diversity of 3 J calls, 3 junction lengths, and 6 distinct CDR3 fingerprints. This second pass prevents one convenient lineage representative from hiding contradictory descendants.

Candidates can be exported as proposed V FASTA after explicit selection. No reference is added to the running assignment database automatically.

## Personalized expressed V-set inference

This is a separate model from the high-specificity missing-V warning. It uses every assigned lineage with a usable V alignment, with exactly one vote per lineage.

### Representative rule

Swig selects the member with the lowest nucleotide mismatch rate in its **current** `v_sequence_alignment` versus `v_germline_alignment`. Ties prefer more aligned V bases and then the earlier AIRR ordinal. There is no random draw, temporal “early” requirement, allele-invariant mask, or circularity correction in this selection step. Lower-SHM observations are used simply because they carry more information about the rearranged germline allele.

### Candidate and likelihood model

Each subject, locus, and called V gene is fitted independently. Candidate database alleles must have the same ungapped length as the current allele and fall within the configured Hamming radius. Recurrent, linked substitution patterns can additionally create same-length novel candidates; insertions, deletions, boundary changes, cross-gene swaps, and D/J hypotheses are not admitted.

For each lineage observation, mutation exposure is estimated from covered positions that are invariant across the local candidate set. Candidate likelihoods use a fixed five-nucleotide-window motif prior:

- WRCY/RGYW AID hot spots have relative mutability 5;
- SYC/GRS cold spots have relative mutability 0.25;
- WA/TW polymerase-eta hot spots have relative mutability 2;
- other contexts have relative mutability 1.

This is a cheap fixed context approximation, **not** the empirical 1,024-entry human-heavy S5F table. A small configurable sequencing-error floor is mixed into every site probability. Allele-discriminating positions and their ±2-nucleotide context halo do not contribute to the observation's mutation-exposure estimate.

### Stepwise allele-set search

Within each expressed gene, Swig begins with the best penalized single allele. It tests an inactive candidate by optimizing the frequency transferred from the current mixture, using only the candidate's per-lineage emission and the cached current mixture. An accepted addition is followed by an exact EM frequency refit. Backward deletion removes components that become redundant.

The objective is log likelihood minus a BIC penalty for mixture-frequency parameters and every learned nucleotide in a novel candidate. The optional extra log-evidence threshold is applied after that penalty. Known and novel candidates therefore use the same likelihood, while data-derived sequences pay for the bases learned from the same repertoire.

The downloadable V FASTA is conservative: inferred active alleles replace only genes actually tested for the selected subject/locus. Untested genes and other loci remain unchanged. Novel records inherit the exact parent's `SWIGMETA` annotation. The download is a candidate reference and must be used in a complete new assignment run; the current analysis is never silently changed.

SHM, possible-missing-V evidence, and personalized inference are separate full-width result sections. They share a post-lineage location because all consume lineage assignments; none modifies another result.

## Prominent incomplete-reference escalation

Swig displays a blocking-style visual warning when either condition holds:

- at least one retained candidate is supported by more than 50 independent lineages; or
- more than five candidates with support at or below 50 independent lineages pass the two evidence stages.

The warning tells the user to rerun with a more complete germline database or use personalized germline discovery with **IgDiscover**. It does not present the candidate FASTA as a validated genotype and does not automatically rerun assignment.

## Literature relationship

The screen is **custom and deliberately conservative**. It shares the problem of subject-specific allele discovery with TIgGER, IgDiscover, and partis germline inference but does not implement their genotype/novel-allele models:

- [TIgGER subject-specific allele identification](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.00129/full)
- [IgDiscover individualized V databases](https://doi.org/10.1038/ncomms13642)
- [partis per-sample germline inference](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007133)

## Limitations

- AIRR alignment or lineage errors can create recurrent artifacts.
- Clonal dependence is reduced, not eliminated, by one donor-lineage unit.
- Expressed repertoires cannot prove genomic presence/absence or copy number.
- The binomial background is simplified and substitutions are not modeled with S5F context probabilities.
- A warning must be validated by dedicated germline inference and, where important, genomic evidence.
- Personalized inference covers expressed same-length V alleles only. A silent allele, unobserved gene, length-changing allele, or cross-gene reassignment cannot be resolved by this implementation.
- The fixed motif likelihood is a computational approximation and may be misspecified for non-human species or non-heavy-chain loci.
