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

### Candidate graph and homologous coordinates

Each subject/locus is processed separately. Same-length V references may compete across gene labels only along explicitly substitution-compatible paths. The affine global alignment scores are match +2, mismatch -3, gap open -5, extension -1. An ungapped path may tie the optimum; a strictly better gapped path excludes that pair. This prevents compensating indels from authorizing raw-coordinate swaps. Every observation scores only directly compatible candidates, even inside a connected component.

SHM calibration uses the published 1,024-context human silent five-mer targeting and substitution tables (HS5F; provenance in `src/shm-model/README.md`). Observed reference sequences are globally aligned to a common longest reference within each subject/locus/V family. Only mapped homologous positions and the same reference-to-alternate event share calibration. The tested physical gene label is excluded. Gapped columns without an anchor mapping do not borrow evidence.

### Novel-haplotype validity test

A proposal is a linked substitution haplotype. Mutation burden for that test excludes precisely its complete proposed site set, plus the uncertain last 12 V positions; it never masks the union of other proposals. Each proposed alternate has an effective SHM hazard and a Gamma random-effect prior learned across other aligned genes. Genes contribute equally to between-gene dispersion; large base counts do not collapse that dispersion to zero. Sparse calibration falls back to a unit-mean exponential rate prior. The prior is evaluated on log rate with its Jacobian included.

The null explains each alternate through its SHM hazard. The alternative adds one lineage-level mixture component for the complete linked allele. Covered site states and mutation burden (rounded to 0.0005) are grouped as sufficient statistics. Candidate-specific other-mutation exposure is held identical under null and alternative. Alternate-specific hazards can approach probability one rather than being capped at a fixed human substitution fraction.

The test charges half log N for the added mixture frequency plus `log C(L,s) + s log 3 + log P` for discovering s nucleotide changes among L positions and P expressed parent sequences. Nested proposals are tested as additions to the most specific accepted subset, with support breaking ties. A second sequence-level conditional test compares weak proposals with better-supported accepted alleles on different reference backbones; inherited backbone bases must be supported too. This is a proposal-validity test, not a claim of a fully integrated posterior over all allele sets.

### Common final emission model and set selection

All admitted known and novel sequences compete using the same categorical mutation model: HS5F alternate-specific rates multiplied by the mean of the leave-gene-out Gamma calibration. A lineage's exposure is anchored to its actual current parent alignment, excluding the uncertain terminal window, and cannot change with the size of the proposed candidate catalogue. The three alternate hazards define a total mutation hazard and a normalized substitution distribution; the sequencing-error floor completes the nucleotide probabilities.

The final mixture does **not** multiply a ratio from the haplotype null into an uncorrected parent emission. That would make scores from different parent alleles incomparable. The joint validity test and the final common emission model are distinct documented stages.

Greedy additions optimize a new mixture frequency, followed by an EM refit. Backward removal refits the remaining components. The set objective charges mixture-frequency BIC and the full sequence-search code length for each novel sequence. The active-allele setting is a generous computational guard per gene label, scaled to the number of observed genes in a component; it is not a diploid genotype or copy-number prior.

### Terminal V uncertainty

AIRR `sequence`, `rev_comp`, and `v_sequence_start` reconstruct raw downstream bases from an aligned anchor 12 nt before the V reference end. Query insertions advance raw query coordinates but never become a second observation at a germline position. Both matching and nonmatching downstream bases enter symmetrically. Without raw sequence/anchor coverage, terminal evidence is missing.

The likelihood sums over V trimming states 0–12 (the last state includes deeper trimming). Trimmed bases follow a learned junction composition, not necessarily uniform A/C/G/T. Junction composition and trimming sufficient statistics come from other genes in the same subject/locus. Trimming uses a fixed broad geometric prior for one E step followed by leave-gene-out aggregation. `d_sequence_start` or `j_sequence_start`, when available, bounds the downstream composition sample. The terminal model is approximate; an unrecovered terminal base must not be counted as an exact allele recovery.

### CLI

`swig-cli personalized-germline --airr sample.processed.airr.tsv.gz --v-reference V.fasta --out personalized`

Use the complete processed AIRR table with numeric `clone_id` or `lineage_id`, not a lineage-viewer export limited to the top displayed lineages. The command writes the dashboard JSON (including proposal evidence), evidence TSV and per-pool V FASTA. Web and CLI call the same TypeScript engine. New inference remains in a lazy Web Worker in the browser; V(D)J assignment remains WASM.

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
- Personalized inference covers expressed substitution-compatible V alleles. Silent genes, unsupported terminal bases and indel-containing alternatives remain unidentifiable.
- The HS5F calibration, log-rate MAP test, fixed burden estimates, approximate boundary background and two-stage search are computational approximations. Paired dropout controls do not establish the truth of stable background candidates.
