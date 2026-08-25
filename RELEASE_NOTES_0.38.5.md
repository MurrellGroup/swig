# Swig 0.38.5

## Personalized expressed V-set inference

- Adds a separate **Infer expressed V allele set** analysis beside the existing high-specificity missing-V screen.
- Uses all assigned lineages without the 10,000-lineage display-summary cap. Exactly one usable observation contributes per subject-by-lineage group.
- Selects the observation with the lowest mismatch rate under its current V assignment; ties prefer more aligned V bases and then the earlier AIRR ordinal. It does not choose a random observation or impose an early-timepoint rule.
- Fits each subject, locus, and called V gene independently. Candidate known alleles must be from the same gene, have the same ungapped length, and fall within the configured SNP radius.
- Can propose recurrent linked same-length substitution haplotypes. Length-changing and cross-gene hypotheses are deliberately excluded.
- Uses a cheap fixed context-dependent SHM approximation with WRCY/RGYW and WA/TW hot spots, SYC/GRS cold spots, and a configurable sequencing-error floor. It is explicitly not an S5F implementation.
- Uses cached per-lineage emissions, greedy candidate addition, exact EM mixture refits, backward removal, and BIC penalties for mixture parameters and bases learned in novel candidates.
- Exports evidence as CSV, TSV, or JSON Lines and a conservative personalized V FASTA. Only genes tested for the selected subject/locus are pruned; all untested genes and other loci remain. The current assignment is never silently changed.
- Personalized-germline settings and results survive session save/restore, and the analysis runs in a lazy cancellable worker so it adds no idle worker or assignment-time overhead.

## Validation

- The complete 254-test suite passes, including regression tests for lowest-current-V-SHM representative selection, one-vote-per-lineage weighting, known two-allele recovery, recurrent linked two-SNP novel-candidate recovery, conservative FASTA generation, context ranking, and interface wording.
- TypeScript checking and the production web/CLI build pass.
- A synthetic 10,000-lineage, 300-nt, two-allele benchmark recovered the exact 50/50 expressed set in about 1.68 seconds on the development machine (about 0.64 seconds ingestion and 1.03 seconds fitting). This is an engineering benchmark, not biological validation.

## Packaging

Browser and CLI version strings are 0.38.5. The source archive includes the GitHub Actions workflow that rebuilds the pinned WebAssembly core and compiles standalone Linux, macOS, and Windows CLI binaries when a `v0.38.5` tag is pushed.
