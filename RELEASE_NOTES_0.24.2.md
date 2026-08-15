# Swig 0.24.2

## Phylogenetic UCA codon posterior

- Replaced the amino-acid product-of-nucleotide-marginals approximation with an exact three-column posterior under the existing V(D)J HMM.
- The codon contraction keeps V/J candidate identity, D identity and position, trimming/recombination path, N-region state, and local UCA-placement hypothesis latent until after the three characters are jointly marginalized.
- Added nucleotide, codon, and amino-acid frequency-logo modes. Amino-acid probabilities are obtained by summing synonymous states from the exact 125-state `A/C/G/T/gap` codon posterior.
- Added a long-form codon posterior TSV and included exact codon vectors in the complete result JSON/session state.
- Bound saved UCA results to the selected alignment reading frame so a posterior cannot silently be displayed in a different phase.

## Alignment and gap semantics

- FastTree still omits columns that are terminally missing in every observed sequence.
- The posterior pass now restores the full curated alignment width and treats those terminal columns as missing phylogenetic data. This preserves codon phase and permits the germline/recombination prior to infer their UCA states.
- Leading/trailing gaps remain missing per tip; only internal gaps use the fifth character in automatic GTR5 mode. Complete gap codons translate to `-`; partially gapped codons translate to `X`.

## Verification

- Added a correlated two-V-allele regression showing that the exact codon posterior retains coherent allele codons that an independent-site product disperses.
- Added a marginal-consistency test: summing the 125-state codon posterior over either of the other two positions reproduces the corresponding nucleotide posterior to numerical precision.
- Added codon/AA normalization, frame, complete-gap, split-gap, full-width alignment, and UI regression coverage.
