# Swig 0.24.1

## Phylogenetic UCA posterior visualization

- Added a reusable fixed-height probability-logo component with clipped glyph bounds, serif lettering, horizontal scrolling, and self-contained SVG export.
- Added nucleotide UCA posterior logos showing `A/C/G/T/gap` marginal probabilities without entropy scaling.
- Added an optional amino-acid projection using the current alignment reading frame. The view explicitly reports its product-of-nucleotide-marginals approximation; complete gap codons remain `-` and split-gap codons become `X`.
- UCA posterior vectors are normalized after HMM inference and again after local-placement mixing, with a floating-point residual correction so every displayed/exported site sums to one.

## Gap observation semantics

- Leading and trailing gap runs are now treated as missing sequence coverage for each observed tip.
- Only gaps between a tip's first and last observed nucleotide/IUPAC character are exact fifth-state observations.
- Automatic character selection now activates GTR5 only when an internal gap exists; alignments containing terminal padding alone use GTR4.
- Added explicit result warnings reporting terminal missing gaps and internal fifth-state gaps.
- Added regression tests for terminal, internal, forced-GTR4, posterior normalization, codon projection, serif font selection, and clipped logo bounds.
