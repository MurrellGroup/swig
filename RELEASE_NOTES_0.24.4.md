# Swig 0.24.4

## Exact frequency-logo geometry

- Replaced em-box text scaling with measured SVG ink-bound fitting.
- Removed the two-pixel column inset and all deliberate inter-letter spacing.
- Removed vertical font-box padding from every probability rectangle.
- Switched logo glyphs to a bold monospace stack so narrow amino-acid symbols remain visible.
- Descenders, including the tail of `Q`, are included in the fitted bounds rather than clipped.
- The reusable nucleotide, codon, and amino-acid logo still normalizes every stack to one and never applies entropy scaling.

## Phylo-HMM source annotation

- Added a compact HMM annotation surface directly above the UCA posterior logo on the same horizontal grid and scroll position.
- **Best path** shows the Viterbi recombination path at the single best tree placement.
- **Marginalized** shows forward-backward source occupancy mixed over retained attachment and pendant-length hypotheses.
- V, D, and J rows are generated from the phylo-HMM state, never copied from the input AIRR calls.
- Allele rows remain pure template tracks. D register uncertainty is split into separate allele/register rows rather than mixed within an allele.
- N and trimming-boundary rows show character mixtures whose total letter height equals their unnormalized posterior source occupancy.
- Non-negligible allele groups are retained with explicit visualization-only thresholds; omitted subthreshold rows are reported.
- Nucleotide subcolumns remain aligned above nucleotide, codon, and amino-acid logo modes.
- Codon display does not quantize HMM occupancy or recombination boundaries: V/D/J/N transitions remain possible after any nucleotide, including positions one and two within a displayed codon.

## Verification

- Added statistical invariants for Viterbi source coverage, pure allele rows, HMM-derived V/D/J/N tracks, and schema-3 session state.
- Added UI regressions for ink-bound measurement, zero column padding, the `Q` descender safeguard, shared fitted glyphs, and both annotation modes.
