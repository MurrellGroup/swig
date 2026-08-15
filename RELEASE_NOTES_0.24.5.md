# Swig 0.24.5

## Exact painted-contour frequency logos

- Replaced live SVG text measurement with embedded DejaVu Sans Mono Bold outline paths.
- Fits each glyph's literal contour extrema to its assigned probability rectangle, eliminating the vertical gaps left by logical font bounds.
- Preserves true descender extents, including the tail of `Q`, without clipping or adding font-box padding.
- Uses the same path-based fitter for nucleotide, codon, amino-acid, and aligned phylo-HMM source tracks.
- Keeps the application fully static and compatible with GitHub Pages deployment.

## Verification

- Added regression coverage for every supported nucleotide and amino-acid outline, multi-character codon runs, contour geometry, and removal of `getBBox()` measurement.
- Verified the complete automated test suite and production build.
