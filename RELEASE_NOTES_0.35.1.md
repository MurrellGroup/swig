# Swig 0.35.1 release notes

## N-aware lineage MSA

- The Alivibe-compatible WASM nucleotide aligner now gives `N` the full match score against any non-gap residue in its progressive POA and pairwise refinement dynamic programs.
- Gap opening and extension are unchanged. Exact k-mer anchors also remain literal; ambiguity-containing windows fall through to the wildcard-aware DP without adding a slower wildcard seed index.
- Protein mode retains the original literal score, so amino-acid `N` continues to mean asparagine rather than a wildcard.
- The original WASM entry point remains available for literal/protein scoring; Swig's nucleotide routes call a separate nucleotide entry point.
- The bundled Alivibe bridge version is incremented so an already-open editor from an older release cannot accidentally apply nucleotide wildcard scoring to amino-acid mode.

## Compact lineage plots

- The lineage MSA preview can hide sequence identifiers and immediately reassign the removed grid column to aligned residues.
- Every coordinated lineage tree, including the UCA tree, can hide tip identifiers. The SVG sequence matrix moves left instead of retaining the former 208 px label gutter.
- Tip and residue hover text remains complete.

## Verification

- Non-`N` inputs remain byte-identical to the pinned Alivibe JavaScript oracle.
- A deterministic germline-like ambiguity/indel fixture verifies wildcard-`N` placement while preserving every ungapped input sequence.
- UI regressions verify that both name controls collapse their corresponding layout gutters.
