# Swig 0.11.0

## Method D: indel-aware denoising

- Adds a fourth explicit collapse method without changing methods A–C.
- Partitions exact-dereplicated trimmed VDJ sequences by locus and configurable V/J call resolution.
- Uses a complete length-aware `d+1` segment join for edit radii one or two. Candidate generation scales with accepted variants rather than enumerating deletion neighborhoods.
- Verifies candidates with a reusable narrow-band dynamic program in `O(sequence length × edit radius)` memory-bounded work and records substitutions, insertions, and deletions on an exact optimal path.
- Collapses indel-containing children only toward a strictly more abundant accepted parent that meets the configurable parent:child ratio. Substitution-only paths use method C’s sequence-specific Poisson rule.
- Never assigns a read to a distant centroid. Candidate-cap truncation remains explicit, and indel/substitution merge counts are reported separately.
- Preserves summed `duplicate_count` through exports, lineage abundance, and phylogeny tip bubbles.

The bundled 50,000-record indel benchmark (`npm run benchmark:indel-denoise-50k`) completed in 246 ms in the final verification run (about 203,000 input records/s), performed 1,000 exact candidate verifications, and measured an 8 MiB V8 heap increase. Timing and garbage-collector-dependent heap deltas vary between runs; this is a regression benchmark, not an accuracy estimate.

## Tree/alignment corrections

- Branch labels now follow the active sequence mode. Nucleotide mode shows reconstructed nucleotide changes. Amino-acid mode translates each reconstructed parent/child codon, reports one nonsynonymous replacement per changed codon, and omits synonymous changes and unresolved `X` codons.
- Non-contiguous displayed alignment runs, including the all-CDR preset, receive a half-residue blank separator in the live tree/alignment view and exported SVG.
- Lineage alignment cells now occupy the complete row height with no vertical space between adjacent sequences.
- Custom alignment positions and motifs use local edit drafts and are applied only on Enter or blur. Numeric controls throughout results/post-analysis use the same deferred-commit behavior, so transient values such as an empty field or `0.` are no longer resolved while typing.

## Verification

- Adds exhaustive one-edit insertion/deletion boundary coverage, mixed two-edit collapse and abundance-protection tests, exact edit-profile tests, synonymous-mutation suppression tests, and discontinuous-column spacing tests.
- The complete test suite and production TypeScript/Vite build pass in the packaged source.
