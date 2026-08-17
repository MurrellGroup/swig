# Swig 0.31.0 release notes

## Alivibe MSA WebAssembly port

- The lineage alignment view now has a one-click **MSA lineage · Alivibe WASM** action; it runs entirely inside Swig, without an Alivibe round trip, in a dedicated cancellable worker backed by a 90-KiB optimized WASM module.
- The C++ port preserves the pinned MurrellGroup/WebWidgets `refinedMSA` behavior: ordered progressive POA, unique 15-mer seeding, LIS and anchor erosion, strict tie priorities, graph column packing, three 20-column-flank refinement passes, Double-DP pairwise alignment, and the existing amino-acid-to-codon projection.
- The optional Alivibe open/edit/return workflow remains available for manual curation and uses the same WASM core. Standalone Alivibe retains its original JavaScript implementation as a fallback.
- Differential tests compare zero through three refinement passes byte-for-byte on nucleotide, amino-acid, indel, ambiguity, order-permutation, randomized, and low-complexity fixtures. The suite explicitly covers JavaScript's reversed-`substring` behavior for overlapping anchors.
- The MSA task no longer blocks the editor UI and can be cancelled before any aligned state is committed.

On a deterministic 80-row × 360-nt benchmark, the warm WASM core took 59.610 ms versus 2,230.370 ms for the pinned JavaScript implementation (37.42×). At 200 rows × 360 nt it took 216.080 ms versus 6,357.812 ms (29.42×). The benchmark asserts exact output equality before reporting timing.

## Repository attribution cleanup

Duplicate third-party-style license copies and notices for MurrellGroup-owned WebWidgets/Alivibe and SwiftIG components were removed. Source-revision provenance remains. Notices and terms for unrelated external dependencies and data sources remain unchanged.
