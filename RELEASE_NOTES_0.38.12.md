# Swig 0.38.12

## Clipped V junction anchors

A local V alignment can stop before the conserved cysteine because clustered terminal substitutions lower its score. Junction annotation previously required that anchor to lie inside the committed local alignment, so an intact coding sequence could have no CDR3/productivity annotation. The detail panel then mislabeled its empty productivity field as nonproductive.

The annotation layer now checks a short, reference-guided co-linear mapping to the selected allele's annotated cysteine. It requires an intact query cysteine at that projected coordinate, unambiguous compared bases, a tail of at most 18 bases, and no overlap with the selected D/J tract. An anchored affine comparison includes up to 12 retained co-linear bases and alternative query endpoints. Any equally good or better gapped explanation prevents rescue. It does not search for a convenient cysteine elsewhere in the junction or use an in-frame result to select a mapping.

The selected allele calls, segment scores, alignment strings, boundaries, support values and NP decomposition remain unchanged. Stop-codon and V/J frameshift vetoes still apply. The rescue is limited to clipped V cysteine anchors; unresolved J anchors and ambiguous or longer V extensions remain unresolved. An observed cysteine's structural coordinate does not prove that all its bases came from V: junction addition can recreate a deleted codon.

Empty AIRR productivity now appears as **Productivity unresolved** in sequence details, with a neutral table marker. Explicit `F` calls still appear as nonproductive.

## Validation

- The supplied alpaca `Dyl` sequence, using the complete bundled Vicugna pacos IGH reference and AER-R/R-optimized, now has `productive=T`, `cdr3_start=289`, `cdr3_end=339`, and CDR3 `AARFAGIVAGPWDEYNY`. Its V call remains `IGHV3S66*01`, V score 452, and V endpoint 276. Its reverse complement gives the same annotation.
- Regression cases retain stop-codon, junction-frame and V-frameshift vetoes. Deleted/misplaced anchors, unknown clipped bases, and competitive gapped mappings remain unresolved.
- All 290 tests pass; the final strengthened indel test and production build pass. The rebuilt CLI also reproduces the Dyl result.
- A paired before/after comparison used 3,000 simulated alpaca sequences and 3,000 simulated human sequences (seed 765127, no indels or sequencing errors, with SHM, exonuclease trimming and junction additions), plus 2,000 sampled sequences from each preserved macaque full-reference AIRR input (reservoir seed 75612). All runs used AER-R/R-optimized. All non-annotation fields and all previously resolved annotations were identical across the 10,000 sequences.
- The simulation rescued 12 alpaca and 3 human anchors; every rescued coordinate agreed with the simulator's projected germline coordinate. One human rescue became productive. Eleven rescues had a codon partially or wholly recreated by junction addition; all eleven remained nonproductive. These checks validate structural coordinate recovery, not exact germline provenance or biological function.
- Neither macaque sample had a newly rescued anchor; their complete AIRR rows were unchanged. They validate assignment stability, not additional rescue sensitivity.

Reproduce the paired check against an unmodified 0.38.11 binary with:

```sh
node tests/benchmark-junction-anchor.mjs /path/to/0.38.11/swiftig.wasm
```

The small rescue runs only when the original mapped V anchor is missing. Existing assignment kernels and candidate searches are unchanged.
