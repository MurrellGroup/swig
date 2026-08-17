# Swig 0.35.2 release notes

## Frame-preserving Alivibe MSA

- The lineage workbench retains the existing nucleotide **MSA lineage · Alivibe WASM** action and adds **AA MSA · Alivibe WASM**.
- AA MSA translates each ungapped sequence from its AIRR-derived frame, maps every ambiguous codon—including germline codons containing `N`—to `X`, aligns the amino acids, and back-translates each amino-acid gap as one `---` codon.
- The Alivibe-compatible WASM core has a separate amino-acid entry point in which `X` is a symmetric unknown-residue wildcard. Asparagine `N` remains literal and all gap costs are unchanged.
- The bundled Alivibe bridge now selects nucleotide-`N` or amino-acid-`X` scoring from its active mode.

## Productive-row filtering

- **Drop non-productive** reversibly excludes false or unresolved AIRR productivity rows from the current lineage workspace.
- Excluded rows cannot contribute to germline construction, either MSA route, tree metadata, or UCA inference. Repertoire records and lineage assignments are not mutated.
- All-row and productive-only manual alignments are saved independently in sessions.

## Deleted-row UCA integrity

- Alivibe row deletion already removed the sequence from FastTree and the phylogenetic likelihood.
- UCA preparation now also intersects AIRR metadata with the exact retained FASTA names, closing a residual path through which a deleted row could still affect germline-candidate screening or segment-boundary preparation.

## Verification

- Differential tests retain byte-identical literal-scoring behavior against the pinned JavaScript implementation.
- Deterministic fixtures verify nucleotide `N` and amino-acid `X` wildcard placement, literal protein `N`, ambiguous-codon translation, triplet back-projection, productive-row classification, and deleted-row metadata exclusion.
