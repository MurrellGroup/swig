# Swig 0.20.2

## Allele exclusion before initial assignment

- Every applicable locus/segment cell in the initial reference-composition
  matrix now has a visible **Exclude alleles…** action.
- The editor searches exact FASTA identifiers, supports individual toggles and
  bulk exclusion/inclusion of the current search, and reports retained and
  excluded record counts.
- Exclusions are applied independently per V/D/J/C cell while composing the
  FASTA passed to `compileReferences`; SwiftIG therefore never indexes or calls
  removed alleles as primary or alternative hits. Other loci and segments are
  unaffected.
- IMGT, compatible published databases, and loaded custom FASTA all use the
  same exclusion path. Changing a cell's source clears that cell's old
  exclusions rather than transferring names to an unrelated database.
- The resulting exact compiled references and the per-cell exclusion map are
  retained in saved sessions and project checkpoints.

## Chimera-analysis boundary

- Removed the independent allele-exclusion UI, pipeline field, filtering call,
  and saved exclusion state from CHMMAIRRa.
- CHMMAIRRa can build an MSA from the immutable references already used for
  initial assignment, or accept a complete aligned FASTA file. It no longer
  defines a second allele-filtering stage.

Regression coverage verifies that excluding one IGH V allele removes exactly
that record from the combined BCR assignment FASTA, preserves IGK/IGL and J,
and survives session serialization.
