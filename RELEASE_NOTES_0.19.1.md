# Swig 0.19.1

## Alivibe roundtrip and amino-acid frame

- Swig now retains the shared AA reading-frame offset selected in Alivibe. The compact lineage alignment and coordinated tree/alignment viewer translate from that same nucleotide column.
- The injected **Return alignment to Swig** action records Alivibe's frame selector, switches Alivibe to NT before extracting FASTA, and therefore cannot accidentally return the current AA view as a nucleotide alignment.
- Downloaded/imported corrected FASTA has no frame metadata, so Swig infers the phase by minimizing mixed base/gap codons across the MSA, with stop and ambiguity counts as tie-breakers. The resulting frame is visible and manually selectable as nucleotide column 1, 2, or 3.
- Frame selection is saved with manually edited alignments and with the current session, and it propagates to Kabat numbering, region projection, AA parsimony labels, and SVG tree/alignment output.
- Translation semantics remain strict: only `---` in the selected codon phase is displayed as `-`; a mixed base/gap codon remains `X`.
- Regression coverage verifies that a genuine complete gap codon in Alivibe frame 2 remains `-`, while a split codon remains `X`.
