# Germline reference preparation

## Reference composition

Swig composes a locus-by-segment matrix for V, D, J, and C. A cell may come from the bundled IMGT-derived pack, a named online collection, or a local FASTA. Unchanged cells retain their recommended source, so replacing one segment does not implicitly replace the others.

The UI reports the exact number of active V/D/J/C records before a run. Per-cell allele exclusions filter FASTA records by their exact FASTA identifier. Exclusion changes the reference supplied to initial assignment and consequently every later method; CHMMAIRRa does not maintain a second independent allele-exclusion list.

## Parsing and validation

Reference FASTA is parsed into named records, normalized to uppercase nucleotide sequences, and checked for empty/duplicate-invalid content. Segment preparation computes the metadata and indexes needed by SwiftIG in a worker. Online choices fetch only the chosen public FASTA; query data are not uploaded.

The standalone CLI uses the same parsing, normalization, validation, taxonomic template tiers, nucleotide coordinate projection, and conserved-anchor checks for local FASTA supplied through `references.files`. This metadata preparation is enabled by default with `references.prepareMetadata: true`; setting it to `false` deliberately passes the local FASTA through without metadata transfer. The bundled fixed reference pack supplies the coordinate donors, so CLI preparation has no network dependency. Exact FASTA embedded in `references.inline` by Swig Web is already prepared and is retained byte-for-byte rather than compared again.

Exact sequence-identical records may remain distinct labels for assignment output. The later repertoire allele model collapses sequence-identical labels into one unresolved inference node because no nucleotide evidence can distinguish them; this does not alter the initial database.

Reference preparation is transactional at the UI level. A whole-database choice is displayed only after all required cells have downloaded, parsed, and prepared successfully. Cancelling terminates preparation/fetches and retains the previously committed matrix. Per-cell preparation follows the same rule.

## Coordinate preparation

SwiftIG computes local affine alignments against the selected references and reports AIRR coordinates, alignments, CIGAR strings, identities, and region annotations. Swig uses IMGT-style region fields when they can be mapped from the reference annotations; it does not claim to run IMGT/V-QUEST. See the IMGT unique numbering overview: [Lefranc et al., 2003](https://www.imgt.org/textes/PDF/DCI/27_55-77_2003.pdf).

## Literature relationship

This reference matrix and browser preparation layer are **custom Swig infrastructure**. IMGT, KI databases, and uploaded FASTA are data sources, not algorithms implemented by Swig. Personalized germline databases can materially improve assignment, as demonstrated by IgDiscover ([Corcoran et al., 2016](https://doi.org/10.1038/ncomms13642)) and TIgGER ([Gadala-Maria et al., 2015](https://doi.org/10.1073/pnas.1417683112)); selecting or excluding existing references in Swig is not equivalent to either method's novel-allele/genotype inference.

## Limitations

- Reference names and sequences are trusted after format validation; biological curation remains the user's responsibility.
- Public coverage differs by species and locus.
- Removing a true allele can force assignments to a neighbour; retaining hundreds of impossible alleles can increase ambiguity. Repertoire refinement can redistribute evidence only among retained reference nodes.
- Reference source URLs and collection versions should be recorded with exported work when external collections are used.
