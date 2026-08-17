# Third-party notices

Swig source code is distributed under the repository's MIT license. That
license does not replace the terms that apply to bundled third-party data and
components.

## IMGT/GENE-DB reference data

Swig includes a compact reference pack derived from the official
[IMGT/GENE-DB reference-sequence download](https://www.imgt.org/download/GENE-DB/):

- source release: `202632-7`
- source retrieved: `2026-08-08`
- bundled files: `public/references/imgt-202632-7-swig-0.7.json.gz` and the
  identical CLI asset `cli/assets/imgt-reference-pack.json.gz`

Copyright © 1995-2026 IMGT®, the international ImMunoGeneTics information
system®. IMGT®, the international ImMunoGeneTics information system®, and
its logo are registered marks of the Centre National de la Recherche
Scientifique (CNRS).

IMGT data and metadata are provided under the
[Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
The [IMGT terms of use](https://www.imgt.org/about/termsofuse.php) and
[IMGT citation guidance](https://www.imgt.org/about/CitingIMGT.php) provide the
authoritative current terms and requested citations.

Attribution: IMGT®, the international ImMunoGeneTics information system®,
https://www.imgt.org/, Institute of Human Genetics, Université de Montpellier
and CNRS.

Swig modifies the source data when producing its derived pack. The build:

- selects IG and TR V, D, J, and C records and organizes them by species,
  locus, and segment;
- normalizes nucleotide symbols, converts `U` to `T`, and removes IMGT gap
  characters from indexed sequences;
- retains one selected sequence per allele identifier and stores compact
  coordinate metadata derived from IMGT-gapped V records and J motifs;
- joins selected coding constant-region exons for IGH and TR alleles while
  omitting membrane-only and untranslated exons; and
- serializes the result as a compressed, Swig-specific JSON structure.

These modifications are made by the Swig project. IMGT, Université de
Montpellier, and CNRS do not endorse Swig or warrant the modified pack or its
use.

Please cite IMGT/GENE-DB as requested by its maintainers:

> Giudicelli V, Chaume D, Lefranc M-P. IMGT/GENE-DB: a comprehensive database
> for human and mouse immunoglobulin and T cell receptor genes. Nucleic Acids
> Research. 2005;33:D593-D597. https://doi.org/10.1093/nar/gki010

## Other components

Other bundled third-party components retain their own adjacent license or
notice files, including `DEJAVU_FONT_LICENSE.txt` and the notices under
`public/tools/` and `public/references/kimdb-1.1/`.
