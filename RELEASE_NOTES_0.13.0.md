# Swig 0.13.0

## Explicit, composable repertoire populations

Post-analysis now has a separate preview-and-commit selection stage after collapse/chimera filtering and before lineage-scale analyses. Filters compose sequence ID, locus, V/D1/D2/J/C/isotype calls, CDR3 nucleotide or amino-acid substrings, motifs (literal, IUPAC, or regular expression), productivity/QC states, CDR3 lengths, V/J identity, V mutation fraction, and double-D status. Text drafts commit on Enter or blur and do not repeatedly scan while being edited. Positive double-D mode resolves only the sparse evidence ordinals. Lineage assignment, SHM, missing-allele diagnostics, query, and expansion inherit the committed mask.

## Double-D exploration

Supported VDDJ records have a dedicated results tab with D1/D2, score-gain, V–J-span, and sequence-ID filters; a D1→D2 pair summary; and an on-click V–D1–D2–J alignment in nucleotide or amino-acid mode. The sequence explorer also exposes a prominent indexed “only double-D positive” switch. The ordinary AIRR annotation remains unchanged when screening is off, and double-D evidence remains a sparse sidecar.

## SHM and possible missing-allele evidence

The selected population can be summarized by V nucleotide mutation count/rate, V amino-acid replacement count/rate, synonymous codons, CDR1/2 nucleotide rate, or framework nucleotide rate. Duplicate counts weight the figures. SVG outputs include lineage violins, repertoire histograms, lineage abundance/SHM bubbles, V-gene summaries, and optional locus/V/isotype strata; underlying data export separately.

The missing-allele diagnostic uses one least-mutated observation per lineage and V call by default. Recurrent linked substitutions must satisfy configurable support, coverage, fraction, low-SHM, tail-probability, J-diversity, and junction-length-diversity rules. Known reference sequences are excluded and common SHM hotspot contexts are flagged. Results are explicitly warnings/candidate FASTA for follow-up with IgDiscover, TIgGER, partis, or another genotype-aware workflow; Swig never modifies the active germline set automatically.

## Exports and linked sessions

- AIRR, collapse, lineage, CHMMAIRRa, active-population, double-D, SHM, and candidate-evidence tables export as TSV, CSV, or JSON Lines where applicable.
- Repertoire and post-analysis figures export as standalone SVG plus CSV data.
- Lineage alignments export as aligned FASTA, Clustal, relaxed PHYLIP, Stockholm, or NEXUS; trees export as Newick or rooted NEXUS in addition to the annotated SVG.
- **Save session** writes gzip-compressed JSON containing exact composed references, analysis settings, sparse double-D evidence, typed masks/counts/assignments, query state, alignment/tree state, and downstream dashboards. It does not embed the AIRR rows. **Load session** requests the linked AIRR TSV/TSV.gz and verifies its columns, record count, and streaming 128-bit fingerprint before restoring analysis state.

## Scale and verification

Converted tabular outputs use streaming writers when the browser exposes the File System Access API. Double-D-positive selection uses its dedicated IndexedDB index. SHM plot/session retention is bounded globally and per lineage. Missing-allele observations store compact mutation lists and coverage intervals rather than full aligned rows.

Verification includes a complete TypeScript check, production build, the existing 50,000-record AIRR store/export regression, the existing 50,000-record post-analysis regression, and new tests for standard serializers, sparse double-D selection, SHM codon accounting, lineage-aware candidate evidence, and linked-session typed-vector round trips.
