# Swig 0.36.0

## Reusable custom-reference annotation

- Added `swig-cli prepare-reference`, which performs custom V/D/J/C normalization and germline metadata inference without processing reads.
- Preparation writes per-segment FASTAs containing validated `SWIGMETA`, a portable SHA-256-checked `.swig-reference.json` manifest, and a per-allele `.annotation-diagnostics.tsv` report.
- Added `swig-cli --vdj --prepared-reference MANIFEST`; it verifies every prepared FASTA before worker initialization and restores saved IgBLAST J FWR4 offsets.
- Added strict, permissive, and best-guess metadata-transfer modes, explicit V/J same-gene and nearest-template identity controls, a fallback-candidate limit, and `--require-complete`.
- Best-guess disables identity rejection but retains hard coordinate, frame, V-anchor, and J-motif validation. These settings affect only annotation-schema transfer, never V/D/J read mapping or scoring.
- Detailed diagnostics now identify the selected template and identity or count each rejection cause for unresolved alleles.
- Cached immutable template k-mer sets and deferred broad nearest-neighbour ranking until named-gene candidates fail, removing avoidable repeated work without changing strict-mode decisions.

The legacy spelling `--precompute_aux` is accepted as an alias for the new preparation command, but the output is deliberately called a prepared Swig reference rather than an IgBLAST auxiliary file: IgBLAST `.aux` data describes J metadata only, while the bundle can contain V/D/J/C annotations.
