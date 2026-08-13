# Swig 0.19.0

This release adds explicit reference-alignment composition, candidate-reference exports, expanded lineage-tree metadata coloring, constant-aware collapse defaults, and truthful staged denoising progress.

## Reference alignments and missing-allele candidates

- CHMMAIRRa now has an expandable exact-identifier allele exclusion editor for both the run-composed reference FASTA and a loaded aligned FASTA.
- Search, per-allele checkboxes, and bulk include/exclude actions operate before Kalign or MSA validation. At least two retained records are required.
- Pipeline mode exposes the same exclusions up front. Interactive and pipeline selections are retained by portable sessions and project checkpoints.
- Excluded assigned alleles remain explicitly unevaluated/missing; they are never replaced by a nearby reference.
- The missing-V diagnostic now exports all candidates or only selected candidates as FASTA.
- Selected candidates can be appended to a downloadable copy of the exact V reference used by the run. Every original reference remains; candidate headers record diagnostic status, parent, and substitutions. `SWIGMETA` is inherited only from the exact parent record.
- Candidate selection is retained in saved sessions. Constructing an augmented FASTA does not silently change or rerun the completed analysis.

## Constant-aware collapse

- The default exact-collapse key is now the V–J-trimmed nucleotide sequence.
- Exact collapse and denoising separately partition on the normalized top C-gene/isotype call by default. Different called constant genes therefore do not merge, while random constant-tail sequence length is not part of the comparison key.
- Records without a C call form a separate unassigned partition.
- An explicit toggle in interactive and pipeline modes permits collapse while ignoring constant assignments.
- Existing/imported `duplicate_count` remains multiplicity-weighted after every policy.

## Denoising progress and allocation behavior

- The former apparent stall after streamed ingestion is removed. Ingestion, indexed-variant denoising, and representative/multiplicity materialization are now separate progress stages.
- The worker emits bounded variant-level and finalization progress messages while computation is running; the progress bar resets for each real stage rather than reaching 100% before the main algorithm begins.
- Hot candidate sets are reused across variants to reduce large-run allocation and garbage-collection pressure without changing candidate completeness or scoring.
- The included 50k method-D benchmark completed at approximately 191k records/s with an 18 MiB heap delta in the release environment; this synthetic result is not a guarantee for other repertoire structures or browsers.

## Lineage tree metadata

- Multiplicity bubbles can now be colored by sample, original lineage, isotype, constant gene, donor/subject, cohort, timepoint, compartment, V gene, J gene, productivity, or double-D status, or drawn uniformly.
- Sample mode keeps the centrally editable sample-color assignments. Other metadata use stable deterministic categorical colors.
- The active tip-color legend is included inside the downloaded coordinated tree/alignment SVG.

## Validation

- Added regression coverage for exact FASTA exclusions and metadata preservation, selected-candidate FASTA and augmented references, C-aware exact collapse and denoising, worker progress phases, pipeline defaults, and stable categorical colors.
- Full TypeScript, production-build, post-analysis, storage, streaming, germline, alignment, and WebAssembly suites are run before packaging.
