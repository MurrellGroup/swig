# Swig 0.14.0

## Multi-dataset processing

- Upload multiple FASTA, FASTQ, or AIRR datasets in one analysis.
- Assign editable dataset, biological sample, donor/subject, cohort, and timepoint metadata before the run.
- Combine results into one collision-safe AIRR table with `sample_id`, `subject_id`, and explicit `swig_*` study fields.
- Filter the result index and downstream working set by every study field.
- Stratify SHM summaries by sample, donor, cohort, or timepoint.

## Hard biological processing boundaries

- Exact collapse and denoising can be scoped independently to uploaded dataset, biological sample, donor, cohort, or the complete study.
- Sample scope is the default. Technical replicate files collapse together only when they share a sample ID.
- Lineage assignment has a separate scope. The longitudinal preset uses donor scope, allowing timepoints from the same donor to share a lineage without comparing different donors.
- Scope values are embedded into the compact candidate indexes, so prohibited cross-boundary comparisons are not generated.

## Pipeline mode

- Configure repertoire-scale stages before annotation and run them unattended in a fixed order.
- Supported automatic stages are collapse/denoise, CHMMAIRRa exclusion, repertoire selection, lineage assignment, SHM, and lineage-aware missing-allele hints.
- Every stage consumes the cumulative retained-record mask produced by the preceding stage.
- Up-front controls include exact-collapse key, collapse boundary and unresolved policy; CHMMAIRRa model, MSA and threshold; dataset/sample/donor/cohort/timepoint plus call/CDR3/QC selection; lineage boundary, identity, call resolution and ambiguity policy; and the SHM metric.
- The study manifest, pipeline options, masks, collapse multiplicities, and lineage assignments are included in saved-session state.

Single-input analysis remains supported. Multi-dataset execution is sequential by dataset and parallel within each dataset, so the configured WebAssembly worker count and bounded-memory guarantees are not multiplied by the number of inputs.
