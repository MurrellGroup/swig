# Swig 0.30.0 release notes

## Browser-to-CLI parity

- Pipeline configs now include FASTQ expected-error filtering, optional 3′ quality trimming, and exact seeded per-dataset reservoir subsampling. `swig-cli` uses the same streaming parsers, filter, and random sampler as Swig Web.
- The Analyze page exposes **Export CLI config** as soon as Pipeline execution mode is selected; an analysis no longer has to be run first. The existing Results-page export remains available.
- Browser exports explicitly reannotate AIRR inputs, matching browser behavior. Hand-written CLI configs retain the backward-compatible `annotation.airrMode: "preserve"` default and can opt into `"reannotate"`.
- Configs retain the original input/member display name separately from its path. Separate samples selected from one concatenated gzip are represented by exact compressed byte ranges and can be reproduced by the CLI without extracting temporary files.
- Small inputs pasted into Analyze are embedded in the exported JSON so the pre-run config is immediately executable instead of pointing at a nonexistent local text file.
- Completed Results exports now retain advanced collapse/denoising, CHMMAIRRa, selection, allele-remapping, lineage, SHM, and missing-allele settings rather than falling back to pipeline defaults.
- CHMMAIRRa configs can retain detailed Viterbi breakpoint precomputation. The CLI writes the corresponding parent and recombination labels when enabled. Browser-selected unaligned references are aligned with Kalign during config export and embedded as the executable MSA.

## Verification

- Regression tests cover wrapped FASTQ QC before subsampling, exact config export, split concatenated-gzip members, metadata-safe sequence IDs, AIRR preserve/reannotate behavior, browser-accepted AIRR CSV input, and completed interactive-state export.
- The upstream v0.29.4 standalone release workflow is retained unchanged: Linux and macOS artifacts are compiled and smoke-tested on matching runners, while the Windows artifact is cross-compiled and published.
