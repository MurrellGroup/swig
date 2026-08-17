# Standalone CLI, portable pipeline configuration, and lazy lineage studies

**Evidence label:** direct reuse of Swig analysis kernels plus custom execution/interchange engineering.

## Runtime and distribution

`swig-cli` is compiled into one operating-system executable. The release binary embeds the Bun runtime, the same `swiftig.wasm` annotation core used by Swig Web, the browser-WASI compatibility layer, the built-in reference pack, the CLI worker, and the shared TypeScript post-analysis kernels. The person running an analysis does not install Node, npm, Bun, Python, a container engine, or a separate WASM runtime.

Release automation builds Linux x86-64, Linux ARM64, macOS Intel, macOS Apple silicon, and Windows x86-64 binaries. It executes the end-to-end smoke analysis against the actual Linux and macOS artifacts on matching runners; the Windows artifact is cross-compiled on Linux and is not executed by that workflow. Baseline x86-64 targets avoid requiring AVX2. Bun is a build-time packaging choice, not a distinct scientific implementation: annotation and post-analysis call the repository's shared code.

The repository also emits an unpackaged JavaScript CLI during developer tests. That artifact exists to make local regression tests cheap; it is not the documented end-user installation path.

## Configuration and defaults

The configuration is JSON with `schema: 1`. Every block is optional when a default is acceptable. `swig-cli init` writes a complete example. Swig Web exposes **Export CLI config** both before execution in **Analyze → Execution mode → Pipeline** and after execution in Results. It writes:

- exact V/D/J/C FASTA used by the browser run;
- species and receptor/locus scope;
- FASTQ expected-error filtering, 3′ trimming, exact per-dataset subsample size, and base random seed;
- assignment strategy, calling profile, AIRR preserve/reannotate policy, strand, identity floor, worker count, and Double-D settings;
- enabled state and parameters for allele refinement/reassignment, collapse or denoising, CHMMAIRRa, repertoire selection, lineage assignment, SHM, and possible-missing-V evidence;
- input paths and formats, concatenated-gzip member ranges when applicable, plus dataset, biological-sample, donor/subject, cohort, timepoint, and compartment labels; and
- output naming.

Input paths may be edited or replaced for batch use. `subjectId` is the donor grouping key: different files/samples with the same nonempty value are in the same donor. Collapse, allele-pooling, and lineage scopes remain separate explicit settings, exactly as in the browser. `annotation.workers: 0` (the hand-written/default value) selects up to eight host threads; a positive integer is an exact worker count. If no config is supplied, one input may be annotated with `--sample`, `--donor`, `--dataset`, and `--out` flags.

The fixed order is annotation → optional allele fit and policy reassignment → collapse/denoising → optional chimera filter → optional explicit selection → lineage assignment → SHM → possible-missing-V evidence. The CLI deliberately does not infer lineage phylogenies or UCAs because those analyses require choosing and inspecting a lineage MSA.

FASTA, wrapped FASTA, FASTQ, wrapped FASTQ, AIRR TSV/CSV, and gzip-compressed versions are streamed in bounded batches. Annotation uses a bounded worker pool. FASTQ filtering and seeded reservoir sampling use the same shared functions as the browser. A `gzipRange` identifies one complete gzip member by end-exclusive compressed byte offsets, allowing two config inputs to address separate members of one source file. Small data pasted into Analyze are stored in that exported input's `inline` field, making the config self-contained instead of referring to a nonexistent `pasted-sequences.txt`; file inputs remain path based. AIRR outputs include collision-safe dataset/sample/donor fields and retain the source sequence ID; call reassignment occurs before call-dependent collapse and lineage grouping. The resolved config and a machine-readable run summary are always emitted.

`annotation.airrMode` is deliberately explicit. Browser exports use `reannotate`, matching browser upload: the AIRR `sequence` is sent through SwiftIG again and old assignments are not trusted. The default for a hand-written/legacy CLI config is `preserve`, which treats existing AIRR calls as the annotation result. Double-D screening with AIRR requires reannotation because the sparse evidence calculation needs the assigner pass.

## Audited browser-to-CLI option coverage

| Browser pipeline area | Portable config and CLI behavior |
|---|---|
| Input and study structure | Format, source/member name, path or pasted inline content, optional gzip-member range, dataset/sample/donor/cohort/timepoint/compartment, and study design are retained. |
| FASTQ preprocessing | Expected-error threshold, Phred offset, all 3′ trim settings, per-dataset reservoir size, and base seed run through shared browser code. |
| References | Exact composed V/D/J/C FASTA is embedded, including alternative databases, local replacements, and allele exclusions already applied. |
| Assignment | Workers, batch size/default, AER/RIAT-MP/standard strategy, calling profile, identity floor, strand, AIRR mode, and every Double-D threshold are executed. |
| Allele refinement | Model, pooling boundary, segments, weighting, every advanced prior/kernel/inference option, and best-versus-confidence reassignment policy are executed before collapse. |
| Collapse/denoising | Exact key and every FAD/conservative/indel parameter, unresolved policy, scope, candidate cap, and constant-call partition are executed. |
| CHMMAIRRa | Segment, BW/DB/auto choice, exact MSA, prior/mutation states/switch, DFR floor, filter threshold, unevaluated policy, and detailed Viterbi labels are executed. Selected unaligned browser references are aligned with Kalign and embedded during export. |
| Repertoire selection | The complete interactive selection schema—not only the compact pre-run fields—is retained by Results export and applied by the CLI. |
| Lineage, SHM, missing-V evidence | Lineage scope/call policy/identity/productivity/candidate cap, SHM metric, and every two-pass missing-V threshold are executed. |
| Outputs | Directory/prefix, annotated and retained AIRR, method summaries, and optional lazy lineage-study bundle are controlled by config. |

The audit intentionally excludes display-only state, repertoire plots, ad hoc interactive queries, and lineage phylogenetics/UCA. Those are not unattended pipeline stages. The browser still exports the config from Results after they have been used, but the CLI stops at non-phylogenetic repertoire and lineage outputs.

## Lazy lineage-study bundle

When lineage assignment is enabled, the CLI writes two linked files:

1. an uncompressed AIRR TSV sorted by numeric lineage ID; and
2. a compressed lineage-study manifest containing exact references/options, dataset metadata, the bounded lineage summary table, optional lineage SHM mean/q95, and byte ranges into that AIRR file.

Every range stores its row count and SHA-256. The manifest also stores the complete AIRR size and SHA-256. Swig Web's **Load lineage study** first reads only the small manifest and checks the selected linked AIRR size. It can then search/filter the summary table without importing the AIRR rows. When a lineage is selected, the browser reads exactly that byte slice, verifies its row count and range SHA-256, and constructs a temporary result store containing only those rows. The normal lineage alignment, FastTree, and phylogenetic-UCA workbench then operates on that small store.

This route deliberately avoids the full-session metadata-index rebuild. It does not make arbitrary whole-repertoire filtering available: changing to another lineage discards the previous temporary store and reads the new range. The manifest contains at most the same 10,000 largest lineage summaries retained by the normal interactive explorer; assignment itself still covered every eligible record when the CLI ran.

The linked AIRR must be uncompressed because byte offsets into gzip streams are not independently seekable without an additional block index. It is an analysis interchange file, not a replacement for the ordinary final AIRR export.

## Failure and reproducibility behavior

- A malformed or unsupported config is normalized only where a documented default exists; missing input or unusable references fail the run.
- CHMMAIRRa requires an aligned reference MSA. Browser export builds and embeds it before download when selected references are unaligned; an already aligned explicit MSA is retained byte-for-byte.
- Output is written to a dedicated directory; the resolved config records all inherited defaults.
- A lineage range is rejected if its byte count, row count, or SHA-256 does not agree with the manifest.
- The standalone binary has no network dependency during analysis unless a future explicit input option says otherwise.

## Engineering references and relationship

- [Bun standalone executables](https://bun.sh/docs/bundler/executables) provide the runtime embedding, asset embedding, worker bundling, and per-platform compilation used for distribution. Bun does not define any biological method.
- [AIRR Rearrangement schema](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) supplies interoperable field conventions. The JSON config, lineage-sort order, byte-range manifest, range hashes, and lazy browser store are Swig-specific.
