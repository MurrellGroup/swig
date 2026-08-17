# Standalone CLI, portable pipeline configuration, and lazy lineage studies

**Evidence label:** direct reuse of Swig analysis kernels plus custom execution/interchange engineering.

## Runtime and distribution

`swig-cli` is compiled into one operating-system executable. The release binary embeds the Bun runtime, the same `swiftig.wasm` annotation core used by Swig Web, the browser-WASI compatibility layer, the built-in reference pack, the CLI worker, and the shared TypeScript post-analysis kernels. The person running an analysis does not install Node, npm, Bun, Python, a container engine, or a separate WASM runtime.

Release automation builds Linux x86-64, Linux ARM64, macOS Intel, macOS Apple silicon, and Windows x86-64 binaries. It executes the end-to-end smoke analysis against the actual Linux and macOS artifacts on matching runners; the Windows artifact is cross-compiled on Linux and is not executed by that workflow. Baseline x86-64 targets avoid requiring AVX2. Bun is a build-time packaging choice, not a distinct scientific implementation: annotation and post-analysis call the repository's shared code.

The repository also emits an unpackaged JavaScript CLI during developer tests. That artifact exists to make local regression tests cheap; it is not the documented end-user installation path.

## Configuration and defaults

The configuration is JSON with `schema: 1`. Every analysis block is optional when a default is acceptable; `output.directory` is deliberately required unless `--out` supplies it. `swig-cli init` writes a complete example. Swig Web exposes **Export CLI config** both before execution in **Analyze → Execution mode → Pipeline** and after execution in Results. It writes:

- exact V/D/J/C FASTA used by the browser run;
- species and receptor/locus scope;
- FASTQ expected-error filtering, 3′ trimming, exact per-dataset subsample size, and base random seed;
- assignment strategy, calling profile, AIRR preserve/reannotate policy, strand, identity floor, worker count, and Double-D settings;
- enabled state and parameters for allele refinement/reassignment, collapse or denoising, CHMMAIRRa, repertoire selection, lineage assignment, SHM, and possible-missing-V evidence;
- input paths and formats, concatenated-gzip member ranges when applicable, plus dataset, biological-sample, donor/subject, cohort, timepoint, and compartment labels; and
- output naming.

Input paths may be edited or replaced for batch use. `subjectId` is the donor grouping key: different files/samples with the same nonempty value are in the same donor. Collapse, allele-pooling, and lineage scopes remain separate explicit settings, exactly as in the browser. `annotation.workers: 0` (the hand-written/default value) selects up to eight host threads; a positive integer is an exact worker count. `--workers N` is applied after reading JSON and therefore always overrides `annotation.workers`. If no config is supplied, one input may be annotated with `--sample`, `--donor`, `--dataset`, and `--out` flags. Every pipeline run must name an output directory explicitly through `output.directory` or `--out`; normalization no longer turns an omitted destination into a silent working-directory write.

Local replacement FASTA in `references.files` is prepared before SwiftIG indexing by default. The CLI calls the same germline preprocessor and species-first taxonomic template tiers as Swig Web, using the embedded fixed reference pack to transfer and validate V FWR/CDR boundaries and J frame/CDR3-anchor metadata. D/C files receive the same parsing, locus/segment checks, normalization, and metadata validation. No reference or read is sent over the network. Set `references.prepareMetadata` to `false` only when the file must be passed through without metadata transfer:

```json
{
  "references": {
    "species": "Homo sapiens",
    "scope": "IGH",
    "prepareMetadata": true,
    "files": {"V": "custom-V.fasta", "J": "custom-J.fasta"}
  }
}
```

Swig Web exports already-prepared exact FASTA through `references.inline`; those strings are not reprocessed. If both an inline and file value exist for one segment, the file remains the replacement and follows the `prepareMetadata` setting. Preparation counts and warnings are written to stderr and the run summary, while the resolved config records the effective flag.

The fixed order is annotation → optional allele fit and policy reassignment → collapse/denoising → optional chimera filter → optional explicit selection → lineage assignment → SHM → possible-missing-V evidence. The CLI deliberately does not infer lineage phylogenies or UCAs because those analyses require choosing and inspecting a lineage MSA.

FASTA, wrapped FASTA, FASTQ, wrapped FASTQ, AIRR TSV/CSV, and gzip-compressed versions are streamed in bounded batches. Annotation uses a bounded worker pool. The annotated AIRR header is written after the first completed input batch and each subsequent batch body is appended in input order; preserved AIRR inputs therefore may not introduce a new column after that header has been committed. The retained table is serialized in bounded row chunks after global policy/filter decisions instead of constructing one repertoire-sized string. FASTQ filtering and seeded reservoir sampling use the same shared functions as the browser. A `gzipRange` identifies one complete gzip member by end-exclusive compressed byte offsets, allowing two config inputs to address separate members of one source file. Small data pasted into Analyze are stored in that exported input's `inline` field, making the config self-contained instead of referring to a nonexistent `pasted-sequences.txt`; file inputs remain path based. AIRR outputs include collision-safe dataset/sample/donor fields and retain the source sequence ID; call reassignment occurs before call-dependent collapse and lineage grouping. The resolved config and a machine-readable run summary are always emitted. Repertoire-wide enabled methods still retain the compact/full row state those methods actually consume; streaming the file writer does not turn a global allele, deduplication, chimera, or lineage calculation into an online algorithm.

`annotation.airrMode` is deliberately explicit. Browser exports use `reannotate`, matching browser upload: the AIRR `sequence` is sent through SwiftIG again and old assignments are not trusted. The default for a hand-written/legacy CLI config is `preserve`, which treats existing AIRR calls as the annotation result. Double-D screening with AIRR requires reannotation because the sparse evidence calculation needs the assigner pass.

## Assignment-only `--vdj` route

`swig-cli --vdj` is a separate execution route for high-volume assignment. It does not call `annotateAirrBatch`, construct study manifests, parse results into post-analysis records, retain repertoire rows, or run allele pooling/collapse/chimera/selection/lineage/SHM/missing-allele code. The query parser supplies one bounded FASTA batch to the same worker pool, completed promises are consumed in submission order, and each AIRR body is written immediately after its single header. `-out`/`--out` is mandatory; `-out -` is an explicit streamed stdout destination. `-query -` reads FASTA from stdin. `-outfmt` currently accepts only `19` because other IgBLAST renderers are not SwiftIG output formats.

The route accepts the common IgBLAST names `-query`, `-germline_db_V`, `-germline_db_D`, `-germline_db_J`, `-c_region_db`, `-out`, `-outfmt`, `-organism`, `-ig_seqtype`, `-domain_system`, `-strand`, `-num_threads`, `-custom_internal_data`, `-auxiliary_data`, and `-d_frame_data`. Source germline arguments are nucleotide FASTA (plain or gzip), not BLAST database prefixes. `--workers` has higher precedence than `-num_threads`. SwiftIG-specific `--batch-records`, `--minimum-identity`, `--assigner`, and `--calling-profile` remain available. A small set of D/J search controls maps to the already exported SwiftIG tuning ABI; unsupported BLAST scoring/search/report options fail rather than being silently ignored.

Reference/annotation behavior is intentionally three-way:

| Invocation | Metadata used by SwiftIG |
|---|---|
| no annotation file and no `--swigannots` | Headers are canonicalized, IMGT `.` gaps (and any `-` gaps) are removed, valid IUPAC bases are retained, and existing `SWIGMETA` is removed. Segment assignment/alignment/endpoints are emitted; FWR/CDR/junction fields remain empty. |
| `-custom_internal_data` and/or `-auxiliary_data` | V `.ndm.imgt` entries are exact-identifier matches with 1-based inclusive FWR/CDR intervals converted to 0-based half-open SwiftIG bounds. J `.aux` frame/CDR3-stop fields remain 0-based. Its optional fifth field is used to map the FWR4 end before each batch is written. A custom V file must cover every selected V identifier, matching IgBLAST's documented custom-data requirement. If J auxiliary data are given without custom V data, the embedded fixed species pack supplies validated V metadata in place of IgBLAST's installed organism directory. |
| `--swigannots` | The same species-first, progressively broadened and anchor-validated preprocessor used for a Swig Web upload. This is mutually exclusive with the two V/J IgBLAST annotation inputs. |

`-d_frame_data` can additionally write exact 0-based D frame-one starts. Auxiliary/internal chain labels are checked against the locus inferred from each germline identifier. Duplicate `.aux` identifiers follow IgBLAST's map behavior (the last record wins); the official human file itself contains a duplicate. The route does not claim to reproduce IgBLAST's BLAST search, result ranking, or non-AIRR formatting. It reuses IgBLAST's public file/interface conventions to initialize SwiftIG metadata.

## Audited browser-to-CLI option coverage

| Browser pipeline area | Portable config and CLI behavior |
|---|---|
| Input and study structure | Format, source/member name, path or pasted inline content, optional gzip-member range, dataset/sample/donor/cohort/timepoint/compartment, and study design are retained. |
| FASTQ preprocessing | Expected-error threshold, Phred offset, all 3′ trim settings, per-dataset reservoir size, and base seed run through shared browser code. |
| References | Exact browser-composed V/D/J/C FASTA is embedded unchanged. Hand-written `references.files` replacements default to the same local metadata preparation used by the browser; `references.prepareMetadata: false` is the explicit opt-out. |
| Assignment | Workers, batch size/default, AER/RIAT-MP/standard strategy, calling profile, identity floor, strand, AIRR mode, and every Double-D threshold are executed. |
| Allele refinement | Model, pooling boundary, segments, weighting, every advanced prior/kernel/inference option, and best-versus-confidence reassignment policy are executed before collapse. |
| Collapse/denoising | Exact key and every FAD/conservative/indel parameter, unresolved policy, scope, candidate cap, and constant-call partition are executed. |
| CHMMAIRRa | Segment, BW/DB/auto choice, exact MSA, prior/mutation states/switch, DFR floor, filter threshold, unevaluated policy, and detailed Viterbi labels are executed. Selected unaligned browser references are aligned with Kalign and embedded during export. |
| Repertoire selection | The complete interactive selection schema—not only the compact pre-run fields—is retained by Results export and applied by the CLI. |
| Lineage, SHM, missing-V evidence | Lineage scope/call policy/identity/productivity/candidate cap, SHM metric, and every two-pass missing-V threshold are executed. |
| Outputs | An explicit directory/prefix, batch-streamed annotated AIRR, chunk-written retained AIRR, method summaries, and optional lazy lineage-study bundle are controlled by config; `--workers` overrides JSON. |

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
- With `references.prepareMetadata: true`, malformed identifiers, duplicate names, locus/segment mismatches, invalid IUPAC symbols, and invalid embedded `SWIGMETA` fail before SwiftIG indexing. A V/J record for which no transfer validates is retained without invented metadata and reported as a warning.
- CHMMAIRRa requires an aligned reference MSA. Browser export builds and embeds it before download when selected references are unaligned; an already aligned explicit MSA is retained byte-for-byte.
- Pipeline output requires a dedicated explicit directory; assignment-only output requires an explicit file or `-` stdout destination. The resolved pipeline config records all inherited defaults.
- A lineage range is rejected if its byte count, row count, or SHA-256 does not agree with the manifest.
- The standalone binary has no network dependency during analysis unless a future explicit input option says otherwise.

## Engineering references and relationship

- [Bun standalone executables](https://bun.sh/docs/bundler/executables) provide the runtime embedding, asset embedding, worker bundling, and per-platform compilation used for distribution. Bun does not define any biological method.
- [AIRR Rearrangement schema](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) supplies interoperable field conventions. The JSON config, lineage-sort order, byte-range manifest, range hashes, and lazy browser store are Swig-specific.
- [NCBI's IgBLAST setup documentation](https://ncbi.github.io/igblast/cook/How-to-set-up.html) defines the roles of `.ndm.imgt`, `.aux`, `-custom_internal_data`, `-auxiliary_data`, and `-d_frame_data`; the [official 1.22.0 release source/package](https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/1.22.0/) supplies their exact coordinate comments and parser behavior. Swig consumes those metadata conventions but does not implement the IgBLAST search engine.
