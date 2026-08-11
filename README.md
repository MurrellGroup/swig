# Swig — GitHub Pages edition

Swig is a completely static browser interface for [SwiftIG](https://github.com/MurrellGroup/swiftig). The SwiftIG C++20 annotation core runs as WebAssembly in a bounded pool of workers; sequence data and uploaded germline sets never leave the browser.

There is no backend, database, authentication layer, or platform-specific hosting configuration in this repository.

## Deploy on GitHub Pages

1. Create an empty GitHub repository and place the contents of this directory at its root.
2. Push it to the repository's `main` branch.
3. Open **Settings → Pages** in GitHub.
4. Under **Build and deployment**, choose **GitHub Actions** as the source.
5. The included `Deploy Swig to GitHub Pages` workflow will test, build, and publish the application.

The Vite configuration derives the correct base path from `GITHUB_REPOSITORY`, so both forms work without editing source:

- `https://USERNAME.github.io/` from a `USERNAME.github.io` repository
- `https://USERNAME.github.io/REPOSITORY/` from a normal project repository

## Local development

Node.js 22.13 or newer is required.

```bash
npm ci
npm run dev
```

Production verification:

```bash
npm test
npm run build
npm run preview
```

The generated static site is written to `dist/`. It contains the JavaScript/CSS bundles, `swiftig.wasm`, the compressed reference pack, and an `index.html`; it can be served by any ordinary static file host.

## User flow

Swig is deliberately split into three views:

1. **Overview** explains the biological and computational workflow before asking for data.
2. **Analyze** accepts a file or pasted records, configures species/receptor/locus and independent V/D/J/C substitutions, and shows measured read, reference-indexing, annotation, and result-indexing progress.
3. **Results** downloads the complete AIRR TSV or explores the local result index. Runs of one to three records open directly into evidence; larger runs stay paged and filterable.

Selecting a record retrieves its full AIRR row and exposes its rearrangement map, framework/CDR regions, junction evidence, every populated AIRR field, and separate query-to-germline V, D, J, and optional C alignments. The alignment viewer switches between nucleotide and amino-acid modes. Protein alignments are translated and calculated only for the selected record.

## Large-run behavior

- Input is read incrementally from `File.stream()`. Gzip is decompressed through `DecompressionStream`; FASTA, multiline FASTQ, and AIRR rows are parsed without making a decompressed whole-file string.
- A coordinator schedules one active batch per independent WASM worker. Results may compute out of order but are committed in input order. Compute may run one batch ahead of storage, with a strict high-water mark of twice the worker count—not the number of input records.
- In Auto mode, known-large runs offer a save location before annotation and write the AIRR TSV incrementally while computation is running. Otherwise AIRR batches are gzip-compressed in IndexedDB. A service worker exposes those batches as a streaming download without assembling a full result Blob.
- The table renders 50 records at a time. Exact locus, productivity, and V/D/J allele filters use browser-local indexes; ID, CDR3 substring, identity, length, and QC filters scan the narrowed candidate set on demand and are cancellable.
- Additional filters cover V/D/J identity floors, CDR3 amino-acid length, D/CDR3 presence, receptor locus, and productivity.
- Full alignments and translated evidence are retrieved or calculated only after a row is opened.
- Every input still receives one lightweight IndexedDB record for filtering. Full AIRR rows and alignments remain in chunked storage and are read only for the selected query.
- Direct-to-disk output uses the File System Access API (Chromium-family browsers in a secure context). Other modern browsers use compressed IndexedDB plus the streaming download service worker. A conventional Blob is only the small-result fallback.

## Supported biology and formats

- FASTA, FASTQ, and AIRR Rearrangement TSV input, optionally gzip-compressed
- AIRR Rearrangement TSV output
- BCR: `IGH`, `IGK`, and `IGL`
- TCR: `TRA`, `TRB`, `TRD`, and `TRG`
- D-segment search for `IGH`, `TRB`, and `TRD`
- independent V, D, J, and optional C germline FASTA replacement
- 64 complete species/strain reference sets derived from IMGT/GENE-DB release `202632-7`

Swig runs SwiftIG, not IgBLAST. Receptor behavior is determined from the selected locus and locus-bearing gene identifiers rather than an IgBLAST `ig_seqtype` setting. Custom germline identifiers should contain one of the supported locus names above.

## Rebuild the WebAssembly core

Install [WASI SDK 25](https://github.com/WebAssembly/wasi-sdk/releases/tag/wasi-sdk-25) or a compatible newer SDK, then run:

```bash
WASI_SDK=/path/to/wasi-sdk npm run build:wasm
npm test
```

The build uses `-O3`, LTO, WebAssembly SIMD, bulk-memory operations, section garbage collection, and a pinned Binaryen `wasm-opt -O4` convergence pass. The GitHub Pages workflow downloads the pinned SDK and rebuilds this optimized core before testing and publishing. The vendored SwiftIG source and browser ABI are under `wasm/`. The upstream MIT license is retained in `LICENSE` and `wasm/LICENSE.swiftig`.

## Performance and stress verification

`npm run benchmark:50k` runs 50,000 annotations through the same one-batch-per-worker scheduling bound used by the browser. On the included nine-logical-CPU build environment, the deliberately heavy all-full-length human IGH profile (543 germline alleles, both strands, complete AIRR formatting and transfer) completed at 3,125 reads/s with eight scalar WASM workers. This is not the same workload as SwiftIG's published 5,242 reads/s native benchmark, which mixes short and full IGH/IGK/IGL reads and uses 396 NCBI alleles; do not treat the two rates as a Web/native ratio.

The normal test suite also pushes 50,000 gzip FASTA records through the incremental parser and 50,000 synthetic AIRR records through chunk compression, IndexedDB indexing, filtering, detail retrieval, and batchwise export. It asserts a 1,000-record parser high-water mark. Scaling from 50,000 to one million increases durable IndexedDB/file data and lightweight index rows, but not the number of resident query or AIRR batches.

## Update the IMGT reference pack

After downloading the official all-species ungapped IMGT/GENE-DB FASTA:

```bash
npm run build:references -- \
  /path/to/IMGTGENEDB-ReferenceSequences.fasta \
  RELEASE_ID \
  YYYY-MM-DD \
  public/references/imgt-RELEASE_ID.json.gz
```

Update the filename in `src/reference-pack.ts` when the release ID changes. The bundled file is a compact, segment-organized derivative; [IMGT attribution and terms](https://www.imgt.org/about/termsofuse.php) apply.

## Scientific scope

SwiftIG 0.2 is research software. Benchmark study-critical calls against IgBLAST or another validated workflow for the organism, assay, read length, somatic-hypermutation regime, and germline set relevant to the analysis. Public germline coverage is not literally available for every animal species; arbitrary or partial new sets can be supplied with the per-segment upload controls.
