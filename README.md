# Swig — GitHub Pages edition

Swig is a completely static browser interface for [SwiftIG](https://github.com/MurrellGroup/swiftig). The SwiftIG C++20 annotation core runs as WebAssembly inside a Web Worker; sequence data and uploaded germline sets never leave the browser.

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

- SwiftIG runs in a Web Worker and annotates bounded batches.
- AIRR batches are acknowledged only after they have been committed to IndexedDB, providing backpressure between WebAssembly and browser storage.
- The table renders 50 records at a time. Exact locus, productivity, allele, frame, stop-codon, completeness, and orientation queries use browser-local indexes; ID and CDR3 substring searches scan candidates on demand and are cancellable.
- Additional filters cover V/D/J identity floors, CDR3 amino-acid length, D/CDR3 presence, receptor locus, and productivity.
- Full alignments and translated evidence are retrieved or calculated only after a row is opened.
- On browsers with the File System Access API, large AIRR downloads are streamed chunk by chunk to the chosen file. Other browsers use a conventional Blob download.

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

The vendored SwiftIG source and browser ABI are under `wasm/`. The upstream MIT license is retained in `LICENSE` and `wasm/LICENSE.swiftig`.

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
