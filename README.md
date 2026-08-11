# Swig — GitHub Pages edition

Swig is a static browser interface for [SwiftIG](https://github.com/MurrellGroup/swiftig). The SwiftIG C++20 annotation core runs as WebAssembly in a bounded worker pool. Query sequences, uploaded germline sets, and AIRR results are processed in the browser and are not transmitted by Swig. Selecting an online KI collection downloads only the chosen germline FASTA from its provider.

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
3. **Results** downloads the complete AIRR TSV, opens a repertoire dashboard, or explores the local result index. Runs of one to three records open directly at the record detail; larger runs remain paged and filterable.

Selecting a record scrolls directly to its detail panel. A shared coordinate view layers the query, IMGT FWR/CDR intervals, and mapped V/D/J/C hits. The same panel provides the stitched V(D)J alignment, individual query-to-germline alignments, junction calls, and every populated AIRR field. Nucleotide and amino-acid modes share one biologically derived rearrangement frame; frames are not optimized independently per segment. Exact co-optimal calls and sparse near-tied candidates are shown with score, identity, and coordinates. Full alternate alignment strings are not duplicated into every AIRR row.

The repertoire dashboard accumulates summaries while AIRR batches are committed, so opening it does not rescan a million-row output. It includes customizable V/D/J/C or isotype frequency bars, CDR3-length and V-identity distributions, and a V–J pairing bubble matrix. Every rendered figure downloads as a standalone SVG. Ambiguous comma-separated calls can be counted by first call or split fractionally across co-optimal alleles.

## Large-run behavior

- Input is read incrementally from `File.stream()`. Gzip is decompressed through `DecompressionStream`; FASTA, multiline FASTQ, and AIRR rows are parsed without making a decompressed whole-file string.
- Optional seeded random subsampling uses exact reservoir sampling. It scans the complete stream, retains only `k` records in memory, preserves their original order, and annotates/outputs only that uniform sample.
- A coordinator schedules one active batch per independent WASM worker. Results may compute out of order but are committed in input order. Compute may run one batch ahead of storage, with a strict high-water mark of twice the worker count—not the number of input records.
- In Auto mode, known-large runs offer a save location before annotation and write the AIRR TSV incrementally while computation is running. Otherwise AIRR batches are gzip-compressed in IndexedDB. A service worker exposes those batches as a streaming download without assembling a full result Blob.
- The table renders 50 records at a time. Exact locus, productivity, V/D/J/C allele, and isotype filters use browser-local indexes; ID, CDR3 substring, identity, length, and QC filters scan the narrowed candidate set on demand and are cancellable.
- Additional filters cover V/D/J/C identity floors, CDR3 amino-acid length, D/CDR3 presence, receptor locus, and productivity.
- Full alignment rows are retrieved only after a record is opened.
- Every input still receives one lightweight IndexedDB record for filtering. Full AIRR rows and alignments remain in chunked storage and are read only for the selected query.
- Direct-to-disk output uses the File System Access API (Chromium-family browsers in a secure context). Other modern browsers use compressed IndexedDB plus the streaming download service worker. A conventional Blob is only the small-result fallback.

## Supported biology and formats

- FASTA, FASTQ, and AIRR Rearrangement TSV input, optionally gzip-compressed
- AIRR Rearrangement TSV output
- BCR: `IGH`, `IGK`, and `IGL`
- TCR: `TRA`, `TRB`, `TRD`, and `TRG`
- D-segment search for `IGH`, `TRB`, and `TRD`
- independent V, D, J, and optional C germline FASTA replacement
- constant-region calls require at least 30 aligned nucleotides; the derived `isotype` AIRR extension additionally requires at least 65% identity
- 64 species/strain reference sets derived from the IMGT-gapped all-species file in IMGT/GENE-DB release `202632-7`

Swig runs SwiftIG, not IgBLAST. Receptor behavior is determined from the selected IG/TR locus and locus-bearing gene identifiers rather than an IgBLAST `ig_seqtype` setting. The bundled pack retains C-REGION records where the IMGT source provides them; a matching C FASTA can be supplied independently.

## Germline collections and preprocessing

The default pack is derived from IMGT/GENE-DB. The chain-specific selector also exposes the following KI resources:

| Collection | Scope | Delivery |
| --- | --- | --- |
| KIARVA | human IGH V/D/J | fetched from the KIARVA API when selected; KIARVA data are published under CC BY-NC 4.0 |
| KI human TCR database | human TRA, TRB, TRD, and TRG V/D/J as applicable | fetched from the publisher's HTTPS FASTA endpoint when selected |
| KIMDB 1.1 | rhesus and cynomolgus macaque IGH V/D/J | bundled unchanged because the publisher endpoint is HTTP-only and cannot be fetched by an HTTPS page |

Each collection can replace all of its segments or only V, D, or J. Custom FASTA uses the same in-browser preprocessing worker. The preprocessor:

- validates identifiers, locus/segment consistency, nucleotide symbols, duplicate names, and embedded `SWIGMETA` coordinates;
- reads exact FWR/CDR boundaries from an IMGT-gapped V record when present;
- otherwise aligns each V allele to a fully delineated, locus-matched relative, beginning with the selected species and progressively broadening the taxonomic search, then projects the IMGT intervals through that nucleotide alignment;
- resolves J frame and CDR3-stop metadata by homologous anchor transfer when possible and re-verifies the mapped conserved F/W–G motif in the submitted sequence;
- ignores functional, ORF, and pseudogene labels when selecting coordinate donors;
- leaves a record without region or junction metadata if it does not span the necessary interval or the mapped anchor cannot be verified. It does not substitute a motif-derived or fixed-coordinate default silently.

No BLAST database files or IgBLAST `edit_imgt_file.pl`/`makeblastdb` step is required because SwiftIG indexes normalized FASTA directly. The IMGT gap layout is used only to derive coordinate metadata before the sequence is ungapped for SwiftIG's index.

## Rebuild the WebAssembly core

Install [WASI SDK 25](https://github.com/WebAssembly/wasi-sdk/releases/tag/wasi-sdk-25) or a compatible newer SDK, then run:

```bash
WASI_SDK=/path/to/wasi-sdk npm run build:wasm
npm test
```

The build uses `-O3`, LTO, WebAssembly SIMD, bulk-memory operations, section garbage collection, and a pinned Binaryen `wasm-opt -O4` convergence pass. The GitHub Pages workflow downloads the pinned SDK and rebuilds this optimized core before testing and publishing. The vendored SwiftIG source and browser ABI are under `wasm/`. The upstream MIT license is retained in `LICENSE` and `wasm/LICENSE.swiftig`.

## Performance and stress verification

`npm run benchmark:50k` runs 50,000 annotations through the same one-batch-per-worker scheduling bound used by the browser. The deliberately heavy profile uses all full-length human IGH references, both strands, complete AIRR formatting and transfer. This is not the same workload as SwiftIG's published 5,242 reads/s native benchmark, which mixes short and full IGH/IGK/IGL reads and uses 396 NCBI alleles; do not treat the two rates as a Web/native ratio.

The normal test suite also pushes 50,000 gzip FASTA records through the incremental parser and 50,000 synthetic AIRR records through chunk compression, IndexedDB indexing, filtering, detail retrieval, and batchwise export. It asserts a 1,000-record parser high-water mark. Scaling from 50,000 to one million increases durable IndexedDB/file data and lightweight index rows, but not the number of resident query or AIRR batches.

## Update the IMGT reference pack

After downloading the official IMGT-gapped all-species IMGT/GENE-DB FASTA:

```bash
npm run build:references -- \
  /path/to/IMGTGENEDB-ReferenceSequences.fasta \
  RELEASE_ID \
  YYYY-MM-DD \
  public/references/imgt-RELEASE_ID.json.gz
```

Update the filename in `src/reference-pack.ts` when the release ID changes. The bundled file is a compact, segment-organized derivative; [IMGT attribution and terms](https://www.imgt.org/about/termsofuse.php) apply.

## Scientific scope

SwiftIG/Swig 0.5 is research software. Benchmark study-critical calls against IgBLAST or another validated workflow for the organism, assay, read length, somatic-hypermutation regime, and germline set relevant to the analysis. Public germline coverage is not available for every animal species; arbitrary or partial new sets can be supplied with the per-segment upload controls.
