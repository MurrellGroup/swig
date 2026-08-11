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
2. **Analyze** accepts a file or pasted records, configures species/receptor/locus and a composable locus-by-segment reference matrix, and shows measured read, reference-indexing, annotation, and result-indexing progress.
3. **Results** downloads the complete AIRR TSV, opens a repertoire dashboard, explores the local result index, or runs opt-in post-analysis. Runs of one to three records open directly at the record detail; larger runs remain paged and filterable.

Selecting a record scrolls directly to its detail panel. A shared coordinate view layers the query, IMGT FWR/CDR intervals, and mapped V/D/J/C hits. The same panel provides the stitched V(D)J alignment, individual query-to-germline alignments, junction calls, and every populated AIRR field. Nucleotide and amino-acid modes share one biologically derived rearrangement frame; frames are not optimized independently per segment. Exact co-optimal calls and sparse near-tied candidates are shown with score, identity, and coordinates. Full alternate alignment strings are not duplicated into every AIRR row.

The repertoire dashboard accumulates summaries while AIRR batches are committed, so opening it does not rescan a million-row output. It includes customizable V/D/J/C or isotype frequency bars, CDR3-length and V-identity distributions, and a V–J pairing bubble matrix. Every rendered figure downloads as a standalone SVG. Ambiguous comma-separated calls can be counted by first call or split fractionally across co-optimal alleles.

## Post-analysis

The **Post-analysis** results tab is opt-in and keeps the full AIRR rows in chunked storage. It provides:

- **Explicit cumulative working set:** deduplication and CHMMAIRRa exclusion are separate applyable stages. Before committing either filter, the interface shows its threshold/action and retained/excluded counts. CHMMAIRRa, lineage assignment, initial sequence queries, and single-linkage expansion consume the current working set. Reset restores all assigned records; the underlying AIRR result is not deleted.
- **Exact collapse or denoising with abundance:** method A performs exact deduplication by full input sequence, VDJ-aligned sequence, locus + CDR3 nucleotide, or locus + V/J calls + CDR3. Method B reproduces FAD's corrected 6-mer distance and abundance/Poisson template rule, including the published nearest-centroid assignment, while partitioning by locus/V/J and replacing dense all-pairs scans with a complete radius index plus an exact VP-tree nearest-neighbor index. Method C is a conservative experimental alternative: it generates exact Hamming candidates with a `d+1` block index, applies a sequence-specific Poisson error model, retains isolated singletons, and never forces a read onto a distant centroid. Method D extends C to one- or two-edit Levenshtein neighborhoods: a complete length-aware `d+1` segment join generates candidates, an allocation-bounded banded dynamic program reports the exact substitution/insertion/deletion path, low-abundance indel paths use a configurable parent:child abundance ratio, and substitution-only paths retain C's Poisson rule. Method D never performs a distant-centroid assignment and reports indel and substitution merges separately. Call resolution, ambiguity handling, error rate, alpha, parent abundance, ambiguous-base policy, and method-specific radii are configurable. Trimmed VDJ sequence is streamed into a two-bit chunked arena; temporary profiles exist for one V/J partition at a time. Existing and newly collapsed `duplicate_count` values are summed and propagated into lineage weights, exports, and phylogeny bubbles.
- **CHMMAIRRa:** a browser implementation of the CHMMera structured HMM can be run after assignment on V (default) or J. IG defaults to the per-reference Baum–Welch model; TCR defaults to the discretized Bayesian model with a fixed 0.005 mutation-rate state. The default posterior threshold is 0.95, chimera prior 0.05, and minimum distance-from-reference is one. The selected run's V or J FASTA can be aligned with Kalign WASM, or an equal-length aligned FASTA MSA can be uploaded. Missing reference alleles are reported and never silently replaced. Clicking a high-posterior record reruns a detailed Viterbi trace on demand and opens an SVG-downloadable alignment with standard nucleotide coloring or an identity-only highlighter: red means query matches parent A only, blue parent B only, and gray both/neither. Tile colors never use the Viterbi state path; the inferred breakpoint remains a dashed labeled line. Columns that are gaps in the query and both displayed parents are hidden only in the visualization. Results export as a separate TSV rather than mutating the original assignment table.
- **Lineage assignment:** records are partitioned by locus, compatible V/J assignments, and exact CDR3 nucleotide length. Default gene-level ambiguity handling accepts any overlapping call, then performs single-linkage clustering at 85% CDR3 nucleotide identity. A pigeonhole `d+1` block index generates candidates, exact normalized Hamming distance verifies every retained edge, and union–find constructs components. Call resolution, ambiguity policy, productivity, identity, and the pathological-bucket comparison cap are configurable. Assignments can be exported as an AIRR table with `clone_id`; records excluded by the selected criteria retain an empty value.
- **Sequence-directed retrieval:** one or more CDR3 nucleotide, CDR3 amino-acid, or full VDJ-aligned queries can be searched by exact match, substring, bounded Hamming distance, banded edit distance, or a compact 7-mer MinHash estimate. Constraints can be manual or inferred independently for every complete nucleotide seed by running the same SwiftIG caller with the run's composed references, identity floor, strand setting, and ambiguity-aware V/J output; inferred constraints are displayed and remain overrideable. **Single-linkage expand set** starts from the current matches and follows CDR3 nucleotide edges to a fixed point or a user-visible result cap, so a candidate lineage can be recovered without assigning every lineage first.
- **Lineage figures and inspection:** abundance distribution, largest components, and V/J use are rendered as customizable SVG figures. The largest lineages are clickable; only the selected component's AIRR rows are decompressed.
- **On-demand alignment and phylogeny:** selected lineage members can use the existing AIRR coordinate strings as a fast approximate view, Kalign 3.3.1 WASM for nucleotide MSA, or a codon-aware path that aligns translated sequences and projects gaps back as triplets. The compact alignment preview has no inter-sequence vertical padding. Nucleotide and amino-acid views use standard Alivibe palettes. An N-masked reconstructed germline is included. The current alignment can be opened in Alivibe for manual/codon-preserving correction and returned directly on the same GitHub Pages origin. Cross-origin use prefers Alivibe's full-alignment download because its Copy action can contain only the current selection. Swig accepts a correction only when all original identifiers and ungapped sequences are present, fingerprints the accepted MSA, and invalidates any prior tree. Before every request, the complete current MSA is explicitly rewritten into the WASM filesystem; FastTree 2.1.11's double-precision build then consumes those exact bytes. The UI records the command, dimensions, source and fingerprint, and downloads both the named MSA and exact numeric-label input. Newick serialization writes an explicit length for every non-root edge, including `:0` for zero-length terminal branches. Germline rerooting places the root exactly at the N-masked germline (`0` root-edge length) and assigns the complete original connecting length to the ingroup side, preserving pairwise patristic distances. The complete germline-rooted FastTree resolution is shown by default; internal branches at or below `1e-8` substitutions/site can optionally be collapsed because they are at the FastTreeDbl numerical floor. The untouched raw FastTree Newick also remains inspectable and downloadable. Phylogram mode preserves retained branch lengths, and topology-only cladogram mode is explicit.
- **Coordinated tree/alignment workbench:** the exact named FastTree input is drawn to the right of its tree with leaf rows vertically locked and dotted tip-to-sequence connectors. Alignment rows have no vertical inter-sequence padding. Fractional substitution/site branch lengths are normalized to the observed root-to-tip maximum, so the deepest tip uses the complete tree span. Width and tip-spacing sliders reach zero without a hidden geometry floor; Newick order, large-clade-first ladderization, and small-clade-first ladderization are selectable. Tip bubble area reflects `duplicate_count`. Users can select FWR/CDR presets, arbitrary alignment or Kabat ranges, or only motif-matched columns; every non-contiguous displayed run receives a half-residue blank separator. Custom-region and motif drafts are applied on Enter or blur so keystrokes do not rebuild the SVG. Cells can use standard nucleotide/amino-acid palettes or motif colors; tree width, row height, and residue width are independently scalable; the complete composition can enter full screen and export as one SVG.
- **Ancestral mutation mapping and numbering:** an equal-cost nucleotide parsimony pass reconstructs internal nodes on demand. Known reconstructed-germline bases constrain the UCA; germline `N` is treated as an unknown A/C/G/T state and inferred from descendants, while gaps remain an explicit indel state. Branch labels follow the active viewer mode: nucleotide mode shows nucleotide changes, while amino-acid mode compares reconstructed parent/child codons, emits one amino-acid replacement per changed codon, and omits every synonymous nucleotide change and unresolved `X` codon. Labels can be toggled, are restricted to the displayed region, and are centered along the branches where changes occur. IG amino-acid views can switch from alignment coordinates to Kabat numbering computed locally with the MIT-licensed [Immunum](https://github.com/ENPICOM/immunum) WASM engine; multiple lineage members vote on each aligned Kabat column. A terminal partial codon or stop codon produces a visible reliability warning rather than disabling Kabat for the complete alignment. Kabat is not offered for TCR chains because the scheme is defined for IGH/IGK/IGL.

Kalign and FastTree are provided through the [bioWASM Aioli runtime](https://biowasm.com/). Their program modules are downloaded on first use. Sequence data are mounted in the worker's browser-local virtual filesystem and are not uploaded to bioWASM.

## Large-run behavior

- Input is read incrementally from `File.stream()`. Gzip is decompressed through `DecompressionStream`; FASTA, multiline FASTQ, and AIRR rows are parsed without making a decompressed whole-file string.
- Optional seeded random subsampling uses exact reservoir sampling. It scans the complete stream, retains only `k` records in memory, preserves their original order, and annotates/outputs only that uniform sample.
- A coordinator schedules one active batch per independent WASM worker. Results may compute out of order but are committed in input order. Compute may run one batch ahead of storage, with a strict high-water mark of twice the worker count—not the number of input records.
- In Auto mode, known-large runs offer a save location before annotation and write the AIRR TSV incrementally while computation is running. Otherwise AIRR batches are gzip-compressed in IndexedDB. A service worker exposes those batches as a streaming download without assembling a full result Blob.
- The table renders 50 records at a time. Exact locus, productivity, V/D/J/C allele, and isotype filters use browser-local indexes; ID, CDR3 substring, identity, length, and QC filters scan the narrowed candidate set on demand and are cancellable.
- Additional filters cover V/D/J/C identity floors, CDR3 amino-acid length, D/CDR3 presence, receptor locus, and productivity.
- Full alignment rows are retrieved only after a record is opened.
- Post-analysis scans chunked AIRR rows into a dedicated worker and retains compact sequence fingerprints, call/CDR3 keys, and typed assignment arrays—not full AIRR rows or a million-row DOM. Denoising adds a two-bit sequence arena and processes one V/J partition's neighbor profiles at a time. The 32 MB-per-million VDJ MinHash array is allocated only if full VDJ sketch search is requested.
- CHMMAIRRa uses a bounded batch per worker and stores posterior/DFR outputs in 6 bytes per input record plus a bounded high-posterior list. Its TSV writer rescans sequence IDs batchwise.
- Lineage alignment and trees are never precomputed repertoire-wide. At most the first 500 records of an opened lineage are retrieved, and the alignment limit defaults to 200.
- Every input still receives one lightweight IndexedDB record for filtering. Full AIRR rows and alignments remain in chunked storage and are read only for the selected query.
- Direct-to-disk output uses the File System Access API (Chromium-family browsers in a secure context). Other modern browsers use compressed IndexedDB plus the streaming download service worker. A conventional Blob is only the small-result fallback.

## Supported biology and formats

- FASTA, FASTQ, and AIRR Rearrangement TSV input, optionally gzip-compressed
- AIRR Rearrangement TSV output
- BCR: `IGH`, `IGK`, and `IGL`
- TCR: `TRA`, `TRB`, `TRD`, and `TRG`
- D-segment search for `IGH`, `TRB`, and `TRD`
- independent V, D, J, and optional C source selection or FASTA upload for every active locus
- constant-region calls require at least 30 aligned nucleotides; the derived `isotype` AIRR extension additionally requires at least 65% identity
- 64 species/strain reference sets derived from the IMGT-gapped all-species file in IMGT/GENE-DB release `202632-7`

Swig runs SwiftIG, not IgBLAST. Receptor behavior is determined from the selected IG/TR locus and locus-bearing gene identifiers rather than an IgBLAST `ig_seqtype` setting. The bundled pack retains light-chain `C-REGION` records and assembles IMGT's exon-level IGH/TCR constant records into one secreted/coding reference per allele; membrane-only and UTR exons are not mixed into that path. A matching C FASTA can still be supplied independently.

## Germline databases and preprocessing

The **Database** selector always defaults to the bundled IMGT/GENE-DB pack. Alternative entries are shown only for compatible species; for example, no KI entry is shown for cat. Applying a database is a preset operation, not an all-or-nothing switch: it changes only the locus/segment cells supplied by that resource and leaves every other cell unchanged. The compatible alternatives are:

| Collection | Scope | Delivery |
| --- | --- | --- |
| KIARVA | human IGH V/D/J | fetched from the KIARVA API when selected; KIARVA data are published under CC BY-NC 4.0 |
| KI human TCR database | human TRA, TRB, TRD, and TRG V/D/J as applicable | fetched from the publisher's HTTPS FASTA endpoint when selected |
| KIMDB 1.1 | rhesus and cynomolgus macaque IGH V/D/J | bundled unchanged because the publisher endpoint is HTTP-only and cannot be fetched by an HTTPS page |

The reference matrix has loci as rows and V/D/J/C as columns. Every applicable cell independently selects IMGT, a compatible published database, or an uploaded FASTA. Thus a combined BCR run can use KIARVA for IGH V/D/J while retaining IMGT for IGK/IGL, and a single cell such as IGH J can then be replaced again without changing IGH V or D. Custom and remotely loaded FASTA use the same in-browser preprocessing worker. The preprocessor:

- validates identifiers, locus/segment consistency, nucleotide symbols, duplicate names, and embedded `SWIGMETA` coordinates;
- reads exact FWR/CDR boundaries from an IMGT-gapped V record when present and anchors FWR3 on the nearest frame-consistent conserved V cysteine, accommodating species-specific IMGT gap counts;
- otherwise aligns each V allele to a fully delineated, locus-matched relative, beginning with the selected species and progressively broadening the taxonomic search, then projects the IMGT intervals through that nucleotide alignment and re-verifies the mapped V-cysteine anchor;
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

`npm run benchmark:denoise-50k` runs 50,000 synthetic assigned reads through the FAD-compatible streaming, V/J partition, exact-radius, and nearest-centroid path and reports throughput, verified candidate comparisons, and V8 heap delta.

`npm run benchmark:indel-denoise-50k` runs the same 50,000-record shape through method D with one-base insertion children. In the final verification run it retained 1,000 representatives across 100 V/J partitions, verified 1,000 candidates in 246 ms (about 203,000 input records/s), and increased measured V8 heap use by 8 MiB. This synthetic benchmark is intended to detect indexing or memory regressions, not to estimate biological correction accuracy; elapsed time and garbage-collector-dependent heap deltas vary between runs.

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

SwiftIG/Swig 0.11.0 is research software. Benchmark study-critical calls and post-analysis thresholds against an independently validated workflow for the organism, assay, read length, somatic-hypermutation regime, and germline set relevant to the analysis. The 85% CDR3 identity lineage threshold and method D's abundance ratio are starting values, not universal biological constants. Public germline coverage is not available for every animal species; arbitrary or partial new sets can be supplied with the per-segment upload controls.
