# Swig 0.33.0 release notes

## Streaming CLI output

- Pipeline runs now require an explicit output directory in `output.directory` or `--out`.
- The annotated AIRR table is written as worker batches complete in input order. The retained table is serialized in bounded chunks after enabled repertoire-wide decisions, avoiding repertoire-sized output strings.
- `--workers N` is applied after JSON loading and therefore overrides `annotation.workers`; `--workers 0` selects the bounded host default.

## Low-overhead `swig-cli --vdj`

- Added a separate assignment-only path with IgBLAST-style `-query`, germline, annotation, thread, strand, format, and output option names.
- `-out` is mandatory and AIRR outfmt 19 is streamed directly from bounded SwiftIG batches. No study metadata, AIRR row objects, downstream repertoire analysis, summaries, or lineage state are constructed.
- The default path canonicalizes/ungaps source germline FASTA but deliberately leaves FWR/CDR/junction fields unannotated.
- `--swigannots` enables the browser-equivalent germline metadata inference route.
- IgBLAST `-custom_internal_data`, `-auxiliary_data`, and `-d_frame_data` are parsed with their real coordinate conventions. Exact identifiers and chain/locus consistency are checked, and the optional fifth J auxiliary field trims trailing non-coding bases from FWR4 within each bounded output batch.
- Supplying only J auxiliary data uses the embedded fixed species pack for the V-domain role normally served by IgBLAST's installed internal-data directory.
- Germline arguments accept source FASTA, not binary `makeblastdb` prefixes. Unsupported BLAST search/rendering flags fail explicitly; this remains SwiftIG assignment and AIRR output, not an IgBLAST implementation.

## Verification

- Added unit tests for IMGT-gap removal, exact V internal-data coordinates, J auxiliary coordinates/FWR4 offsets, and D frame files.
- Added end-to-end CLI tests for explicit-output enforcement, command-line worker precedence, annotation-free VDJ output, IgBLAST-data annotations, Swig-inferred annotations, mandatory VDJ output, and absence of downstream artifacts.
