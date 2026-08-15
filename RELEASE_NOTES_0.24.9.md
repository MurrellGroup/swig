# Swig 0.24.9

## Concatenated gzip imports

- Detects and validates concatenated gzip members during dataset upload.
- When every member begins at an independent FASTA, FASTQ, or AIRR record boundary, asks whether to import the members as separate editable samples or merge them into one sample.
- Infers useful member names from gzip metadata or SRA-style first-record identifiers when available.
- Decodes merged members sequentially through the existing bounded streaming parser, avoiding the browser `DecompressionStream` error after the first member.
- Keeps block-compressed members that continue a record merged automatically instead of incorrectly presenting blocks as biological samples.
