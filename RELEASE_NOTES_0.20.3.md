# Swig 0.20.3

## AIRR-anchored lineage frame

- Newly generated reference-anchored lineage alignments no longer choose an AA
  frame by scoring alignment gaps. The selected germline-template member
  supplies `sequence_frame` and `v_sequence_start`; `v_germline_start` padding
  is then projected into the displayed alignment columns.
- Nucleotide Kalign output derives its displayed phase from that same named
  anchor after alignment. Codon-projected Kalign output starts at column one by
  construction. Gap-based inference is retained only for legacy inputs without
  frame metadata.
- Alivibe's frame selector remains enabled in both NT and AA modes. An explicit
  change is returned to Swig and retained with manually edited alignments,
  saved sessions, tree state, and SVG mutation coordinates.
- Corrected FASTA imported without separate frame metadata keeps the current
  explicit phase instead of silently redefining it from newly moved gaps.

## Sequence-level CDR3 context

- The paged AIRR table shows CDR3 AA with nucleotide fallback and, when both
  exist, retains the nucleotide sequence as a second compact line.
- The Double-D explorer restores CDR3 to its clickable rows, adds CDR3 search,
  and shows CDR3 beside the selected VDDJ alignment. Values are hydrated from
  the compact main index rather than duplicated into the sparse sidecar.
- Clickable CHMMAIRRa candidates and selected-lineage member controls now show
  their sequence identifier and CDR3 before opening the detailed alignment.

Regression tests cover AIRR phase propagation through reference padding,
Alivibe frame-control availability in NT mode, and CDR3 presence in all
sequence-level click targets.
