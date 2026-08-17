# Swig 0.35.3 release notes

## AA back-projection correction

- Fixed an empty-padding error that inserted a synthetic leading `---` into every frame-0 AA-guided alignment row.
- The same correction removes a synthetic terminal `---` when no residual partial codon exists.
- Genuine 1–2 nucleotide pre-frame and terminal fragments remain represented as padded partial-codon blocks, while every internal amino-acid alignment gap still projects to exactly `---`.
- A deterministic regression now requires `ATGGGGCCC` aligned as `M-GP` to project to exactly `ATG---GGGCCC`, with neither a leading nor trailing gap codon.

## Alivibe edit recovery

- The bridge ABI remains strict for running MSA: a v3 editor cannot invoke the v4 string-valued nucleotide/amino-acid scoring protocol.
- The unchanged nucleotide snapshot schema is now recovery-compatible from bridge v3, allowing an already-open older editor to return existing manual nucleotide edits safely.
- Newly opened editor URLs include bridge-version and pinned-source cache keys, preventing stale editor HTML from being paired with a newer Swig app.
