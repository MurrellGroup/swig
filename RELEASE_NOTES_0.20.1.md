# Swig 0.20.1

## Alivibe nucleotide round trip

- Replaced the external-editor timing/DOM shim with a pinned, bundled Alivibe
  copy and a versioned same-origin bridge.
- Direct return now forces Alivibe into NT mode and snapshots every ordered row
  from the same `state.viewSequences` state used by Alivibe's NT canvas and full
  NT FASTA export. It cannot return the AA display, a canvas selection, or stale
  clipboard text.
- Alivibe controls are attached before its optional Aioli/WASM initialization,
  eliminating the interval in which the editor functions existed but the NT
  and frame controls were inert.
- Every open editor is bound to the originating lineage set, input alignment
  fingerprint, and bundled Alivibe revision. Swig refuses stale or mismatched
  returns rather than applying them to another lineage.
- The AA reading-frame offset remains separate metadata. Nucleotide characters,
  gap positions, row order, and retained identifiers must survive the direct
  return byte-for-byte after canonical FASTA serialization; any return that
  would require silent normalization is rejected.
- Downloaded nucleotide FASTA remains an explicit fallback. Direct return no
  longer reads or writes the system clipboard.

Regression tests cover exact NT rows and gaps, AA/stale/inconsistent snapshot
rejection, full-row export semantics, pinned source revision, and control
initialization order.
