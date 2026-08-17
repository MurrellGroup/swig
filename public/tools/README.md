# Bundled Alivibe

`alivibe.html`, `nw.js`, `phylotools.js`, and `frameclean.js` are pinned from
MurrellGroup/WebWidgets revision
`cbcd02719dd0a5f1f05d3127666f00e8579f2423`.

Swig adds a small, versioned same-origin bridge to `alivibe.html`. The bridge
loads nucleotide FASTA and snapshots every ordered row from Alivibe's
`state.viewSequences` after calling `setMode('NT')`; this is the same state used
by the nucleotide canvas and the full nucleotide FASTA export. It does not read
the system clipboard or selected canvas cells.

For Swig-hosted editors, the bridge also injects a cancellable worker running
`public/alivibe-msa.wasm`. That module is a behavior-preserving C++/WASM port of
the pinned `refinedMSA` route in `nw.js`; standalone Alivibe keeps the original
JavaScript function as its fallback.
