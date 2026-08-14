# Swig 0.22.0

This release reorganizes the interface without removing analysis capabilities.

- V(D)J assignment is one continuous action page: data, study metadata, biological search space, assignment/input options, optional pipeline stages, and launch are no longer split across setup tabs.
- Strong-default parameters use progressive disclosure. Biologically consequential choices—species, receptor/locus, reference database, collapse method and boundary, unresolved-record policy, lineage boundary/threshold, and double-D mode—remain directly editable.
- KI collections are the context-aware defaults where available: KIARVA for human IGH, the KI human TCR database for human TR loci, and KIMDB for macaque IGH. IMGT remains the fallback and supplies loci/segments absent from the selected KI collection.
- Species, dataset/sample metadata, and V/D/J/C calls use searchable inputs rather than long enumerated menus.
- V/D/J/C result filters accept a gene or allele and independently control whether ambiguous multi-hit assignments containing that target are included. Exact comma-separated multi-hit call-set matching remains available.
- Gene/allele tokens are stored in browser-local multi-entry indexes so targeted filters remain indexed at repertoire scale.
- Double-D exploration moved under Sequences. The redundant top-level Double-D results tab was removed; D1/D2 evidence filters, tables, complete VDDJ alignments, and full-record inspection remain available.
- Repertoire population/figure settings are editable on every figure page instead of living in a separate settings panel.
- Collapse/denoising, CHMMAIRRa, repertoire selection, lineage assignment, SHM, and missing-allele controls were regrouped into primary biological choices and expandable calibrated parameters.
- Every ordinary labeled control receives mouseover help; setting-style buttons retain inline explanations and hover text.
- Large sample lists in the lineage explorer are searched and bounded rather than rendered as an unbounded checkbox wall.

All sequence analysis remains browser-local. Selecting a remotely hosted reference collection requests that reference from its named provider.
