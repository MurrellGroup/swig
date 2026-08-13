# Swig 0.17.0

- Uses the longitudinal / compartmental study preset by default. Exact collapse remains sample-bounded; lineage assignment is donor-bounded.
- Adds optional compartment/tissue metadata to AIRR rows, result filtering, pipeline selection, saved sessions, and post-analysis strata.
- Keeps the lineage-assignment module open after assignment so its table remains the immediate next analysis surface.
- Adds lineage-table filters for sample breadth, selected-sample presence, locus, V/J text, abundance, unique members, CDR3 length, and SHM thresholds.
- Adds sortable lineage SHM mean, maximum, and upper-95%-quantile columns.
- Adds an SVG-exportable lineage-by-sample SHM heatmap with editable sample order and input/timepoint/compartment/sample-ID ordering presets.
- Changes missing-allele result wording from “validated” to candidate/evidence terminology.
- Revises interface copy toward descriptive scientific language.
