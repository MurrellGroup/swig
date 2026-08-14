# Swig 0.22.1

This corrective release repairs the interaction and layout regressions in 0.22.0 without removing analysis options.

- Species is again a finite native dropdown. Macaque species/strain entries remain directly selectable without first typing into a search field.
- KI remains the context-aware default where available. Selecting a species, receptor, database preset, or per-segment source updates the displayed choice immediately; reference preprocessing continues asynchronously and stale requests cannot overwrite a newer selection.
- Advanced controls are compact buttons on their parent panel. Expanded controls appear as one attached surface beneath that button instead of as an empty standalone card.
- Sequence filters and the records table now share one workspace. Sample, donor, cohort, timepoint, compartment, locus, isotype, and call filters use observed-value selectors; high-cardinality V/D/J/C selectors are searchable and accept typed genes or alleles.
- Ambiguous multi-call matching is explicit for V/D/J/C filters. The selector presents individual gene/allele tokens rather than every comma-separated multi-hit combination.
- Double-D evidence is demoted to an optional rare-event filter and a secondary sequence explorer rather than occupying a primary results tab.
- The per-record viewer restores its dark scientific surface, bounded responsive grids, and horizontal scrolling for irreducibly wide alignments and evidence tables.
- Collapse/denoising controls use responsive columns, and the constant-gene policy no longer escapes its panel.

All analysis remains local to the browser. Selecting a remotely hosted reference collection requests that reference from its named provider.
