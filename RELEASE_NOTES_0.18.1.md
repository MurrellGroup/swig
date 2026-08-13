# Swig 0.18.1

## Result navigation

- Repertoire, Sequences, Double-D, and Post-analysis are now sticky floating tabs at the top of the results viewport.
- Every opened view remains mounted, preserving filters, expanded records, internal scroll containers, and other interactive state.
- Swig records a separate page-scroll position for each view and restores it when the user switches back.

## Double-D lineage filtering

- Lineage summaries now retain exact counts of active unique members with supported Double-D evidence and their multiplicity-weighted abundance.
- The lineage explorer can show any lineage, lineages containing at least one Double-D-positive member, lineages whose active members are all positive, or lineages with no positive member.
- The committed repertoire-level Double-D selection continues to constrain lineage assignment itself; the new control filters the already-assigned lineage table by membership.
- Restored sessions recompute these counts from the linked Double-D sidecar and saved lineage assignments.
