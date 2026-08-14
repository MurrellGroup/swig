# Swig 0.21.0

## Contextual scientific workspace

- Analysis setup now uses a narrow left rail for Data, References, Assignment,
  Pipeline, and Run. Selecting a section replaces the main workspace while the
  state of every other section remains intact.
- Each Results tab has its own contextual tools. Repertoire figures separate
  settings, gene use, distributions, V–J pairs, and sample composition;
  Sequences separates filters, the paged record table, and record detail;
  Double-D separates supported calls from the selected VDDJ alignment.
- Post-analysis now exposes Overview, Collapse, Chimera, Selection, Lineages,
  Diagnostics, Workbench, and Query in one contextual rail. Only one module is
  displayed at a time. Guided stage transitions select the next relevant
  module, while lineage assignment remains open for immediate inspection.
- Opening a sequence, double-D record, query lineage, or lineage table entry
  activates its detail workspace directly rather than updating content below
  the current scroll position.

## Visual system

- The interface uses a quieter neutral-green scientific palette, flat surfaces,
  compact type, restrained borders, and substantially reduced shadows and
  decorative treatments.
- Top-level Results tabs remain sticky. Context rails are independently sticky
  on desktop and become compact horizontally scrollable section bars on narrow
  screens.
- Landing, setup, result summaries, plots, and post-analysis panels use the same
  spacing, controls, and typography without changing analysis behavior or
  output formats.

## Compatibility

- This release changes navigation and presentation only. V(D)J assignment,
  reference composition, FASTQ filtering, post-analysis algorithms, saved
  sessions, project directories, AIRR output, and SVG/data exports retain their
  v0.20.4 behavior.
