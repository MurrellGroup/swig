# Swig 0.24.8

## UCA placement and posterior views

- Added a default fixed-alignment V/J nucleotide-mixture screen over multiple interior points on every observed-tree edge. The cheap screen only chooses starting edges; every retained point and every refinement step is evaluated under the complete recombination HMM.
- Added controls for screen mode, screen density, full-HMM edge breadth (`0` means all edges), full-HMM attachment/pendant grids, refinement rounds, and local marginalization breadth.
- Added an observed-tree placement surface showing every point used in local marginalization at its exact branch fraction, its UCA pendant length, raw full-HMM ΔLL, `exp(ΔLL)` red-to-blue color, and posterior weight.
- Added AIRR CDR1/CDR2/CDR3 bands below the UCA logo numbering.
- Added complete-canvas and exact-current-viewport SVG downloads for the HMM multi-track view.

## HMM-source tracks

- Marginalized mode now collapses every non-template and unresolved-boundary route into one NT nucleotide-mixture track at the top, followed by V, D, and J rows ordered within segment by posterior-mass location.
- Best-path/Viterbi mode retains its spatial V → NT → D-block → NT → J organization.
- Same-allele route/register alternatives remain visualization-only mixtures; inference and serialized raw tracks are unchanged.

## Repertoire allele pooling

- Moved optional pooling/reassignment before collapse, chimera filtering, and selection in both automatic and interactive workflows. Fits use the complete assigned input; downstream call-dependent partitions use the applied policy-selected calls.
- Made donor/subject the explicit safe default, combining samples for one donor while never crossing participants. Cohort and whole-study scopes are labeled cross-donor overrides.
- Kept the before/after display as hard one-record/one-allele assignment counts, with explicit every-MAP and confidence-gated (80% default) policies and vanished-allele filtering.
- Added a current-policy surviving-reference FASTA with a configurable minimum post-reassignment count (default 0).
- Constrained advanced sequence/kernel surfaces so expanded settings cannot force the page wider, and removed horizontal padding from compact sequence rows.

## Terminal alignment endpoints

- Corrected local-alignment endpoint materialization so a terminal run of V 5′ or J 3′ somatic substitutions is not misreported as missing read/reference sequence.
- Candidate discovery and dynamic programming are unchanged; the selected-hit endpoint extension is linear only in omitted terminal bases.
- Tested all 34 supplied monoclonal sequences against the KIARVA human IGH database: V query/reference starts are now 1 for every record, while complete J 3′ endpoints remain complete. Added a synthetic terminal-SHM WebAssembly regression test.

## Verification

- Rebuilt the optimized WebAssembly core and verified TypeScript, focused UCA/allele tests, the full automated suite, and the production GitHub Pages build.
