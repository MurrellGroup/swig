# V(D)J assignment and AIRR annotation

## Core execution

SwiftIG's C++20 core runs as WebAssembly in a bounded worker pool. Each worker receives the same composed germline references and processes independent record batches. Swig preserves input order when committing AIRR batches. Search strand, minimum accepted alignment identity, assignment strategy, and calling profile are explicit run settings.

For each oriented query, SwiftIG uses exact nucleotide seed indexes to retrieve candidates, exact affine local alignment for selected candidates, then reconstructs a consistent V-(D)-J rearrangement and AIRR coordinate/alignment fields. D and J evidence is evaluated in the V/J-bounded region. Co-optimal calls are emitted as comma-separated names; sparse near-tied candidates retain score/identity/coordinate evidence without duplicating full alignment strings into every AIRR row.

The minimum identity control is an acceptance floor for candidate alignments, not a global read-trimming heuristic. The `sequence` field remains the user's input sequence (or the chosen reverse complement); V/J-aligned and junction fields are derived outputs.

## V assignment strategies

The strategy changes V candidate retrieval/refinement only. D and J use the selected calling profile in every strategy.

### AER—Adaptive Exact Refinement (default)

AER begins from the complete V-allele seed index. It exactly aligns the leading candidates and increases exact affine-alignment depth only when the leading 9-mer vote ranking is ambiguous: within 5% relative vote count or 8 weighted votes, up to 16 candidates. It never substitutes a propagated approximate score for the final retained exact alignment.

### RIAT-MP

RIAT-MP groups close V alleles into root-indexed trees. It aligns up to three representative roots and propagates sparse descendant differences without performing a full descendant alignment. When the provisional winner contains an indel, it tests at most two root traceback geometries within four raw-score units, with a 1,024-state traceback cap. This is a Swig/SwiftIG algorithm, not a literature package.

### Standard SwiftIG

The standard strategy exactly aligns the three leading strong-seed V candidates. Weak-seed and seedless cases retain a larger safety pool. It is fixed-depth rather than ambiguity-adaptive.

## Calling profiles

The profile is orthogonal to V strategy.

- **Truth-optimized (default):** the SwiftIG settings selected on the supplied simulated human-IGH truth data.
- **IgBLAST-agreement:** D scoring `+2/−4/−11/−1`, minimum 5-nt exact run and 3 candidates; J scoring `+2/−4/−13/−1` and 2 candidates. These values were selected for agreement with supplied IgBLAST calls.
- **IgBLAST-balanced:** the agreement settings plus removal of a D call only when its strongest support is exactly five consecutive matches and `j_sequence_start − v_sequence_end ≤ 11`. It maximized tuning-set IgBLAST agreement subject to mean first-call and ambiguity-aware V/D/J truth accuracy exceeding IgBLAST on that simulation.

Neither agreement profile is an IgBLAST implementation. IgBLAST uses its own BLAST-based germline search and annotation logic; see Ye et al., [IgBLAST](https://doi.org/10.1093/nar/gkt382). The profile names describe calibration objectives on one supplied simulation, not held-out validation or universal superiority.

## Output and record detail

Swig writes an AIRR rearrangement table plus `swig_*` study/provenance fields. The interactive record view is a rendering of those committed fields: query and germline coordinates, FWR/CDR intervals, stitched V(D)J alignment, individual segment alignments, junction decomposition, and retained candidate evidence. It does not rerun assignment when opened.

The schema follows the AIRR Rearrangement representation so downstream tools can consume the table: [AIRR Community standardized representations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6173121/). Swig adds namespaced columns where the AIRR core schema has no equivalent. AIRR compatibility does not imply that SwiftIG is another AIRR-producing caller such as IgBLAST or IMGT/V-QUEST.

## Performance and cancellation

Reference indexes are built once per worker. Records are processed in bounded batches, and results are appended to browser storage or a user-selected writable stream. Cancellation aborts input reading, terminates the worker pool, aborts unfinalized output, and clears the temporary result store; no partial run is opened as a completed analysis.

## Limitations

- AER/RIAT-MP/profile calibration used simulated human IGH and is not held-out validation.
- Short D segments and highly mutated reads can remain intrinsically ambiguous.
- Co-optimal labels describe equal scores under this caller, not equal biological posterior probabilities.
- C calls depend on whether the input actually covers constant sequence.
- Important study-specific calls should be benchmarked against an independent workflow and suitable truth/control material.
