# V(D)J assignment and AIRR annotation

## Core execution

SwiftIG's C++20 core runs as WebAssembly in a bounded worker pool. Each worker receives the same composed germline references and processes independent record batches. Swig preserves input order when committing AIRR batches. Search strand, minimum accepted alignment identity, assignment strategy, and calling profile are explicit run settings.

For each oriented query, SwiftIG uses exact nucleotide seed indexes to retrieve candidates, exact affine local alignment for selected candidates, then reconstructs a consistent V-(D)-J rearrangement and AIRR coordinate/alignment fields. D and J evidence is evaluated in the V/J-bounded region. Co-optimal calls are emitted as comma-separated names; sparse near-tied candidates retain score/identity/coordinate evidence without duplicating full alignment strings into every AIRR row.

The minimum identity control is an acceptance floor for candidate alignments, not a global read-trimming heuristic. The `sequence` field remains the user's input sequence (or the chosen reverse complement); V/J-aligned and junction fields are derived outputs.

## V assignment strategies

The strategy changes V candidate retrieval/refinement only. D and J use the selected calling profile in every strategy.

Swig Web starts all three front-page workflows with RIAT-MP. The command-line config and `swig-cli --vdj` retain AER as their default for backward compatibility. Either surface can select any strategy explicitly, and a browser-exported config records the web selection.

### AER—Adaptive Exact Refinement

AER begins from the complete V-allele seed index. It exactly aligns the leading candidates and increases exact affine-alignment depth only when the leading 9-mer vote ranking is ambiguous: within 5% relative vote count or 8 weighted votes, up to 16 candidates. It never substitutes a propagated approximate score for the final retained exact alignment.

### RIAT-MP (Swig Web default)

RIAT-MP groups close V alleles into root-indexed trees. It aligns up to three representative roots and propagates sparse descendant differences without performing a full descendant alignment. When the provisional winner contains an indel, it tests at most two root traceback geometries within four raw-score units, with a 1,024-state traceback cap. This is a Swig/SwiftIG algorithm, not a literature package.

### Standard SwiftIG

The standard strategy exactly aligns the three leading strong-seed V candidates. Weak-seed and seedless cases retain a larger safety pool. It is fixed-depth rather than ambiguity-adaptive.

## Calling profiles

The profile is orthogonal to V strategy.

- **Truth-optimized (default):** the SwiftIG settings selected on the supplied simulated human-IGH truth data.
- **IgBLAST-agreement:** D scoring `+2/−4/−11/−1`, minimum 5-nt exact run and 3 candidates; J scoring `+2/−4/−13/−1` and 2 candidates. These values were selected for agreement with supplied IgBLAST calls.
- **IgBLAST-balanced:** the agreement settings plus removal of a D call only when its strongest support is exactly five consecutive matches and `j_sequence_start − v_sequence_end ≤ 11`. It maximized tuning-set IgBLAST agreement subject to mean first-call and ambiguity-aware V/D/J truth accuracy exceeding IgBLAST on that simulation.

Neither agreement profile is an IgBLAST implementation. IgBLAST uses its own BLAST-based germline search and annotation logic; see Ye et al., [IgBLAST](https://doi.org/10.1093/nar/gkt382). The profile names describe calibration objectives on one supplied simulation, not held-out validation or universal superiority.

## Segment support values

For every selected hit whose scoring tuple has been calibrated, SwiftIG writes the AIRR fields `v_support`, `d_support`, `j_support`, and `c_support`. These are BLAST-form expectation values: the expected number of chance local alignments with score at least as large as the reported raw segment score. Smaller is stronger. They are neither posterior probabilities nor repertoire-level evidence that an allele is present.

For raw score \(S\), SwiftIG evaluates

\[
E = K m' n' e^{-\lambda S}, \qquad
m' = m-\ell, \quad n'=n-N\ell,
\]

where \(m\) is the searched query length, \(n\) and \(N\) are the total nucleotide length and record count of the supplied segment database, and \(\ell\) is the standard finite-length adjustment fitted with \(\alpha\) and \(\beta\). D uses the actual V/J-bounded junction span searched; V, J, and C use the complete oriented query. Candidate indexing does not shrink the reference search space: support remains relative to every supplied allele in that segment database.

The \(\lambda,K,\alpha,\beta\) constants were fitted offline by deterministic uniform-DNA simulation of SwiftIG's exact affine local-alignment recurrence for all shipped truth-optimized and IgBLAST-agreement scoring tuples. Runtime does no simulation, realignment, or database scan. The reference totals are cached once per worker; each retained hit adds only constant-time arithmetic. An uncalibrated custom/benchmark-only scoring tuple produces an empty support field rather than borrowing invalid constants.

The mathematical form follows Karlin and Altschul's [local-alignment statistics](https://doi.org/10.1073/pnas.87.6.2264) and the gapped-search treatment used by BLAST ([Altschul et al. 1997](https://doi.org/10.1093/nar/25.17.3389)). The numbers are calibrated to SwiftIG's scores and searched spans. They are therefore interpretable E-values but are **not expected to be numerically identical to IgBLAST E-values**, because IgBLAST uses different scoring/search pipelines, databases, query ranges, and hit construction. Calibration and held-out null checks are recorded in [`BENCHMARK_AIRR_SUPPORT_0.34.0.md`](../../BENCHMARK_AIRR_SUPPORT_0.34.0.md).

## Output and record detail

Swig writes an AIRR rearrangement table plus `swig_*` study/provenance fields. The interactive record view is a rendering of those committed fields: query and germline coordinates, FWR/CDR intervals, stitched V(D)J alignment, individual segment alignments, segment E-values, junction decomposition, and retained candidate evidence. It does not rerun assignment when opened.

The schema follows the AIRR Rearrangement representation so downstream tools can consume the table: [AIRR Community standardized representations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6173121/). Swig adds namespaced columns where the AIRR core schema has no equivalent. AIRR compatibility does not imply that SwiftIG is another AIRR-producing caller such as IgBLAST or IMGT/V-QUEST.

## Performance and cancellation

Reference indexes are built once per worker. Records are processed in bounded batches, and results are appended to browser storage or a user-selected writable stream. Cancellation aborts input reading, terminates the worker pool, aborts unfinalized output, and clears the temporary result store; no partial run is opened as a completed analysis.

## Limitations

- AER/RIAT-MP/profile calibration used simulated human IGH and is not held-out validation.
- Short D segments and highly mutated reads can remain intrinsically ambiguous.
- Co-optimal labels describe equal scores under this caller, not equal biological posterior probabilities.
- C calls depend on whether the input actually covers constant sequence.
- Important study-specific calls should be benchmarked against an independent workflow and suitable truth/control material.
