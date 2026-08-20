# V(D)J assignment and AIRR annotation

## Core execution

SwiftIG's C++20 core runs as WebAssembly in a bounded worker pool. Each worker receives the same composed germline references and processes independent record batches. Swig preserves input order when committing AIRR batches. Search strand, minimum accepted alignment identity, assignment strategy, and calling profile are explicit run settings.

For each oriented query, SwiftIG uses exact nucleotide seed indexes to retrieve candidates, exact affine local alignment for selected candidates, then reconstructs a consistent V-(D)-J rearrangement and AIRR coordinate/alignment fields. The ordinary RIAT-MP, AER, and standard paths evaluate D in the V/J-bounded region. Experimental AER-R begins there too, but can expand an uncertain junction to mapped conserved V/J anchors before comparing complete non-overlapping partitions. Co-optimal calls are emitted as comma-separated names; sparse near-tied candidates retain score/identity/coordinate evidence without duplicating full alignment strings into every AIRR row.

V and J are first optimized independently, but their junction-facing endpoints are not treated as irrevocable. If the ordinary D call has less than 10 consecutive exact bases and a V 3′ or J 5′ affine gap lies within 24 query bases of the junction-facing edge, SwiftIG checks the reclaimed boundary span for an exact 10-nt D seed. Only then does it perform a second D refinement. It clips the competing V/J alignment at each candidate D boundary, recomputes the retained affine scores, and accepts the split only when the non-overlapping V+D+J score is strictly greater than the original V+(weak D)+J score. This custom rescue prevents a gappy terminal V/J optimum from hiding a long exact D tract. Strong ordinary D calls and endpoints without a junction-facing gap never execute the additional search.

The minimum identity control is an acceptance floor for candidate alignments, not a global read-trimming heuristic. The `sequence` field remains the user's input sequence (or the chosen reverse complement); V/J-aligned and junction fields are derived outputs.

## Assignment strategies

RIAT-MP, ordinary AER, and standard SwiftIG change V candidate retrieval/refinement only. AER-R is a separately selectable experimental AER derivative that also changes how uncertain V/D/J boundaries and candidates are adjudicated. Every strategy uses the same selected calling-profile scores and minimum evidence settings.

RIAT-MP is the default for all three front-page workflows, hand-written/normalized CLI configs, and the direct `swig-cli --vdj` route. Either surface can select any strategy explicitly, and a browser-exported config records the web selection.

### AER—Adaptive Exact Refinement

AER begins from the complete V-allele seed index. It exactly aligns the leading candidates and increases exact affine-alignment depth only when the leading 9-mer vote ranking is ambiguous: within 5% relative vote count or 8 weighted votes, up to 16 candidates. It never substitutes a propagated approximate score for the final retained exact alignment.

The production AER kernel evaluates independent affine recurrences in four exact SIMD128 lanes. A score-only pass over a wide safety pool omits trace storage but uses the same recurrence; full DP/traceback continues through every raw-score group capable of entering the retained top-N set. The scalar affine implementation remains the executable reference path used by equivalence tests.

### AER-R—experimental robust joint caller

AER-R leaves ordinary AER intact and is selected explicitly as `aer_robust` in portable config or `--assigner aer_robust` in the CLI. It reuses AER's exact seed candidates, affine recurrence, calling-profile scores, AIRR writer, and clean V/J-bounded D search. Its additional decisions are:

- compare explicit non-overlapping V–D–J partitions, charging the complete affine-score loss when V or J must be clipped;
- retain the nominal V/J-bounded D location and, only when it is absent, boundary-truncated, or weak beside a junction-facing gap, add a conserved-anchor-bounded D search as a distinct location hypothesis;
- map the final V cysteine and J F/W–G anchors from reference through the alignment, with a database-length-derived fallback only when metadata or alignment coverage cannot map an anchor;
- retain exact score ties at a configured candidate cutoff and add exact retained-D-substring allele labels that a top-N list alone could hide;
- if a strongly V-supported orientation has no valid V/J partition, exactly align the complete small J database; the opposite/weak orientation does not pay for that rescue;
- optimize a non-overlapping J/C split when independent local hits overlap, retaining it only when the scored constant contribution exceeds the lost J score; and
- when the separate double-D screen is enabled, permit a shorter seed rescue that must ungapped-extend back to at least the configured ordinary seed length before entering the unchanged pair-gain and pseudo-tandem tests.

Anchor projection and the maximum D-reference length are cached or evaluated lazily. Clean geometry therefore performs the same number of D searches as ordinary AER. The development simulator contains skewed allele use, V/D/J deletion, P/N addition, occasional tandem D, hotspot-weighted SHM, short indels, base-call errors, ambiguous bases, partial reads, flanks, reverse complements, and exact base provenance. It and the compile-time all-reference oracle live under `tests/`; neither is linked into the runtime WebAssembly. Results and limitations are in [`BENCHMARK_AER_ROBUST_0.37.6.md`](../../BENCHMARK_AER_ROBUST_0.37.6.md).

### RIAT-MP (Swig default)

RIAT-MP groups close V alleles into root-indexed trees. It aligns up to three representative roots and propagates sparse descendant differences without performing a full descendant alignment. When the provisional winner contains an indel, it tests at most two root traceback geometries within four raw-score units, with a 1,024-state traceback cap. This is a Swig/SwiftIG algorithm, not a literature package.

Production root alignments use the same exact four-lane SIMD kernel. Tree child links are precomputed once, and alternative partial tracebacks use parent-linked arena nodes; neither changes the path priority, traceback cap, or sparse descendant score. The complete 0.37.0 profile and reference-equivalence protocol are in [`BENCHMARK_SWIFTIG_KERNELS_0.37.0.md`](../../BENCHMARK_SWIFTIG_KERNELS_0.37.0.md).

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

## Reference metadata versus read assignment

FWR/CDR coordinate transfer for a custom germline is a reference-preparation operation, not part of read-to-allele mapping. Strict, permissive, and best-guess metadata matching select an annotated template from which to project V-region boundaries or J frame/CDR3-anchor information into `SWIGMETA`. They do not alter the custom allele's nucleotide sequence, add or remove reference alleles, or enter the V/D/J seed, affine-alignment, score, or call-selection calculations.

After assignment, the chosen allele's `SWIGMETA` is projected through its already-selected read alignment to create FWR/CDR sequences, junction/CDR3 boundaries, reading frame, translation, and productivity. Thus a permissive metadata choice can change those annotation-derived fields and analyses that subsequently filter/group on them, but cannot change the underlying V/D/J alignment or call. `swig-cli prepare-reference` persists this schema once; the associated diagnostics report the transfer template, identity, structural checks, and unresolved reason for every allele.

## Output and record detail

Swig writes an AIRR rearrangement table plus `swig_*` study/provenance fields. The interactive record view is a rendering of those committed fields: query and germline coordinates, FWR/CDR intervals, stitched V(D)J alignment, individual segment alignments, segment E-values, junction decomposition, and retained candidate evidence. It does not rerun assignment when opened.

The schema follows the AIRR Rearrangement representation so downstream tools can consume the table: [AIRR Community standardized representations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6173121/). Swig adds namespaced columns where the AIRR core schema has no equivalent. AIRR compatibility does not imply that SwiftIG is another AIRR-producing caller such as IgBLAST or IMGT/V-QUEST.

## Performance and cancellation

Reference indexes are built once per worker. Records are processed in bounded batches, and results are appended to browser storage or a user-selected writable stream. Cancellation aborts input reading, terminates the worker pool, aborts unfinalized output, and clears the temporary result store; no partial run is opened as a completed analysis.

AER/RIAT production indexes use direct-address CSR k-mer hits and reusable generation-stamped vote workspaces; their original unordered-map candidate implementation is retained for byte-equivalence tests. AIRR production output is appended directly with byte-identical numeric formatting, while the original ostream writer is likewise retained as a reference. Those kernel optimizations do not change scoring, search breadth, band choice, tie handling, or output fields. AER-R's explicitly documented search and tie differences are method changes, not kernel optimizations.

## Limitations

- AER/RIAT-MP/profile calibration and the AER-R development audit use simulated material and are not held-out validation.
- The AER-R simulator is realistic enough to expose boundary and pruning failures, but it is not a claim that its SHM, indel, trimming, or junction distributions match every repertoire.
- AER-R improved retained-D-tract recovery on the audited difficult cohort but also increased calls on truth records retaining fewer than six D bases; it remains experimental and opt-in.
- Short D segments and highly mutated reads can remain intrinsically ambiguous.
- Co-optimal labels describe equal scores under this caller, not equal biological posterior probabilities.
- C calls depend on whether the input actually covers constant sequence.
- Important study-specific calls should be benchmarked against an independent workflow and suitable truth/control material.
