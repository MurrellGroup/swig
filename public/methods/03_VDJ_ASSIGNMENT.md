# V(D)J assignment and AIRR annotation

## Core execution

SwiftIG's C++20 core runs as WebAssembly in a bounded worker pool. Each worker receives the same composed germline references and processes independent record batches. Swig preserves input order when committing AIRR batches. Search strand, minimum accepted alignment identity, assignment strategy, and calling profile are explicit run settings.

For each oriented query, SwiftIG uses exact nucleotide seed indexes to retrieve candidates, exact affine local alignment for selected candidates, then reconstructs a consistent V-(D)-J rearrangement and AIRR coordinate/alignment fields. The ordinary RIAT-MP, AER, and standard paths evaluate D in the V/J-bounded region. Experimental AER-R begins there too, but can expand an uncertain junction to mapped conserved V/J anchors before comparing complete non-overlapping partitions. Co-optimal calls are emitted as comma-separated names; sparse near-tied candidates retain score/identity/coordinate evidence without duplicating full alignment strings into every AIRR row.

V and J are first optimized independently, but their junction-facing endpoints are not treated as irrevocable. If the ordinary D call has less than 10 consecutive exact bases and a V 3′ or J 5′ affine gap lies within 24 query bases of the junction-facing edge, SwiftIG checks the reclaimed boundary span for an exact 10-nt D seed. Only then does it perform a second D refinement. It clips the competing V/J alignment at each candidate D boundary, recomputes the retained affine scores, and accepts the split only when the non-overlapping V+D+J score is strictly greater than the original V+(weak D)+J score. This custom rescue prevents a gappy terminal V/J optimum from hiding a long exact D tract. Strong ordinary D calls and endpoints without a junction-facing gap never execute the additional search.

The minimum identity control is an acceptance floor for candidate alignments, not a global read-trimming heuristic. The `sequence` field remains the user's input sequence (or the chosen reverse complement); V/J-aligned and junction fields are derived outputs.

## Assignment strategies

RIAT-MP, ordinary AER, and standard SwiftIG change V candidate retrieval/refinement only. AER-R is a separately selectable experimental AER derivative that also changes how uncertain V/D/J boundaries and candidates are adjudicated. Every strategy uses the selected calling-profile scores except that R-optimized and Sensitive-D are intentionally restricted to AER-R because they also change AER-R's endpoint and D-presence decisions. Ordinary exact-run evidence floors remain unchanged; AER-R alone has the conservative distributed-D alternative described below.

AER-R with R-optimized is the default for all three front-page workflows, hand-written/normalized CLI configs, and the direct `swig-cli --vdj` route. Either surface can select any strategy explicitly, and a browser-exported config records the web selection.

### AER—Adaptive Exact Refinement

AER begins from the complete V-allele seed index. It exactly aligns the leading candidates and increases exact affine-alignment depth only when the leading 9-mer vote ranking is ambiguous: within 5% relative vote count or 8 weighted votes, up to 16 candidates. It never substitutes a propagated approximate score for the final retained exact alignment.

The production AER kernel evaluates independent affine recurrences in four exact SIMD128 lanes. A score-only pass over a wide safety pool omits trace storage but uses the same recurrence; full DP/traceback continues through every raw-score group capable of entering the retained top-N set. The scalar affine implementation remains the executable reference path used by equivalence tests.

### AER-R—experimental robust joint caller

AER-R leaves ordinary AER intact and is selected explicitly as `aer_robust` in portable config or `--assigner aer_robust` in the CLI. It reuses AER's exact seed candidates, affine recurrence, calling-profile scores, AIRR writer, and clean V/J-bounded D search. Its additional decisions are:

- compare explicit non-overlapping V–D–J partitions, charging the complete affine-score loss when V or J must be clipped;
- keep independently optimized V/J alignments immutable while candidates are evaluated; for an overlapping D candidate, compare retain-D/clip-V-or-J with retain-V/J/clip-D as separate complete partitions and commit boundaries only from the winner. The symmetric allocation must improve its matching and global partition strictly, except for the separately calibrated one-base V tie rule in R-optimized profiles;
- retain the nominal V/J-bounded D location and, only when it is absent, boundary-truncated, or weak beside a junction-facing gap, add a conserved-anchor-bounded D search as a distinct location hypothesis;
- map the final V cysteine and J F/W–G anchors from reference through the alignment, with a database-length-derived fallback only when metadata or alignment coverage cannot map an anchor;
- retain exact score ties at a configured candidate cutoff and add exact retained-D-substring allele labels that a top-N list alone could hide;
- when an already-scored D alignment fails the calling profile's consecutive-exact-run floor, retain it only if it has at least 16 aligned substitution columns, at least 14 nucleotide matches, at least 80% identity across those columns, and score at least 19 under the fixed `+2/−3/−13/−1` D evidence tuple; this lets distributed SHM evidence compete in the same joint partition without weakening the ordinary short-D rule;
- if a strongly V-supported orientation has no valid V/J partition, exactly align the complete small J database; the opposite/weak orientation does not pay for that rescue;
- call C only from the post-J slice, beginning at the single J/C-overlap base permitted by the reference convention; a C candidate can never clip or otherwise revise the selected J alignment; and
- when the separate double-D screen is enabled, permit a shorter seed rescue that must ungapped-extend back to at least the configured ordinary seed length before entering the unchanged pair-gain and pseudo-tandem tests.

For R-optimized profiles, AER-R also repairs a junction-facing local-alignment reset when the selected V identity is at most 0.88: a co-linear five-base V suffix or J prefix is materialized only when it has total substitution score zero under the selected `+2/−3` tuple, no shorter junction-adjacent prefix has positive score, and it lies entirely inside sequence not occupied by the actual selected segments. This decision uses the selected D boundary (or the opposite V/J boundary when D is absent), never an unselected D hypothesis. The five-base/identity condition improved held-out V-end and J-start losses in the supplied IgG-like simulation and left the lower-SHM J distribution unchanged.

The distributed rule is not gated on an empty D call. AER-R first traces the same leading strong-seed D candidates as ordinary AER. If none is a convincing long alignment—even when one or more short D calls exist—it applies exact four-lane SIMD score screening to the deferred safety pool through score 16, the conservative lower bound for any alignment capable of passing the aggregate rule under either shipped D-scoring tuple. If that still yields no convincing long D, AER-R forces the 3-mer tier across the complete locus-matched D set, scores only new gene/diagonal hypotheses, and traces only score-capable candidates. Every traced alignment that actually passes the aggregate rule is retained for joint V–D–J comparison even below the generic top-N cutoff.

This two-stage quality trigger avoids both brittle cases: “no D” is not a special boundary, and a short seed-driven call cannot exhaust the candidate/traceback budget. A deterministic adversarial regression inserts fifty perfect 7-nt D decoys ahead of the supplied mutated 17-nt KIMDB match; AER-R still selects the longer distributed match under the legacy truth/agreement profiles and R-optimized. An unconditional wide D scan slowed clean reads by about 44% and was rejected. A simpler eight-traceback retry added no correct calls and promoted many short chance matches, so it is also not shipped. Anchor projection and the maximum D-reference length are cached or evaluated lazily; convincing ordinary D calls pay none of the complete-set fallback.

The development simulator contains skewed allele use, V/D/J deletion, P/N addition, occasional tandem D, hotspot-weighted SHM, short indels, base-call errors, ambiguous bases, partial reads, flanks, reverse complements, fully trimmed-D controls, and exact base provenance. It and the compile-time all-reference candidate oracle live under `tests/`; neither is linked into the runtime WebAssembly. The oracle removes V/D/J candidate-count pruning but deliberately retains the production evidence rules, so it diagnoses candidate search rather than defining an unconstrained biological optimum. Results and limitations are in [`BENCHMARK_AER_ROBUST_0.37.7.md`](../../BENCHMARK_AER_ROBUST_0.37.7.md).

### RIAT-MP

RIAT-MP groups close V alleles into root-indexed trees. It aligns up to three representative roots and propagates sparse descendant differences without performing a full descendant alignment. When the provisional winner contains an indel, it tests at most two root traceback geometries within four raw-score units, with a 1,024-state traceback cap. This is a Swig/SwiftIG algorithm, not a literature package.

Production root alignments use the same exact four-lane SIMD kernel. Tree child links are precomputed once, and alternative partial tracebacks use parent-linked arena nodes; neither changes the path priority, traceback cap, or sparse descendant score. The complete 0.37.0 profile and reference-equivalence protocol are in [`BENCHMARK_SWIFTIG_KERNELS_0.37.0.md`](../../BENCHMARK_SWIFTIG_KERNELS_0.37.0.md).

### Standard SwiftIG

The standard strategy exactly aligns the three leading strong-seed V candidates. Weak-seed and seedless cases retain a larger safety pool. It is fixed-depth rather than ambiguity-adaptive.

## Calling profiles

The three legacy profiles are orthogonal to V strategy. R-optimized and Sensitive-D are AER-R-only and are rejected with every other strategy.

- **Truth-optimized:** the SwiftIG settings selected on the supplied simulated human-IGH truth data. Its ordinary D floor is a 6-nt exact run.
- **IgBLAST-agreement:** D scoring `+2/−4/−11/−1`, minimum 5-nt exact run and 3 candidates; J scoring `+2/−4/−13/−1` and 2 candidates. These values were selected for agreement with supplied IgBLAST calls.
- **IgBLAST-balanced:** the agreement settings plus removal of a D call only when its strongest support is exactly five consecutive matches and `j_sequence_start − v_sequence_end ≤ 11`. It maximized tuning-set IgBLAST agreement subject to mean first-call and ambiguity-aware V/D/J truth accuracy exceeding IgBLAST on that simulation.
- **R-optimized (experimental, AER-R only):** V `+2/−3/−9/−1`; D `+2/−3/−13/−1`, minimum 4-nt exact run, 2 candidates, and a 12-point cost for introducing D into the joint V–D–J partition. The cost is relaxed to 10 when the selected D hypothesis has raw score at least 20 or its exact ungapped tract occurs in at least two distinct locus-matched D template sequences. J uses `+2/−3/−17/−2`, 2 candidates. D alleles aligned to the same query span and within one raw alignment-score point of the selected D are emitted as a comma-separated uncertainty set. The 4-nt admission floor was promoted from 0.38.2 Sensitive-D because it improved fair D accuracy in both supplied simulations. The joint threshold was retained by uniform-set Brier score across the same simulations; it is not a posterior probability.
- **Sensitive-D (experimental, AER-R only):** the complete R-optimized scoring tuple, 4-nt exact-run floor, candidate depth, and ambiguity rule, but with a fixed 10-point joint D-presence cost instead of R-optimized's evidence-conditioned 12-to-10 cost. Candidate admission still requires an exact four-base run or the conservative distributed-D alternative, and the D-bearing complete partition must still beat the no-D partition. A bare four-base exact seed scores 8 under the D tuple and cannot pay the joint cost by itself. The profile is intended to expose the best supported D hypothesis in more weak-evidence junctions, not to force a D call on every read.

R-optimized's V gap-open cost prevents a local V alignment from cheaply opening a junction-facing gap and then consuming downstream N/D bases, while the restored −3 mismatch cost retains mutated high-SHM edge tracts. The D cost corrects the otherwise free advantage obtained by maximizing short matches across many alleles and junction positions. Its two-point relaxation prevents that correction from suppressing a strong single-template alignment or a weak tract supported by independent exact templates. Sensitive-D applies that relaxed cost to every already-admitted D candidate, which changes only the final joint D/no-D decision. Template-class support is cached per existing D hit and scans the small D reference only until two distinct sequences are found; it performs no additional alignment or candidate pass. The same-span, one-point D uncertainty merge scans already-ranked hits and changes labels only, not the selected boundary.

Calibration used the supplied 20,000-record lower-SHM and 20,000-record IgG-like KIMDB macaque simulations with the deduplicated union of bundled rhesus and cynomolgus V/D/J references. Calls were scored with a normalized multiclass Brier loss after collapsing identical full-reference sequences; a comma-separated report was interpreted as a uniform distribution over its unique sequence classes, and absence of D was an explicit class. Plain accuracy also reports the equivalent truth probability (`1/K` for `K` predicted sequence classes when truth is present). Boundary accuracy was evaluated separately, so a wider ambiguity list could not hide a bad endpoint. Tandem-D truths were excluded from the single-D primary call score. See the original [`BENCHMARK_R_OPTIMIZED_0.38.0.md`](../../BENCHMARK_R_OPTIMIZED_0.38.0.md), the cross-regime correction in [`RELEASE_NOTES_0.38.1.md`](../../RELEASE_NOTES_0.38.1.md), the 0.38.2 [`boundary audit`](../../BENCHMARK_SENSITIVE_D_0.38.2.md), and the 0.38.3 [`Sensitive-D decision audit`](../../BENCHMARK_SENSITIVE_D_0.38.3.md).

For every profile, AER-R's distributed-D rule is an alternative to—not a reduction of—the consecutive exact-run floor. It uses one fixed aggregate evidence tuple so changing calling profile does not silently change that rescue boundary. Ordinary AER, RIAT-MP, and standard assignment retain only their profile-specific exact-run rule.

Neither agreement profile is an IgBLAST implementation. IgBLAST uses its own BLAST-based germline search and annotation logic; see Ye et al., [IgBLAST](https://doi.org/10.1093/nar/gkt382). The profile names describe calibration objectives on one supplied simulation, not held-out validation or universal superiority.

## Segment support values

For every selected hit whose scoring tuple has been calibrated, SwiftIG writes the AIRR fields `v_support`, `d_support`, `j_support`, and `c_support`. These are BLAST-form expectation values: the expected number of chance local alignments with score at least as large as the reported raw segment score. Smaller is stronger. They are neither posterior probabilities nor repertoire-level evidence that an allele is present.

For raw score \(S\), SwiftIG evaluates

\[
E = K m' n' e^{-\lambda S}, \qquad
m' = m-\ell, \quad n'=n-N\ell,
\]

where \(m\) is the searched query length, \(n\) and \(N\) are the total nucleotide length and record count of the supplied segment database, and \(\ell\) is the standard finite-length adjustment fitted with \(\alpha\) and \(\beta\). D uses the actual V/J-bounded junction span searched; C uses only the post-J slice that was searched; V and J use the complete oriented query. Candidate indexing does not shrink the reference search space: support remains relative to every supplied allele in that segment database.

## Constant-region calling and isotype labels

C annotation is downstream of V(D)J assignment. Once J has been selected, SwiftIG searches from `j_sequence_end - 1`, allowing the one J-derived base that some C references use to complete the joining codon. It does not align long C references against the upstream V(D)J portion and does not give a C hypothesis permission to shorten J. Within the post-J slice, local alignment may stop before a mismatching primer, adapter, low-quality tail, or sequence outside the supplied C reference; downstream disagreement therefore does not invalidate an otherwise contiguous C match.

A reported C call always requires at least 90% identity over its local alignment and calibrated chance evidence. The ordinary local route requires at least 30 aligned query nucleotides and `c_support <= 1e-5`. A separate short-amplicon route is permitted only for a gapless hit near both expected origins (query offset at most 8 nt from the post-J search start and reference offset at most 3 nt from C base 1); it uses the number of locus-matched C references and allowed offsets as its search space and requires anchored `c_support <= 1e-4`. Twelve nucleotides is only the candidate-generation floor: the score/database-size gate determines whether any particular short tract is sufficient. Thus a small C database can support a shorter exact tract than a large one without pretending that 30 nt is a biological minimum.

These C defaults are independent of the general V/D/J minimum-identity control. An additional **covered C-prefix identity** gate is available but off by default. When enabled, it counts every base from inferred C-reference base 1 through the end of the retained local tract, including leading mismatches and indels that the positive-scoring local HSP could omit. A hit below the user-selected percentage is discarded. C is annotation-only: neither a provisional nor an accepted C candidate contributes to V/D/J boundary, call, strand, or orientation selection, so failing this optional gate clears only C fields. AIRR output reports the measured value as `c_prefix_identity`.

The bundled reference builder assembles one secreted/coding path per IMGT C allele and excludes membrane-only and UTR exons. Consequently, the caller can distinguish only the constant genes/alleles represented in the selected C FASTA and only when the sequenced tract covers discriminating bases. Co-optimal C names remain comma-separated in AIRR output. The derived browser `isotype` label accepts the same long evidence route or a short supported route (`c_support <= 1e-4`); ambiguity across IgG or IgA subclasses is conservatively collapsed to `IgG` or `IgA` rather than choosing the first subclass arbitrarily. Secreted-versus-membrane isoform assignment requires a reference FASTA that explicitly contains the relevant alternative exon paths and reads that cover them; it is not inferred from an unobserved tail.

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

- AER/RIAT-MP/profile calibration and the AER-R development audit use simulated material and are not external validation. R-optimized has a deterministic internal test partition, but it remains the same supplied simulator/database regime.
- The AER-R simulator is realistic enough to expose boundary and pruning failures, but it is not a claim that its SHM, indel, trimming, or junction distributions match every repertoire.
- AER-R improved retained-D-tract recovery on the audited difficult cohort but its broader joint-boundary behavior still calls many truth records retaining fewer than six D bases; it remains experimental and opt-in. The complete 0.37.7 change left that rate unchanged on the primary cohort and added one call among 5,000 fully trimmed-D controls. The large pre-existing short-exact-match call rate in the latter remains a separate model-selection limitation.
- R-optimized and Sensitive-D were calibrated on two supplied simulator regimes but are not external biological validation. Very short retained D tracts remain difficult to distinguish from chance junction matches. Sensitive-D deliberately accepts a materially larger no-D false-positive rate and an approximately 0.2-percentage-point fair D-accuracy loss on each supplied simulation in exchange for higher D detection; use it only when that recovery/specificity trade-off is appropriate.
- Short D segments and highly mutated reads can remain intrinsically ambiguous.
- Co-optimal labels describe equal scores under this caller, not equal biological posterior probabilities.
- C calls depend on sufficient calibrated post-J evidence and whether the selected C reference contains the relevant gene/isoform path. Very short reads may still be insufficient, especially with a large C database; a fixed minimum biological length is not asserted.
- Important study-specific calls should be benchmarked against an independent workflow and suitable truth/control material.

## Assignment progress

Multi-file progress uses read counts only when all inputs have known counts. Otherwise every file uses the same expanded-file-size workload proxy; an unknown count is never assigned a one-read weight. Gzip footer sizes are read without a counting pass, with a rough compressed-size expansion fallback for large/wrapped sizes.

Within a file, parsed character progress is discounted by the fraction of parsed records actually acknowledged by the AIRR store. Compressed read-ahead and dispatched batches are not completed assignments. End-of-input switches to the exact committed/parsed count. Subsampling separates the scan from selected-read commits. Index finalization reserves the last percentage; 100% is shown only when results are ready. These are approximate work fractions, not elapsed-time predictions.
