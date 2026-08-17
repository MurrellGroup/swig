# SwiftIG AIRR support calibration — Swig 0.34.0

## Scope and interpretation

Swig 0.34.0 adds `v_support`, `d_support`, `j_support`, and `c_support` to AIRR output from the shared SwiftIG WebAssembly core. The same implementation therefore serves the browser, the full pipeline CLI, and the streaming `swig-cli --vdj` route.

Each populated value is an expectation value for the selected SwiftIG local alignment:

\[
E(S)=K(m-\ell)(n-N\ell)e^{-\lambda S}.
\]

It estimates the number of chance local alignments scoring at least \(S\) against the supplied segment reference. It is not an allele posterior or a probability that the call is correct. It uses a BLAST-form model, but it is not copied from IgBLAST: SwiftIG's affine costs do not match a published NCBI parameter row, and its candidate search and segment query spans differ from IgBLAST's. Reporting IgBLAST's constants against SwiftIG raw scores would therefore create precise-looking but invalid numbers.

## Reproducible calibration

`scripts/calibrate-airr-support.cpp` implements a score-only copy of SwiftIG's exact affine local-alignment recurrence. With fixed seed `0x535749474556414c`, it generated independent uniform-DNA query/reference pairs at 15 matrix dimensions from 20 × 20 through 512 × 320. There were 121,480 null alignments per scoring tuple. `scripts/fit-airr-support.py` fitted the extreme-value model and finite-size adjustment; one orientation of each rectangular matrix was fitted and the independently simulated transpose was held out. The scripts require a C++20 compiler and SciPy only for release-time calibration. Neither is used by Swig Web or the distributed CLI.

The embedded constants are:

| SwiftIG use | Match / mismatch / gap-open / gap-extend | λ | K | α | β |
|---|---:|---:|---:|---:|---:|
| V and C, all shipped profiles | `2 / −3 / −5 / −1` | 0.615091627 | 0.657543880 | 2.292298840 | −9.322674236 |
| D, truth-optimized | `2 / −3 / −13 / −1` | 0.633731437 | 0.382178144 | 0.383131231 | 0.317757739 |
| J, truth-optimized | `2 / −3 / −17 / −2` | 0.633731437 | 0.381864932 | 0.379080053 | 0.341396552 |
| D, IgBLAST-agreement/balanced | `2 / −4 / −11 / −1` | 0.645352028 | 0.280867619 | 0.230455685 | −3.339067764 |
| J, IgBLAST-agreement/balanced | `2 / −4 / −13 / −1` | 0.645457955 | 0.281099152 | 0.231233106 | −3.350080245 |

The fitted gapped \(\lambda\) was constrained not to exceed the corresponding ungapped theoretical limit. Across 45,740 held-out alignments per tuple, root-mean-square error in log expected exceedance count over central null thresholds was 0.16 for truth D/J, 0.36 for agreement D/J, and 0.24 for V/C. In held-out upper tails with empirical survival from 0.1% to 20%, median predicted/empirical expected-count ratios were 0.98, 0.99, and 0.77 respectively. Discrete short-alignment scores make individual thresholds stepwise; these checks support an approximate E-value interpretation, not false exactness in the extreme unobserved tail.

Runtime uses the actual record count and total nucleotide length of each supplied V/D/J/C reference. These totals are cached when a worker indexes the references. V/J/C use the complete oriented query length; D uses the V/J-bounded query span actually searched. An unknown scoring tuple leaves support empty.

## Runtime check

The repository's 50,000-read human-IGH WASM benchmark used 543 alleles, four workers, and 1,000-record batches. Two unchanged-core runs took 26.741 and 27.568 seconds. Three support-enabled runs took 26.146, 28.234, and 27.255 seconds. Their means were 27.155 and 27.212 seconds, a 0.21% difference within run-to-run noise; no additional alignment or database scan occurs. AIRR output grew from 149.5 to 151.4 MiB (1.27%) because of the four fields. Final peak RSS was 330.2 MiB versus 321.9–327.9 MiB in the two baseline runs; allocator high-water measurements vary, and the implementation adds no repertoire-sized state.

## Relationship to prior work

The expectation-value form follows Karlin and Altschul, [*Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes*](https://doi.org/10.1073/pnas.87.6.2264), and the gapped BLAST framework of Altschul et al., [*Gapped BLAST and PSI-BLAST*](https://doi.org/10.1093/nar/25.17.3389). IgBLAST applies BLAST machinery to immunoglobulin and T-cell receptor queries (Ye et al., [IgBLAST](https://doi.org/10.1093/nar/gkt382)). SwiftIG borrows the statistical interpretation and AIRR field convention, not IgBLAST's complete search implementation or numerical parameter tables.

## Limitations

- The null assumes independent equiprobable A/C/G/T; composition-biased reads or references can be less well calibrated.
- Extremely small reported values extrapolate the fitted extreme-value tail beyond directly simulatable frequencies.
- Co-optimal comma-separated alleles share the selected alignment's value; support does not apportion probability among them.
- Support is tied to a segment alignment score. It does not include V–J pairing bonuses, D exact-run gates, annotation metadata, productivity, or repertoire-level allele evidence.
- Numerical comparisons across tools require identical scoring, databases, query ranges, and search semantics. Those conditions do not generally hold between SwiftIG and IgBLAST.
