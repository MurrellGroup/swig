# Swig 0.38.8

Personalized germline discovery now distinguishes an unmutated lineage component, flexible gene-specific SHM and local coordinated mutation events from inherited linked haplotypes. Terminal proposals must also reject flexible trimming/junction explanations. Partial-coverage competition and stable log likelihoods fix concrete false-call and rare-recovery defects.

The 5% novel-fraction floor is removed; linked support can include additional SHM. The candidate budget is exposed in Web and CLI (`--max-candidates`). Discovery evidence is reported separately from final set-selection gain. No hard diploid allele cap is imposed.

Across two independently processed macaques and two six-allele panels per animal, 22/24 targets are recovered exactly and two at all identifiable interior sites. Original-panel additional calls fall from 10/8 to 5/3. Remaining calls are unvalidated; their persistence with full KIMDB does not establish specificity. See BENCHMARK_PERSONALIZED_GERMLINE_0.38.8.md for results, limitations and replay data.

AER-R/R-optimized assignment defaults and multi-file progress estimates from 0.38.7 are retained. Full tests (271) and production build pass. This is a complete source release with generated CLI, WASM assets and GitHub Actions binary-release workflow.
