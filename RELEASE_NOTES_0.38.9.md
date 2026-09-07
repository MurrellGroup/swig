# Swig 0.38.9

Adds a separate **Joint inherited / SHM model** box alongside the existing personalized-reference method. Either or both can be run; their workers, results, downloads and saved-session fields are independent. The existing inference algorithm is preserved.

The new experimental method uses one normalized latent-exposure likelihood for inherited and somatic sequence hypotheses, a four-base CTMC observation kernel, raw V-end marginalization, training-only proposals, and held-out split-likelihood evidence. A concave-fitting bound prevents unfinished null optimization from inflating evidence. It reports evidence rather than exporting a supposedly certain personalized genotype.

This is **not yet a superior germline caller**. Validation shows substantial loss of rare-allele sensitivity. The benchmark reports nomination, positive evidence and the multiple-testing evidence rule separately. Rejecting more hypotheses is not presented as improved accuracy.

The new CLI command is `swig-cli joint-germline`; `swig-cli personalized-germline` retains the existing engine. AER-R/R-optimized assignment defaults and multi-file progress fixes remain in place.

See `public/methods/10_JOINT_INHERITED_SHM_MODEL.md` for the model, assumptions and interpretation, and `BENCHMARK_JOINT_GERMLINE_0.38.9.md` for the empirical comparison. This complete source release includes generated CLI, WASM assets, both methods, tests, replay inputs and the GitHub Actions release workflow.
