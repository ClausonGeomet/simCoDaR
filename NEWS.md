# simCoDaR 0.1.0.9000

* Initial development API: five functions for explicit hypothetical domain
  compositions and deterministic ideal mineral chemistry.
* Support-specific ILR draws preserve declared structural zeros. Fraction
  bounds use a finite rejection budget; explicit seeds restore caller RNG state.
* Separate sample IDs/domain labels, synthetic reporting tables, mineral truth
  and metadata, with mineral mass %, theoretical oxide wt% and elemental ppm.
* Physical/elemental mass audits and oxygen-equivalent reporting differences are
  retained; reporting equivalents are never renormalized.
* Standalone 300-sample two-domain example, reproducible PNG/SVG figures and
  expanded usage/assumption documentation. Scenarios remain uncalibrated with
  no instrument model; scenarios remain illustrative despite numerical API checks.
* GPL-3 source attribution preserved in inst/PROVENANCE.md.

* Add empirical `fit_assay_simulator()` and `simulate_assays()` for Dirichlet,
  full-covariance logistic normal, fixed-df full-scale ILR Student t and multistart
  full-covariance ILR mixtures with training-only BIC. Preserve five existing APIs.
* Explicit measurement/structural-support contracts, scoped RNG, convergence and
  inverse-floor diagnostics; no total model, censoring replacement or raw-row fallback.
* Add `validate_assay_simulator()` for positive ungrouped heldout common-ILR scores,
  energy/logratio/tail/covariance diagnostics, explicit failures and conditional MCSE.
* Bundle the pinned CC BY 4.0 IODP U1418 data with attribution, extraction notices
  and a five-split installed reproduction companion. Historical numeric outputs are
  unavailable; new model rankings are empirical observations, not acceptance criteria.
* Add executable measured-assays vignette, dependence and raincloud PNG/SVG figures,
  two-workflow README and metadata-only six-source catalogue with qualified rights,
  measurement limits and corrected Cannington Pb143 numeric/10 nondetect counts.
