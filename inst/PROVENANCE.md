# Source provenance and modifications

simCoDaR is a GPL-3 adaptation of selected reusable source from CoDA
(source package declaration: GeoAssaySim).
Immutable source commit: `c2184f44994a7c3f2d4f36ec74c749e0267a4127`.
Source tree: `63e9eb5438c9b2a73b53385032ff7507a2582652`.
Source DESCRIPTION declares GPL-3 and attributes **Research Team**.
Inspected source history attributes **Matt Clauson**. Both origin declarations
are retained here; simCoDaR modifications are by Matt Clauson.
No separate source LICENSE/COPYING/NOTICE or conflicting notices were present
in the selected committed modules. No datasets or third-party assets were copied in the initial code recovery.

| Original relative source | Selected functions | Destination | Modifications |
| --- | --- | --- | --- |
| R/composition_study.R | .ilr_basis, .to_ilr, .from_ilr | R/utils-ilr.R | Extract only normalized Helmert ILR helpers; validate dimensions and finite intermediates; retain stabilized inverse floor -700. |
| R/scenario.R | new_assay_scenario, simulate_assay_scenario, .validate_composition_bounds | R/scenario.R | Remove totals, grades and templates; harden names, scalar inputs, closure, support and bounds; safe normalization; revalidate objects; scoped seeds with withr; keyed rows and assumptions. |
| R/mineral_conversion.R | mineral_assay_library, .mineral_entry, constants, .assay_formula_mass, .validate_mineral_library, .assay_oxide_molar_masses, minerals_to_assay | R/mineral_conversion.R | Retain atom counts, atomic weights, oxidation states, water/CO2 and oxygen audit; harden library/input names, types, finite masses and aligned output rows. |
| R/mineral_conversion.R | simulate_mineral_assay_scenario | R/mineral_scenario.R | Separate observed synthetic reporting blocks from full ideal truth; add domains, units, RNG policy, provenance and schema metadata; propagate draw budget. |

Classification of the initial recovery: the selected functions are reusable
numerical, explicit scenario, ideal chemistry and thin integration code.
Fitted models were excluded from that initial recovery and are added only in
the scoped empirical recovery below. Deposit templates, research analysis,
causal/process modules, upstream tests, launchers and data remain excluded.
There is no geological calibration or instrument measurement model.
All code remains GPL-3; no relicensing or source-side modification occurred.

## Empirical fitting recovery (Milestone 1)

Selected numerical families and group/support integration are additionally
recovered from immutable commit
`989577aaa449d28ac554d806ef57daa47122069c`, with the historical implementation
at `c2184f44994a7c3f2d4f36ec74c749e0267a4127` inspected for continuity.
The same GPL-3, Research Team origin declaration and Matt Clauson history
attribution apply. Modifications are by Matt Clauson. Original source was
inspected read-only; no source-side execution or modifications occurred.

| Original relative source | Selected functions | Destination | Modifications |
| --- | --- | --- | --- |
| R/composition_study.R | .mvnorm_logpdf/draw, .mvt_logpdf/draw, .logsumexp_rows, .dirichlet_logpdf/draw, .mixture_logpdf/draw | R/utils-empirical-numerics.R | Retain Cholesky whitening, lgamma/log1p, categorical component draw, stable row logsumexp. Reject individual zero gamma draws. Keep n-by-k density dimensions for single observations. |
| R/composition_study.R | .regularize_covariance, .fit_dirichlet, .fit_logistic_normal, .fit_student_t, .fit_gaussian_mixture | R/utils-empirical-numerics.R | Preserve moment initialization, bounded log-alpha optimization, full covariance/scale, declared ridge defaults, fixed-df weighted updates and multistart EM. Add analytical Dirichlet gradient and score, bound-saturation refusal, ridge/eigen/rank/condition diagnostics, explicit convergence and degeneracy checks. |
| R/simulator.R | .assay_pattern_key, .fit_assay_pattern, .draw_assay_pattern, fit_assay_simulator, simulate_assays | R/empirical.R, R/utils-empirical-contract.R | Preserve fitted group/support hierarchy and training BIC. New additive package signatures use scoped default seeds 42/123; explicitly validate common units/basis, names, groups, supports, controls and serializable fit parameters. Remove empirical row fallback, stored observed rows/training totals and independent amount generation. |
| R/composition_study.R | run_composition_study Jacobian scoring expression only | R/utils-empirical-density.R | Internal common ILR-volume density for all four families; explicit simplex conversion on active supports, discrete support/group mass and impossible-support -Inf. No study coordinator or evaluation data copied. |

Numerical changes are deliberate safeguards rather than new statistical
families. Mixture starts reaching iteration exhaustion or low component mass
are failures; only converged admissible starts can compete. Likelihood is
recomputed at updated states, and material decreases beyond the declared
relative tolerance fail with ridge involvement reported. Student exhaustion
also fails. Covariance regularization is applied once at its declared amount,
never escalated secretly. Gaussian covariance retains the source n-1 sample
convention; Student scale and df/(df-2) covariance are labelled separately.
BIC uses training rows only and smallest-k ties. Optim's Dirichlet diagnostic
provides function/gradient evaluation counts; an iteration count is unavailable.

The -700 inverse ILR floor remains unchanged in numerical meaning; the new
API reports affected cells and rows because it modifies extreme tails.
Selected closure is a subcomposition, not whole-rock mass closure. The fit
retains aggregate counts and necessary user group labels; it has no privacy
or anonymization guarantee. No total-composition model, censoring model,
unit conversion, physical mineral calibration or empirical dataset is added.
Source launchers, study orchestration, causal/process modules, private data,
research records and infrastructure remain excluded. Existing five scenario
and ideal-chemistry APIs and their numerical helpers are unchanged.

## Public data and bounded diagnostic addition (Milestone 2)

The separately licensed pinned IODP CSV is now bundled; attribution, selection
changes and byte identity are in extdata/IODP-LICENSE-PROVENANCE.md. Those data
are CC BY 4.0, not GPL-3. No other empirical dataset is included.

R/validate-empirical.R is a new bounded diagnostic interface by Matt Clauson,
using the recovered density/sampler unchanged. It reports common ILR scores,
energy, logratio quantiles/tails, covariance error, named-part correlation,
reference cell tails, draw validity and conditional draw MCSE. It does not fit,
tune or establish scientific validation. The installed companion example
freezes public split/config manifests before fitting. Historical numeric
outputs, primary start and RNG environment are unavailable: its results are
a new isolated-package reproduction, not recovered original numeric evidence.
