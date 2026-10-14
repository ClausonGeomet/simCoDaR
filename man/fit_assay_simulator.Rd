\name{fit_assay_simulator}
\alias{fit_assay_simulator}
\title{Fit an Empirical Simulator for Selected Assay Subcompositions}
\usage{
fit_assay_simulator(data, parts,
  model = c("logistic_normal", "ilr_student_t", "ilr_mixture", "dirichlet"),
  group_col = NULL, seed = 42, zero_mode = c("error", "structural"),
  structural_parts = NULL, min_pattern_n = 30L, mixture_components = 1:3,
  student_df = 5, unit_label = NULL, part_units = NULL,
  measurement_basis = "unspecified", control = list())
}
\arguments{
\item{data}{Nonempty data frame with complete finite numeric selected parts.}
\item{parts}{Distinct selected column names in desired output order.}
\item{model}{One of the four fitted families listed above.}
\item{group_col}{Optional separate complete geological group column.}
\item{seed}{Integer from zero to .Machine$integer.max-1, or NULL to consume caller RNG.}
\item{zero_mode}{Reject zeros, or permit explicitly declared structural absence.}
\item{structural_parts}{Distinct selected names whose zeros are confirmed structural.}
\item{min_pattern_n}{Minimum multivariate support count, also subject to p+2 rows for p ILR dimensions.}
\item{mixture_components}{Positive integer candidate counts selected by training BIC.}
\item{student_df}{Fixed predeclared Student degrees of freedom greater than two.}
\item{unit_label}{Descriptive common unit, without conversion.}
\item{part_units}{Optional named common unit for every selected part.}
\item{measurement_basis}{One of unspecified, oxide_mass, element_mass, mineral_mass, mineral_volume, abstract_composition.}
\item{control}{Named iteration/tolerance overrides described below.}
}
\description{Recover four compositional fitting methods on a caller-confirmed common,
nonoverlapping measurement basis. Selected amounts are closed to one. No
whole-rock mass balance, total model or censoring replacement is inferred.}
\details{
Controls are dirichlet_maxit=2000, student_max_iter=500,
mixture_max_iter=300, mixture_starts=5 and tolerance=1e-7. Iteration counts
must be positive integers and tolerance must be in (0, 1e-3]. Sparse or
identical multivariate supports and nonconvergence fail without empirical
row-resampling fallback. Singleton supports are deterministic.

Dirichlet fits positive alpha by bounded log-alpha maximum likelihood;
bound saturation fails. Distinct-part Dirichlet covariances are negative.
Logistic normal uses full sample covariance with denominator n-1 plus ridge,
so this is not the exact unpenalized covariance MLE. Student sigma is a full
scale matrix, with covariance df/(df-2) times sigma. The fixed df is never
selected on evaluation data. Mixture uses multiple starts, full component
covariance and minimum training BIC; failures are reported in diagnostics.
Declared ridges are not escalated silently. Each support records its ordered
normalized Helmert ILR basis and convergence/regularization diagnostics.

Supplied seeds are scoped with withr, restoring caller RNG even on errors.
NULL consumes the caller RNG. Recognized overlapping elemental/oxide or
iron-equivalent names and mixed recognized units are rejected. Arbitrary
names cannot prove chemical nonoverlap; callers must resolve unidentified
species. Unspecified basis records an unverified caller assumption. No
unit conversion or mineral-volume to mass inference is performed.

The object retains parameters and aggregate counts, not observed rows,
identifiers or totals. Necessary user group labels are retained; there is
no anonymization or differential privacy guarantee. No total-composition
dependence is modelled. Zeros from censoring, rounding or missingness require
separate treatment and cannot be automatically declared structural.
}
\value{A serializable geo_assay_simulator list with schema version, ordered
parts/groups/supports, fitted parameters and basis, source summary, controls,
RNG/version metadata and assumptions.}
\examples{
d <- data.frame(a = exp(seq(-1, 1, length.out = 40)), b = 1)
f <- fit_assay_simulator(d, c("a", "b"), measurement_basis = "abstract_composition")
simulate_assays(f, 3)
}
\seealso{\code{\link{simulate_assays}}}
