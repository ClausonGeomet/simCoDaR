\name{validate_assay_simulator}
\alias{validate_assay_simulator}
\title{Diagnose a fitted simulator against heldout positive assays}
\usage{
validate_assay_simulator(fit, heldout, reference = NULL,
  n_draws = 71L, replicates = 20L, seed = 20261014L)
}
\arguments{
\item{fit}{Valid ungrouped full positive support object from \code{fit_assay_simulator}.}
\item{heldout}{Complete finite positive same-basis assay data frame, at least two rows. Extra valid columns are ignored.}
\item{reference}{Optional positive same-basis training assays, at least two rows, used only for tail bounds.}
\item{n_draws}{Integer rows per draw, at least two.}
\item{replicates}{Integer repetitions, at least two.}
\item{seed}{Known nonnegative integer; seed plus replicates must be less than \code{.Machine$integer.max}. NULL is rejected. Repetition r uses seed+r with caller RNG restoration on success and error.}
}
\description{Evaluates an already fitted simulator without fitting, selecting or tuning it. Grouped fits and structural faces require a separately specified measure and weighting protocol and are rejected. Units and nonoverlap remain caller responsibilities.}
\details{
All pointwise log scores use Lebesgue volume in the fitted ordered orthonormal Helmert ILR basis; larger mean score is better. Dirichlet simplex density is converted by subtracting log J, where log J = -log(D)/2 - sum(log(x)).

Energy is 2 mean pairwise observed/generated ILR distances minus within-observed and within-generated mean distances, including diagonal zeros. Lower is better. All unordered logratios use type 7 quantiles at .05, .25, .5, .75 and .95; discrepancy is absolute quantile difference. Covariance error is relative Frobenius distance, with sample denominator n-1 and denominator floor 1e-10. Closed Fe2O3--MgO correlation is reported when named parts exist; constant-part correlations are unavailable with an explicit metadata reason. Tail fractions use strict inequalities below/above heldout .05/.95 bounds. Reference cell-tail fractions use per-part closed .01/.99 bounds; without reference they are unavailable.

Sampler clipping counts, invalid-draw count, cell validity counts and maximum closure residual are reported without corrections or discarded draws. Deterministic input/scoring errors stop; each sampling/diagnostic error is recorded and remaining repetitions proceed. For each scalar generated/discrepancy metric, MCSE is sd divided by sqrt(successful repetitions), unavailable below two successes. Any missing repetitions mark the estimate incomplete. Unavailable metrics also count as missing, with reasons in metadata. Draw MCSE measures conditional simulation noise, not uncertainty about a field population or scientific validation.
}
\value{A plain list with:
\describe{
\item{scores}{Data frame of row_index and log_score in heldout input order.}
\item{score_summary}{n, mean_log_score, measure and basis description.}
\item{replicate_metrics}{Long data frame: replicate, seed, metric, part_a, part_b, probability, observed, generated, discrepancy. Unused identifiers/values are NA. Energy, relative covariance error and aggregate quantile error are generated scalar statistics; their observed targets are zero.}
\item{mc_summary}{Metric identifiers plus statistic (generated or discrepancy), n_success, n_failure, mean, sd, mcse and incomplete. Structurally unused discrepancies are omitted.}
\item{failures}{Data frame: replicate, seed, message.}
\item{metadata}{Configuration, parts, units/basis, seeds, RNG kind, versions, unavailable reasons and limitations.}
}
No raw assay tables, calls, paths or identifiers are retained. No anonymity or differential privacy is promised.}
\seealso{\code{\link{fit_assay_simulator}}, \code{\link{simulate_assays}}}
