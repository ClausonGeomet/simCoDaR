\name{new_assay_scenario}
\alias{new_assay_scenario}
\title{Construct an Uncalibrated Composition Scenario}
\usage{
new_assay_scenario(endmembers, weights = NULL, logratio_sd = 0.25,
  structural_parts = NULL, composition_bounds = NULL)
}
\arguments{
\item{endmembers}{Finite nonnegative numeric matrix or data frame with unique named domain rows and at least two uniquely named part columns. Positive row totals are required. Parts must share a mass unit; mixing percent and ppm is unsupported. Rows are closed internally as ILR centers, not arithmetic mean guarantees.}
\item{weights}{Positive domain weights in row order, or named for every domain. Normalized safely to probabilities; NULL gives equal weights.}
\item{logratio_sd}{Nonnegative finite ILR standard deviation: an unnamed scalar or a vector for every domain. Named vectors must match all domain names exactly.}
\item{structural_parts}{NULL or distinct part names. Every zero cell requires its part declared here. Zeros remain exact; no pseudocount or censoring model is used.}
\item{composition_bounds}{NULL or a named list of part ranges c(lower, upper) in closed fractions within [0,1]. Centers must lie inside bounds. Support-specific closure must be feasible. Zero-width active bounds with positive spread and multiple active parts are rejected.}
}
\description{Define hypothetical domain compositions and their sampling probabilities, log-ratio spread, structural absence and fraction bounds. This is an explicit scenario, with no fitted geological calibration.}
\details{Normalized positive target parts are geometric centers in isometric log-ratio coordinates, not guaranteed arithmetic sample means. Isotropic normal variation is applied only to the active support in each domain. A part declared structural remains zero only where its target is zero; it can be positive in another domain. Domain labels and weights define independent categorical sampling, not fixed domain counts or spatial structure. No censored-zero or instrument model is supplied.}
\value{An \code{assay_scenario} list retaining endmembers, facies (domain names), parts, normalized weights, named logratio_sd, structural_parts, composition_bounds (lower/upper matrix), and calibrated = FALSE.}
\seealso{\code{\link{simulate_assay_scenario}}, \code{\link{mineral_assay_library}}, \code{\link{minerals_to_assay}}, \code{\link{simulate_mineral_assay_scenario}}}
