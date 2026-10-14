\name{simulate_assay_scenario}
\alias{simulate_assay_scenario}
\title{Draw Synthetic Domain Compositions}
\usage{
simulate_assay_scenario(spec, n, seed = NULL, max_draw_attempts = 1000L)
}
\arguments{
\item{spec}{An assay_scenario constructed by new_assay_scenario. Fields and labels are revalidated on every draw.}
\item{n}{Numeric integer scalar between 1 and .Machine$integer.max.}
\item{seed}{NULL or numeric integer scalar from 0 to .Machine$integer.max - 1. Explicit seeds restore caller RNG state even on error; NULL consumes it normally. RNGkind is unchanged. Cross-R-version bitwise identity is not promised.}
\item{max_draw_attempts}{Positive integer scalar rejection-round budget per domain. Exhaustion stops with a domain/bounds message without relaxing bounds.}
}
\description{Draw independent domain labels and closed compositions from an explicit hypothetical scenario. No deposit calibration or measurement model is supplied.}
\value{A list with sample_id (sample_1 through sample_n), composition (numeric data frame of closed fractions, rownames equal to IDs), facies (domain labels), and assumptions (validated scenario, seed, RNGkind, rejection rate and RNG policy).}
\details{Within each sampled domain the positive support uses a normalized Helmert isometric log-ratio basis and independent normal coordinate variation. Structural zeros remain absent. Single-part support and zero spread are deterministic. Bounds use inclusive comparisons with 1e-12 slack; fraction closure is checked to 1e-10.}
\seealso{\code{\link{new_assay_scenario}}, \code{\link{mineral_assay_library}}, \code{\link{minerals_to_assay}}, \code{\link{simulate_mineral_assay_scenario}}}
