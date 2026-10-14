\name{simulate_assays}
\alias{simulate_assays}
\title{Draw Synthetic Selected-Part Proportions from a Fitted Simulator}
\usage{simulate_assays(fit, n, seed = 123, group = NULL, output = "proportion")}
\arguments{
\item{fit}{A valid geo_assay_simulator from fit_assay_simulator.}
\item{n}{Positive integer number of rows.}
\item{seed}{Integer from zero to .Machine$integer.max-1, or NULL to consume caller RNG.}
\item{group}{Optional fitted group label, length one or n. Omission samples fitted group frequencies.}
\item{output}{Only proportion is supported; total generation is not included.}
}
\description{Revalidate the fitted schema and draw from full fitted family
parameters on declared active supports, preserving structural absence.}
\details{Supplied seeds restore the caller RNG on success and error via withr.
Inverse ILR subtracts each row maximum and floors shifted log values at -700;
reported row/cell floor counts identify modification of extreme tails.
Invalid gamma entries or nonfinite draws fail rather than receive pseudocounts.
Structural supports and groups are sampled from their fitted discrete
probabilities. Parts are returned in the original caller order.}
\value{A data frame of closed fractions with synthetic row names and an
optional group column after parts. Attributes assumptions and diagnostics
record limits, seed policy, inverse floor events and invalid draw count
(zero for a successfully returned sample).}
\seealso{\code{\link{fit_assay_simulator}}}
