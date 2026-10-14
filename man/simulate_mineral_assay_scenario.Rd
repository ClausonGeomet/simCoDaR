\name{simulate_mineral_assay_scenario}
\alias{simulate_mineral_assay_scenario}
\title{Draw Minerals and Separate Reporting Features from Truth}
\usage{
simulate_mineral_assay_scenario(spec, n, seed = NULL,
  library = mineral_assay_library(),
  iron_reporting = c("separate", "total_Fe2O3"), max_draw_attempts = 1000L)
}
\arguments{
\item{spec}{An assay_scenario whose parts all name minerals in library.}
\item{n}{Positive numeric integer scalar number of samples.}
\item{seed}{NULL to consume caller RNG, or a scoped integer seed in [0, .Machine$integer.max - 1].}
\item{library}{Ideal mineral formula library; see mineral_assay_library.}
\item{iron_reporting}{Separate physical iron oxides or total Fe2O3 reporting equivalent.}
\item{max_draw_attempts}{Positive numeric integer scalar rejection-round budget passed to the composition draw.}
}
\description{Draw one hypothetical mineral composition dataset and convert it into theoretical reporting tables, retaining domain labels, mineral truth and reproducibility metadata. No geological calibration or instrument measurement model is supplied.}
\details{Samples are independent conditional on domain. The observed field contains deterministic ideal reporting features, with no added measurement noise. Mineral fractions close to one, while reported oxide equivalents are not renormalized. Water and carbon dioxide equivalents and the oxygen-gain audit are available in \code{truth$assay}; see \code{\link{minerals_to_assay}}.

Explicit seeds restore the caller's RNG state on success or error, including when initially absent. A NULL seed consumes the caller's stream. RNG kind is unchanged and recorded. Reproducibility depends on package/R versions and RNG kind; cross-version bitwise identity is not promised.}
\value{A list with sample_id, domains (sample_id and domain columns), observed, truth and metadata. Observed contains xrd_phase_percent (ideal phase mass percent, not diffraction intensities), xrf_oxide_wt_percent (ideal reported oxide equivalents), and icp_element_ppm (ideal elemental ppm by mass). Truth contains mineral_fraction and the full assay converter result. Metadata records schema/package versions, units, seed/RNGkind/policy, iron reporting, source attribution, scenario assumptions, calibrated = FALSE, measurement_model = "none", and truth_equals_ideal_reporting = TRUE. All sample tables are aligned by sample ID.}
\examples{
target <- rbind(a = c(Hematite = 0.7, Quartz = 0.3),
                b = c(Hematite = 0.3, Quartz = 0.7))
x <- simulate_mineral_assay_scenario(new_assay_scenario(target), 10, seed = 42)
x$observed$icp_element_ppm
}
\seealso{\code{\link{new_assay_scenario}}, \code{\link{simulate_assay_scenario}}, \code{\link{mineral_assay_library}}, \code{\link{minerals_to_assay}}}
