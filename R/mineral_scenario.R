# Adapted wrapper from CoDA R/mineral_conversion.R at
# c2184f44994a7c3f2d4f36ec74c749e0267a4127; GPL-3.
# Origin: Research Team; history: Matt Clauson. See inst/PROVENANCE.md.
#' Draw ideal mineral fractions and separate reporting features from truth.
#' @export
simulate_mineral_assay_scenario <- function(spec, n, seed = NULL,
    library = mineral_assay_library(), iron_reporting = c("separate", "total_Fe2O3"),
    max_draw_attempts = 1000L) {
  spec <- .revalidate_scenario(spec)
  iron_reporting <- match.arg(iron_reporting)
  .validate_mineral_library(library)
  if (!all(spec$parts %in% names(library))) stop("Every scenario part must name a mineral in library.")
  draw <- simulate_assay_scenario(spec, n, seed, max_draw_attempts)
  assay <- minerals_to_assay(draw$composition, library, "fraction", iron_reporting)
  list(sample_id = draw$sample_id,
    domains = data.frame(sample_id = draw$sample_id, domain = draw$facies,
                         row.names = draw$sample_id),
    observed = list(xrd_phase_percent = 100 * draw$composition,
      xrf_oxide_wt_percent = assay$oxides_reported, icp_element_ppm = assay$element_ppm),
    truth = list(mineral_fraction = draw$composition, assay = assay),
    metadata = list(schema_version = "0.1", package_version = "0.1.0.9000",
      seed = draw$assumptions$seed, RNGkind = draw$assumptions$RNGkind,
      rng_policy = draw$assumptions$rng_policy,
      units = c(xrd_phase_percent = "mass percent", xrf_oxide_wt_percent = "wt percent",
                icp_element_ppm = "ppm by mass"),
      iron_reporting = iron_reporting,
      source_attribution = list(origin = "Research Team", history_author = "Matt Clauson",
        commit = "c2184f44994a7c3f2d4f36ec74c749e0267a4127", license = "GPL-3"),
      scenario_assumptions = draw$assumptions, calibrated = FALSE,
      measurement_model = "none", truth_equals_ideal_reporting = TRUE))
}
