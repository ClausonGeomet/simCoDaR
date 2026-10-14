\name{mineral_assay_library}
\alias{mineral_assay_library}
\title{Ideal Fixed Formula Mineral Library}
\usage{
mineral_assay_library()
}
\description{Return ideal fixed-formula minerals for deterministic mass-to-chemistry conversion. Formula labels are descriptive; explicit atom counts drive the calculations. No solid-solution fit or geological calibration is supplied.}
\value{A named list of ideal minerals; each entry has a formula label and a named positive numeric atoms vector. Explicit atom counts drive calculations; labels are not parsed. Fe2/Fe3 and Mn2/Mn3/Mn4 specify oxidation states. Variable solid solutions are not fitted.}
\seealso{\code{\link{new_assay_scenario}}, \code{\link{simulate_assay_scenario}}, \code{\link{minerals_to_assay}}, \code{\link{simulate_mineral_assay_scenario}}}
