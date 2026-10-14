\name{minerals_to_assay}
\alias{minerals_to_assay}
\title{Convert Ideal Mineral Mass Proportions to Theoretical Assays}
\usage{
minerals_to_assay(mineral_proportions, library = mineral_assay_library(),
  basis = c("fraction", "percent"),
  iron_reporting = c("separate", "total_Fe2O3"), closure_tolerance = 1e-6)
}
\arguments{
\item{mineral_proportions}{Finite nonnegative named numeric vector or numeric matrix/data frame with unique supported mineral column names. Rows must sum to the declared basis; there is no implicit normalization. Column order is respected and row names are retained.}
\item{library}{Named list of ideal formula entries with positive supported atom counts. Oxygen/volatile stoichiometry must be compatible with conventional oxides, water and carbon dioxide; oxygen-free sulfide/native formulas are supported.}
\item{basis}{Input mass fraction (sum 1) or mass percent (sum 100).}
\item{iron_reporting}{Separate FeO and Fe2O3, or total iron as Fe2O3 equivalent. The latter adds oxygen-equivalent mass; reported totals can exceed 100 percent.}
\item{closure_tolerance}{Finite nonnegative numeric scalar maximum input row-sum deviation in declared basis units.}
}
\description{Use ideal fixed stoichiometry to convert already closed mineral mass proportions into theoretical elemental and oxide reporting tables. These are bookkeeping outputs, with no instrument response or error model.}
\details{All tables have one row per input sample. The library supplies formula masses and oxidation-state atom counts. Elements include oxygen and merge Fe/Mn valence labels; one wt percent equals 10000 ppm by mass. Water and carbon dioxide equivalents remain separate from oxide tables. Physical components include physical oxides, water, carbon dioxide and nonoxide mineral mass.

The \code{mass_balance} data frame records \code{mineral_mass_wt_pct}, \code{element_mass_wt_pct}, \code{physical_oxides_plus_H2O_wt_pct}, \code{CO2_equiv_wt_pct}, \code{nonoxide_mineral_mass_wt_pct}, \code{physical_components_wt_pct}, \code{reported_components_wt_pct} and \code{reporting_oxygen_gain_wt_pct}. The last is reported component mass minus 100 wt percent, not measurement error. Oxide equivalents for nonoxide cations and total-iron reporting can produce totals above 100. Never force these equivalents to close or mix them with elemental ppm.}
\value{A list of aligned numeric data frames: elements (elemental wt percent), element_ppm (elements times 10000), oxides_physical (excluding nonoxide minerals), oxides_reported (cation equivalents), water (H2O equivalent), volatiles (H2O and CO2 equivalents), unreported_elements, iron_valence, and mass_balance; plus iron_reporting and assumptions. All chemistry except element_ppm is wt percent. Reported equivalents are never renormalized. Elemental and physical component mass are audited to 1e-6 wt percent against input mineral mass.}
\seealso{\code{\link{new_assay_scenario}}, \code{\link{simulate_assay_scenario}}, \code{\link{mineral_assay_library}}, \code{\link{simulate_mineral_assay_scenario}}}
