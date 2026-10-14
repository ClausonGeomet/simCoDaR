# Adapted from CoDA R/mineral_conversion.R at c2184f44994a7c3f2d4f36ec74c749e0267a4127; GPL-3.
# Origin declaration: Research Team; history attribution: Matt Clauson.
# See inst/PROVENANCE.md for modifications.
# Ideal formula-mineral to bulk element/oxide assay conversion. Formula strings
# are labels; explicit atom counts drive all calculations.

.assay_atomic_weights <- c(H = 1.008, O = 15.999, Na = 22.98977,
                           Mg = 24.305, Al = 26.98154, Si = 28.085,
                           P = 30.97376, S = 32.06, K = 39.0983, Ca = 40.078,
                           Ti = 47.867, Cr = 51.9961, Mn = 54.93804,
                           Fe = 55.845, Ni = 58.6934, Cu = 63.546,
                           Zn = 65.38, As = 74.92160, Ag = 107.8682,
                           Pb = 207.2, Au = 196.96657, C = 12.011)

.assay_oxide_formula <- list(
  Na2O = c(Na = 2, O = 1), MgO = c(Mg = 1, O = 1),
  Al2O3 = c(Al = 2, O = 3), SiO2 = c(Si = 1, O = 2),
  P2O5 = c(P = 2, O = 5), K2O = c(K = 2, O = 1),
  CaO = c(Ca = 1, O = 1), TiO2 = c(Ti = 1, O = 2),
  Cr2O3 = c(Cr = 2, O = 3), MnO = c(Mn = 1, O = 1),
  Mn2O3 = c(Mn = 2, O = 3), MnO2 = c(Mn = 1, O = 2),
  FeO = c(Fe = 1, O = 1), Fe2O3 = c(Fe = 2, O = 3),
  NiO = c(Ni = 1, O = 1), CuO = c(Cu = 1, O = 1),
  ZnO = c(Zn = 1, O = 1), PbO = c(Pb = 1, O = 1),
  H2O = c(H = 2, O = 1), CO2 = c(C = 1, O = 2)
)

.assay_oxide_cation <- c(Na2O = "Na", MgO = "Mg", Al2O3 = "Al",
                          SiO2 = "Si", P2O5 = "P", K2O = "K",
                          CaO = "Ca", TiO2 = "Ti", Cr2O3 = "Cr",
                          MnO = "Mn2", Mn2O3 = "Mn3", MnO2 = "Mn4",
                          FeO = "Fe2", Fe2O3 = "Fe3", NiO = "Ni",
                          CuO = "Cu", ZnO = "Zn", PbO = "Pb",
                          H2O = "H", CO2 = "C")

.mineral_entry <- function(formula, atoms) list(formula = formula, atoms = atoms)

#' Ideal fixed-formula minerals supported by the assay converter.
#'
#' Each entry is `list(formula = <label>, atoms = <named numeric vector>)`.
#' Iron atom names `Fe2` and `Fe3` denote oxidation states; ordinary `Fe` is
#' deliberately unsupported. Variable-composition solid solutions require a
#' user-supplied idealized stoichiometry.
#' @return A named list of ideal mineral formula entries.
mineral_assay_library <- function() {
  list(
    Hematite = .mineral_entry("Fe2O3", c(Fe3 = 2, O = 3)),
    Goethite = .mineral_entry("FeO(OH)", c(Fe3 = 1, O = 2, H = 1)),
    Magnetite = .mineral_entry("Fe3O4 = Fe(II)Fe(III)2O4",
                                c(Fe2 = 1, Fe3 = 2, O = 4)),
    Quartz = .mineral_entry("SiO2", c(Si = 1, O = 2)),
    Kaolinite = .mineral_entry("Al2Si2O5(OH)4",
                                c(Al = 2, Si = 2, O = 9, H = 4)),
    Forsterite = .mineral_entry("Mg2SiO4", c(Mg = 2, Si = 1, O = 4)),
    Fayalite = .mineral_entry("Fe2SiO4", c(Fe2 = 2, Si = 1, O = 4)),
    Enstatite = .mineral_entry("MgSiO3", c(Mg = 1, Si = 1, O = 3)),
    Lizardite = .mineral_entry("Mg3Si2O5(OH)4",
                                c(Mg = 3, Si = 2, O = 9, H = 4)),
    Brucite = .mineral_entry("Mg(OH)2", c(Mg = 1, O = 2, H = 2)),
    Chromite = .mineral_entry("FeCr2O4", c(Fe2 = 1, Cr = 2, O = 4)),
    Bunsenite = .mineral_entry("NiO", c(Ni = 1, O = 1)),
    Pyrolusite = .mineral_entry("MnO2", c(Mn4 = 1, O = 2)),
    Manganite = .mineral_entry("MnO(OH)", c(Mn3 = 1, O = 2, H = 1)),
    Rhodochrosite = .mineral_entry("MnCO3", c(Mn2 = 1, C = 1, O = 3)),
    Galena = .mineral_entry("PbS", c(Pb = 1, S = 1)),
    Sphalerite = .mineral_entry("ZnS", c(Zn = 1, S = 1)),
    Pyrite = .mineral_entry("FeS2", c(Fe2 = 1, S = 2)),
    Chalcopyrite = .mineral_entry("CuFeS2", c(Cu = 1, Fe3 = 1, S = 2)),
    Chalcocite = .mineral_entry("Cu2S", c(Cu = 2, S = 1)),
    Arsenopyrite = .mineral_entry("FeAsS", c(Fe2 = 1, As = 1, S = 1)),
    Acanthite = .mineral_entry("Ag2S", c(Ag = 2, S = 1)),
    NativeGold = .mineral_entry("Au", c(Au = 1)),
    NativeSilver = .mineral_entry("Ag", c(Ag = 1)),
    Calcite = .mineral_entry("CaCO3", c(Ca = 1, C = 1, O = 3)),
    Cerussite = .mineral_entry("PbCO3", c(Pb = 1, C = 1, O = 3)),
    Smithsonite = .mineral_entry("ZnCO3", c(Zn = 1, C = 1, O = 3))
  )
}

.assay_formula_mass <- function(atoms) {
  sum(atoms * .assay_atomic_weights[names(atoms)])
}

.validate_mineral_library <- function(library) {
  if (!is.list(library) || !length(library) || is.null(names(library)) ||
      !.valid_names(names(library), length(library))) {
    stop("library must be a named list of unique mineral entries.")
  }
  supported <- c(setdiff(names(.assay_atomic_weights), c("Fe", "Mn")),
                 "Fe2", "Fe3", "Mn2", "Mn3", "Mn4")
  elements <- unique(c("O", "H", "S", "As", "Au", "Ag",
                       unname(.assay_oxide_cation)))
  atomic <- matrix(0, nrow = length(library), ncol = length(elements),
                   dimnames = list(names(library), elements))
  masses <- stats::setNames(numeric(length(library)), names(library))
  compatible <- stats::setNames(logical(length(library)), names(library))
  for (i in seq_along(library)) {
    entry <- library[[i]]
    if (!is.list(entry) || !.valid_names(names(entry), length(entry)) ||
        !all(c("formula", "atoms") %in% names(entry)) ||
        !is.character(entry$formula) || !is.null(dim(entry$formula)) ||
        length(entry$formula) != 1L || is.na(entry$formula) ||
        !nzchar(trimws(entry$formula)) ||
        !is.numeric(entry$atoms) || !is.null(dim(entry$atoms)) ||
        !length(entry$atoms) || !.valid_names(names(entry$atoms), length(entry$atoms)) || any(!is.finite(entry$atoms)) ||
        any(entry$atoms <= 0) || !all(names(entry$atoms) %in% supported)) {
      stop("Invalid mineral entry '", names(library)[i],
           "': supply a formula label and positive named supported atom counts; use Fe2/Fe3 for iron.")
    }
    atoms <- entry$atoms
    atomic[i, names(atoms)] <- atoms
    base_names <- ifelse(names(atoms) %in% c("Fe2", "Fe3"), "Fe",
                         ifelse(names(atoms) %in% c("Mn2", "Mn3", "Mn4"),
                                "Mn", names(atoms)))
    masses[i] <- sum(atoms * .assay_atomic_weights[base_names])
    if (!is.finite(masses[i]) || masses[i] <= 0)
      stop("Mineral formula mass must be positive and finite.")
    # The chosen conventional oxides plus H2O must account for every oxygen.
    required_oxygen <- 0
    for (oxide in names(.assay_oxide_cation)) {
      cation <- .assay_oxide_cation[[oxide]]
      required_oxygen <- required_oxygen + atomic[i, cation] /
        .assay_oxide_formula[[oxide]][1L] * .assay_oxide_formula[[oxide]]["O"]
    }
    compatible[i] <- abs(atomic[i, "O"] - required_oxygen) <= 1e-8 &&
      !any(atomic[i, c("S", "As", "Au", "Ag")] > 0)
    if (!compatible[i] && (atomic[i, "O"] > 0 || atomic[i, "H"] > 0 ||
                            atomic[i, "C"] > 0)) {
      stop("Mineral '", names(library)[i],
           "' has unsupported oxygen/volatile stoichiometry. Only oxide/hydroxide/carbonate formulas or oxygen-free sulfide/native formulas are supported.")
    }
  }
  list(atoms = atomic, molar_mass = masses, oxide_compatible = compatible)
}

.assay_oxide_molar_masses <- function() {
  vapply(.assay_oxide_formula, .assay_formula_mass, numeric(1))
}

#' Convert ideal mineral mass fractions to theoretical bulk assays.
#'
#' @param mineral_proportions Named numeric vector or numeric matrix/data frame
#'   with mineral names as columns. Rows must already sum to 1 for `fraction`
#'   basis or 100 for `percent` basis; no implicit renormalization occurs.
#' @param library Named list of entries with `formula` labels and explicit
#'   positive named `atoms` stoichiometry; defaults to `mineral_assay_library()`.
#' @param basis Input mass basis, `fraction` or `percent`.
#' @param iron_reporting `separate` retains physical FeO and Fe2O3. `total_Fe2O3`
#'   reports all Fe as Fe2O3 equivalent; this adds oxygen-equivalent mass and
#'   can produce oxide + H2O totals above 100 percent.
#' @param closure_tolerance Maximum allowed row-sum deviation on input basis.
#' @return A list of elemental wt%, physical oxide wt%, reported oxide wt%,
#'   theoretical H2O-equivalent wt%, Fe(II)/Fe(III) wt%, and mass-balance audit.
#'   No output is renormalized.
minerals_to_assay <- function(mineral_proportions,
                               library = mineral_assay_library(),
                               basis = c("fraction", "percent"),
                               iron_reporting = c("separate", "total_Fe2O3"),
                               closure_tolerance = 1e-6) {
  if (!is.null(dim(mineral_proportions)) &&
      !is.matrix(mineral_proportions) && !is.data.frame(mineral_proportions))
    stop("mineral_proportions must be a vector, matrix or data frame.")
  basis <- match.arg(basis)
  iron_reporting <- match.arg(iron_reporting)
  checked <- .validate_mineral_library(library)
  if (!(is.numeric(mineral_proportions) || is.data.frame(mineral_proportions)) ||
      (is.data.frame(mineral_proportions) &&
       !all(vapply(mineral_proportions, is.numeric, logical(1)))))
    stop("mineral_proportions must be numeric.")
  x <- if (is.numeric(mineral_proportions) && is.null(dim(mineral_proportions))) {
    matrix(mineral_proportions, nrow = 1L,
           dimnames = list(NULL, names(mineral_proportions)))
  } else as.matrix(mineral_proportions)
  if (!is.numeric(x) || length(dim(x)) != 2L ||
      !.valid_names(colnames(x), ncol(x)) ||
      !nrow(x) || !ncol(x) || any(!is.finite(x)) || any(x < 0) ||
      !all(colnames(x) %in% names(library))) {
    stop("mineral_proportions must be a nonnegative numeric vector/matrix with supported mineral names.")
  }
  target <- if (basis == "fraction") 1 else 100
  if (!is.numeric(closure_tolerance) || !is.null(dim(closure_tolerance)) ||
      length(closure_tolerance) != 1L || !is.finite(closure_tolerance) ||
      closure_tolerance < 0 ||
      any(abs(rowSums(x) - target) > closure_tolerance)) {
    stop("Mineral mass proportions do not sum to the declared basis; no renormalization is performed.")
  }
  fractions <- x / target
  atoms <- checked$atoms[colnames(x), , drop = FALSE]
  molar_mass <- checked$molar_mass[colnames(x)]
  mineral_moles <- sweep(fractions * 100, 2L, molar_mass, "/")
  atom_moles <- mineral_moles %*% atoms
  oxide_moles <- do.call(cbind, lapply(names(.assay_oxide_cation), function(oxide) {
    cation <- .assay_oxide_cation[[oxide]]
    atom_moles[, cation] / .assay_oxide_formula[[oxide]][1L]
  }))
  colnames(oxide_moles) <- names(.assay_oxide_cation)
  oxide_mass <- sweep(oxide_moles, 2L, .assay_oxide_molar_masses(), "*")
  volatile_names <- c("H2O", "CO2")
  report_names <- setdiff(colnames(oxide_mass), volatile_names)
  reported <- oxide_mass[, report_names, drop = FALSE]
  compatible_moles <- mineral_moles
  compatible_moles[, !checked$oxide_compatible[colnames(x)]] <- 0
  compatible_atom_moles <- compatible_moles %*% atoms
  compatible_oxide_moles <- do.call(cbind, lapply(names(.assay_oxide_cation),
    function(oxide) {
      cation <- .assay_oxide_cation[[oxide]]
      compatible_atom_moles[, cation] / .assay_oxide_formula[[oxide]][1L]
    }))
  colnames(compatible_oxide_moles) <- names(.assay_oxide_cation)
  compatible_oxide_mass <- sweep(compatible_oxide_moles, 2L,
                                 .assay_oxide_molar_masses(), "*")
  physical <- compatible_oxide_mass[, report_names, drop = FALSE]
  water <- compatible_oxide_mass[, "H2O"]
  carbon_dioxide <- compatible_oxide_mass[, "CO2"]
  nonoxide_mass <- rowSums(fractions[, !checked$oxide_compatible[colnames(x)],
                                     drop = FALSE]) * 100
  if (iron_reporting == "total_Fe2O3") {
    feo_to_fe2o3 <- .assay_oxide_molar_masses()["Fe2O3"] /
      (2 * .assay_oxide_molar_masses()["FeO"])
    reported[, "Fe2O3"] <- reported[, "Fe2O3"] +
      reported[, "FeO"] * feo_to_fe2o3
    reported <- reported[, setdiff(colnames(reported), "FeO"), drop = FALSE]
    colnames(reported)[colnames(reported) == "Fe2O3"] <- "Fe2O3_total"
  }
  base_elements <- ifelse(colnames(atom_moles) %in% c("Fe2", "Fe3"), "Fe",
                          ifelse(colnames(atom_moles) %in%
                            c("Mn2", "Mn3", "Mn4"), "Mn", colnames(atom_moles)))
  element_mass <- sweep(atom_moles, 2L,
                        .assay_atomic_weights[base_elements], "*")
  iron_valence <- cbind(FeII_wt_pct = element_mass[, "Fe2"],
                        FeIII_wt_pct = element_mass[, "Fe3"])
  elements <- element_mass[, !colnames(element_mass) %in%
                             c("Fe2", "Fe3", "Mn2", "Mn3", "Mn4"),
                           drop = FALSE]
  elements <- cbind(elements, Fe = rowSums(iron_valence),
                    Mn = rowSums(element_mass[, c("Mn2", "Mn3", "Mn4"),
                                               drop = FALSE]))
  # Keep only present elements, apart from Fe which is a useful reported field.
  keep <- colSums(elements) > 0
  elements <- elements[, keep, drop = FALSE]
  unreported_names <- intersect(c("S", "As", "Au", "Ag"), colnames(elements))
  unreported <- if (length(unreported_names))
    elements[, unreported_names, drop = FALSE] else
    matrix(0, nrow = nrow(x), ncol = 0L, dimnames = list(rownames(x), NULL))
  reported_total <- rowSums(reported) + water + carbon_dioxide +
    rowSums(unreported)
  audit <- data.frame(
    mineral_mass_wt_pct = rowSums(fractions) * 100,
    element_mass_wt_pct = rowSums(elements),
    physical_oxides_plus_H2O_wt_pct = rowSums(physical) + water,
    CO2_equiv_wt_pct = carbon_dioxide,
    nonoxide_mineral_mass_wt_pct = nonoxide_mass,
    physical_components_wt_pct = rowSums(physical) + water +
      carbon_dioxide + nonoxide_mass,
    reported_components_wt_pct = reported_total,
    reporting_oxygen_gain_wt_pct = reported_total - 100
  )
  if (any(!is.finite(as.matrix(audit))) || any(!is.finite(elements)) ||
      any(!is.finite(reported)) || any(!is.finite(elements * 10000)) ||
      max(abs(audit$element_mass_wt_pct - audit$mineral_mass_wt_pct)) > 1e-6 ||
      max(abs(audit$physical_components_wt_pct -
              audit$mineral_mass_wt_pct)) > 1e-6) {
    stop("Internal mineral-to-assay mass balance failed.")
  }
  rownames(elements) <- rownames(physical) <- rownames(reported) <-
    rownames(iron_valence) <- rownames(audit) <- rownames(x)
  list(elements = as.data.frame(elements),
       element_ppm = as.data.frame(elements * 10000),
       oxides_physical = as.data.frame(physical),
       oxides_reported = as.data.frame(reported),
       water = data.frame(H2O_equiv_wt_pct = water, row.names = rownames(x)),
       volatiles = data.frame(H2O_equiv_wt_pct = water,
                              CO2_equiv_wt_pct = carbon_dioxide,
                              row.names = rownames(x)),
       unreported_elements = as.data.frame(unreported),
       iron_valence = as.data.frame(iron_valence), mass_balance = audit,
       iron_reporting = iron_reporting,
       assumptions = list(ideal_fixed_stoichiometry = TRUE,
                          hydrogen_as_H2O_equivalent = TRUE,
                          oxides_physical_exclude_nonoxide_minerals = TRUE,
                          oxides_reported_are_cation_equivalents = TRUE,
                          no_unmeasured_volatiles = TRUE,
                          input_basis = basis))
}
