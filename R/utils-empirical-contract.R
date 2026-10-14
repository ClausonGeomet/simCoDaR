# Input and serialization contracts for the recovered empirical API.
.empirical_integer <- function(x, name, minimum = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x > .Machine$integer.max - 1 || x != floor(x))
    stop(name, " must be an integer >= ", minimum, " and < .Machine$integer.max.")
  as.integer(x)
}
.empirical_seed <- function(seed) {
  if (is.null(seed)) NULL else .empirical_integer(seed, "seed", 0)
}
.empirical_names <- function(x, name, minimum = 1L) {
  if (!is.character(x) || length(x) < minimum || anyNA(x) ||
      any(!nzchar(trimws(x))) || anyDuplicated(x))
    stop(name, " must contain distinct nonempty names.")
  x
}
.empirical_control <- function(control) {
  defaults <- list(dirichlet_maxit = 2000L, student_max_iter = 500L,
                   mixture_max_iter = 300L, mixture_starts = 5L, tolerance = 1e-7)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% names(defaults))))) stop("Unknown or unnamed control override.")
  for (key in names(control)) {
    value <- control[[key]]
    if (key == "tolerance") {
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value <= 0 || value > 1e-3) stop("tolerance must be in (0, 1e-3].")
      defaults[[key]] <- value
    } else defaults[[key]] <- .empirical_integer(value, key)
  }
  defaults
}
.empirical_close <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || any(!is.finite(x)) || any(x < 0))
    stop("Assay parts must be complete, finite numeric nonnegative values; censoring requires separate treatment.")
  maximum <- apply(x, 1L, max)
  if (any(maximum <= 0)) stop("Every row must have positive selected-part total.")
  scaled <- x / maximum
  totals <- rowSums(scaled)
  # Reject unrepresentable selected-part totals; closure itself remains overflow-safe.
  if (any(maximum > .Machine$double.xmax / totals)) stop("Selected-part row total is nonfinite.")
  out <- scaled / totals
  if (any(out[x > 0] <= 0)) stop("Positive assay parts underflowed during closure.")
  out
}
.empirical_units <- function(parts, unit_label, part_units, measurement_basis) {
  basis <- match.arg(measurement_basis, c("unspecified", "oxide_mass", "element_mass",
    "mineral_mass", "mineral_volume", "abstract_composition"))
  if (!is.null(unit_label)) .empirical_names(unit_label, "unit_label")
  if (!is.null(unit_label) && length(unit_label) != 1L) stop("unit_label must be one string.")
  if (!is.null(part_units)) {
    if (!is.character(part_units) || length(part_units) != length(parts) ||
        anyNA(part_units) || any(!nzchar(trimws(part_units)))) stop("part_units must give one unit per selected part.")
    .empirical_names(names(part_units), "part_units names")
    if (!setequal(names(part_units), parts)) stop("part_units must name every selected part exactly once.")
    part_units <- part_units[parts]
    if (length(unique(part_units)) != 1L) stop("Mixed part units are unsupported; no silent unit conversion.")
    if (!is.null(unit_label) && unit_label != part_units[1L]) stop("unit_label and part_units disagree.")
    unit_label <- unname(part_units[1L])
  }
  # Standard converter suffixes and total-equivalent iron labels.
  chemistry <- sub("_(wt_pct|ppm|fraction|vol_pct)$", "", parts)
  total <- grepl("(^total_Fe|^Fe_total|^Fe2O3_total|^FeO_total|^FeOT$|^Fe2O3T$)", chemistry)
  stripped <- sub("_total$", "", chemistry)
  atoms <- .assay_oxide_formula
  for (element in names(.assay_atomic_weights)) {
    oxides <- names(atoms)[vapply(atoms, function(a) element %in% names(a), logical(1))]
    # O itself is not an independently nonoverlapping oxide component.
    if (element %in% stripped && any(stripped %in% oxides))
      stop("Overlapping elemental/oxide reporting for ", element, "; resolve the measurement basis.")
  }
  if ("Fe" %in% stripped && any(stripped %in% c("FeII", "FeIII")))
    stop("Overlapping total elemental iron and separate valence reporting.")
  if (any(total) && (sum(total) > 1L || any(!total & stripped %in% c("Fe", "FeO", "Fe2O3", "FeII", "FeIII"))))
    stop("Overlapping total iron and separate iron reporting.")
  if (anyDuplicated(stripped)) stop("Duplicate chemical species in selected reporting columns.")
  suffix_units <- ifelse(grepl("_wt_pct$", parts), "wt%",
    ifelse(grepl("_ppm$", parts), "ppm", ifelse(grepl("_vol_pct$", parts), "vol%",
      ifelse(grepl("_fraction$", parts), "fraction", NA_character_))))
  detected <- unique(suffix_units[!is.na(suffix_units)])
  if (length(detected) > 1L) stop("Mixed units in recognized reporting suffixes.")
  canonical <- function(u) {
    if (is.null(u)) return(NULL)
    u <- tolower(gsub("[ _]", "", u))
    if (u %in% c("wt%", "wtpct", "weightpercent", "mass%")) "wt%" else
      if (u %in% c("vol%", "volpct", "volumepercent")) "vol%" else u
  }
  declared <- canonical(unit_label)
  if (length(detected) && !is.null(declared) && detected != declared)
    stop("Declared unit disagrees with recognized column suffix.")
  effective <- if (is.null(declared) && length(detected)) detected else declared
  if (basis == "mineral_volume" && !is.null(effective) && effective %in% c("wt%", "ppm"))
    stop("Mass units cannot declare mineral_volume.")
  if (basis %in% c("oxide_mass", "element_mass", "mineral_mass") &&
      identical(effective, "vol%")) stop("Volume units cannot declare a mass basis.")
  if (basis == "oxide_mass" && any(stripped %in% setdiff(names(.assay_atomic_weights), "O")))
    stop("Recognized elemental column conflicts with oxide_mass.")
  if (basis == "element_mass" && any(stripped %in% names(atoms)))
    stop("Recognized oxide column conflicts with element_mass.")
  list(measurement_basis = basis, unit_label = unit_label, part_units = part_units,
       recognized_suffix_unit = if (length(detected)) detected else NULL,
       verified = FALSE, conversion = "none")
}
.empirical_labels <- function(x) {
  if (!(is.character(x) || is.numeric(x) || is.factor(x)) || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x)))) stop("Groups must be complete finite numeric or nonempty character labels.")
  x <- as.character(x)
  if (any(!nzchar(trimws(x)))) stop("Group labels must be nonempty.")
  x
}
.assay_pattern_key <- function(positive) paste0(as.integer(positive), collapse = "")

.empirical_probabilities <- function(x, labels, name) {
  if (!is.numeric(x) || length(x) != length(labels) || any(!is.finite(x)) ||
      any(x <= 0) || abs(sum(x) - 1) > 1e-10 || !identical(names(x), labels))
    stop("Invalid fitted ", name, ".")
}
.empirical_spd <- function(s, p) {
  if (!is.matrix(s) || !is.numeric(s) || !identical(dim(s), c(p, p)) ||
      any(!is.finite(s)) || max(abs(s - t(s))) > 1e-10 * max(1, max(abs(s))))
    stop("Invalid fitted covariance/scale.")
  tryCatch(chol(s), error = function(e) stop("Invalid fitted covariance/scale: not positive definite."))
  invisible(TRUE)
}
.validate_empirical_fit <- function(fit) {
  if (!inherits(fit, "geo_assay_simulator") || !is.list(fit) ||
      !identical(fit$schema_version, 1L)) stop("Invalid geo_assay_simulator schema.")
  parts <- .empirical_names(fit$parts, "fitted parts", 2L)
  models <- c("logistic_normal", "ilr_student_t", "ilr_mixture", "dirichlet")
  if (length(fit$model) != 1L || !fit$model %in% models) stop("Invalid fitted model.")
  .empirical_seed(fit$seed)
  .empirical_control(fit$control)
  .empirical_integer(fit$fitted_rows, "fitted_rows")
  minimum <- .empirical_integer(fit$min_pattern_n, "fitted min_pattern_n", 2L)
  if (!is.list(fit$units)) stop("Invalid fitted unit metadata.")
  .empirical_units(parts, fit$units$unit_label, fit$units$part_units,
                   fit$units$measurement_basis)
  if (!is.numeric(fit$student_df) || length(fit$student_df) != 1L ||
      !is.finite(fit$student_df) || fit$student_df <= 2) stop("Invalid fitted student_df.")
  ks <- fit$mixture_components
  if (!is.integer(ks) || !length(ks) || anyNA(ks) || any(ks < 1L) ||
      !identical(ks, sort(unique(ks)))) stop("Invalid fitted mixture candidates.")
  if (!fit$zero_mode %in% c("error", "structural")) stop("Invalid fitted zero policy.")
  if (!is.null(fit$structural_parts)) {
    .empirical_names(fit$structural_parts, "fitted structural_parts")
    if (!all(fit$structural_parts %in% parts)) stop("Invalid structural parts.")
  }
  if (!is.null(fit$group_col)) {
    .empirical_names(fit$group_col, "fitted group_col")
    if (length(fit$group_col) != 1L || fit$group_col %in% parts) stop("Invalid fitted group column.")
  }
  if (!is.list(fit$groups) || !length(fit$groups)) stop("Invalid fitted groups.")
  labels <- .empirical_names(names(fit$groups), "fitted group names")
  if (is.null(fit$group_col) && !identical(labels, ".all")) stop("Invalid ungrouped fit.")
  .empirical_probabilities(fit$group_probabilities, labels, "group probabilities")
  counts <- numeric(length(labels))
  for (i in seq_along(labels)) {
    g <- fit$groups[[i]]
    counts[i] <- .empirical_integer(g$count, "group count")
    if (!is.list(g$patterns) || !length(g$patterns)) stop("Invalid fitted patterns.")
    keys <- .empirical_names(names(g$patterns), "support names")
    .empirical_probabilities(g$pattern_probabilities, keys, "support probabilities")
    pattern_counts <- numeric(length(keys))
    for (j in seq_along(keys)) {
      h <- g$patterns[[j]]
      pattern_counts[j] <- .empirical_integer(h$count, "pattern count")
      active <- .empirical_names(h$active_parts, "active parts")
      if (!identical(active, parts[parts %in% active]) ||
          keys[j] != .assay_pattern_key(parts %in% active)) stop("Invalid ordered support.")
      absent <- setdiff(parts, active)
      if (length(absent) && (fit$zero_mode != "structural" ||
          !all(absent %in% fit$structural_parts))) stop("Support violates fitted zero policy.")
      p <- as.integer(length(active) - 1L)
      if (p == 0L) {
        if (!identical(h$kind, "singleton") || !is.null(h$model)) stop("Invalid singleton fit.")
        next
      }
      if (h$count < max(minimum, p + 2L)) stop("Invalid sparse fitted support.")
      if (!identical(h$kind, fit$model)) stop("Inconsistent pattern model.")
      b <- h$basis
      if (!is.matrix(b) || !identical(dim(b), c(p, p + 1L)) || any(!is.finite(b)) ||
          !identical(colnames(b), active) || max(abs(rowSums(b))) > 1e-10 ||
          max(abs(tcrossprod(b) - diag(p))) > 1e-10) stop("Invalid fitted ILR basis.")
      m <- h$model
      if (!is.list(m) || !identical(m$convergence, 0L) ||
          length(m$loglik) != 1L || !is.finite(m$loglik)) stop("Invalid unconverged model parameters.")
      location <- function(mu) {
        if (!is.numeric(mu) || length(mu) != p || any(!is.finite(mu))) stop("Invalid fitted location.")
      }
      if (h$kind == "dirichlet") {
        if (!is.numeric(m$alpha) || length(m$alpha) != p + 1L ||
            any(!is.finite(m$alpha)) || any(m$alpha <= 0) || !identical(names(m$alpha), active) ||
            any(abs(log(m$alpha)) >= 16 - 1e-5)) stop("Invalid fitted alpha.")
      } else if (h$kind == "ilr_mixture") {
        k <- length(m$weights)
        if (!k || !is.numeric(m$weights) || any(!is.finite(m$weights)) ||
            any(m$weights <= 0) || abs(sum(m$weights) - 1) > 1e-10 ||
            !is.matrix(m$means) || !identical(dim(m$means), c(as.integer(k), p)) ||
            any(!is.finite(m$means)) || !is.list(m$covariances) || length(m$covariances) != k)
          stop("Invalid fitted mixture parameters.")
        if (!identical(as.integer(k), m$selected_components) || !k %in% ks ||
            h$count < max(20, 3 * k * (p + 1))) stop("Invalid fitted mixture selection.")
        lapply(m$covariances, .empirical_spd, p = p)
      } else {
        location(m$mu)
        .empirical_spd(m$sigma, p)
        if (h$kind == "ilr_student_t" && (!is.numeric(m$df) || length(m$df) != 1L ||
            !is.finite(m$df) || m$df <= 2 || !identical(m$df, fit$student_df))) stop("Invalid fitted t df.")
        if (h$kind == "ilr_student_t") {
          .empirical_spd(m$covariance, p)
          if (max(abs(m$covariance - m$sigma * m$df / (m$df - 2))) >
              1e-10 * max(1, max(abs(m$covariance)))) stop("Inconsistent t scale/covariance.")
        }
      }
    }
    if (sum(pattern_counts) != g$count ||
        max(abs(pattern_counts / g$count - g$pattern_probabilities)) > 1e-10)
      stop("Inconsistent support counts/probabilities.")
  }
  if (sum(counts) != fit$fitted_rows ||
      max(abs(counts / fit$fitted_rows - fit$group_probabilities)) > 1e-10)
    stop("Inconsistent group counts/probabilities.")
  invisible(TRUE)
}
