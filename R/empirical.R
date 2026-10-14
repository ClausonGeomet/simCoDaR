# Adapted from GPL-3 CoDA R/simulator.R and R/composition_study.R.
# Immutable source mapping, origin attribution and changes: inst/PROVENANCE.md.

.fit_assay_pattern <- function(x, model, mixture_components, student_df, control) {
  if (ncol(x) == 1L) return(list(kind = "singleton", model = NULL,
    diagnostics = list(status = "deterministic", measure = "point_mass")))
  basis <- .ilr_basis(ncol(x))
  colnames(basis) <- colnames(x)
  rownames(basis) <- paste0("ilr", seq_len(nrow(basis)))
  z <- .to_ilr(x, basis)
  if (any(!is.finite(z)) || sum(apply(z, 2L, stats::var)) <= 0)
    stop("Multivariate support is unidentifiable: no logratio dispersion.")
  if (model == "dirichlet") {
    m <- .fit_dirichlet(x, control$dirichlet_maxit)
    # Store common ILR-volume likelihood alongside original simplex likelihood.
    m$ilr_loglik <- m$loglik - sum(-log(ncol(x)) / 2 - rowSums(log(x)))
  } else if (model == "logistic_normal") m <- .fit_logistic_normal(z) else
    if (model == "ilr_student_t") m <- .fit_student_t(z, student_df,
      control$student_max_iter, control$tolerance) else {
      candidates <- lapply(mixture_components, function(k) tryCatch(
        .fit_gaussian_mixture(z, k, control$mixture_starts, control$mixture_max_iter,
          control$tolerance), error = function(e) list(failure = conditionMessage(e))))
      names(candidates) <- as.character(mixture_components)
      bic <- vapply(seq_along(candidates), function(i) {
        f <- candidates[[i]]
        if (!is.null(f$failure)) return(Inf)
        k <- mixture_components[i]
        p <- ncol(z)
        -2 * f$loglik + ((k - 1) + k * p + k * p * (p + 1) / 2) * log(nrow(z))
      }, numeric(1))
      names(bic) <- names(candidates)
      if (!any(is.finite(bic))) stop("No candidate mixture converged: ",
        paste(vapply(candidates, function(f) f$failure, character(1)), collapse = "; "))
      chosen <- which.min(bic)
      m <- candidates[[chosen]]
      m$selected_components <- mixture_components[chosen]
      m$bic <- bic
      m$diagnostics$candidates <- lapply(seq_along(candidates), function(i) list(
        k = mixture_components[i], bic = bic[i], failure = candidates[[i]]$failure,
        diagnostics = candidates[[i]]$diagnostics))
      m$diagnostics$selection <- "minimum_training_BIC_smallest_k_tie"
    }
  list(kind = model, model = m, basis = basis,
       diagnostics = list(measure = "ilr_volume_on_active_support", fit = m$diagnostics))
}

#' Fit an empirical simulator for selected assay subcompositions
#'
#' Fit four recovered compositional families on a caller-confirmed common,
#' nonoverlapping measurement basis. Selected amounts are closed to one;
#' no whole-rock mass balance, total model or censoring replacement is inferred.
#' A supplied seed is scoped with withr and restores the caller RNG on errors.
#' The serializable object retains fitted parameters and aggregate counts, not
#' observed rows or identifiers. Necessary group labels are retained.
#' @param data Complete data frame; selected numeric columns must be finite.
#' @param parts Distinct selected column names in desired output order.
#' @param model One of logistic_normal, ilr_student_t, ilr_mixture, dirichlet.
#' @param group_col Optional separate complete geological group column.
#' @param seed Integer seed or NULL to consume the caller RNG.
#' @param zero_mode Reject zeros (error) or explicitly declared structural zeros.
#' @param structural_parts Selected names allowed to be structurally absent.
#' @param min_pattern_n Minimum rows for multivariate support; also needs p+2.
#' @param mixture_components Positive integer candidate counts; training BIC only.
#' @param student_df Fixed predeclared t degrees of freedom greater than two.
#' @param unit_label Descriptive common unit; never performs a conversion.
#' @param part_units Optional named common unit for every selected part.
#' @param measurement_basis Unspecified, oxide_mass, element_mass, mineral_mass,
#'   mineral_volume, or abstract_composition.
#' @param control Named iteration/tolerance overrides; see Details.
#' @details Controls: dirichlet_maxit=2000, student_max_iter=500,
#'   mixture_max_iter=300, mixture_starts=5, tolerance=1e-7. Sparse or identical
#'   multivariate supports and nonconvergence fail without row-resampling fallback.
#'   Gaussian covariance uses n-1 plus ridge. Student sigma is scale, with
#'   covariance df/(df-2)*sigma. Mixture uses full covariance and multiple starts.
#'   Arbitrary names cannot prove absence of overlapping chemical species;
#'   callers must resolve unidentified overlaps. Group labels are not anonymized.
#' @return A geo_assay_simulator list with ordered supports, parameters,
#'   basis, diagnostics, controls, RNG metadata and assumptions.
#' @export
fit_assay_simulator <- function(data, parts,
    model = c("logistic_normal", "ilr_student_t", "ilr_mixture", "dirichlet"),
    group_col = NULL, seed = 42, zero_mode = c("error", "structural"),
    structural_parts = NULL, min_pattern_n = 30L, mixture_components = 1:3,
    student_df = 5, unit_label = NULL, part_units = NULL,
    measurement_basis = "unspecified", control = list()) {
  model <- match.arg(model)
  zero_mode <- match.arg(zero_mode)
  seed <- .empirical_seed(seed)
  control <- .empirical_control(control)
  if (!is.data.frame(data) || !nrow(data))
    stop("data must be a nonempty data frame with unique column names.")
  .empirical_names(names(data), "data column names")
  parts <- .empirical_names(parts, "parts", 2L)
  if (!all(parts %in% names(data)) || !all(vapply(data[parts], is.numeric, logical(1))))
    stop("parts must name existing numeric columns.")
  min_pattern_n <- .empirical_integer(min_pattern_n, "min_pattern_n", 2L)
  if (!is.numeric(mixture_components) || !length(mixture_components) ||
      any(!is.finite(mixture_components)) || any(mixture_components < 1) ||
      any(mixture_components >= .Machine$integer.max) ||
      any(mixture_components != floor(mixture_components)))
    stop("mixture_components must contain positive integers without truncation.")
  mixture_components <- sort(unique(as.integer(mixture_components)))
  if (!is.numeric(student_df) || length(student_df) != 1L || !is.finite(student_df) || student_df <= 2)
    stop("student_df must be a finite scalar greater than two.")
  units <- .empirical_units(parts, unit_label, part_units, measurement_basis)
  x <- .empirical_close(as.matrix(data[parts]))
  if (!is.null(structural_parts)) {
    structural_parts <- .empirical_names(structural_parts, "structural_parts")
    if (!all(structural_parts %in% parts)) stop("structural_parts must name selected parts.")
  }
  if (zero_mode == "error" && any(x == 0))
    stop("Zero assay parts found: declare structural absence explicitly; censoring/rounding needs separate treatment.")
  if (any(x == 0) && (is.null(structural_parts) ||
      any(x[, setdiff(parts, structural_parts), drop = FALSE] == 0)))
    stop("Zeros occur outside declared structural_parts.")
  if (is.null(group_col)) labels <- rep(".all", nrow(x)) else {
    .empirical_names(group_col, "group_col")
    if (length(group_col) != 1L || !group_col %in% names(data) || group_col %in% parts)
      stop("group_col must name one separate complete column.")
    labels <- .empirical_labels(data[[group_col]])
  }
  fit_core <- function() {
    # First-occurrence order is explicit and stable, including numeric labels.
    group_names <- unique(labels)
    groups <- lapply(group_names, function(label) {
      subset <- x[labels == label, , drop = FALSE]
      keys <- apply(subset > 0, 1L, .assay_pattern_key)
      support_names <- unique(keys)
      patterns <- lapply(support_names, function(key) {
        observed <- subset[keys == key, , drop = FALSE]
        active <- parts[observed[1L, ] > 0]
        observed <- observed[, active, drop = FALSE]
        count <- nrow(observed)
        contextual <- function(e) stop("Group '", label, "', support ", key,
          " (n=", count, "): ", conditionMessage(e), call. = FALSE)
        tryCatch({
          if (length(active) > 1L && count < max(min_pattern_n, length(active) + 1L))
            stop("Unsupported sparse multivariate support; no empirical fallback.")
          fitted <- .fit_assay_pattern(observed, model, mixture_components, student_df, control)
          c(fitted, list(count = count, active_parts = active))
        }, error = contextual)
      })
      names(patterns) <- support_names
      list(count = nrow(subset), patterns = patterns,
        pattern_probabilities = stats::setNames(vapply(patterns, function(h) h$count,
          numeric(1)) / nrow(subset), support_names))
    })
    names(groups) <- group_names
    structure(list(schema_version = 1L, parts = parts, model = model,
      group_col = group_col, groups = groups,
      group_probabilities = stats::setNames(vapply(groups, function(g) g$count,
        numeric(1)) / nrow(x), group_names), fitted_rows = nrow(x), seed = seed,
      rng = list(kind = RNGkind(), policy = if (is.null(seed)) "consume_caller" else "scoped_restore",
        R_version = as.character(getRversion()),
        package_version = as.character(getNamespaceVersion("simCoDaR")),
        withr_version = as.character(getNamespaceVersion(asNamespace("withr")))),
      units = units, zero_mode = zero_mode, structural_parts = structural_parts,
      min_pattern_n = min_pattern_n, student_df = student_df,
      mixture_components = mixture_components, control = control,
      source_summary = list(rows = nrow(x), selected_parts = length(parts),
        groups = length(groups), totals_retained = FALSE, observed_rows_retained = FALSE),
      assumptions = list(selected_part_subcomposition = TRUE, whole_rock_mass_closure = FALSE,
        total_composition_dependence_modelled = FALSE, common_nonoverlapping_basis_confirmed_by_caller = TRUE,
        chemical_overlap_detection_limited_to_recognized_names = TRUE,
        units_converted = FALSE, censoring_modelled = FALSE,
        group_labels_retained = !is.null(group_col), privacy_guarantee = "none",
        density_measure = "ilr_volume_on_active_support_plus_discrete_support_and_group")),
      class = "geo_assay_simulator")
  }
  if (is.null(seed)) fit_core() else withr::with_seed(seed, fit_core())
}

.draw_assay_pattern <- function(pattern, n) {
  if (pattern$kind == "singleton") return(list(x = matrix(1, n, 1L), floored_rows = 0L, floored_cells = 0L))
  m <- pattern$model
  if (pattern$kind == "dirichlet") return(list(x = .dirichlet_draw(n, m), floored_rows = 0L, floored_cells = 0L))
  z <- switch(pattern$kind,
    logistic_normal = .mvnorm_draw(n, m$mu, m$sigma),
    ilr_student_t = .mvt_draw(n, m), ilr_mixture = .mixture_draw(n, m))
  if (any(!is.finite(z))) stop("Nonfinite ILR draw.")
  logx <- z %*% pattern$basis
  if (any(!is.finite(logx))) stop("Nonfinite inverse ILR draw.")
  shifted <- sweep(logx, 1L, apply(logx, 1L, max), "-")
  if (any(!is.finite(shifted))) stop("Nonfinite stabilized inverse ILR draw.")
  clipped <- shifted < -700
  x <- exp(pmax(shifted, -700))
  list(x = x / rowSums(x), floored_rows = sum(rowSums(clipped) > 0), floored_cells = sum(clipped))
}

#' Draw synthetic selected-part proportions from a fitted assay simulator
#' @param fit A valid object from fit_assay_simulator.
#' @param n Positive integer number of rows.
#' @param seed Integer seed or NULL to consume the caller RNG.
#' @param group Optional fitted group label, length one or n.
#' @param output Only proportion is supported; total generation is not included.
#' @return Data frame in caller part order, with optional group column,
#'   synthetic row names and assumptions/diagnostics attributes. Diagnostics
#'   count inverse ILR floor events, which modify extreme tails.
#' @export
simulate_assays <- function(fit, n, seed = 123, group = NULL, output = "proportion") {
  .validate_empirical_fit(fit)
  n <- .empirical_integer(n, "n")
  seed <- .empirical_seed(seed)
  if (!identical(output, "proportion")) stop("Only output='proportion' is supported; total generation is not included.")
  if (!is.null(group)) {
    if (is.null(fit$group_col)) stop("group requires a grouped fit.")
    if (!length(group) %in% c(1L, n)) stop("group must have length one or n.")
    group <- .empirical_labels(group)
    if (!all(group %in% names(fit$groups))) stop("Unknown fitted group label.")
  }
  generate <- function() {
    available <- names(fit$groups)
    labels <- if (!is.null(group)) rep(group, length.out = n) else
      if (is.null(fit$group_col)) rep(".all", n) else
        available[sample.int(length(available), n, replace = TRUE, prob = fit$group_probabilities)]
    out <- matrix(0, n, length(fit$parts), dimnames = list(NULL, fit$parts))
    floor_rows <- floor_cells <- 0L
    for (label in unique(labels)) {
      rows <- which(labels == label)
      g <- fit$groups[[label]]
      keys <- names(g$patterns)
      chosen <- keys[sample.int(length(keys), length(rows), replace = TRUE, prob = g$pattern_probabilities)]
      for (key in unique(chosen)) {
        where <- rows[chosen == key]
        h <- g$patterns[[key]]
        draw <- .draw_assay_pattern(h, length(where))
        if (any(!is.finite(draw$x)) || any(draw$x <= 0)) stop("Invalid active-support draw.")
        out[where, h$active_parts] <- draw$x
        floor_rows <- floor_rows + draw$floored_rows
        floor_cells <- floor_cells + draw$floored_cells
      }
    }
    if (any(!is.finite(out)) || any(out < 0) || max(abs(rowSums(out) - 1)) > 1e-10)
      stop("Simulation produced invalid compositions.")
    result <- as.data.frame(out, check.names = FALSE)
    if (!is.null(fit$group_col)) result[[fit$group_col]] <- labels
    attr(result, "assumptions") <- c(fit$assumptions, list(output = "proportion", unit = "closed_fraction"))
    attr(result, "diagnostics") <- list(inverse_floor_rows = floor_rows, inverse_floor_cells = floor_cells,
      invalid_draw_count = 0L, seed = seed, RNGkind = RNGkind(),
      rng_policy = if (is.null(seed)) "consume_caller" else "scoped_restore")
    result
  }
  if (is.null(seed)) generate() else withr::with_seed(seed, generate())
}
