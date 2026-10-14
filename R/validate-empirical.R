# New bounded diagnostics built on the recovered density and sampler.
.validation_amounts <- function(data, parts, label) {
  if (!is.data.frame(data) || nrow(data) < 2L)
    stop(label, " must be a data frame with at least two rows.")
  .empirical_names(names(data), paste(label, "column names"))
  if (!all(parts %in% names(data)) || !all(vapply(data[parts], is.numeric, logical(1))))
    stop(label, " must contain the fitted numeric parts.")
  x <- .empirical_close(as.matrix(data[parts]))
  if (any(x <= 0)) stop(label, " must be fully positive; cross-support validation is not specified.")
  x
}

.validation_logratios <- function(x, pairs) {
  log(x[, pairs[1L, ], drop = FALSE]) - log(x[, pairs[2L, ], drop = FALSE])
}

.validation_distances <- function(a, b) {
  # Direct differences avoid cancellation for nearly identical ILR vectors.
  out <- matrix(0, nrow(a), nrow(b))
  for (j in seq_len(ncol(a))) out <- out + outer(a[, j], b[, j], "-")^2
  sqrt(out)
}

.validation_metrics <- function(x, z, y, basis, pairs, reference_bounds,
                                observed_distance, diagnostics) {
  w <- .to_ilr(y, basis)
  rows <- list()
  add <- function(metric, observed, generated, discrepancy = NA_real_,
                  part_a = NA_character_, part_b = NA_character_, probability = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(metric = metric, part_a = part_a,
      part_b = part_b, probability = probability, observed = observed,
      generated = generated, discrepancy = discrepancy, stringsAsFactors = FALSE)
  }
  energy <- 2 * mean(.validation_distances(z, w)) - observed_distance -
    mean(.validation_distances(w, w))
  add("ilr_energy", 0, energy)
  cx <- stats::cov(z)
  cy <- stats::cov(w)
  add("ilr_covariance_relative_error", 0,
      sqrt(sum((cx - cy)^2)) / max(sqrt(sum(cx^2)), 1e-10))
  lx <- .validation_logratios(x, pairs)
  ly <- .validation_logratios(y, pairs)
  probs <- c(.05, .25, .5, .75, .95)
  quantile_error <- numeric(ncol(lx))
  for (j in seq_len(ncol(lx))) {
    qx <- stats::quantile(lx[, j], probs, type = 7, names = FALSE)
    qy <- stats::quantile(ly[, j], probs, type = 7, names = FALSE)
    a <- colnames(x)[pairs[1L, j]]
    b <- colnames(x)[pairs[2L, j]]
    for (q in seq_along(probs)) add("logratio_quantile", qx[q], qy[q],
      abs(qx[q] - qy[q]), a, b, probs[q])
    quantile_error[j] <- mean(abs(qx - qy))
    add("logratio_lower_tail", mean(lx[, j] < qx[1L]), mean(ly[, j] < qx[1L]),
      abs(mean(lx[, j] < qx[1L]) - mean(ly[, j] < qx[1L])), a, b, .05)
    add("logratio_upper_tail", mean(lx[, j] > qx[5L]), mean(ly[, j] > qx[5L]),
      abs(mean(lx[, j] > qx[5L]) - mean(ly[, j] > qx[5L])), a, b, .95)
  }
  add("mean_absolute_logratio_quantile_error", 0, mean(quantile_error))
  correlation <- function(v) {
    if (any(vapply(as.data.frame(v), stats::sd, numeric(1)) == 0)) return(NA_real_)
    stats::cor(v[, 1L], v[, 2L])
  }
  if (all(c("Fe2O3", "MgO") %in% colnames(x))) {
    a <- correlation(x[, c("Fe2O3", "MgO"), drop = FALSE])
    b <- correlation(y[, c("Fe2O3", "MgO"), drop = FALSE])
    add("closed_part_correlation", a, b, abs(a - b), "Fe2O3", "MgO")
  }
  outside <- function(v) {
    mean(sweep(v, 2L, reference_bounds[1L, ], "<") |
         sweep(v, 2L, reference_bounds[2L, ], ">"))
  }
  if (is.null(reference_bounds)) add("reference_outside_cell_fraction", NA_real_, NA_real_) else {
    a <- outside(x)
    b <- outside(y)
    add("reference_outside_cell_fraction", a, b, abs(a - b))
  }
  add("maximum_closure_residual", max(abs(rowSums(x) - 1)), max(abs(rowSums(y) - 1)))
  add("nonfinite_cell_count", sum(!is.finite(x)), sum(!is.finite(y)))
  add("negative_cell_count", sum(x < 0), sum(y < 0))
  add("zero_cell_count", sum(x == 0), sum(y == 0))
  for (key in c("inverse_floor_rows", "inverse_floor_cells", "invalid_draw_count"))
    add(key, NA_real_, diagnostics[[key]])
  do.call(rbind, rows)
}

.validation_mc_summary <- function(metrics, template, replicates) {
  out <- list()
  for (i in seq_len(nrow(template))) {
    key <- template[i, c("metric", "part_a", "part_b", "probability")]
    same <- metrics$metric == key$metric
    for (name in c("part_a", "part_b", "probability"))
      same <- same & if (is.na(key[[name]])) is.na(metrics[[name]]) else
        !is.na(metrics[[name]]) & metrics[[name]] == key[[name]]
    for (statistic in c("generated", "discrepancy")) {
      # Retain specified discrepancies even if their inputs are unavailable.
      discrepancy_metrics <- c("logratio_quantile", "logratio_lower_tail",
        "logratio_upper_tail", "closed_part_correlation", "reference_outside_cell_fraction")
      if (statistic == "discrepancy" && !key$metric %in% discrepancy_metrics) next
      values <- metrics[[statistic]][same]
      values <- values[is.finite(values)]
      n <- length(values)
      s <- if (n >= 2L) stats::sd(values) else NA_real_
      out[[length(out) + 1L]] <- cbind(key, data.frame(statistic = statistic,
        n_success = n, n_failure = replicates - n,
        mean = if (n) mean(values) else NA_real_, sd = s,
        mcse = if (n >= 2L) s / sqrt(n) else NA_real_,
        incomplete = n < replicates, stringsAsFactors = FALSE))
    }
  }
  do.call(rbind, out)
}

#' Diagnose a fitted simulator against heldout positive assays
#'
#' Scores an already fitted, ungrouped full-support simulator in a common ILR
#' measure and compares repeated synthetic draws with heldout assays. This
#' function never fits, selects or tunes a model and makes no validation claim.
#' @param fit A valid ungrouped full-support geo_assay_simulator.
#' @param heldout Complete positive same-basis assays with at least two rows.
#' @param reference Optional positive training assays for reference tail bounds.
#' @param n_draws Integer number of generated rows per repetition, at least two.
#' @param replicates Integer number of repetitions, at least two.
#' @param seed Nonnegative integer base seed; repetition r uses seed+r.
#' @return A plain list containing scores, score_summary, replicate_metrics,
#'   mc_summary, failures and metadata. No assay tables or identifiers are kept.
#' @details Scores use fitted ordered Helmert ILR Lebesgue volume, including the
#'   Dirichlet Jacobian. Energy is the ILR V-statistic including diagonal zeros.
#'   Logratio quantiles use type 7; covariance uses denominator n-1. MCSE is
#'   sd/sqrt(successful repetitions), unavailable below two successes. Failures
#'   are explicit and summaries with missing repetitions are incomplete.
#'   Reference tails count cells outside per-part closed .01/.99 bounds.
#'   Supplied seeds restore caller RNG even on error. Units and nonoverlap are
#'   caller responsibilities. Draw MCSE is conditional simulation noise, not
#'   field uncertainty. Group/support weighting requires a separate protocol.
#' @export
validate_assay_simulator <- function(fit, heldout, reference = NULL,
    n_draws = 71L, replicates = 20L, seed = 20261014L) {
  .validate_empirical_fit(fit)
  patterns <- fit$groups[[1L]]$patterns
  if (!is.null(fit$group_col) || length(patterns) != 1L ||
      !identical(patterns[[1L]]$active_parts, fit$parts) ||
      patterns[[1L]]$kind == "singleton")
    stop("Validation requires an ungrouped full positive support; cross-support/group weighting is not specified.")
  n_draws <- .empirical_integer(n_draws, "n_draws", 2L)
  replicates <- .empirical_integer(replicates, "replicates", 2L)
  if (is.null(seed)) stop("seed must be a known nonnegative integer, not NULL.")
  seed <- .empirical_seed(seed)
  if (as.double(seed) + replicates > .Machine$integer.max - 1)
    stop("seed + replicates must be < .Machine$integer.max.")
  x <- .validation_amounts(heldout, fit$parts, "heldout")
  ref <- if (is.null(reference)) NULL else .validation_amounts(reference, fit$parts, "reference")
  basis <- patterns[[1L]]$basis
  z <- .to_ilr(x, basis)
  scores <- .assay_log_density(fit, heldout, measure = "ilr")
  if (any(!is.finite(scores))) stop("Nonfinite heldout score.")
  bounds <- if (is.null(ref)) NULL else apply(ref, 2L, stats::quantile,
    probs = c(.01, .99), type = 7, names = FALSE)
  pairs <- t(which(upper.tri(matrix(FALSE, length(fit$parts), length(fit$parts))),
    arr.ind = TRUE))
  distance <- mean(.validation_distances(z, z))
  empty_diagnostics <- list(inverse_floor_rows = 0L, inverse_floor_cells = 0L, invalid_draw_count = 0L)
  template <- .validation_metrics(x, z, x, basis, pairs, bounds, distance, empty_diagnostics)
  metrics <- cbind(data.frame(replicate = integer(), seed = integer()), template[FALSE, ])
  failures <- data.frame(replicate = integer(), seed = integer(), message = character())
  for (r in seq_len(replicates)) {
    s <- seed + r
    result <- tryCatch(withr::with_seed(s, {
      draw <- simulate_assays(fit, n_draws, seed = s)
      y <- as.matrix(draw[fit$parts])
      .validation_metrics(x, z, y, basis, pairs, bounds, distance, attr(draw, "diagnostics"))
    }), error = function(e) e)
    if (inherits(result, "error")) failures <- rbind(failures,
      data.frame(replicate = r, seed = s, message = conditionMessage(result))) else
        metrics <- rbind(metrics, cbind(replicate = r, seed = s, result))
  }
  unavailable <- c(if (is.null(ref)) "Reference cell-tail metric unavailable: no reference supplied.",
    if (!all(c("Fe2O3", "MgO") %in% fit$parts)) "Fe2O3-MgO correlation unavailable: named parts absent.",
    if (any(template$metric == "closed_part_correlation" & !is.finite(template$observed)) ||
      any(metrics$metric == "closed_part_correlation" &
        (!is.finite(metrics$observed) | !is.finite(metrics$generated))))
      "Undefined closed-part correlations unavailable: at least one part is constant.")
  list(scores = data.frame(row_index = seq_len(nrow(x)), log_score = as.numeric(scores)),
    score_summary = list(n = nrow(x), mean_log_score = mean(scores),
      measure = "ilr", basis = "fitted_ordered_orthonormal_Helmert_ILR_Lebesgue_volume"),
    replicate_metrics = metrics,
    mc_summary = .validation_mc_summary(metrics, template, replicates), failures = failures,
    metadata = list(configuration = list(n_draws = n_draws, replicates = replicates,
      seed = seed, replicate_seeds = seed + seq_len(replicates),
      reference_rows = if (is.null(ref)) 0L else nrow(ref), quantile_type = 7L,
      quantile_probabilities = c(.05, .25, .5, .75, .95),
      reference_probabilities = c(.01, .99), covariance_denominator = "n-1"),
      parts = fit$parts, units = fit$units, RNGkind = RNGkind(),
      R_version = as.character(getRversion()),
      package_version = as.character(getNamespaceVersion("simCoDaR")),
      withr_version = as.character(getNamespaceVersion(asNamespace("withr"))),
      rng_policy = "each_seed_plus_r_scoped_restore", unavailable = unavailable,
      limitations = c("Positive ungrouped full support only; units/nonoverlap confirmed by caller.",
        "MCSE is conditional draw noise, not field uncertainty or scientific validation.",
        "No fitting, tuning, total generation, anonymity or differential privacy.")))
}
