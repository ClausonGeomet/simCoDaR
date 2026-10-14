# Adapted from CoDA R/scenario.R at
# c2184f44994a7c3f2d4f36ec74c749e0267a4127; GPL-3.
# Origin: Research Team; history: Matt Clauson. See inst/PROVENANCE.md.
.bound_slack <- 1e-12
.closure_slack <- 1e-10
.valid_names <- function(x, n) {
  is.character(x) && length(x) == n && !anyNA(x) &&
    all(nzchar(trimws(x))) && !anyDuplicated(x)
}
.integer_scalar <- function(x, label, lower = 1, upper = .Machine$integer.max) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != 1L ||
      !is.finite(x) || x < lower || x > upper || x != floor(x))
    stop(label, " must be an integer scalar in [", lower, ", ", upper, "].")
  as.integer(x)
}
.close_safe <- function(x) {
  scaled <- x / max(x)
  out <- scaled / sum(scaled)
  if (any(!is.finite(out)) || any(out[x > 0] <= 0))
    stop("Normalization cannot represent the positive support; rescale inputs.")
  out
}
.domain_vector <- function(x, facies, label, scalar = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
      any(!is.finite(x)) || any(x < 0) || (label == "weights" && any(x <= 0)))
    stop(label, " must contain finite ", if (label == "weights") "positive" else "nonnegative", " numbers.")
  if (!is.null(names(x))) {
    if (!.valid_names(names(x), length(facies)) || !setequal(names(x), facies))
      stop("Named ", label, " must match every domain exactly.")
    x <- x[facies]
  } else if (scalar && length(x) == 1L) x <- rep(x, length(facies))
  if (length(x) != length(facies)) stop(label, " must match the number of domains.")
  stats::setNames(x, facies)
}
.validate_composition_bounds <- function(bounds, target, sd) {
  parts <- colnames(target)
  out <- rbind(lower = rep(0, length(parts)), upper = rep(1, length(parts)))
  colnames(out) <- parts
  if (!is.null(bounds)) {
    if (!is.list(bounds) || !.valid_names(names(bounds), length(bounds)) ||
        !all(names(bounds) %in% parts))
      stop("composition_bounds must be a named list of selected part ranges.")
    for (part in names(bounds)) {
      r <- bounds[[part]]
      if (!is.numeric(r) || !is.null(dim(r)) || length(r) != 2L ||
          any(!is.finite(r)) || r[1L] < 0 || r[2L] > 1 || r[1L] > r[2L])
        stop("Each bound must be a finite [lower, upper] within [0, 1].")
      out[, part] <- r
    }
  }
  for (i in seq_len(nrow(target))) {
    active <- target[i, ] > 0
    if (sum(out[1L, active]) > 1 + .bound_slack ||
        sum(out[2L, active]) < 1 - .bound_slack || any(out[1L, !active] > 0))
      stop("Bounds are incompatible with domain support: ", rownames(target)[i])
    if (any(target[i, ] < out[1L, ] - .bound_slack) ||
        any(target[i, ] > out[2L, ] + .bound_slack))
      stop("Domain target lies outside bounds: ", rownames(target)[i])
    if (sd[i] > 0 && sum(active) > 1L &&
        any(out[1L, active] == out[2L, active]))
      stop("Zero-width active bounds have zero acceptance probability with positive spread.")
  }
  out
}
#' Construct an explicit, uncalibrated domain composition scenario.
#' @export
new_assay_scenario <- function(endmembers, weights = NULL, logratio_sd = 0.25,
                               structural_parts = NULL, composition_bounds = NULL) {
  if (!(is.matrix(endmembers) || is.data.frame(endmembers)) ||
      (is.data.frame(endmembers) && !all(vapply(endmembers, is.numeric, logical(1)))))
    stop("endmembers must be a numeric matrix or data frame.")
  target <- as.matrix(endmembers)
  if (!is.numeric(target) || nrow(target) < 1L || ncol(target) < 2L ||
      !.valid_names(rownames(target), nrow(target)) ||
      !.valid_names(colnames(target), ncol(target)) ||
      any(!is.finite(target)) || any(target < 0) || any(rowSums(target > 0) == 0))
    stop("endmembers needs unique domain/part names and finite nonnegative targets with positive totals.")
  parts <- colnames(target)
  facies <- rownames(target)
  if (!is.null(structural_parts) &&
      (!is.character(structural_parts) || !is.null(dim(structural_parts)) ||
       anyNA(structural_parts) ||
       anyDuplicated(structural_parts) || !all(structural_parts %in% parts)))
    stop("structural_parts must contain distinct selected part names.")
  if (!all(parts[colSums(target == 0) > 0] %in% structural_parts))
    stop("Every zero target requires its part in structural_parts; censored zeros are unsupported.")
  closed <- target
  for (i in seq_len(nrow(target))) closed[i, ] <- .close_safe(target[i, ])
  if (is.null(weights)) weights <- rep(1, length(facies))
  weights <- .close_safe(.domain_vector(weights, facies, "weights"))
  sd <- .domain_vector(logratio_sd, facies, "logratio_sd", scalar = TRUE)
  structure(list(endmembers = target, facies = facies, parts = parts,
    weights = weights, logratio_sd = sd, structural_parts = structural_parts,
    composition_bounds = .validate_composition_bounds(composition_bounds, closed, sd),
    calibrated = FALSE), class = "assay_scenario")
}
.revalidate_scenario <- function(spec) {
  fields <- c("endmembers", "facies", "parts", "weights", "logratio_sd",
              "structural_parts", "composition_bounds", "calibrated")
  if (!inherits(spec, "assay_scenario") || !is.list(spec) ||
      !.valid_names(names(spec), length(fields)) || !setequal(names(spec), fields) ||
      !identical(spec$calibrated, FALSE)) stop("spec must be an uncalibrated assay_scenario with supported fields.")
  b <- spec$composition_bounds
  if (!is.matrix(b) || !is.numeric(b) || nrow(b) != 2L ||
      !identical(rownames(b), c("lower", "upper")) ||
      !identical(colnames(b), spec$parts)) stop("Invalid scenario bounds matrix.")
  rebuilt <- new_assay_scenario(spec$endmembers, spec$weights, spec$logratio_sd,
    spec$structural_parts, stats::setNames(lapply(seq_len(ncol(b)), function(j) b[, j]), colnames(b)))
  if (!identical(spec$facies, rebuilt$facies) || !identical(spec$parts, rebuilt$parts))
    stop("Scenario domain/part labels do not match endmembers.")
  rebuilt
}
#' Draw closed synthetic compositions with support-specific ILR variation.
#' @export
simulate_assay_scenario <- function(spec, n, seed = NULL, max_draw_attempts = 1000L) {
  spec <- .revalidate_scenario(spec)
  n <- .integer_scalar(n, "n")
  max_draw_attempts <- .integer_scalar(max_draw_attempts, "max_draw_attempts")
  if (!is.null(seed)) seed <- .integer_scalar(seed, "seed", 0, .Machine$integer.max - 1)
  rng_kind <- RNGkind()
  generate <- function() {
    facies <- sample(spec$facies, n, replace = TRUE, prob = spec$weights)
    ids <- paste0("sample_", seq_len(n))
    out <- matrix(0, n, length(spec$parts), dimnames = list(ids, spec$parts))
    attempted <- rejected <- 0
    for (g in spec$facies) {
      rows <- which(facies == g)
      if (!length(rows)) next
      center <- .close_safe(spec$endmembers[g, ])
      active <- which(center > 0)
      sd <- spec$logratio_sd[g]
      if (length(active) == 1L || sd == 0) {
        out[rows, ] <- matrix(rep(center, each = length(rows)), nrow = length(rows))
        attempted <- attempted + length(rows)
        next
      }
      basis <- .ilr_basis(length(active))
      mu <- as.numeric(.to_ilr(matrix(center[active], 1L), basis))
      pending <- rows
      for (attempt in seq_len(max_draw_attempts)) {
        z <- matrix(stats::rnorm(length(pending) * length(mu), sd = sd), nrow = length(pending))
        candidate <- .from_ilr(sweep(z, 2L, mu, "+"), basis)
        valid <- rowSums(sweep(candidate, 2L, spec$composition_bounds[1L, active] - .bound_slack, "<") |
          sweep(candidate, 2L, spec$composition_bounds[2L, active] + .bound_slack, ">")) == 0
        attempted <- attempted + length(pending)
        rejected <- rejected + sum(!valid)
        out[pending[valid], active] <- candidate[valid, , drop = FALSE]
        pending <- pending[!valid]
        if (!length(pending)) break
      }
      if (length(pending)) stop("Bounds exhausted max_draw_attempts in domain '", g,
                                "'; relax bounds or reduce logratio_sd.")
    }
    if (any(!is.finite(out)) || any(out < 0) || any(abs(rowSums(out) - 1) > .closure_slack))
      stop("Scenario produced invalid compositions.")
    list(sample_id = ids, composition = as.data.frame(out, check.names = FALSE),
      facies = facies, assumptions = list(spec = spec, seed = seed, RNGkind = rng_kind,
        composition_rejection_rate = rejected / attempted,
        validation = "illustrative_unvalidated_scenario", calibrated = FALSE,
        rng_policy = if (is.null(seed)) "consume_caller_rng" else "scoped_seed_restore_caller_rng"))
  }
  if (is.null(seed)) generate() else withr::with_seed(seed, generate())
}
