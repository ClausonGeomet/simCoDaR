# Common measure on each declared active support, plus discrete support mass.
# Dirichlet's original coordinate-simplex density uses dx_1 ... dx_(d-1).
# ILR volume: log f_z = log f_x - log J,
# where log J = -log(d)/2 - sum(log x). Singleton support uses point mass.
.assay_log_density <- function(fit, data, group = NULL,
    measure = c("ilr", "simplex_coordinate", "simplex"), conditional_group = TRUE) {
  .validate_empirical_fit(fit)
  measure <- match.arg(measure)
  if (measure == "simplex") measure <- "simplex_coordinate"
  if (!is.logical(conditional_group) || length(conditional_group) != 1L || is.na(conditional_group))
    stop("conditional_group must be TRUE or FALSE.")
  if (!is.data.frame(data)) stop("Density data must be a data frame.")
  .empirical_names(names(data), "Density data column names")
  if (!all(fit$parts %in% names(data)) || !all(vapply(data[fit$parts], is.numeric, logical(1))))
    stop("Density data must contain unique fitted numeric parts.")
  x <- .empirical_close(as.matrix(data[fit$parts]))
  n <- nrow(x)
  labels <- NULL
  if (!is.null(group)) {
    if (is.null(fit$group_col)) stop("group requires a grouped fit.")
    if (!length(group) %in% c(1L, n)) stop("group must have length one or number of observations.")
    labels <- rep(.empirical_labels(group), length.out = n)
  } else if (!is.null(fit$group_col) && fit$group_col %in% names(data)) {
    labels <- .empirical_labels(data[[fit$group_col]])
  }
  if (!is.null(fit$group_col) && is.null(labels)) stop("Grouped density requires group labels.")
  if (!is.null(labels) && !all(labels %in% names(fit$groups))) stop("Unknown density group.")
  keys <- apply(x > 0, 1L, .assay_pattern_key)
  score_group <- function(label) {
    g <- fit$groups[[label]]
    value <- rep(-Inf, n)
    for (key in intersect(unique(keys), names(g$patterns))) {
      rows <- which(keys == key)
      h <- g$patterns[[key]]
      active <- x[rows, h$active_parts, drop = FALSE]
      if (h$kind == "singleton") logpdf <- rep(0, length(rows)) else {
        z <- .to_ilr(active, h$basis)
        m <- h$model
        jacobian <- -log(ncol(active)) / 2 - rowSums(log(active))
        logpdf <- switch(h$kind,
          dirichlet = .dirichlet_logpdf(active, m) - jacobian,
          logistic_normal = .mvnorm_logpdf(z, m$mu, m$sigma),
          ilr_student_t = .mvt_logpdf(z, m$mu, m$sigma, m$df),
          ilr_mixture = .mixture_logpdf(z, m))
        if (measure == "simplex_coordinate") logpdf <- logpdf + jacobian
        if (any(!is.finite(logpdf))) stop("Nonfinite active-support density.")
      }
      value[rows] <- logpdf + log(g$pattern_probabilities[[key]])
    }
    if (!conditional_group) value <- value + log(fit$group_probabilities[[label]])
    value
  }
  if (!is.null(labels)) {
    out <- rep(-Inf, n)
    for (label in unique(labels)) {
      rows <- which(labels == label)
      out[rows] <- score_group(label)[rows]
    }
  } else {
    out <- score_group(".all")
  }
  attr(out, "measure") <- paste0(measure, "_volume_on_active_support_plus_discrete_mass")
  attr(out, "support") <- keys
  attr(out, "conditional_group") <- conditional_group
  attr(out, "basis") <- "fitted_ordered_Helmert_on_each_active_support"
  out
}
