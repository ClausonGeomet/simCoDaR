# Selected GPL-3 numerical code recovered from CoDA R/composition_study.R
# at 989577aaa449d28ac554d806ef57daa47122069c (historical c2184f4).
# Origin: Research Team; history and modifications: Matt Clauson.
# See inst/PROVENANCE.md.
.mvnorm_logpdf <- function(z, mu, sigma) {
  z <- as.matrix(z)
  p <- ncol(z)
  u <- chol(sigma)
  centered <- sweep(z, 2L, mu, "-")
  whitened <- t(forwardsolve(t(u), t(centered)))
  -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(u))) + rowSums(whitened^2))
}

.mvnorm_draw <- function(n, mu, sigma) {
  noise <- matrix(stats::rnorm(n * length(mu)), nrow = n)
  sweep(noise %*% chol(sigma), 2L, mu, "+")
}

# Fixed-df multivariate Student t in ilr coordinates. Keeping df fixed avoids
# an unstable extra tail parameter in small assay datasets. Sigma is the t
# scale matrix (its covariance is nu/(nu-2) * Sigma), not the covariance.
.mvt_logpdf <- function(z, mu, sigma, df) {
  z <- as.matrix(z)
  p <- ncol(z)
  u <- chol(sigma)
  centered <- sweep(z, 2L, mu, "-")
  whitened <- t(forwardsolve(t(u), t(centered)))
  delta <- rowSums(whitened^2)
  lgamma((df + p) / 2) - lgamma(df / 2) -
    (p * log(df * pi) + 2 * sum(log(diag(u)))) / 2 -
    (df + p) * log1p(delta / df) / 2
}

.mvt_draw <- function(n, fit) {
  noise <- .mvnorm_draw(n, rep(0, length(fit$mu)), fit$sigma)
  scale <- sqrt(fit$df / stats::rchisq(n, df = fit$df))
  sweep(noise * scale, 2L, fit$mu, "+")
}

.logsumexp_rows <- function(x) {
  peak <- apply(x, 1L, max)
  peak + log(rowSums(exp(sweep(x, 1L, peak, "-"))))
}

.dirichlet_logpdf <- function(x, fit) {
  a <- fit$alpha
  lgamma(sum(a)) - sum(lgamma(a)) + as.vector(log(x) %*% (a - 1))
}

.dirichlet_draw <- function(n, fit) {
  x <- matrix(stats::rgamma(n * length(fit$alpha),
                            shape = rep(fit$alpha, each = n)), nrow = n)
  totals <- rowSums(x)
  if (any(!is.finite(x)) || any(x <= 0) || any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Dirichlet gamma draws underflowed or were non-finite; increase concentration parameters.")
  }
  x / totals
}

.mixture_logpdf <- function(z, fit) {
  z <- as.matrix(z)
  terms <- vapply(seq_along(fit$weights), function(j) {
    log(fit$weights[j]) + .mvnorm_logpdf(z, fit$means[j, ], fit$covariances[[j]])
  }, numeric(nrow(z)))
  # Preserve n by k even for one observation or one component.
  terms <- matrix(terms, nrow = nrow(z), ncol = length(fit$weights))
  .logsumexp_rows(terms)
}

.mixture_draw <- function(n, fit) {
  component <- sample.int(length(fit$weights), n, replace = TRUE, prob = fit$weights)
  p <- ncol(fit$means)
  out <- matrix(NA_real_, n, p)
  for (j in seq_along(fit$weights)) {
    where <- which(component == j)
    if (length(where)) out[where, ] <- .mvnorm_draw(length(where), fit$means[j, ], fit$covariances[[j]])
  }
  out
}


# Regularize once; never silently increase a declared ridge.
.regularize_covariance <- function(s, ridge = 1e-6) {
  s <- as.matrix(s)
  if (!is.numeric(s) || nrow(s) != ncol(s) || !ncol(s) ||
      any(!is.finite(s))) stop("Nonfinite or malformed covariance.")
  s <- (s + t(s)) / 2
  before <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
  amount <- ridge * max(mean(diag(s)), 1e-8)
  out <- s + diag(amount, ncol(s))
  after <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(out)) || min(after) <= 0)
    stop("Declared covariance ridge did not yield a positive definite matrix.")
  chol(out)
  attr(out, "regularization") <- list(
    ridge = ridge, added = amount, eigenvalues_before = before,
    eigenvalues_after = after,
    rank_before = sum(before > max(abs(before)) * ncol(s) * .Machine$double.eps),
    condition_before = if (min(before) <= 0) Inf else max(before) / min(before),
    condition_after = max(after) / min(after))
  out
}

.fit_dirichlet <- function(x, maxit = 2000L) {
  mu <- colMeans(x)
  moment_total <- mu * (1 - mu) / pmax(apply(x, 2L, stats::var), 1e-12) - 1
  moment_total <- stats::median(moment_total[is.finite(moment_total) & moment_total > 0])
  if (!is.finite(moment_total)) moment_total <- 10
  start <- log(pmax(mu * min(moment_total, 1e5), 1e-5))
  logsum <- colSums(log(x))
  n <- nrow(x)
  objective <- function(a) {
    alpha <- exp(a)
    -(n * (lgamma(sum(alpha)) - sum(lgamma(alpha))) + sum((alpha - 1) * logsum))
  }
  score <- function(a) {
    alpha <- exp(a)
    n * (digamma(sum(alpha)) - digamma(alpha)) + logsum
  }
  opt <- stats::optim(start, objective, gr = function(a) -exp(a) * score(a),
                      method = "L-BFGS-B", lower = rep(-16, ncol(x)),
                      upper = rep(16, ncol(x)), control = list(maxit = maxit))
  if (opt$convergence != 0 || !is.finite(opt$value))
    stop("Dirichlet maximum-likelihood fit failed to converge: code ",
         opt$convergence, "; ", opt$message, "; objective ", opt$value)
  distance <- min(16 - abs(opt$par))
  if (distance <= 1e-5)
    stop("Dirichlet solution saturated log-alpha bounds; unrestricted MLE unsupported.")
  list(alpha = stats::setNames(exp(opt$par), colnames(x)),
       loglik = -opt$value, convergence = opt$convergence,
       diagnostics = list(status = "converged", message = opt$message,
         evaluations = opt$counts, iterations = NA_integer_,
         iteration_count_available = FALSE, score = score(opt$par),
         log_alpha_bounds = c(-16, 16), minimum_bound_distance = distance))
}

.fit_logistic_normal <- function(z) {
  mu <- colMeans(z)
  sigma <- .regularize_covariance(stats::cov(z), ridge = 1e-6)
  list(mu = mu, sigma = sigma, loglik = sum(.mvnorm_logpdf(z, mu, sigma)),
       convergence = 0L,
       diagnostics = list(status = "analytic_regularized", covariance_denominator = nrow(z) - 1L,
                          regularization = attr(sigma, "regularization")))
}

.check_likelihood_update <- function(value, previous, tolerance, label, iteration) {
  if (!is.finite(value)) stop(label, " nonfinite likelihood at iteration ", iteration)
  if (is.finite(previous) && value < previous - tolerance * (1 + abs(previous)))
    stop(label, " material likelihood decrease with declared ridge at iteration ",
         iteration, "; previous ", previous, "; updated ", value)
  is.finite(previous) && abs(value - previous) < tolerance * (1 + abs(previous))
}

.fit_student_t <- function(z, df = 5, max_iter = 500L, tolerance = 1e-7) {
  n <- nrow(z)
  p <- ncol(z)
  mu <- apply(z, 2L, stats::median)
  sigma <- .regularize_covariance(stats::cov(z), ridge = 1e-5)
  initial <- attr(sigma, "regularization")
  previous <- -Inf
  for (iter in seq_len(max_iter)) {
    centered <- sweep(z, 2L, mu, "-")
    delta <- rowSums(t(forwardsolve(t(chol(sigma)), t(centered)))^2)
    weights <- (df + p) / (df + delta)
    next_mu <- colSums(z * weights) / sum(weights)
    residual <- sweep(z, 2L, next_mu, "-")
    next_sigma <- .regularize_covariance(crossprod(residual * sqrt(weights)) / n, 1e-6)
    value <- sum(.mvt_logpdf(z, next_mu, next_sigma, df))
    likelihood_converged <- .check_likelihood_update(value, previous, tolerance, "Student-t", iter)
    # A nearly flat likelihood alone can stop while location/scale still move.
    # Require both parameter updates to settle at the same declared tolerance.
    location_change <- max(abs(next_mu - mu) / (1 + abs(mu)))
    scale_change <- max(abs(next_sigma - sigma)) / max(abs(sigma))
    parameter_change <- max(location_change, scale_change)
    converged <- likelihood_converged && parameter_change < tolerance
    mu <- next_mu
    sigma <- next_sigma
    if (converged) return(list(mu = mu, sigma = sigma, df = df,
      covariance = df / (df - 2) * sigma, loglik = value, iterations = iter,
      convergence = 0L, diagnostics = list(status = "converged", tolerance = tolerance,
        convergence_criterion = "likelihood_and_relative_location_scale_updates",
        location_change = location_change, scale_change = scale_change,
        scale_is_not_covariance = TRUE, initial_regularization = initial,
        regularization = attr(sigma, "regularization"))))
    previous <- value
  }
  stop("Student-t failed to converge: iteration ", max_iter, "; last likelihood ", previous,
       "; last relative parameter change ", parameter_change)
}

.check_mixture_dispersion <- function(z, membership, component, stage) {
  support <- z[membership > 0, , drop = FALSE]
  # Repeats do not identify a component scale. Exact variation is sufficient:
  # rank-deficient but varying support may still use the declared ridge.
  if (!nrow(support) || !any(sweep(support, 2L, support[1L, ], "-") != 0))
    stop("Degenerate ", stage, " component ", component,
         ": no within-component dispersion (repeated-point support).")
}

.fit_gaussian_mixture <- function(z, k, starts = 5L, max_iter = 300L, tolerance = 1e-7) {
  n <- nrow(z)
  p <- ncol(z)
  if (n < max(20L, 3 * k * (p + 1))) stop("Too few training rows for mixture k=", k)
  if (k == 1L) {
    fit <- .fit_logistic_normal(z)
    return(list(weights = 1, means = matrix(fit$mu, nrow = 1L),
      covariances = list(fit$sigma), loglik = fit$loglik, convergence = 0L,
      diagnostics = fit$diagnostics))
  }
  base_cov <- .regularize_covariance(stats::cov(z), 1e-5)
  best <- NULL
  attempts <- vector("list", starts)
  for (attempt in seq_len(starts)) {
    iteration <- 0L
    last_value <- NA_real_
    fit <- tryCatch({
      centers <- z[sample.int(n, k), , drop = FALSE]
      cluster <- max.col(-vapply(seq_len(k), function(j)
        rowSums(sweep(z, 2L, centers[j, ], "-")^2), numeric(n)), ties.method = "first")
      if (length(unique(cluster)) < k) stop("Degenerate initial allocation")
      for (j in seq_len(k)) .check_mixture_dispersion(z, cluster == j, j, "initial")
      weights <- tabulate(cluster, nbins = k) / n
      means <- do.call(rbind, lapply(seq_len(k), function(j)
        colMeans(z[cluster == j, , drop = FALSE])))
      covariances <- replicate(k, base_cov, simplify = FALSE)
      previous <- -Inf
      converged <- FALSE
      for (iter in seq_len(max_iter)) {
        iteration <- iter
        terms <- vapply(seq_len(k), function(j)
          log(weights[j]) + .mvnorm_logpdf(z, means[j, ], covariances[[j]]), numeric(n))
        norms <- .logsumexp_rows(terms)
        responsibility <- exp(terms - norms)
        mass <- colSums(responsibility)
        if (any(!is.finite(mass)) || any(mass < p + 2)) stop("Degenerate component mass < p+2")
        weights <- mass / n
        for (j in seq_len(k)) {
          .check_mixture_dispersion(z, responsibility[, j], j, "updated")
          means[j, ] <- colSums(z * responsibility[, j]) / mass[j]
          residual <- sweep(z, 2L, means[j, ], "-")
          covariances[[j]] <- .regularize_covariance(
            crossprod(residual * sqrt(responsibility[, j])) / mass[j], 1e-4)
        }
        current <- list(weights = weights, means = means, covariances = covariances)
        last_value <- sum(.mixture_logpdf(z, current))
        converged <- .check_likelihood_update(last_value, previous, tolerance, "Mixture", iter)
        if (converged) break
        previous <- last_value
      }
      if (!converged) stop("Iteration exhaustion")
      # Check the final updated responsibilities as well, before admitting a start.
      final_terms <- vapply(seq_len(k), function(j)
        log(weights[j]) + .mvnorm_logpdf(z, means[j, ], covariances[[j]]), numeric(n))
      final_responsibility <- exp(final_terms - .logsumexp_rows(final_terms))
      final_mass <- colSums(final_responsibility)
      if (any(!is.finite(final_mass)) || any(final_mass < p + 2)) stop("Degenerate final component mass")
      for (j in seq_len(k)) .check_mixture_dispersion(z, final_responsibility[, j], j, "final")
      c(current, list(loglik = last_value, iterations = iter, convergence = 0L,
        diagnostics = list(status = "converged", effective_mass = final_mass,
          initial_regularization = attr(base_cov, "regularization"),
          regularization = lapply(covariances, attr, "regularization"))))
    }, error = function(e) list(failure = conditionMessage(e)))
    attempts[[attempt]] <- list(start = attempt, iterations = iteration, loglik = last_value,
      status = if (is.null(fit$failure)) "converged" else "failed", reason = fit$failure)
    if (is.null(fit$failure) && (is.null(best) || fit$loglik > best$loglik)) best <- fit
  }
  if (is.null(best)) stop("All Gaussian mixture starts failed: ",
    paste(vapply(attempts, function(a) paste0(a$reason, " (iteration ", a$iterations,
      ", likelihood ", a$loglik, ")"), character(1)), collapse = "; "))
  best$diagnostics$starts <- attempts
  best
}
