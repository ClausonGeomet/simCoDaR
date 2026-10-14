test_that("density aliases share a measure and grouped scores require labels", {
  d <- withr::with_seed(824, {
    x <- matrix(rgamma(450, shape = rep(c(3, 6, 9), each = 150)), 150, 3)
    data.frame(a = x[, 1], b = x[, 2], c = x[, 3],
               rock = rep(c("first", "second"), c(60, 90)))
  })
  for (family in c("dirichlet", "logistic_normal", "ilr_student_t", "ilr_mixture")) {
    f <- fit_assay_simulator(d, c("a", "b", "c"), family,
      group_col = "rock", mixture_components = 1L)
    density <- simCoDaR:::.assay_log_density
    conditional <- density(f, d)
    explicit <- density(f, d, group = d$rock)
    expect_identical(conditional, explicit)
    joint <- density(f, d, conditional_group = FALSE)
    expect_equal(as.numeric(joint - conditional),
      log(unname(f$group_probabilities[d$rock])), tolerance = 1e-12)
    simplex <- density(f, d, measure = "simplex_coordinate")
    expect_identical(simplex, density(f, d, measure = "simplex"))
    expect_identical(attr(simplex, "measure"),
      "simplex_coordinate_volume_on_active_support_plus_discrete_mass")
    x <- as.matrix(d[c("a", "b", "c")])
    x <- x / rowSums(x)
    expect_equal(as.numeric(simplex - conditional),
      -log(3) / 2 - rowSums(log(x)), tolerance = 1e-12)
    for (condition in c(TRUE, FALSE)) {
      expect_error(density(f, d[c("a", "b", "c")], conditional_group = condition),
        "requires group labels")
      unknown <- d
      unknown$rock[1] <- "unknown"
      expect_error(density(f, unknown, conditional_group = condition), "Unknown")
      unknown$rock[1] <- NA_character_
      expect_error(density(f, unknown, conditional_group = condition), "complete")
    }
    # Explicit labels override an absent or invalid observed label column.
    expect_identical(density(f, unknown, group = d$rock), conditional)
    plain <- fit_assay_simulator(d, c("a", "b", "c"), family, mixture_components = 1L)
    expect_true(all(is.finite(density(plain, d[c("a", "b", "c")]))))
    expect_error(density(plain, d, group = "first"), "grouped fit")
  }
})

test_that("complete input names are validated without inspecting unused values", {
  d <- data.frame(a = exp(seq(-1, 1, length.out = 40)), b = 1,
                  unused = rep(NA_character_, 40))
  expect_s3_class(fit_assay_simulator(d, c("a", "b")), "geo_assay_simulator")
  f <- fit_assay_simulator(d, c("a", "b"))
  for (bad in c("", " ", NA_character_, "a")) {
    malformed <- d
    names(malformed)[3] <- bad
    expect_error(fit_assay_simulator(malformed, c("a", "b")), "column names")
    expect_error(simCoDaR:::.assay_log_density(f, malformed), "column names")
  }
})

test_that("Student parameters settle against an independently optimized density", {
  z <- withr::with_seed(925, {
    normal <- matrix(rnorm(1200), 600, 2) %*% chol(matrix(c(1, .45, .45, 1.6), 2))
    sweep(normal * sqrt(5 / rchisq(600, 5)), 2, c(.8, -.4), "+")
  })
  # Parameterize positive scale by a lower Cholesky factor and optimize the
  # fixed-df likelihood directly, independently of the weighted-update fitter.
  decode <- function(theta) {
    l <- matrix(c(exp(theta[3]), theta[4], 0, exp(theta[5])), 2)
    list(mu = theta[1:2], sigma = tcrossprod(l))
  }
  objective <- function(theta) {
    m <- decode(theta)
    if (any(!is.finite(m$sigma))) return(Inf)
    residual <- sweep(z, 2, m$mu, "-")
    inverse <- tryCatch(solve(m$sigma), error = function(e) NULL)
    if (is.null(inverse)) return(Inf)
    delta <- rowSums((residual %*% inverse) * residual)
    -sum(lgamma(3.5) - lgamma(2.5) - log(5 * pi) -
      as.numeric(determinant(m$sigma, logarithm = TRUE)$modulus) / 2 -
      3.5 * log1p(delta / 5))
  }
  reference <- optim(c(.5, 0, 0, 0, 0), objective, method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 1000))
  expect_identical(reference$convergence, 0L)
  expected <- decode(reference$par)
  fitted <- simCoDaR:::.fit_student_t(z)
  expect_lt(max(abs(fitted$mu - expected$mu)), 1e-4)
  expect_lt(max(abs(fitted$sigma - expected$sigma)) / max(abs(expected$sigma)), 1e-4)
  expect_equal(fitted$covariance, 5 / 3 * fitted$sigma)
  expect_lt(fitted$diagnostics$location_change, 1e-7)
  expect_lt(fitted$diagnostics$scale_change, 1e-7)
  expect_error(simCoDaR:::.fit_student_t(z, max_iter = 1),
    "last relative parameter change")
})

test_that("repeated-point mixtures fail while varying deficient support remains eligible", {
  density_data <- data.frame(a = rep(c(1, 100), each = 50), b = 1)
  withr::with_seed(632, {
    before <- .Random.seed
    expect_error(fit_assay_simulator(density_data, c("a", "b"), "ilr_mixture",
      mixture_components = 2L), "Degenerate.*component|Degenerate initial allocation")
    expect_identical(.Random.seed, before)
  })
  admissible <- fit_assay_simulator(density_data, c("a", "b"), "ilr_mixture",
    mixture_components = 1:2)
  m <- admissible$groups[[1]]$patterns[[1]]$model
  expect_identical(m$selected_components, 1L)
  expect_true(is.infinite(m$bic[["2"]]))
  expect_match(m$diagnostics$candidates[[2]]$failure, "Degenerate")
  # Repeated rows remain permitted when each regime really varies, even if
  # its full ILR covariance is rank deficient and needs the existing ridge.
  z <- cbind(rep(c(seq(-4, -2, length.out = 30), seq(2, 4, length.out = 30)), each = 2), 0)
  varied <- as.data.frame(simCoDaR:::.from_ilr(z, simCoDaR:::.ilr_basis(3)))
  names(varied) <- c("a", "b", "c")
  fit <- fit_assay_simulator(varied, names(varied), "ilr_mixture", mixture_components = 2L)
  component <- fit$groups[[1]]$patterns[[1]]$model
  expect_identical(component$selected_components, 2L)
  expect_true(all(vapply(component$diagnostics$regularization,
    function(r) r$rank_before == 1L && r$added > 0, logical(1))))
  expect_true(all(is.finite(simCoDaR:::.assay_log_density(fit, varied))))
})
