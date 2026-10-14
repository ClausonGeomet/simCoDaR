test_that("recovered families return closed draws and restore RNG", {
  data <- withr::with_seed(731, {
    x <- matrix(rgamma(360, shape = rep(c(3, 5, 7), each = 120)), 120, 3)
    as.data.frame(x)
  })
  names(data) <- c("a", "b", "c")
  withr::with_seed(219, {
    before <- .Random.seed
    for (model in c("dirichlet", "logistic_normal", "ilr_student_t", "ilr_mixture")) {
      fit <- fit_assay_simulator(data, names(data), model,
        mixture_components = 1L, measurement_basis = "abstract_composition")
      expect_identical(.Random.seed, before)
      out <- simulate_assays(fit, 12)
      expect_identical(.Random.seed, before)
      expect_equal(rowSums(out), rep(1, 12), tolerance = 1e-12)
      expect_true(all(as.matrix(out) > 0))
      expect_identical(out, simulate_assays(fit, 12))
      expect_false("training_totals" %in% names(fit))
      expect_false("observed" %in% names(fit$groups[[1]]$patterns[[1]]))
    }
  })
})

test_that("structural singleton supports stay absent without a row fallback", {
  d <- data.frame(a = c(1, 0), b = c(0, 1), g = c("first", "second"))
  f <- fit_assay_simulator(d, c("a", "b"), group_col = "g",
    zero_mode = "structural", structural_parts = c("a", "b"))
  x <- simulate_assays(f, 3, group = "first")
  expect_identical(x$a, rep(1, 3))
  expect_identical(x$b, rep(0, 3))
  expect_identical(names(x), c("a", "b", "g"))
  expect_error(fit_assay_simulator(d, c("a", "b")), "Zero")
  expect_error(fit_assay_simulator(data.frame(a = c(1, 2), b = 1), c("a", "b")), "sparse")
})

test_that("common-measure Jacobian and analytic covariance match contracts", {
  d <- data.frame(a = exp(seq(-2, 2, length.out = 40)), b = 1)
  f <- fit_assay_simulator(d, c("a", "b"))
  h <- f$groups[[1]]$patterns[[1]]
  x <- as.matrix(d) / rowSums(d)
  z <- log(x) %*% t(h$basis)
  expected <- var(as.numeric(z)) * (1 + 1e-6)
  expect_equal(as.numeric(h$model$sigma), expected)
  ilr <- simCoDaR:::.assay_log_density(f, d)
  simplex <- simCoDaR:::.assay_log_density(f, d, measure = "simplex")
  expect_equal(as.numeric(simplex - ilr), -log(2) / 2 - rowSums(log(x)))
  f$groups[[1]]$patterns[[1]]$model$sigma[1, 1] <- -1
  expect_error(simulate_assays(f, 3), "positive definite")
})

test_that("unsafe units, overlap and iteration exhaustion fail explicitly", {
  d <- withr::with_seed(18, data.frame(a = exp(rnorm(60)), b = exp(rnorm(60))))
  expect_error(fit_assay_simulator(d, c("a", "b"), part_units = c(a = "ppm", b = "wt%")), "Mixed")
  names(d) <- c("Fe_wt_pct", "Fe2O3_wt_pct")
  expect_error(fit_assay_simulator(d, names(d)), "Overlapping")
  names(d) <- c("a", "b")
  withr::with_seed(888, {
    before <- .Random.seed
    expect_error(fit_assay_simulator(d, names(d), "ilr_student_t",
      control = list(student_max_iter = 1)), "failed to converge")
    expect_identical(.Random.seed, before)
  })
  expect_error(fit_assay_simulator(d, names(d), mixture_components = 1.5), "integers")
  expect_error(fit_assay_simulator(d, names(d), seed = -1), "seed")
})

test_that("multicomponent mixtures preserve single-row density dimensions", {
  d <- withr::with_seed(400, {
    z <- matrix(rnorm(400, sd = .3), 200, 2) + c(rep(-2, 100), rep(2, 100))
    as.data.frame(simCoDaR:::.from_ilr(z, simCoDaR:::.ilr_basis(3)))
  })
  names(d) <- c("a", "b", "c")
  f <- fit_assay_simulator(d, names(d), "ilr_mixture", mixture_components = 2:3)
  h <- f$groups[[1]]$patterns[[1]]
  expect_identical(h$model$selected_components, 2L)
  score <- simCoDaR:::.assay_log_density(f, d[1, , drop = FALSE])
  expect_length(score, 1L)
  expect_true(is.finite(score))
  expect_equal(sum(simulate_assays(f, 1)), 1)
})
