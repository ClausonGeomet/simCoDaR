test_that("positive heldout diagnostics preserve RNG and common scoring", {
  data <- withr::with_seed(732L, {
    z <- matrix(rnorm(240), 80, 3)
    x <- exp(cbind(z, 0))
    as.data.frame(x / rowSums(x))
  })
  names(data) <- c("Fe2O3", "MgO", "SiO2", "CaO")
  fit <- fit_assay_simulator(data[1:60, ], names(data), "logistic_normal")
  set.seed(431L)
  before <- .Random.seed
  result <- validate_assay_simulator(fit, data[61:80, ], data[1:60, ],
    n_draws = 12L, replicates = 3L, seed = 100L)
  expect_identical(.Random.seed, before)
  expect_equal(result$scores$log_score,
    as.numeric(simCoDaR:::.assay_log_density(fit, data[61:80, ])))
  expect_equal(result$score_summary$mean_log_score, mean(result$scores$log_score))
  expect_identical(result$score_summary$measure, "ilr")
  expect_equal(nrow(result$failures), 0L)
  expect_identical(unique(result$replicate_metrics$seed), 101:103)
  energy <- subset(result$replicate_metrics, metric == "ilr_energy")$generated
  summary <- subset(result$mc_summary, metric == "ilr_energy")
  expect_equal(summary$mean, mean(energy))
  expect_equal(summary$mcse, sd(energy) / sqrt(3))
  expect_true(all(!summary$incomplete))
  expect_equal(validate_assay_simulator(fit, data[61:80, ], data[1:60, ],
    12L, 3L, 100L), result)
  # Reconstruct a single diagnostic independently from the declared draw seed.
  y <- as.matrix(simulate_assays(fit, 12, seed = 101L))
  b <- fit$groups[[1]]$patterns[[1]]$basis
  x <- as.matrix(data[61:80, ])
  z <- log(x) %*% t(b)
  w <- log(y) %*% t(b)
  cross <- vapply(seq_len(nrow(z)), function(i)
    mean(sqrt(rowSums(sweep(w, 2, z[i, ], "-")^2))), numeric(1))
  within <- function(a) sum(as.matrix(dist(a))) / nrow(a)^2
  expect_equal(energy[1], 2 * mean(cross) - within(z) - within(w))
  cov_error <- sqrt(sum((cov(z) - cov(w))^2)) / sqrt(sum(cov(z)^2))
  expect_equal(subset(result$replicate_metrics,
    replicate == 1 & metric == "ilr_covariance_relative_error")$generated, cov_error)
  row <- subset(result$replicate_metrics, replicate == 1 & metric == "logratio_quantile" &
    part_a == "Fe2O3" & part_b == "MgO" & probability == .5)
  expect_equal(row$generated, unname(quantile(log(y[, 1] / y[, 2]), .5, type = 7)))
  expect_equal(row$discrepancy, abs(row$generated - row$observed))
  expect_false(any(c("heldout", "reference", "draw", "data", "call") %in% names(result)))
})

test_that("validator restrictions and unavailable metrics are explicit", {
  d <- withr::with_seed(23, data.frame(A = exp(rnorm(60)), B = exp(rnorm(60))))
  f <- fit_assay_simulator(d, c("A", "B"))
  out <- validate_assay_simulator(f, d[1:5, ], n_draws = 3, replicates = 2)
  expect_true(any(grepl("no reference", out$metadata$unavailable)))
  expect_true(all(is.na(subset(out$mc_summary,
    metric == "reference_outside_cell_fraction")$mcse)))
  expect_error(validate_assay_simulator(f, d, seed = NULL), "not NULL")
  expect_error(validate_assay_simulator(f, d, seed = .Machine$integer.max - 2,
    replicates = 3), "seed \\+ replicates")
  expect_error(validate_assay_simulator(f, d, n_draws = 1), "n_draws")
  expect_error(validate_assay_simulator(f, d[1, ]), "two rows")
  bad <- d
  bad$A[1] <- 0
  expect_error(validate_assay_simulator(f, bad), "positive")
  expect_error(validate_assay_simulator(f, d, reference = bad), "positive")
  d$g <- rep(c("a", "b"), each = 30)
  grouped <- fit_assay_simulator(d, c("A", "B"), group_col = "g")
  expect_error(validate_assay_simulator(grouped, d), "ungrouped")
  d$g <- NULL
  d$A[1:30] <- 0
  structural <- fit_assay_simulator(d, c("A", "B"), zero_mode = "structural",
    structural_parts = "A")
  expect_error(validate_assay_simulator(structural, d), "positive support")
})

test_that("replicate sampler failures remain visible with no invented MCSE", {
  d <- withr::with_seed(10, data.frame(A = exp(rnorm(40)), B = exp(rnorm(40))))
  f <- fit_assay_simulator(d, c("A", "B"))
  original <- simCoDaR::simulate_assays
  testthat::local_mocked_bindings(simulate_assays = function(fit, n, seed, ...) {
    if (seed != 102L) stop("declared smoke failure")
    original(fit, n, seed, ...)
  }, .package = "simCoDaR")
  set.seed(222L)
  before <- .Random.seed
  out <- validate_assay_simulator(f, d, n_draws = 4, replicates = 3, seed = 100)
  expect_identical(.Random.seed, before)
  expect_equal(out$failures$replicate, c(1L, 3L))
  expect_true(all(grepl("declared smoke failure", out$failures$message)))
  e <- subset(out$mc_summary, metric == "ilr_energy")
  expect_equal(e$n_success, 1L)
  expect_equal(e$n_failure, 2L)
  expect_true(e$incomplete)
  expect_true(is.na(e$mcse))
})

test_that("constant correlations and entirely failed repetitions stay unavailable", {
  d <- withr::with_seed(63, data.frame(Fe2O3 = exp(rnorm(40)), MgO = exp(rnorm(40))))
  f <- fit_assay_simulator(d, names(d))
  heldout <- data.frame(Fe2O3 = rep(1, 3), MgO = rep(2, 3))
  out <- validate_assay_simulator(f, heldout, n_draws = 3, replicates = 2)
  expect_true(any(grepl("constant", out$metadata$unavailable)))
  expect_true(all(is.na(subset(out$replicate_metrics,
    metric == "closed_part_correlation")$observed)))
  expect_true(all(is.na(subset(out$mc_summary,
    metric == "closed_part_correlation" & statistic == "discrepancy")$mcse)))
  testthat::local_mocked_bindings(simulate_assays = function(...) stop("all failed"),
    .package = "simCoDaR")
  failed <- validate_assay_simulator(f, heldout, n_draws = 3, replicates = 2)
  expect_equal(nrow(failed$failures), 2)
  expect_equal(nrow(failed$replicate_metrics), 0)
  expect_true(all(failed$mc_summary$n_failure == 2))
  expect_true(all(failed$mc_summary$n_success == 0))
  expect_true(all(is.na(failed$mc_summary$mcse)))
  expect_true(any(grepl("constant", failed$metadata$unavailable)))
})
