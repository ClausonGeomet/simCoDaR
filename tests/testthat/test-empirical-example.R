test_that("installed example freezes all five public manifests without fits", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("examples", "iodp-empirical.R", package = "simCoDaR",
    mustWork = TRUE), envir = env)
  output <- tempfile("iodp-freeze-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  set.seed(24)
  before <- .Random.seed
  config <- env$run_iodp_example(output, "freeze")
  expect_identical(.Random.seed, before)
  expect_equal(config$student_df, 5)
  expect_identical(config$mixture_components, 2:3)
  expect_identical(config$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection"))
  expect_equal(nrow(config$diagnostic_seeds), 400)
  expect_false(any(grepl("fit.rds$", list.files(output))))
  ids <- character()
  for (f in 0:4) {
    manifest <- read.csv(file.path(output, paste0("split", f, "-manifest.csv")))
    expect_equal(sum(manifest$role == "heldout"), 71)
    expect_equal(sum(manifest$role == "train"), 213)
    expect_equal(nrow(manifest), 284)
    if (f > 0) ids <- c(ids, manifest$sample_id[manifest$role == "heldout"])
  }
  expect_equal(length(unique(ids)), 284)
  expect_error(env$run_iodp_example(output, "freeze"), "empty output")
})

test_that("installed example smoke writes one fit and bounded diagnostics", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("examples", "iodp-empirical.R", package = "simCoDaR",
    mustWork = TRUE), envir = env)
  output <- tempfile("iodp-smoke-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  set.seed(817)
  before <- .Random.seed
  result <- env$run_iodp_example(output, "smoke")
  expect_identical(.Random.seed, before)
  expect_equal(nrow(result$failures), 0)
  expect_equal(nrow(result$summary), 1)
  expect_identical(result$summary$model, "logistic_normal")
  expect_equal(length(list.files(output, pattern = "-manifest.csv$")), 5)
  evaluation <- readRDS(file.path(output, "split0-logistic_normal-validation.rds"))
  expect_equal(evaluation$metadata$configuration$n_draws, 12L)
  expect_equal(evaluation$metadata$configuration$replicates, 2L)
  expect_equal(evaluation$metadata$configuration$replicate_seeds, 20261215:20261216)
  expect_true(is.finite(result$summary$mean_log_score))
})
