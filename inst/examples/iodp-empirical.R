# Installed-package reproduction. Run with:
# Rscript <installed script> [output_directory] [full|smoke|freeze]
# Default output is a fresh temporary directory. full runs 20 fits sequentially;
# smoke runs only primary logistic_normal with 2 x 12 draws, not study evidence.
# No download, private source access, fitting retries or control changes.
run_iodp_example <- function(output_dir = tempfile("simCoDaR-iodp-"),
                             mode = c("full", "smoke", "freeze")) {
  mode <- match.arg(mode)
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) stop("output_dir must be one nonempty path.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create output directory.")
  # Refuse accidental overwrite of previously frozen evidence.
  if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("Use an empty output directory to preserve previous evidence.")
  path <- system.file("extdata", "iodp_u1418_major_trace.csv", package = "simCoDaR",
                      mustWork = TRUE)
  expected_md5 <- "14f5a0bdb8e19c18165986fb88ea9271"
  if (unname(tools::md5sum(path)) != expected_md5) stop("Bundled CSV byte identity changed.")
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  parts <- c("SiO2", "TiO2", "Al2O3", "Fe2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5")
  if (nrow(data) != 284L || ncol(data) != 32L || anyDuplicated(data$sample_id) ||
      any(!is.finite(data$depth_m)) || any(as.matrix(data[parts]) <= 0))
    stop("Pinned public schema/rows/positive oxide contract failed.")
  models <- c("dirichlet", "logistic_normal", "ilr_student_t", "ilr_mixture")
  scoped <- function(seed, code) withr::with_seed(seed, code,
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection")
  ordered <- order(data$depth_m, seq_len(nrow(data)))
  start <- scoped(20260924L, sample.int(214L, 1L))
  primary <- ordered[seq.int(start, length.out = 71L)]
  block <- floor((seq_len(284L) - 1L) * 4 / 284) + 1L
  holdouts <- c(list(primary), lapply(1:4, function(f) ordered[block == f]))
  names(holdouts) <- paste0("split", 0:4)
  controls <- list(dirichlet_maxit = 2000L, student_max_iter = 500L,
    mixture_max_iter = 300L, mixture_starts = 5L, tolerance = 1e-7)
  schedule <- expand.grid(replicate = 1:20, model_index = 1:4, split = 0:4)
  schedule$model <- models[schedule$model_index]
  schedule$seed <- 20261014L + 1000L * schedule$split +
    100L * schedule$model_index + schedule$replicate
  config <- list(evidence_label = "new isolated-package reproduction",
    original_numeric_outputs = "unavailable; historical primary start and RNG environment unavailable",
    mode = mode, source = list(doi = "10.2204/iodp.proc.341.206.2020",
      data_file = "iodp_u1418_major_trace.csv", licence = "CC BY 4.0",
      sha256 = "06b31a651aba7df1dc8b3c3f379018c4b4793fe778ce4920a4d039831eb062a3",
      md5 = expected_md5, code_source_commit = "989577aaa449d28ac554d806ef57daa47122069c",
      accepted_core_commit = "1e65abdff02672cb6d6f625124ec049c22ee2c4d"),
    example_md5 = unname(tools::md5sum(system.file("examples", "iodp-empirical.R",
      package = "simCoDaR", mustWork = TRUE))),
    package = "simCoDaR", package_version = as.character(packageVersion("simCoDaR")),
    R_version = as.character(getRversion()), withr_version = as.character(packageVersion("withr")),
    RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
    parts = parts, measurement_basis = "oxide_mass", unit_label = "wt%",
    zero_mode = "error", controls = controls, min_pattern_n = 30L,
    student_df = 5, mixture_components = 2:3, models = models,
    fit_seed = 20260924L, split_seed = 20260924L, primary_start = start,
    n_draws = 71L, replicates = 20L, diagnostic_seeds = schedule,
    bootstrap_seed = 20261914L, bootstrap_repetitions = c(400L, 0L, 0L, 0L, 0L),
    bootstrap_quantiles = c(.025, .975), quantile_type = 7L,
    smoke_configuration = list(splits = 0L, models = "logistic_normal", n_draws = 12L,
      replicates = 2L, bootstrap_repetitions = 0L),
    limitations = c("Selected oxide subcomposition, not whole-rock mass balance.",
      "Nearby core depths are dependent; iid paired bootstrap is descriptive.",
      "Depth folds partition holdouts but have overlapping training sets.",
      "MCSE is conditional draw noise, not field uncertainty.",
      "Mixture components are statistical, not inferred facies; no external ore calibration."))
  # All configuration, source identity and five manifests exist BEFORE any fit.
  saveRDS(config, file.path(output_dir, "config.rds"))
  write.csv(schedule, file.path(output_dir, "diagnostic-seeds.csv"), row.names = FALSE)
  sorted_position <- match(seq_len(284L), ordered)
  for (f in 0:4) {
    test <- holdouts[[f + 1L]]
    role <- ifelse(seq_len(284L) %in% test, "heldout", "train")
    lower <- ave(data$depth_m, role, FUN = min)
    upper <- ave(data$depth_m, role, FUN = max)
    manifest <- data.frame(split = f, original_row = seq_len(284L),
      sorted_position = sorted_position, sample_id = data$sample_id,
      depth_m = data$depth_m, role = role, role_depth_min = lower,
      role_depth_max = upper, primary_start = if (f == 0L) start else NA_integer_)
    write.csv(manifest, file.path(output_dir, paste0("split", f, "-manifest.csv")),
              row.names = FALSE)
  }
  bootstrap_indices <- scoped(20261914L, matrix(sample.int(71L, 71L * 400L,
    replace = TRUE), nrow = 71L))
  saveRDS(bootstrap_indices, file.path(output_dir, "primary-bootstrap-indices.rds"))
  if (mode == "freeze") return(invisible(config))
  summaries <- list()
  failure_rows <- list()
  paired <- list()
  claims <- list()
  splits <- if (mode == "smoke") 0L else 0:4
  model_indices <- if (mode == "smoke") 2L else 1:4
  for (f in splits) {
    test <- holdouts[[f + 1L]]
    train <- setdiff(seq_len(284L), test)
    # Original row order for all four models and pointwise scoring.
    test <- sort(test)
    score_columns <- list()
    energy_means <- numeric()
    for (m in model_indices) {
      model <- models[m]
      prefix <- paste0("split", f, "-", model)
      fit <- tryCatch(scoped(20260924L, simCoDaR::fit_assay_simulator(
        data[train, ], parts, model, seed = 20260924L, student_df = 5,
        mixture_components = 2:3, unit_label = "wt%", measurement_basis = "oxide_mass",
        zero_mode = "error", min_pattern_n = 30L, control = controls)), error = function(e) e)
      if (inherits(fit, "error")) {
        failure_rows[[length(failure_rows) + 1L]] <- data.frame(split = f, model = model,
          stage = "fit", message = conditionMessage(fit))
        saveRDS(list(status = "failure", message = conditionMessage(fit)),
                file.path(output_dir, paste0(prefix, "-fit-failure.rds")))
        next
      }
      saveRDS(fit, file.path(output_dir, paste0(prefix, "-fit.rds")))
      base_seed <- 20261014L + 1000L * f + 100L * m
      evaluation <- tryCatch(scoped(base_seed, simCoDaR::validate_assay_simulator(
        fit, data[test, ], reference = data[train, ],
        n_draws = if (mode == "smoke") 12L else 71L,
        replicates = if (mode == "smoke") 2L else 20L, seed = base_seed)), error = function(e) e)
      if (inherits(evaluation, "error")) {
        failure_rows[[length(failure_rows) + 1L]] <- data.frame(split = f, model = model,
          stage = "validation", message = conditionMessage(evaluation))
        next
      }
      saveRDS(evaluation, file.path(output_dir, paste0(prefix, "-validation.rds")))
      write.csv(evaluation$scores, file.path(output_dir, paste0(prefix, "-scores.csv")), row.names = FALSE)
      write.csv(evaluation$mc_summary, file.path(output_dir, paste0(prefix, "-mc-summary.csv")), row.names = FALSE)
      write.csv(evaluation$failures, file.path(output_dir, paste0(prefix, "-draw-failures.csv")), row.names = FALSE)
      h <- fit$groups[[1L]]$patterns[[1L]]
      p <- length(parts) - 1L
      k <- if (model == "ilr_mixture") h$model$selected_components else NA_integer_
      q <- switch(model, dirichlet = length(parts), logistic_normal = p + p * (p + 1) / 2,
        ilr_student_t = p + p * (p + 1) / 2,
        ilr_mixture = (k - 1) + k * p + k * p * (p + 1) / 2)
      ll <- if (model == "dirichlet") h$model$ilr_loglik else h$model$loglik
      summaries[[length(summaries) + 1L]] <- data.frame(split = f, model = model,
        evidence_label = if (mode == "smoke") "operational smoke only" else "new isolated-package reproduction",
        mean_log_score = evaluation$score_summary$mean_log_score,
        training_ilr_loglik = ll, training_common_measure_bic = -2 * ll + q * log(length(train)),
        selected_components = k, draw_failures = nrow(evaluation$failures))
      score_columns[[model]] <- evaluation$scores$log_score
      energy_means[model] <- subset(evaluation$mc_summary,
        metric == "ilr_energy" & statistic == "generated")$mean
    }
    if (length(score_columns) >= 2L) {
      pairs <- combn(names(score_columns), 2L)
      for (j in seq_len(ncol(pairs))) {
        a <- pairs[1L, j]
        b <- pairs[2L, j]
        delta <- score_columns[[a]] - score_columns[[b]]
        ci <- if (f == 0L && mode == "full") quantile(colMeans(matrix(delta[bootstrap_indices],
          nrow = 71L)), c(.025, .975), type = 7, names = FALSE) else c(NA_real_, NA_real_)
        paired[[length(paired) + 1L]] <- data.frame(split = f, model_a = a, model_b = b,
          mean_difference = mean(delta), descriptive_iid_lower = ci[1L],
          descriptive_iid_upper = ci[2L], bootstrap_repetitions = if (f == 0L && mode == "full") 400L else 0L)
      }
    }
    complete <- length(score_columns) == 4L
    score_mean <- vapply(score_columns, mean, numeric(1))
    claims[[length(claims) + 1L]] <- list(split = f, all_four_models_scored = complete,
      actual_score_order = names(sort(score_mean, decreasing = TRUE)),
      actual_energy_order = names(sort(energy_means)),
      source_Student_best_score_supported = if (complete) names(which.max(score_mean)) == "ilr_student_t" else NA,
      source_Dirichlet_underperformance_supported = if (complete) names(which.min(score_mean)) == "dirichlet" else NA,
      source_normal_lowest_energy_supported = if (complete && all(is.finite(energy_means)))
        names(which.min(energy_means)) == "logistic_normal" else NA,
      Dirichlet_population_correlation_limit = "Distinct-part population covariance is negative; positive finite-draw correlation can occur by chance.")
  }
  summary <- if (length(summaries)) do.call(rbind, summaries) else data.frame()
  if (nrow(summary)) summary$score_rank_among_successes <- ave(-summary$mean_log_score,
    summary$split, FUN = function(x) rank(x, ties.method = "min"))
  write.csv(summary, file.path(output_dir, "model-summary.csv"), row.names = FALSE)
  failures <- if (length(failure_rows)) do.call(rbind, failure_rows) else
    data.frame(split = integer(), model = character(), stage = character(), message = character())
  write.csv(failures, file.path(output_dir, "fit-validation-failures.csv"), row.names = FALSE)
  if (length(paired)) write.csv(do.call(rbind, paired), file.path(output_dir, "paired-score-differences.csv"), row.names = FALSE)
  saveRDS(claims, file.path(output_dir, "source-claim-comparison.rds"))
  invisible(list(config = config, summary = summary, failures = failures))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args)) args[1L] else tempfile("simCoDaR-iodp-")
  mode <- if (length(args) >= 2L) args[2L] else "full"
  run_iodp_example(output_dir, mode)
  message("Example artifacts: ", output_dir)
}
