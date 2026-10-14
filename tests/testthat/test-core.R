contract_target <- function() {
  matrix(c(8, 2, 3, 7), nrow = 2, byrow = TRUE,
         dimnames = list(c("silica", "iron"), c("Quartz", "Hematite")))
}

test_that("seeded draws are reproducible, closed and labelled", {
  spec <- new_assay_scenario(contract_target())
  a <- simulate_assay_scenario(spec, 24, seed = 41)
  b <- simulate_assay_scenario(spec, 24, seed = 41)
  c <- simulate_assay_scenario(spec, 24, seed = 42)
  expect_identical(a, b)
  expect_false(identical(a$composition, c$composition))
  expect_type(a$facies, "character")
  expect_length(a$facies, 24)
  expect_true(all(a$facies %in% rownames(contract_target())))
  expect_identical(a$sample_id, paste0("sample_", seq_len(24)))
  expect_identical(dim(a$composition), c(24L, 2L))
  expect_identical(rownames(a$composition), a$sample_id)
  expect_identical(names(a$composition), colnames(contract_target()))
  expect_true(all(vapply(a$composition, is.numeric, logical(1))))
  expect_true(all(is.finite(as.matrix(a$composition))))
  expect_true(all(as.matrix(a$composition) >= 0))
  expect_equal(rowSums(a$composition), rep(1, 24), ignore_attr = TRUE,
               tolerance = 1e-10)
  one <- simulate_assay_scenario(spec, 1, seed = 0)
  expect_identical(dim(one$composition), c(1L, 2L))
  expect_identical(one$sample_id, "sample_1")
  expect_equal(sum(one$composition), 1, tolerance = 1e-10)
})

test_that("explicit seeds restore present and absent caller RNG state", {
  withr::local_preserve_seed()
  spec <- new_assay_scenario(contract_target())
  set.seed(109)
  before <- .Random.seed
  kind <- RNGkind()
  simulate_assay_scenario(spec, 4, seed = 12)
  expect_identical(.Random.seed, before)
  simulate_mineral_assay_scenario(spec, 4, seed = 12)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)

  narrow <- new_assay_scenario(contract_target()[1, , drop = FALSE],
    composition_bounds = list(Quartz = c(0.8 - 1e-10, 0.8 + 1e-10)))
  expect_error(simulate_assay_scenario(narrow, 4, seed = 12,
    max_draw_attempts = 1), "Bounds exhausted.*silica")
  expect_identical(.Random.seed, before)
  rm(".Random.seed", envir = .GlobalEnv)
  simulate_assay_scenario(spec, 4, seed = 12)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_error(simulate_assay_scenario(narrow, 4, seed = 12,
    max_draw_attempts = 1), "Bounds exhausted")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
})

test_that("NULL seeds consume and continue the caller stream", {
  withr::local_preserve_seed()
  spec <- new_assay_scenario(contract_target())
  set.seed(105)
  before <- .Random.seed
  a <- simulate_assay_scenario(spec, 8)
  after <- .Random.seed
  b <- simulate_assay_scenario(spec, 8)
  expect_false(identical(before, after))
  expect_false(identical(after, .Random.seed))
  expect_false(identical(a$composition, b$composition))
  set.seed(105)
  expect_identical(simulate_assay_scenario(spec, 8), a)
  expect_identical(simulate_assay_scenario(spec, 8), b)
  expect_null(a$assumptions$seed)
  expect_identical(a$assumptions$rng_policy, "consume_caller_rng")
})

test_that("zero spread and structural absence respect domain targets", {
  target <- contract_target()
  spec <- new_assay_scenario(target,
    weights = c(iron = 3, silica = 1),
    logratio_sd = c(iron = 0, silica = 0),
    composition_bounds = list(Quartz = c(0.3, 0.8)))
  expect_equal(unname(spec$weights[c("silica", "iron")]), c(0.25, 0.75))
  draw <- simulate_assay_scenario(spec, 16, seed = 10)
  expected <- target[draw$facies, , drop = FALSE] / 10
  expect_equal(unname(as.matrix(draw$composition)), unname(expected))

  target["silica", ] <- c(10, 0)
  spec <- new_assay_scenario(target, structural_parts = "Hematite")
  draw <- simulate_assay_scenario(spec, 32, seed = 11)
  absent <- draw$facies == "silica"
  expect_true(any(absent))
  expect_true(any(!absent))
  expect_true(all(draw$composition$Hematite[absent] == 0))
  expect_true(all(draw$composition$Quartz[absent] == 1))
  expect_true(all(draw$composition$Hematite[!absent] > 0))
  bounded <- new_assay_scenario(contract_target(),
    composition_bounds = list(Quartz = c(0.2, 0.85)))
  draw <- simulate_assay_scenario(bounded, 32, seed = 11)
  expect_true(all(draw$composition$Quartz >= 0.2 - 1e-12))
  expect_true(all(draw$composition$Quartz <= 0.85 + 1e-12))
})

test_that("invalid counts, seeds and attempt budgets are rejected", {
  spec <- new_assay_scenario(contract_target())
  bad_counts <- list(TRUE, NA_real_, Inf, c(1, 2), 1.5, 0, -1,
                     .Machine$integer.max + 1)
  for (x in bad_counts) {
    expect_error(simulate_assay_scenario(spec, x), "n must be")
    expect_error(simulate_assay_scenario(spec, 1, max_draw_attempts = x),
                 "max_draw_attempts must be")
  }
  for (x in list(TRUE, NA_real_, Inf, c(1, 2), 1.5, -1,
                 .Machine$integer.max)) {
    expect_error(simulate_assay_scenario(spec, 1, seed = x), "seed must be")
  }
})

test_that("constructor rejects malformed and impossible scenarios", {
  target <- contract_target()
  for (x in list(NA_real_, Inf, -1)) {
    bad <- target
    bad[1, 1] <- x
    expect_error(new_assay_scenario(bad))
  }
  for (labels in list(c("", "iron"), c(NA_character_, "iron"),
                      c("iron", "iron"))) {
    bad <- target
    rownames(bad) <- labels
    expect_error(new_assay_scenario(bad))
    bad <- target
    colnames(bad) <- labels
    expect_error(new_assay_scenario(bad))
  }
  bad <- target
  bad[1, ] <- 0
  expect_error(new_assay_scenario(bad))
  bad <- target
  bad[1, 2] <- 0
  expect_error(new_assay_scenario(bad), "structural")
  for (parts in list("unknown", c("Quartz", "Quartz"), NA_character_)) {
    expect_error(new_assay_scenario(target, structural_parts = parts))
  }
  for (w in list(c(0, 1), c(1, Inf), c(1, NA), c(1, -1),
                c(other = 1, iron = 1), c(iron = 1, iron = 1))) {
    expect_error(new_assay_scenario(target, weights = w))
  }
  for (s in list(-1, Inf, NA_real_, c(other = 0, iron = 0),
                c(iron = 0, iron = 0))) {
    expect_error(new_assay_scenario(target, logratio_sd = s))
  }
  for (bounds in list(list(Quartz = c(0.9, 1)),
                      list(Quartz = c(0.8, 0.2)),
                      list(Quartz = c(-0.1, 1)),
                      list(unknown = c(0, 1)),
                      list(Quartz = c(0.8, 0.8)))) {
    expect_error(new_assay_scenario(target, composition_bounds = bounds))
  }
  expect_error(new_assay_scenario(bad, structural_parts = "Hematite",
    composition_bounds = list(Hematite = c(0.1, 1))))
})

test_that("modified and forged scenarios are revalidated at drawing", {
  spec <- new_assay_scenario(contract_target())
  expect_error(simulate_assay_scenario(structure(list(),
    class = "assay_scenario"), 1))
  for (field in c("weights", "logratio_sd", "composition_bounds")) {
    bad <- spec
    bad[[field]][] <- NA_real_
    expect_error(simulate_assay_scenario(bad, 1))
  }
  for (field in c("facies", "parts")) {
    bad <- spec
    bad[[field]][1] <- "forged"
    expect_error(simulate_assay_scenario(bad, 1))
  }
  bad <- spec
  bad$endmembers[1, 1] <- -1
  expect_error(simulate_assay_scenario(bad, 1))
  bad <- spec
  bad$unsupported_option <- TRUE
  expect_error(simulate_assay_scenario(bad, 1))
  expect_error(simulate_mineral_assay_scenario(new_assay_scenario(
    matrix(c(1, 1), 1, dimnames = list("domain", c("Quartz", "unknown")))), 1))
})

test_that("ideal quartz gives analytical oxide and elemental units", {
  quartz <- minerals_to_assay(c(Quartz = 1))
  # SiO2 formula mass from standard atomic weights, independent of library code.
  silicon_percent <- 100 * 28.0855 / (28.0855 + 2 * 15.999)
  expect_equal(quartz$oxides_physical$SiO2, 100, tolerance = 1e-8)
  expect_equal(quartz$oxides_reported$SiO2, 100, tolerance = 1e-8)
  expect_equal(quartz$elements$Si, silicon_percent, tolerance = 1e-4)
  expect_equal(quartz$elements$O, 100 - silicon_percent, tolerance = 1e-4)
  expect_equal(quartz$element_ppm$Si, 10000 * quartz$elements$Si)
  expect_equal(sum(quartz$elements), 100, tolerance = 1e-6)
  expect_equal(quartz$water$H2O_equiv_wt_pct, 0)
  expect_equal(quartz$mass_balance$physical_components_wt_pct, 100,
               tolerance = 1e-6)
})

test_that("basis conversion preserves chemistry and ideal mineral mass", {
  library <- mineral_assay_library()
  expect_true(all(c("Quartz", "Hematite", "Magnetite", "Calcite", "Pyrite") %in%
                  names(library)))
  x <- diag(length(library))
  dimnames(x) <- list(names(library), names(library))
  assay <- minerals_to_assay(x)
  expect_equal(unname(rowSums(assay$elements)), rep(100, nrow(x)),
               tolerance = 1e-6)
  expect_equal(assay$mass_balance$physical_components_wt_pct,
               rep(100, nrow(x)), tolerance = 1e-6)
  expect_equal(assay$mass_balance$element_mass_wt_pct,
               assay$mass_balance$mineral_mass_wt_pct, tolerance = 1e-6)
  expect_equal(assay$element_ppm, 10000 * assay$elements)
  percent <- minerals_to_assay(100 * x, basis = "percent")
  for (block in names(assay)[vapply(assay, is.data.frame, logical(1))]) {
    expect_equal(percent[[block]], assay[[block]])
    expect_identical(rownames(assay[[block]]), rownames(x))
    expect_true(all(is.finite(as.matrix(assay[[block]]))))
  }
  mix <- rbind(first = c(Quartz = 0.25, Magnetite = 0.75),
               second = c(Quartz = 0.9, Magnetite = 0.1))
  separate <- minerals_to_assay(mix)
  total <- minerals_to_assay(mix, iron_reporting = "total_Fe2O3")
  expect_equal(separate$elements, total$elements)
  expect_equal(separate$oxides_physical, total$oxides_physical)
  expect_true(all(total$mass_balance$reported_components_wt_pct > 100))
  expect_equal(total$mass_balance$reporting_oxygen_gain_wt_pct,
               total$mass_balance$reported_components_wt_pct - 100)
})

test_that("chemistry rejects implicit closure, invalid inputs and libraries", {
  for (x in list(c(Quartz = 0.5), c(Quartz = 100), c(Quartz = -1),
                 c(Quartz = NA_real_), c(Quartz = Inf), c(unknown = 1),
                 c(Quartz = 0.5, Quartz = 0.5), c(1), c(Quartz = "1"))) {
    expect_error(minerals_to_assay(x))
  }
  expect_error(minerals_to_assay(c(Quartz = 1), basis = "percent"))
  for (tolerance in list(-1, NA_real_, Inf, c(0, 1), TRUE)) {
    expect_error(minerals_to_assay(c(Quartz = 1), closure_tolerance = tolerance))
  }
  for (atoms in list(numeric(), c(Si = -1), c(Si = Inf), c(Si = NA_real_),
                     c(unknown = 1), c(Si = "1"), c(Si = 1, Si = 1))) {
    library <- list(Quartz = list(formula = "SiO2", atoms = atoms))
    expect_error(minerals_to_assay(c(Quartz = 1), library = library))
  }
})

test_that("integrated outputs separate observed reporting from latent truth", {
  spec <- new_assay_scenario(contract_target())
  for (n in c(1L, 12L)) {
    result <- simulate_mineral_assay_scenario(spec, n, seed = 43)
    draw <- simulate_assay_scenario(spec, n, seed = 43)
    expect_identical(result$sample_id, draw$sample_id)
    expect_identical(result$domains$sample_id, draw$sample_id)
    expect_identical(result$domains$domain, draw$facies)
    expect_identical(result$truth$mineral_fraction, draw$composition)
    expect_equal(result$observed$xrd_phase_percent, 100 * draw$composition)
    expect_identical(result$observed$xrf_oxide_wt_percent,
                     result$truth$assay$oxides_reported)
    expect_identical(result$observed$icp_element_ppm,
                     result$truth$assay$element_ppm)
    for (block in result$observed) {
      expect_s3_class(block, "data.frame")
      expect_equal(nrow(block), n)
      expect_identical(rownames(block), result$sample_id)
      expect_true(all(vapply(block, is.numeric, logical(1))))
      expect_true(all(is.finite(as.matrix(block))))
      expect_true(all(as.matrix(block) >= 0))
    }
    expect_false("mineral_fraction" %in% names(result$observed))
    expect_identical(result$metadata$schema_version, "0.1")
    expect_false(result$metadata$calibrated)
    expect_identical(result$metadata$measurement_model, "none")
    expect_true(result$metadata$truth_equals_ideal_reporting)
    expect_identical(result$metadata$units,
      c(xrd_phase_percent = "mass percent", xrf_oxide_wt_percent = "wt percent",
        icp_element_ppm = "ppm by mass"))
  }
})
