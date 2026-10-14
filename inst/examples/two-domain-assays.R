# simCoDaR 0.1.0.9000 — hypothetical two-domain teaching example; GPL-3.
# Run: Rscript inst/examples/two-domain-assays.R [output-directory]
# With no directory, print summaries only. Requires installed simCoDaR and
# optional ggplot2, ggdist, dplyr, tidyr, purrr and tibble packages.
library(simCoDaR)

targets <- rbind(
  hematite_rich = c(Hematite = 0.65, Quartz = 0.25, Kaolinite = 0.10),
  quartz_rich = c(Hematite = 0.20, Quartz = 0.65, Kaolinite = 0.15)
)
scenario <- new_assay_scenario(
  targets, weights = c(hematite_rich = 1, quartz_rich = 1),
  logratio_sd = c(hematite_rich = 0.30, quartz_rich = 0.25)
)
x <- simulate_mineral_assay_scenario(scenario, n = 300, seed = 42)
phase <- x$observed$xrd_phase_percent
oxide <- x$observed$xrf_oxide_wt_percent
element <- x$observed$icp_element_ppm
domain <- factor(x$domains$domain, levels = rownames(targets))

cat("Hypothetical ideal chemistry; n = 300; seed = 42\n")
cat("Domain counts:\n")
print(table(domain))
cat("Mineral closure (mass percent):\n")
print(range(rowSums(phase)))
cat("Maximum mineral fraction closure error:\n")
print(max(abs(rowSums(x$truth$mineral_fraction) - 1)))
summary_data <- tibble::tibble(
  Hematite_mass_percent = phase$Hematite,
  Quartz_mass_percent = phase$Quartz,
  SiO2_wt_percent = oxide$SiO2,
  Fe_ppm = element$Fe
)
cat("Arithmetic sample means by domain (targets are geometric centers):\n")
presentation <- tibble::add_column(summary_data, sample_id = x$sample_id,
                                  domain = domain, .before = 1)
means <- dplyr::summarise(dplyr::group_by(presentation, domain),
  dplyr::across(dplyr::all_of(names(summary_data)), mean), .groups = "drop")
print(as.data.frame(means), digits = 5)
cat("Pooled Pearson correlations (closure and domain mixing affect these):\n")
print(round(cor(summary_data), 3))
cat("Within-domain Pearson correlations:\n")
print(purrr::map(split(summary_data, domain), function(z) round(cor(z), 3)))
cat("First three theoretical oxide wt% / element ppm rows:\n")
print(head(oxide[c("Fe2O3", "SiO2", "Al2O3")], 3), digits = 5)
print(head(element[c("Fe", "Si", "Al")], 3), digits = 5)

# Plot only generated package output, with each mass unit on its own scale.
# No fitted deposit, observed measurements, instrument response or error model.
draw_hero <- function() {
  colors <- c(hematite_rich = "#007C91", quartz_rich = "#C75B2A")
  labels <- c(hematite_rich = "Hematite-rich", quartz_rich = "Quartz-rich")
  feature_labels <- c(Hematite_mass_percent = "Hematite (mineral mass %)",
                      Quartz_mass_percent = "Quartz (mineral mass %)",
                      SiO2_wt_percent = "Theoretical SiO2 (wt%)",
                      Fe_ppm = "Theoretical Fe (ppm by mass)")
  long <- tidyr::pivot_longer(presentation, dplyr::all_of(names(summary_data)),
                             names_to = "feature", values_to = "value")
  long$feature <- factor(long$feature, levels = names(feature_labels),
                         labels = feature_labels)
  ggplot2::ggplot(long, ggplot2::aes(x = domain, y = value, fill = domain)) +
    ggdist::stat_halfeye(orientation = "vertical", side = "right",
      normalize = "panels", show_slab = TRUE, show_interval = FALSE,
      show_point = FALSE, trim = TRUE, scale = .5, justification = 0,
      position = ggplot2::position_nudge(x = .12), alpha = .65,
      slab_colour = NA) +
    ggplot2::geom_point(ggplot2::aes(color = domain), size = .8, alpha = .55,
      position = ggplot2::position_jitter(width = .06, height = 0,
                                          seed = 20261014), show.legend = FALSE) +
    ggplot2::facet_wrap(~ feature, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors, labels = labels, name = "Domain") +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = labels, expand = ggplot2::expansion(add = c(.3, .7))) +
    ggplot2::scale_y_continuous(labels = function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL,
      title = "simCoDaR | mineral mixtures to ideal assays",
      subtitle = "Two hypothetical domains / 300 synthetic samples / simulation seed 42",
      caption = paste("Half-clouds: smoothed marginal density, thickness normalized within each panel; separate value scales.",
        "Points: individual draws, categorical-only jitter (seed 20261014); values unchanged.",
        "Ideal fixed formulas; independent ILR draws; no geological calibration or measurement model.", sep = "\n")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top", panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9),
      plot.background = ggplot2::element_rect(fill = "#FAFAF7", color = NA),
      plot.margin = ggplot2::margin(12, 18, 12, 12))
}
render_hero <- function(path, format) {
  if (format == "png") {
    png(path, width = 2400, height = 1400, res = 200, type = "cairo")
  } else {
    svg(path, width = 12, height = 7)
  }
  on.exit(dev.off())
  withr::with_preserve_seed(print(draw_hero()))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one output directory.")
if (length(args) == 1L) {
  dir.create(args[1L], recursive = TRUE, showWarnings = FALSE)
  for (format in c("png", "svg")) {
    path <- file.path(args[1L], paste0("two-domain-assays.", format))
    render_hero(path, format)
    cat("Figure written:", path, "\n")
  }
}
