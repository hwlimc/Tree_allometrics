# ---- setup ----
current_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]

  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[[length(frame_files)]], mustWork = TRUE)))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  getwd()
}

script_dir <- current_script_dir()
project_root <- if (basename(script_dir) == "scripts") {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  normalizePath(script_dir, mustWork = TRUE)
}

setwd(project_root)
source(file.path(script_dir, "00_plotting_functions.R"))

# ---- directories ----
resolve_input_dir <- function(args, base_dir) {
  if (length(args) == 0) {
    if (file.exists(file.path(base_dir, "03_loo_distribution_comparison_by_structure.txt"))) {
      return(normalizePath(base_dir, mustWork = TRUE))
    }

    stop("No BEF distribution-comparison output found in ", base_dir)
  }

  arg <- args[[1]]
  candidates <- unique(c(
    arg,
    file.path(base_dir, arg),
    file.path(project_root, arg)
  ))
  candidates <- candidates[dir.exists(candidates)]

  if (length(candidates) == 0) {
    stop("Distribution-comparison output directory not found: ", arg)
  }

  normalizePath(candidates[[1]], mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
base_input_dir <- file.path(project_root, "processed_data", "bef_bayes_dist_comparision")
input_dir <- resolve_input_dir(args, base_input_dir)
figure_dir <- file.path(project_root, "figures", "Bayes_dist_comparison")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_output <- function(file) {
  path <- file.path(input_dir, file)

  if (!file.exists(path)) {
    stop("Missing required distribution-comparison output: ", path)
  }

  read_tsv_table(path)
}

# ---- data ----
mcmc <- read_output("01_mcmc_diagnostics.txt")
ppcheck <- read_output("02_ppcheck_observed_yrep_summary.txt")
loo_response <- read_output("03_loo_metrics_by_response.txt")
loo_structure <- read_output("03_loo_distribution_comparison_by_structure.txt")

family_levels <- intersect(
  c("gamma", "lognormal", "student", "gaussian"),
  unique(c(mcmc$family, ppcheck$family, loo_response$family, loo_structure$family))
)
family_cols <- setNames(
  c("#1B9E77", "#D95F02", "#7570B3", "#666666")[seq_along(family_levels)],
  family_levels
)

structure_levels <- unique(loo_structure$model_structure)

prepare_family_data <- function(d) {
  d$model_structure <- factor(d$model_structure, levels = structure_levels)
  d$family <- factor(d$family, levels = family_levels)
  d
}

plot_family_metric_by_structure <- function(
  d,
  file,
  metric,
  y_lab,
  main,
  hline = NULL,
  warning_col = NULL,
  rank_axis = FALSE
) {
  d <- prepare_family_data(d)
  structures <- levels(d$model_structure)
  families <- family_levels[family_levels %in% as.character(unique(d$family))]
  offsets <- setNames(seq(-0.25, 0.25, length.out = length(families)), families)
  y <- d[[metric]]
  ylim <- if (rank_axis) {
    c(max(y, na.rm = TRUE) + 0.4, 0.6)
  } else {
    plot_range(y)
  }

  with_pdf_device(file, width = 12, height = 7, {
    graphics::par(mar = c(7, 5, 4, 1))
    graphics::plot(
      NA,
      NA,
      xlim = c(0.5, length(structures) + 0.5),
      ylim = ylim,
      xaxt = "n",
      xlab = "",
      ylab = y_lab,
      main = main
    )
    graphics::axis(1, at = seq_along(structures), labels = structures, las = 2, cex.axis = 0.75)
    graphics::grid(nx = NA, ny = NULL, col = "gray88")

    if (!is.null(hline)) {
      graphics::abline(h = hline, col = "gray35", lty = 2)
    }

    for (family in families) {
      dd <- d[as.character(d$family) == family, , drop = FALSE]
      x <- match(as.character(dd$model_structure), structures) + offsets[[family]]
      warning_flag <- if (!is.null(warning_col) && warning_col %in% names(dd)) {
        as.logical(dd[[warning_col]])
      } else {
        rep(FALSE, nrow(dd))
      }

      graphics::points(
        x[!warning_flag],
        dd[[metric]][!warning_flag],
        pch = 16,
        col = family_cols[[family]],
        cex = 1.1
      )

      if (any(warning_flag, na.rm = TRUE)) {
        graphics::points(
          x[warning_flag],
          dd[[metric]][warning_flag],
          pch = 17,
          col = family_cols[[family]],
          cex = 1.2
        )
      }
    }

    graphics::legend(
      "bottomleft",
      legend = c(families, "flagged"),
      col = c(unname(family_cols[families]), "black"),
      pch = c(rep(16, length(families)), 17),
      bty = "n",
      cex = 0.85
    )
  })
}

plot_response_elpd_diff_pdf <- function(d, file) {
  d <- prepare_family_data(d)
  responses <- unique(d$response)
  structures <- levels(d$model_structure)
  families <- family_levels[family_levels %in% as.character(unique(d$family))]
  offsets <- setNames(seq(-0.25, 0.25, length.out = length(families)), families)

  with_pdf_device(file, width = 12, height = 7, {
    for (response in responses) {
      dd_resp <- d[d$response == response, , drop = FALSE]
      graphics::par(mar = c(7, 5, 4, 1))
      graphics::plot(
        NA,
        NA,
        xlim = c(0.5, length(structures) + 0.5),
        ylim = plot_range(dd_resp$elpd_diff_from_structure_response_best),
        xaxt = "n",
        xlab = "",
        ylab = "ELPD difference from best",
        main = paste("Response-level LOO distribution comparison:", response)
      )
      graphics::axis(1, at = seq_along(structures), labels = structures, las = 2, cex.axis = 0.75)
      graphics::grid(nx = NA, ny = NULL, col = "gray88")
      graphics::abline(h = 0, col = "gray35", lty = 2)

      for (family in families) {
        dd <- dd_resp[as.character(dd_resp$family) == family, , drop = FALSE]
        x <- match(as.character(dd$model_structure), structures) + offsets[[family]]
        warning_flag <- as.logical(dd$pareto_warning)
        graphics::points(x[!warning_flag], dd$elpd_diff_from_structure_response_best[!warning_flag], pch = 16, col = family_cols[[family]], cex = 1.1)
        graphics::points(x[warning_flag], dd$elpd_diff_from_structure_response_best[warning_flag], pch = 17, col = family_cols[[family]], cex = 1.2)
      }

      graphics::legend(
        "bottomleft",
        legend = c(families, "Pareto k > 0.7"),
        col = c(unname(family_cols[families]), "black"),
        pch = c(rep(16, length(families)), 17),
        bty = "n",
        cex = 0.85
      )
    }
  })
}

plot_ppcheck_rmse_pdf <- function(d, file) {
  d <- prepare_family_data(d)
  responses <- unique(d$response)
  structures <- levels(d$model_structure)
  families <- family_levels[family_levels %in% as.character(unique(d$family))]
  offsets <- setNames(seq(-0.25, 0.25, length.out = length(families)), families)

  with_pdf_device(file, width = 12, height = 7, {
    for (response in responses) {
      dd_resp <- d[d$response == response, , drop = FALSE]
      graphics::par(mar = c(7, 5, 4, 1))
      graphics::plot(
        NA,
        NA,
        xlim = c(0.5, length(structures) + 0.5),
        ylim = plot_range(dd_resp$rmse),
        xaxt = "n",
        xlab = "",
        ylab = "Posterior predictive RMSE",
        main = paste("Posterior predictive RMSE:", response)
      )
      graphics::axis(1, at = seq_along(structures), labels = structures, las = 2, cex.axis = 0.75)
      graphics::grid(nx = NA, ny = NULL, col = "gray88")

      for (family in families) {
        dd <- dd_resp[as.character(dd_resp$family) == family, , drop = FALSE]
        x <- match(as.character(dd$model_structure), structures) + offsets[[family]]
        graphics::points(x, dd$rmse, pch = 16, col = family_cols[[family]], cex = 1.1)
      }

      graphics::legend(
        "topright",
        legend = families,
        col = unname(family_cols[families]),
        pch = 16,
        bty = "n",
        cex = 0.85
      )
    }
  })
}

# ---- figure-generation ----
generated_figures <- c(
  "01_mcmc_max_rhat_by_structure.pdf",
  "01_mcmc_divergences_by_structure.pdf",
  "02_ppcheck_rmse_by_response.pdf",
  "03_loo_elpd_diff_by_structure.pdf",
  "03_loo_rank_by_structure.pdf",
  "03_loo_elpd_diff_by_response.pdf"
)

mcmc$mcmc_flag <- !as.logical(mcmc$passes_mcmc_screen)

plot_family_metric_by_structure(
  mcmc,
  file.path(figure_dir, "01_mcmc_max_rhat_by_structure.pdf"),
  metric = "max_rhat",
  y_lab = "Max Rhat",
  main = "MCMC max Rhat by model structure and distribution",
  hline = 1.01,
  warning_col = "mcmc_flag"
)

plot_family_metric_by_structure(
  mcmc,
  file.path(figure_dir, "01_mcmc_divergences_by_structure.pdf"),
  metric = "divergences",
  y_lab = "Divergences",
  main = "MCMC divergences by model structure and distribution",
  hline = 0,
  warning_col = "mcmc_flag"
)

plot_ppcheck_rmse_pdf(
  ppcheck,
  file.path(figure_dir, "02_ppcheck_rmse_by_response.pdf")
)

plot_family_metric_by_structure(
  loo_structure,
  file.path(figure_dir, "03_loo_elpd_diff_by_structure.pdf"),
  metric = "elpd_diff_from_structure_best",
  y_lab = "ELPD difference from best",
  main = "LOO distribution comparison by model structure",
  hline = 0,
  warning_col = "pareto_warning"
)

plot_family_metric_by_structure(
  loo_structure,
  file.path(figure_dir, "03_loo_rank_by_structure.pdf"),
  metric = "rank_within_structure",
  y_lab = "Rank within model structure",
  main = "LOO distribution rank by model structure",
  warning_col = "pareto_warning",
  rank_axis = TRUE
)

plot_response_elpd_diff_pdf(
  loo_response,
  file.path(figure_dir, "03_loo_elpd_diff_by_response.pdf")
)

# ---- figure-annotation ----
figure_dictionary <- data.frame(
  file = generated_figures,
  source_data = c(
    "01_mcmc_diagnostics.txt",
    "01_mcmc_diagnostics.txt",
    "02_ppcheck_observed_yrep_summary.txt",
    "03_loo_distribution_comparison_by_structure.txt",
    "03_loo_distribution_comparison_by_structure.txt",
    "03_loo_metrics_by_response.txt"
  ),
  purpose = c(
    "Max Rhat by model structure and distribution; triangles mark failed MCMC screen.",
    "Divergence counts by model structure and distribution; triangles mark failed MCMC screen.",
    "Posterior predictive RMSE by response, model structure, and distribution.",
    "Main LOO comparison: ELPD difference from the best distribution within each model structure.",
    "LOO rank of each distribution within each model structure.",
    "Response-level LOO ELPD difference from best distribution within each model structure and response."
  ),
  stringsAsFactors = FALSE
)
write_tsv_table(figure_dictionary, file.path(figure_dir, "00_figure_dictionary.txt"))

writeLines(
  c(
    "BEF Bayesian distribution-comparison figure outputs",
    "",
    paste("Generated:", as.character(Sys.time())),
    paste("Source output directory:", input_dir),
    paste("Figure output directory:", figure_dir),
    "",
    "Figure notes:",
    "- Figures are rendered from numeric .txt outputs written by scripts/02_BEF_Bayes_dist_testing.R.",
    "- Triangles indicate either failed MCMC screen or Pareto k warnings, depending on the figure.",
    "- Figure file descriptions are in 00_figure_dictionary.txt."
  ),
  file.path(figure_dir, "README_FIGURES.txt")
)

writeLines(
  c(
    paste("Figure run generated at:", as.character(Sys.time())),
    paste("Project root:", project_root),
    paste("Source output directory:", input_dir),
    paste("Figure output directory:", figure_dir),
    "",
    "Figure files:",
    generated_figures,
    "",
    "Figure metadata files:",
    c("00_figure_dictionary.txt", "README_FIGURES.txt", "FIGURE_RUN_SUMMARY.txt")
  ),
  file.path(figure_dir, "FIGURE_RUN_SUMMARY.txt")
)

figure_dir
