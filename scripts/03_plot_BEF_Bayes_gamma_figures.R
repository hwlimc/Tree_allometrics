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

# ---- output-directory ----
resolve_output_dir <- function(args, base_dir) {
  if (length(args) == 0) {
    if (file.exists(file.path(base_dir, "00_model_info.txt"))) {
      return(normalizePath(base_dir, mustWork = TRUE))
    }

    dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    dirs <- dirs[
      file.exists(file.path(dirs, "00_model_info.txt")) |
        file.exists(file.path(dirs, "00_model_manifest.txt"))
    ]

    if (length(dirs) == 0) {
      stop("No numeric BEF gamma output directories found in ", base_dir)
    }

    return(normalizePath(dirs[order(basename(dirs), decreasing = TRUE)][[1]], mustWork = TRUE))
  }

  arg <- args[[1]]
  candidates <- unique(c(
    arg,
    file.path(base_dir, arg),
    file.path(project_root, arg)
  ))
  candidates <- candidates[dir.exists(candidates)]

  if (length(candidates) == 0) {
    stop("Output directory not found: ", arg)
  }

  normalizePath(candidates[[1]], mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
base_output_dir <- file.path(project_root, "processed_data", "bef_bayes_gamma")
out_dir <- resolve_output_dir(args, base_output_dir)
figure_dir <- file.path(project_root, "figures", "Bayes_gamma")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_output <- function(file) {
  path <- file.path(out_dir, file)

  if (!file.exists(path)) {
    stop("Missing required numeric output: ", path)
  }

  read_tsv_table(path)
}

# ---- data ----

model_info <- read_output("00_model_info.txt")
step2_yrep <- read_output("02_ppcheck_obs_vs_yrep.txt")
step5_observed_fit <- read_output("05_observed_vs_predicted.txt")

make_mean_residuals_long <- function(d, mean_groups = c("h1", "h2")) {
  base_cols <- c(
    "model", "model_structure", "hierarchy", "k_depth", "family",
    "response", "row_id", "h1", "h2", "x", "observed"
  )

  out <- lapply(mean_groups, function(mean_group) {
    mean_col <- paste0("mean_observed_", mean_group)
    residual_col <- paste0("residual_from_", mean_group, "_mean")

    if (!all(c(mean_group, mean_col, residual_col) %in% names(d))) {
      return(NULL)
    }

    data.frame(
      d[, base_cols, drop = FALSE],
      mean_group = mean_group,
      mean_level = d[[mean_group]],
      mean_observed = d[[mean_col]],
      residual_from_mean = d[[residual_col]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

step5_mean_residuals <- make_mean_residuals_long(step5_observed_fit)

model_names <- model_info$model
response_names <- unique(step2_yrep$response)
preferred_figure_responses <- c("befa.st", "beft.st")
figure_response_names <- preferred_figure_responses[
  preferred_figure_responses %in% unique(step5_observed_fit$response)
]

if (length(figure_response_names) == 0) {
  figure_response_names <- unique(step5_observed_fit$response)
}

response_labels <- response_label(figure_response_names)
response_labels[figure_response_names == "befa.st"] <- "BEFA"
response_labels[figure_response_names == "beft.st"] <- "BEFT"
response_plot_styles <- make_response_plot_styles(figure_response_names, labels = response_labels)

# ---- ppc-density-from-models ----
response_info_for_fit <- function(fit) {
  data_names <- names(fit$data)

  if (!all(c("y1m1", "y2") %in% data_names)) {
    stop(
      "Gamma fits should contain response columns y1m1 and y2. Found: ",
      paste(data_names, collapse = ", ")
    )
  }

  data.frame(
    resp = c("y1m1", "y2"),
    response = c("befa.st", "befr.st"),
    transform = c("shift_y1m1", "identity"),
    stringsAsFactors = FALSE
  )
}

resp_label <- function(resp_info) {
  resp_info$response[[1]]
}

backtransform <- function(resp_info, x) {
  switch(
    resp_info$transform[[1]],
    shift_y1m1 = x + 1,
    identity = x,
    stop("Unknown response transform: ", resp_info$transform[[1]])
  )
}

plot_ppc_density_one <- function(fit, model, resp_info, ndraws = 50) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  obs <- backtransform(resp_info, dat[[resp]])

  yrep <- brms::posterior_predict(
    fit,
    newdata = dat,
    resp = resp,
    ndraws = ndraws
  )
  yrep <- backtransform(resp_info, yrep)

  obs_den <- stats::density(obs, na.rm = TRUE)
  yrep_den <- lapply(seq_len(nrow(yrep)), function(i) stats::density(yrep[i, ], na.rm = TRUE))
  y_lim <- range(c(obs_den$y, unlist(lapply(yrep_den, `[[`, "y"))), na.rm = TRUE)

  graphics::plot(
    obs_den,
    lwd = 2,
    col = "black",
    ylim = y_lim,
    xlab = resp_label(resp_info),
    main = paste(resp_label(resp_info), model, sep = " | ")
  )

  for (i in seq_along(yrep_den)) {
    graphics::lines(yrep_den[[i]], col = grDevices::rgb(1, 0, 0, 0.15))
  }

  graphics::lines(obs_den, lwd = 2, col = "black")
  graphics::legend(
    "topright",
    legend = c("Observed", "yrep"),
    col = c("black", "red"),
    lwd = 2,
    bty = "n"
  )
}

plot_ppc_density_pdf <- function(fits, model_info, file, responses, ndraws = 50) {
  layout_dim <- model_panel_layout(length(fits))

  with_pdf_device(file, width = 12, height = model_panel_pdf_height(length(fits)), {
    for (response_name in responses) {
      graphics::par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (i in seq_along(fits)) {
        resp_info <- response_info_for_fit(fits[[i]])
        resp_info <- resp_info[resp_info$response == response_name, , drop = FALSE]

        if (nrow(resp_info) == 0) {
          empty_panel(paste("No response:", response_name, model_info$model[[i]]))
          next
        }

        plot_ppc_density_one(fits[[i]], model_info$model[[i]], resp_info, ndraws = ndraws)
      }

      graphics::mtext(paste("Density PPC:", response_name, "| gamma"), outer = TRUE, cex = 1.2)
    }
  })
}

load_model_info_models <- function(model_info) {
  model_files <- model_info$model_file
  model_files <- ifelse(file.exists(model_files), model_files, file.path(project_root, model_files))
  missing_files <- model_files[!file.exists(model_files)]

  if (length(missing_files) > 0) {
    stop("Missing model files for PPC density plot:\n", paste(missing_files, collapse = "\n"))
  }

  setNames(lapply(model_files, readRDS), model_info$model)
}

# ---- figure-generation ----
generated_figures <- character(0)
skipped_figures <- character(0)

make_ppc_density <- tolower(Sys.getenv("BEF_MAKE_PPC_DENSITY", "true")) %in% c("1", "true", "yes")

if (make_ppc_density) {
  ppc_file <- "03_ppcheck_density_observed_yrep_gamma.pdf"

  if (requireNamespace("brms", quietly = TRUE)) {
    fits <- load_model_info_models(model_info)
    plot_ppc_density_pdf(
      fits,
      model_info,
      file.path(figure_dir, ppc_file),
      responses = response_names
    )
    generated_figures <- c(generated_figures, ppc_file)
  } else {
    skipped_figures <- c(skipped_figures, paste(ppc_file, "(requires the brms package)"))
  }
}

plot_observed_vs_predicted_pdf(
  step5_observed_fit,
  file.path(figure_dir, "04_observed_vs_predicted.pdf"),
  responses = unique(step5_observed_fit$response),
  models = model_names
)
generated_figures <- c(generated_figures, "04_observed_vs_predicted.pdf")

plot_grouped_residual_pdf(
  step5_observed_fit,
  file.path(figure_dir, "05_residual_vs_predicted_h1_h2_befa_beft.pdf"),
  models = model_names,
  responses = figure_response_names,
  styles = response_plot_styles,
  x_var = "predicted",
  x_lab = "Predicted"
)
generated_figures <- c(generated_figures, "05_residual_vs_predicted_h1_h2_befa_beft.pdf")

plot_grouped_residual_pdf(
  step5_observed_fit,
  file.path(figure_dir, "05_residual_vs_rsd_h1_h2_befa_beft.pdf"),
  models = model_names,
  responses = figure_response_names,
  styles = response_plot_styles,
  x_var = "x",
  x_lab = "Relative stand density"
)
generated_figures <- c(generated_figures, "05_residual_vs_rsd_h1_h2_befa_beft.pdf")

presentation_residual_files <- plot_grouped_residual_pdfs_by_model(
  step5_observed_fit,
  out_dir = figure_dir,
  file_prefix = "05_residual_diagnostics_h1_h2_befa_beft_",
  models = model_names,
  responses = figure_response_names,
  styles = response_plot_styles
)
generated_figures <- c(generated_figures, unname(presentation_residual_files))

mean_residual_groups <- intersect(c("h1", "h2"), unique(step5_mean_residuals$mean_group))

for (mean_group in mean_residual_groups) {
  file <- paste0("05_residual_vs_rsd_mean_", mean_group, "_befa_beft.pdf")
  plot_grouped_mean_residual_pdf(
    step5_mean_residuals,
    file.path(figure_dir, file),
    models = model_names,
    responses = figure_response_names,
    styles = response_plot_styles,
    mean_group = mean_group
  )
  generated_figures <- c(generated_figures, file)
}

presentation_mean_residual_files <- plot_grouped_mean_residual_pdfs_by_model(
  step5_mean_residuals,
  out_dir = figure_dir,
  file_prefix = "05_residual_vs_rsd_mean_h1_h2_befa_beft_",
  models = model_names,
  responses = figure_response_names,
  styles = response_plot_styles,
  mean_groups = mean_residual_groups
)
generated_figures <- c(generated_figures, unname(presentation_mean_residual_files))

# ---- figure-annotation ----
figure_dictionary <- data.frame(
  file = c(
    "03_ppcheck_density_observed_yrep_gamma.pdf",
    "04_observed_vs_predicted.pdf",
    "05_residual_vs_predicted_h1_h2_befa_beft.pdf",
    "05_residual_vs_rsd_h1_h2_befa_beft.pdf",
    "05_residual_diagnostics_h1_h2_befa_beft_*.pdf",
    "05_residual_vs_rsd_mean_h1_befa_beft.pdf",
    "05_residual_vs_rsd_mean_h2_befa_beft.pdf",
    "05_residual_vs_rsd_mean_h1_h2_befa_beft_*.pdf"
  ),
  source_data = c(
    "00_model_info.txt plus model .rds files",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt",
    "05_observed_vs_predicted.txt"
  ),
  purpose = c(
    "Density posterior predictive checks from observed values and new posterior predictive draws.",
    "Observed versus posterior fitted medians for calibration checks.",
    "Residual versus predicted panels by h1 with befa.st and beft.st overlaid.",
    "Residual versus relative stand density panels by h1 with befa.st and beft.st overlaid.",
    "One-model residual diagnostic PDFs for presentation.",
    "Residuals from h1 observed means versus relative stand density.",
    "Residuals from h2 observed means versus relative stand density.",
    "One-model residual-from-mean PDFs for presentation."
  ),
  stringsAsFactors = FALSE
)
write_tsv_table(figure_dictionary, file.path(figure_dir, "00_figure_dictionary.txt"))

writeLines(
  c(
    "BEF Bayesian gamma figure outputs",
    "",
    paste("Generated:", as.character(Sys.time())),
    paste("Source output directory:", out_dir),
    paste("Figure output directory:", figure_dir),
    "",
    "Figure notes:",
    "- Figures are rendered from numeric .txt outputs written by scripts/02_BEF_Bayes_gamma.R.",
    "- The PPC density figure is the only figure that uses model .rds files directly.",
    "- Set BEF_MAKE_PPC_DENSITY=false to skip the PPC density figure when only table-based figures are needed.",
    "- Figure file descriptions are in 00_figure_dictionary.txt."
  ),
  file.path(figure_dir, "README_FIGURES.txt")
)

figure_metadata_files <- c("00_figure_dictionary.txt", "README_FIGURES.txt", "FIGURE_RUN_SUMMARY.txt")

writeLines(
  c(
    paste("Figure run generated at:", as.character(Sys.time())),
    paste("Project root:", project_root),
    paste("Source output directory:", out_dir),
    paste("Figure output directory:", figure_dir),
    "",
    "Figure files:",
    sort(unique(generated_figures)),
    "",
    "Figure metadata files:",
    figure_metadata_files,
    if (length(skipped_figures) > 0) c("", "Skipped figures:", skipped_figures) else character(0)
  ),
  file.path(figure_dir, "FIGURE_RUN_SUMMARY.txt")
)

figure_dir
